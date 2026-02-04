import warnings
import sys
import os
import io
import sqlite3
import numpy as np
import pandas as pd

# Suppress RDKit deprecation warnings from stderr
warnings.filterwarnings('ignore')
os.environ['PYTHONWARNINGS'] = 'ignore'

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report, roc_auc_score, accuracy_score
import joblib

# Disable RDKit logging completely
RDLogger.DisableLog('rdApp.*')


def load_structural_alerts(db_path):
    query = """
    SELECT 
        sa.alert_id,
        sa.alert_name,
        sa.smarts,
        sas.set_name,
        sas.priority
    FROM structural_alerts sa
    INNER JOIN structural_alert_sets sas ON sa.alert_set_id = sas.alert_set_id
    """
    conn = sqlite3.connect(db_path)
    alerts_df = pd.read_sql_query(query, conn)
    conn.close()
    print(f"Loaded {len(alerts_df)} structural alerts")
    return alerts_df


def load_molecules_with_alerts(db_path):
    query = """
    SELECT 
        cs.molregno,
        cs.canonical_smiles,
        COUNT(DISTINCT csa.alert_id) as alert_count,
        GROUP_CONCAT(DISTINCT sa.alert_name) as alert_names
    FROM compound_structures cs
    LEFT JOIN compound_structural_alerts csa ON cs.molregno = csa.molregno
    LEFT JOIN structural_alerts sa ON csa.alert_id = sa.alert_id
    WHERE cs.canonical_smiles IS NOT NULL
    GROUP BY cs.molregno, cs.canonical_smiles
    LIMIT 300000
    """
    conn = sqlite3.connect(db_path)
    df = pd.read_sql_query(query, conn)
    conn.close()
    
    df['has_alerts'] = (df['alert_count'] > 0).astype(int)
    
    print(f"Loaded {len(df)} molecules")
    print(f"With alerts: {df['has_alerts'].sum()} ({100*df['has_alerts'].mean():.1f}%)")
    print(f"Without alerts: {(1-df['has_alerts']).sum()} ({100*(1-df['has_alerts'].mean()):.1f}%)")
    
    return df


def load_drug_warnings(db_path):
    query = """
    SELECT DISTINCT
        cs.molregno,
        cs.canonical_smiles,
        dw.warning_type,
        dw.warning_class
    FROM compound_structures cs
    INNER JOIN drug_warning dw ON cs.molregno = dw.molregno
    WHERE cs.canonical_smiles IS NOT NULL
    """
    conn = sqlite3.connect(db_path)
    df = pd.read_sql_query(query, conn)
    conn.close()
    print(f"Loaded {len(df)} drug warning records")
    return df


def clean_data(df):
    initial_count = len(df)
    
    df = df.dropna(subset=['canonical_smiles'])
    
    valid_smiles = []
    for idx, row in df.iterrows():
        mol = Chem.MolFromSmiles(row['canonical_smiles'])
        if mol is not None:
            valid_smiles.append(idx)
    
    df = df.loc[valid_smiles].copy()
    print(f"Valid molecules: {len(df)} (removed {initial_count - len(df)})")
    
    return df


def compute_fingerprints(smiles_list, radius=2, n_bits=1024):
    fingerprints = []
    valid_indices = []
    
    for idx, smiles in enumerate(smiles_list):
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            # Redirect stderr to suppress C++ warnings
            old_stderr = sys.stderr
            sys.stderr = io.StringIO()
            try:
                with warnings.catch_warnings():
                    warnings.simplefilter('ignore')
                    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=n_bits)
                fingerprints.append(np.array(fp))
                valid_indices.append(idx)
            finally:
                sys.stderr = old_stderr
    
    X = np.array(fingerprints)
    print(f"Generated {len(fingerprints)} fingerprints")
    return X, valid_indices


def train_toxicity_classifier(X_train, y_train):
    print(f"Training toxicity classifier...")
    print(f"Training samples: {len(X_train)}")
    print(f"Positive class: {y_train.sum()} ({100*y_train.mean():.1f}%)")
    
    model = RandomForestClassifier(
        n_estimators=100,
        random_state=42,
        n_jobs=-1,
        max_depth=15,
        min_samples_split=5,
        class_weight='balanced'
    )
    
    model.fit(X_train, y_train)
    print("Toxicity classifier training completed")
    return model


def evaluate_toxicity_model(model, X_test, y_test):
    print(f"Evaluating on {len(X_test)} test samples...")
    
    y_pred = model.predict(X_test)
    y_prob = model.predict_proba(X_test)[:, 1]
    
    accuracy = accuracy_score(y_test, y_pred)
    roc_auc = roc_auc_score(y_test, y_prob)
    
    print("=" * 60)
    print("Toxicity Classifier Results")
    print("=" * 60)
    print(f"Accuracy: {accuracy:.4f}")
    print(f"ROC-AUC: {roc_auc:.4f}")
    print("-" * 60)
    print(classification_report(y_test, y_pred, target_names=['Safe', 'Toxic']))
    print("=" * 60)
    
    return {
        'accuracy': accuracy,
        'roc_auc': roc_auc,
        'y_pred': y_pred,
        'y_prob': y_prob
    }


def load_checkpoint(checkpoint_path):
    """Load training checkpoint if it exists"""
    if os.path.exists(checkpoint_path):
        checkpoint = joblib.load(checkpoint_path)
        print(f"\n✓ Checkpoint loaded from {checkpoint_path}")
        print(f"  Resuming from step: {checkpoint.get('step', 'unknown')}")
        return checkpoint
    return None


def save_checkpoint(checkpoint_path, step, data):
    """Save training checkpoint"""
    os.makedirs(os.path.dirname(checkpoint_path) or '.', exist_ok=True)
    joblib.dump({'step': step, **data}, checkpoint_path)
    print(f"\n✓ Checkpoint saved: {step}")


def save_toxicity_model(model, alerts_df, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    
    model_path = os.path.join(output_dir, 'toxicity_model.joblib')
    joblib.dump(model, model_path)
    print(f"Model saved to {model_path}")
    
    alerts_path = os.path.join(output_dir, 'structural_alerts.joblib')
    alerts_data = alerts_df[['alert_name', 'smarts', 'set_name']].to_dict('records')
    joblib.dump(alerts_data, alerts_path)
    print(f"Structural alerts saved to {alerts_path}")
    
    config = {
        'model_type': 'toxicity_classifier',
        'fingerprint_type': 'Morgan',
        'fingerprint_radius': 2,
        'fingerprint_bits': 1024,
        'num_alerts': len(alerts_df)
    }
    config_path = os.path.join(output_dir, 'toxicity_config.joblib')
    joblib.dump(config, config_path)
    print(f"Config saved to {config_path}")
    
    return model_path


def main():
    print("=" * 60)
    print("Toxicity Prediction Model Training Pipeline")
    print("=" * 60)
    
    db_path = "chembl_36/chembl_36_sqlite/chembl_36.db"
    output_dir = "trained_models"
    
    print("\n[Step 1/5] Loading structural alerts...")
    alerts_df = load_structural_alerts(db_path)
    
    print("\n[Step 2/5] Loading molecules with alert annotations...")
    df = load_molecules_with_alerts(db_path)
    
    print("\n[Step 3/5] Cleaning and balancing dataset...")
    df_clean = clean_data(df)

    # Balance the dataset (1:1 ratio)
    print("Balancing dataset (1:1 ratio)...")
    toxic_df = df_clean[df_clean['has_alerts'] == 1]
    safe_df = df_clean[df_clean['has_alerts'] == 0]
    
    n_toxic = len(toxic_df)
    n_safe = len(safe_df)
    
    if n_safe > n_toxic:
        safe_df = safe_df.sample(n=n_toxic, random_state=42)
        print(f"Downsampled safe molecules from {n_safe} to {n_toxic}")
    else:
        print("Dataset already balanced or toxic > safe")
        
    df_clean = pd.concat([toxic_df, safe_df]).sample(frac=1, random_state=42).reset_index(drop=True)
    print(f"Balanced dataset size: {len(df_clean)}")
    
    print("\n[Step 4/5] Computing fingerprints...")
    smiles_list = df_clean['canonical_smiles'].tolist()
    X, valid_indices = compute_fingerprints(smiles_list)
    y = df_clean['has_alerts'].iloc[valid_indices].values
    
    print(f"Final dataset: {len(X)} molecules")
    
    print("\n[Step 5/5] Training and evaluating...")
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.2, random_state=42, stratify=y
    )
    print(f"Training: {len(X_train)}, Test: {len(X_test)}")
    
    model = train_toxicity_classifier(X_train, y_train)
    results = evaluate_toxicity_model(model, X_test, y_test)
    
    save_toxicity_model(model, alerts_df, output_dir)
    
    print("\nToxicity model training completed!")
    return model, alerts_df, results


if __name__ == "__main__":
    model, alerts_df, results = main()
