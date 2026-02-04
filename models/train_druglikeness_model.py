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
from rdkit.Chem import AllChem, Descriptors, QED
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor, RandomForestClassifier
from sklearn.metrics import mean_squared_error, r2_score, accuracy_score, classification_report
import joblib

# Disable RDKit logging completely
RDLogger.DisableLog('rdApp.*')


def load_druglikeness_data(db_path):
    query = """
    SELECT 
        cs.molregno,
        cs.canonical_smiles,
        cp.mw_freebase,
        cp.alogp,
        cp.hba,
        cp.hbd,
        cp.psa,
        cp.rtb,
        cp.num_ro5_violations,
        cp.qed_weighted,
        md.max_phase,
        md.therapeutic_flag,
        md.oral,
        md.black_box_warning
    FROM compound_structures cs
    INNER JOIN compound_properties cp ON cs.molregno = cp.molregno
    INNER JOIN molecule_dictionary md ON cs.molregno = md.molregno
    WHERE cs.canonical_smiles IS NOT NULL
        AND cp.mw_freebase IS NOT NULL
        AND cp.qed_weighted IS NOT NULL
    LIMIT 500000
    """
    conn = sqlite3.connect(db_path)
    df = pd.read_sql_query(query, conn)
    conn.close()
    print(f"Loaded {len(df)} compounds with drug-likeness data")
    return df


def clean_druglikeness_data(df):
    initial_count = len(df)
    
    df = df.dropna(subset=['canonical_smiles', 'qed_weighted'])
    
    valid_smiles = []
    for idx, row in df.iterrows():
        mol = Chem.MolFromSmiles(row['canonical_smiles'])
        if mol is not None:
            valid_smiles.append(idx)
    
    df = df.loc[valid_smiles].copy()
    print(f"Valid molecules: {len(df)}")
    
    df = df.drop_duplicates(subset=['canonical_smiles'], keep='first')
    print(f"After deduplication: {len(df)}")
    
    df['is_druglike'] = (
        (df['num_ro5_violations'] <= 1) &
        (df['qed_weighted'] >= 0.3)
    ).astype(int)
    
    df['is_approved'] = (df['max_phase'] == 4).astype(int)
    
    print(f"Drug-like molecules: {df['is_druglike'].sum()} ({100*df['is_druglike'].mean():.1f}%)")
    print(f"Approved drugs: {df['is_approved'].sum()}")
    print(f"Cleaned: {len(df)} records (removed {initial_count - len(df)})")
    
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


def train_qed_regressor(X_train, y_train):
    print("Training QED regressor...")
    
    model = RandomForestRegressor(
        n_estimators=100,
        random_state=42,
        n_jobs=-1,
        max_depth=15
    )
    model.fit(X_train, y_train)
    
    print("QED regressor training completed")
    return model


def train_druglike_classifier(X_train, y_train):
    print("Training drug-likeness classifier...")
    
    model = RandomForestClassifier(
        n_estimators=100,
        random_state=42,
        n_jobs=-1,
        max_depth=15,
        class_weight='balanced'
    )
    model.fit(X_train, y_train)
    
    print("Drug-likeness classifier training completed")
    return model


def evaluate_models(qed_model, druglike_model, X_test, y_qed_test, y_druglike_test):
    print("=" * 60)
    print("Drug-likeness Model Evaluation")
    print("=" * 60)
    
    y_qed_pred = qed_model.predict(X_test)
    rmse = np.sqrt(mean_squared_error(y_qed_test, y_qed_pred))
    r2 = r2_score(y_qed_test, y_qed_pred)
    print(f"\nQED Regressor:")
    print(f"  RMSE: {rmse:.4f}")
    print(f"  R²: {r2:.4f}")
    
    y_druglike_pred = druglike_model.predict(X_test)
    accuracy = accuracy_score(y_druglike_test, y_druglike_pred)
    print(f"\nDrug-likeness Classifier:")
    print(f"  Accuracy: {accuracy:.4f}")
    print(classification_report(y_druglike_test, y_druglike_pred, 
                               target_names=['Not Drug-like', 'Drug-like']))
    print("=" * 60)
    
    return {
        'qed_rmse': rmse,
        'qed_r2': r2,
        'druglike_accuracy': accuracy
    }


def calculate_lipinski_violations(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    
    violations = 0
    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    hbd = Descriptors.NumHDonors(mol)
    hba = Descriptors.NumHAcceptors(mol)
    
    if mw > 500:
        violations += 1
    if logp > 5:
        violations += 1
    if hbd > 5:
        violations += 1
    if hba > 10:
        violations += 1
    
    return {
        'mw': mw,
        'logp': logp,
        'hbd': hbd,
        'hba': hba,
        'violations': violations,
        'is_druglike': violations <= 1
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


def save_druglikeness_models(qed_model, druglike_model, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    
    qed_path = os.path.join(output_dir, 'qed_model.joblib')
    joblib.dump(qed_model, qed_path)
    print(f"QED model saved to {qed_path}")
    
    druglike_path = os.path.join(output_dir, 'druglikeness_model.joblib')
    joblib.dump(druglike_model, druglike_path)
    print(f"Drug-likeness model saved to {druglike_path}")
    
    config = {
        'qed_model_type': 'regressor',
        'druglike_model_type': 'classifier',
        'fingerprint_type': 'Morgan',
        'fingerprint_radius': 2,
        'fingerprint_bits': 1024,
        'lipinski_rules': {
            'mw_max': 500,
            'logp_max': 5,
            'hbd_max': 5,
            'hba_max': 10,
            'max_violations': 1
        }
    }
    config_path = os.path.join(output_dir, 'druglikeness_config.joblib')
    joblib.dump(config, config_path)
    print(f"Config saved to {config_path}")
    
    return qed_path, druglike_path


def main():
    print("=" * 60)
    print("Drug-likeness Model Training Pipeline")
    print("=" * 60)
    
    db_path = "chembl_36/chembl_36_sqlite/chembl_36.db"
    output_dir = "trained_models"
    checkpoint_path = "trained_models/.checkpoint_druglikeness.joblib"
    
    # Load checkpoint if exists
    checkpoint = load_checkpoint(checkpoint_path)
    
    if checkpoint and checkpoint.get('step') == 'complete':
        print("\n✓ Training already complete! Returning saved models...")
        qed_model = checkpoint.get('qed_model')
        druglike_model = checkpoint.get('druglike_model')
        results = checkpoint.get('results', {})
        return qed_model, druglike_model, results
    
    print("\n[Step 1/5] Loading drug-likeness data...")
    if checkpoint and checkpoint.get('step') == 'data_loaded':
        df = checkpoint.get('df')
        print("✓ Using cached data from checkpoint")
    else:
        df = load_druglikeness_data(db_path)
        save_checkpoint(checkpoint_path, 'data_loaded', {'df': df})
    
    print("\n[Step 2/5] Cleaning dataset...")
    if checkpoint and checkpoint.get('step') == 'data_cleaned':
        df_clean = checkpoint.get('df_clean')
        print("✓ Using cleaned data from checkpoint")
    else:
        df_clean = clean_druglikeness_data(df)
        save_checkpoint(checkpoint_path, 'data_cleaned', {'df_clean': df_clean})
    
    print("\n[Step 3/5] Computing fingerprints...")
    if checkpoint and checkpoint.get('step') == 'fingerprints_computed':
        X = checkpoint.get('X')
        y_qed = checkpoint.get('y_qed')
        y_druglike = checkpoint.get('y_druglike')
        print("✓ Using cached fingerprints from checkpoint")
    else:
        smiles_list = df_clean['canonical_smiles'].tolist()
        X, valid_indices = compute_fingerprints(smiles_list)
        
        df_valid = df_clean.iloc[valid_indices].copy()
        y_qed = df_valid['qed_weighted'].values
        y_druglike = df_valid['is_druglike'].values
        
        save_checkpoint(checkpoint_path, 'fingerprints_computed', {'X': X, 'y_qed': y_qed, 'y_druglike': y_druglike})
    
    print(f"Final dataset: {len(X)} molecules")
    
    print("\n[Step 4/5] Training models...")
    if checkpoint and checkpoint.get('step') == 'models_trained':
        qed_model = checkpoint.get('qed_model')
        druglike_model = checkpoint.get('druglike_model')
        X_test = checkpoint.get('X_test')
        y_qed_test = checkpoint.get('y_qed_test')
        y_dl_test = checkpoint.get('y_dl_test')
        print("✓ Using cached models from checkpoint")
    else:
        X_train, X_test, y_qed_train, y_qed_test, y_dl_train, y_dl_test = train_test_split(
            X, y_qed, y_druglike, test_size=0.2, random_state=42
        )
        print(f"Training: {len(X_train)}, Test: {len(X_test)}")
        
        qed_model = train_qed_regressor(X_train, y_qed_train)
        druglike_model = train_druglike_classifier(X_train, y_dl_train)
        
        save_checkpoint(checkpoint_path, 'models_trained', {
            'qed_model': qed_model,
            'druglike_model': druglike_model,
            'X_test': X_test,
            'y_qed_test': y_qed_test,
            'y_dl_test': y_dl_test
        })
    
    print("\n[Step 5/5] Evaluating and saving...")
    results = evaluate_models(qed_model, druglike_model, X_test, y_qed_test, y_dl_test)
    save_druglikeness_models(qed_model, druglike_model, output_dir)
    save_checkpoint(checkpoint_path, 'complete', {'qed_model': qed_model, 'druglike_model': druglike_model, 'results': results})
    
    print("\nDrug-likeness model training completed!")
    return qed_model, druglike_model, results


if __name__ == "__main__":
    qed_model, druglike_model, results = main()
