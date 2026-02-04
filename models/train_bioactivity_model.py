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
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor
from sklearn.metrics import mean_squared_error, r2_score, mean_absolute_error
import joblib

# Disable RDKit logging completely
RDLogger.DisableLog('rdApp.*')


def load_bioactivity_data(db_path, limit=50000):
    query = f"""
    SELECT DISTINCT
        md.molregno,
        cs.canonical_smiles,
        act.standard_value,
        act.standard_type,
        td.pref_name as target_name,
        td.chembl_id as target_chembl_id
    FROM molecule_dictionary md
    INNER JOIN compound_structures cs ON md.molregno = cs.molregno
    INNER JOIN activities act ON md.molregno = act.molregno
    INNER JOIN assays a ON act.assay_id = a.assay_id
    INNER JOIN target_dictionary td ON a.tid = td.tid
    WHERE act.standard_type = 'IC50'
        AND act.standard_units = 'nM'
        AND act.standard_value IS NOT NULL
        AND act.standard_value > 0
        AND a.assay_type = 'B'
        AND cs.canonical_smiles IS NOT NULL
    LIMIT {limit}
    """
    conn = sqlite3.connect(db_path)
    df = pd.read_sql_query(query, conn)
    conn.close()
    print(f"Loaded {len(df):,} bioactivity records from database")
    return df


def clean_bioactivity_data(df):
    initial_count = len(df)
    
    df = df.dropna(subset=['canonical_smiles', 'standard_value'])
    print(f"After removing nulls: {len(df):,} records")
    
    valid_smiles_mask = df['canonical_smiles'].apply(lambda x: Chem.MolFromSmiles(x) is not None)
    df = df[valid_smiles_mask].copy()
    print(f"After removing invalid SMILES: {len(df):,} records")
    
    df = df.drop_duplicates(subset=['canonical_smiles'], keep='first')
    print(f"After removing duplicates: {len(df):,} records")
    
    df['pIC50'] = -np.log10(df['standard_value'] * 1e-9)
    
    df = df[(df['pIC50'] >= 0) & (df['pIC50'] <= 15)].copy()
    print(f"After filtering pIC50 range [0, 15]: {len(df):,} records")
    
    removed = initial_count - len(df)
    print(f"Cleaned: {len(df):,} records (removed {removed:,})")
    return df


def compute_morgan_fingerprints(smiles_list, radius=2, n_bits=2048):
    fingerprints = []
    valid_indices = []
    total = len(smiles_list)
    
    for idx, smiles in enumerate(smiles_list):
        if (idx + 1) % max(1, total // 10) == 0 or idx == 0:
            print(f"  Processing: {idx + 1:,}/{total:,} molecules...")
        
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
    print(f"  ✓ Generated {len(fingerprints):,} fingerprints")
    return X, valid_indices


def train_bioactivity_model(X_train, y_train, model_type='random_forest'):
    print(f"Training {model_type} model...")
    print(f"Training set size: {len(X_train)} samples")
    
    if model_type == 'random_forest':
        model = RandomForestRegressor(
            n_estimators=100,
            random_state=42,
            n_jobs=-1,
            max_depth=20,
            min_samples_split=5,
            min_samples_leaf=2
        )
    elif model_type == 'gradient_boosting':
        model = GradientBoostingRegressor(
            n_estimators=100,
            random_state=42,
            max_depth=10,
            min_samples_split=5,
            learning_rate=0.1
        )
    
    model.fit(X_train, y_train)
    print("Model training completed")
    return model


def evaluate_model(model, X_test, y_test):
    print(f"Evaluating model on {len(X_test)} test samples...")
    y_pred = model.predict(X_test)
    
    rmse = np.sqrt(mean_squared_error(y_test, y_pred))
    mae = mean_absolute_error(y_test, y_pred)
    r2 = r2_score(y_test, y_pred)
    
    print("=" * 50)
    print("Bioactivity Model Evaluation Results")
    print("=" * 50)
    print(f"RMSE: {rmse:.4f}")
    print(f"MAE: {mae:.4f}")
    print(f"R² Score: {r2:.4f}")
    print("=" * 50)
    
    return {'rmse': rmse, 'mae': mae, 'r2': r2, 'y_pred': y_pred, 'y_true': y_test}


def save_bioactivity_model(model, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    
    model_path = os.path.join(output_dir, 'bioactivity_model.joblib')
    joblib.dump(model, model_path)
    print(f"Model saved to {model_path}")
    
    config = {
        'model_type': 'pIC50_regressor',
        'fingerprint_type': 'Morgan',
        'fingerprint_radius': 2,
        'fingerprint_bits': 2048,
        'target_column': 'pIC50'
    }
    config_path = os.path.join(output_dir, 'bioactivity_config.joblib')
    joblib.dump(config, config_path)
    print(f"Config saved to {config_path}")
    
    return model_path, config_path


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


def main():
    print("=" * 60)
    print("Bioactivity Model Training Pipeline")
    print("=" * 60)
    
    db_path = "chembl_36/chembl_36_sqlite/chembl_36.db"
    output_dir = "trained_models"
    checkpoint_path = "trained_models/.checkpoint_bioactivity.joblib"
    
    # Load checkpoint if exists
    checkpoint = load_checkpoint(checkpoint_path)
    
    if checkpoint and checkpoint.get('step') == 'complete':
        print("\n✓ Training already complete! Returning saved model...")
        model_path = os.path.join(output_dir, 'bioactivity_model.joblib')
        model = joblib.load(model_path)
        results = checkpoint.get('results', {})
        return model, results
    
    print("\n[Step 1/5] Loading bioactivity data...")
    if checkpoint and checkpoint.get('step') == 'data_loaded':
        df = checkpoint.get('df')
        print("✓ Using cached data from checkpoint")
    else:
        df = load_bioactivity_data(db_path)
        save_checkpoint(checkpoint_path, 'data_loaded', {'df': df})
    
    print("\n[Step 2/5] Cleaning and balancing dataset...")
    if checkpoint and checkpoint.get('step') == 'data_cleaned':
        df_clean = checkpoint.get('df_clean')
        print("✓ Using cleaned data from checkpoint")
    else:
        df_clean = clean_bioactivity_data(df)
        # Balance the continuous distribution
        print("Balancing pIC50 distribution...")
        df_clean['bin'] = pd.cut(df_clean['pIC50'], bins=10, labels=False)
        min_count = df_clean['bin'].value_counts().min()
        median_count = df_clean['bin'].value_counts().median()
        cap = int(max(median_count, 1000))
        
        balanced_dfs = []
        for bin_id in df_clean['bin'].unique():
            bin_df = df_clean[df_clean['bin'] == bin_id]
            if len(bin_df) > cap:
                bin_df = bin_df.sample(n=cap, random_state=42)
            balanced_dfs.append(bin_df)
        
        df_clean = pd.concat(balanced_dfs).sample(frac=1, random_state=42).reset_index(drop=True)
        print(f"Balanced dataset size: {len(df_clean):,}")
        
        print(f"\n✓ Checkpoint saved: data_cleaned")
        save_checkpoint(checkpoint_path, 'data_cleaned', {'df_clean': df_clean})
    
    print("\n[Step 3/5] Computing molecular fingerprints...")
    if checkpoint and checkpoint.get('step') == 'fingerprints_computed':
        X = checkpoint.get('X')
        y = checkpoint.get('y')
        print("✓ Using cached fingerprints from checkpoint")
    else:
        smiles_list = df_clean['canonical_smiles'].tolist()
        X, valid_indices = compute_morgan_fingerprints(smiles_list)
        y = df_clean['pIC50'].iloc[valid_indices].values
        
        print(f"\n✓ Checkpoint saved: fingerprints_computed")
        save_checkpoint(checkpoint_path, 'fingerprints_computed', {'X': X, 'y': y})
    
    print(f"Final dataset: {len(X):,} molecules")
    
    print("\n[Step 4/5] Splitting and training...")
    if checkpoint and checkpoint.get('step') == 'model_trained':
        model = checkpoint.get('model')
        print("✓ Using cached model from checkpoint")
    else:
        X_train, X_test, y_train, y_test = train_test_split(
            X, y, test_size=0.2, random_state=42
        )
        print(f"Training: {len(X_train)}, Test: {len(X_test)}")
        model = train_bioactivity_model(X_train, y_train)
        save_checkpoint(checkpoint_path, 'model_trained', {'model': model, 'X_test': X_test, 'y_test': y_test})
    
    # Get test data from checkpoint if needed
    if not checkpoint or checkpoint.get('step') != 'model_trained':
        X_train, X_test, y_train, y_test = train_test_split(
            X, y, test_size=0.2, random_state=42
        )
    else:
        X_test = checkpoint.get('X_test')
        y_test = checkpoint.get('y_test')
    
    print("\n[Step 5/5] Evaluating model...")
    results = evaluate_model(model, X_test, y_test)
    save_bioactivity_model(model, output_dir)
    save_checkpoint(checkpoint_path, 'complete', {'results': results})
    
    return model, results
    
    print("\n[Step 5/5] Evaluating and saving...")
    results = evaluate_model(model, X_test, y_test)
    save_bioactivity_model(model, output_dir)
    
    print("\nBioactivity model training completed!")
    return model, results


if __name__ == "__main__":
    model, results = main()
