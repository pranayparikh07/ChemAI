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
from rdkit.Chem import AllChem, Descriptors, Lipinski
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor
from sklearn.multioutput import MultiOutputRegressor
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.preprocessing import StandardScaler
import joblib

# Disable RDKit logging completely
RDLogger.DisableLog('rdApp.*')


def load_property_data(db_path):
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
        cp.aromatic_rings,
        cp.heavy_atoms,
        cp.qed_weighted,
        cp.num_ro5_violations
    FROM compound_structures cs
    INNER JOIN compound_properties cp ON cs.molregno = cp.molregno
    WHERE cs.canonical_smiles IS NOT NULL
        AND cp.mw_freebase IS NOT NULL
        AND cp.alogp IS NOT NULL
    LIMIT 500000
    """
    conn = sqlite3.connect(db_path)
    df = pd.read_sql_query(query, conn)
    conn.close()
    print(f"Loaded {len(df)} compounds with properties")
    return df


def clean_property_data(df):
    initial_count = len(df)
    
    df = df.dropna()
    print(f"After removing nulls: {len(df)} records")
    
    valid_smiles = []
    for idx, row in df.iterrows():
        mol = Chem.MolFromSmiles(row['canonical_smiles'])
        if mol is not None:
            valid_smiles.append(idx)
    
    df = df.loc[valid_smiles].copy()
    print(f"After removing invalid SMILES: {len(df)} records")
    
    df = df.drop_duplicates(subset=['canonical_smiles'], keep='first')
    print(f"After removing duplicates: {len(df)} records")
    
    df = df[
        (df['mw_freebase'] > 0) & (df['mw_freebase'] < 1000) &
        (df['alogp'] > -10) & (df['alogp'] < 10) &
        (df['psa'] >= 0) & (df['psa'] < 500)
    ].copy()
    print(f"After filtering outliers: {len(df)} records")
    
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


def train_property_models(X_train, y_train, property_names):
    print(f"Training multi-output property predictor...")
    print(f"Properties: {property_names}")
    print(f"Training samples: {len(X_train)}")
    
    base_model = RandomForestRegressor(
        n_estimators=100,
        random_state=42,
        n_jobs=-1,
        max_depth=15,
        min_samples_split=5
    )
    
    model = MultiOutputRegressor(base_model)
    model.fit(X_train, y_train)
    
    print("Property model training completed")
    return model


def evaluate_property_models(model, X_test, y_test, property_names):
    print(f"Evaluating on {len(X_test)} test samples...")
    y_pred = model.predict(X_test)
    
    results = {}
    print("=" * 60)
    print("Property Prediction Results")
    print("=" * 60)
    print(f"{'Property':<20} {'RMSE':<12} {'R²':<12}")
    print("-" * 60)
    
    for i, prop in enumerate(property_names):
        rmse = np.sqrt(mean_squared_error(y_test[:, i], y_pred[:, i]))
        r2 = r2_score(y_test[:, i], y_pred[:, i])
        results[prop] = {'rmse': rmse, 'r2': r2}
        print(f"{prop:<20} {rmse:<12.4f} {r2:<12.4f}")
    
    print("=" * 60)
    return results


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


def save_property_model(model, scaler, property_names, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    
    model_path = os.path.join(output_dir, 'property_model.joblib')
    joblib.dump(model, model_path)
    print(f"Model saved to {model_path}")
    
    scaler_path = os.path.join(output_dir, 'property_scaler.joblib')
    joblib.dump(scaler, scaler_path)
    print(f"Scaler saved to {scaler_path}")
    
    config = {
        'model_type': 'multi_property_regressor',
        'property_names': property_names,
        'fingerprint_type': 'Morgan',
        'fingerprint_radius': 2,
        'fingerprint_bits': 1024
    }
    config_path = os.path.join(output_dir, 'property_config.joblib')
    joblib.dump(config, config_path)
    print(f"Config saved to {config_path}")
    
    return model_path


def main():
    print("=" * 60)
    print("Property Prediction Model Training Pipeline")
    print("=" * 60)
    
    db_path = "chembl_36/chembl_36_sqlite/chembl_36.db"
    output_dir = "trained_models"
    checkpoint_path = "trained_models/.checkpoint_property.joblib"
    
    property_columns = ['mw_freebase', 'alogp', 'hba', 'hbd', 'psa', 'rtb', 'qed_weighted']
    
    # Load checkpoint if exists
    checkpoint = load_checkpoint(checkpoint_path)
    
    if checkpoint and checkpoint.get('step') == 'complete':
        print("\n✓ Training already complete! Returning saved model...")
        model_path = os.path.join(output_dir, 'property_model.joblib')
        model = joblib.load(model_path)
        results = checkpoint.get('results', {})
        return model, results
    
    print("\n[Step 1/5] Loading property data...")
    if checkpoint and checkpoint.get('step') == 'data_loaded':
        df = checkpoint.get('df')
        print("✓ Using cached data from checkpoint")
    else:
        df = load_property_data(db_path)
        save_checkpoint(checkpoint_path, 'data_loaded', {'df': df})
    
    print("\n[Step 2/5] Cleaning dataset...")
    if checkpoint and checkpoint.get('step') == 'data_cleaned':
        df_clean = checkpoint.get('df_clean')
        print("✓ Using cleaned data from checkpoint")
    else:
        df_clean = clean_property_data(df)
        save_checkpoint(checkpoint_path, 'data_cleaned', {'df_clean': df_clean})
    
    print("\n[Step 3/5] Computing fingerprints...")
    if checkpoint and checkpoint.get('step') == 'fingerprints_computed':
        X = checkpoint.get('X')
        y = checkpoint.get('y')
        y_scaled = checkpoint.get('y_scaled')
        scaler = checkpoint.get('scaler')
        print("✓ Using cached fingerprints from checkpoint")
    else:
        smiles_list = df_clean['canonical_smiles'].tolist()
        X, valid_indices = compute_fingerprints(smiles_list)
        
        df_valid = df_clean.iloc[valid_indices].copy()
        y = df_valid[property_columns].values
        
        scaler = StandardScaler()
        y_scaled = scaler.fit_transform(y)
        
        save_checkpoint(checkpoint_path, 'fingerprints_computed', {'X': X, 'y': y, 'y_scaled': y_scaled, 'scaler': scaler})
    
    print(f"Final dataset: {len(X)} molecules, {len(property_columns)} properties")
    
    print("\n[Step 4/5] Training property models...")
    if checkpoint and checkpoint.get('step') == 'model_trained':
        model = checkpoint.get('model')
        X_test = checkpoint.get('X_test')
        y_test = checkpoint.get('y_test')
        print("✓ Using cached model from checkpoint")
    else:
        X_train, X_test, y_train, y_test = train_test_split(
            X, y_scaled, test_size=0.2, random_state=42
        )
        print(f"Training: {len(X_train)}, Test: {len(X_test)}")
        model = train_property_models(X_train, y_train, property_columns)
        save_checkpoint(checkpoint_path, 'model_trained', {'model': model, 'X_test': X_test, 'y_test': y_test, 'scaler': scaler})
    
    print("\n[Step 5/5] Evaluating and saving...")
    y_test_orig = scaler.inverse_transform(y_test)
    y_pred_scaled = model.predict(X_test)
    y_pred_orig = scaler.inverse_transform(y_pred_scaled)
    
    results = {}
    print("=" * 60)
    print("Property Prediction Results (Original Scale)")
    print("=" * 60)
    print(f"{'Property':<20} {'RMSE':<12} {'R²':<12}")
    print("-" * 60)
    
    for i, prop in enumerate(property_columns):
        rmse = np.sqrt(mean_squared_error(y_test_orig[:, i], y_pred_orig[:, i]))
        r2 = r2_score(y_test_orig[:, i], y_pred_orig[:, i])
        results[prop] = {'rmse': rmse, 'r2': r2}
        print(f"{prop:<20} {rmse:<12.4f} {r2:<12.4f}")
    print("=" * 60)
    
    save_property_model(model, scaler, property_columns, output_dir)
    save_checkpoint(checkpoint_path, 'complete', {'results': results})
    
    print("\nProperty model training completed!")
    return model, results


if __name__ == "__main__":
    model, results = main()
