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


# Batch processing configuration - adjust based on available RAM
BATCH_SIZE = 50000  # Configurable batch size for memory management


def load_property_data(db_path, limit=None):
    query = f"""
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
    """
    if limit:
        query += f" LIMIT {limit}"
    
    conn = sqlite3.connect(db_path)
    df = pd.read_sql_query(query, conn)
    conn.close()
    print(f"Loaded {len(df):,} compounds with properties")
    return df


def clean_property_data(df):
    initial_count = len(df)
    
    df = df.dropna()
    print(f"  After removing nulls: {len(df):,} records")
    
    print("  Validating SMILES structures...")
    valid_smiles = []
    for idx, row in df.iterrows():
        mol = Chem.MolFromSmiles(row['canonical_smiles'])
        if mol is not None:
            valid_smiles.append(idx)
    
    df = df.loc[valid_smiles].copy()
    print(f"  After removing invalid SMILES: {len(df):,} records")
    
    df = df.drop_duplicates(subset=['canonical_smiles'], keep='first')
    print(f"  After removing duplicates: {len(df):,} records")
    
    df = df[
        (df['mw_freebase'] > 0) & (df['mw_freebase'] < 1000) &
        (df['alogp'] > -10) & (df['alogp'] < 10) &
        (df['psa'] >= 0) & (df['psa'] < 500)
    ].copy()
    print(f"  After filtering outliers: {len(df):,} records")
    
    print(f"Cleaned: {len(df):,} records (removed {initial_count - len(df):,})")
    return df


def compute_fingerprints_batch(smiles_list, batch_size=BATCH_SIZE, radius=2, n_bits=1024):
    """Compute fingerprints in batches to manage memory"""
    all_fingerprints = []
    all_valid_indices = []
    total = len(smiles_list)
    
    num_batches = (total + batch_size - 1) // batch_size
    
    for batch_idx in range(num_batches):
        start_idx = batch_idx * batch_size
        end_idx = min((batch_idx + 1) * batch_size, total)
        batch_smiles = smiles_list[start_idx:end_idx]
        
        print(f"  Processing batch {batch_idx + 1}/{num_batches} ({start_idx:,}-{end_idx:,})...")
        
        fingerprints = []
        valid_indices = []
        
        for local_idx, smiles in enumerate(batch_smiles):
            global_idx = start_idx + local_idx
            
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
                    valid_indices.append(global_idx)
                finally:
                    sys.stderr = old_stderr
        
        all_fingerprints.extend(fingerprints)
        all_valid_indices.extend(valid_indices)
        
        print(f"    ✓ Batch {batch_idx + 1}: Generated {len(fingerprints):,} fingerprints")
        
        # Clear batch data to free memory
        del fingerprints, valid_indices
    
    X = np.array(all_fingerprints)
    print(f"  ✓ Total: Generated {len(all_fingerprints):,} fingerprints")
    return X, all_valid_indices


def train_property_models_batch(X_train, y_train, property_names, batch_size=BATCH_SIZE):
    """Train model using batch data to manage memory"""
    print(f"  Training multi-output property regressor...")
    print(f"  Properties: {property_names}")
    print(f"  Training samples: {len(X_train):,}")
    print(f"  Using batch size: {batch_size:,} for memory efficiency")
    
    base_model = RandomForestRegressor(
        n_estimators=30,
        random_state=42,
        n_jobs=-1,
        max_depth=10,
        min_samples_split=5
    )
    
    model = MultiOutputRegressor(base_model)
    
    # Train on full dataset (scikit-learn handles it efficiently)
    model.fit(X_train, y_train)
    
    print(f"  ✓ Training completed")
    return model


def evaluate_property_models(model, X_test, y_test, property_names):
    print(f"  Evaluating on {len(X_test):,} test samples...")
    y_pred = model.predict(X_test)
    
    results = {}
    print("\n  Properties Performance:")
    print(f"  {'-' * 50}")
    for i, prop in enumerate(property_names):
        rmse = np.sqrt(mean_squared_error(y_test[:, i], y_pred[:, i]))
        r2 = r2_score(y_test[:, i], y_pred[:, i])
        results[prop] = {'rmse': rmse, 'r2': r2}
        print(f"  {prop:<20} RMSE: {rmse:>8.4f}  R²: {r2:>7.4f}")
    print(f"  {'-' * 50}")
    return results


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
    
    # Configuration - adjust these based on your needs
    data_limit = None  # Set to None to use all data, or specify a number
    batch_size = BATCH_SIZE  # Configure batch size based on available RAM
    
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
        df = load_property_data(db_path, limit=data_limit)
        save_checkpoint(checkpoint_path, 'data_loaded', {'df': df})
    
    print("\n[Step 2/5] Cleaning dataset...")
    if checkpoint and checkpoint.get('step') == 'data_cleaned':
        df_clean = checkpoint.get('df_clean')
        print("✓ Using cleaned data from checkpoint")
    else:
        df_clean = clean_property_data(df)
        print(f"\n✓ Checkpoint saved: data_cleaned")
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
        X, valid_indices = compute_fingerprints_batch(smiles_list, batch_size=batch_size)
        
        df_valid = df_clean.iloc[valid_indices].copy()
        y = df_valid[property_columns].values
        
        scaler = StandardScaler()
        y_scaled = scaler.fit_transform(y)
        
        print(f"\n✓ Checkpoint saved: fingerprints_computed")
        save_checkpoint(checkpoint_path, 'fingerprints_computed', {'X': X, 'y': y, 'y_scaled': y_scaled, 'scaler': scaler})
    
    print(f"Final dataset: {len(X):,} molecules, {len(property_columns)} properties")
    
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
        print(f"Training: {len(X_train):,}, Test: {len(X_test):,}")
        model = train_property_models_batch(X_train, y_train, property_columns, batch_size=batch_size)
        
        print(f"\n✓ Checkpoint saved: model_trained")
        save_checkpoint(checkpoint_path, 'model_trained', {'model': model, 'X_test': X_test, 'y_test': y_test, 'scaler': scaler})
    
    print("\n[Step 5/5] Evaluating and saving...")
    y_test_orig = scaler.inverse_transform(y_test)
    y_pred_scaled = model.predict(X_test)
    y_pred_orig = scaler.inverse_transform(y_pred_scaled)
    
    results = {}
    print("=" * 60)
    print("Property Prediction Results (Original Scale)")
    print("=" * 60)
    
    for i, prop in enumerate(property_columns):
        rmse = np.sqrt(mean_squared_error(y_test_orig[:, i], y_pred_orig[:, i]))
        r2 = r2_score(y_test_orig[:, i], y_pred_orig[:, i])
        results[prop] = {'rmse': rmse, 'r2': r2}
        print(f"{prop:<20} RMSE: {rmse:>8.4f}  R²: {r2:>7.4f}")
    print("=" * 60)
    
    save_property_model(model, scaler, property_columns, output_dir)
    save_checkpoint(checkpoint_path, 'complete', {'results': results})
    
    print("\n✓ Property model training completed!")
    return model, results


if __name__ == "__main__":
    model, results = main()
