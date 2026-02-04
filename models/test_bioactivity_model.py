import warnings
import os
import io
import sys
import sqlite3
import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
from sklearn.metrics import mean_squared_error, r2_score, mean_absolute_error
import joblib

# Suppress warnings
warnings.filterwarnings('ignore')
os.environ['PYTHONWARNINGS'] = 'ignore'
RDLogger.DisableLog('rdApp.*')


def load_bioactivity_test_data(db_path, limit=5000):
    """Load test data from ChEMBL database"""
    query = f"""
    SELECT DISTINCT
        md.molregno,
        cs.canonical_smiles,
        act.standard_value,
        td.pref_name as target_name
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
    try:
        conn = sqlite3.connect(db_path)
        df = pd.read_sql_query(query, conn)
        conn.close()
        print(f"✓ Loaded {len(df):,} test samples")
        return df
    except Exception as e:
        print(f"✗ Error loading test data: {e}")
        return None


def compute_fingerprints(smiles_list):
    """Generate Morgan fingerprints for SMILES"""
    fingerprints = []
    valid_indices = []
    
    for idx, smiles in enumerate(smiles_list):
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            old_stderr = sys.stderr
            sys.stderr = io.StringIO()
            try:
                fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
                fingerprints.append(np.array(fp))
                valid_indices.append(idx)
            finally:
                sys.stderr = old_stderr
    
    X = np.array(fingerprints)
    return X, valid_indices


def test_bioactivity_model():
    """Test the bioactivity prediction model"""
    print("\n" + "=" * 70)
    print("BIOACTIVITY MODEL TESTING")
    print("=" * 70)
    
    # Load model
    model_path = "trained_models/bioactivity_model.joblib"
    
    if not os.path.exists(model_path):
        print(f"✗ Model file not found at {model_path}")
        print("Note: Bioactivity model may not be trained yet.")
        return None
    
    try:
        model = joblib.load(model_path)
        print(f"✓ Loaded bioactivity model from {model_path}")
    except Exception as e:
        print(f"✗ Error loading model: {e}")
        return None
    
    # Load test data
    db_path = "chembl_36/chembl_36_sqlite/chembl_36.db"
    df = load_bioactivity_test_data(db_path, limit=5000)
    
    if df is None or len(df) == 0:
        print("✗ No test data available")
        return None
    
    # Calculate pIC50
    df['pIC50'] = -np.log10(df['standard_value'] * 1e-9)
    df = df[(df['pIC50'] >= 0) & (df['pIC50'] <= 15)].copy()
    
    # Compute fingerprints
    print("\nComputing fingerprints...")
    X, valid_indices = compute_fingerprints(df['canonical_smiles'].tolist())
    print(f"✓ Generated {len(X):,} valid fingerprints")
    
    # Filter data
    df_valid = df.iloc[valid_indices].reset_index(drop=True)
    y_test = df_valid['pIC50'].values
    
    # Make predictions
    print("\nMaking predictions...")
    y_pred = model.predict(X)
    
    # Calculate metrics
    print("\n" + "-" * 70)
    print("BIOACTIVITY MODEL METRICS")
    print("-" * 70)
    
    rmse = np.sqrt(mean_squared_error(y_test, y_pred))
    mae = mean_absolute_error(y_test, y_pred)
    r2 = r2_score(y_test, y_pred)
    
    print(f"\nTarget: pIC50 (Bioactivity)")
    print(f"  RMSE: {rmse:.4f}")
    print(f"  MAE:  {mae:.4f}")
    print(f"  R²:   {r2:.4f}")
    print(f"  Test Samples: {len(X):,}")
    
    # Additional statistics
    residuals = y_test - y_pred
    print(f"\nResiduals Statistics:")
    print(f"  Mean: {np.mean(residuals):.4f}")
    print(f"  Std:  {np.std(residuals):.4f}")
    print(f"  Min:  {np.min(residuals):.4f}")
    print(f"  Max:  {np.max(residuals):.4f}")
    
    print("\n" + "-" * 70)
    
    return {
        'model_name': 'Bioactivity Model',
        'test_samples': len(X),
        'metrics': {
            'rmse': rmse,
            'mae': mae,
            'r2': r2,
            'residuals_mean': np.mean(residuals),
            'residuals_std': np.std(residuals),
            'residuals_min': np.min(residuals),
            'residuals_max': np.max(residuals),
            'y_true': y_test,
            'y_pred': y_pred
        }
    }


if __name__ == "__main__":
    results = test_bioactivity_model()
    if results:
        print("\n✓ Bioactivity model testing completed successfully!")
