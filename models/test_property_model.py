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


def load_test_data(db_path, limit=5000):
    """Load test data from ChEMBL database"""
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
        cp.qed_weighted
    FROM compound_structures cs
    INNER JOIN compound_properties cp ON cs.molregno = cp.molregno
    WHERE cs.canonical_smiles IS NOT NULL
        AND cp.mw_freebase IS NOT NULL
        AND cp.qed_weighted IS NOT NULL
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
                fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024)
                fingerprints.append(np.array(fp))
                valid_indices.append(idx)
            finally:
                sys.stderr = old_stderr
    
    X = np.array(fingerprints)
    return X, valid_indices


def test_property_model():
    """Test the property prediction model"""
    print("\n" + "=" * 70)
    print("PROPERTY MODEL TESTING")
    print("=" * 70)
    
    # Load model and scaler
    model_path = "trained_models/property_model.joblib"
    scaler_path = "trained_models/property_scaler.joblib"
    
    if not os.path.exists(model_path) or not os.path.exists(scaler_path):
        print("✗ Model files not found!")
        return None
    
    try:
        model = joblib.load(model_path)
        scaler = joblib.load(scaler_path)
        print(f"✓ Loaded property model from {model_path}")
        print(f"✓ Loaded scaler from {scaler_path}")
    except Exception as e:
        print(f"✗ Error loading model: {e}")
        return None
    
    # Load test data
    db_path = "chembl_36/chembl_36_sqlite/chembl_36.db"
    df = load_test_data(db_path, limit=5000)
    
    if df is None or len(df) == 0:
        print("✗ No test data available")
        return None
    
    # Compute fingerprints
    print("\nComputing fingerprints...")
    X, valid_indices = compute_fingerprints(df['canonical_smiles'].tolist())
    print(f"✓ Generated {len(X):,} valid fingerprints")
    
    # Filter data
    df_valid = df.iloc[valid_indices].reset_index(drop=True)
    property_columns = ['mw_freebase', 'alogp', 'hba', 'hbd', 'psa', 'rtb', 'qed_weighted']
    y_test = df_valid[property_columns].values
    
    # Make predictions
    print("\nMaking predictions...")
    y_pred_scaled = model.predict(X)
    y_pred = scaler.inverse_transform(y_pred_scaled)
    
    # Calculate metrics
    print("\n" + "-" * 70)
    print("PROPERTY MODEL METRICS")
    print("-" * 70)
    
    metrics = {}
    for i, prop in enumerate(property_columns):
        rmse = np.sqrt(mean_squared_error(y_test[:, i], y_pred[:, i]))
        mae = mean_absolute_error(y_test[:, i], y_pred[:, i])
        r2 = r2_score(y_test[:, i], y_pred[:, i])
        
        metrics[prop] = {
            'rmse': rmse,
            'mae': mae,
            'r2': r2,
            'y_true': y_test[:, i],
            'y_pred': y_pred[:, i]
        }
        
        print(f"\n{prop:>15}")
        print(f"  RMSE: {rmse:>10.4f}")
        print(f"  MAE:  {mae:>10.4f}")
        print(f"  R²:   {r2:>10.4f}")
    
    print("\n" + "-" * 70)
    
    return {
        'model_name': 'Property Model',
        'test_samples': len(X),
        'properties': property_columns,
        'metrics': metrics,
        'y_true': y_test,
        'y_pred': y_pred
    }


if __name__ == "__main__":
    results = test_property_model()
    if results:
        print("\n✓ Property model testing completed successfully!")
