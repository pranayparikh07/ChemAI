import warnings
import os
import io
import sys
import sqlite3
import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
from sklearn.metrics import accuracy_score, roc_auc_score, classification_report, confusion_matrix
import joblib

# Suppress warnings
warnings.filterwarnings('ignore')
os.environ['PYTHONWARNINGS'] = 'ignore'
RDLogger.DisableLog('rdApp.*')


def load_toxicity_test_data(db_path, limit=5000):
    """Load test data from ChEMBL database"""
    query = f"""
    SELECT 
        cs.molregno,
        cs.canonical_smiles,
        COUNT(DISTINCT csa.alert_id) as alert_count
    FROM compound_structures cs
    LEFT JOIN compound_structural_alerts csa ON cs.molregno = csa.molregno
    WHERE cs.canonical_smiles IS NOT NULL
    GROUP BY cs.molregno, cs.canonical_smiles
    LIMIT {limit}
    """
    try:
        conn = sqlite3.connect(db_path)
        df = pd.read_sql_query(query, conn)
        conn.close()
        df['has_alerts'] = (df['alert_count'] > 0).astype(int)
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


def test_toxicity_model():
    """Test the toxicity prediction model"""
    print("\n" + "=" * 70)
    print("TOXICITY MODEL TESTING")
    print("=" * 70)
    
    # Load models
    model_path = "trained_models/toxicity_model.joblib"
    
    if not os.path.exists(model_path):
        print(f"✗ Toxicity model file not found at {model_path}")
        return None
    
    try:
        model = joblib.load(model_path)
        print(f"✓ Loaded toxicity model from {model_path}")
    except Exception as e:
        print(f"✗ Error loading model: {e}")
        return None
    
    # Load test data
    db_path = "chembl_36/chembl_36_sqlite/chembl_36.db"
    df = load_toxicity_test_data(db_path, limit=5000)
    
    if df is None or len(df) == 0:
        print("✗ No test data available")
        return None
    
    # Compute fingerprints
    print("\nComputing fingerprints...")
    X, valid_indices = compute_fingerprints(df['canonical_smiles'].tolist())
    print(f"✓ Generated {len(X):,} valid fingerprints")
    
    # Filter data
    df_valid = df.iloc[valid_indices].reset_index(drop=True)
    y_true = df_valid['has_alerts'].values
    
    # Make predictions
    print("\nMaking predictions...")
    y_pred = model.predict(X)
    y_prob = model.predict_proba(X)[:, 1]
    
    # Calculate metrics
    print("\n" + "-" * 70)
    print("TOXICITY MODEL METRICS")
    print("-" * 70)
    
    accuracy = accuracy_score(y_true, y_pred)
    roc_auc = roc_auc_score(y_true, y_prob)
    conf_matrix = confusion_matrix(y_true, y_pred)
    
    print(f"\nAccuracy: {accuracy:.4f}")
    print(f"ROC-AUC: {roc_auc:.4f}")
    print(f"\nConfusion Matrix:")
    print(f"  True Negatives:  {conf_matrix[0, 0]}")
    print(f"  False Positives: {conf_matrix[0, 1]}")
    print(f"  False Negatives: {conf_matrix[1, 0]}")
    print(f"  True Positives:  {conf_matrix[1, 1]}")
    
    print("\n" + "-" * 70)
    print("Classification Report:")
    print("-" * 70)
    print(classification_report(y_true, y_pred, target_names=['Safe', 'Toxic']))
    print("-" * 70)
    
    return {
        'test_samples': len(X),
        'metrics': {
            'accuracy': accuracy,
            'roc_auc': roc_auc,
            'confusion_matrix': conf_matrix
        },
        'y_true': y_true,
        'y_pred': y_pred,
        'y_prob': y_prob
    }


if __name__ == "__main__":
    results = test_toxicity_model()
    if results:
        print(f"\n✓ Toxicity model test completed successfully!")
    else:
        print("\n✗ Toxicity model test failed!")
