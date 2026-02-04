import warnings
import os
import io
import sys
import sqlite3
import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, Descriptors
from sklearn.metrics import mean_squared_error, r2_score, mean_absolute_error, accuracy_score, classification_report, confusion_matrix
import joblib

# Suppress warnings
warnings.filterwarnings('ignore')
os.environ['PYTHONWARNINGS'] = 'ignore'
RDLogger.DisableLog('rdApp.*')


def load_druglikeness_test_data(db_path, limit=5000):
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
        cp.num_ro5_violations,
        cp.qed_weighted,
        md.max_phase
    FROM compound_structures cs
    INNER JOIN compound_properties cp ON cs.molregno = cp.molregno
    INNER JOIN molecule_dictionary md ON cs.molregno = md.molregno
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
    print(f"Fingerprint shape: {X.shape} (should be (N, 1024))")
    return X, valid_indices


def calculate_drug_likeness_label(row):
    """Create drug-likeness label based on Lipinski's rule of 5"""
    violations = 0
    if row['mw_freebase'] > 500:
        violations += 1
    if row['alogp'] > 5:
        violations += 1
    if row['hbd'] > 5:
        violations += 1
    if row['hba'] > 10:
        violations += 1
    
    return 1 if violations <= 1 else 0


def test_druglikeness_models():
    """Test the QED and drug-likeness models"""
    print("\n" + "=" * 70)
    print("DRUGLIKENESS MODEL TESTING")
    print("=" * 70)
    
    # Load models
    qed_model_path = "trained_models/qed_model.joblib"
    druglike_model_path = "trained_models/druglikeness_model.joblib"
    
    if not os.path.exists(qed_model_path):
        print(f"✗ QED model file not found at {qed_model_path}")
        return None
    
    if not os.path.exists(druglike_model_path):
        print(f"✗ Drug-likeness model file not found at {druglike_model_path}")
        return None
    
    try:
        qed_model = joblib.load(qed_model_path)
        druglike_model = joblib.load(druglike_model_path)
        print(f"✓ Loaded QED model from {qed_model_path}")
        print(f"✓ Loaded drug-likeness model from {druglike_model_path}")
    except Exception as e:
        print(f"✗ Error loading models: {e}")
        return None
    
    # Load test data
    db_path = "chembl_36/chembl_36_sqlite/chembl_36.db"
    df = load_druglikeness_test_data(db_path, limit=5000)
    
    if df is None or len(df) == 0:
        print("✗ No test data available")
        return None
    
    # Compute fingerprints
    print("\nComputing fingerprints...")
    X, valid_indices = compute_fingerprints(df['canonical_smiles'].tolist())
    print(f"✓ Generated {len(X):,} valid fingerprints")
    
    # Filter data
    df_valid = df.iloc[valid_indices].reset_index(drop=True)
    
    # Create labels
    y_qed = df_valid['qed_weighted'].values
    y_druglike = df_valid.apply(calculate_drug_likeness_label, axis=1).values
    
    # Make predictions
    print("\nMaking predictions...")
    y_qed_pred = qed_model.predict(X)
    y_druglike_pred = druglike_model.predict(X)
    
    # Calculate QED metrics
    print("\n" + "-" * 70)
    print("QED REGRESSOR METRICS")
    print("-" * 70)
    
    qed_rmse = np.sqrt(mean_squared_error(y_qed, y_qed_pred))
    qed_mae = mean_absolute_error(y_qed, y_qed_pred)
    qed_r2 = r2_score(y_qed, y_qed_pred)
    
    print(f"RMSE: {qed_rmse:.4f}")
    print(f"MAE:  {qed_mae:.4f}")
    print(f"R²:   {qed_r2:.4f}")
    print(f"Test Samples: {len(X):,}")
    
    # Calculate drug-likeness metrics
    print("\n" + "-" * 70)
    print("DRUG-LIKENESS CLASSIFIER METRICS")
    print("-" * 70)
    
    accuracy = accuracy_score(y_druglike, y_druglike_pred)
    cm = confusion_matrix(y_druglike, y_druglike_pred)
    
    print(f"Accuracy: {accuracy:.4f}")
    print(f"Test Samples: {len(X):,}")
    print("\nConfusion Matrix:")
    print(f"  True Negatives:  {cm[0, 0]}")
    print(f"  False Positives: {cm[0, 1]}")
    print(f"  False Negatives: {cm[1, 0]}")
    print(f"  True Positives:  {cm[1, 1]}")
    
    print("\nClassification Report:")
    print(classification_report(y_druglike, y_druglike_pred, 
                               target_names=['Not Drug-like', 'Drug-like']))
    
    print("-" * 70)
    
    return {
        'model_name': 'Druglikeness Models',
        'test_samples': len(X),
        'qed_metrics': {
            'rmse': qed_rmse,
            'mae': qed_mae,
            'r2': qed_r2,
            'y_true': y_qed,
            'y_pred': y_qed_pred
        },
        'druglike_metrics': {
            'accuracy': accuracy,
            'confusion_matrix': cm,
            'y_true': y_druglike,
            'y_pred': y_druglike_pred
        }
    }


if __name__ == "__main__":
    results = test_druglikeness_models()
    if results:
        print("\n✓ Druglikeness models testing completed successfully!")
