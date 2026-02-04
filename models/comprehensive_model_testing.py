"""
Comprehensive Model Testing - Calculate Accuracy, F1, and R² for all models
"""

import warnings
import os
import io
import sys
import sqlite3
import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
from sklearn.metrics import (
    accuracy_score, f1_score, r2_score, 
    mean_squared_error, mean_absolute_error,
    classification_report, confusion_matrix
)
import joblib

# Suppress warnings
warnings.filterwarnings('ignore')
os.environ['PYTHONWARNINGS'] = 'ignore'
RDLogger.DisableLog('rdApp.*')


class ModelTester:
    """Comprehensive model testing framework"""
    
    def __init__(self, db_path="chembl_36/chembl_36_sqlite/chembl_36.db"):
        self.db_path = db_path
        self.results = {}
    
    def compute_fingerprints(self, smiles_list, n_bits=1024):
        """Generate Morgan fingerprints for SMILES"""
        fingerprints = []
        valid_indices = []
        
        for idx, smiles in enumerate(smiles_list):
            mol = Chem.MolFromSmiles(smiles)
            if mol is not None:
                old_stderr = sys.stderr
                sys.stderr = io.StringIO()
                try:
                    fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=n_bits)
                    fingerprints.append(np.array(fp))
                    valid_indices.append(idx)
                finally:
                    sys.stderr = old_stderr
        
        X = np.array(fingerprints) if fingerprints else None
        return X, valid_indices
    
    # ==================== BIOACTIVITY MODEL ====================
    def test_bioactivity_model(self):
        """Test bioactivity model - Regression"""
        print("\n" + "=" * 80)
        print("BIOACTIVITY MODEL - REGRESSION METRICS")
        print("=" * 80)
        
        model_path = "trained_models/bioactivity_model.joblib"
        
        if not os.path.exists(model_path):
            print(f"✗ Model not found at {model_path}")
            return None
        
        try:
            model = joblib.load(model_path)
            print(f"[+] Loaded bioactivity model")
        except Exception as e:
            print(f"✗ Error loading model: {e}")
            return None
        
        # Load test data
        query = """
        SELECT DISTINCT
            cs.canonical_smiles,
            act.standard_value
        FROM molecule_dictionary md
        INNER JOIN compound_structures cs ON md.molregno = cs.molregno
        INNER JOIN activities act ON md.molregno = act.molregno
        INNER JOIN assays a ON act.assay_id = a.assay_id
        WHERE act.standard_type = 'IC50'
            AND act.standard_units = 'nM'
            AND act.standard_value IS NOT NULL
            AND act.standard_value > 0
            AND a.assay_type = 'B'
            AND cs.canonical_smiles IS NOT NULL
        LIMIT 5000
        """
        
        try:
            conn = sqlite3.connect(self.db_path)
            df = pd.read_sql_query(query, conn)
            conn.close()
        except Exception as e:
            print(f"✗ Error loading test data: {e}")
            return None
        
        # Calculate pIC50
        df['pIC50'] = -np.log10(df['standard_value'] * 1e-9)
        df = df[(df['pIC50'] >= 0) & (df['pIC50'] <= 15)].copy()
        
        # Compute fingerprints
        X, valid_indices = self.compute_fingerprints(df['canonical_smiles'].tolist(), n_bits=2048)
        df_valid = df.iloc[valid_indices].reset_index(drop=True)
        y_true = df_valid['pIC50'].values
        
        # Make predictions
        y_pred = model.predict(X)
        
        # Calculate metrics
        r2 = r2_score(y_true, y_pred)
        rmse = np.sqrt(mean_squared_error(y_true, y_pred))
        mae = mean_absolute_error(y_true, y_pred)
        
        # Note: Accuracy and F1 not applicable for regression
        results = {
            'model_name': 'Bioactivity Model',
            'type': 'Regression (pIC50)',
            'test_samples': len(X),
            'metrics': {
                'R²': r2,
                'RMSE': rmse,
                'MAE': mae,
                'Accuracy': 'N/A (Regression)',
                'F1': 'N/A (Regression)'
            }
        }
        
        print(f"\nModel: Bioactivity Prediction (pIC50)")
        print(f"Test Samples: {len(X):,}")
        print(f"\nMetrics:")
        print(f"  R²:       {r2:.4f}")
        print(f"  RMSE:     {rmse:.4f}")
        print(f"  MAE:      {mae:.4f}")
        print(f"  Accuracy: N/A (Regression Task)")
        print(f"  F1-Score: N/A (Regression Task)")
        
        self.results['Bioactivity'] = results
        return results
    
    # ==================== TOXICITY MODEL ====================
    def test_toxicity_model(self):
        """Test toxicity model - Classification"""
        print("\n" + "=" * 80)
        print("TOXICITY MODEL - CLASSIFICATION METRICS")
        print("=" * 80)
        
        model_path = "trained_models/toxicity_model.joblib"
        
        if not os.path.exists(model_path):
            print(f"✗ Model not found at {model_path}")
            return None
        
        try:
            model = joblib.load(model_path)
            print(f"✓ Loaded toxicity model")
        except Exception as e:
            print(f"✗ Error loading model: {e}")
            return None
        
        # Load test data
        query = """
        SELECT 
            cs.canonical_smiles,
            COUNT(DISTINCT csa.alert_id) as alert_count
        FROM compound_structures cs
        LEFT JOIN compound_structural_alerts csa ON cs.molregno = csa.molregno
        WHERE cs.canonical_smiles IS NOT NULL
        GROUP BY cs.molregno, cs.canonical_smiles
        LIMIT 5000
        """
        
        try:
            conn = sqlite3.connect(self.db_path)
            df = pd.read_sql_query(query, conn)
            conn.close()
        except Exception as e:
            print(f"✗ Error loading test data: {e}")
            return None
        
        # Create labels (has_alerts = 1 if alert_count > 0)
        y_true = (df['alert_count'] > 0).astype(int).values
        
        # Compute fingerprints
        X, valid_indices = self.compute_fingerprints(df['canonical_smiles'].tolist(), n_bits=1024)
        y_true = y_true[valid_indices]
        
        # Make predictions
        y_pred = model.predict(X)
        y_pred_proba = model.predict_proba(X)[:, 1]
        
        # Calculate metrics
        accuracy = accuracy_score(y_true, y_pred)
        f1 = f1_score(y_true, y_pred, zero_division=0)
        # R² not directly applicable for classification, but can use on probabilities
        r2_proba = r2_score(y_true, y_pred_proba)
        
        results = {
            'model_name': 'Toxicity Model',
            'type': 'Classification (Binary)',
            'test_samples': len(X),
            'metrics': {
                'Accuracy': accuracy,
                'F1-Score': f1,
                'R²': r2_proba,
                'Confusion Matrix': confusion_matrix(y_true, y_pred).tolist()
            }
        }
        
        print(f"\nModel: Toxicity Prediction (Structural Alerts)")
        print(f"Test Samples: {len(X):,}")
        print(f"\nMetrics:")
        print(f"  Accuracy:     {accuracy:.4f}")
        print(f"  F1-Score:     {f1:.4f}")
        print(f"  R² (Proba):   {r2_proba:.4f}")
        print(f"\nConfusion Matrix:")
        cm = confusion_matrix(y_true, y_pred)
        print(f"  TN: {cm[0,0]}, FP: {cm[0,1]}")
        print(f"  FN: {cm[1,0]}, TP: {cm[1,1]}")
        
        self.results['Toxicity'] = results
        return results
    
    # ==================== PROPERTY MODEL ====================
    def test_property_model(self):
        """Test property model - Multi-target Regression"""
        print("\n" + "=" * 80)
        print("PROPERTY MODEL - MULTI-TARGET REGRESSION METRICS")
        print("=" * 80)
        
        model_path = "trained_models/property_model.joblib"
        scaler_path = "trained_models/property_scaler.joblib"
        
        if not os.path.exists(model_path) or not os.path.exists(scaler_path):
            print(f"✗ Model files not found")
            return None
        
        try:
            model = joblib.load(model_path)
            scaler = joblib.load(scaler_path)
            print(f"✓ Loaded property model and scaler")
        except Exception as e:
            print(f"✗ Error loading models: {e}")
            return None
        
        # Load test data
        query = """
        SELECT 
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
        LIMIT 5000
        """
        
        try:
            conn = sqlite3.connect(self.db_path)
            df = pd.read_sql_query(query, conn)
            conn.close()
        except Exception as e:
            print(f"✗ Error loading test data: {e}")
            return None
        
        # Compute fingerprints
        X, valid_indices = self.compute_fingerprints(df['canonical_smiles'].tolist(), n_bits=1024)
        df_valid = df.iloc[valid_indices].reset_index(drop=True)
        
        property_columns = ['mw_freebase', 'alogp', 'hba', 'hbd', 'psa', 'rtb', 'qed_weighted']
        y_true = df_valid[property_columns].values
        
        # Make predictions
        y_pred_scaled = model.predict(X)
        y_pred = scaler.inverse_transform(y_pred_scaled)
        
        # Calculate metrics for each property
        results_dict = {
            'model_name': 'Property Model',
            'type': 'Multi-target Regression (7 properties)',
            'test_samples': len(X),
            'metrics': {}
        }
        
        print(f"\nModel: Molecular Property Prediction")
        print(f"Test Samples: {len(X):,}")
        print(f"Properties: {len(property_columns)}")
        print("\n" + "-" * 80)
        
        for i, prop in enumerate(property_columns):
            r2 = r2_score(y_true[:, i], y_pred[:, i])
            rmse = np.sqrt(mean_squared_error(y_true[:, i], y_pred[:, i]))
            mae = mean_absolute_error(y_true[:, i], y_pred[:, i])
            
            results_dict['metrics'][prop] = {
                'R²': r2,
                'RMSE': rmse,
                'MAE': mae,
                'Accuracy': 'N/A',
                'F1': 'N/A'
            }
            
            print(f"\n{prop.upper()}")
            print(f"  R²:       {r2:.4f}")
            print(f"  RMSE:     {rmse:.4f}")
            print(f"  MAE:      {mae:.4f}")
        
        self.results['Property'] = results_dict
        return results_dict
    
    # ==================== DRUGLIKENESS MODEL ====================
    def test_druglikeness_model(self):
        """Test druglikeness model - Classification"""
        print("\n" + "=" * 80)
        print("DRUGLIKENESS MODEL - CLASSIFICATION METRICS")
        print("=" * 80)
        
        qed_model_path = "trained_models/qed_model.joblib"
        druglike_model_path = "trained_models/druglikeness_model.joblib"
        
        if not os.path.exists(qed_model_path) or not os.path.exists(druglike_model_path):
            print(f"✗ Model files not found")
            return None
        
        try:
            qed_model = joblib.load(qed_model_path)
            druglike_model = joblib.load(druglike_model_path)
            print(f"✓ Loaded QED and druglikeness models")
        except Exception as e:
            print(f"✗ Error loading models: {e}")
            return None
        
        # Load test data
        query = """
        SELECT 
            cs.canonical_smiles,
            cp.mw_freebase,
            cp.alogp,
            cp.hba,
            cp.hbd,
            cp.psa,
            cp.rtb,
            cp.qed_weighted,
            cp.num_ro5_violations
        FROM compound_structures cs
        INNER JOIN compound_properties cp ON cs.molregno = cp.molregno
        WHERE cs.canonical_smiles IS NOT NULL
            AND cp.mw_freebase IS NOT NULL
            AND cp.qed_weighted IS NOT NULL
        LIMIT 5000
        """
        
        try:
            conn = sqlite3.connect(self.db_path)
            df = pd.read_sql_query(query, conn)
            conn.close()
        except Exception as e:
            print(f"✗ Error loading test data: {e}")
            return None
        
        # Compute fingerprints
        X, valid_indices = self.compute_fingerprints(df['canonical_smiles'].tolist(), n_bits=1024)
        df_valid = df.iloc[valid_indices].reset_index(drop=True)
        
        # QED Regression
        y_qed_true = df_valid['qed_weighted'].values
        y_qed_pred = qed_model.predict(X)
        
        # Drug-likeness Classification
        y_druglike_true = (df_valid['num_ro5_violations'] <= 1).astype(int).values
        y_druglike_pred = druglike_model.predict(X)
        y_druglike_proba = druglike_model.predict_proba(X)[:, 1]
        
        # Calculate QED metrics
        qed_r2 = r2_score(y_qed_true, y_qed_pred)
        qed_rmse = np.sqrt(mean_squared_error(y_qed_true, y_qed_pred))
        qed_mae = mean_absolute_error(y_qed_true, y_qed_pred)
        
        # Calculate Drug-likeness metrics
        druglike_accuracy = accuracy_score(y_druglike_true, y_druglike_pred)
        druglike_f1 = f1_score(y_druglike_true, y_druglike_pred, zero_division=0)
        druglike_r2 = r2_score(y_druglike_true, y_druglike_proba)
        
        results = {
            'model_name': 'Druglikeness Model',
            'type': 'Mixed (Regression + Classification)',
            'test_samples': len(X),
            'metrics': {
                'QED (Regression)': {
                    'R²': qed_r2,
                    'RMSE': qed_rmse,
                    'MAE': qed_mae,
                    'Accuracy': 'N/A',
                    'F1': 'N/A'
                },
                'Drug-likeness (Classification)': {
                    'Accuracy': druglike_accuracy,
                    'F1-Score': druglike_f1,
                    'R²': druglike_r2,
                    'Confusion Matrix': confusion_matrix(y_druglike_true, y_druglike_pred).tolist()
                }
            }
        }
        
        print(f"\nModel 1: QED Prediction (Regression)")
        print(f"  R²:       {qed_r2:.4f}")
        print(f"  RMSE:     {qed_rmse:.4f}")
        print(f"  MAE:      {qed_mae:.4f}")
        
        print(f"\nModel 2: Drug-likeness Classification (Lipinski's Rule)")
        print(f"  Accuracy:     {druglike_accuracy:.4f}")
        print(f"  F1-Score:     {druglike_f1:.4f}")
        print(f"  R² (Proba):   {druglike_r2:.4f}")
        print(f"\nConfusion Matrix:")
        cm = confusion_matrix(y_druglike_true, y_druglike_pred)
        print(f"  TN: {cm[0,0]}, FP: {cm[0,1]}")
        print(f"  FN: {cm[1,0]}, TP: {cm[1,1]}")
        
        self.results['Druglikeness'] = results
        return results
    
    def generate_summary_report(self):
        """Generate comprehensive summary report"""
        print("\n" + "=" * 80)
        print("COMPREHENSIVE MODEL TESTING SUMMARY")
        print("=" * 80)
        
        summary_data = []
        
        for model_name, result in self.results.items():
            if result is None:
                continue
            
            print(f"\n{model_name.upper()} MODEL")
            print(f"  Type: {result['type']}")
            print(f"  Test Samples: {result['test_samples']:,}")
            
            if isinstance(result['metrics'], dict):
                if model_name == 'Property':
                    print(f"  Properties: {len(result['metrics'])}")
                    for prop, metrics in result['metrics'].items():
                        print(f"    {prop:15} | R²: {metrics.get('R²', 'N/A'):.4f} | " + 
                              f"F1: {metrics.get('F1', 'N/A')} | " +
                              f"Accuracy: {metrics.get('Accuracy', 'N/A')}")
                        summary_data.append({
                            'Model': model_name,
                            'Target': prop,
                            'Accuracy': metrics.get('Accuracy', 'N/A'),
                            'F1-Score': metrics.get('F1', 'N/A'),
                            'R²': f"{metrics.get('R²', 0):.4f}"
                        })
                elif model_name == 'Druglikeness':
                    # QED
                    print(f"  QED (Regression):")
                    print(f"    R²: {result['metrics']['QED (Regression)'].get('R²', 'N/A'):.4f}")
                    print(f"    F1: {result['metrics']['QED (Regression)'].get('F1', 'N/A')}")
                    print(f"    Accuracy: {result['metrics']['QED (Regression)'].get('Accuracy', 'N/A')}")
                    
                    # Drug-likeness
                    print(f"  Drug-likeness (Classification):")
                    print(f"    Accuracy: {result['metrics']['Drug-likeness (Classification)'].get('Accuracy', 'N/A'):.4f}")
                    print(f"    F1-Score: {result['metrics']['Drug-likeness (Classification)'].get('F1-Score', 'N/A'):.4f}")
                    print(f"    R²: {result['metrics']['Drug-likeness (Classification)'].get('R²', 'N/A'):.4f}")
                    
                    summary_data.append({
                        'Model': model_name,
                        'Target': 'QED',
                        'Accuracy': 'N/A',
                        'F1-Score': 'N/A',
                        'R²': f"{result['metrics']['QED (Regression)'].get('R²', 0):.4f}"
                    })
                    summary_data.append({
                        'Model': model_name,
                        'Target': 'Drug-likeness',
                        'Accuracy': f"{result['metrics']['Drug-likeness (Classification)'].get('Accuracy', 0):.4f}",
                        'F1-Score': f"{result['metrics']['Drug-likeness (Classification)'].get('F1-Score', 0):.4f}",
                        'R²': f"{result['metrics']['Drug-likeness (Classification)'].get('R²', 0):.4f}"
                    })
                else:
                    acc = result['metrics'].get('Accuracy', 'N/A')
                    f1 = result['metrics'].get('F1-Score', 'N/A')
                    r2 = result['metrics'].get('R²', 'N/A')
                    
                    if isinstance(acc, (int, float)):
                        print(f"  Accuracy: {acc:.4f}")
                    else:
                        print(f"  Accuracy: {acc}")
                    
                    if isinstance(f1, (int, float)):
                        print(f"  F1-Score: {f1:.4f}")
                    else:
                        print(f"  F1-Score: {f1}")
                        
                    if isinstance(r2, (int, float)):
                        print(f"  R²: {r2:.4f}")
                    else:
                        print(f"  R²: {r2}")
                    
                    acc_str = f"{acc:.4f}" if isinstance(acc, (int, float)) else str(acc)
                    f1_str = f"{f1:.4f}" if isinstance(f1, (int, float)) else str(f1)
                    r2_str = f"{r2:.4f}" if isinstance(r2, (int, float)) else str(r2)
                    
                    summary_data.append({
                        'Model': model_name,
                        'Target': model_name,
                        'Accuracy': acc_str,
                        'F1-Score': f1_str,
                        'R²': r2_str
                    })
        
        # Create summary DataFrame
        if summary_data:
            summary_df = pd.DataFrame(summary_data)
            print("\n" + "=" * 80)
            print("FINAL SUMMARY TABLE")
            print("=" * 80)
            print(summary_df.to_string(index=False))
        
        print("\n" + "=" * 80)


def main():
    """Main testing function"""
    tester = ModelTester()
    
    # Test all models
    tester.test_bioactivity_model()
    tester.test_toxicity_model()
    tester.test_property_model()
    tester.test_druglikeness_model()
    
    # Generate summary
    tester.generate_summary_report()


if __name__ == "__main__":
    main()
