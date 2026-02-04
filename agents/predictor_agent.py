import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
import joblib
import os
import warnings

# Suppress RDKit and scikit-learn warnings
warnings.filterwarnings('ignore')
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

# Suppress sklearn parallel warning
import sklearn.utils.parallel
original_warn = warnings.warn
def custom_warn(*args, **kwargs):
    if 'sklearn.utils.parallel.delayed' not in str(args):
        original_warn(*args, **kwargs)
warnings.warn = custom_warn


class PredictorAgent:
    def __init__(self, models_dir='trained_models'):
        self.models_dir = models_dir
        self.models = {}
        self.configs = {}
        self._load_models()
    
    def _load_models(self):
        print("PredictorAgent: Loading trained models...")
        
        model_files = {
            'bioactivity': 'bioactivity_model.joblib',
            'property': 'property_model.joblib',
            'toxicity': 'toxicity_model.joblib',
            'qed': 'qed_model.joblib',
            'druglikeness': 'druglikeness_model.joblib'
        }
        
        config_files = {
            'bioactivity': 'bioactivity_config.joblib',
            'property': 'property_config.joblib',
            'toxicity': 'toxicity_config.joblib',
            'druglikeness': 'druglikeness_config.joblib'
        }
        
        for name, filename in model_files.items():
            path = os.path.join(self.models_dir, filename)
            if os.path.exists(path):
                self.models[name] = joblib.load(path)
                print(f"  ✓ Loaded {name} model")
            else:
                print(f"  ✗ {name} model not found at {path}")
        
        for name, filename in config_files.items():
            path = os.path.join(self.models_dir, filename)
            if os.path.exists(path):
                self.configs[name] = joblib.load(path)
        
        scaler_path = os.path.join(self.models_dir, 'property_scaler.joblib')
        if os.path.exists(scaler_path):
            self.property_scaler = joblib.load(scaler_path)
        else:
            self.property_scaler = None
        
        alerts_path = os.path.join(self.models_dir, 'structural_alerts.joblib')
        if os.path.exists(alerts_path):
            self.structural_alerts = joblib.load(alerts_path)
        else:
            self.structural_alerts = []
        
        print(f"PredictorAgent: {len(self.models)} models loaded")
    
    def _smiles_to_fingerprint(self, smiles, radius=2, n_bits=2048):
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=n_bits)
        return np.array(fp).reshape(1, -1)
    
    def _smiles_to_fingerprint_1024(self, smiles):
        return self._smiles_to_fingerprint(smiles, radius=2, n_bits=1024)
    
    def predict_bioactivity(self, smiles):
        if 'bioactivity' not in self.models:
            return {'error': 'Bioactivity model not loaded'}
        
        fp = self._smiles_to_fingerprint(smiles, n_bits=2048)
        if fp is None:
            return {'error': 'Invalid SMILES'}
        
        pIC50 = self.models['bioactivity'].predict(fp)[0]
        IC50_nM = 10 ** (9 - pIC50)
        
        if pIC50 >= 8:
            potency = 'Very High'
        elif pIC50 >= 7:
            potency = 'High'
        elif pIC50 >= 6:
            potency = 'Moderate'
        elif pIC50 >= 5:
            potency = 'Low'
        else:
            potency = 'Very Low'
        
        return {
            'pIC50': round(pIC50, 3),
            'IC50_nM': round(IC50_nM, 2),
            'potency_class': potency
        }
    
    def predict_properties(self, smiles):
        if 'property' not in self.models:
            return {'error': 'Property model not loaded'}
        
        fp = self._smiles_to_fingerprint_1024(smiles)
        if fp is None:
            return {'error': 'Invalid SMILES'}
        
        pred_scaled = self.models['property'].predict(fp)
        
        if self.property_scaler is not None:
            pred = self.property_scaler.inverse_transform(pred_scaled)[0]
        else:
            pred = pred_scaled[0]
        
        property_names = self.configs.get('property', {}).get(
            'property_names', 
            ['mw_freebase', 'alogp', 'hba', 'hbd', 'psa', 'rtb', 'qed_weighted']
        )
        
        return {name: round(float(pred[i]), 3) for i, name in enumerate(property_names)}
    
    def predict_toxicity(self, smiles):
        if 'toxicity' not in self.models:
            return {'error': 'Toxicity model not loaded'}
        
        fp = self._smiles_to_fingerprint_1024(smiles)
        if fp is None:
            return {'error': 'Invalid SMILES'}
        
        has_alerts = self.models['toxicity'].predict(fp)[0]
        alert_prob = self.models['toxicity'].predict_proba(fp)[0][1]
        
        matched_alerts = []
        mol = Chem.MolFromSmiles(smiles)
        if mol and self.structural_alerts:
            for alert in self.structural_alerts:
                try:
                    pattern = Chem.MolFromSmarts(alert['smarts'])
                    if pattern and mol.HasSubstructMatch(pattern):
                        matched_alerts.append(alert['alert_name'])
                except:
                    # Skip invalid SMARTS patterns
                    pass
        
        if alert_prob < 0.3:
            risk = 'Low'
        elif alert_prob < 0.6:
            risk = 'Medium'
        else:
            risk = 'High'
        
        return {
            'has_structural_alerts': bool(has_alerts),
            'toxicity_probability': round(alert_prob, 3),
            'risk_level': risk,
            'matched_alerts': matched_alerts[:5]
        }
    
    def predict_druglikeness(self, smiles):
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {'error': 'Invalid SMILES'}
        
        from rdkit.Chem import Descriptors, QED
        
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        hbd = Descriptors.NumHDonors(mol)
        hba = Descriptors.NumHAcceptors(mol)
        psa = Descriptors.TPSA(mol)
        rotatable_bonds = Descriptors.NumRotatableBonds(mol)
        
        qed_score = QED.qed(mol)
        
        violations = 0
        violation_details = []
        if mw > 500:
            violations += 1
            violation_details.append(f'MW={mw:.1f}>500')
        if logp > 5:
            violations += 1
            violation_details.append(f'LogP={logp:.2f}>5')
        if hbd > 5:
            violations += 1
            violation_details.append(f'HBD={hbd}>5')
        if hba > 10:
            violations += 1
            violation_details.append(f'HBA={hba}>10')
        
        is_druglike = violations <= 1 and qed_score >= 0.3
        
        if 'qed' in self.models:
            fp = self._smiles_to_fingerprint_1024(smiles)
            if fp is not None:
                predicted_qed = self.models['qed'].predict(fp)[0]
            else:
                predicted_qed = qed_score
        else:
            predicted_qed = qed_score
        
        return {
            'qed_score': round(qed_score, 3),
            'predicted_qed': round(predicted_qed, 3),
            'lipinski_violations': violations,
            'violation_details': violation_details,
            'is_druglike': is_druglike,
            'descriptors': {
                'molecular_weight': round(mw, 2),
                'logP': round(logp, 2),
                'hbd': hbd,
                'hba': hba,
                'psa': round(psa, 2),
                'rotatable_bonds': rotatable_bonds
            }
        }
    
    def predict_all(self, smiles):
        results = {
            'smiles': smiles,
            'valid': Chem.MolFromSmiles(smiles) is not None
        }
        
        if not results['valid']:
            results['error'] = 'Invalid SMILES'
            return results
        
        results['bioactivity'] = self.predict_bioactivity(smiles)
        results['properties'] = self.predict_properties(smiles)
        results['toxicity'] = self.predict_toxicity(smiles)
        results['druglikeness'] = self.predict_druglikeness(smiles)
        
        return results
    
    def batch_predict(self, smiles_list):
        return [self.predict_all(smiles) for smiles in smiles_list]


if __name__ == "__main__":
    agent = PredictorAgent()
    
    test_smiles = [
        "CC(=O)Oc1ccccc1C(=O)O",
        "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
        "CC(C)Cc1ccc(cc1)C(C)C(=O)O"
    ]
    
    print("\n" + "=" * 60)
    print("Predictor Agent Test")
    print("=" * 60)
    
    for smiles in test_smiles:
        print(f"\nSMILES: {smiles}")
        results = agent.predict_all(smiles)
        print(f"  Bioactivity: pIC50={results['bioactivity'].get('pIC50', 'N/A')}")
        print(f"  Druglikeness: QED={results['druglikeness'].get('qed_score', 'N/A')}, Drug-like={results['druglikeness'].get('is_druglike', 'N/A')}")
        print(f"  Toxicity: Risk={results['toxicity'].get('risk_level', 'N/A')}")
