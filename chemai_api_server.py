#!/usr/bin/env python
"""
ChemAI Production API Server
Drug Discovery AI System - Production Level
"""

from flask import Flask, request, jsonify
from flask_cors import CORS
import joblib
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
import logging
import os
import json
from datetime import datetime
from pathlib import Path
import traceback

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('logs/chemai_server.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# Initialize Flask app
app = Flask(__name__)
CORS(app)

# Configuration
CONFIG = {
    'MODEL_DIR': 'trained_models',
    'MAX_MOLECULES': 100,
    'FINGERPRINT_RADIUS': 2,
    'FINGERPRINT_BITS': 1024,
}

class ModelManager:
    """Manages all ML models"""
    
    def __init__(self):
        self.models = {}
        self.config = {}
        self.load_models()
    
    def load_models(self):
        """Load all trained models"""
        try:
            logger.info("Loading models...")
            
            # Property model
            self.models['property'] = joblib.load(
                os.path.join(CONFIG['MODEL_DIR'], 'property_model.joblib')
            )
            self.models['property_scaler'] = joblib.load(
                os.path.join(CONFIG['MODEL_DIR'], 'property_scaler.joblib')
            )
            
            # QED model
            self.models['qed'] = joblib.load(
                os.path.join(CONFIG['MODEL_DIR'], 'qed_model.joblib')
            )
            
            # Druglikeness model
            self.models['druglikeness'] = joblib.load(
                os.path.join(CONFIG['MODEL_DIR'], 'druglikeness_model.joblib')
            )
            
            # Bioactivity model
            self.models['bioactivity'] = joblib.load(
                os.path.join(CONFIG['MODEL_DIR'], 'bioactivity_model.joblib')
            )
            
            # Toxicity model
            if os.path.exists(os.path.join(CONFIG['MODEL_DIR'], 'toxicity_model.joblib')):
                self.models['toxicity'] = joblib.load(
                    os.path.join(CONFIG['MODEL_DIR'], 'toxicity_model.joblib')
                )
            
            logger.info(f"✓ Successfully loaded {len(self.models)} models")
            return True
        except Exception as e:
            logger.error(f"✗ Failed to load models: {e}")
            return False
    
    def compute_fingerprint(self, smiles):
        """Compute Morgan fingerprint for a SMILES string"""
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return None
            
            fp = AllChem.GetMorganFingerprintAsBitVect(
                mol, 
                CONFIG['FINGERPRINT_RADIUS'],
                nBits=CONFIG['FINGERPRINT_BITS']
            )
            return np.array(fp)
        except Exception as e:
            logger.warning(f"Failed to compute fingerprint for {smiles}: {e}")
            return None
    
    def predict_properties(self, fingerprint):
        """Predict molecular properties"""
        try:
            properties = {}
            
            # Reshape for prediction
            X = fingerprint.reshape(1, -1)
            
            # Scale features if scaler available
            if 'property_scaler' in self.models:
                X = self.models['property_scaler'].transform(X)
            
            # Predict properties
            predictions = self.models['property'].predict(X)[0]
            property_names = ['mw_freebase', 'alogp', 'hba', 'hbd', 'psa', 'rtb', 'qed_weighted']
            
            for name, value in zip(property_names, predictions):
                properties[name] = float(value)
            
            return properties
        except Exception as e:
            logger.error(f"Error predicting properties: {e}")
            return None
    
    def predict_qed(self, fingerprint):
        """Predict QED score"""
        try:
            X = fingerprint.reshape(1, -1)
            qed_score = float(self.models['qed'].predict(X)[0])
            return qed_score
        except Exception as e:
            logger.error(f"Error predicting QED: {e}")
            return None
    
    def predict_druglikeness(self, fingerprint):
        """Predict drug-likeness (Lipinski Rule of 5)"""
        try:
            X = fingerprint.reshape(1, -1)
            prediction = self.models['druglikeness'].predict(X)[0]
            probability = float(self.models['druglikeness'].predict_proba(X)[0][1])
            
            return {
                'is_druglike': int(prediction),
                'probability': probability
            }
        except Exception as e:
            logger.error(f"Error predicting druglikeness: {e}")
            return None
    
    def predict_bioactivity(self, fingerprint):
        """Predict bioactivity (pIC50)"""
        try:
            X = fingerprint.reshape(1, -1)
            pic50 = float(self.models['bioactivity'].predict(X)[0])
            return pic50
        except Exception as e:
            logger.error(f"Error predicting bioactivity: {e}")
            return None
    
    def predict_toxicity(self, fingerprint):
        """Predict toxicity risk"""
        try:
            if 'toxicity' not in self.models:
                return None
            
            X = fingerprint.reshape(1, -1)
            prediction = self.models['toxicity'].predict(X)[0]
            probability = float(self.models['toxicity'].predict_proba(X)[0][1])
            
            return {
                'is_toxic': int(prediction),
                'probability': probability
            }
        except Exception as e:
            logger.error(f"Error predicting toxicity: {e}")
            return None


# Initialize models
model_manager = ModelManager()


@app.route('/health', methods=['GET'])
def health_check():
    """Health check endpoint"""
    return jsonify({
        'status': 'healthy',
        'timestamp': datetime.now().isoformat(),
        'models_loaded': len(model_manager.models)
    })


@app.route('/api/v1/predict', methods=['POST'])
def predict():
    """
    Predict all properties for a molecule
    
    Request:
    {
        "smiles": "CC(C)Cc1ccc(cc1)C(C)C(O)=O",
        "include_toxicity": true
    }
    """
    try:
        data = request.get_json()
        
        if not data or 'smiles' not in data:
            return jsonify({'error': 'Missing SMILES string'}), 400
        
        smiles = data['smiles']
        include_toxicity = data.get('include_toxicity', True)
        
        logger.info(f"Prediction request for: {smiles}")
        
        # Validate SMILES
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return jsonify({'error': 'Invalid SMILES string'}), 400
        
        # Compute fingerprint
        fp = model_manager.compute_fingerprint(smiles)
        if fp is None:
            return jsonify({'error': 'Failed to compute fingerprint'}), 400
        
        # Get predictions
        result = {
            'smiles': smiles,
            'timestamp': datetime.now().isoformat(),
            'predictions': {
                'properties': model_manager.predict_properties(fp),
                'qed': model_manager.predict_qed(fp),
                'druglikeness': model_manager.predict_druglikeness(fp),
                'bioactivity': model_manager.predict_bioactivity(fp),
            }
        }
        
        if include_toxicity:
            result['predictions']['toxicity'] = model_manager.predict_toxicity(fp)
        
        logger.info(f"✓ Prediction successful for {smiles}")
        return jsonify(result), 200
        
    except Exception as e:
        logger.error(f"✗ Prediction error: {e}\n{traceback.format_exc()}")
        return jsonify({'error': str(e)}), 500


@app.route('/api/v1/batch_predict', methods=['POST'])
def batch_predict():
    """
    Batch predict for multiple molecules
    
    Request:
    {
        "smiles_list": ["SMILES1", "SMILES2", ...],
        "include_toxicity": true
    }
    """
    try:
        data = request.get_json()
        
        if not data or 'smiles_list' not in data:
            return jsonify({'error': 'Missing smiles_list'}), 400
        
        smiles_list = data['smiles_list']
        include_toxicity = data.get('include_toxicity', True)
        
        if len(smiles_list) > CONFIG['MAX_MOLECULES']:
            return jsonify({'error': f'Maximum {CONFIG["MAX_MOLECULES"]} molecules allowed'}), 400
        
        logger.info(f"Batch prediction request for {len(smiles_list)} molecules")
        
        results = []
        for smiles in smiles_list:
            fp = model_manager.compute_fingerprint(smiles)
            if fp is None:
                results.append({'smiles': smiles, 'error': 'Invalid SMILES'})
                continue
            
            result = {
                'smiles': smiles,
                'predictions': {
                    'properties': model_manager.predict_properties(fp),
                    'qed': model_manager.predict_qed(fp),
                    'druglikeness': model_manager.predict_druglikeness(fp),
                    'bioactivity': model_manager.predict_bioactivity(fp),
                }
            }
            
            if include_toxicity:
                result['predictions']['toxicity'] = model_manager.predict_toxicity(fp)
            
            results.append(result)
        
        logger.info(f"✓ Batch prediction completed for {len(results)} molecules")
        return jsonify({'results': results, 'count': len(results)}), 200
        
    except Exception as e:
        logger.error(f"✗ Batch prediction error: {e}\n{traceback.format_exc()}")
        return jsonify({'error': str(e)}), 500


@app.route('/api/v1/models', methods=['GET'])
def get_models_info():
    """Get information about loaded models"""
    return jsonify({
        'models_loaded': list(model_manager.models.keys()),
        'count': len(model_manager.models),
        'config': CONFIG
    }), 200


@app.errorhandler(404)
def not_found(error):
    """Handle 404 errors"""
    return jsonify({'error': 'Endpoint not found'}), 404


@app.errorhandler(500)
def internal_error(error):
    """Handle 500 errors"""
    logger.error(f"Internal server error: {error}")
    return jsonify({'error': 'Internal server error'}), 500


def main():
    """Start the server"""
    # Create logs directory
    os.makedirs('logs', exist_ok=True)
    
    logger.info("=" * 70)
    logger.info("ChemAI Production API Server Starting")
    logger.info("=" * 70)
    logger.info(f"Models loaded: {len(model_manager.models)}")
    logger.info(f"Configuration: {CONFIG}")
    logger.info("=" * 70)
    
    # Run Flask app
    app.run(
        host='0.0.0.0',
        port=5000,
        debug=False,
        threaded=True
    )


if __name__ == '__main__':
    main()
