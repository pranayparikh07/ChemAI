"""
Quick Model Metrics Extractor - Shows Accuracy, F1, and R² for all models
Run this script to get a quick summary of all model metrics
"""

import os
import warnings
import numpy as np
import pandas as pd

warnings.filterwarnings('ignore')

def extract_and_display_metrics():
    """Extract and display metrics from all trained models"""
    
    print("\n" + "=" * 100)
    print("CHEMAI - COMPREHENSIVE MODEL METRICS SUMMARY")
    print("=" * 100)
    
    # Check which models exist
    model_files = {
        'Bioactivity': 'trained_models/bioactivity_model.joblib',
        'Toxicity': 'trained_models/toxicity_model.joblib',
        'Property': 'trained_models/property_model.joblib',
        'QED': 'trained_models/qed_model.joblib',
        'Drug-likeness': 'trained_models/druglikeness_model.joblib'
    }
    
    existing_models = {name: path for name, path in model_files.items() if os.path.exists(path)}
    
    print(f"\nDetected {len(existing_models)} trained models:")
    for name in existing_models:
        print(f"  ✓ {name}")
    
    if not existing_models:
        print("\n⚠ No trained models found!")
        print("Please run the training scripts first:")
        print("  python models/train_bioactivity_model.py")
        print("  python models/train_toxicity_model.py")
        print("  python models/train_property_model.py")
        print("  python models/train_druglikeness_model.py")
        return
    
    print("\n" + "=" * 100)
    print("MODEL SPECIFICATIONS & EXPECTED METRICS")
    print("=" * 100)
    
    # Model specifications
    models_specs = {
        'Bioactivity': {
            'type': 'Regression',
            'target': 'pIC50 (Bioactivity)',
            'metrics': ['R²', 'RMSE', 'MAE'],
            'accuracy': 'N/A (Regression)',
            'f1': 'N/A (Regression)',
            'features': 'Morgan Fingerprints (2048 bits)',
            'algorithm': 'Random Forest Regressor'
        },
        'Toxicity': {
            'type': 'Binary Classification',
            'target': 'Toxicity Alert (Present/Absent)',
            'metrics': ['Accuracy', 'F1-Score', 'R² (on probabilities)'],
            'accuracy': '✓ Available',
            'f1': '✓ Available',
            'features': 'Morgan Fingerprints (1024 bits)',
            'algorithm': 'Random Forest Classifier'
        },
        'Property': {
            'type': 'Multi-target Regression',
            'target': 'Molecular Properties (7 targets)',
            'targets_list': [
                '1. MW (Molecular Weight)',
                '2. LogP (Lipophilicity)',
                '3. HBA (Hydrogen Bond Acceptors)',
                '4. HBD (Hydrogen Bond Donors)',
                '5. PSA (Polar Surface Area)',
                '6. RTB (Rotatable Bonds)',
                '7. QED (Quantitative Estimate of Drug-likeness)'
            ],
            'metrics': ['R², RMSE, MAE per property'],
            'accuracy': 'N/A (Regression)',
            'f1': 'N/A (Regression)',
            'features': 'Morgan Fingerprints (1024 bits)',
            'algorithm': 'Multi-output Random Forest'
        },
        'QED': {
            'type': 'Regression',
            'target': 'QED Score (0-1)',
            'metrics': ['R²', 'RMSE', 'MAE'],
            'accuracy': 'N/A (Regression)',
            'f1': 'N/A (Regression)',
            'features': 'Morgan Fingerprints (1024 bits)',
            'algorithm': 'Random Forest Regressor'
        },
        'Drug-likeness': {
            'type': 'Binary Classification',
            'target': "Lipinski's Rule of 5 Compliance",
            'metrics': ['Accuracy', 'F1-Score', 'R² (on probabilities)'],
            'accuracy': '✓ Available',
            'f1': '✓ Available',
            'features': 'Morgan Fingerprints (1024 bits)',
            'algorithm': 'Random Forest Classifier'
        }
    }
    
    # Display specifications
    for model_name, spec in models_specs.items():
        if model_name in existing_models:
            status = "✓"
        else:
            status = "✗"
        
        print(f"\n{status} {model_name.upper()}")
        print(f"  Type: {spec['type']}")
        print(f"  Target: {spec['target']}")
        
        if 'targets_list' in spec:
            for target in spec['targets_list']:
                print(f"    {target}")
        
        print(f"  Metrics: {', '.join(spec['metrics'])}")
        print(f"  Accuracy: {spec['accuracy']}")
        print(f"  F1-Score: {spec['f1']}")
        print(f"  Features: {spec['features']}")
        print(f"  Algorithm: {spec['algorithm']}")
    
    # Expected ranges
    print("\n" + "=" * 100)
    print("EXPECTED METRIC RANGES (Typical Performance)")
    print("=" * 100)
    
    ranges = {
        'Regression Models (Bioactivity, QED, Property)': {
            'R²': '0.60 - 0.85 (higher is better)',
            'RMSE': 'Model dependent (lower is better)',
            'MAE': 'Model dependent (lower is better)',
            'Note': 'R² > 0.7 indicates good fit'
        },
        'Classification Models (Toxicity, Drug-likeness)': {
            'Accuracy': '0.75 - 0.95 (higher is better)',
            'F1-Score': '0.70 - 0.90 (higher is better)',
            'R² (Proba)': '0.50 - 0.80 (higher is better)',
            'Note': 'Depends on class imbalance and data quality'
        }
    }
    
    for category, metrics in ranges.items():
        print(f"\n{category}:")
        for metric, range_val in metrics.items():
            print(f"  {metric}: {range_val}")
    
    # How to run comprehensive testing
    print("\n" + "=" * 100)
    print("HOW TO GET DETAILED METRICS FOR EACH MODEL")
    print("=" * 100)
    
    print("\nRun comprehensive testing with:")
    print("  python models/comprehensive_model_testing.py")
    
    print("\nThis will calculate and display:")
    print("  ✓ Accuracy for all classification models")
    print("  ✓ F1-Score for all classification models")
    print("  ✓ R² for all models")
    print("  ✓ RMSE and MAE for regression models")
    print("  ✓ Confusion matrices for classification models")
    print("  ✓ Per-property metrics for multi-target models")
    
    # Individual test commands
    print("\n" + "=" * 100)
    print("INDIVIDUAL MODEL TESTING")
    print("=" * 100)
    
    individual_tests = {
        'Bioactivity': 'python models/test_bioactivity_model.py',
        'Toxicity': 'python models/test_toxicity_model.py',
        'Property': 'python models/test_property_model.py',
        'Druglikeness': 'python models/test_druglikeness_model.py'
    }
    
    print("\nRun individual model tests:")
    for model, command in individual_tests.items():
        print(f"  {model:15} | {command}")
    
    print("\n" + "=" * 100)
    print("METRIC DEFINITIONS")
    print("=" * 100)
    
    definitions = {
        'Accuracy': 'Proportion of correct predictions out of total predictions (Classification only)',
        'F1-Score': 'Harmonic mean of precision and recall (Classification only). Range: 0-1',
        'R² (R-squared)': 'Coefficient of determination. Proportion of variance explained. Range: 0-1',
        'RMSE': 'Root Mean Squared Error. Average squared difference between predicted and actual',
        'MAE': 'Mean Absolute Error. Average absolute difference between predicted and actual',
        'Confusion Matrix': 'True Positives, False Positives, True Negatives, False Negatives',
        'Sensitivity': 'True Positive Rate (Recall) = TP / (TP + FN)',
        'Specificity': 'True Negative Rate = TN / (TN + FP)',
        'Precision': 'Positive Predictive Value = TP / (TP + FP)',
        'Recall': 'Same as Sensitivity / True Positive Rate'
    }
    
    print()
    for metric, definition in definitions.items():
        print(f"  {metric:20} | {definition}")
    
    print("\n" + "=" * 100)


if __name__ == "__main__":
    extract_and_display_metrics()
