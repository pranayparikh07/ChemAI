import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from models.train_bioactivity_model import main as train_bioactivity
from models.train_property_model import main as train_properties
from models.train_toxicity_model import main as train_toxicity
from models.train_druglikeness_model import main as train_druglikeness


def main():
    print("=" * 70)
    print("ChemAI - Complete Model Training Pipeline")
    print("Agentic AI-Based Molecular Design and Optimization System")
    print("=" * 70)
    
    start_time = time.time()
    results = {}
    
    print("\n" + "=" * 70)
    print("[1/4] TRAINING BIOACTIVITY MODEL (pIC50 Prediction)")
    print("=" * 70)
    try:
        bio_model, bio_results = train_bioactivity()
        results['bioactivity'] = {'status': 'success', 'metrics': bio_results}
        print("✓ Bioactivity model trained successfully")
    except Exception as e:
        results['bioactivity'] = {'status': 'failed', 'error': str(e)}
        print(f"✗ Bioactivity model training failed: {e}")
    
    print("\n" + "=" * 70)
    print("[2/4] TRAINING PROPERTY PREDICTION MODEL")
    print("=" * 70)
    try:
        prop_model, prop_results = train_properties()
        results['properties'] = {'status': 'success', 'metrics': prop_results}
        print("✓ Property model trained successfully")
    except Exception as e:
        results['properties'] = {'status': 'failed', 'error': str(e)}
        print(f"✗ Property model training failed: {e}")
    
    print("\n" + "=" * 70)
    print("[3/4] TRAINING TOXICITY PREDICTION MODEL")
    print("=" * 70)
    try:
        tox_model, alerts_df, tox_results = train_toxicity()
        results['toxicity'] = {'status': 'success', 'metrics': tox_results}
        print("✓ Toxicity model trained successfully")
    except Exception as e:
        results['toxicity'] = {'status': 'failed', 'error': str(e)}
        print(f"✗ Toxicity model training failed: {e}")
    
    print("\n" + "=" * 70)
    print("[4/4] TRAINING DRUG-LIKENESS MODEL")
    print("=" * 70)
    try:
        qed_model, druglike_model, dl_results = train_druglikeness()
        results['druglikeness'] = {'status': 'success', 'metrics': dl_results}
        print("✓ Drug-likeness model trained successfully")
    except Exception as e:
        results['druglikeness'] = {'status': 'failed', 'error': str(e)}
        print(f"✗ Drug-likeness model training failed: {e}")
    
    elapsed_time = time.time() - start_time
    
    print("\n" + "=" * 70)
    print("TRAINING SUMMARY")
    print("=" * 70)
    print(f"\nTotal training time: {elapsed_time/60:.1f} minutes")
    print("\nModel Status:")
    for model_name, result in results.items():
        status = "✓" if result['status'] == 'success' else "✗"
        print(f"  {status} {model_name}: {result['status']}")
    
    print("\nTrained models saved to: trained_models/")
    print("  - bioactivity_model.joblib (pIC50 prediction)")
    print("  - property_model.joblib (molecular properties)")
    print("  - toxicity_model.joblib (structural alerts)")
    print("  - qed_model.joblib (QED score)")
    print("  - druglikeness_model.joblib (drug-likeness)")
    
    print("\n" + "=" * 70)
    print("Model training completed! Ready for agent system.")
    print("=" * 70)
    
    return results


if __name__ == "__main__":
    results = main()
