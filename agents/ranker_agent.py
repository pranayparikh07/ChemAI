import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
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


class RankerAgent:
    def __init__(self, predictor_agent=None):
        self.predictor = predictor_agent
        
        self.ranking_criteria = {
            'potency': {'weight': 0.30, 'higher_better': True},
            'druglikeness': {'weight': 0.25, 'higher_better': True},
            'safety': {'weight': 0.25, 'higher_better': True},
            'novelty': {'weight': 0.10, 'higher_better': True},
            'diversity': {'weight': 0.10, 'higher_better': True}
        }
        
        print("RankerAgent initialized")
    
    def set_predictor(self, predictor_agent):
        self.predictor = predictor_agent
    
    def set_criteria_weights(self, weights):
        for key, value in weights.items():
            if key in self.ranking_criteria:
                self.ranking_criteria[key]['weight'] = value
        
        total = sum(c['weight'] for c in self.ranking_criteria.values())
        for key in self.ranking_criteria:
            self.ranking_criteria[key]['weight'] /= total
    
    def score_molecule(self, smiles, known_molecules=None):
        if self.predictor is None:
            return {'error': 'Predictor agent not set'}
        
        predictions = self.predictor.predict_all(smiles)
        if not predictions.get('valid'):
            return {'smiles': smiles, 'total_score': 0, 'valid': False}
        
        scores = {}
        
        bioactivity = predictions.get('bioactivity', {})
        pIC50 = bioactivity.get('pIC50', 5)
        scores['potency'] = min(pIC50 / 9.0, 1.0)
        
        druglikeness = predictions.get('druglikeness', {})
        qed = druglikeness.get('qed_score', 0.5)
        violations = druglikeness.get('lipinski_violations', 0)
        is_druglike = druglikeness.get('is_druglike', False)
        scores['druglikeness'] = qed * (1 if is_druglike else 0.5)
        
        toxicity = predictions.get('toxicity', {})
        tox_prob = toxicity.get('toxicity_probability', 0.5)
        has_alerts = toxicity.get('has_structural_alerts', False)
        scores['safety'] = (1 - tox_prob) * (0.7 if has_alerts else 1.0)
        
        if known_molecules:
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024)
                similarities = []
                for known_smiles in known_molecules:
                    known_mol = Chem.MolFromSmiles(known_smiles)
                    if known_mol:
                        known_fp = AllChem.GetMorganFingerprintAsBitVect(known_mol, 2, nBits=1024)
                        sim = DataStructs.TanimotoSimilarity(fp, known_fp)
                        similarities.append(sim)
                max_sim = max(similarities) if similarities else 0
                scores['novelty'] = 1 - max_sim
            else:
                scores['novelty'] = 0
        else:
            scores['novelty'] = 0.5
        
        scores['diversity'] = 0.5
        
        total_score = sum(
            self.ranking_criteria[k]['weight'] * scores.get(k, 0)
            for k in self.ranking_criteria
        )
        
        return {
            'smiles': smiles,
            'total_score': round(total_score, 4),
            'component_scores': {k: round(v, 4) for k, v in scores.items()},
            'predictions': predictions,
            'valid': True
        }
    
    def rank_molecules(self, smiles_list, known_molecules=None, top_k=None):
        scored_molecules = []
        
        for smiles in smiles_list:
            result = self.score_molecule(smiles, known_molecules)
            if result.get('valid'):
                scored_molecules.append(result)
        
        scored_molecules.sort(key=lambda x: x['total_score'], reverse=True)
        
        for i, mol in enumerate(scored_molecules):
            mol['rank'] = i + 1
        
        if top_k:
            scored_molecules = scored_molecules[:top_k]
        
        return scored_molecules
    
    def calculate_diversity(self, smiles_list):
        fingerprints = []
        valid_smiles = []
        
        for smiles in smiles_list:
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024)
                fingerprints.append(fp)
                valid_smiles.append(smiles)
        
        if len(fingerprints) < 2:
            return {'diversity_score': 0, 'pairwise_similarities': []}
        
        similarities = []
        for i in range(len(fingerprints)):
            for j in range(i + 1, len(fingerprints)):
                sim = DataStructs.TanimotoSimilarity(fingerprints[i], fingerprints[j])
                similarities.append({
                    'mol1': valid_smiles[i],
                    'mol2': valid_smiles[j],
                    'similarity': round(sim, 4)
                })
        
        avg_similarity = np.mean([s['similarity'] for s in similarities])
        diversity_score = 1 - avg_similarity
        
        return {
            'diversity_score': round(diversity_score, 4),
            'average_similarity': round(avg_similarity, 4),
            'num_molecules': len(valid_smiles),
            'num_pairs': len(similarities)
        }
    
    def select_diverse_top_k(self, smiles_list, k=10, diversity_weight=0.3):
        scored = self.rank_molecules(smiles_list)
        
        if len(scored) <= k:
            return scored
        
        fingerprints = {}
        for mol in scored:
            m = Chem.MolFromSmiles(mol['smiles'])
            if m:
                fingerprints[mol['smiles']] = AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=1024)
        
        selected = [scored[0]]
        remaining = scored[1:]
        
        while len(selected) < k and remaining:
            best_idx = 0
            best_combined_score = 0
            
            for i, candidate in enumerate(remaining):
                if candidate['smiles'] not in fingerprints:
                    continue
                
                cand_fp = fingerprints[candidate['smiles']]
                min_dist = 1.0
                
                for sel in selected:
                    if sel['smiles'] in fingerprints:
                        sim = DataStructs.TanimotoSimilarity(cand_fp, fingerprints[sel['smiles']])
                        min_dist = min(min_dist, 1 - sim)
                
                combined = (1 - diversity_weight) * candidate['total_score'] + diversity_weight * min_dist
                
                if combined > best_combined_score:
                    best_combined_score = combined
                    best_idx = i
            
            selected.append(remaining.pop(best_idx))
        
        for i, mol in enumerate(selected):
            mol['diverse_rank'] = i + 1
        
        return selected
    
    def generate_report(self, ranked_molecules, report_name="Candidate Molecules Report"):
        report = []
        report.append("=" * 70)
        report.append(f"{report_name}")
        report.append("=" * 70)
        report.append(f"\nTotal candidates: {len(ranked_molecules)}")
        
        diversity = self.calculate_diversity([m['smiles'] for m in ranked_molecules])
        report.append(f"Set diversity: {diversity['diversity_score']:.4f}")
        
        report.append("\n" + "-" * 70)
        report.append("TOP CANDIDATES")
        report.append("-" * 70)
        
        for mol in ranked_molecules[:10]:
            report.append(f"\nRank #{mol.get('rank', 'N/A')}: {mol['smiles']}")
            report.append(f"  Total Score: {mol['total_score']:.4f}")
            report.append(f"  Component Scores:")
            for comp, score in mol.get('component_scores', {}).items():
                report.append(f"    - {comp}: {score:.4f}")
            
            preds = mol.get('predictions', {})
            bio = preds.get('bioactivity', {})
            drug = preds.get('druglikeness', {})
            tox = preds.get('toxicity', {})
            
            report.append(f"  Key Metrics:")
            report.append(f"    - pIC50: {bio.get('pIC50', 'N/A')}")
            report.append(f"    - QED: {drug.get('qed_score', 'N/A')}")
            report.append(f"    - Drug-like: {drug.get('is_druglike', 'N/A')}")
            report.append(f"    - Toxicity Risk: {tox.get('risk_level', 'N/A')}")
        
        report.append("\n" + "=" * 70)
        
        return "\n".join(report)


if __name__ == "__main__":
    from agents.predictor_agent import PredictorAgent
    
    predictor = PredictorAgent()
    ranker = RankerAgent(predictor)
    
    print("\n" + "=" * 60)
    print("Ranker Agent Test")
    print("=" * 60)
    
    test_molecules = [
        "CC(=O)Oc1ccccc1C(=O)O",
        "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
        "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
        "CC(C)NCC(O)c1ccc(O)c(O)c1",
        "COc1ccc2[nH]cc(CCNC(C)=O)c2c1"
    ]
    
    print("\n--- Ranking molecules ---")
    ranked = ranker.rank_molecules(test_molecules, top_k=5)
    
    for mol in ranked:
        print(f"\nRank #{mol['rank']}: {mol['smiles'][:40]}...")
        print(f"  Score: {mol['total_score']:.4f}")
    
    print("\n--- Diversity analysis ---")
    diversity = ranker.calculate_diversity(test_molecules)
    print(f"Diversity score: {diversity['diversity_score']:.4f}")
    
    print("\n--- Full report ---")
    report = ranker.generate_report(ranked)
    print(report[:2000])
