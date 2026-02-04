import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, QED
from rdkit import DataStructs
import random
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


class OptimizerAgent:
    def __init__(self, predictor_agent=None):
        self.predictor = predictor_agent
        
        self.default_weights = {
            'bioactivity': 0.35,
            'druglikeness': 0.25,
            'toxicity': 0.25,
            'novelty': 0.15
        }
        
        self.optimization_history = []
        print("OptimizerAgent initialized")
    
    def set_predictor(self, predictor_agent):
        self.predictor = predictor_agent
    
    def calculate_multi_objective_score(self, smiles, weights=None, reference_fps=None):
        if self.predictor is None:
            return {'error': 'Predictor agent not set'}
        
        weights = weights or self.default_weights
        
        predictions = self.predictor.predict_all(smiles)
        if not predictions.get('valid', False):
            return {'smiles': smiles, 'total_score': 0, 'valid': False}
        
        scores = {}
        
        bioactivity = predictions.get('bioactivity', {})
        pIC50 = bioactivity.get('pIC50', 5)
        scores['bioactivity'] = min(pIC50 / 10.0, 1.0)
        
        druglikeness = predictions.get('druglikeness', {})
        qed = druglikeness.get('qed_score', 0.5)
        violations = druglikeness.get('lipinski_violations', 0)
        scores['druglikeness'] = qed * (1 - 0.2 * violations)
        
        toxicity = predictions.get('toxicity', {})
        tox_prob = toxicity.get('toxicity_probability', 0.5)
        scores['toxicity'] = 1 - tox_prob
        
        if reference_fps:
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024)
                max_sim = max(DataStructs.TanimotoSimilarity(fp, ref_fp) for ref_fp in reference_fps)
                scores['novelty'] = 1 - max_sim
            else:
                scores['novelty'] = 0
        else:
            scores['novelty'] = 0.5
        
        total_score = sum(weights.get(k, 0) * v for k, v in scores.items())
        
        return {
            'smiles': smiles,
            'total_score': round(total_score, 4),
            'component_scores': {k: round(v, 4) for k, v in scores.items()},
            'predictions': predictions,
            'valid': True
        }
    
    def optimize_molecule(self, smiles, num_iterations=50, num_children=10, 
                          mutation_rate=0.3, weights=None):
        from agents.generator_agent import GeneratorAgent
        
        generator = GeneratorAgent()
        weights = weights or self.default_weights
        
        population = [smiles]
        best_score = 0
        best_molecule = smiles
        history = []
        
        initial_result = self.calculate_multi_objective_score(smiles, weights)
        if initial_result.get('valid'):
            best_score = initial_result['total_score']
            history.append({
                'iteration': 0,
                'best_score': best_score,
                'best_smiles': smiles
            })
        
        print(f"Starting optimization from: {smiles}")
        print(f"Initial score: {best_score:.4f}")
        
        for iteration in range(1, num_iterations + 1):
            new_population = []
            
            for parent in population[:5]:
                for _ in range(num_children):
                    if random.random() < mutation_rate:
                        child = generator.mutate_smiles(parent, num_mutations=random.randint(1, 2))
                    else:
                        child = parent
                    
                    if child:
                        new_population.append(child)
            
            if len(population) >= 2:
                for _ in range(num_children // 2):
                    p1, p2 = random.sample(population[:10], 2)
                    child = generator.crossover(p1, p2)
                    if child:
                        new_population.append(child)
            
            scored = []
            for mol_smiles in new_population:
                result = self.calculate_multi_objective_score(mol_smiles, weights)
                if result.get('valid'):
                    scored.append(result)
            
            scored.sort(key=lambda x: x['total_score'], reverse=True)
            
            population = [s['smiles'] for s in scored[:20]]
            
            if scored and scored[0]['total_score'] > best_score:
                best_score = scored[0]['total_score']
                best_molecule = scored[0]['smiles']
            
            if iteration % 10 == 0:
                print(f"Iteration {iteration}: Best score = {best_score:.4f}")
            
            history.append({
                'iteration': iteration,
                'best_score': best_score,
                'best_smiles': best_molecule,
                'population_size': len(population)
            })
        
        final_result = self.calculate_multi_objective_score(best_molecule, weights)
        
        print(f"\nOptimization complete!")
        print(f"Final score: {best_score:.4f}")
        print(f"Best molecule: {best_molecule}")
        
        return {
            'original': smiles,
            'optimized': best_molecule,
            'original_score': initial_result.get('total_score', 0),
            'optimized_score': best_score,
            'improvement': best_score - initial_result.get('total_score', 0),
            'history': history,
            'final_predictions': final_result
        }
    
    def pareto_optimization(self, molecules, objectives=None):
        objectives = objectives or ['bioactivity', 'druglikeness', 'toxicity']
        
        scored_molecules = []
        for smiles in molecules:
            result = self.calculate_multi_objective_score(smiles)
            if result.get('valid'):
                result['objectives'] = {obj: result['component_scores'].get(obj, 0) 
                                        for obj in objectives}
                scored_molecules.append(result)
        
        pareto_front = []
        
        for mol in scored_molecules:
            is_dominated = False
            
            for other in scored_molecules:
                if mol['smiles'] == other['smiles']:
                    continue
                
                all_worse_or_equal = all(
                    mol['objectives'][obj] <= other['objectives'][obj] 
                    for obj in objectives
                )
                any_strictly_worse = any(
                    mol['objectives'][obj] < other['objectives'][obj] 
                    for obj in objectives
                )
                
                if all_worse_or_equal and any_strictly_worse:
                    is_dominated = True
                    break
            
            if not is_dominated:
                pareto_front.append(mol)
        
        return pareto_front
    
    def suggest_improvements(self, smiles):
        if self.predictor is None:
            return {'error': 'Predictor agent not set'}
        
        predictions = self.predictor.predict_all(smiles)
        if not predictions.get('valid'):
            return {'error': 'Invalid molecule'}
        
        suggestions = []
        
        druglikeness = predictions.get('druglikeness', {})
        violations = druglikeness.get('violation_details', [])
        if violations:
            for v in violations:
                if 'MW' in v:
                    suggestions.append("Consider reducing molecular weight by removing heavy atoms or side chains")
                if 'LogP' in v:
                    suggestions.append("Add polar groups (OH, NH2) to reduce LogP")
                if 'HBD' in v:
                    suggestions.append("Reduce hydrogen bond donors")
                if 'HBA' in v:
                    suggestions.append("Reduce hydrogen bond acceptors")
        
        toxicity = predictions.get('toxicity', {})
        if toxicity.get('risk_level') in ['Medium', 'High']:
            alerts = toxicity.get('matched_alerts', [])
            suggestions.append(f"Address structural alerts: {', '.join(alerts[:3])}")
        
        bioactivity = predictions.get('bioactivity', {})
        if bioactivity.get('pIC50', 0) < 6:
            suggestions.append("Consider adding aromatic rings or hydrogen bond acceptors to improve binding")
        
        qed = druglikeness.get('qed_score', 0)
        if qed < 0.5:
            suggestions.append("Overall drug-likeness is low - consider a more balanced molecular profile")
        
        return {
            'smiles': smiles,
            'current_scores': predictions,
            'suggestions': suggestions
        }


if __name__ == "__main__":
    from agents.predictor_agent import PredictorAgent
    
    predictor = PredictorAgent()
    optimizer = OptimizerAgent(predictor)
    
    print("\n" + "=" * 60)
    print("Optimizer Agent Test")
    print("=" * 60)
    
    test_smiles = "CC(=O)Oc1ccccc1C(=O)O"
    
    print("\n--- Multi-objective scoring ---")
    result = optimizer.calculate_multi_objective_score(test_smiles)
    print(f"Total score: {result['total_score']}")
    print(f"Component scores: {result['component_scores']}")
    
    print("\n--- Improvement suggestions ---")
    suggestions = optimizer.suggest_improvements(test_smiles)
    for s in suggestions.get('suggestions', []):
        print(f"  • {s}")
