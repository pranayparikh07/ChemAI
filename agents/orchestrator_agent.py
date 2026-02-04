import os
import sys
import json
from datetime import datetime

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from agents.predictor_agent import PredictorAgent
from agents.generator_agent import GeneratorAgent
from agents.optimizer_agent import OptimizerAgent
from agents.ranker_agent import RankerAgent
from web_dashboard import init_dashboard, add_molecule_to_dashboard, finalize_dashboard


class OrchestratorAgent:
    def __init__(self, models_dir='trained_models'):
        print("=" * 70)
        print("Initializing Orchestrator Agent")
        print("=" * 70)
        
        self.predictor = PredictorAgent(models_dir)
        self.generator = GeneratorAgent()
        self.optimizer = OptimizerAgent(self.predictor)
        self.ranker = RankerAgent(self.predictor)
        
        self.workflow_history = []
        self.best_candidates = []
        
        print("All agents initialized and connected")
        print("=" * 70)
    
    def run_discovery_pipeline(self, seed_molecules, config=None):
        config = config or {}
        num_generations = config.get('num_generations', 3)
        molecules_per_generation = config.get('molecules_per_generation', 50)
        optimization_iterations = config.get('optimization_iterations', 30)
        top_k = config.get('top_k', 10)
        diversity_weight = config.get('diversity_weight', 0.3)
        
        # Initialize dashboard
        dashboard = init_dashboard()
        
        print("\n" + "=" * 70)
        print("DRUG DISCOVERY PIPELINE")
        print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        print("=" * 70)
        
        print(f"\nConfiguration:")
        print(f"  - Seed molecules: {len(seed_molecules)}")
        print(f"  - Generations: {num_generations}")
        print(f"  - Molecules per generation: {molecules_per_generation}")
        print(f"  - Top candidates to select: {top_k}")
        
        self.generator.set_seed_molecules(seed_molecules)
        
        all_candidates = list(seed_molecules)
        generation_results = []
        molecule_counter = 0
        
        for gen in range(num_generations):
            dashboard.set_generation(gen + 1)
            
            print(f"\n{'='*50}")
            print(f"GENERATION {gen + 1}/{num_generations}")
            print(f"{'='*50}")
            
            print("\n[1] Generating new molecules...")
            generated = self.generator.generate_diverse_set(
                num_molecules=molecules_per_generation,
                min_diversity=0.3
            )
            print(f"    Generated {len(generated)} molecules")
            
            # Track each generated molecule
            for smiles in generated:
                molecule_counter += 1
                mol_name = f"G{gen+1}_M{molecule_counter}"
                predictions = self.predictor.predict_all(smiles)
                add_molecule_to_dashboard(
                    smiles=smiles,
                    name=mol_name,
                    predictions=predictions,
                    generation=gen + 1
                )
            
            print("\n[2] Predicting properties...")
            predictions = self.predictor.batch_predict(generated)
            valid_molecules = [p['smiles'] for p in predictions if p.get('valid')]
            print(f"    Valid molecules: {len(valid_molecules)}")
            
            print("\n[3] Ranking molecules...")
            ranked = self.ranker.rank_molecules(valid_molecules, known_molecules=seed_molecules)
            
            print("\n[4] Selecting top candidates for optimization...")
            top_candidates = ranked[:10]
            
            print("\n[5] Optimizing top candidates...")
            optimized = []
            for i, candidate in enumerate(top_candidates[:3]):
                print(f"\n    Optimizing candidate {i+1}/3: {candidate['smiles'][:30]}...")
                opt_result = self.optimizer.optimize_molecule(
                    candidate['smiles'],
                    num_iterations=optimization_iterations,
                    num_children=5
                )
                optimized.append(opt_result)
                print(f"    Improvement: {opt_result['improvement']:.4f}")
            
            optimized_smiles = [o['optimized'] for o in optimized]
            all_candidates.extend(valid_molecules)
            all_candidates.extend(optimized_smiles)
            
            self.generator.set_seed_molecules(
                [r['smiles'] for r in ranked[:20]]
            )
            
            gen_result = {
                'generation': gen + 1,
                'generated': len(generated),
                'valid': len(valid_molecules),
                'top_score': ranked[0]['total_score'] if ranked else 0,
                'optimized': len(optimized)
            }
            generation_results.append(gen_result)
            
            print(f"\n    Generation {gen + 1} Summary:")
            print(f"    - Best score: {gen_result['top_score']:.4f}")
        
        print("\n" + "=" * 70)
        print("FINAL RANKING")
        print("=" * 70)
        
        all_candidates = list(set(all_candidates))
        print(f"\nTotal unique candidates: {len(all_candidates)}")
        
        final_ranked = self.ranker.select_diverse_top_k(
            all_candidates, k=top_k, diversity_weight=diversity_weight
        )
        
        self.best_candidates = final_ranked
        
        report = self.ranker.generate_report(final_ranked, "Final Discovery Results")
        print(report)
        
        pipeline_result = {
            'timestamp': datetime.now().isoformat(),
            'config': config,
            'seed_molecules': seed_molecules,
            'generation_results': generation_results,
            'total_candidates_evaluated': len(all_candidates),
            'top_candidates': final_ranked,
            'report': report
        }
        
        self.workflow_history.append(pipeline_result)
        
        # Finalize and save dashboard
        finalize_dashboard()
        
        return pipeline_result
    
    def quick_evaluate(self, smiles):
        print(f"\nQuick evaluation: {smiles}")
        
        dashboard = init_dashboard()
        
        predictions = self.predictor.predict_all(smiles)
        score = self.ranker.score_molecule(smiles)
        suggestions = self.optimizer.suggest_improvements(smiles)
        
        result = {
            'smiles': smiles,
            'valid': predictions.get('valid', False),
            'predictions': predictions,
            'ranking_score': score.get('total_score', 0),
            'component_scores': score.get('component_scores', {}),
            'improvement_suggestions': suggestions.get('suggestions', [])
        }
        
        # Add to dashboard
        if result['valid']:
            add_molecule_to_dashboard(
                smiles=smiles,
                name=f"Evaluated_Compound",
                predictions=predictions.get('druglikeness', {}),
                generation=1
            )
        
        print(f"\n  Valid: {result['valid']}")
        if result['valid']:
            print(f"  Ranking Score: {result['ranking_score']:.4f}")
            print(f"  pIC50: {predictions.get('bioactivity', {}).get('pIC50', 'N/A')}")
            print(f"  QED: {predictions.get('druglikeness', {}).get('qed_score', 'N/A')}")
            print(f"  Drug-like: {predictions.get('druglikeness', {}).get('is_druglike', 'N/A')}")
            print(f"  Toxicity Risk: {predictions.get('toxicity', {}).get('risk_level', 'N/A')}")
            
            if result['improvement_suggestions']:
                print(f"\n  Suggestions:")
                for s in result['improvement_suggestions'][:3]:
                    print(f"    • {s}")
        
        return result
    
    def compare_molecules(self, smiles_list):
        print("\n" + "=" * 70)
        print("MOLECULE COMPARISON")
        print("=" * 70)
        
        results = []
        for smiles in smiles_list:
            result = self.quick_evaluate(smiles)
            results.append(result)
        
        ranked = self.ranker.rank_molecules(smiles_list)
        
        print("\n" + "-" * 70)
        print("COMPARISON SUMMARY")
        print("-" * 70)
        print(f"{'Rank':<6} {'Score':<10} {'pIC50':<8} {'QED':<8} {'Drug-like':<12} {'SMILES':<30}")
        print("-" * 70)
        
        for mol in ranked:
            preds = mol.get('predictions', {})
            bio = preds.get('bioactivity', {})
            drug = preds.get('druglikeness', {})
            
            print(f"{mol['rank']:<6} {mol['total_score']:<10.4f} "
                  f"{bio.get('pIC50', 'N/A'):<8} {drug.get('qed_score', 'N/A'):<8} "
                  f"{str(drug.get('is_druglike', 'N/A')):<12} {mol['smiles'][:30]}")
        
        return ranked
    
    def optimize_lead(self, smiles, iterations=50):
        print("\n" + "=" * 70)
        print("LEAD OPTIMIZATION")
        print("=" * 70)
        
        print(f"Lead compound: {smiles}")
        
        initial = self.quick_evaluate(smiles)
        
        result = self.optimizer.optimize_molecule(
            smiles,
            num_iterations=iterations,
            num_children=10
        )
        
        print("\n" + "-" * 70)
        print("OPTIMIZATION RESULTS")
        print("-" * 70)
        print(f"Original: {result['original']}")
        print(f"Optimized: {result['optimized']}")
        print(f"Score improvement: {result['original_score']:.4f} → {result['optimized_score']:.4f}")
        print(f"Delta: +{result['improvement']:.4f}")
        
        return result
    
    def save_results(self, filepath):
        results = {
            'workflow_history': self.workflow_history,
            'best_candidates': [
                {
                    'smiles': c['smiles'],
                    'score': c['total_score'],
                    'rank': c.get('rank', 0)
                }
                for c in self.best_candidates
            ]
        }
        
        with open(filepath, 'w') as f:
            json.dump(results, f, indent=2, default=str)
        
        print(f"Results saved to {filepath}")


def main():
    orchestrator = OrchestratorAgent()
    
    seed_molecules = [
        "CC(=O)Oc1ccccc1C(=O)O",
        "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
        "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
        "CC(C)NCC(O)c1ccc(O)c(O)c1",
        "COc1ccc2[nH]cc(CCNC(C)=O)c2c1"
    ]
    
    print("\n" + "=" * 70)
    print("ChemAI - Agentic Molecular Design System")
    print("=" * 70)
    
    print("\n[Demo 1] Quick evaluation of a single molecule")
    orchestrator.quick_evaluate(seed_molecules[0])
    
    print("\n[Demo 2] Comparing multiple molecules")
    orchestrator.compare_molecules(seed_molecules[:3])
    
    print("\n[Demo 3] Full discovery pipeline (reduced for demo)")
    results = orchestrator.run_discovery_pipeline(
        seed_molecules,
        config={
            'num_generations': 2,
            'molecules_per_generation': 20,
            'optimization_iterations': 20,
            'top_k': 5
        }
    )
    
    orchestrator.save_results('discovery_results.json')
    
    print("\n" + "=" * 70)
    print("Demo completed!")
    print("=" * 70)


if __name__ == "__main__":
    main()
