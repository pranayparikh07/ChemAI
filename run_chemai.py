import os
import sys
import argparse

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from agents.orchestrator_agent import OrchestratorAgent


def main():
    parser = argparse.ArgumentParser(
        description='ChemAI - Agentic AI-Based Molecular Design and Optimization System'
    )
    
    parser.add_argument('--mode', type=str, default='demo',
                        choices=['demo', 'discover', 'evaluate', 'optimize', 'compare'],
                        help='Operation mode')
    parser.add_argument('--smiles', type=str, nargs='+',
                        help='SMILES string(s) for evaluation/optimization')
    parser.add_argument('--seeds', type=str, nargs='+',
                        help='Seed molecules for discovery pipeline')
    parser.add_argument('--generations', type=int, default=3,
                        help='Number of generations for discovery')
    parser.add_argument('--top-k', type=int, default=10,
                        help='Number of top candidates to return')
    parser.add_argument('--output', type=str, default='results.json',
                        help='Output file for results')
    parser.add_argument('--models-dir', type=str, default='trained_models',
                        help='Directory containing trained models')
    
    args = parser.parse_args()
    
    print("=" * 70)
    print("ChemAI - Agentic AI-Based Molecular Design and Optimization System")
    print("=" * 70)
    print(f"\nMode: {args.mode}")
    
    orchestrator = OrchestratorAgent(args.models_dir)
    
    if args.mode == 'demo':
        seeds = [
            "CC(=O)Oc1ccccc1C(=O)O",
            "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
            "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
            "CC(C)NCC(O)c1ccc(O)c(O)c1",
            "COc1ccc2[nH]cc(CCNC(C)=O)c2c1"
        ]
        
        print("\n[1] Quick molecule evaluation")
        orchestrator.quick_evaluate(seeds[0])
        
        print("\n[2] Molecule comparison")
        orchestrator.compare_molecules(seeds[:3])
        
        print("\n[3] Mini discovery pipeline")
        results = orchestrator.run_discovery_pipeline(
            seeds,
            config={
                'num_generations': 2,
                'molecules_per_generation': 20,
                'optimization_iterations': 20,
                'top_k': 5
            }
        )
        
        orchestrator.save_results(args.output)
    
    elif args.mode == 'evaluate':
        if not args.smiles:
            print("Error: --smiles required for evaluate mode")
            return
        
        for smiles in args.smiles:
            orchestrator.quick_evaluate(smiles)
    
    elif args.mode == 'compare':
        if not args.smiles or len(args.smiles) < 2:
            print("Error: At least 2 --smiles required for compare mode")
            return
        
        orchestrator.compare_molecules(args.smiles)
    
    elif args.mode == 'optimize':
        if not args.smiles:
            print("Error: --smiles required for optimize mode")
            return
        
        for smiles in args.smiles:
            orchestrator.optimize_lead(smiles)
    
    elif args.mode == 'discover':
        seeds = args.seeds or args.smiles
        if not seeds:
            print("Error: --seeds or --smiles required for discover mode")
            return
        
        results = orchestrator.run_discovery_pipeline(
            seeds,
            config={
                'num_generations': args.generations,
                'molecules_per_generation': 50,
                'optimization_iterations': 30,
                'top_k': args.top_k
            }
        )
        
        orchestrator.save_results(args.output)
    
    print("\n" + "=" * 70)
    print("Completed!")
    print("=" * 70)


if __name__ == "__main__":
    main()
