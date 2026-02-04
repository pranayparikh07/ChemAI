# ChemAI - Agentic Molecular Discovery System

## Project Overview
ChemAI is an intelligent, automated system designed to discover potential drug-like compounds using AI. It uses a network of software agents (Orchestrator, Generator, Predictor, Ranker, Optimizer) to evolve and evaluate molecules.

## System Architecture

### 1. Agents
- **OrchestratorAgent**: Manages the entire pipeline and coordinates other agents.
- **GeneratorAgent**: Creates new molecules using genetic algorithms (mutation/crossover) from seed compounds.
- **PredictorAgent**: Uses trained ML models to predict bioactivity, toxicity, and properties.
- **RankerAgent**: Scores molecules based on multiple criteria (QED, pIC50, Safety).
- **OptimizerAgent**: Iteratively improves a specific molecule to fix flaws (e.g., toxicity).

### 2. Machine Learning Models
The system relies on four key models located in `models/`:
- **Bioactivity Model**: Predicts pIC50 (potency) against targets.
- **Property Model**: Predicts MW, LogP, HBD, HBA, PSA.
- **Toxicity Model**: a Classifier that predicts structural alerts/toxicity.
- **Drug-Likeness Model**: Predicts QED (Quantitative Estimation of Drug-likeness) scores.

## Setup & Configuration for Laptop (i5, 8GB+)

To ensure the system runs smoothly on your hardware, the training scripts have been optimized:
- **Dataset Limit**: Capped at 50,000 records per model (down from 300k+).
- **Model Complexity**: Random Forest models use `n_estimators=30` and `max_depth=10` (lighter memory usage).
- **Balancing**: Datasets are automatically balanced to ensure fair learning even with smaller data.

## How to Run

### Step 1: Prepare the Models (One-time setup)
When you are ready to train the models, run the master script. This will take ~5-10 minutes on your system.
```bash
python train_all_models.py
```

### Step 2: Run the Discovery Pipeline
Once models are trained, you can run the full agent system.
```bash
python run_chemai.py --mode discover --smiles "CC(=O)Oc1ccccc1C(=O)O" --generations 3
```

### Step 3: Quick Evaluation (Demo)
To see the system in action without Training (if you just want to test imports):
```bash
python run_chemai.py --mode demo
```

## Directory Structure
- `agents/`: Source code for AI agents.
- `models/`: Training scripts for ML models.
- `trained_models/`: Where .joblib model files will be saved.
- `chembl_36/`: Database directory.
