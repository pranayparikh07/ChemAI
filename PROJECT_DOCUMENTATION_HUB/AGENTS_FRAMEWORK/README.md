# 🤖 AGENTS FRAMEWORK

**Agent Design, Implementation, and Integration**

---

## 📂 Files in This Folder

| File | Purpose | Location |
|------|---------|----------|
| `orchestrator_agent.py` | Main coordination agent | `agents/` |
| `generator_agent.py` | Molecule generation | `agents/` |
| `optimizer_agent.py` | RL optimization | `agents/` |
| `predictor_agent.py` | Property prediction | `agents/` |
| `ranker_agent.py` | Candidate ranking | `agents/` |

---

## 🔗 Agent Integration Map

```
┌─────────────────────────────────────┐
│   Orchestrator Agent (Shreya)       │
│   - Coordinates all agents          │
│   - Makes decisions                 │
│   - Logs experiments                │
└──────────────┬──────────────────────┘
        │ ┌────┴────┬────────┬─────────┐
        │ │         │        │         │
        ▼ ▼         ▼        ▼         ▼
   Generator  Optimizer  Predictor  Ranker
   (Pranay)  (Vishwa)    (All)     (All)
```

---

## 🎯 Agent Responsibilities

### **1. Orchestrator Agent (Shreya's Main Work)**

**Location**: `agents/orchestrator_agent.py`

**Key Methods**:
```python
class OrchestratorAgent:
    def __init__(self, models_dir='trained_models')
    def run_discovery_pipeline(self, seed_molecules, config)
    def make_autonomous_decision(self, context)
    def log_experiment(self, decision, data)
    def execute_generation_phase(self, config)
    def execute_optimization_phase(self, candidates, config)
    def execute_ranking_phase(self, molecules)
    def evaluate_convergence(self, scores_history)
    def save_experiment_results(self, filepath)
```

**Inputs**: Configuration, seed molecules
**Outputs**: Best candidates, experiment logs

---

### **2. Generator Agent (Pranay's Main Work)**

**Location**: `agents/generator_agent.py`

**Key Methods**:
```python
class GeneratorAgent:
    def __init__(self, seed_molecules=None)
    def generate_molecules(self, num_molecules=100)
    def apply_mutation(self, smiles, operation)
    def generate_with_constraints(self, seed, constraints)
    def calculate_validity(self, molecules)
    def calculate_novelty(self, molecules, reference_set)
    def calculate_diversity(self, molecules)
```

**Inputs**: Seed molecules, number to generate
**Outputs**: Novel SMILES strings, metrics

---

### **3. Optimizer Agent (Vishwa's Main Work)**

**Location**: `agents/optimizer_agent.py`

**Key Methods**:
```python
class OptimizerAgent:
    def __init__(self, predictor_agent=None)
    def optimize_molecule(self, smiles, config)
    def calculate_reward(self, smiles, weights)
    def apply_rl_policy(self, molecule, state)
    def run_rl_training(self, molecules, iterations)
    def compare_random_vs_rl(self, test_set, iterations)
```

**Inputs**: Molecules, optimization parameters
**Outputs**: Optimized molecules, reward scores

---

### **4. Predictor Agent (Property Evaluation)**

**Location**: `agents/predictor_agent.py`

**Key Methods**:
```python
class PredictorAgent:
    def __init__(self, models_dir)
    def predict_all(self, smiles)
    def predict_bioactivity(self, smiles)
    def predict_toxicity(self, smiles)
    def predict_druglikeness(self, smiles)
    def predict_properties(self, smiles)
    def calculate_multi_objective_score(self, smiles, weights)
```

**Inputs**: SMILES strings
**Outputs**: Multi-objective scores

---

### **5. Ranker Agent (Candidate Selection)**

**Location**: `agents/ranker_agent.py`

**Key Methods**:
```python
class RankerAgent:
    def __init__(self, predictor_agent)
    def rank_molecules(self, molecules, top_k=10)
    def calculate_diversity_score(self, molecules)
    def filter_by_criteria(self, molecules, criteria)
    def select_diverse_subset(self, molecules, num_select)
```

**Inputs**: Candidate molecules with scores
**Outputs**: Top-ranked, diverse candidates

---

## 🔄 Typical Execution Flow

```python
# Initialize
orchestrator = OrchestratorAgent()
predictor = orchestrator.predictor
generator = orchestrator.generator
optimizer = orchestrator.optimizer
ranker = orchestrator.ranker

# Configuration
config = {
    'num_generations': 5,
    'molecules_per_generation': 100,
    'optimization_iterations': 30,
    'top_k': 10
}

# Run pipeline
results = orchestrator.run_discovery_pipeline(
    seed_molecules=['SMILES1', 'SMILES2'],
    config=config
)

# Results contain:
# - best_candidates: Top 10 molecules
# - scores: All evaluation scores
# - experiment_log: All decisions & reasoning
# - metrics: Validity, novelty, diversity
```

---

## 📋 Agent Communication Protocol

### **Generator → Predictor**
```
Generator produces: SMILES list
↓
Predictor evaluates each SMILES
↓
Returns: Scores (bioactivity, toxicity, etc.)
```

### **Predictor → Optimizer**
```
Predictor scores: {smiles, score_dict}
↓
Optimizer selects promising molecules
↓
Applies RL modifications
↓
Returns: Improved molecules
```

### **All Agents → Orchestrator**
```
Each agent logs its output
↓
Orchestrator collects data
↓
Makes decision: Generate? Optimize? Rank?
↓
Logs to Neo4j for traceability
```

---

## ⚙️ Configuration Parameters

```python
config = {
    # Generation
    'num_generations': 5,                # How many times to iterate
    'molecules_per_generation': 100,     # How many to generate each time
    
    # Optimization
    'optimization_iterations': 30,       # RL training steps
    'lr': 0.001,                        # Learning rate
    'batch_size': 32,                   # Batch size
    
    # Ranking & Selection
    'top_k': 10,                        # Top molecules to keep
    'diversity_weight': 0.3,            # How much to favor diversity
    
    # Convergence
    'convergence_threshold': 0.001,     # When to stop
    'early_stopping': True,             # Stop if not improving
    'patience': 5,                      # Iterations before stop
    
    # Logging
    'verbose': True,                    # Print progress
    'save_results': True,               # Save to file
    'results_path': './results/'
}
```

---

## 🧪 Testing Agents Individually

### **Test Generator**
```python
generator = GeneratorAgent(seed_molecules=['CC(C)Cc1ccc(cc1)C(C)C(=O)O'])
molecules = generator.generate_molecules(num_molecules=50)
print(f"Generated: {len(molecules)} molecules")
print(f"Validity: {generator.calculate_validity(molecules)}")
```

### **Test Optimizer**
```python
optimizer = OptimizerAgent(predictor)
optimized = optimizer.optimize_molecule('SMILES_HERE', config)
print(f"Original score: {original_score}")
print(f"Optimized score: {optimized['score']}")
```

### **Test Full Pipeline**
```python
orchestrator = OrchestratorAgent()
results = orchestrator.run_discovery_pipeline(['seed1', 'seed2'])
print(results['best_candidates'])
```

---

## 🎛️ How to Extend Agents

### **Add New Agent**
1. Create: `agents/myagent_agent.py`
2. Inherit from base if needed
3. Implement required methods
4. Register in `orchestrator_agent.py`

### **Modify Generator** (Pranay)
- File: `agents/generator_agent.py`
- Add new mutation operations
- Implement GNN logic
- Add constraint checking

### **Modify Optimizer** (Vishwa)
- File: `agents/optimizer_agent.py`
- Update reward function
- Implement RL algorithm
- Add environment interaction

### **Modify Orchestrator** (Shreya)
- File: `agents/orchestrator_agent.py`
- Update decision logic
- Add new logging
- Implement causal analysis

---

## 📊 Output and Logging

### **Experiment Results Format**
```python
{
    'best_candidates': [
        {'smiles': '...', 'score': 0.92, 'source': 'generation_2'},
        {'smiles': '...', 'score': 0.88, 'source': 'optimization_3'},
        ...
    ],
    'metrics': {
        'validity': 0.95,
        'novelty': 0.87,
        'diversity': 0.73,
        'avg_improvement': 0.15
    },
    'decisions_log': [
        {
            'generation': 1,
            'decision': 'GENERATE',
            'input': [...],
            'output': [...],
            'timestamp': '2026-02-04 10:30:45',
            'reason': 'Initial generation phase'
        },
        ...
    ]
}
```

---

## 🔗 Integration with Neo4j

```python
# Shreya's code stores agent decisions in Neo4j
from graph_db.graph_loader import GraphLoader

loader = GraphLoader(neo4j_uri, username, password)
loader.create_experiment_node(results)
loader.link_molecules(results['best_candidates'])
loader.add_decision_trace(results['decisions_log'])
```

---

## 📚 Related Documentation

- `../../ARCHITECTURE_AND_DESIGN/` - System design
- `../../ML_MODELS/` - Predictor models
- `../../DATABASE_AND_GRAPHS/` - Neo4j integration
- `../../TESTING_VALIDATION/` - Agent tests

---

**Master Index**: Go back to `../../INDEX.md`
