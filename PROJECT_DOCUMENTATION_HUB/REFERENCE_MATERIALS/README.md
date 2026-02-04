# 📚 REFERENCE MATERIALS

**Knowledge Base, Quick References, and Lookup Guides**

---

## 📂 Files in This Folder

| File | Purpose |
|------|---------|
| `METRICS_GUIDE.py` | How to calculate metrics |
| `METRICS_REFERENCE.md` | Metric definitions |
| `MODEL_METRICS_QUICK_VIEW.txt` | Quick reference |
| `QUICK_REFERENCE.md` | Command cheat sheet |
| `FILE_INDEX.md` | Complete file listing |

---

## 📊 Metrics Reference

### **Classification Metrics**

#### **Accuracy**
- Definition: (TP + TN) / (TP + TN + FP + FN)
- Range: 0-1 (0-100%)
- Interpretation: Percentage of correct predictions
- Use Case: Toxicity classification, drug-likeness
- **Good Value**: > 0.85

#### **Precision**
- Definition: TP / (TP + FP)
- Range: 0-1
- Interpretation: Of positive predictions, how many were correct?
- **Good Value**: > 0.80

#### **Recall (Sensitivity)**
- Definition: TP / (TP + FN)
- Range: 0-1
- Interpretation: Of actual positives, how many did we find?
- **Good Value**: > 0.80

#### **F1-Score**
- Definition: 2 × (Precision × Recall) / (Precision + Recall)
- Range: 0-1
- Interpretation: Harmonic mean of precision & recall
- **Good Value**: > 0.85

#### **ROC-AUC**
- Definition: Area under Receiver Operating Characteristic curve
- Range: 0-1
- Interpretation: Model's ability to distinguish between classes
- **Good Value**: > 0.90

---

### **Regression Metrics**

#### **R² Score (Coefficient of Determination)**
- Definition: 1 - (SS_res / SS_tot)
- Range: -∞ to 1
- Interpretation: Variance explained by model
- Use Case: Bioactivity prediction, property modeling
- **Good Value**: > 0.80

#### **Root Mean Squared Error (RMSE)**
- Definition: √(Σ(y_pred - y_true)² / n)
- Range: 0 to ∞
- Interpretation: Average prediction error
- Units: Same as target variable
- **Good Value**: < 1.0

#### **Mean Absolute Error (MAE)**
- Definition: Σ|y_pred - y_true| / n
- Range: 0 to ∞
- Interpretation: Average absolute deviation
- Units: Same as target variable
- **Good Value**: < 0.8

#### **Spearman Correlation**
- Definition: Rank-based correlation coefficient
- Range: -1 to 1
- Interpretation: Monotonic relationship strength
- **Good Value**: > 0.80

---

## 🧬 Molecular Descriptors

### **Common Descriptors**

| Descriptor | Abbr | Range | Significance |
|-----------|------|-------|-------------|
| Molecular Weight | MW | 100-500 | Drug-like range |
| LogP | LogP | -1 to 5 | Lipophilicity |
| QED | QED | 0-1 | Drug-likeness |
| TPSA | TPSA | 0-120 | Polar surface area |
| H-bond Donors | HBD | 0-5 | H-bonding capacity |
| H-bond Acceptors | HBA | 0-10 | H-bonding capacity |
| Rotatable Bonds | RotBonds | 0-15 | Flexibility |
| Aromatic Rings | AromRings | 0-6 | Aromaticity |

---

## 🔬 Chemistry Rules

### **Lipinski's Rule of Five**
- Molecular weight: ≤ 500 Da
- LogP: ≤ 5
- H-bond donors: ≤ 5
- H-bond acceptors: ≤ 10

### **Veber's Rule**
- Rotatable bonds: ≤ 10
- TPSA: 20-130 Ų

### **Lead-likeness**
- Molecular weight: ≤ 350 Da
- LogP: ≤ 3
- H-bond donors: ≤ 3
- H-bond acceptors: ≤ 3-4

---

## 🔧 Quick Command Reference

### **Training**
```bash
# Train all models
python train_all_models.py

# Train specific model
python models/train_bioactivity_model.py

# Train with custom parameters
python models/train_bioactivity_model.py --epochs 100 --batch_size 64
```

### **Testing**
```bash
# Run all tests
python run_tests.py

# Run specific test file
pytest models/test_bioactivity_model.py -v

# Run with coverage report
pytest --cov=models --cov=agents
```

### **Execution**
```bash
# Run orchestrator
python run_chemai.py

# Run with specific config
python run_chemai.py --config config.json

# Generate molecules
python -c "from agents.generator_agent import GeneratorAgent; \
  g = GeneratorAgent(); \
  mols = g.generate_molecules(100); \
  print(f'Generated {len(mols)} molecules')"
```

### **Dashboard**
```bash
# Start web dashboard
python web_dashboard.py

# Start API server
python chemai_api_server.py

# Load to Neo4j
python TEAM_SHREYA/load_to_neo4j.py
```

---

## 📋 File Organization Guide

### **By Purpose**

**Main Scripts**
- `run_chemai.py` - Orchestration
- `web_dashboard.py` - Web UI
- `chemai_api_server.py` - REST API

**Training Scripts** (`models/`)
- `train_*.py` - Train each model
- `train_all_models.py` - Train all

**Testing Scripts** (`models/`)
- `test_*.py` - Test each model
- `run_tests.py` - Run all tests
- `comprehensive_model_testing.py` - Full suite

**Agent Scripts** (`agents/`)
- `orchestrator_agent.py` - Coordination
- `generator_agent.py` - Generation
- `optimizer_agent.py` - Optimization
- `predictor_agent.py` - Prediction
- `ranker_agent.py` - Ranking

---

## 🎯 Common Workflows

### **Workflow 1: Train a New Model**
```bash
# 1. Prepare data
# 2. Create training script: models/train_newmodel.py
# 3. Create test script: models/test_newmodel.py
# 4. Run training
python models/train_newmodel.py
# 5. Run tests
pytest models/test_newmodel.py -v
# 6. Integrate with predictor
# Edit: agents/predictor_agent.py
```

### **Workflow 2: Run a Discovery Pipeline**
```bash
# 1. Setup seed molecules
seeds = ['SMILES1', 'SMILES2', ...]

# 2. Create config
config = {
    'num_generations': 5,
    'molecules_per_generation': 100,
    'optimization_iterations': 30
}

# 3. Run orchestrator
python run_chemai.py

# 4. View results on dashboard
# Open: http://localhost:5000

# 5. Query results in Neo4j
# Open: http://localhost:7474/browser/
```

### **Workflow 3: Debug an Agent**
```bash
# 1. Start Python REPL
python

# 2. Import and test agent
from agents.generator_agent import GeneratorAgent
gen = GeneratorAgent()
mols = gen.generate_molecules(10)

# 3. Check outputs
print(len(mols))
print(gen.calculate_validity(mols))

# 4. Debug specific issue
# Add print statements
# Run tests: pytest models/test_*.py -v --pdb
```

---

## 📊 Data Format Reference

### **SMILES String**
- Simplified Molecular Input Line Entry System
- Linear notation for chemical structures
- Example: `CC(C)Cc1ccc(cc1)C(C)C(=O)O` (Ibuprofen)

### **CSV Format for Loading**

**molecules_kg.csv**:
```csv
id,smiles,name,MW,LogP,QED,is_drug_like,source,confidence
mol_001,CC(C)Cc1ccc(cc1)C(C)C(=O)O,Ibuprofen,206.28,3.97,0.83,true,ChEMBL,0.95
```

**bioactivity_edges.csv**:
```csv
source,target,pic50,activity_type,confidence,source
mol_001,prot_001,5.2,IC50,0.92,ChEMBL
```

---

## 🔗 Cross-Reference

| Need | Go To |
|------|-------|
| How to train models | ML_MODELS/README.md |
| How to run agents | AGENTS_FRAMEWORK/README.md |
| How to test code | TESTING_VALIDATION/README.md |
| How to use Neo4j | DATABASE_AND_GRAPHS/README.md |
| How to setup project | SETUP_AND_QUICKSTART/README.md |
| How to understand architecture | ARCHITECTURE_AND_DESIGN/README.md |

---

## 🎓 Learning Resources

### **Python & Chemistry**
- RDKit documentation: https://www.rdkit.org/docs/
- Scikit-learn: https://scikit-learn.org/
- PyTorch: https://pytorch.org/

### **Graph Databases**
- Neo4j documentation: https://neo4j.com/docs/
- Cypher query language: https://neo4j.com/docs/cypher-manual/

### **Cheminformatics**
- ChEMBL database: https://www.ebi.ac.uk/chembl/
- PubChem: https://pubchem.ncbi.nlm.nih.gov/
- SMILES format: https://www.daylight.com/smiles/

---

## 💡 Tips & Tricks

### **Performance Tips**
- Cache molecule descriptors
- Use batch processing for predictions
- Index Neo4j queries
- Profile code with cProfile

### **Debugging Tips**
- Use Python debugger: `import pdb; pdb.set_trace()`
- Add verbose logging
- Use pytest with `--pdb` flag
- Test smaller datasets first

### **Optimization Tips**
- Parallelize expensive operations
- Use GPU for deep learning
- Minimize Neo4j queries
- Cache frequently accessed data

---

**Master Index**: Go back to `../../INDEX.md`
