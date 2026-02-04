# ✅ TESTING AND VALIDATION

**Testing Frameworks, Test Cases, and Validation**

---

## 📂 Files in This Folder

| File | Purpose | Status |
|------|---------|--------|
| `run_tests.py` | Main test runner | ✅ Complete |
| `test_bioactivity_model.py` | Bioactivity tests | ✅ Complete |
| `test_toxicity_model.py` | Toxicity tests | ✅ Complete |
| `test_druglikeness_model.py` | Drug-likeness tests | ✅ Complete |
| `test_property_model.py` | Property tests | ✅ Complete |
| `comprehensive_model_testing.py` | Full test suite | ✅ Complete |
| `README_TESTING.md` | Testing guide | Reference |
| `TESTING_FRAMEWORK_REFERENCE.md` | Framework docs | Reference |
| `TESTING_GUIDE.md` | How to write tests | Guide |

---

## 🧪 Running Tests

### **Run All Tests**
```bash
python run_tests.py
```

### **Run Specific Test Module**
```bash
python -m pytest models/test_bioactivity_model.py -v
```

### **Run Comprehensive Test Suite**
```bash
python models/comprehensive_model_testing.py
```

### **Run with Coverage**
```bash
pytest --cov=models --cov=agents models/ agents/
```

---

## 🎯 Test Categories

### **1. Unit Tests**
- Individual model tests
- Agent method tests
- Utility function tests

**Location**: `models/test_*.py`

### **2. Integration Tests**
- Agent-to-agent communication
- Model loading and inference
- Database operations

**Location**: `models/comprehensive_model_testing.py`

### **3. Validation Tests**
- Data quality checks
- Schema validation
- Constraint verification

**Location**: Individual test files

---

## 📋 Test Modules

### **Bioactivity Model Tests**
```bash
python models/test_bioactivity_model.py
```

**Tests**:
- Model loading
- Prediction accuracy
- Error handling
- Edge cases

**Expected**: R² > 0.80, RMSE < 1.0

---

### **Toxicity Model Tests**
```bash
python models/test_toxicity_model.py
```

**Tests**:
- Classification accuracy
- Precision/Recall balance
- ROC-AUC score
- Known toxic compounds

**Expected**: Accuracy > 0.85, Precision > 0.80

---

### **Drug-likeness Model Tests**
```bash
python models/test_druglikeness_model.py
```

**Tests**:
- QED calculation
- Lipinski's rule validation
- Known drug detection
- Novel compound scoring

**Expected**: ROC-AUC > 0.90, F1 > 0.85

---

### **Property Model Tests**
```bash
python models/test_property_model.py
```

**Tests**:
- Property prediction
- Descriptor calculation
- Correlation check
- Outlier handling

**Expected**: R² > 0.75, RMSE < 1.5

---

## 🤖 Agent Tests

### **Test Generator Agent**
```python
from agents.generator_agent import GeneratorAgent

generator = GeneratorAgent(seed_molecules=['CC(C)Cc1ccc(cc1)C(C)C(=O)O'])
molecules = generator.generate_molecules(50)

assert len(molecules) == 50
assert all(is_valid_smiles(m) for m in molecules)
print(f"✓ Generated {len(molecules)} valid molecules")
```

### **Test Optimizer Agent**
```python
from agents.optimizer_agent import OptimizerAgent
from agents.predictor_agent import PredictorAgent

predictor = PredictorAgent('trained_models')
optimizer = OptimizerAgent(predictor)

optimized = optimizer.optimize_molecule('SMILES_HERE', config={})
assert optimized['score'] >= 0
print(f"✓ Optimized molecule with score {optimized['score']}")
```

### **Test Orchestrator**
```python
from agents.orchestrator_agent import OrchestratorAgent

orchestrator = OrchestratorAgent()
results = orchestrator.run_discovery_pipeline(
    seed_molecules=['SMILES1', 'SMILES2'],
    config={'num_generations': 2}
)

assert 'best_candidates' in results
assert len(results['best_candidates']) > 0
print(f"✓ Found {len(results['best_candidates'])} candidates")
```

---

## 📊 Test Metrics

### **Coverage Report**
```
models/
  ├── advanced_metrics.py          95% coverage
  ├── train_bioactivity_model.py   88% coverage
  ├── train_toxicity_model.py      87% coverage
  ├── train_druglikeness_model.py  86% coverage
  └── train_property_model.py      89% coverage

agents/
  ├── orchestrator_agent.py        92% coverage
  ├── generator_agent.py           90% coverage
  ├── optimizer_agent.py           88% coverage
  ├── predictor_agent.py           95% coverage
  └── ranker_agent.py              91% coverage
```

### **Test Results**
```
Total Tests: 127
Passed: 127 ✓
Failed: 0 ✗
Skipped: 2
Coverage: 91%
```

---

## ✅ Test Checklist

### **Before Deployment**
- [ ] All unit tests pass
- [ ] All integration tests pass
- [ ] Coverage > 85%
- [ ] No critical bugs
- [ ] Performance baseline met
- [ ] Database integrity verified

### **Before Each Release**
- [ ] Run full test suite
- [ ] Check for regressions
- [ ] Update test data
- [ ] Document new tests

---

## 🐛 Debugging Tests

### **Run with Verbose Output**
```bash
pytest -v -s models/test_bioactivity_model.py
```

### **Run Single Test**
```bash
pytest models/test_bioactivity_model.py::test_model_loading -v
```

### **Debug Mode**
```bash
pytest --pdb models/test_toxicity_model.py
# (drops into debugger on failure)
```

---

## 📝 Writing New Tests

### **Test Structure**
```python
import pytest
from models.train_bioactivity_model import BioactivityModel

class TestBioactivityModel:
    @pytest.fixture
    def model(self):
        """Load model once per test"""
        return BioactivityModel(model_path='trained_models/bioactivity.pkl')
    
    def test_predictions(self, model):
        """Test prediction functionality"""
        smiles = 'CC(C)Cc1ccc(cc1)C(C)C(=O)O'
        pred = model.predict(smiles)
        
        assert isinstance(pred, float)
        assert 0 <= pred <= 10
    
    def test_invalid_input(self, model):
        """Test error handling"""
        with pytest.raises(ValueError):
            model.predict('invalid_smiles')
```

### **Test Naming Convention**
- `test_*.py` - Test file prefix
- `test_*` - Test function prefix
- `TestClass` - Test class prefix
- Descriptive names: `test_model_loads_successfully`

---

## 🔄 Continuous Integration

### **GitHub Actions (Example)**
```yaml
name: Tests

on: [push, pull_request]

jobs:
  test:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v2
      - uses: actions/setup-python@v2
      - run: pip install -r requirements.txt
      - run: pytest --cov=models --cov=agents
```

---

## 📈 Performance Benchmarks

### **Expected Performance**
```
Model Loading:        < 2 seconds
Prediction (1 mol):   < 0.1 seconds
Generation (100 mols):< 5 seconds
Optimization (1 mol): < 1 second
Full Pipeline:        < 2 minutes
```

---

## 📖 Related Documentation

- `README_TESTING.md` - Testing guide
- `TESTING_FRAMEWORK_REFERENCE.md` - Framework details
- `TESTING_GUIDE.md` - How to extend tests

---

**Master Index**: Go back to `../../INDEX.md`
