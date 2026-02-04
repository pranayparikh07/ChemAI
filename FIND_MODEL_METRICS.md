# Where to Find Model Metrics - Accuracy, F1, and R²

## Files in This Directory

### 📊 Main Reference Documents
1. **[METRICS_REFERENCE.md](METRICS_REFERENCE.md)** ⭐ START HERE
   - Complete definitions of Accuracy, F1-Score, and R²
   - Which metrics apply to each model
   - Expected performance ranges
   - Metric definitions and formulas

2. **[MODELS_METRICS_SUMMARY.md](MODELS_METRICS_SUMMARY.md)** ⭐ DETAILED SUMMARY
   - Quick reference table of all models
   - Detailed metrics for each model
   - How to read the test output
   - Performance benchmarks

### 🧪 Testing Scripts

#### Comprehensive Testing (Recommended)
- **File:** `models/comprehensive_model_testing.py` 
- **Run:** `python models/comprehensive_model_testing.py`
- **Output:** Shows Accuracy, F1, and R² for all models at once
- **Takes:** 5-15 minutes depending on database access

#### Individual Model Tests
1. **Bioactivity Model**
   - File: `models/test_bioactivity_model.py`
   - Metrics: R², RMSE, MAE
   - Run: `python models/test_bioactivity_model.py`

2. **Toxicity Model**
   - File: `models/test_toxicity_model.py`
   - Metrics: Accuracy, F1-Score, ROC-AUC
   - Run: `python models/test_toxicity_model.py`

3. **Property Model**
   - File: `models/test_property_model.py`
   - Metrics: R², RMSE, MAE (per property)
   - Run: `python models/test_property_model.py`

4. **Drug-likeness Model**
   - File: `models/test_druglikeness_model.py`
   - Metrics: Accuracy, F1-Score, R² (QED + Drug-likeness)
   - Run: `python models/test_druglikeness_model.py`

#### Quick Reference
- **File:** `METRICS_GUIDE.py`
- **Run:** `python METRICS_GUIDE.py`
- **Output:** Shows model specifications and expected ranges

---

## Quick Answer: Which Metric Goes With Which Model?

### ✓ Accuracy Available For:
- ☑ Toxicity Model (Binary Classification)
- ☑ Drug-likeness Model (Binary Classification)
- ✗ Bioactivity Model (Regression)
- ✗ Property Model (Regression)
- ✗ QED Model (Regression)

### ✓ F1-Score Available For:
- ☑ Toxicity Model (Binary Classification)
- ☑ Drug-likeness Model (Binary Classification)
- ✗ Bioactivity Model (Regression)
- ✗ Property Model (Regression)
- ✗ QED Model (Regression)

### ✓ R² Available For:
- ☑ Bioactivity Model (Regression)
- ☑ Toxicity Model (on probabilities)
- ☑ Property Model (per property, 7 total)
- ☑ QED Model (Regression)
- ☑ Drug-likeness Model (on probabilities)

---

## Code Snippets to Calculate These Metrics

### In Python - Calculate Accuracy
```python
from sklearn.metrics import accuracy_score
accuracy = accuracy_score(y_true, y_pred)
print(f"Accuracy: {accuracy:.4f}")
```

### In Python - Calculate F1-Score
```python
from sklearn.metrics import f1_score
f1 = f1_score(y_true, y_pred)
print(f"F1-Score: {f1:.4f}")
```

### In Python - Calculate R²
```python
from sklearn.metrics import r2_score
r2 = r2_score(y_true, y_pred)
print(f"R²: {r2:.4f}")
```

### All Three Together
```python
from sklearn.metrics import accuracy_score, f1_score, r2_score

# For classification
accuracy = accuracy_score(y_true, y_pred)
f1 = f1_score(y_true, y_pred)
r2 = r2_score(y_true, y_pred_proba)  # usually on probabilities

print(f"Accuracy: {accuracy:.4f}")
print(f"F1-Score: {f1:.4f}")
print(f"R²:       {r2:.4f}")
```

---

## Example Expected Output

```
================================================================================
BIOACTIVITY MODEL - REGRESSION METRICS
================================================================================
Model: Bioactivity Prediction (pIC50)
Test Samples: 4,523

Metrics:
  R²:       0.7234
  RMSE:     0.6543
  MAE:      0.4892
  Accuracy: N/A (Regression Task)
  F1-Score: N/A (Regression Task)

================================================================================
TOXICITY MODEL - CLASSIFICATION METRICS
================================================================================
Model: Toxicity Prediction (Structural Alerts)
Test Samples: 4,789

Metrics:
  Accuracy:     0.8734
  F1-Score:     0.8421
  R² (Proba):   0.6892

================================================================================
PROPERTY MODEL - MULTI-TARGET REGRESSION METRICS
================================================================================
Model: Molecular Property Prediction
Test Samples: 4,456
Properties: 7

mw_freebase
  R²:       0.7821
  RMSE:     25.3456
  MAE:      15.2341

alogp
  R²:       0.6542
  RMSE:     0.7234
  MAE:      0.5123

[... more properties ...]

================================================================================
DRUGLIKENESS MODEL - CLASSIFICATION METRICS
================================================================================
Model 1: QED Prediction (Regression)
  R²:       0.7456
  RMSE:     0.0834
  MAE:      0.0612

Model 2: Drug-likeness Classification (Lipinski's Rule)
  Accuracy:     0.8934
  F1-Score:     0.8756
  R² (Proba):   0.7123
```

---

## Interpretation Guide

### R² Interpretation
- **R² = 1.0** → Perfect predictions
- **R² > 0.7** → Good model, explains 70%+ of variance
- **R² = 0.5-0.7** → Moderate model
- **R² < 0.3** → Poor model, not recommended

### Accuracy Interpretation
- **Accuracy > 0.9** → Excellent
- **Accuracy 0.8-0.9** → Very Good
- **Accuracy 0.7-0.8** → Good
- **Accuracy < 0.7** → Needs improvement

### F1-Score Interpretation
- **F1 > 0.9** → Excellent balance
- **F1 0.8-0.9** → Very Good
- **F1 0.7-0.8** → Good
- **F1 < 0.7** → Needs improvement

---

## Summary Table

| Want to Know... | Check File | Run Command |
|-----------------|-----------|------------|
| Which metrics for each model | METRICS_REFERENCE.md | - |
| Complete metrics summary | MODELS_METRICS_SUMMARY.md | - |
| Actual metric values | models/comprehensive_model_testing.py | `python models/comprehensive_model_testing.py` |
| Bioactivity R², RMSE, MAE | models/test_bioactivity_model.py | `python models/test_bioactivity_model.py` |
| Toxicity Accuracy & F1 | models/test_toxicity_model.py | `python models/test_toxicity_model.py` |
| Property R² per attribute | models/test_property_model.py | `python models/test_property_model.py` |
| Drug-likeness Accuracy & F1 | models/test_druglikeness_model.py | `python models/test_druglikeness_model.py` |
| Quick reference | METRICS_GUIDE.py | `python METRICS_GUIDE.py` |

---

## Files Structure
```
ChemAI/
├── METRICS_REFERENCE.md          ← Definitions and formulas
├── MODELS_METRICS_SUMMARY.md     ← Summary table and details
├── FIND_MODEL_METRICS.md         ← This file
├── METRICS_GUIDE.py              ← Quick reference script
└── models/
    ├── comprehensive_model_testing.py  ← Run all tests at once
    ├── test_bioactivity_model.py       ← Bioactivity testing
    ├── test_toxicity_model.py          ← Toxicity testing
    ├── test_property_model.py          ← Property testing
    └── test_druglikeness_model.py      ← Drug-likeness testing
```

---

**Last Updated:** 2026-01-31  
**For Questions:** Refer to METRICS_REFERENCE.md for detailed definitions
