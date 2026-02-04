# ChemAI Models - Metrics Summary Table

## Quick Reference Table

| Model | Type | Accuracy | F1-Score | R² | Target(s) | Features |
|-------|------|----------|----------|-----|-----------|----------|
| **Bioactivity** | Regression | N/A | N/A | ✓ | pIC50 (Bioactivity) | Morgan FP (2048) |
| **Toxicity** | Classification | ✓ | ✓ | ✓ | Structural Alerts | Morgan FP (1024) |
| **Property** | Multi-Regression | N/A | N/A | ✓×7 | 7 Molecular Properties | Morgan FP (1024) |
| **QED** | Regression | N/A | N/A | ✓ | QED Score (0-1) | Morgan FP (1024) |
| **Drug-likeness** | Classification | ✓ | ✓ | ✓ | Lipinski's Rule of 5 | Morgan FP (1024) |

---

## Detailed Metrics for Each Model

### 1. BIOACTIVITY MODEL
```
Task Type:          Regression
Target Variable:    pIC50 (Bioactivity)
Output:             Continuous value (0-15)
Accuracy:           N/A
F1-Score:           N/A
R²:                 Primary metric
Additional:         RMSE, MAE
Test Samples:       ~5,000
```

### 2. TOXICITY MODEL
```
Task Type:          Binary Classification
Target Variable:    Structural Alerts (0=Safe, 1=Toxic)
Output:             Binary (0 or 1)
Accuracy:           Available (typical: 0.80-0.92)
F1-Score:           Available (typical: 0.75-0.88)
R²:                 Available (on probabilities)
Additional:         Confusion Matrix, ROC-AUC
Test Samples:       ~5,000
```

### 3. PROPERTY MODEL
```
Task Type:          Multi-target Regression
Target Variables:   7 properties
  - Molecular Weight (MW)
  - Lipophilicity (LogP)
  - H-Bond Acceptors (HBA)
  - H-Bond Donors (HBD)
  - Polar Surface Area (PSA)
  - Rotatable Bonds (RTB)
  - QED Score

Accuracy:           N/A (one per property)
F1-Score:           N/A (one per property)
R²:                 Available (one per property)
Additional:         RMSE, MAE per property
Test Samples:       ~5,000
```

### 4. QED MODEL
```
Task Type:          Regression
Target Variable:    QED Score (0-1 range)
Output:             Continuous value
Accuracy:           N/A
F1-Score:           N/A
R²:                 Primary metric
Additional:         RMSE, MAE
Test Samples:       ~5,000
```

### 5. DRUG-LIKENESS MODEL
```
Task Type:          Binary Classification
Target Variable:    Lipinski Compliance (0=Violates, 1=Complies)
Output:             Binary (0 or 1)
Accuracy:           Available (typical: 0.85-0.95)
F1-Score:           Available (typical: 0.80-0.92)
R²:                 Available (on probabilities)
Additional:         Confusion Matrix
Test Samples:       ~5,000
```

---

## Metric Calculations

### For Regression Models (Bioactivity, Property, QED)

**R² Score:**
```
R² = 1 - (SS_residual / SS_total)
where:
  SS_residual = Σ(y_actual - y_predicted)²
  SS_total = Σ(y_actual - y_mean)²

Range: 0 to 1 (higher is better)
Interpretation:
  R² = 1.0  → Perfect fit
  R² > 0.7  → Good fit
  R² > 0.5  → Acceptable fit
  R² < 0.3  → Poor fit
```

**RMSE:**
```
RMSE = √(Σ(y_actual - y_predicted)² / n)

Same units as target variable
Lower is better
```

**MAE:**
```
MAE = Σ|y_actual - y_predicted| / n

Same units as target variable
Lower is better
```

---

### For Classification Models (Toxicity, Drug-likeness)

**Accuracy:**
```
Accuracy = (TP + TN) / (TP + TN + FP + FN)

Range: 0 to 1 (higher is better)
Caution: Can be misleading with imbalanced data
```

**F1-Score:**
```
F1 = 2 × (Precision × Recall) / (Precision + Recall)

where:
  Precision = TP / (TP + FP)  [accuracy of positive predictions]
  Recall = TP / (TP + FN)     [sensitivity/true positive rate]

Range: 0 to 1 (higher is better)
Better than Accuracy when data is imbalanced
```

**Confusion Matrix Terms:**
```
TP (True Positive):   Model predicted 1, actual is 1 ✓
TN (True Negative):   Model predicted 0, actual is 0 ✓
FP (False Positive):  Model predicted 1, actual is 0 ✗
FN (False Negative):  Model predicted 0, actual is 1 ✗
```

---

## How to Read the Comprehensive Test Output

When you run `python models/comprehensive_model_testing.py`, you'll see:

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

Confusion Matrix:
  TN: 2145, FP: 312
  FN: 287, TP: 1945
```

---

## Performance Benchmarks

### Typical Expected Results

| Model | Metric | Expected Range | Notes |
|-------|--------|-----------------|-------|
| Bioactivity | R² | 0.65 - 0.80 | Depends on data quality |
| Toxicity | Accuracy | 0.80 - 0.92 | Binary classification task |
| Toxicity | F1-Score | 0.75 - 0.88 | Balanced performance |
| Property | R² (avg) | 0.60 - 0.78 | Varies per property |
| QED | R² | 0.70 - 0.85 | Easier target |
| Drug-likeness | Accuracy | 0.85 - 0.95 | Well-defined rules |
| Drug-likeness | F1-Score | 0.80 - 0.92 | Good discrimination |

---

## Key Takeaways

1. **Regression models** focus on **R²** (with RMSE/MAE for magnitude)
2. **Classification models** focus on **Accuracy & F1-Score**
3. **F1-Score > Accuracy** when classes are imbalanced
4. **All metrics are reported** for comparison and validation
5. **Test data** comes from ChEMBL database (~5,000 samples each)
6. **Higher values** = Better performance (for all metrics shown)

---

## Testing Commands

```bash
# Option 1: Quick metrics overview
python METRICS_GUIDE.py

# Option 2: Comprehensive testing of all models
python models/comprehensive_model_testing.py

# Option 3: Individual model testing
python models/test_bioactivity_model.py
python models/test_toxicity_model.py
python models/test_property_model.py
python models/test_druglikeness_model.py
```

---

**Document Version:** 1.0  
**Last Updated:** 2026-01-31  
**ChemAI Project**
