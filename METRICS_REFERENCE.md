# ChemAI Model Metrics - Accuracy, F1, and R²

## Quick Summary

This document outlines which metrics are calculated for each model in ChemAI.

---

## Model Metrics Overview

### 1. **Bioactivity Model** (Regression)
**Target:** pIC50 (Bioactivity Prediction)
- **Accuracy:** N/A (Regression task)
- **F1-Score:** N/A (Regression task)
- **R²:** ✓ Available (Primary metric)
- **Other metrics:** RMSE, MAE

### 2. **Toxicity Model** (Binary Classification)
**Target:** Structural Alert Detection (Toxic/Safe)
- **Accuracy:** ✓ Available
- **F1-Score:** ✓ Available
- **R²:** ✓ Available (on prediction probabilities)
- **Other metrics:** Confusion Matrix, Precision, Recall, ROC-AUC

### 3. **Property Model** (Multi-target Regression)
**Targets:** 7 Molecular Properties
1. MW (Molecular Weight)
2. LogP (Lipophilicity)
3. HBA (Hydrogen Bond Acceptors)
4. HBD (Hydrogen Bond Donors)
5. PSA (Polar Surface Area)
6. RTB (Rotatable Bonds)
7. QED (Drug-likeness Score)

- **Accuracy:** N/A (Regression tasks)
- **F1-Score:** N/A (Regression tasks)
- **R²:** ✓ Available (per property)
- **Other metrics:** RMSE, MAE (per property)

### 4. **QED Model** (Regression)
**Target:** QED Score (0-1 scale)
- **Accuracy:** N/A (Regression task)
- **F1-Score:** N/A (Regression task)
- **R²:** ✓ Available
- **Other metrics:** RMSE, MAE

### 5. **Drug-likeness Model** (Binary Classification)
**Target:** Lipinski's Rule of 5 Compliance
- **Accuracy:** ✓ Available
- **F1-Score:** ✓ Available
- **R²:** ✓ Available (on prediction probabilities)
- **Other metrics:** Confusion Matrix, Precision, Recall

---

## How to Extract Metrics

### Option 1: Quick Reference
```bash
python METRICS_GUIDE.py
```
Shows all model specifications and expected metric ranges.

### Option 2: Comprehensive Testing
```bash
python models/comprehensive_model_testing.py
```
Runs all models on test data and calculates:
- **Accuracy** (for classification models)
- **F1-Score** (for classification models)
- **R²** (for all models)
- Additional metrics (RMSE, MAE, confusion matrices)

### Option 3: Individual Model Testing
```bash
# Bioactivity
python models/test_bioactivity_model.py

# Toxicity
python models/test_toxicity_model.py

# Property
python models/test_property_model.py

# Drug-likeness
python models/test_druglikeness_model.py
```

---

## Metric Definitions

### **Accuracy**
- **Definition:** (TP + TN) / Total
- **Range:** 0 to 1 (higher is better)
- **When to use:** Classification tasks
- **Caution:** Can be misleading with imbalanced data

### **F1-Score**
- **Definition:** 2 × (Precision × Recall) / (Precision + Recall)
- **Range:** 0 to 1 (higher is better)
- **When to use:** Classification tasks, especially with imbalanced data
- **Better than:** Accuracy alone because it considers both false positives and false negatives

### **R² (Coefficient of Determination)**
- **Definition:** 1 - (SS_res / SS_tot) where SS_res = Σ(y_true - y_pred)², SS_tot = Σ(y_true - y_mean)²
- **Range:** 0 to 1 (higher is better)
- **When to use:** Regression tasks (and can be used on classification probabilities)
- **Interpretation:** 
  - R² = 1: Perfect prediction
  - R² = 0.7+: Good model
  - R² = 0.5-0.7: Moderate model
  - R² < 0.3: Poor model

### **RMSE (Root Mean Squared Error)**
- **Definition:** √(Σ(y_true - y_pred)² / n)
- **Units:** Same as target variable
- **Range:** 0 to ∞ (lower is better)
- **When to use:** Regression tasks
- **Properties:** Penalizes larger errors more heavily

### **MAE (Mean Absolute Error)**
- **Definition:** Σ|y_true - y_pred| / n
- **Units:** Same as target variable
- **Range:** 0 to ∞ (lower is better)
- **When to use:** Regression tasks
- **Properties:** Equal penalty for all errors

---

## Classification Model Confusion Matrix

```
                 Predicted Negative    Predicted Positive
Actual Negative         TN                    FP
Actual Positive         FN                    TP
```

### Derived Metrics:
- **Sensitivity (Recall)** = TP / (TP + FN) — Ability to find all positives
- **Specificity** = TN / (TN + FP) — Ability to find all negatives
- **Precision** = TP / (TP + FP) — Accuracy of positive predictions

---

## Expected Performance Ranges

### Regression Models (Bioactivity, Property, QED)
- **R²:** 0.65 - 0.85 (depends on target complexity)
- **RMSE:** Model-dependent
- **MAE:** Model-dependent

### Classification Models (Toxicity, Drug-likeness)
- **Accuracy:** 0.75 - 0.95
- **F1-Score:** 0.70 - 0.90
- **R² (Proba):** 0.50 - 0.85

---

## Interpreting Results

### If Accuracy is High but F1 is Low
- **Problem:** Class imbalance or model biased toward majority class
- **Solution:** Check F1-score, use stratified splits, adjust class weights

### If R² is Low for Regression
- **Problem:** Poor model fit or difficult prediction task
- **Solution:** Improve features, try different algorithms, check for outliers

### If Metrics are Different Between Train and Test
- **Problem:** Potential overfitting
- **Solution:** Add regularization, more training data, simpler model

---

## File Locations

- **Trained Models:** `trained_models/`
- **Test Scripts:** `models/test_*.py`
- **Comprehensive Test:** `models/comprehensive_model_testing.py`
- **This Guide:** `METRICS_GUIDE.py`

---

## Notes

1. **Regression models** (Bioactivity, Property, QED): Use **R²** as primary metric
2. **Classification models** (Toxicity, Drug-likeness): Use **Accuracy** and **F1-Score**
3. **F1-Score preferred over Accuracy** when data is imbalanced
4. **All models** report R² for comparison purposes
5. Test data is loaded from ChEMBL database when running tests

---

Last Updated: 2026-01-31
