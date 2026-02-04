# ChemAI Model Testing Suite

Complete automated testing framework for all trained ChemAI models with HTML report generation.

## Overview

The testing suite provides comprehensive evaluation of all trained models:
- **Property Model** - Predicts 7 molecular properties (MW, LogP, HBA, HBD, PSA, RTB, QED)
- **Druglikeness Models** - QED regressor & Lipinski classifier
- **Bioactivity Model** - Predicts pIC50 bioactivity scores

## Files

### Main Test Orchestrator
- **`test_all_models.py`** - Runs all tests and generates HTML report
  - Executes all model tests in sequence
  - Aggregates results from individual test modules
  - Generates beautiful HTML report with detailed metrics
  - Located in root directory: `d:\ChemAI\test_all_models.py`

### Individual Model Tests
Located in `models/` directory:

1. **`test_property_model.py`**
   - Tests multi-output property prediction
   - Metrics: RMSE, MAE, R² for each property
   - Sample properties: MW, LogP, HBA, HBD, PSA, RTB, QED

2. **`test_druglikeness_model.py`**
   - Tests QED score regression
   - Tests Lipinski drug-likeness classification
   - Metrics: RMSE/MAE/R² for QED, Accuracy/Confusion matrix for classification

3. **`test_bioactivity_model.py`**
   - Tests pIC50 bioactivity prediction
   - Metrics: RMSE, MAE, R², Residuals statistics

## Running Tests

### Option 1: Run All Tests (Recommended)
```bash
cd d:\ChemAI
python test_all_models.py
```

This will:
1. Execute all three model tests sequentially
2. Collect comprehensive metrics from each model
3. Generate `model_test_report.html` with visualized results
4. Display summary in console

### Option 2: Run Individual Tests
```bash
# Test property model only
cd d:\ChemAI\models
python test_property_model.py

# Test druglikeness models only
python test_druglikeness_model.py

# Test bioactivity model only
python test_bioactivity_model.py
```

## Output

### Console Output
- Real-time test execution status
- Metrics for each model displayed as they complete
- Summary statistics for each model

### HTML Report
Generated file: `model_test_report.html`

Features:
- **Professional Dashboard** - Purple gradient theme with modern UI
- **Model-by-Model Results** - Detailed sections for each model
- **Comprehensive Metrics** - All performance indicators displayed
- **Visual Cards** - Color-coded metric cards for quick assessment
- **Detailed Tables** - Full breakdown of results per property/metric
- **Status Badges** - Pass/Fail indicators for each model
- **Summary Section** - Overall testing statistics

#### Report Sections:
1. **Header** - Title, timestamp, general info
2. **Property Model** - All 7 property predictions with R², RMSE, MAE
3. **Druglikeness Models** - QED regressor metrics + classification accuracy
4. **Bioactivity Model** - pIC50 prediction with residuals analysis
5. **Summary** - Overall test results and coverage
6. **Footer** - Report metadata

## Test Data

Tests use real ChEMBL database molecules:
- **Property Model**: 5,000 molecules with computed properties
- **Druglikeness Models**: 5,000 molecules with QED/Lipinski properties
- **Bioactivity Model**: 5,000 IC50 bioactivity measurements

## Metrics Explained

### Regression Metrics
- **RMSE** (Root Mean Square Error) - Lower is better
- **MAE** (Mean Absolute Error) - Average absolute deviation
- **R²** (Coefficient of Determination) - 0-1 scale, higher is better

### Classification Metrics
- **Accuracy** - Percentage of correct predictions
- **Confusion Matrix** - TP, TN, FP, FN breakdown
- **Classification Report** - Precision, recall, F1-score per class

### Residuals Analysis
- **Mean** - Should be close to 0
- **Std Dev** - Spread of errors
- **Min/Max** - Error range

## Troubleshooting

### Models not found
```
✗ Model file not found at trained_models/property_model.joblib
```
→ Run `train_all_models.py` first to train models

### Database connection error
```
✗ Error loading test data
```
→ Ensure ChEMBL database exists at `chembl_36/chembl_36_sqlite/chembl_36.db`
→ Run `scan_and_clean.py` to extract database if needed

### Module import errors
```
ModuleNotFoundError: No module named 'rdkit'
```
→ Install requirements: `pip install -r requirements.txt`

## Generated Files

After running tests:
```
d:\ChemAI\
├── model_test_report.html          ← Main HTML report (OPEN IN BROWSER)
├── test_all_models.py              ← Test orchestrator
└── models/
    ├── test_property_model.py       ← Individual test modules
    ├── test_druglikeness_model.py
    └── test_bioactivity_model.py
```

## Viewing the Report

1. Run `python test_all_models.py`
2. Open `model_test_report.html` in any web browser
3. Review detailed metrics and visualizations

## Performance Benchmarks

Expected Results (approximate):
- **Property Model R²**: 0.6-0.8 (varies by property)
- **QED Regressor R²**: 0.7-0.85
- **Drug-likeness Accuracy**: 85-95%
- **Bioactivity R²**: 0.5-0.75

## Testing Architecture

```
test_all_models.py (Main Orchestrator)
├── test_property_model.py
│   ├── Load model + scaler
│   ├── Get test data from ChEMBL
│   ├── Compute Morgan fingerprints
│   └── Evaluate on 7 properties
│
├── test_druglikeness_model.py
│   ├── Load QED + classifier models
│   ├── Get test data from ChEMBL
│   ├── Compute Morgan fingerprints
│   ├── Evaluate QED regressor
│   └── Evaluate Lipinski classifier
│
└── test_bioactivity_model.py
    ├── Load bioactivity model
    ├── Get bioactivity test data
    ├── Compute Morgan fingerprints
    └── Evaluate pIC50 prediction

Generate HTML Report
└── model_test_report.html (Beautiful Dashboard)
```

## Requirements

- Python 3.8+
- Dependencies in `requirements.txt`:
  - pandas, numpy
  - rdkit
  - scikit-learn
  - joblib
  - sqlite3

## Notes

- Tests run sequentially to avoid memory issues
- Each test loads ~5,000 molecules from ChEMBL
- Total execution time: 2-5 minutes (depends on system)
- Report is self-contained HTML (can be shared/emailed)
- All metrics calculated on unseen test data
