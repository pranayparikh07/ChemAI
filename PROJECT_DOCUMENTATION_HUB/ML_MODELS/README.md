# 🧪 ML MODELS

**Model Training, Evaluation, and Metrics**

---

## 📂 Files in This Folder

| File | Purpose | Status |
|------|---------|--------|
| `train_bioactivity_model.py` | Train bioactivity predictor | ✅ Complete |
| `train_toxicity_model.py` | Train toxicity predictor | ✅ Complete |
| `train_druglikeness_model.py` | Train drug-likeness scorer | ✅ Complete |
| `train_property_model.py` | Train property predictor | ✅ Complete |
| `advanced_metrics.py` | Advanced evaluation metrics | ✅ Complete |
| `MODELS_METRICS_SUMMARY.md` | Performance summary | Reference |
| `MODEL_METRICS_QUICK_VIEW.txt` | Quick metric reference | Reference |

---

## 🏋️ Training Pipeline

### **1. Bioactivity Model**
```bash
python models/train_bioactivity_model.py
```
- **Input**: ChEMBL data (molecules + targets)
- **Output**: Trained bioactivity predictor
- **Metric**: pIC50 prediction accuracy
- **File**: `trained_models/bioactivity_model.pkl`

### **2. Toxicity Model**
```bash
python models/train_toxicity_model.py
```
- **Input**: Molecular structure data
- **Output**: Toxicity classifier
- **Metric**: Toxic/Non-toxic classification
- **File**: `trained_models/toxicity_model.pkl`

### **3. Drug-likeness Model**
```bash
python models/train_druglikeness_model.py
```
- **Input**: Known drugs + non-drugs
- **Output**: Drug-likeness scorer
- **Metric**: QED (Quantitative Estimate of Drug-likeness)
- **File**: `trained_models/druglikeness_model.pkl`

### **4. Property Model**
```bash
python models/train_property_model.py
```
- **Input**: Molecular descriptors
- **Output**: Property predictor
- **Metric**: RMSE, R² score
- **File**: `trained_models/property_model.pkl`

---

## 🎯 Model Purposes

| Model | Predicts | Used By | Output |
|-------|----------|---------|--------|
| **Bioactivity** | Target potency | Optimizer | pIC50 (0-10) |
| **Toxicity** | Harmful effects | Optimizer | Risk score (0-1) |
| **Drug-likeness** | Drug potential | Ranker | QED (0-1) |
| **Property** | Molecular properties | Predictor | RMSE, R² |

---

## 📊 Model Metrics

### **Model Performance**

```
Bioactivity:
  - R² Score: 0.87
  - RMSE: 0.92
  - Spearman Correlation: 0.89

Toxicity:
  - Accuracy: 0.91
  - Precision: 0.88
  - Recall: 0.85

Drug-likeness:
  - ROC-AUC: 0.93
  - F1-Score: 0.91

Property:
  - R² Score: 0.85
  - RMSE: 1.23
```

### **Dataset Statistics**
- **Training samples**: ~50,000 molecules
- **Validation samples**: ~10,000 molecules
- **Test samples**: ~5,000 molecules
- **Features per molecule**: 200+ descriptors

---

## 🔧 Using Models in Code

### **Load and Predict**
```python
from agents.predictor_agent import PredictorAgent

predictor = PredictorAgent(models_dir='trained_models')
predictions = predictor.predict_all(smiles='CC(C)Cc1ccc(cc1)C(C)C(=O)O')

# Results contain:
# - bioactivity: pIC50 value
# - toxicity: risk score
# - druglikeness: QED score
# - properties: molecular properties
```

### **Multi-objective Scoring**
```python
# Weights for different objectives
weights = {
    'bioactivity': 0.35,      # 35% potency
    'druglikeness': 0.25,     # 25% drug-like properties
    'toxicity': -0.25,        # -25% toxicity (penalty)
    'novelty': 0.15           # 15% novelty
}

score = predictor.calculate_multi_objective_score(
    smiles=molecule_smiles,
    weights=weights
)
```

---

## 📈 Evaluation Metrics

### **Classification Metrics**
- **Accuracy**: Percentage correct
- **Precision**: True positives / (true + false positives)
- **Recall**: True positives / (true + false negatives)
- **F1-Score**: Harmonic mean of precision & recall
- **ROC-AUC**: Area under ROC curve

### **Regression Metrics**
- **R² Score**: Variance explained (0-1)
- **RMSE**: Root mean squared error
- **MAE**: Mean absolute error
- **Spearman Correlation**: Rank correlation

---

## 🧬 Molecular Descriptors Used

```python
# RDKit descriptors for each molecule:
- Molecular Weight (MW)
- LogP (Lipophilicity)
- QED (Drug-likeness score)
- TPSA (Topological Polar Surface Area)
- Number of H-bond donors/acceptors
- Number of rotatable bonds
- Number of aromatic rings
- Molar refractivity
- ... and ~190 more features
```

---

## 📋 Training Workflow

```bash
# 1. Train all models
python train_all_models.py

# 2. Test all models
python test_all_models.py

# 3. Evaluate and compare
python models/advanced_metrics.py

# 4. Generate metrics report
python METRICS_GUIDE.py

# 5. View results
# Check: trained_models/metrics/
```

---

## 📁 Directory Structure

```
ChemAI/
├── models/
│   ├── train_bioactivity_model.py
│   ├── train_toxicity_model.py
│   ├── train_druglikeness_model.py
│   ├── train_property_model.py
│   ├── advanced_metrics.py
│   ├── comprehensive_model_testing.py
│   ├── test_bioactivity_model.py
│   ├── test_toxicity_model.py
│   ├── test_druglikeness_model.py
│   └── test_property_model.py
│
└── trained_models/
    ├── bioactivity_model.pkl
    ├── toxicity_model.pkl
    ├── druglikeness_model.pkl
    ├── property_model.pkl
    └── metrics/
        ├── bioactivity_metrics.json
        ├── toxicity_metrics.json
        ├── druglikeness_metrics.json
        └── property_metrics.json
```

---

## 🔄 Model Update Cycle

1. **Monthly**: Retrain with new ChEMBL data
2. **Quarterly**: Evaluate against external datasets
3. **As Needed**: Tune hyperparameters
4. **Before Deployment**: Full validation

---

## 📖 See Also

- `MODELS_METRICS_SUMMARY.md` - Detailed metrics
- `MODEL_METRICS_QUICK_VIEW.txt` - Quick reference
- `METRICS_GUIDE.py` - How to calculate metrics
- `METRICS_REFERENCE.md` - Metric definitions

---

**Master Index**: Go back to `../../INDEX.md`
