# 🧬 ChemAI Project Documentation Hub

**Master Index of All Project Documentation and Resources**

---

## 📚 Folder Structure Overview

```
PROJECT_DOCUMENTATION_HUB/
├── SETUP_AND_QUICKSTART/           → Getting started and installation
├── TEAM_ROLES_AND_ASSIGNMENTS/     → Team responsibilities and deliverables
├── ARCHITECTURE_AND_DESIGN/        → System design and technical roadmap
├── ML_MODELS/                      → Model training and evaluation
├── AGENTS_FRAMEWORK/               → Agent orchestration and design
├── DATABASE_AND_GRAPHS/            → Neo4j and knowledge graph docs
├── TESTING_VALIDATION/             → Testing frameworks and guides
├── DEPLOYMENT_AND_DASHBOARDS/      → Dashboards and deployment
└── REFERENCE_MATERIALS/            → Knowledge base and references
```

---

## 🚀 Quick Navigation

### **1️⃣ SETUP_AND_QUICKSTART** 
*Start here if you're new to the project*

- **_START_HERE_FIRST.md** - Initial project overview
- **START_HERE.md** - Quick installation guide
- **QUICKSTART_WEB_DASHBOARD.txt** - How to run the web interface
- **README_CHEMAI.md** - Project description
- **requirements.txt** - Python dependencies (linked)

**Purpose**: Get the environment running in minutes

---

### **2️⃣ TEAM_ROLES_AND_ASSIGNMENTS**
*Who does what and deliverables*

- **README_TEAM_ORGANIZATION.txt** - Team structure
- **TEAM_STATUS_AND_REMAINING_WORK.txt** - Current progress and TODOs
- **SHREYA_KNOWLEDGE_GRAPH_COMPLETION.md** - Shreya's focus area
- **Team Subfolders**: TEAM_SHREYA/, TEAM_VISHWA/, TEAM_PRANAY/

**Team Assignments**:
- 🔵 **Shreya**: Agent Orchestration & Autonomous Experiment Loop
- 🟢 **Vishwa**: Reinforcement Learning for Molecular Optimization
- 🟡 **Pranay**: Graph Neural Network-based Molecular Generator

---

### **3️⃣ ARCHITECTURE_AND_DESIGN**
*Technical design and system architecture*

- **TECHNICAL_ROADMAP.md** - Long-term technical plan
- **PROJECT_PROGRESS_ANALYSIS.md** - Implementation status
- **IMPLEMENTATION_SUMMARY.txt** - What's been built
- **VISUAL_FLOW_DIAGRAM.txt** - System flow visualization
- **COMPREHENSIVE_MONTH1_REPORT.txt** - Month 1 progress

**What's Inside**: 
- System design patterns
- Data flow architecture
- Integration points between components

---

### **4️⃣ ML_MODELS**
*Machine learning model documentation*

- **train_bioactivity_model.py** - Bioactivity model training
- **train_toxicity_model.py** - Toxicity prediction model
- **train_druglikeness_model.py** - Drug-likeness scoring
- **train_property_model.py** - Molecular property prediction
- **advanced_metrics.py** - Evaluation metrics
- **MODELS_METRICS_SUMMARY.md** - Model performance metrics
- **MODEL_METRICS_QUICK_VIEW.txt** - Quick reference
- **FIND_MODEL_METRICS.md** - How to extract metrics

**Locations**: 
- Source: `d:\ChemAI\models/`
- Trained: `d:\ChemAI\trained_models/`

---

### **5️⃣ AGENTS_FRAMEWORK**
*Agent orchestration and coordination*

- **orchestrator_agent.py** - Main coordination logic
- **generator_agent.py** - Molecule generation (Pranay's work)
- **optimizer_agent.py** - RL optimization (Vishwa's work)
- **predictor_agent.py** - Property prediction
- **ranker_agent.py** - Molecule ranking

**Location**: `d:\ChemAI\agents/`

**Key Functions**:
- `run_discovery_pipeline()` - Main autonomous loop
- `generate_molecules()` - Generation step
- `optimize_molecule()` - Optimization step
- `rank_candidates()` - Ranking step

---

### **6️⃣ DATABASE_AND_GRAPHS**
*Neo4j and knowledge graph documentation*

- **graph_schema.py** - Graph database schema
- **graph_loader.py** - Data loading utilities
- **graph_reasoning.py** - Query and reasoning logic
- **graph_algorithms.py** - Graph algorithms
- **load_to_neo4j.py** - Neo4j ingestion script
- **README.md** - Graph DB guide
- **QUICKSTART.py** - Quick start example

**Location**: `d:\ChemAI\graph_db/`

**ChEMBL Database**: `d:\ChemAI\chembl_36/`

---

### **7️⃣ TESTING_VALIDATION**
*Testing frameworks and validation*

- **run_tests.py** - Main test runner
- **test_bioactivity_model.py** - Bioactivity tests
- **test_toxicity_model.py** - Toxicity tests
- **test_druglikeness_model.py** - Drug-likeness tests
- **test_property_model.py** - Property tests
- **comprehensive_model_testing.py** - Full test suite
- **README_TESTING.md** - Testing guide
- **TESTING_FRAMEWORK_REFERENCE.md** - Test framework docs
- **TESTING_GUIDE.md** - How to write tests

**Location**: `d:\ChemAI\models/` and `d:\ChemAI/`

---

### **8️⃣ DEPLOYMENT_AND_DASHBOARDS**
*Running dashboards and deployment*

- **web_dashboard.py** - Web dashboard server
- **QUICKSTART_WEB_DASHBOARD.txt** - Dashboard guide
- **DASHBOARD_CSS_DESIGN.txt** - UI styling reference
- **models_list.html** - HTML templates
- **VISUAL_STATUS_SUMMARY.md** - Status visualization

**How to Run**:
```bash
python web_dashboard.py
# Visit: http://localhost:5000
```

---

### **9️⃣ REFERENCE_MATERIALS**
*Knowledge base and reference docs*

- **METRICS_GUIDE.py** - Metrics calculation guide
- **METRICS_REFERENCE.md** - Metric definitions
- **MODEL_METRICS_QUICK_VIEW.txt** - Quick metric reference
- **QUICK_REFERENCE.md** - Quick lookup guide
- **FILE_INDEX.md** - Complete file listing

---

## 📊 Key Entry Points

### **To Start Development**
1. Read: `SETUP_AND_QUICKSTART/_START_HERE_FIRST.md`
2. Install: Follow `requirements.txt`
3. Review: `TEAM_ROLES_AND_ASSIGNMENTS/README_TEAM_ORGANIZATION.txt`

### **To Understand Architecture**
1. Read: `ARCHITECTURE_AND_DESIGN/TECHNICAL_ROADMAP.md`
2. Visualize: `ARCHITECTURE_AND_DESIGN/VISUAL_FLOW_DIAGRAM.txt`
3. Review: `AGENTS_FRAMEWORK/orchestrator_agent.py`

### **To Run Experiments**
1. See: `AGENTS_FRAMEWORK/orchestrator_agent.py`
2. Run: `run_chemai.py` (main entry point)
3. Monitor: Web dashboard at `http://localhost:5000`

### **To Add New Models**
1. Create: `models/train_yourmodel.py`
2. Test: `models/test_yourmodel.py`
3. Update: Run `test_all_models.py`

---

## 🔗 File Cross-References

| Task | Read This | Then This |
|------|-----------|-----------|
| Setup project | SETUP/START_HERE | SETUP/QUICKSTART_WEB_DASHBOARD |
| Understand team roles | TEAM/README_TEAM_ORG | TEAM/TEAM_STATUS |
| Learn architecture | ARCHITECTURE/TECHNICAL_ROADMAP | AGENTS/orchestrator_agent.py |
| Train models | ML_MODELS/train_bioactivity | TESTING/comprehensive_testing |
| Run dashboard | DEPLOYMENT/QUICKSTART_DASHBOARD | web_dashboard.py |
| Query graph DB | DATABASE/README | DATABASE/QUICKSTART |

---

## 📝 File Descriptions Summary

**Python Scripts** (`d:\ChemAI\`)
- `run_chemai.py` - Main orchestration runner
- `train_all_models.py` - Train all models
- `test_all_models.py` - Test all models
- `chemai_api_server.py` - REST API server
- `web_dashboard.py` - Web UI

**Documentation** (Various `.md` and `.txt` files)
- 45+ documentation files covering every aspect
- See subfolders for organized access

**Trained Models** (`d:\ChemAI\trained_models/`)
- Saved model weights and configurations

**Data** (`d:\ChemAI\data/`)
- Datasets and CSVs for training

---

## 🎯 Next Steps by Role

### **Shreya (Orchestration)**
→ `TEAM_ROLES/SHREYA_KNOWLEDGE_GRAPH_COMPLETION.md`
→ `AGENTS_FRAMEWORK/`

### **Vishwa (RL Optimization)**
→ `TEAM_ROLES/TEAM_VISHWA/` (if exists)
→ `AGENTS_FRAMEWORK/optimizer_agent.py`

### **Pranay (GNN Generator)**
→ `TEAM_ROLES/TEAM_PRANAY/` (if exists)
→ `AGENTS_FRAMEWORK/generator_agent.py`

---

## ❓ Common Questions

**Q: How do I run the project?**
A: See `SETUP_AND_QUICKSTART/` folder

**Q: Where are the trained models?**
A: In `d:\ChemAI\trained_models/` (referenced in ML_MODELS folder)

**Q: How do I understand the agent orchestration?**
A: Read `AGENTS_FRAMEWORK/` and `ARCHITECTURE_AND_DESIGN/`

**Q: How do I contribute a new model?**
A: Follow pattern in `ML_MODELS/` and add tests in `TESTING_VALIDATION/`

**Q: How is data organized?**
A: See `DATABASE_AND_GRAPHS/` for Neo4j, `d:\ChemAI\data/` for CSVs

---

## 📞 Quick Links

- **Main Entry**: `run_chemai.py`
- **Web Dashboard**: http://localhost:5000 (after running `web_dashboard.py`)
- **Neo4j Browser**: http://localhost:7474/browser/ (if Neo4j running)
- **API Server**: http://localhost:8000 (after running `chemai_api_server.py`)

---

**Last Updated**: February 4, 2026
**Hub Version**: 1.0
**Total Documentation**: 45+ files organized in 10 categories
