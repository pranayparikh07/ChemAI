# 🚀 SETUP AND QUICKSTART

**Get ChemAI running in minutes**

---

## 📂 Files in This Folder

| File | Purpose | Status |
|------|---------|--------|
| `_START_HERE_FIRST.md` | Project overview & first steps | Reference |
| `START_HERE.md` | Installation & setup guide | Setup Guide |
| `QUICKSTART_WEB_DASHBOARD.txt` | How to launch dashboard | Tutorial |
| `README_CHEMAI.md` | Full project description | Reference |
| `requirements.txt` | Python dependencies | Required |

---

## ⚡ Quick Setup (5 minutes)

### 1. **Install Dependencies**
```bash
cd d:\ChemAI
pip install -r requirements.txt
```

### 2. **Run Main Orchestrator**
```bash
python run_chemai.py
```

### 3. **Launch Dashboard**
```bash
python web_dashboard.py
# Open: http://localhost:5000
```

---

## 📖 Step-by-Step

1. **Read**: `_START_HERE_FIRST.md` (overview)
2. **Read**: `START_HERE.md` (installation)
3. **Read**: `QUICKSTART_WEB_DASHBOARD.txt` (dashboard)
4. **Run**: `python run_chemai.py`
5. **Visit**: `http://localhost:5000`

---

## 🔧 Common Commands

```bash
# Train all models
python train_all_models.py

# Test all models
python test_all_models.py

# Run orchestrator
python run_chemai.py

# Start API server
python chemai_api_server.py

# Start web dashboard
python web_dashboard.py

# Load to Neo4j
python TEAM_SHREYA/load_to_neo4j.py
```

---

## 📍 Main Entry Points

- **`run_chemai.py`** - Main orchestration (Shreya)
- **`web_dashboard.py`** - Web interface
- **`chemai_api_server.py`** - REST API

---

## ✅ Checklist

- [ ] Install requirements.txt
- [ ] Run run_chemai.py once
- [ ] Access dashboard at localhost:5000
- [ ] Train models (or use pre-trained)
- [ ] Read TEAM_ROLES_AND_ASSIGNMENTS to understand structure

---

**Next**: Go to `TEAM_ROLES_AND_ASSIGNMENTS/` to understand team structure
