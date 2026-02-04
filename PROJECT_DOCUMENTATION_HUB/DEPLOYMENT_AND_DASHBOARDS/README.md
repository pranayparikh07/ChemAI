# 📊 DEPLOYMENT AND DASHBOARDS

**Web Dashboard, API Deployment, and Monitoring**

---

## 📂 Files in This Folder

| File | Purpose |
|------|---------|
| `web_dashboard.py` | Flask-based web interface |
| `QUICKSTART_WEB_DASHBOARD.txt` | Dashboard setup guide |
| `DASHBOARD_CSS_DESIGN.txt` | UI styling reference |
| `models_list.html` | HTML templates |

---

## 🚀 Web Dashboard

### **Start Dashboard**
```bash
python web_dashboard.py
# Open: http://localhost:5000
```

### **Features**
- 📊 Real-time molecule visualization
- 📈 Performance metrics display
- 🧬 Interactive molecule explorer
- 📋 Experiment history
- 📊 Statistical summaries

---

## 🖥️ Dashboard Pages

### **1. Home Page**
- Project overview
- Quick statistics
- Recent experiments
- Quick access buttons

**URL**: `http://localhost:5000/`

### **2. Molecules Explorer**
- Browse all molecules
- Filter by properties
- View chemical structures
- Search by SMILES

**URL**: `http://localhost:5000/molecules`

### **3. Models Dashboard**
- Model performance metrics
- Training status
- Prediction accuracy
- Comparison charts

**URL**: `http://localhost:5000/models`

### **4. Experiments**
- List all experiments
- View experiment details
- Download results
- Compare experiments

**URL**: `http://localhost:5000/experiments`

### **5. Neo4j Graph Browser**
- Interactive graph visualization
- Molecule-protein networks
- Similarity networks
- Query interface

**URL**: (Integrated or link to localhost:7474)

---

## ⚙️ Configuration

### **Dashboard Settings**
```python
# In web_dashboard.py
app = Flask(__name__)
app.config['UPLOAD_FOLDER'] = 'uploads/'
app.config['MAX_CONTENT_LENGTH'] = 16 * 1024 * 1024  # 16MB max

# Database
DB_URI = 'sqlite:///chemai.db'
NEO4J_URI = 'bolt://localhost:7687'

# Flask settings
DEBUG = True
PORT = 5000
HOST = '0.0.0.0'
```

---

## 🎨 UI Components

### **Molecule Viewer**
```html
<!-- Display molecular structure -->
<canvas id="molecule-canvas"></canvas>

<!-- Molecule properties -->
<div class="molecule-properties">
  <p>SMILES: <code>CC(C)Cc1ccc...</code></p>
  <p>Bioactivity: <strong>7.5</strong></p>
  <p>Toxicity: <strong>0.2</strong></p>
  <p>Drug-likeness: <strong>0.87</strong></p>
</div>
```

### **Metrics Chart**
```html
<div id="metrics-chart"></div>
<script>
  // Chart.js for visualization
  var ctx = document.getElementById('metrics-chart').getContext('2d');
  var chart = new Chart(ctx, {
    type: 'bar',
    data: {...}
  });
</script>
```

### **Experiment Timeline**
```html
<div class="timeline">
  <div class="event">
    <p class="title">Generation Phase</p>
    <p class="time">10:30 AM</p>
    <p class="detail">Generated 100 molecules</p>
  </div>
  ...
</div>
```

---

## 📡 REST API Endpoints

### **Start API Server**
```bash
python chemai_api_server.py
# Base URL: http://localhost:8000
```

### **Available Endpoints**

#### **Molecules**
```
GET    /api/molecules              # List all
GET    /api/molecules/<id>         # Get one
POST   /api/molecules              # Create new
PUT    /api/molecules/<id>         # Update
DELETE /api/molecules/<id>         # Delete
```

#### **Predictions**
```
POST   /api/predict               # Predict properties
POST   /api/predict/batch         # Batch prediction
GET    /api/predict/status/<job>  # Check status
```

#### **Experiments**
```
GET    /api/experiments           # List experiments
GET    /api/experiments/<id>      # Get details
POST   /api/experiments           # Create new
GET    /api/experiments/<id>/results  # Get results
```

---

## 🔌 API Usage Examples

### **Get Molecule Info**
```bash
curl http://localhost:8000/api/molecules/mol_001
```

**Response**:
```json
{
  "id": "mol_001",
  "smiles": "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
  "name": "Ibuprofen",
  "properties": {
    "mw": 206.28,
    "logp": 3.97,
    "qed": 0.83
  }
}
```

### **Predict Properties**
```bash
curl -X POST http://localhost:8000/api/predict \
  -H "Content-Type: application/json" \
  -d '{"smiles": "CC(C)Cc1ccc(cc1)C(C)C(=O)O"}'
```

**Response**:
```json
{
  "smiles": "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
  "predictions": {
    "bioactivity": 7.5,
    "toxicity": 0.2,
    "druglikeness": 0.87
  }
}
```

---

## 📈 Performance Monitoring

### **Metrics Dashboard**
- Request latency
- API throughput
- Model inference time
- Database query time
- Cache hit rate

### **Logging**
```python
import logging

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('logs/app.log'),
        logging.StreamHandler()
    ]
)
```

---

## 🐳 Docker Deployment

### **Dockerfile**
```dockerfile
FROM python:3.9-slim

WORKDIR /app

COPY requirements.txt .
RUN pip install -r requirements.txt

COPY . .

EXPOSE 5000 8000

CMD ["python", "web_dashboard.py"]
```

### **Docker Compose**
```yaml
version: '3.8'

services:
  chemai:
    build: .
    ports:
      - "5000:5000"
      - "8000:8000"
    environment:
      - FLASK_ENV=production
      - NEO4J_URI=bolt://neo4j:7687
    depends_on:
      - neo4j

  neo4j:
    image: neo4j:latest
    ports:
      - "7474:7474"
      - "7687:7687"
    environment:
      - NEO4J_AUTH=neo4j/password
```

### **Run with Docker**
```bash
docker-compose up -d
# Access at: http://localhost:5000
```

---

## 🔐 Security

### **HTTPS Configuration**
```python
from flask import Flask
from flask_talisman import Talisman

app = Flask(__name__)
Talisman(app)  # Enables HTTPS security headers
```

### **Authentication**
```python
from flask_login import LoginManager, login_required

login_manager = LoginManager()
login_manager.init_app(app)

@app.route('/experiments')
@login_required
def experiments():
    return render_template('experiments.html')
```

---

## 📊 Monitoring Dashboard

### **Check Dashboard Status**
```bash
# Logs
tail -f logs/app.log

# System resources
watch -n 1 'free -h && df -h'

# API metrics
curl http://localhost:8000/metrics
```

---

## 🚨 Troubleshooting

### **Port Already in Use**
```bash
# Find process using port 5000
lsof -i :5000

# Kill it
kill -9 <PID>
```

### **Database Connection Error**
```python
# Check Neo4j connection
python -c "from neo4j import GraphDatabase; \
  driver = GraphDatabase.driver('bolt://localhost:7687'); \
  print('Connected' if driver.verify_connectivity() else 'Failed')"
```

### **Model Loading Error**
```bash
# Verify models exist
ls -la trained_models/

# Test model import
python -c "from agents.predictor_agent import PredictorAgent; \
  p = PredictorAgent(); \
  print('Models loaded successfully')"
```

---

## 📊 Sample Screenshots

### **Dashboard Overview**
```
┌────────────────────────────────────┐
│    ChemAI Dashboard                │
├────────────────────────────────────┤
│ Stats:                             │
│ • Total Molecules: 5,234          │
│ • Recent Experiments: 12           │
│ • Active Agents: 5                │
│                                    │
│ [Molecules] [Models] [Experiments] │
└────────────────────────────────────┘
```

---

## 📖 Related Files

- `QUICKSTART_WEB_DASHBOARD.txt` - Setup guide
- `DASHBOARD_CSS_DESIGN.txt` - Styling reference
- `models_list.html` - HTML templates

---

**Master Index**: Go back to `../../INDEX.md`
