# 📊 DATABASE AND GRAPHS

**Neo4j Knowledge Graph and Data Management**

---

## 📂 Files in This Folder

| File | Purpose | Location |
|------|---------|----------|
| `graph_schema.py` | Neo4j schema definition | `graph_db/` |
| `graph_loader.py` | Data loading utilities | `graph_db/` |
| `graph_reasoning.py` | Query and reasoning logic | `graph_db/` |
| `graph_algorithms.py` | Graph algorithms | `graph_db/` |
| `load_to_neo4j.py` | Neo4j ingestion script | `TEAM_SHREYA/` |
| `README.md` | Graph DB guide | `graph_db/` |
| `QUICKSTART.py` | Quick start example | `graph_db/` |

---

## 🗄️ Neo4j Database Setup

### **Start Neo4j**
```bash
# Using Docker
docker run -d \
  -p 7474:7474 \
  -p 7687:7687 \
  -e NEO4J_AUTH=neo4j/password \
  neo4j:latest

# Then access at: http://localhost:7474/browser/
```

### **Load Data to Neo4j**
```bash
cd TEAM_SHREYA
python load_to_neo4j.py
# Enter URI: bolt://localhost:7687
# Enter username: neo4j
# Enter password: [your password]
```

---

## 📐 Graph Schema

### **Nodes**

#### **Molecule Node**
```cypher
:Molecule {
  id: "mol_001",          # Unique identifier
  smiles: "CC(C)Cc1ccc...", # SMILES string
  name: "Ibuprofen",      # Drug name
  mw: 206.28,             # Molecular weight
  logp: 3.97,             # Lipophilicity
  qed: 0.83,              # Drug-likeness
  is_drug_like: true,     # Boolean flag
  source: "ChEMBL",       # Data origin
  confidence: 0.95        # Data confidence
}
```

#### **Protein Node**
```cypher
:Protein {
  id: "prot_001",
  uniprot: "P23458",      # UniProt ID
  name: "Cyclooxygenase-2",
  gene: "PTGS2",
  organism: "Homo sapiens",
  function: "Prostaglandin synthesis",
  source: "ChEMBL",
  confidence: 0.95
}
```

#### **Experiment Node**
```cypher
:Experiment {
  id: "exp_001",
  timestamp: "2026-02-04 10:30:45",
  parameters: {
    num_generations: 5,
    molecules_per_generation: 100
  },
  status: "completed"
}
```

### **Relationships**

#### **TARGETS (Molecule → Protein)**
```cypher
(m:Molecule)-[:TARGETS {
  pic50: 5.2,             # pIC50 value (potency)
  activity_type: "IC50",
  confidence: 0.92,
  source: "ChEMBL"
}]->(p:Protein)
```

#### **SIMILAR_TO (Molecule → Molecule)**
```cypher
(m1:Molecule)-[:SIMILAR_TO {
  similarity: 0.82,       # Tanimoto similarity
  method: "Tanimoto",
  confidence: 0.95,
  source: "Computed"
}]->(m2:Molecule)
```

#### **GENERATED_IN (Molecule → Experiment)**
```cypher
(m:Molecule)-[:GENERATED_IN {
  generation_number: 2,
  score: 0.89,
  source: "Agent"
}]->(e:Experiment)
```

---

## 🔍 Important Cypher Queries

### **Find Active Molecules**
```cypher
MATCH (m:Molecule)-[r:TARGETS]->(p:Protein {name: "Cyclooxygenase-2"})
WHERE r.pic50 > 5.0
RETURN m.name, r.pic50, r.confidence
ORDER BY r.pic50 DESC
LIMIT 10
```

### **Find Similar Molecules**
```cypher
MATCH (m1:Molecule {name: "Ibuprofen"})-[r:SIMILAR_TO]->(m2:Molecule)
WHERE r.similarity > 0.8
RETURN m2.name, m2.smiles, r.similarity
ORDER BY r.similarity DESC
```

### **Analyze Generation Performance**
```cypher
MATCH (m:Molecule)-[r:GENERATED_IN]->(e:Experiment)
WITH e.id as experiment_id, 
     AVG(r.score) as avg_score,
     MAX(r.score) as best_score,
     COUNT(m) as num_molecules
RETURN experiment_id, avg_score, best_score, num_molecules
```

### **Drug-like Molecule Discovery**
```cypher
MATCH (m:Molecule)
WHERE m.qed > 0.7 
  AND m.mw < 500 
  AND m.logp < 5
RETURN m.name, m.qed, m.mw, m.logp
ORDER BY m.qed DESC
LIMIT 20
```

---

## 💾 Database Operations

### **Create Indexes**
```cypher
CREATE INDEX molecule_smiles FOR (m:Molecule) ON (m.smiles);
CREATE INDEX protein_name FOR (p:Protein) ON (p.name);
CREATE INDEX molecule_qed FOR (m:Molecule) ON (m.qed);
```

### **Create Constraints**
```cypher
CREATE CONSTRAINT molecule_id FOR (m:Molecule) REQUIRE m.id IS UNIQUE;
CREATE CONSTRAINT protein_id FOR (p:Protein) REQUIRE p.id IS UNIQUE;
```

### **Backup Database**
```bash
neo4j-admin database dump neo4j backup.dump
```

### **Restore Database**
```bash
neo4j-admin database load neo4j backup.dump --overwrite
```

---

## 📊 Database Statistics

### **Query to Check Stats**
```cypher
MATCH (n) RETURN labels(n) as Node_Type, count(*) as Count
UNION ALL
MATCH ()-[r]->() RETURN type(r) as Relationship_Type, count(*) as Count
```

### **Expected Sizes**
- **Molecules**: ~50,000+
- **Proteins**: ~5,000+
- **TARGETS relationships**: ~500,000+
- **SIMILAR_TO relationships**: ~100,000+

---

## 🔗 Loading Process

### **Step 1: Clear Database (Optional)**
```python
loader.clear_database()
```

### **Step 2: Create Constraints**
```python
loader.create_constraints()
```

### **Step 3: Load Molecules**
```python
loader.load_molecules('data/molecules_kg.csv')
```

### **Step 4: Load Proteins**
```python
loader.load_proteins('data/proteins_kg.csv')
```

### **Step 5: Load Relationships**
```python
loader.load_bioactivity_edges('data/bioactivity_edges.csv')
loader.load_similarity_edges('data/similarity_edges.csv')
```

### **Step 6: View Statistics**
```python
loader.get_statistics()
```

---

## 🧬 ChEMBL Database

### **Location**
```
d:\ChemAI\chembl_36\chembl_36_sqlite\
```

### **Access ChEMBL Data**
```python
import sqlite3

conn = sqlite3.connect('chembl_36/chembl_36_sqlite/chembl_36.db')
cursor = conn.cursor()

# Get molecules
cursor.execute("SELECT DISTINCT chembl_id, canonical_smiles FROM compound_structures LIMIT 10")
molecules = cursor.fetchall()

# Get bioactivity data
cursor.execute("""
  SELECT * FROM bioassay_result 
  WHERE standard_value > 0 LIMIT 10
""")
bioactivity = cursor.fetchall()
```

---

## 📈 Performance Optimization

### **For Large Queries**
1. Use indexes (already created)
2. Limit result sets
3. Use APOC procedures for complex operations

### **Example: Efficient Similarity Search**
```cypher
MATCH (m1:Molecule {id: "mol_001"})-[r:SIMILAR_TO]->(m2:Molecule)
WHERE r.similarity > 0.7
RETURN m2 LIMIT 100
```

---

## 🛠️ Maintenance Tasks

### **Weekly**
- Backup database
- Check query performance
- Monitor disk usage

### **Monthly**
- Rebuild indexes
- Analyze query logs
- Update statistics

### **Quarterly**
- Full database optimization
- Update molecule data
- Archive old experiments

---

## 📖 Related Files

- `graph_schema.py` - Detailed schema definition
- `QUICKSTART.py` - Quick start example
- `graph_algorithms.py` - Advanced operations

---

## 🔗 Integration with Agents

### **How Shreya Uses Neo4j**

```python
from graph_db.graph_loader import GraphLoader

loader = GraphLoader(uri, username, password)

# After each generation
loader.create_experiment_node(experiment_data)
loader.add_molecules(generated_molecules)
loader.create_targets_edges(bioactivity_data)
loader.log_decision(decision_trace)

# For analysis
top_molecules = loader.query_drug_like_molecules(qed_threshold=0.7)
```

---

**Master Index**: Go back to `../../INDEX.md`
