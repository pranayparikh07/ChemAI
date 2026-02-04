# Knowledge Graph Implementation Guide

## Overview
This guide provides step-by-step instructions to implement the ChemAI Knowledge Graph for drug discovery reasoning and prediction explanation.

## Architecture

```
ChemAI Knowledge Graph System
├── Graph Schema (graph_schema.py)
│   ├── Node Types: 6 (Molecule, Protein, Pathway, Disease, Gene, SideEffect)
│   ├── Edge Types: 10 (TARGETS, INHIBITS, AGONIST_OF, etc.)
│   └── Queries: 6 predefined reasoning queries
│
├── Data Layer (graph_loader.py)
│   ├── ChEMBL Loading: 50K molecules + targets
│   ├── Interaction Loading: Drug-target interactions
│   └── Similarity: Compound similarity edges
│
├── Reasoning Layer (graph_algorithms.py)
│   ├── Shortest Path: Drug repurposing
│   ├── Centrality: Target importance ranking
│   ├── Similarity: Compound/target analogs
│   └── Polypharmacology: Drug synergies
│
└── Explanation Layer (graph_reasoning.py)
    ├── Predict Explanation: Why does drug hit target?
    ├── Activity Explanation: Why is compound active?
    ├── Repurposing Hypothesis: New indications
    └── Combination Hypothesis: Drug synergies
```

## Prerequisites

### 1. Install Neo4j
**Option A: Local Installation**
```bash
# Download from https://neo4j.com/download/
# Or use Docker
docker run -p 7687:7687 -p 7474:7474 \
  -e NEO4J_AUTH=neo4j/your_password \
  neo4j:latest
```

**Option B: Neo4j AuraDB (Managed Service)**
```
Visit: https://neo4j.com/cloud/aura/
Create free instance
Copy connection URI and credentials
```

### 2. Python Dependencies
```bash
pip install neo4j==5.x
pip install pandas
pip install numpy
pip install rdkit
```

## Step 1: Schema Initialization

### Code
```python
from graph_schema import NODE_TYPES, EDGE_TYPES, INIT_CYPHER_COMMANDS

# Connect to Neo4j
from neo4j import GraphDatabase

driver = GraphDatabase.driver("neo4j://localhost:7687", 
                              auth=("neo4j", "password"))

# Initialize graph
with driver.session() as session:
    for command in INIT_CYPHER_COMMANDS:
        session.run(command)
```

### What This Creates
- 6 node types with properties
- 10 relationship types with causality attributes
- 5 unique constraints (SMILES, UniProt, Disease ID, Gene ID, Pathway ID)
- 5 indexes for fast lookups

### Validation
```cypher
# Check nodes created
MATCH (n) RETURN labels(n), count(n)

# Check relationships
MATCH ()-[r]->() RETURN type(r), count(r)

# Check constraints
CALL db.constraints()
```

## Step 2: Data Loading

### Basic Loading
```python
from graph_loader import ChEMBLGraphLoader

loader = ChEMBLGraphLoader()
loader.connect_neo4j("neo4j://localhost:7687", "neo4j", "password")
loader.init_graph()

# Load data
loader.load_molecules(50000)      # ~5-10 minutes
loader.load_proteins(10000)       # ~2-3 minutes
loader.load_interactions(50000)   # ~10-15 minutes
loader.create_similarity_edges(0.7)  # ~5 minutes

# Get statistics
stats = loader.get_statistics()
print(f"Loaded: {stats}")
```

### Data Volumes
| Entity | Count | Notes |
|--------|-------|-------|
| Molecules | 50,000 | With SMILES and properties |
| Proteins | 10,000 | Drug targets |
| Interactions | 50,000+ | Drug-target with pIC50 |
| Similarity Edges | Variable | Tanimoto > 0.7 |

### Loading Pipeline Time
- Phase 1 (Schema): < 1 minute
- Phase 2 (Molecules): 5-10 minutes
- Phase 3 (Proteins): 2-3 minutes
- Phase 4 (Interactions): 10-15 minutes
- Phase 5 (Similarities): 5 minutes
- **Total: ~25-35 minutes**

## Step 3: Algorithms Implementation

### Available Algorithms

#### A. Shortest Path (Drug Repurposing)
```python
from graph_algorithms import GraphAlgorithms

algo = GraphAlgorithms(driver)

# Find path from drug to disease
paths = algo.find_shortest_path(
    start_node="CHEMBL123456",      # Known drug
    end_node="disease_123",         # New disease
    max_depth=5
)
```

**Use Case**: Identify if a known drug could treat a new disease via indirect targets

#### B. Centrality Ranking (Target Importance)
```python
# Find hub proteins
degree_centrality = algo.compute_degree_centrality("Protein")
# Returns: [(target_id, name, degree), ...]

# Rank targets for a disease
ranked = algo.rank_targets_by_importance("Alzheimer's disease")
# Returns: Importance-ranked protein list
```

**Use Case**: Identify the best targets to drug in a disease

#### C. Similarity Search (Compound SAR)
```python
# Find structural analogs
similar = algo.find_similar_compounds(
    molecule_id="SMILES_STRING",
    min_similarity=0.7,
    limit=20
)
```

**Use Case**: Structure-activity relationship analysis

#### D. Polypharmacology (Drug Synergies)
```python
# Find synergistic combinations
synergies = algo.find_drug_synergies("SMILES_STRING")
```

**Use Case**: Identify molecules that could work well together

## Step 4: Reasoning & Explanation

### Prediction Explanation
```python
from graph_reasoning import ReasoningEngine

engine = ReasoningEngine(driver)

# Explain a prediction
explanation = engine.explain_drug_target_prediction(
    molecule_id="SMILES_or_ChEMBL_ID",
    target_id="Protein_ChEMBL_ID"
)

print(explanation['evidence'])        # Supporting evidence
print(explanation['confidence'])      # Confidence score
print(explanation['mechanistic_reasoning'])  # Explanation text
```

**Output Example**:
```python
{
    'evidence': [
        {
            'type': 'structural_analogy',
            'description': 'Similar to known binder CHEMBL456 (similarity: 0.82, pIC50: 7.5)'
        },
        {
            'type': 'pathway_connectivity',
            'description': 'Participates in PI3K-AKT pathway (interconnection: 12)'
        }
    ],
    'confidence': 0.82,
    'mechanistic_reasoning': '2 pieces of evidence...'
}
```

### Repurposing Hypothesis
```python
hypothesis = engine.generate_repurposing_hypothesis("CHEMBL123456")

print(hypothesis['known_indications'])      # Current uses
print(hypothesis['repurposing_candidates']) # Potential new uses
```

**Output Example**:
```python
{
    'known_indications': [
        {'disease': 'Cancer', 'target': 'EGFR', 'strength': 8.2}
    ],
    'repurposing_candidates': [
        {
            'disease': 'Fibrosis',
            'targets': ['FGFR1', 'FGFR2'],
            'pathway': 'FGF-FGFR signaling'
        }
    ]
}
```

### Combination Therapy
```python
hypothesis = engine.generate_combination_therapy_hypothesis([
    "SMILES_1",
    "SMILES_2"
])
```

## Step 5: Integration with ML Pipeline

### Enhance Predictions
```python
# ML model predicts
ml_prediction = {
    'molecule_id': 'SMILES',
    'target_id': 'Protein_ID',
    'predicted_pIC50': 7.5,
    'confidence': 0.72
}

# Validate with graph
validation = engine.validate_prediction_with_graph(ml_prediction)

# Adjusted confidence
print(validation['final_confidence'])  # May be 0.85+ with graph support
print(validation['graph_evidence'])    # Why model prediction is credible
```

## Performance Considerations

### Graph Size
| Nodes | Edges | Storage | Query Time |
|-------|-------|---------|-----------|
| 60K | 100K | ~500MB | < 1sec |
| 200K | 1M | ~2GB | < 2sec |
| 1M | 5M | ~10GB | < 5sec |

### Optimization Tips
1. **Index creation**: Already done in schema initialization
2. **Batch processing**: Load in batches of 1000
3. **Query optimization**: Use MATCH with WHERE filters first
4. **Read replicas**: For large-scale queries, use Neo4j Enterprise

### Scaling to Full ChEMBL (2M+ molecules)
```python
# Adjust batch sizes
loader.load_molecules(2000000)  # Full ChEMBL
# Will take 4-8 hours
# Requires Neo4j Enterprise for best performance
```

## Troubleshooting

### 1. Connection Issues
```python
# Test connection
try:
    driver.verify_connectivity()
    print("✓ Connected")
except Exception as e:
    print(f"✗ Error: {e}")
```

### 2. Data Loading Slow
```python
# Check database statistics
driver.session().run("CALL dbms.queryJmx('*:*') YIELD attributes RETURN attributes")

# Increase batch sizes
loader = ChEMBLGraphLoader()
# Modify batch_size=5000 in load_molecules()
```

### 3. Out of Memory
```bash
# Increase Neo4j heap
export NEO4J_HEAP_MEMORY=8g
./bin/neo4j start
```

### 4. Similarity Computation Slow
```python
# Use fingerprint indexing instead of full computation
# Or limit to subset
loader.create_similarity_edges(0.8, limit=10000)
```

## Testing

### Unit Tests
```python
# Test schema
def test_schema():
    with driver.session() as session:
        result = session.run("MATCH (n:Molecule) RETURN count(n)")
        assert result.single()['count'] > 0

# Test algorithms
def test_algorithms():
    algo = GraphAlgorithms(driver)
    results = algo.compute_degree_centrality()
    assert len(results) > 0

# Test reasoning
def test_reasoning():
    engine = ReasoningEngine(driver)
    stats = engine.get_graph_statistics()
    assert stats['molecules'] > 0
```

### Integration Tests
```python
# End-to-end prediction validation
prediction = {
    'molecule_id': 'CHEMBL1',
    'target_id': 'CHEMBL2',
    'predicted_pIC50': 7.5,
    'confidence': 0.75
}

validation = engine.validate_prediction_with_graph(prediction)
assert validation['final_confidence'] > prediction['confidence']
```

## Next Steps

1. **Phase 1** (DONE): Schema initialization ✓
2. **Phase 2** (NOW): Data loading → Complete first with subset
3. **Phase 3**: Graph algorithms → Implement per use case
4. **Phase 4**: Graph embeddings (Node2Vec, DeepWalk) → For similarity
5. **Phase 5**: Causal reasoning → Advanced mechanistic insights

## Estimated Completion

| Phase | Tasks | Time | Status |
|-------|-------|------|--------|
| 1 | Schema, Init | 1 hour | ✓ DONE |
| 2 | Data loading | 1-2 hours | → IN PROGRESS |
| 3 | Algorithms | 2-3 hours | NEXT |
| 4 | Embeddings | 2-3 hours | FUTURE |
| 5 | Causal + Testing | 2-3 hours | FUTURE |
| **TOTAL** | **Full system** | **9-14 hours** | **50% COMPLETE** |

## Resources

- **Neo4j Documentation**: https://neo4j.com/docs/
- **Cypher Query Language**: https://neo4j.com/docs/cypher-manual/
- **RDKit Cheminformatics**: https://rdkit.org/
- **Drug Discovery Graphs**: https://academic.oup.com/bib/article/21/3/855/5826748

---

**Created for ChemAI Project - Shreya's Knowledge Graph Work (50% Complete)**
