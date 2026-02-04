# ChemAI - Technical Implementation Roadmap

## Overview
This document provides specific technical guidance for completing the ChemAI project's remaining work.

---

## SHREYA: Knowledge Graph Construction (IMMEDIATE PRIORITY)

### Phase 1: Technology Selection (Day 1)

#### Option A: Neo4j (Recommended for drug discovery)
**Pros:**
- Most used for biological networks
- Excellent for pattern discovery
- Built-in graph algorithms
- Easy to query with Cypher
- Good Python support

**Cons:**
- Requires installation (but can use Neo4j Aura cloud for free)
- Learning curve for Cypher

**Implementation**:
```python
# pseudo-code
from neo4j import GraphDatabase

driver = GraphDatabase.driver("neo4j://localhost:7687", auth=("neo4j", "password"))

with driver.session() as session:
    session.run(
        "CREATE (m:Molecule {smiles: $smiles, name: $name})-"
        "[r:TARGETS]->(p:Protein {name: $target})",
        smiles="CC(=O)Oc1ccccc1C(=O)O",
        name="Aspirin",
        target="COX1"
    )
```

#### Option B: RDF/OWL (Best for formal reasoning)
**Pros:**
- Semantic web standards
- Formal ontologies available (ChEMBL-RDF exists)
- SPARQL queries
- Reasoning engines available

**Cons:**
- Steeper learning curve
- Slower for large graphs
- Overkill if not doing heavy reasoning

#### Option C: Property Graph (Middle ground)
**Pros:**
- Simpler than RDF
- More flexible than relational
- Good performance

**Cons:**
- Less standardized

**Recommendation**: **USE NEO4J**
- De facto standard for biological networks
- Excellent community support
- ChEMBL already has Neo4j export available

---

### Phase 2: Graph Schema Design (Days 2-3)

#### Node Types (Required)
```
Node Types & Properties:
1. Molecule
   - smiles (string) - Primary identifier
   - iupac_name (string)
   - molecular_weight (float)
   - logp (float)
   - qed_score (float)
   - pIC50 (float) - Optional, for bioactive molecules
   - chembl_id (string)

2. Protein/Target
   - uniprot_id (string) - Primary identifier
   - gene_name (string)
   - protein_name (string)
   - protein_class (string)
   - species (string)
   - chembl_id (string)

3. Pathway
   - pathway_id (string)
   - pathway_name (string)
   - source (string) - KEGG, Reactome, etc.
   - description (string)

4. Disease
   - disease_id (string)
   - disease_name (string)
   - mesh_id (string)
   - icd_code (string)

5. Gene
   - gene_id (string)
   - gene_name (string)
   - entrez_id (string)
   - uniprot_id (string)

6. SideEffect
   - side_effect_id (string)
   - side_effect_name (string)
   - meddra_term (string)
   - severity (string)

7. Compound (Alternative to Molecule, for full ChEMBL)
   - chembl_id (string)
   - preferred_name (string)
   - pref_type (string)
   - max_phase (int)
```

#### Relationship Types (Required)
```
TARGETS
  FROM: Molecule
  TO: Protein
  PROPERTIES:
    - pIC50 (float)
    - ic50_value (float)
    - affinity_type (string) - Ki, Kd, IC50, etc.
    - confidence_score (float)
    - pubmed_id (string)

INHIBITS / AGONIST / ANTAGONIST
  FROM: Molecule
  TO: Protein
  Similar structure to TARGETS

PART_OF (Gene -> Pathway)
  FROM: Gene
  TO: Pathway
  PROPERTIES:
    - role (string)

TARGETS_DISEASE
  FROM: Protein
  TO: Disease
  PROPERTIES:
    - evidence (string)

CAUSES
  FROM: Molecule
  TO: SideEffect
  PROPERTIES:
    - frequency (float)
    - severity (string)
    - confidence (float)

ENCODES
  FROM: Gene
  TO: Protein

PARTICIPATES_IN
  FROM: Pathway
  TO: Disease

SIMILAR_TO
  FROM: Molecule
  TO: Molecule
  PROPERTIES:
    - similarity_score (float)
    - fingerprint_type (string)
```

#### Cypher Query Examples (for Phase 3)
```
# Find direct drug targets for a molecule
MATCH (m:Molecule {smiles: "CC(=O)Oc1ccccc1C(=O)O"})-[r:TARGETS]->(p:Protein)
RETURN p.gene_name, r.pIC50, r.confidence_score
ORDER BY r.confidence_score DESC

# Find drug repurposing opportunities (diseases a drug could treat)
MATCH (m:Molecule {smiles: "..."})-[r:TARGETS]->(p:Protein)-[r2:TARGETS_DISEASE]->(d:Disease)
RETURN DISTINCT d.disease_name, count(p) as num_targets
ORDER BY num_targets DESC

# Find off-targets and side effects
MATCH (m:Molecule)-[r:CAUSES]->(se:SideEffect)
WHERE r.frequency > 0.05
RETURN m.smiles, se.side_effect_name, r.frequency
ORDER BY r.frequency DESC

# Find similar compounds with known MOA
MATCH (m1:Molecule {smiles: "..."})-[s:SIMILAR_TO {similarity_score: {$gte: 0.7}}]->(m2:Molecule),
      (m2)-[t:TARGETS]->(p:Protein)-[d:TARGETS_DISEASE]->(disease:Disease)
RETURN DISTINCT m2.smiles, p.gene_name, disease.disease_name
```

---

### Phase 3: Data Loading Pipeline (Days 4-7)

#### Script: Load ChEMBL Data

```python
# File: graph_db/load_chembl_graph.py

import sqlite3
import pandas as pd
from neo4j import GraphDatabase
from tqdm import tqdm

class ChEMBLGraphLoader:
    def __init__(self, db_path="chembl_36/chembl_36_sqlite/chembl_36.db", 
                 neo4j_uri="neo4j://localhost:7687", 
                 neo4j_user="neo4j", neo4j_pass="password"):
        self.conn = sqlite3.connect(db_path)
        self.driver = GraphDatabase.driver(neo4j_uri, auth=(neo4j_user, neo4j_pass))
        
    def load_molecules(self, batch_size=1000):
        """Load molecules from ChEMBL"""
        query = """
        SELECT cd.chembl_id, cd.pref_name, cs.canonical_smiles, 
               cd.max_phase
        FROM compound_structures cs
        JOIN molecule_dictionary cd ON cs.molregno = cd.molregno
        WHERE cs.canonical_smiles IS NOT NULL
        LIMIT 50000
        """
        
        df = pd.read_sql_query(query, self.conn)
        
        for i in tqdm(range(0, len(df), batch_size)):
            batch = df.iloc[i:i+batch_size]
            self._create_molecule_batch(batch)
    
    def _create_molecule_batch(self, batch):
        """Create molecule nodes in Neo4j"""
        with self.driver.session() as session:
            for _, row in batch.iterrows():
                session.run(
                    """
                    CREATE (m:Molecule {
                        smiles: $smiles,
                        chembl_id: $chembl_id,
                        name: $name,
                        max_phase: $max_phase
                    })
                    """,
                    smiles=row['canonical_smiles'],
                    chembl_id=row['chembl_id'],
                    name=row['pref_name'],
                    max_phase=row['max_phase']
                )
    
    def load_targets(self, batch_size=500):
        """Load protein targets from ChEMBL"""
        query = """
        SELECT DISTINCT td.chembl_id, td.pref_name, td.target_type
        FROM target_dictionary td
        WHERE td.pref_name IS NOT NULL
        LIMIT 10000
        """
        
        df = pd.read_sql_query(query, self.conn)
        
        for i in tqdm(range(0, len(df), batch_size)):
            batch = df.iloc[i:i+batch_size]
            self._create_target_batch(batch)
    
    def _create_target_batch(self, batch):
        """Create protein nodes"""
        with self.driver.session() as session:
            for _, row in batch.iterrows():
                session.run(
                    """
                    CREATE (p:Protein {
                        chembl_id: $chembl_id,
                        name: $name,
                        type: $type
                    })
                    """,
                    chembl_id=row['chembl_id'],
                    name=row['pref_name'],
                    type=row['target_type']
                )
    
    def load_interactions(self, batch_size=1000):
        """Load drug-target interactions"""
        query = """
        SELECT mc.chembl_id as molecule_id, 
               td.chembl_id as target_id,
               a.standard_value, a.standard_units, a.standard_type,
               a.activity_comment
        FROM activities a
        JOIN molecule_dictionary mc ON a.molregno = mc.molregno
        JOIN target_dictionary td ON a.tid = td.tid
        WHERE a.standard_value IS NOT NULL 
        AND a.standard_type = 'IC50'
        AND a.standard_units = 'nM'
        LIMIT 50000
        """
        
        df = pd.read_sql_query(query, self.conn)
        
        for i in tqdm(range(0, len(df), batch_size)):
            batch = df.iloc[i:i+batch_size]
            self._create_interaction_batch(batch)
    
    def _create_interaction_batch(self, batch):
        """Create relationships between molecules and targets"""
        with self.driver.session() as session:
            for _, row in batch.iterrows():
                # Compute pIC50
                pIC50 = 9 - np.log10(row['standard_value'])
                
                session.run(
                    """
                    MATCH (m:Molecule {chembl_id: $mol_id})
                    MATCH (p:Protein {chembl_id: $target_id})
                    CREATE (m)-[r:TARGETS {
                        ic50_nM: $ic50,
                        pIC50: $pIC50,
                        type: $type
                    }]->(p)
                    """,
                    mol_id=row['molecule_id'],
                    target_id=row['target_id'],
                    ic50=row['standard_value'],
                    pIC50=pIC50,
                    type=row['standard_type']
                )
    
    def close(self):
        self.driver.close()
        self.conn.close()

# Usage
loader = ChEMBLGraphLoader()
loader.load_molecules()
loader.load_targets()
loader.load_interactions()
loader.close()
```

---

### Phase 4: Graph Algorithms (Days 8-12)

#### Script: Graph Algorithms

```python
# File: graph_db/graph_algorithms.py

from neo4j import GraphDatabase
import networkx as nx

class GraphAnalyzer:
    def __init__(self, neo4j_uri, neo4j_user, neo4j_pass):
        self.driver = GraphDatabase.driver(neo4j_uri, auth=(neo4j_user, neo4j_pass))
    
    def find_drug_repurposing_opportunities(self, molecule_smiles):
        """
        Find diseases a molecule could treat based on target similarity
        """
        with self.driver.session() as session:
            result = session.run("""
                MATCH (m:Molecule {smiles: $smiles})-[r:TARGETS]->(p1:Protein),
                      (p1)<-[r2:TARGETS]-(m2:Molecule)-[r3:TARGETS]->(p2:Protein),
                      (p2)-[r4:TARGETS_DISEASE]->(d:Disease)
                WHERE r.pIC50 > 6 AND r3.pIC50 > 6
                RETURN DISTINCT d.disease_name, count(DISTINCT p2) as num_common_targets,
                       avg(r3.pIC50) as avg_affinity
                ORDER BY num_common_targets DESC
                LIMIT 10
                """, smiles=molecule_smiles)
            
            return [dict(record) for record in result]
    
    def find_similar_active_compounds(self, molecule_smiles, similarity_threshold=0.7):
        """
        Find similar compounds with known bioactivity
        """
        with self.driver.session() as session:
            result = session.run("""
                MATCH (m:Molecule {smiles: $smiles})-[s:SIMILAR_TO]-(m2:Molecule),
                      (m2)-[t:TARGETS]->(p:Protein)
                WHERE s.similarity_score >= $threshold
                RETURN DISTINCT m2.smiles, m2.name, 
                       collect(p.gene_name) as targets,
                       max(t.pIC50) as max_pIC50
                ORDER BY max_pIC50 DESC
                """, smiles=molecule_smiles, threshold=similarity_threshold)
            
            return [dict(record) for record in result]
    
    def predict_toxicity_risk(self, molecule_smiles):
        """
        Predict potential toxicity based on similar compounds' side effects
        """
        with self.driver.session() as session:
            result = session.run("""
                MATCH (m:Molecule {smiles: $smiles})-[s:SIMILAR_TO]-(m2:Molecule),
                      (m2)-[c:CAUSES]->(se:SideEffect)
                WHERE s.similarity_score >= 0.7
                RETURN se.side_effect_name, count(*) as frequency,
                       avg(c.frequency) as avg_frequency
                ORDER BY frequency DESC
                LIMIT 10
                """, smiles=molecule_smiles)
            
            return [dict(record) for record in result]
    
    def rank_targets_by_importance(self):
        """
        Rank protein targets by betweenness centrality
        """
        with self.driver.session() as session:
            result = session.run("""
                MATCH (p:Protein)
                RETURN p.gene_name, p.chembl_id,
                       size((p)<-[r:TARGETS]-()) as num_ligands,
                       size((p)-[r:TARGETS_DISEASE]->()) as num_diseases
                ORDER BY num_ligands DESC, num_diseases DESC
                LIMIT 20
                """)
            
            return [dict(record) for record in result]
    
    def close(self):
        self.driver.close()
```

---

### Phase 5: Graph Embeddings (Days 13-14)

#### Install Required Package
```bash
pip install pyvis node2vec
```

#### Script: Graph Embeddings

```python
# File: graph_db/graph_embeddings.py

import numpy as np
from node2vec import Node2Vec
from sklearn.metrics.pairwise import cosine_similarity

class GraphEmbedder:
    def __init__(self, neo4j_uri, neo4j_user, neo4j_pass, embedding_dim=128):
        self.driver = GraphDatabase.driver(neo4j_uri, auth=(neo4j_user, neo4j_pass))
        self.embedding_dim = embedding_dim
        self.embeddings = {}
        self.node_id_map = {}
    
    def export_graph_to_networkx(self):
        """Export Neo4j graph to NetworkX"""
        import networkx as nx
        
        with self.driver.session() as session:
            # Get all nodes
            nodes_result = session.run("MATCH (n) RETURN id(n) as id, labels(n)[0] as label")
            # Get all edges
            edges_result = session.run("MATCH (a)--(b) RETURN id(a) as a_id, id(b) as b_id")
            
            G = nx.Graph()
            
            for record in nodes_result:
                G.add_node(record['id'], label=record['label'])
            
            for record in edges_result:
                G.add_edge(record['a_id'], record['b_id'])
            
            return G
    
    def compute_embeddings(self):
        """Compute Node2Vec embeddings"""
        G = self.export_graph_to_networkx()
        
        # Train Node2Vec
        node2vec = Node2Vec(G, dimensions=self.embedding_dim, walks_per_node=10, 
                           walk_length=80, workers=4)
        model = node2vec.fit(window=10, min_count=1, batch_words=4)
        
        # Store embeddings
        for node in G.nodes():
            self.embeddings[node] = model.wv[str(node)]
        
        return self.embeddings
    
    def find_similar_molecules(self, molecule_smiles, top_k=10):
        """Find similar molecules using embeddings"""
        with self.driver.session() as session:
            # Get molecule node ID
            result = session.run(
                "MATCH (m:Molecule {smiles: $smiles}) RETURN id(m) as node_id",
                smiles=molecule_smiles
            )
            
            mol_node_id = result.single()['node_id']
            mol_embedding = self.embeddings.get(mol_node_id)
            
            if mol_embedding is None:
                return []
            
            # Compute similarities
            similarities = {}
            for node_id, embedding in self.embeddings.items():
                sim = cosine_similarity([mol_embedding], [embedding])[0][0]
                similarities[node_id] = sim
            
            # Get top-k similar molecules
            top_nodes = sorted(similarities.items(), key=lambda x: x[1], reverse=True)[1:top_k+1]
            
            node_ids = [node_id for node_id, _ in top_nodes]
            node_ids_str = ', '.join(str(nid) for nid in node_ids)
            
            result = session.run(f"""
                MATCH (m:Molecule) WHERE id(m) IN [{node_ids_str}]
                RETURN m.smiles, m.name, m.max_phase
                """)
            
            return [(record['m.smiles'], record['m.name'], record['m.max_phase']) 
                    for record in result]
    
    def close(self):
        self.driver.close()
```

---

## VISHWA: Advanced Features (Weeks 2-3)

### Graph Feature Engineering

#### Script: Compute Graph Features

```python
# File: models/graph_feature_engineer.py

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
import numpy as np

class GraphFeatureEngineer:
    @staticmethod
    def get_node_features(smiles):
        """Get atom-level features for GNN"""
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        
        node_features = []
        for atom in mol.GetAtoms():
            features = [
                atom.GetAtomicNum(),                    # Atomic number
                atom.GetDegree(),                       # Number of neighbors
                atom.GetTotalDegree(),                  # Including hydrogens
                atom.GetExplicitValence(),              # Explicit valence
                atom.GetImplicitValence(),              # Implicit valence
                int(atom.GetIsAromatic()),              # Aromaticity
                int(atom.GetIsInRing()),                # Ring membership
                atom.GetFormalCharge(),                 # Formal charge
                atom.GetNumExplicitHs(),                # Explicit hydrogens
                atom.GetTotalNumHs(),                   # Total hydrogens
            ]
            node_features.append(features)
        
        return np.array(node_features, dtype=np.float32)
    
    @staticmethod
    def get_edge_features(smiles):
        """Get bond-level features for GNN"""
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None, None
        
        edges = []
        edge_features = []
        
        for bond in mol.GetBonds():
            i = bond.GetBeginAtomIdx()
            j = bond.GetEndAtomIdx()
            edges.append([i, j])
            edges.append([j, i])  # Bidirectional
            
            bond_type = int(bond.GetBondType())
            is_aromatic = int(bond.GetIsAromatic())
            is_conjugated = int(bond.GetIsConjugated())
            is_ring = int(bond.IsInRing())
            
            features = [bond_type, is_aromatic, is_conjugated, is_ring]
            edge_features.append(features)
            edge_features.append(features)  # Bidirectional
        
        return np.array(edges, dtype=np.int32), np.array(edge_features, dtype=np.float32)
    
    @staticmethod
    def get_topological_descriptors(smiles):
        """Compute topological descriptors"""
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        
        descriptors = {
            'Wiener_Index': Descriptors.Ipc(mol),  # Connectivity index
            'Randic_Index': Descriptors.RadiusOfGyration(mol),
            'TPSA_Per_Atom': Descriptors.TPSA(mol) / mol.GetNumAtoms(),
            'Num_Rings': Chem.GetSSSR(mol).__len__(),
            'Num_Aromatic_Rings': len([x for x in Chem.GetSSSR(mol) 
                                       if all(mol.GetAtomWithIdx(i).GetIsAromatic() 
                                              for i in x)]),
        }
        
        return descriptors
```

---

## PRANAY: Research Documentation (Week 1-2)

### Document 1: Literature Review Template

```markdown
# Literature Review: AI-Driven Molecular Discovery

## 1. Molecular Generation Methods
### 1.1 Variational Autoencoders (VAE)
- Gomez-Bombarelli et al. (2018) - Automatic chemical design using a data-driven continuous...
- [5+ more papers]

### 1.2 Generative Adversarial Networks (GAN)
- De Cao & Kipf (2018) - MolGAN
- [5+ more papers]

### 1.3 Reinforcement Learning
- Popova et al. (2018) - Deep reinforcement learning for de novo drug discovery
- [5+ more papers]

### 1.4 Genetic Algorithms & Evolutionary Methods
- [5+ papers]

## 2. Agentic Systems in Chemistry
- [10+ papers on multi-agent systems, autonomy levels, etc.]

## 3. Knowledge Graphs in Drug Discovery
- [10+ papers on KG construction, reasoning, etc.]

## 4. Causal Reasoning in Chemistry
- [5+ papers on causal inference, do-calculus applications]

## 5. Gaps & Opportunities
- Lack of formal causal reasoning ✓
- Limited autonomous decision loops ✓
- Insufficient knowledge integration ✓
```

---

## Implementation Timeline

```
WEEK 1 (Shreya - Critical Path):
  Mon-Tue: Neo4j setup + Graph schema design
  Wed-Thu: ChEMBL data loading pipeline
  Fri: Basic queries & verification

WEEK 2 (Shreya + Pranay):
  Mon-Tue: Graph algorithms implementation
  Wed-Thu: Literature review (Pranay begins)
  Fri: Graph algorithms testing

WEEK 3 (Shreya + Vishwa):
  Mon-Tue: Graph embeddings
  Wed-Thu: Advanced features (Vishwa)
  Fri: Integration testing

WEEK 4 (Vishwa + Pranay):
  Mon-Tue: Multi-dataset support
  Wed-Thu: Research documentation completion
  Fri: All documentation review

WEEK 5 (All):
  Full system integration & validation
```

---

## Success Checklist

**Shreya - Knowledge Graph:**
- [ ] Neo4j instance running
- [ ] 50K+ molecules in graph
- [ ] 10K+ proteins in graph
- [ ] 100K+ drug-target interactions
- [ ] Graph queries working
- [ ] Embeddings computed
- [ ] Similarity search functional

**Vishwa - Features:**
- [ ] Graph node features computed
- [ ] Graph edge features computed
- [ ] Topological descriptors available
- [ ] QM9 dataset loaded
- [ ] ZINC dataset loaded
- [ ] Feature importance ranking done

**Pranay - Documentation:**
- [ ] Literature review complete (50+ citations)
- [ ] Hypotheses formally stated
- [ ] Gap analysis written
- [ ] Experimental design documented
- [ ] Expected performance targets defined

