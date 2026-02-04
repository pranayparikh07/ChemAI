"""
Knowledge Graph Schema Definition for ChemAI
Defines all node types and relationships for drug discovery reasoning
"""

# ============================================================================
# NODE TYPES
# ============================================================================

NODE_TYPES = {
    'Molecule': {
        'properties': {
            'smiles': 'string',              # Canonical SMILES
            'chembl_id': 'string',            # ChEMBL identifier
            'name': 'string',                 # Chemical name
            'molecular_weight': 'float',      # MW
            'logp': 'float',                  # Lipophilicity
            'qed_score': 'float',             # Drug-likeness
            'max_phase': 'int',               # Development phase
            'pIC50': 'float',                 # Bioactivity (optional)
        },
        'unique_id': 'smiles',
        'indexes': ['chembl_id', 'smiles']
    },
    
    'Protein': {
        'properties': {
            'chembl_id': 'string',            # ChEMBL target ID
            'gene_name': 'string',            # Gene name
            'uniprot_id': 'string',           # UniProt identifier
            'protein_name': 'string',         # Full protein name
            'protein_class': 'string',        # Classification (kinase, GPCR, etc.)
            'species': 'string',              # Organism
            'description': 'string',          # Protein description
        },
        'unique_id': 'uniprot_id',
        'indexes': ['chembl_id', 'gene_name', 'uniprot_id']
    },
    
    'Pathway': {
        'properties': {
            'pathway_id': 'string',           # KEGG or Reactome ID
            'pathway_name': 'string',         # Human-readable name
            'source': 'string',               # KEGG, Reactome, WikiPathways
            'description': 'string',          # Pathway description
        },
        'unique_id': 'pathway_id',
        'indexes': ['pathway_id']
    },
    
    'Disease': {
        'properties': {
            'disease_id': 'string',           # Disease identifier
            'disease_name': 'string',         # Disease name
            'mesh_id': 'string',              # MeSH identifier
            'icd_code': 'string',             # ICD code
            'description': 'string',          # Disease description
        },
        'unique_id': 'disease_id',
        'indexes': ['disease_id', 'mesh_id']
    },
    
    'Gene': {
        'properties': {
            'gene_id': 'string',              # Gene identifier
            'gene_name': 'string',            # Gene symbol
            'entrez_id': 'string',            # NCBI Entrez ID
            'uniprot_id': 'string',           # UniProt mapping
            'chromosome': 'string',           # Chromosome location
        },
        'unique_id': 'gene_id',
        'indexes': ['gene_name', 'entrez_id']
    },
    
    'SideEffect': {
        'properties': {
            'side_effect_id': 'string',       # Identifier
            'side_effect_name': 'string',     # Side effect name
            'meddra_term': 'string',          # MedDRA term
            'severity': 'string',             # Mild, Moderate, Severe
            'frequency': 'float',             # Occurrence frequency
        },
        'unique_id': 'side_effect_id',
        'indexes': ['side_effect_id']
    }
}

# ============================================================================
# EDGE TYPES (RELATIONSHIPS)
# ============================================================================

EDGE_TYPES = {
    'TARGETS': {
        'from': 'Molecule',
        'to': 'Protein',
        'properties': {
            'pIC50': 'float',                 # Potency measure
            'ic50_nM': 'float',               # IC50 in nanoMolar
            'affinity_type': 'string',        # Ki, Kd, IC50, etc.
            'confidence_score': 'float',      # Data quality (0-1)
            'pubmed_id': 'string',            # Reference
            'measurement_type': 'string',     # Binding, Activity, etc.
        },
        'direction': 'directed'
    },
    
    'INHIBITS': {
        'from': 'Molecule',
        'to': 'Protein',
        'properties': {
            'inhibition_type': 'string',      # Competitive, non-competitive
            'mechanism': 'string',            # MOA
            'potency': 'float',
            'selectivity': 'float',           # Selectivity vs other targets
        },
        'direction': 'directed'
    },
    
    'AGONIST_OF': {
        'from': 'Molecule',
        'to': 'Protein',
        'properties': {
            'ec50': 'float',                  # Activation potency
            'efficacy': 'float',              # Max response
        },
        'direction': 'directed'
    },
    
    'ANTAGONIST_OF': {
        'from': 'Molecule',
        'to': 'Protein',
        'properties': {
            'ic50': 'float',
            'selectivity': 'float',
        },
        'direction': 'directed'
    },
    
    'PART_OF': {
        'from': 'Gene',
        'to': 'Pathway',
        'properties': {
            'role': 'string',                 # Enzyme, substrate, regulator
            'evidence': 'string',             # Literature support
        },
        'direction': 'directed'
    },
    
    'ENCODES': {
        'from': 'Gene',
        'to': 'Protein',
        'properties': {
            'evidence': 'string',
            'isoform': 'string',              # Protein isoform
        },
        'direction': 'directed'
    },
    
    'TARGETS_DISEASE': {
        'from': 'Protein',
        'to': 'Disease',
        'properties': {
            'evidence': 'string',             # Research/clinical evidence
            'association_strength': 'float',  # 0-1
            'mechanism': 'string',            # How protein relates to disease
        },
        'direction': 'directed'
    },
    
    'PARTICIPATES_IN': {
        'from': 'Pathway',
        'to': 'Disease',
        'properties': {
            'significance': 'float',          # 0-1
            'evidence': 'string',
        },
        'direction': 'directed'
    },
    
    'CAUSES': {
        'from': 'Molecule',
        'to': 'SideEffect',
        'properties': {
            'frequency': 'float',             # Occurrence frequency (0-1)
            'severity': 'string',             # Mild, Moderate, Severe
            'confidence': 'float',            # Confidence in association
            'population': 'string',           # Patient population
        },
        'direction': 'directed'
    },
    
    'SIMILAR_TO': {
        'from': 'Molecule',
        'to': 'Molecule',
        'properties': {
            'similarity_score': 'float',      # 0-1
            'fingerprint_type': 'string',     # Morgan, ECFP, etc.
            'tanimoto_similarity': 'float',   # Tanimoto coefficient
        },
        'direction': 'undirected'
    }
}

# ============================================================================
# GRAPH QUERIES (CYPHER)
# ============================================================================

GRAPH_QUERIES = {
    'find_drug_targets': """
        MATCH (m:Molecule {smiles: $smiles})-[r:TARGETS]->(p:Protein)
        WHERE r.pIC50 > 6
        RETURN p.gene_name as gene, r.pIC50 as potency, r.confidence_score as confidence
        ORDER BY potency DESC
    """,
    
    'find_repurposing_opportunities': """
        MATCH (m1:Molecule {smiles: $smiles})-[r1:TARGETS]->(p1:Protein),
              (p1)-[r2:TARGETS_DISEASE]->(d:Disease)
        WHERE r1.pIC50 > 6
        RETURN DISTINCT d.disease_name as disease, p1.gene_name as target, r1.pIC50 as potency
        ORDER BY potency DESC
    """,
    
    'find_similar_active_compounds': """
        MATCH (m:Molecule {smiles: $smiles})-[s:SIMILAR_TO {similarity_score: {$gte: $threshold}}]->(m2:Molecule),
              (m2)-[t:TARGETS]->(p:Protein)
        WHERE t.pIC50 > 6
        RETURN DISTINCT m2.smiles, m2.name, collect(p.gene_name) as targets, max(t.pIC50) as max_potency
        ORDER BY max_potency DESC
    """,
    
    'predict_side_effects': """
        MATCH (m:Molecule {smiles: $smiles})-[s:SIMILAR_TO]-(m2:Molecule),
              (m2)-[c:CAUSES]->(se:SideEffect)
        WHERE s.similarity_score >= 0.7
        RETURN se.side_effect_name, count(*) as frequency, avg(c.frequency) as avg_frequency
        ORDER BY frequency DESC
    """,
    
    'rank_targets_by_importance': """
        MATCH (p:Protein)
        RETURN p.gene_name, 
               size((p)<-[r:TARGETS]-()) as num_ligands,
               size((p)-[r:TARGETS_DISEASE]->()) as num_diseases,
               (size((p)<-[r:TARGETS]-()) * size((p)-[r:TARGETS_DISEASE]-())) as importance_score
        ORDER BY importance_score DESC
        LIMIT 20
    """,
    
    'find_drug_synergies': """
        MATCH (m1:Molecule)-[t1:TARGETS]->(p1:Protein)-[:PART_OF]->(path:Pathway),
              (m2:Molecule)-[t2:TARGETS]->(p2:Protein)-[:PART_OF]->(path)
        WHERE m1.smiles < m2.smiles AND t1.pIC50 > 6 AND t2.pIC50 > 6
        RETURN m1.smiles as mol1, m2.smiles as mol2, path.pathway_name, 
               collect(DISTINCT p1.gene_name) as targets1, 
               collect(DISTINCT p2.gene_name) as targets2
    """
}

# ============================================================================
# INDEXES AND CONSTRAINTS
# ============================================================================

INDEXES = [
    # Molecule indexes
    ("Molecule", "smiles", "UNIQUE"),
    ("Molecule", "chembl_id", ""),
    
    # Protein indexes
    ("Protein", "uniprot_id", "UNIQUE"),
    ("Protein", "chembl_id", ""),
    ("Protein", "gene_name", ""),
    
    # Disease indexes
    ("Disease", "disease_id", "UNIQUE"),
    ("Disease", "mesh_id", ""),
    
    # Pathway indexes
    ("Pathway", "pathway_id", "UNIQUE"),
    
    # Gene indexes
    ("Gene", "gene_id", "UNIQUE"),
    ("Gene", "gene_name", ""),
]

# ============================================================================
# INITIALIZATION SCRIPT
# ============================================================================

INIT_CYPHER_COMMANDS = [
    # Create constraints
    "CREATE CONSTRAINT molecule_smiles IF NOT EXISTS ON (m:Molecule) ASSERT m.smiles IS UNIQUE",
    "CREATE CONSTRAINT protein_uniprot IF NOT EXISTS ON (p:Protein) ASSERT p.uniprot_id IS UNIQUE",
    "CREATE CONSTRAINT disease_id IF NOT EXISTS ON (d:Disease) ASSERT d.disease_id IS UNIQUE",
    "CREATE CONSTRAINT gene_id IF NOT EXISTS ON (g:Gene) ASSERT g.gene_id IS UNIQUE",
    "CREATE CONSTRAINT pathway_id IF NOT EXISTS ON (p:Pathway) ASSERT p.pathway_id IS UNIQUE",
    
    # Create indexes
    "CREATE INDEX molecule_chembl IF NOT EXISTS FOR (m:Molecule) ON (m.chembl_id)",
    "CREATE INDEX protein_chembl IF NOT EXISTS FOR (p:Protein) ON (p.chembl_id)",
    "CREATE INDEX protein_gene IF NOT EXISTS FOR (p:Protein) ON (p.gene_name)",
    "CREATE INDEX disease_mesh IF NOT EXISTS FOR (d:Disease) ON (d.mesh_id)",
    "CREATE INDEX gene_name IF NOT EXISTS FOR (g:Gene) ON (g.gene_name)",
]

if __name__ == "__main__":
    print("Knowledge Graph Schema Definition")
    print(f"Node Types: {len(NODE_TYPES)}")
    print(f"Edge Types: {len(EDGE_TYPES)}")
    print(f"Predefined Queries: {len(GRAPH_QUERIES)}")
    for node_type in NODE_TYPES:
        print(f"  • {node_type}: {len(NODE_TYPES[node_type]['properties'])} properties")
