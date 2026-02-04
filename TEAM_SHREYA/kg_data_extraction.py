#!/usr/bin/env python3
"""
Knowledge Graph Data Extraction Scripts
Extract molecules, proteins, bioactivity, and relationships from multiple sources
"""

import pandas as pd
import sqlite3
import json
from pathlib import Path


def extract_molecules_from_chembl(db_path, limit=50000):
    """Extract molecules from ChEMBL database"""
    query = """
    SELECT DISTINCT
        md.molregno as molecule_id,
        cs.canonical_smiles as smiles,
        md.pref_name as name,
        'ChEMBL' as source
    FROM molecule_dictionary md
    INNER JOIN compound_structures cs ON md.molregno = cs.molregno
    LIMIT ?
    """
    conn = sqlite3.connect(db_path)
    df = pd.read_sql_query(query, conn, params=(limit,))
    conn.close()
    return df


def extract_bioactivity_from_chembl(db_path):
    """Extract bioactivity measurements (IC50, pIC50)"""
    query = """
    SELECT DISTINCT
        md.molregno as molecule_id,
        t.tid as protein_id,
        t.pref_name as protein_name,
        act.standard_type,
        act.standard_value,
        act.standard_units,
        a.assay_type,
        'ChEMBL' as source
    FROM molecule_dictionary md
    INNER JOIN activities act ON md.molregno = act.molregno
    INNER JOIN assays a ON act.assay_id = a.assay_id
    INNER JOIN target_dictionary t ON a.tid = t.tid
    WHERE act.standard_type IN ('IC50', 'EC50', 'Ki', 'Kd')
        AND act.standard_value IS NOT NULL
        AND act.standard_units = 'nM'
    """
    conn = sqlite3.connect(db_path)
    df = pd.read_sql_query(query, conn)
    conn.close()
    
    # Calculate pIC50
    import numpy as np
    df['pic50'] = -np.log10(df['standard_value'] * 1e-9)
    
    return df


def extract_targets_from_chembl(db_path):
    """Extract protein target information"""
    query = """
    SELECT DISTINCT
        t.tid as protein_id,
        t.pref_name as name,
        t.target_type,
        'ChEMBL' as source
    FROM target_dictionary t
    WHERE t.target_type IN ('SINGLE PROTEIN', 'PROTEIN FAMILY', 'CHIMERIC PROTEIN')
    """
    conn = sqlite3.connect(db_path)
    df = pd.read_sql_query(query, conn)
    conn.close()
    return df


def create_molecule_nodes(molecules_df, properties_df=None):
    """Create molecule nodes for KG"""
    nodes = []
    
    for idx, row in molecules_df.iterrows():
        node = {
            'id': f"mol_{row['molecule_id']}",
            'type': 'Molecule',
            'properties': {
                'smiles': row['smiles'],
                'name': row['name'],
                'source': row['source'],
                'source_id': row['molecule_id'],
            }
        }
        
        # Add properties if available
        if properties_df is not None and 'smiles' in properties_df.columns:
            prop_row = properties_df[properties_df['smiles'] == row['smiles']]
            if not prop_row.empty:
                node['properties'].update({
                    'MW': prop_row.iloc[0].get('MW'),
                    'LogP': prop_row.iloc[0].get('LogP'),
                    'QED': prop_row.iloc[0].get('QED'),
                    'is_drug_like': prop_row.iloc[0].get('drug_like'),
                })
        
        nodes.append(node)
    
    return nodes


def create_protein_nodes(proteins_df):
    """Create protein nodes for KG"""
    nodes = []
    
    for idx, row in proteins_df.iterrows():
        node = {
            'id': f"prot_{row['protein_id']}",
            'type': 'Protein',
            'properties': {
                'name': row['name'],
                'target_type': row['target_type'],
                'source': row['source'],
                'source_id': row['protein_id'],
            }
        }
        nodes.append(node)
    
    return nodes


def create_bioactivity_edges(bioactivity_df):
    """Create molecule-TARGETS-protein edges"""
    edges = []
    
    for idx, row in bioactivity_df.iterrows():
        edge = {
            'source': f"mol_{row['molecule_id']}",
            'target': f"prot_{row['protein_id']}",
            'type': 'TARGETS',
            'properties': {
                'pic50': row['pic50'],
                'activity_type': row['standard_type'],
                'activity_value': row['standard_value'],
                'assay_type': row['assay_type'],
                'source': row['source'],
                'confidence': 0.85,
            }
        }
        edges.append(edge)
    
    return edges


def create_similarity_edges(molecules_df):
    """Create molecule-SIMILAR_TO-molecule edges"""
    from rdkit import Chem
    from rdkit.Chem import AllChem
    
    edges = []
    smiles_list = molecules_df['smiles'].tolist()
    
    # Compute fingerprints
    fps = []
    for smiles in smiles_list:
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
            fps.append(fp)
        else:
            fps.append(None)
    
    # Find similar pairs (simplified: sample pairs)
    import numpy as np
    sampled_indices = np.random.choice(len(smiles_list), min(500, len(smiles_list)), replace=False)
    
    for i in sampled_indices:
        if fps[i] is None:
            continue
        
        # Find top-k similar
        similarities = []
        for j, fp_j in enumerate(fps):
            if i != j and fp_j is not None:
                similarity = Chem.DataStructs.TanimotoSimilarity(fps[i], fp_j)
                if similarity > 0.7:
                    similarities.append((j, similarity))
        
        for j, sim in sorted(similarities, key=lambda x: x[1], reverse=True)[:3]:
            edge = {
                'source': f"mol_{molecules_df.iloc[i]['molecule_id']}",
                'target': f"mol_{molecules_df.iloc[j]['molecule_id']}",
                'type': 'SIMILAR_TO',
                'properties': {
                    'similarity': sim,
                    'method': 'Morgan_Tanimoto',
                }
            }
            edges.append(edge)
    
    return edges


def save_nodes_to_csv(nodes, output_file):
    """Save nodes to CSV for graph loading"""
    data = []
    for node in nodes:
        row = {'id': node['id'], 'type': node['type']}
        row.update(node['properties'])
        data.append(row)
    
    df = pd.DataFrame(data)
    df.to_csv(output_file, index=False)
    print(f"Saved {len(nodes)} nodes to {output_file}")


def save_edges_to_csv(edges, output_file):
    """Save edges to CSV for graph loading"""
    data = []
    for edge in edges:
        row = {
            'source': edge['source'],
            'target': edge['target'],
            'type': edge['type'],
        }
        row.update(edge['properties'])
        data.append(row)
    
    df = pd.DataFrame(data)
    df.to_csv(output_file, index=False)
    print(f"Saved {len(edges)} edges to {output_file}")


if __name__ == "__main__":
    # Extract data
    molecules_df = extract_molecules_from_chembl('chembl_36/chembl_36_sqlite/chembl_36.db', limit=50000)
    proteins_df = extract_targets_from_chembl('chembl_36/chembl_36_sqlite/chembl_36.db')
    bioactivity_df = extract_bioactivity_from_chembl('chembl_36/chembl_36_sqlite/chembl_36.db')
    
    # Create nodes
    mol_nodes = create_molecule_nodes(molecules_df)
    prot_nodes = create_protein_nodes(proteins_df)
    
    # Create edges
    bio_edges = create_bioactivity_edges(bioactivity_df)
    sim_edges = create_similarity_edges(molecules_df)
    
    # Save
    Path('d:\\ChemAI\\TEAM_SHREYA\\data').mkdir(parents=True, exist_ok=True)
    save_nodes_to_csv(mol_nodes, 'd:\\ChemAI\\TEAM_SHREYA\\data\\molecules_kg.csv')
    save_nodes_to_csv(prot_nodes, 'd:\\ChemAI\\TEAM_SHREYA\\data\\proteins_kg.csv')
    save_edges_to_csv(bio_edges, 'd:\\ChemAI\\TEAM_SHREYA\\data\\bioactivity_edges.csv')
    save_edges_to_csv(sim_edges, 'd:\\ChemAI\\TEAM_SHREYA\\data\\similarity_edges.csv')
