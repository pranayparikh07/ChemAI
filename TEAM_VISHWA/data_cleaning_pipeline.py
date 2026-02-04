#!/usr/bin/env python3
\"\"\"
Data Cleaning Pipeline for ChemAI
Removes duplicates, invalid SMILES, and standardizes molecular data
\"\"\"

import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
import sqlite3


def load_chembl_data(db_path, limit=50000):
    \"\"\"Load SMILES data from ChEMBL database\"\"\"
    query = """
    SELECT DISTINCT 
        molregno,
        canonical_smiles as smiles
    FROM compound_structures
    WHERE canonical_smiles IS NOT NULL
    LIMIT ?
    """
    conn = sqlite3.connect(db_path)
    df = pd.read_sql_query(query, conn, params=(limit,))
    conn.close()
    print(f"Loaded {len(df)} molecules from ChEMBL")
    return df


def validate_smiles(smiles):
    \"\"\"Check if SMILES string is valid\"\"\"
    try:
        mol = Chem.MolFromSmiles(smiles)
        return mol is not None
    except:
        return False


def canonicalize_smiles(smiles):
    \"\"\"Canonicalize SMILES string\"\"\"
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        return Chem.MolToSmiles(mol)
    except:
        return None


def remove_salts_and_solvents(smiles):
    \"\"\"Remove salt and solvent molecules, keep largest fragment\"\"\"
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        
        # Split on disconnected components
        parts = Chem.MolToSmiles(mol).split('.')
        if len(parts) > 1:
            # Keep largest fragment
            largest = max(parts, key=len)
            return largest
        return Chem.MolToSmiles(mol)
    except:
        return None


def filter_by_molecular_weight(smiles, min_mw=100, max_mw=800):
    \"\"\"Keep only molecules within MW range\"\"\"
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False
        mw = Descriptors.MolWt(mol)
        return min_mw <= mw <= max_mw
    except:
        return False


def filter_by_element_composition(smiles):
    \"\"\"Keep only organic molecules with common elements\"\"\"
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False
        
        # Check that all atoms are organic
        allowed = {'C', 'H', 'N', 'O', 'S', 'F', 'Cl', 'Br', 'I', 'P'}
        for atom in mol.GetAtoms():
            if atom.GetSymbol() not in allowed:
                return False
        return True
    except:
        return False


def clean_dataset(input_df):
    \"\"\"Apply all cleaning steps\"\"\"
    df = input_df.copy()
    
    print(f"Starting: {len(df)} molecules")
    
    # Step 1: Remove duplicates
    initial_count = len(df)
    df = df.drop_duplicates(subset=['smiles'], keep='first')
    print(f"After removing duplicates: {len(df)} molecules (-{initial_count - len(df)})")
    
    # Step 2: Validate and canonicalize SMILES
    initial_count = len(df)
    df['valid'] = df['smiles'].apply(validate_smiles)
    df = df[df['valid']].copy()
    df['canonical_smiles'] = df['smiles'].apply(canonicalize_smiles)
    print(f"After SMILES validation: {len(df)} molecules (-{initial_count - len(df)})")
    
    # Step 3: Remove salts/solvents
    initial_count = len(df)
    df['processed_smiles'] = df['canonical_smiles'].apply(remove_salts_and_solvents)
    df = df[df['processed_smiles'].notna()].copy()
    print(f"After salt/solvent removal: {len(df)} molecules (-{initial_count - len(df)})")
    
    # Step 4: Filter by molecular weight
    initial_count = len(df)
    df = df[df['processed_smiles'].apply(filter_by_molecular_weight)].copy()
    print(f"After MW filter (100-800): {len(df)} molecules (-{initial_count - len(df)})")
    
    # Step 5: Filter by element composition
    initial_count = len(df)
    df = df[df['processed_smiles'].apply(filter_by_element_composition)].copy()
    print(f"After element filter: {len(df)} molecules (-{initial_count - len(df)})")
    
    # Step 6: Remove remaining duplicates
    initial_count = len(df)
    df = df.drop_duplicates(subset=['processed_smiles'], keep='first')
    print(f"After final duplicate removal: {len(df)} molecules (-{initial_count - len(df)})")
    
    # Cleanup
    df = df.drop(columns=['smiles', 'valid', 'canonical_smiles'])
    df = df.rename(columns={'processed_smiles': 'smiles'})
    
    return df


def save_cleaned_data(df, output_file):
    \"\"\"Save cleaned data to CSV\"\"\"
    df.to_csv(output_file, index=False)
    print(f"Saved {len(df)} molecules to {output_file}")


if __name__ == "__main__":
    # Load data
    df = load_chembl_data('chembl_36/chembl_36_sqlite/chembl_36.db', limit=100000)
    
    # Clean data
    df_clean = clean_dataset(df)
    
    # Save results
    save_cleaned_data(df_clean, 'd:\\\\ChemAI\\\\TEAM_VISHWA\\\\data\\\\clean_molecules.csv')
    
    print(f"\n✓ Cleaning complete!")
    print(f"Final dataset: {len(df_clean)} molecules")
