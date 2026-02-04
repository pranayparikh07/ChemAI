import sqlite3
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski
import os

def scan_and_clean_database(db_path, output_file="clean_candidates.csv", limit=100000):
    print(f"Scanning database: {db_path}")
    print(f"Loading up to {limit} molecules...")
    
    # Query to fetch molecules and properties
    # We prioritize those with high QED (drug-likeness)
    query = f"""
    SELECT 
        cs.molregno,
        cs.canonical_smiles,
        cp.mw_freebase,
        cp.alogp,
        cp.hba,
        cp.hbd,
        cp.psa,
        cp.qed_weighted,
        cp.num_ro5_violations
    FROM compound_structures cs
    JOIN compound_properties cp ON cs.molregno = cp.molregno
    WHERE cs.canonical_smiles IS NOT NULL
      AND cp.mw_freebase IS NOT NULL
    ORDER BY cp.qed_weighted DESC
    LIMIT {limit}
    """
    
    conn = sqlite3.connect(db_path)
    try:
        df = pd.read_sql_query(query, conn)
    finally:
        conn.close()
        
    print(f"Loaded {len(df)} raw records.")
    
    # Cleaning steps
    print("Cleaning data...")
    
    # 1. Drop duplicates
    df = df.drop_duplicates(subset=['canonical_smiles'])
    print(f"After removing duplicates: {len(df)}")
    
    # 2. Filter for Drug-Likeness (Lipinski's Rule of 5)
    # MW <= 500, LogP <= 5, HBD <= 5, HBA <= 10
    # Note: ChEMBL already has num_ro5_violations, we can use that or being strict.
    # Let's be strict but allow small deviations for 'candidates'
    
    df_clean = df[
        (df['mw_freebase'] <= 500) &
        (df['alogp'] <= 5) &
        (df['hbd'] <= 5) &
        (df['hba'] <= 10)
    ].copy()
    
    print(f"After Lipinski filters: {len(df_clean)}")
    
    # 3. Validate SMILES (ensure RDKit can parse them)
    valid_indices = []
    for idx, row in df_clean.iterrows():
        if Chem.MolFromSmiles(row['canonical_smiles']):
            valid_indices.append(idx)
            
    df_final = df_clean.loc[valid_indices].copy()
    print(f"After RDKit validation: {len(df_final)}")
    
    # Rank by QED (already sorted by query, but re-sort to be sure)
    df_final = df_final.sort_values(by='qed_weighted', ascending=False)
    
    # Save
    print(f"Saving top candidates to {output_file}...")
    df_final.to_csv(output_file, index=False)
    print("Done.")
    
    # Show top 5
    print("\nTop 5 Drug-Like Candidates:")
    print(df_final[['canonical_smiles', 'qed_weighted', 'mw_freebase']].head(5))

if __name__ == "__main__":
    db_path = "chembl_36/chembl_36_sqlite/chembl_36.db"
    scan_and_clean_database(db_path)
