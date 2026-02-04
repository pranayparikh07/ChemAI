"""
ChEMBL to Neo4j Graph Loader
Loads molecular data from ChEMBL SQLite into Neo4j Knowledge Graph
"""

import sqlite3
import pandas as pd
from typing import List, Dict
import numpy as np

try:
    from neo4j import GraphDatabase
    NEO4J_AVAILABLE = True
except ImportError:
    NEO4J_AVAILABLE = False
    print("Warning: neo4j package not installed. Install with: pip install neo4j")


class ChEMBLGraphLoader:
    """Loads ChEMBL data into Neo4j graph database"""
    
    def __init__(self, chembl_db_path: str = "chembl_36/chembl_36_sqlite/chembl_36.db"):
        """
        Initialize loader
        
        Args:
            chembl_db_path: Path to ChEMBL SQLite database
        """
        self.chembl_db_path = chembl_db_path
        self.conn = sqlite3.connect(chembl_db_path)
        self.driver = None
        self.stats = {
            'molecules': 0,
            'proteins': 0,
            'interactions': 0,
            'errors': 0
        }
    
    def connect_neo4j(self, uri: str, user: str, password: str):
        """
        Connect to Neo4j instance
        
        Args:
            uri: Neo4j connection URI (e.g., 'neo4j://localhost:7687')
            user: Neo4j username
            password: Neo4j password
        """
        if not NEO4J_AVAILABLE:
            print("Error: neo4j package not installed")
            return False
        
        try:
            self.driver = GraphDatabase.driver(uri, auth=(user, password))
            self.driver.verify_connectivity()
            print(f"✓ Connected to Neo4j at {uri}")
            return True
        except Exception as e:
            print(f"✗ Failed to connect to Neo4j: {e}")
            return False
    
    def init_graph(self):
        """Initialize graph database with constraints and indexes"""
        if not self.driver:
            print("Error: Not connected to Neo4j")
            return False
        
        constraints = [
            "CREATE CONSTRAINT molecule_smiles IF NOT EXISTS ON (m:Molecule) ASSERT m.smiles IS UNIQUE",
            "CREATE CONSTRAINT protein_uniprot IF NOT EXISTS ON (p:Protein) ASSERT p.uniprot_id IS UNIQUE",
            "CREATE CONSTRAINT disease_id IF NOT EXISTS ON (d:Disease) ASSERT d.disease_id IS UNIQUE",
            "CREATE CONSTRAINT gene_id IF NOT EXISTS ON (g:Gene) ASSERT g.gene_id IS UNIQUE",
            "CREATE CONSTRAINT pathway_id IF NOT EXISTS ON (p:Pathway) ASSERT p.pathway_id IS UNIQUE",
        ]
        
        indexes = [
            "CREATE INDEX molecule_chembl IF NOT EXISTS FOR (m:Molecule) ON (m.chembl_id)",
            "CREATE INDEX protein_chembl IF NOT EXISTS FOR (p:Protein) ON (p.chembl_id)",
            "CREATE INDEX protein_gene IF NOT EXISTS FOR (p:Protein) ON (p.gene_name)",
            "CREATE INDEX disease_mesh IF NOT EXISTS FOR (d:Disease) ON (d.mesh_id)",
            "CREATE INDEX gene_name IF NOT EXISTS FOR (g:Gene) ON (g.gene_name)",
        ]
        
        with self.driver.session() as session:
            for constraint in constraints:
                try:
                    session.run(constraint)
                except Exception as e:
                    pass  # Constraint already exists
            
            for index in indexes:
                try:
                    session.run(index)
                except Exception as e:
                    pass  # Index already exists
        
        print("✓ Graph initialized with constraints and indexes")
        return True
    
    def load_molecules(self, limit: int = 50000) -> int:
        """
        Load molecules from ChEMBL to Neo4j
        
        Args:
            limit: Maximum number of molecules to load
            
        Returns:
            Number of molecules loaded
        """
        if not self.driver:
            print("Error: Not connected to Neo4j")
            return 0
        
        print(f"Loading molecules (limit: {limit})...")
        
        query = """
        SELECT cd.chembl_id, cd.pref_name, cs.canonical_smiles, 
               cd.max_phase
        FROM compound_structures cs
        JOIN molecule_dictionary cd ON cs.molregno = cd.molregno
        WHERE cs.canonical_smiles IS NOT NULL
        LIMIT ?
        """
        
        df = pd.read_sql_query(query, self.conn, params=(limit,))
        
        batch_size = 1000
        loaded = 0
        
        with self.driver.session() as session:
            for i in range(0, len(df), batch_size):
                batch = df.iloc[i:i+batch_size]
                
                for _, row in batch.iterrows():
                    try:
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
                            name=row['pref_name'] if pd.notna(row['pref_name']) else 'Unknown',
                            max_phase=int(row['max_phase']) if pd.notna(row['max_phase']) else 0
                        )
                        loaded += 1
                    except Exception as e:
                        self.stats['errors'] += 1
                
                if (i // batch_size + 1) % 10 == 0:
                    print(f"  Loaded {min(i + batch_size, len(df))}/{len(df)} molecules...")
        
        self.stats['molecules'] = loaded
        print(f"✓ Loaded {loaded} molecules")
        return loaded
    
    def load_proteins(self, limit: int = 10000) -> int:
        """
        Load protein targets from ChEMBL to Neo4j
        
        Args:
            limit: Maximum number of proteins to load
            
        Returns:
            Number of proteins loaded
        """
        if not self.driver:
            print("Error: Not connected to Neo4j")
            return 0
        
        print(f"Loading proteins (limit: {limit})...")
        
        query = """
        SELECT DISTINCT td.chembl_id, td.pref_name, td.target_type
        FROM target_dictionary td
        WHERE td.pref_name IS NOT NULL
        LIMIT ?
        """
        
        df = pd.read_sql_query(query, self.conn, params=(limit,))
        
        batch_size = 500
        loaded = 0
        
        with self.driver.session() as session:
            for i in range(0, len(df), batch_size):
                batch = df.iloc[i:i+batch_size]
                
                for _, row in batch.iterrows():
                    try:
                        session.run(
                            """
                            CREATE (p:Protein {
                                chembl_id: $chembl_id,
                                name: $name,
                                protein_class: $protein_class,
                                uniprot_id: $uniprot_id,
                                gene_name: $gene_name
                            })
                            """,
                            chembl_id=row['chembl_id'],
                            name=row['pref_name'],
                            protein_class=row['target_type'] if pd.notna(row['target_type']) else 'Unknown',
                            uniprot_id=row['chembl_id'],  # Placeholder
                            gene_name=row['pref_name'].split()[0] if pd.notna(row['pref_name']) else 'Unknown'
                        )
                        loaded += 1
                    except Exception as e:
                        self.stats['errors'] += 1
                
                if (i // batch_size + 1) % 5 == 0:
                    print(f"  Loaded {min(i + batch_size, len(df))}/{len(df)} proteins...")
        
        self.stats['proteins'] = loaded
        print(f"✓ Loaded {loaded} proteins")
        return loaded
    
    def load_interactions(self, limit: int = 50000) -> int:
        """
        Load drug-target interactions from ChEMBL to Neo4j
        
        Args:
            limit: Maximum number of interactions to load
            
        Returns:
            Number of interactions loaded
        """
        if not self.driver:
            print("Error: Not connected to Neo4j")
            return 0
        
        print(f"Loading interactions (limit: {limit})...")
        
        query = """
        SELECT mc.chembl_id as molecule_id, mc.pref_name as mol_name,
               td.chembl_id as target_id, td.pref_name as target_name,
               a.standard_value, a.standard_units, a.standard_type,
               cs.canonical_smiles
        FROM activities a
        JOIN molecule_dictionary mc ON a.molregno = mc.molregno
        JOIN target_dictionary td ON a.tid = td.tid
        JOIN compound_structures cs ON a.molregno = cs.molregno
        WHERE a.standard_value IS NOT NULL 
        AND a.standard_type = 'IC50'
        AND a.standard_units = 'nM'
        AND cs.canonical_smiles IS NOT NULL
        LIMIT ?
        """
        
        df = pd.read_sql_query(query, self.conn, params=(limit,))
        
        batch_size = 1000
        loaded = 0
        
        with self.driver.session() as session:
            for i in range(0, len(df), batch_size):
                batch = df.iloc[i:i+batch_size]
                
                for _, row in batch.iterrows():
                    try:
                        # Compute pIC50
                        pIC50 = 9 - np.log10(max(float(row['standard_value']), 0.1))
                        
                        session.run(
                            """
                            MATCH (m:Molecule {smiles: $smiles})
                            MATCH (p:Protein {chembl_id: $target_id})
                            CREATE (m)-[r:TARGETS {
                                ic50_nM: $ic50,
                                pIC50: $pIC50,
                                type: $type,
                                confidence_score: 0.8
                            }]->(p)
                            """,
                            smiles=row['canonical_smiles'],
                            target_id=row['target_id'],
                            ic50=float(row['standard_value']),
                            pIC50=pIC50,
                            type=row['standard_type']
                        )
                        loaded += 1
                    except Exception as e:
                        self.stats['errors'] += 1
                
                if (i // batch_size + 1) % 10 == 0:
                    print(f"  Loaded {min(i + batch_size, len(df))}/{len(df)} interactions...")
        
        self.stats['interactions'] = loaded
        print(f"✓ Loaded {loaded} interactions")
        return loaded
    
    def create_similarity_edges(self, similarity_threshold: float = 0.7) -> int:
        """
        Create SIMILAR_TO edges between molecules
        
        Args:
            similarity_threshold: Minimum Tanimoto similarity score
            
        Returns:
            Number of similarity edges created
        """
        if not self.driver:
            print("Error: Not connected to Neo4j")
            return 0
        
        print(f"Creating similarity edges (threshold: {similarity_threshold})...")
        
        # Get all molecules with fingerprints
        query = """
        SELECT m1.smiles as smiles1, m2.smiles as smiles2
        FROM (
            SELECT DISTINCT smiles FROM Molecule LIMIT 1000
        ) as m1
        CROSS JOIN (
            SELECT DISTINCT smiles FROM Molecule LIMIT 1000
        ) as m2
        WHERE m1.smiles < m2.smiles
        """
        
        # This is simplified - in practice, compute Tanimoto similarity using RDKit
        try:
            from rdkit import Chem
            from rdkit.Chem import AllChem
            from rdkit import DataStructs
            
            # Get molecules from graph
            with self.driver.session() as session:
                result = session.run("""
                    MATCH (m:Molecule) 
                    RETURN m.smiles as smiles 
                    LIMIT 100
                """)
                
                smiles_list = [record['smiles'] for record in result]
                mols = [Chem.MolFromSmiles(s) for s in smiles_list if s]
                fps = [AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=2048) for m in mols if m]
                
                # Create edges for similar molecules
                created = 0
                for i in range(len(fps)):
                    for j in range(i+1, len(fps)):
                        sim = DataStructs.TanimotoSimilarity(fps[i], fps[j])
                        if sim >= similarity_threshold:
                            try:
                                session.run("""
                                    MATCH (m1:Molecule {smiles: $smiles1})
                                    MATCH (m2:Molecule {smiles: $smiles2})
                                    CREATE (m1)-[r:SIMILAR_TO {
                                        similarity_score: $sim,
                                        fingerprint_type: 'Morgan'
                                    }]->(m2)
                                """,
                                    smiles1=smiles_list[i],
                                    smiles2=smiles_list[j],
                                    sim=float(sim)
                                )
                                created += 1
                            except:
                                pass
                
                print(f"✓ Created {created} similarity edges")
                return created
        
        except ImportError:
            print("⚠ RDKit not available for similarity computation")
            return 0
    
    def get_statistics(self) -> Dict:
        """Get loading statistics"""
        return self.stats
    
    def close(self):
        """Close database connections"""
        if self.conn:
            self.conn.close()
        if self.driver:
            self.driver.close()
    
    def run_full_load(self, neo4j_uri: str, neo4j_user: str, neo4j_password: str):
        """Run complete loading pipeline"""
        print("="*70)
        print("ChEMBL to Neo4j Full Loading Pipeline")
        print("="*70)
        
        # Connect to Neo4j
        if not self.connect_neo4j(neo4j_uri, neo4j_user, neo4j_password):
            print("Failed to connect to Neo4j")
            return False
        
        # Initialize graph
        self.init_graph()
        
        # Load data
        self.load_molecules(50000)
        self.load_proteins(10000)
        self.load_interactions(50000)
        self.create_similarity_edges(0.7)
        
        # Print statistics
        print("\n" + "="*70)
        print("Loading Complete")
        print("="*70)
        print(f"Molecules loaded: {self.stats['molecules']}")
        print(f"Proteins loaded: {self.stats['proteins']}")
        print(f"Interactions loaded: {self.stats['interactions']}")
        print(f"Errors: {self.stats['errors']}")
        
        return True


if __name__ == "__main__":
    # Example usage
    print("ChEMBL Graph Loader")
    print("Usage: python graph_loader.py")
    print("\nTo load ChEMBL data into Neo4j:")
    print("  loader = ChEMBLGraphLoader()")
    print("  loader.run_full_load('neo4j://localhost:7687', 'neo4j', 'password')")
