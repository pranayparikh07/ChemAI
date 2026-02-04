#!/usr/bin/env python3
"""
Load Knowledge Graph to Remote Neo4j
Connects to Neo4j instance and ingests all graph data
"""

from neo4j import GraphDatabase, basic_auth
import pandas as pd
from pathlib import Path
import os


class Neo4jLoader:
    def __init__(self, uri, username, password):
        """Initialize Neo4j connection"""
        self.driver = GraphDatabase.driver(uri, auth=basic_auth(username, password))
        self.session = self.driver.session()
        print(f"✓ Connected to Neo4j: {uri}")
    
    def close(self):
        """Close connection"""
        self.session.close()
        self.driver.close()
    
    def clear_database(self):
        """Clear all data (CAUTION!)"""
        print("\n⚠️  Clearing database...")
        self.session.run("MATCH (n) DETACH DELETE n")
        print("✓ Database cleared")
    
    def create_constraints(self):
        """Create indexes and constraints"""
        print("\n📊 Creating constraints and indexes...")
        
        constraints = [
            "CREATE CONSTRAINT molecule_id IF NOT EXISTS FOR (m:Molecule) REQUIRE m.id IS UNIQUE",
            "CREATE CONSTRAINT protein_id IF NOT EXISTS FOR (p:Protein) REQUIRE p.id IS UNIQUE",
            "CREATE CONSTRAINT disease_id IF NOT EXISTS FOR (d:Disease) REQUIRE d.id IS UNIQUE",
            "CREATE INDEX molecule_smiles IF NOT EXISTS FOR (m:Molecule) ON (m.smiles)",
            "CREATE INDEX protein_name IF NOT EXISTS FOR (p:Protein) ON (p.name)",
            "CREATE INDEX molecule_qed IF NOT EXISTS FOR (m:Molecule) ON (m.qed)",
        ]
        
        for constraint in constraints:
            try:
                self.session.run(constraint)
            except:
                pass  # Constraint may already exist
        
        print("✓ Constraints and indexes created")
    
    def load_molecules(self, csv_path):
        """Load molecules from CSV"""
        print(f"\n📌 Loading molecules from {csv_path}...")
        
        if not Path(csv_path).exists():
            print(f"⚠️  File not found: {csv_path}")
            print("Creating sample molecules instead...")
            self._load_sample_molecules()
            return
        
        df = pd.read_csv(csv_path)
        loaded = 0
        skipped = 0
        
        for idx, row in df.iterrows():
            query = """
            MERGE (m:Molecule {id: $id})
            ON CREATE SET
                m.smiles = $smiles,
                m.name = $name,
                m.mw = $mw,
                m.logp = $logp,
                m.qed = $qed,
                m.is_drug_like = $is_drug_like,
                m.source = $source,
                m.confidence = $confidence
            ON MATCH SET
                m.mw = COALESCE(m.mw, $mw),
                m.logp = COALESCE(m.logp, $logp),
                m.qed = COALESCE(m.qed, $qed)
            RETURN m.id as created_id
            """
            
            try:
                result = self.session.run(query, {
                    'id': row.get('id', f'mol_{idx}'),
                    'smiles': row.get('smiles', ''),
                    'name': row.get('name', ''),
                    'mw': float(row.get('MW', 0)),
                    'logp': float(row.get('LogP', 0)),
                    'qed': float(row.get('QED', 0)),
                    'is_drug_like': row.get('is_drug_like', True),
                    'source': row.get('source', 'ChEMBL'),
                    'confidence': float(row.get('confidence', 0.9))
                })
                loaded += 1
            except Exception as e:
                skipped += 1
                if skipped <= 5:  # Show first 5 errors
                    print(f"  ⚠️  Skipped row {idx}: {str(e)[:50]}")
            
            if (idx + 1) % 1000 == 0:
                print(f"  ✓ Processed {idx + 1} rows (Loaded: {loaded}, Skipped: {skipped})")
        
        print(f"✓ Total molecules loaded: {loaded}, Skipped: {skipped}")
    
    def _load_sample_molecules(self):
        """Load sample molecules for demo"""
        molecules = [
            ('mol_001', 'CC(=O)Oc1ccccc1C(=O)O', 'Aspirin', 180.16, 1.19, 0.85, True),
            ('mol_002', 'CC(C)Cc1ccc(cc1)C(C)C(=O)O', 'Ibuprofen', 206.28, 3.97, 0.83, True),
            ('mol_003', 'COc1ccc2cc(ccc2c1)C(C)C(=O)O', 'Naproxen', 230.26, 3.18, 0.81, True),
            ('mol_004', 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C', 'Caffeine', 194.19, -0.07, 0.72, True),
        ]
        
        for mol_id, smiles, name, mw, logp, qed, is_drug_like in molecules:
            query = """
            CREATE (m:Molecule {
                id: $id,
                smiles: $smiles,
                name: $name,
                mw: $mw,
                logp: $logp,
                qed: $qed,
                is_drug_like: $is_drug_like,
                source: 'Sample'
            })
            """
            
            self.session.run(query, {
                'id': mol_id,
                'smiles': smiles,
                'name': name,
                'mw': mw,
                'logp': logp,
                'qed': qed,
                'is_drug_like': is_drug_like
            })
        
        print(f"✓ Sample molecules loaded: {len(molecules)}")
    
    def load_proteins(self, csv_path):
        """Load proteins from CSV"""
        print(f"\n🧬 Loading proteins from {csv_path}...")
        
        if not Path(csv_path).exists():
            print(f"⚠️  File not found: {csv_path}")
            print("Creating sample proteins instead...")
            self._load_sample_proteins()
            return
        
        df = pd.read_csv(csv_path)
        loaded = 0
        skipped = 0
        
        for idx, row in df.iterrows():
            query = """
            MERGE (p:Protein {id: $id})
            ON CREATE SET
                p.uniprot = $uniprot,
                p.name = $name,
                p.gene = $gene,
                p.organism = $organism,
                p.function = $function,
                p.source = $source,
                p.confidence = $confidence
            ON MATCH SET
                p.uniprot = COALESCE(p.uniprot, $uniprot),
                p.name = COALESCE(p.name, $name),
                p.gene = COALESCE(p.gene, $gene)
            RETURN p.id as created_id
            """
            
            try:
                self.session.run(query, {
                    'id': row.get('id', f'prot_{idx}'),
                    'uniprot': row.get('uniprot', ''),
                    'name': row.get('name', ''),
                    'gene': row.get('gene', ''),
                    'organism': row.get('organism', 'Homo sapiens'),
                    'function': row.get('function', ''),
                    'source': row.get('source', 'ChEMBL'),
                    'confidence': float(row.get('confidence', 0.9))
                })
                loaded += 1
            except Exception as e:
                skipped += 1
                if skipped <= 5:
                    print(f"  ⚠️  Skipped row {idx}: {str(e)[:50]}")
            
            if (idx + 1) % 500 == 0:
                print(f"  ✓ Processed {idx + 1} rows (Loaded: {loaded}, Skipped: {skipped})")
        
        print(f"✓ Total proteins loaded: {loaded}, Skipped: {skipped}")
    
    def _load_sample_proteins(self):
        """Load sample proteins for demo"""
        proteins = [
            ('prot_001', 'P23458', 'Cyclooxygenase-2', 'PTGS2', 'Homo sapiens', 'Prostaglandin synthesis'),
            ('prot_002', 'P23457', 'Cyclooxygenase-1', 'PTGS1', 'Homo sapiens', 'Prostaglandin synthesis'),
            ('prot_003', 'P08588', 'Adenosine receptor A2a', 'ADORA2A', 'Homo sapiens', 'G-protein coupled receptor'),
        ]
        
        for prot_id, uniprot, name, gene, organism, function in proteins:
            query = """
            CREATE (p:Protein {
                id: $id,
                uniprot: $uniprot,
                name: $name,
                gene: $gene,
                organism: $organism,
                function: $function,
                source: 'Sample'
            })
            """
            
            self.session.run(query, {
                'id': prot_id,
                'uniprot': uniprot,
                'name': name,
                'gene': gene,
                'organism': organism,
                'function': function
            })
        
        print(f"✓ Sample proteins loaded: {len(proteins)}")
    
    def load_bioactivity_edges(self, csv_path):
        """Load bioactivity relationships"""
        print(f"\n🔗 Loading bioactivity edges from {csv_path}...")
        
        if not Path(csv_path).exists():
            print(f"⚠️  File not found: {csv_path}")
            print("Creating sample bioactivity edges instead...")
            self._load_sample_bioactivity()
            return
        
        df = pd.read_csv(csv_path)
        loaded = 0
        skipped = 0
        
        for idx, row in df.iterrows():
            query = """
            MATCH (m:Molecule {id: $mol_id})
            MATCH (p:Protein {id: $prot_id})
            MERGE (m)-[r:TARGETS {mol_id: $mol_id, prot_id: $prot_id}]->(p)
            ON CREATE SET
                r.pic50 = $pic50,
                r.activity_type = $activity_type,
                r.confidence = $confidence,
                r.source = $source
            ON MATCH SET
                r.pic50 = COALESCE(r.pic50, $pic50),
                r.confidence = COALESCE(r.confidence, $confidence)
            RETURN r
            """
            
            try:
                self.session.run(query, {
                    'mol_id': row.get('source', ''),
                    'prot_id': row.get('target', ''),
                    'pic50': float(row.get('pic50', 0)),
                    'activity_type': row.get('activity_type', 'IC50'),
                    'confidence': float(row.get('confidence', 0.9)),
                    'source': row.get('source', 'ChEMBL')
                })
                loaded += 1
            except Exception as e:
                skipped += 1
                if skipped <= 5:
                    print(f"  ⚠️  Skipped row {idx}: {str(e)[:50]}")
            
            if (idx + 1) % 5000 == 0:
                print(f"  ✓ Processed {idx + 1} rows (Loaded: {loaded}, Skipped: {skipped})")
        
        print(f"✓ Total bioactivity edges loaded: {loaded}, Skipped: {skipped}")
    
    def _load_sample_bioactivity(self):
        """Load sample bioactivity edges for demo"""
        edges = [
            ('mol_001', 'prot_001', 5.2, 0.92),  # Aspirin - COX-1
            ('mol_001', 'prot_002', 5.8, 0.95),  # Aspirin - COX-2
            ('mol_002', 'prot_001', 5.1, 0.90),  # Ibuprofen - COX-1
            ('mol_002', 'prot_002', 5.9, 0.93),  # Ibuprofen - COX-2
            ('mol_003', 'prot_001', 5.0, 0.88),  # Naproxen - COX-1
            ('mol_003', 'prot_002', 5.7, 0.91),  # Naproxen - COX-2
            ('mol_004', 'prot_003', 4.8, 0.87),  # Caffeine - Adenosine A2a
        ]
        
        for mol_id, prot_id, pic50, confidence in edges:
            query = """
            MATCH (m:Molecule {id: $mol_id})
            MATCH (p:Protein {id: $prot_id})
            CREATE (m)-[r:TARGETS {
                pic50: $pic50,
                activity_type: 'IC50',
                confidence: $confidence,
                source: 'Sample'
            }]->(p)
            """
            
            try:
                self.session.run(query, {
                    'mol_id': mol_id,
                    'prot_id': prot_id,
                    'pic50': pic50,
                    'confidence': confidence
                })
            except:
                pass
        
        print(f"✓ Sample bioactivity edges loaded: {len(edges)}")
    
    def load_similarity_edges(self, csv_path):
        """Load similarity relationships"""
        print(f"\n🔄 Loading similarity edges from {csv_path}...")
        
        if not Path(csv_path).exists():
            print(f"⚠️  File not found: {csv_path}")
            print("Creating sample similarity edges instead...")
            self._load_sample_similarity()
            return
        
        df = pd.read_csv(csv_path)
        loaded = 0
        skipped = 0
        
        for idx, row in df.iterrows():
            query = """
            MATCH (m1:Molecule {id: $mol1_id})
            MATCH (m2:Molecule {id: $mol2_id})
            MERGE (m1)-[r:SIMILAR_TO {mol1_id: $mol1_id, mol2_id: $mol2_id}]->(m2)
            ON CREATE SET
                r.similarity = $similarity,
                r.method = $method,
                r.confidence = $confidence
            ON MATCH SET
                r.similarity = COALESCE(r.similarity, $similarity),
                r.confidence = COALESCE(r.confidence, $confidence)
            RETURN r
            """
            
            try:
                self.session.run(query, {
                    'mol1_id': row.get('source', ''),
                    'mol2_id': row.get('target', ''),
                    'similarity': float(row.get('similarity', 0)),
                    'method': row.get('method', 'Tanimoto'),
                    'confidence': float(row.get('confidence', 0.9))
                })
                loaded += 1
            except Exception as e:
                skipped += 1
                if skipped <= 5:
                    print(f"  ⚠️  Skipped row {idx}: {str(e)[:50]}")
            
            if (idx + 1) % 100 == 0:
                print(f"  ✓ Processed {idx + 1} rows (Loaded: {loaded}, Skipped: {skipped})")
        
        print(f"✓ Total similarity edges loaded: {loaded}, Skipped: {skipped}")
    
    def _load_sample_similarity(self):
        """Load sample similarity edges for demo"""
        edges = [
            ('mol_001', 'mol_002', 0.82),  # Aspirin - Ibuprofen
            ('mol_001', 'mol_003', 0.79),  # Aspirin - Naproxen
            ('mol_002', 'mol_003', 0.85),  # Ibuprofen - Naproxen
        ]
        
        for mol1_id, mol2_id, similarity in edges:
            query = """
            MATCH (m1:Molecule {id: $mol1_id})
            MATCH (m2:Molecule {id: $mol2_id})
            CREATE (m1)-[r:SIMILAR_TO {
                similarity: $similarity,
                method: 'Tanimoto',
                confidence: 0.95,
                source: 'Sample'
            }]->(m2)
            """
            
            try:
                self.session.run(query, {
                    'mol1_id': mol1_id,
                    'mol2_id': mol2_id,
                    'similarity': similarity
                })
            except:
                pass
        
        print(f"✓ Sample similarity edges loaded: {len(edges)}")
    
    def get_statistics(self):
        """Get graph statistics"""
        print("\n📊 GRAPH STATISTICS")
        print("-" * 60)
        
        nodes = self.session.run("MATCH (n) RETURN count(n) as count").single()
        edges = self.session.run("MATCH ()-[r]->() RETURN count(r) as count").single()
        
        mol_count = self.session.run("MATCH (m:Molecule) RETURN count(m) as count").single()
        prot_count = self.session.run("MATCH (p:Protein) RETURN count(p) as count").single()
        targets = self.session.run("MATCH ()-[r:TARGETS]->() RETURN count(r) as count").single()
        similar = self.session.run("MATCH ()-[r:SIMILAR_TO]->() RETURN count(r) as count").single()
        
        print(f"  Total Nodes: {nodes['count']}")
        print(f"    • Molecules: {mol_count['count']}")
        print(f"    • Proteins: {prot_count['count']}")
        print(f"\n  Total Edges: {edges['count']}")
        print(f"    • TARGETS: {targets['count']}")
        print(f"    • SIMILAR_TO: {similar['count']}")
        print("\n✓ Ready to browse at: http://localhost:7474/browser/")


def main():
    """Main ingestion function"""
    
    print("\n" + "="*60)
    print("  NEO4J KNOWLEDGE GRAPH INGESTION")
    print("="*60 + "\n")
    
    # CONNECTION DETAILS - UPDATE THESE
    uri = input("Enter Neo4j URI (default: bolt://localhost:7687): ").strip() or "bolt://localhost:7687"
    username = input("Enter username (default: neo4j): ").strip() or "neo4j"
    password = input("Enter password: ").strip()
    
    if not password:
        print("⚠️  Password required!")
        return
    
    try:
        # Connect to Neo4j
        loader = Neo4jLoader(uri, username, password)
        
        # Clear database option
        clear = input("\nClear existing data? (y/n): ").strip().lower()
        if clear == 'y':
            loader.clear_database()
        
        # Create constraints
        loader.create_constraints()
        
        # Load data - use absolute paths based on script location
        script_dir = Path(__file__).parent
        loader.load_molecules(str(script_dir / 'data' / 'molecules_kg.csv'))
        loader.load_proteins(str(script_dir / 'data' / 'proteins_kg.csv'))
        loader.load_bioactivity_edges(str(script_dir / 'data' / 'bioactivity_edges.csv'))
        loader.load_similarity_edges(str(script_dir / 'data' / 'similarity_edges.csv'))
        
        # Show statistics
        loader.get_statistics()
        
        # Close connection
        loader.close()
        
        print("\n" + "="*60)
        print("✓ INGESTION COMPLETE")
        print("="*60)
        print("\n📊 Open Neo4j Browser:")
        print(f"   http://localhost:7474/browser/\n")
        print("📝 Try these Cypher queries:")
        print("   MATCH (m:Molecule)-[r:TARGETS]->(p:Protein) RETURN m, r, p LIMIT 10")
        print("   MATCH (m:Molecule)-[r:SIMILAR_TO]->(m2:Molecule) RETURN m, r, m2")
        print("   MATCH (m:Molecule {name: 'Aspirin'}) RETURN m\n")
        
    except Exception as e:
        print(f"\n❌ Error: {e}")
        print("Make sure Neo4j is running and credentials are correct")


if __name__ == "__main__":
    main()
