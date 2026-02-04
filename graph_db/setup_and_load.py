#!/usr/bin/env python
"""
Complete Neo4j Setup and Data Loading Script
Guides through Neo4j setup and loads ChEMBL data into knowledge graph
"""
        
import sys
import time
from pathlib import Path


def print_header(text):
    """Print formatted header"""
    print("\n" + "="*70)
    print(f"  {text}")
    print("="*70)


def print_section(text):
    """Print formatted section"""
    print(f"\n► {text}")
    print("-" * 70)


def check_neo4j_running(uri, user, password):
    """Check if Neo4j instance is running"""
    try:
        from neo4j import GraphDatabase
        driver = GraphDatabase.driver(uri, auth=(user, password))
        driver.verify_connectivity()
        driver.close()
        return True
    except Exception as e:
        print(f"✗ Neo4j not accessible: {e}")
        return False


def setup_neo4j_instructions():
    """Display Neo4j setup instructions"""
    print_header("NEO4J SETUP REQUIRED")
    
    print_section("Option 1: Local Installation (Recommended for Development)")
    print("""
    1. Download Neo4j Desktop from https://neo4j.com/download/
    2. Install and launch Neo4j Desktop
    3. Create a new database:
       - Click "New" button
       - Enter name: "chemai"
       - Choose version: 5.x (latest LTS)
    4. Start the database
    5. Set a password (default user: neo4j)
    6. Note the connection details (usually localhost:7687)
    """)
    
    print_section("Option 2: Docker (Fastest Setup)")
    print("""
    1. Install Docker Desktop
    2. Run this command:
       
       docker run -p 7687:7687 -p 7474:7474 \\
         -e NEO4J_AUTH=neo4j/my_password \\
         neo4j:latest
    
    3. Wait for "Ready to accept connections"
    4. Connection: neo4j://localhost:7687
    """)
    
    print_section("Option 3: Neo4j AuraDB (Cloud - Free Tier)")
    print("""
    1. Go to https://neo4j.com/cloud/aura/
    2. Create free account
    3. Create instance (takes ~2 minutes)
    4. Copy connection URI and credentials
    """)
    
    print_section("Verify Neo4j is Running")
    print("""
    After starting Neo4j, verify connection with:
    
    python setup_and_load.py verify <uri> <user> <password>
    
    Example:
    python setup_and_load.py verify neo4j://localhost:7687 neo4j my_password
    """)


def verify_neo4j(uri, user, password):
    """Verify Neo4j connection"""
    print_header("VERIFYING NEO4J CONNECTION")
    
    print_section(f"Connecting to {uri}")
    
    try:
        from neo4j import GraphDatabase
        
        driver = GraphDatabase.driver(uri, auth=(user, password))
        driver.verify_connectivity()
        
        print(f"✓ Successfully connected to Neo4j!")
        print(f"  URI: {uri}")
        print(f"  User: {user}")
        
        # Get database info
        with driver.session() as session:
            result = session.run("CALL dbms.components() YIELD name, versions, edition")
            for record in result:
                print(f"  Neo4j {record['versions'][0]} ({record['edition']})")
        
        driver.close()
        return True
    
    except ImportError:
        print("✗ neo4j package not installed")
        print("  Run: pip install neo4j")
        return False
    except Exception as e:
        print(f"✗ Connection failed: {e}")
        print(f"\nTroubleshooting:")
        print(f"  1. Verify Neo4j is running")
        print(f"  2. Check URI format (e.g., neo4j://localhost:7687)")
        print(f"  3. Verify username and password")
        print(f"  4. Check firewall settings")
        return False


def load_data(uri, user, password):
    """Load ChEMBL data into Neo4j"""
    print_header("LOADING CHEMBL DATA INTO NEO4J")
    
    try:
        from graph_loader import ChEMBLGraphLoader
    except ImportError as e:
        print(f"✗ Failed to import ChEMBLGraphLoader: {e}")
        return False
    
    try:
        print_section("Initializing loader")
        loader = ChEMBLGraphLoader()
        print("✓ ChEMBLGraphLoader instantiated")
        
        print_section("Connecting to Neo4j")
        if not loader.connect_neo4j(uri, user, password):
            print("✗ Failed to connect to Neo4j")
            return False
        
        print_section("Initializing graph database")
        loader.init_graph()
        print("✓ Constraints and indexes created")
        
        print_section("Loading molecules (50K)")
        print("  This will take 5-10 minutes...")
        start = time.time()
        loader.load_molecules(50000)
        elapsed = time.time() - start
        print(f"✓ Molecules loaded in {elapsed:.1f} seconds")
        
        print_section("Loading proteins (10K)")
        print("  This will take 2-3 minutes...")
        start = time.time()
        loader.load_proteins(10000)
        elapsed = time.time() - start
        print(f"✓ Proteins loaded in {elapsed:.1f} seconds")
        
        print_section("Loading interactions (50K+)")
        print("  This will take 10-15 minutes...")
        start = time.time()
        loader.load_interactions(50000)
        elapsed = time.time() - start
        print(f"✓ Interactions loaded in {elapsed:.1f} seconds")
        
        print_section("Creating similarity edges")
        print("  This will take 5-10 minutes...")
        start = time.time()
        loader.create_similarity_edges(0.7)
        elapsed = time.time() - start
        print(f"✓ Similarity edges created in {elapsed:.1f} seconds")
        
        stats = loader.get_statistics()
        
        print_section("Data Loading Complete!")
        print(f"✓ Molecules loaded:      {stats['molecules']:,}")
        print(f"✓ Proteins loaded:       {stats['proteins']:,}")
        print(f"✓ Interactions loaded:   {stats['interactions']:,}")
        print(f"✓ Errors encountered:    {stats['errors']}")
        
        loader.close()
        return True
    
    except Exception as e:
        print(f"✗ Error during loading: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_algorithms(uri, user, password):
    """Test algorithms on loaded data"""
    print_header("TESTING GRAPH ALGORITHMS")
    
    try:
        from neo4j import GraphDatabase
        from graph_algorithms import GraphAlgorithms
        
        driver = GraphDatabase.driver(uri, auth=(user, password))
        algo = GraphAlgorithms(driver)
        
        print_section("Testing centrality ranking")
        print("  Computing protein degree centrality...")
        centrality = algo.compute_degree_centrality("Protein")
        if centrality:
            print(f"✓ Found {len(centrality)} proteins")
            for i, (node_id, name, degree) in enumerate(centrality[:5]):
                print(f"  {i+1}. {name}: degree={degree}")
        
        print_section("Testing similarity search")
        print("  Searching for similar compounds...")
        # Get a sample molecule from the graph
        with driver.session() as session:
            result = session.run("MATCH (m:Molecule) RETURN m.smiles LIMIT 1")
            records = list(result)
            if records:
                smiles = records[0]['m.smiles']
                similar = algo.find_similar_compounds(smiles, min_similarity=0.7, limit=5)
                print(f"✓ Found {len(similar)} similar compounds for {smiles[:30]}...")
                for i, comp in enumerate(similar[:3]):
                    print(f"  {i+1}. {comp.get('name', 'Unknown')}: similarity={comp.get('similarity', 0):.2f}")
        
        algo.close()
        return True
    
    except ImportError as e:
        print(f"⚠ Test skipped (missing dependencies): {e}")
        return True
    except Exception as e:
        print(f"⚠ Test error: {e}")
        return True


def main():
    """Main execution"""
    print_header("CHEMAI KNOWLEDGE GRAPH - NEO4J SETUP & DATA LOADING")
    
    if len(sys.argv) < 2:
        print_section("Usage")
        print("""
Examples:
  1. Setup instructions:
     python setup_and_load.py setup
  
  2. Verify Neo4j connection:
     python setup_and_load.py verify neo4j://localhost:7687 neo4j password
  
  3. Load data:
     python setup_and_load.py load neo4j://localhost:7687 neo4j password
  
  4. Run complete workflow (setup -> verify -> load -> test):
     python setup_and_load.py full neo4j://localhost:7687 neo4j password

Default (no args): Shows setup instructions
        """)
        setup_neo4j_instructions()
        return
    
    command = sys.argv[1].lower()
    
    if command == "setup":
        setup_neo4j_instructions()
    
    elif command == "verify":
        if len(sys.argv) < 5:
            print("✗ Missing arguments")
            print("  Usage: python setup_and_load.py verify <uri> <user> <password>")
            print("  Example: python setup_and_load.py verify neo4j://localhost:7687 neo4j password")
            return
        
        uri = sys.argv[2]
        user = sys.argv[3]
        password = sys.argv[4]
        
        if verify_neo4j(uri, user, password):
            print("\n✓ Neo4j is ready for data loading!")
            print("  Next: python setup_and_load.py load <uri> <user> <password>")
        else:
            print("\n✗ Neo4j setup required. Run: python setup_and_load.py setup")
    
    elif command == "load":
        if len(sys.argv) < 5:
            print("✗ Missing arguments")
            print("  Usage: python setup_and_load.py load <uri> <user> <password>")
            return
        
        uri = sys.argv[2]
        user = sys.argv[3]
        password = sys.argv[4]
        
        print_section("Pre-flight checks")
        if not verify_neo4j(uri, user, password):
            return
        
        if load_data(uri, user, password):
            print("\n✓ Data loading complete!")
            print("  Data is ready for analysis and reasoning")
        else:
            print("\n✗ Data loading failed")
    
    elif command == "full":
        if len(sys.argv) < 5:
            print("✗ Missing arguments")
            print("  Usage: python setup_and_load.py full <uri> <user> <password>")
            return
        
        uri = sys.argv[2]
        user = sys.argv[3]
        password = sys.argv[4]
        
        print_section("Pre-flight checks")
        if not verify_neo4j(uri, user, password):
            return
        
        if load_data(uri, user, password):
            if test_algorithms(uri, user, password):
                print_header("✓ COMPLETE SUCCESS")
                print("\nYour ChemAI Knowledge Graph is ready!")
                print("\nNext steps:")
                print("  1. Use graph algorithms for drug discovery")
                print("  2. Integrate with ML prediction pipeline")
                print("  3. Generate hypotheses and explanations")
                print("\nDocumentation:")
                print("  • See: IMPLEMENTATION_GUIDE.md")
                print("  • See: README.md")
            else:
                print("\n⚠ Data loaded but algorithm testing had issues")
        else:
            print("\n✗ Workflow failed at data loading")
    
    else:
        print(f"✗ Unknown command: {command}")
        print("  Available: setup, verify, load, full")


if __name__ == "__main__":
    main()
