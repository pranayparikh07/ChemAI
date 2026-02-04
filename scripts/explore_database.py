import sqlite3
import pandas as pd

db_path = "chembl_36/chembl_36_sqlite/chembl_36.db"

conn = sqlite3.connect(db_path)
cursor = conn.cursor()

cursor.execute("SELECT name FROM sqlite_master WHERE type='table' ORDER BY name;")
tables = cursor.fetchall()

print("=" * 80)
print(f"ChEMBL v36 Database Schema Explorer")
print(f"Database: {db_path}")
print("=" * 80)
print(f"\nTotal tables: {len(tables)}\n")

for table in tables:
    table_name = table[0]
    cursor.execute(f"PRAGMA table_info({table_name})")
    columns = cursor.fetchall()
    
    cursor.execute(f"SELECT COUNT(*) FROM {table_name}")
    row_count = cursor.fetchone()[0]
    
    print("-" * 80)
    print(f"TABLE: {table_name} ({row_count:,} rows)")
    print("-" * 80)
    print(f"{'Column Name':<40} {'Type':<15} {'Nullable':<10} {'PK'}")
    print("-" * 80)
    for col in columns:
        cid, name, dtype, notnull, default, pk = col
        nullable = "NO" if notnull else "YES"
        is_pk = "YES" if pk else ""
        print(f"{name:<40} {dtype:<15} {nullable:<10} {is_pk}")
    print()

key_tables = ['molecule_dictionary', 'compound_structures', 'activities', 'assays', 'target_dictionary']

print("\n" + "=" * 80)
print("SAMPLE DATA FROM KEY TABLES")
print("=" * 80)

for table_name in key_tables:
    try:
        df = pd.read_sql_query(f"SELECT * FROM {table_name} LIMIT 3", conn)
        print(f"\n--- {table_name} (first 3 rows) ---")
        print(df.to_string())
        print()
    except Exception as e:
        print(f"Could not read {table_name}: {e}")

print("\n" + "=" * 80)
print("SEARCHING FOR SMILES COLUMNS")
print("=" * 80)

for table in tables:
    table_name = table[0]
    cursor.execute(f"PRAGMA table_info({table_name})")
    columns = cursor.fetchall()
    for col in columns:
        col_name = col[1].lower()
        if 'smiles' in col_name or 'structure' in col_name or 'mol' in col_name:
            print(f"  {table_name}.{col[1]} ({col[2]})")

conn.close()
print("\nDone!")
