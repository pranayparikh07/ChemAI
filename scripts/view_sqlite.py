#!/usr/bin/env python3
"""Simple SQLite database viewer for ChEMBL database"""

import sqlite3
import sys

db_path = "chembl_36/chembl_36_sqlite/chembl_36.db"

try:
    conn = sqlite3.connect(db_path)
    conn.row_factory = sqlite3.Row
    cursor = conn.cursor()
    
    print("=" * 100)
    print(f"ChEMBL v36 SQLite Database Viewer")
    print(f"Database: {db_path}")
    print("=" * 100)
    
    # Get all tables
    cursor.execute("SELECT name FROM sqlite_master WHERE type='table' ORDER BY name;")
    tables = cursor.fetchall()
    
    print(f"\n📊 Total Tables: {len(tables)}\n")
    
    for table in tables:
        table_name = table[0]
        
        # Get row count
        cursor.execute(f"SELECT COUNT(*) FROM {table_name}")
        row_count = cursor.fetchone()[0]
        
        # Get columns
        cursor.execute(f"PRAGMA table_info({table_name})")
        columns = cursor.fetchall()
        col_names = [col[1] for col in columns]
        
        print(f"\n{'─' * 100}")
        print(f"📋 TABLE: {table_name} ({row_count:,} rows)")
        print(f"{'─' * 100}")
        print(f"Columns ({len(columns)}): {', '.join(col_names)}")
        
        # Show first 3 rows
        try:
            cursor.execute(f"SELECT * FROM {table_name} LIMIT 3")
            rows = cursor.fetchall()
            
            if rows:
                print(f"\nFirst {min(3, len(rows))} rows:")
                print("─" * 100)
                
                for i, row in enumerate(rows, 1):
                    print(f"\nRow {i}:")
                    for col_name, value in zip(col_names, row):
                        if value is not None:
                            value_str = str(value)[:80]  # Limit string length
                            print(f"  {col_name:<30} = {value_str}")
        except Exception as e:
            print(f"Error reading data: {e}")
        
        print()
    
    conn.close()
    print("=" * 100)
    print("✅ Database view complete!")
    
except FileNotFoundError:
    print(f"❌ Database file not found: {db_path}")
    sys.exit(1)
except Exception as e:
    print(f"❌ Error: {e}")
    sys.exit(1)
