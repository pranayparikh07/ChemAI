#!/usr/bin/env python3
"""
Decompress data files after cloning the repository

Run this script to decompress all .csv.gz files:
    python decompress_data.py
"""

import gzip
from pathlib import Path

def decompress_csv():
    """Decompress all CSV files"""
    
    gz_files = sorted(Path('.').rglob('*.csv.gz'))
    
    if not gz_files:
        print("❌ No .csv.gz files found")
        return
    
    print(f"📦 Found {len(gz_files)} compressed files. Decompressing...\n")
    
    for gz_file in gz_files:
        output_file = gz_file.with_suffix('')
        file_size_mb = gz_file.stat().st_size / (1024 * 1024)
        
        print(f"Decompressing {gz_file.name} ({file_size_mb:.1f} MB)...", end=' ')
        
        try:
            with gzip.open(gz_file, 'rb') as f_in:
                with open(output_file, 'wb') as f_out:
                    f_out.write(f_in.read())
            
            expanded_size = output_file.stat().st_size / (1024 * 1024)
            print(f"✓ ({expanded_size:.1f} MB)")
        except Exception as e:
            print(f"❌ Error: {e}")

if __name__ == '__main__':
    decompress_csv()
    print("\n✅ Decompression complete!")
