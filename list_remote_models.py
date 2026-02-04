#!/usr/bin/env python
"""Check remote server and download models"""

import urllib.request
import re

url = "http://192.168.137.4:5500/models"

try:
    response = urllib.request.urlopen(url)
    content = response.read().decode()
    
    # Extract all href links
    files = re.findall(r'href="([^"]+)"', content)
    
    print("Files available on remote server:")
    print("=" * 70)
    
    joblib_files = [f for f in files if '.joblib' in f and not f.startswith('?')]
    
    if joblib_files:
        for f in joblib_files:
            print(f"  - {f}")
    else:
        print("No .joblib files found")
        print("\nAll files found:")
        for f in files:
            if not f.startswith('?') and f.strip():
                print(f"  - {f}")
                
except Exception as e:
    print(f"Error: {e}")
