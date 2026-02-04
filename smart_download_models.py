#!/usr/bin/env python
"""Advanced model downloader with multiple format support"""

import os
import urllib.request
import urllib.error
from pathlib import Path

# Try different URL formats
URLS_TO_TRY = [
    "http://192.168.137.4:5500/models",
    "http://192.168.137.4:5500/trained_models",
    "http://192.168.137.4:5500",
    "\\\\192.168.137.4\\models",
    "\\\\192.168.137.4\\shared\\models",
]

LOCAL_DIR = "trained_models"

MODELS = [
    "property_model.joblib",
    "property_scaler.joblib",
    "qed_model.joblib",
    "druglikeness_model.joblib",
    "bioactivity_model.joblib",
    "toxicity_model.joblib",
]

def test_url(url):
    """Test if URL is accessible"""
    try:
        response = urllib.request.urlopen(f"{url}/property_model.joblib", timeout=3)
        return True
    except:
        return False

def download_from_http(base_url):
    """Download models from HTTP server"""
    print(f"\nTrying HTTP: {base_url}")
    print("=" * 70)
    
    os.makedirs(LOCAL_DIR, exist_ok=True)
    downloaded = 0
    
    for model in MODELS:
        url = f"{base_url}/{model}"
        local_path = os.path.join(LOCAL_DIR, model)
        
        try:
            print(f"Downloading {model}...", end=" ")
            urllib.request.urlretrieve(url, local_path, timeout=10)
            size = os.path.getsize(local_path)
            print(f"✓ ({size:,} bytes)")
            downloaded += 1
        except Exception as e:
            print(f"✗")
    
    return downloaded

def copy_from_network_share(unc_path):
    """Copy models from network share"""
    print(f"\nTrying Network Share: {unc_path}")
    print("=" * 70)
    
    import shutil
    os.makedirs(LOCAL_DIR, exist_ok=True)
    downloaded = 0
    
    for model in MODELS:
        source = os.path.join(unc_path, model)
        dest = os.path.join(LOCAL_DIR, model)
        
        try:
            print(f"Copying {model}...", end=" ")
            if os.path.exists(source):
                shutil.copy2(source, dest)
                size = os.path.getsize(dest)
                print(f"✓ ({size:,} bytes)")
                downloaded += 1
            else:
                print(f"✗ File not found")
        except Exception as e:
            print(f"✗ {e}")
    
    return downloaded

def main():
    print("=" * 70)
    print("ChemAI Model Downloader")
    print("=" * 70)
    
    # Try each URL
    for url in URLS_TO_TRY:
        if url.startswith("http"):
            downloaded = download_from_http(url)
            if downloaded > 0:
                print(f"\n✓ Successfully downloaded {downloaded} models!")
                return
        else:
            downloaded = copy_from_network_share(url)
            if downloaded > 0:
                print(f"\n✓ Successfully copied {downloaded} models!")
                return
    
    print("\n✗ Failed to download models from any source")
    print("\nPlease:")
    print("1. Verify the server IP and port: 192.168.137.4:5500")
    print("2. Check if the models folder path is correct")
    print("3. Or provide the network share path (e.g., \\\\PC_NAME\\shared\\models)")

if __name__ == "__main__":
    main()
