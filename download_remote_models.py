#!/usr/bin/env python
"""Download model files from remote server"""

import os
import urllib.request
import urllib.error
from pathlib import Path

# Remote server URL
REMOTE_URL = "http://192.168.137.4:5500/models"
LOCAL_DIR = "trained_models"

# Models to download
MODELS = [
    "property_model.joblib",
    "property_scaler.joblib",
    "property_config.joblib",
    "qed_model.joblib",
    "druglikeness_model.joblib",
    "druglikeness_config.joblib",
    "bioactivity_model.joblib",
    "bioactivity_config.joblib",
    "toxicity_model.joblib",
    "toxicity_config.joblib",
    "structural_alerts.joblib",
    ".checkpoint_property.joblib",
    ".checkpoint_druglikeness.joblib",
    ".checkpoint_bioactivity.joblib"
]

def download_model(filename):
    """Download a single model file"""
    # Handle filenames that start with dot by properly encoding them
    encoded_filename = filename.replace('.', '%2E') if filename.startswith('.') else filename
    remote_path = f"{REMOTE_URL}/{encoded_filename}"
    local_path = os.path.join(LOCAL_DIR, filename)
    
    try:
        print(f"Downloading {filename}...", end=" ")
        urllib.request.urlretrieve(remote_path, local_path)
        file_size = os.path.getsize(local_path)
        print(f"✓ ({file_size:,} bytes)")
        return True
    except urllib.error.HTTPError as e:
        if e.code == 404:
            print(f"✗ Not found (404)")
        else:
            print(f"✗ HTTP {e.code}")
        return False
    except Exception as e:
        print(f"✗ Error: {e}")
        return False

def main():
    """Download all models"""
    # Create directory if needed
    os.makedirs(LOCAL_DIR, exist_ok=True)
    
    print("=" * 70)
    print("Downloading Models from Remote Server")
    print(f"Source: {REMOTE_URL}")
    print(f"Destination: {LOCAL_DIR}")
    print("=" * 70)
    
    downloaded = 0
    failed = 0
    
    for model in MODELS:
        if download_model(model):
            downloaded += 1
        else:
            failed += 1
    
    print("=" * 70)
    print(f"Downloaded: {downloaded}/{len(MODELS)}")
    print(f"Failed: {failed}/{len(MODELS)}")
    print("=" * 70)

if __name__ == "__main__":
    main()
