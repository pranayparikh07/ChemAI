#!/usr/bin/env python3
"""
Data Loaders for ChemAI
Create train/validation/test splits and data loading utilities
"""

import pandas as pd
import numpy as np
import pickle
from sklearn.model_selection import train_test_split
from pathlib import Path


def create_data_splits(feature_matrix_file, test_size=0.15, val_size=0.15, random_state=42):
    """Create train/val/test splits"""
    df = pd.read_csv(feature_matrix_file)
    
    # Get fingerprint columns
    fp_cols = [col for col in df.columns if col.startswith('FP_')]
    
    # Separate SMILES from features
    smiles = df['smiles'].values
    X = df[fp_cols].values
    
    # First split: train + val vs test
    X_temp, X_test, smiles_temp, smiles_test = train_test_split(
        X, smiles, test_size=test_size, random_state=random_state
    )
    
    # Second split: train vs val
    val_size_adj = val_size / (1 - test_size)
    X_train, X_val, smiles_train, smiles_val = train_test_split(
        X_temp, smiles_temp, test_size=val_size_adj, random_state=random_state
    )
    
    print(f"Dataset splits created:")
    print(f"  Training:   {len(X_train)} molecules ({100*len(X_train)/len(X):.1f}%)")
    print(f"  Validation: {len(X_val)} molecules ({100*len(X_val)/len(X):.1f}%)")
    print(f"  Testing:    {len(X_test)} molecules ({100*len(X_test)/len(X):.1f}%)")
    
    return {
        'X_train': X_train, 'smiles_train': smiles_train,
        'X_val': X_val, 'smiles_val': smiles_val,
        'X_test': X_test, 'smiles_test': smiles_test,
    }


def save_splits(splits, output_dir):
    """Save splits to pickle files"""
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    for key, value in splits.items():
        if isinstance(value, np.ndarray):
            with open(f'{output_dir}/{key}.pkl', 'wb') as f:
                pickle.dump(value, f)
    
    print(f"✓ Splits saved to {output_dir}")


class DataLoader:
    """Utility for loading data batches"""
    
    def __init__(self, X, smiles, batch_size=32):
        self.X = X
        self.smiles = smiles
        self.batch_size = batch_size
        self.n_samples = len(X)
    
    def __iter__(self):
        indices = np.random.permutation(self.n_samples)
        for i in range(0, self.n_samples, self.batch_size):
            batch_indices = indices[i:i+self.batch_size]
            yield self.X[batch_indices], self.smiles[batch_indices]


if __name__ == "__main__":
    splits = create_data_splits('d:\\ChemAI\\TEAM_VISHWA\\data\\feature_matrix.csv')
    save_splits(splits, 'd:\\ChemAI\\TEAM_VISHWA\\data\\splits')
