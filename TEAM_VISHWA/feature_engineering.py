#!/usr/bin/env python3
"""
Feature Engineering for ChemAI
Extract molecular properties, fingerprints, and create feature matrix
"""

import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, QED
import pickle


def extract_properties(smiles):
    """Extract 15+ molecular properties from SMILES"""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        
        props = {
            'MW': Descriptors.MolWt(mol),
            'LogP': Descriptors.MolLogP(mol),
            'HBD': Descriptors.NumHDonors(mol),
            'HBA': Descriptors.NumHAcceptors(mol),
            'PSA': Descriptors.TPSA(mol),
            'RotBonds': Descriptors.NumRotatableBonds(mol),
            'AromRings': Descriptors.NumAromaticRings(mol),
            'Rings': Descriptors.RingCount(mol),
            'MaxRingSize': Chem.GetSSSR(mol).__len__(),
            'QED': QED.qed(mol),
            'LabuteASA': Descriptors.LabuteASA(mol),
            'TPSA': Descriptors.TPSA(mol),
            'NumSatCycles': Descriptors.NumSaturatedCycles(mol),
            'NumAliphaticCycles': Descriptors.NumAliphaticCycles(mol),
            'FractionCsp3': Descriptors.FractionCsp3(mol),
        }
        return props
    except:
        return None


def generate_fingerprints(smiles):
    """Generate Morgan fingerprints"""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        
        # Morgan fingerprint (2048 bits, radius 2)
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
        return list(fp)
    except:
        return None


def check_lipinski_compliance(props):
    """Check Lipinski's Rule of 5 compliance"""
    if props is None:
        return False
    
    violations = 0
    if props['MW'] > 500:
        violations += 1
    if props['LogP'] > 5:
        violations += 1
    if props['HBD'] > 5:
        violations += 1
    if props['HBA'] > 10:
        violations += 1
    
    return violations <= 1


def build_feature_matrix(input_file, output_file):
    """Build complete feature matrix"""
    print(f"Reading {input_file}...")
    df = pd.read_csv(input_file)
    
    print(f"Extracting properties for {len(df)} molecules...")
    df['properties'] = df['smiles'].apply(extract_properties)
    
    # Filter out failed extractions
    df = df[df['properties'].notna()].copy()
    
    print(f"Extracting fingerprints...")
    df['fingerprints'] = df['smiles'].apply(generate_fingerprints)
    df = df[df['fingerprints'].notna()].copy()
    
    # Expand properties into columns
    print(f"Expanding properties...")
    props_df = pd.DataFrame(df['properties'].tolist())
    df = pd.concat([df[['smiles']], props_df, 
                    pd.DataFrame(df['fingerprints'].tolist(), 
                                columns=[f'FP_{i}' for i in range(2048)])], axis=1)
    
    # Add Lipinski compliance
    df['lipinski_compliant'] = props_df.apply(check_lipinski_compliance, axis=1)
    df['drug_like'] = (props_df['QED'] > 0.6) & (df['lipinski_compliant'])
    
    print(f"Saving feature matrix...")
    df.to_csv(output_file, index=False)
    
    print(f"\n✓ Feature matrix created!")
    print(f"Molecules: {len(df)}")
    print(f"Features: {len(df.columns)}")
    print(f"Drug-like molecules: {df['drug_like'].sum()} ({100*df['drug_like'].sum()/len(df):.1f}%)")
    
    return df


if __name__ == "__main__":
    df = build_feature_matrix(
        'd:\\ChemAI\\TEAM_VISHWA\\data\\clean_molecules.csv',
        'd:\\ChemAI\\TEAM_VISHWA\\data\\feature_matrix.csv'
    )
