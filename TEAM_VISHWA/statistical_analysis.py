#!/usr/bin/env python3
"""
Statistical Analysis of Molecular Dataset
Generates comprehensive analysis and visualizations
"""

import pandas as pd
import numpy as np
import json
from pathlib import Path


def statistical_analysis(feature_matrix_file):
    """Perform comprehensive statistical analysis"""
    df = pd.read_csv(feature_matrix_file)
    
    # Select property columns (exclude SMILES and fingerprints)
    property_cols = ['MW', 'LogP', 'HBD', 'HBA', 'PSA', 'RotBonds', 'AromRings', 
                     'Rings', 'QED', 'LabuteASA', 'FractionCsp3']
    
    stats = {}
    for col in property_cols:
        if col in df.columns:
            stats[col] = {
                'mean': float(df[col].mean()),
                'median': float(df[col].median()),
                'std': float(df[col].std()),
                'min': float(df[col].min()),
                'max': float(df[col].max()),
                '25%': float(df[col].quantile(0.25)),
                '75%': float(df[col].quantile(0.75)),
                'count': int(df[col].count()),
            }
    
    # Lipinski compliance
    stats['lipinski_compliance'] = {
        'compliant': int(df['lipinski_compliant'].sum()),
        'non_compliant': int((~df['lipinski_compliant']).sum()),
        'compliance_rate': float(df['lipinski_compliant'].mean()),
    }
    
    # Drug-likeness
    stats['drug_like'] = {
        'drug_like': int(df['drug_like'].sum()),
        'non_drug_like': int((~df['drug_like']).sum()),
        'drug_like_rate': float(df['drug_like'].mean()),
    }
    
    return stats


def generate_report(stats, output_file):
    """Generate human-readable report"""
    report = []
    report.append("=" * 80)
    report.append("STATISTICAL ANALYSIS OF MOLECULAR DATASET")
    report.append("=" * 80)
    report.append("")
    
    report.append("PROPERTY STATISTICS")
    report.append("-" * 80)
    for prop, values in stats.items():
        if prop not in ['lipinski_compliance', 'drug_like']:
            report.append(f"\n{prop}:")
            report.append(f"  Mean:   {values['mean']:.2f}")
            report.append(f"  Median: {values['median']:.2f}")
            report.append(f"  Std:    {values['std']:.2f}")
            report.append(f"  Range:  {values['min']:.2f} - {values['max']:.2f}")
            report.append(f"  Q1-Q3:  {values['25%']:.2f} - {values['75%']:.2f}")
            report.append(f"  Count:  {values['count']}")
    
    report.append("\n" + "=" * 80)
    report.append("DRUG-LIKENESS METRICS")
    report.append("=" * 80)
    
    lip = stats['lipinski_compliance']
    report.append(f"\nLipinski's Rule of 5 Compliance:")
    report.append(f"  Compliant:     {lip['compliant']} ({100*lip['compliance_rate']:.1f}%)")
    report.append(f"  Non-compliant: {lip['non_compliant']} ({100*(1-lip['compliance_rate']):.1f}%)")
    
    drug = stats['drug_like']
    report.append(f"\nDrug-like Classification (QED > 0.6 AND Lipinski):")
    report.append(f"  Drug-like:     {drug['drug_like']} ({100*drug['drug_like_rate']:.1f}%)")
    report.append(f"  Non-drug-like: {drug['non_drug_like']} ({100*(1-drug['drug_like_rate']):.1f}%)")
    
    report_text = "\n".join(report)
    
    with open(output_file, 'w') as f:
        f.write(report_text)
    
    print(report_text)


if __name__ == "__main__":
    stats = statistical_analysis('d:\\ChemAI\\TEAM_VISHWA\\data\\feature_matrix.csv')
    generate_report(stats, 'd:\\ChemAI\\TEAM_VISHWA\\statistical_analysis.txt')
    
    # Save as JSON
    with open('d:\\ChemAI\\TEAM_VISHWA\\statistics.json', 'w') as f:
        json.dump(stats, f, indent=2)
