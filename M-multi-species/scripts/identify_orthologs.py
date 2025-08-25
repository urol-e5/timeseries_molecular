#!/usr/bin/env python3
"""
Simple script to identify reciprocal best hits from BLAST results.
This demonstrates the orthology identification approach.
"""

import pandas as pd
import sys
import os

def load_blast_results(blast_file):
    """Load BLAST results with proper column names."""
    columns = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 
               'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'qlen', 'slen']
    
    if not os.path.exists(blast_file):
        print(f"Warning: {blast_file} not found")
        return pd.DataFrame(columns=columns)
    
    try:
        df = pd.read_csv(blast_file, sep='\t', header=None, names=columns)
        return df
    except Exception as e:
        print(f"Error reading {blast_file}: {e}")
        return pd.DataFrame(columns=columns)

def find_reciprocal_best_hits(blast1, blast2, min_identity=30, min_coverage=50):
    """Find reciprocal best hits between two BLAST results."""
    
    # Filter and get best hits for each query
    def get_best_hits(df):
        if df.empty:
            return df
        
        # Calculate coverage
        df['qcoverage'] = (df['qend'] - df['qstart'] + 1) / df['qlen'] * 100
        df['scoverage'] = (df['send'] - df['sstart'] + 1) / df['slen'] * 100
        
        # Filter by identity and coverage
        filtered = df[
            (df['pident'] >= min_identity) & 
            ((df['qcoverage'] >= min_coverage) | (df['scoverage'] >= min_coverage))
        ]
        
        if filtered.empty:
            return filtered
        
        # Get best hit per query (highest bitscore)
        best_hits = filtered.loc[filtered.groupby('qseqid')['bitscore'].idxmax()]
        return best_hits
    
    best1 = get_best_hits(blast1)
    best2 = get_best_hits(blast2)
    
    if best1.empty or best2.empty:
        print("No hits meet criteria")
        return pd.DataFrame()
    
    # Find reciprocal best hits
    # Merge where query1 -> subject1 and subject1 -> query1
    rbh = pd.merge(
        best1[['qseqid', 'sseqid', 'pident', 'evalue', 'bitscore']],
        best2[['qseqid', 'sseqid', 'pident', 'evalue', 'bitscore']],
        left_on=['qseqid', 'sseqid'],
        right_on=['sseqid', 'qseqid'],
        suffixes=('_1', '_2')
    )
    
    if rbh.empty:
        print("No reciprocal best hits found")
        return rbh
    
    # Clean up column names
    rbh = rbh.rename(columns={
        'qseqid_1': 'gene1',
        'sseqid_1': 'gene2',
        'pident_1': 'identity1',
        'pident_2': 'identity2',
        'evalue_1': 'evalue1',
        'evalue_2': 'evalue2',
        'bitscore_1': 'bitscore1',
        'bitscore_2': 'bitscore2'
    })
    
    rbh['avg_identity'] = (rbh['identity1'] + rbh['identity2']) / 2
    rbh['min_evalue'] = rbh[['evalue1', 'evalue2']].min(axis=1)
    rbh['max_bitscore'] = rbh[['bitscore1', 'bitscore2']].max(axis=1)
    
    return rbh[['gene1', 'gene2', 'avg_identity', 'min_evalue', 'max_bitscore']]

def main():
    # Paths to BLAST results
    output_dir = "/home/runner/work/timeseries_molecular/timeseries_molecular/M-multi-species/output/11-orthology-analysis"
    
    # Load test BLAST results
    print("Loading BLAST results...")
    apul_vs_peve = load_blast_results(f"{output_dir}/test_Apul_vs_Peve.blastp")
    peve_vs_apul = load_blast_results(f"{output_dir}/test_Peve_vs_Apul.blastp")
    
    print(f"Apul vs Peve hits: {len(apul_vs_peve)}")
    print(f"Peve vs Apul hits: {len(peve_vs_apul)}")
    
    # Find reciprocal best hits
    print("\nFinding reciprocal best hits...")
    rbh = find_reciprocal_best_hits(apul_vs_peve, peve_vs_apul)
    
    print(f"Reciprocal best hits found: {len(rbh)}")
    
    if not rbh.empty:
        print("\nReciprocal Best Hits:")
        print(rbh.to_string(index=False))
        
        # Save results
        rbh.to_csv(f"{output_dir}/test_apul_peve_orthologs.csv", index=False)
        print(f"\nResults saved to {output_dir}/test_apul_peve_orthologs.csv")
        
        # Summary statistics
        print(f"\nSummary Statistics:")
        print(f"Average identity: {rbh['avg_identity'].mean():.2f}%")
        print(f"Range: {rbh['avg_identity'].min():.2f}% - {rbh['avg_identity'].max():.2f}%")
        print(f"Min E-value: {rbh['min_evalue'].min():.2e}")
    else:
        print("No orthologous pairs found in test dataset")

if __name__ == "__main__":
    main()