#!/usr/bin/env python3
"""
Comprehensive orthology analysis for three coral species.
This script identifies orthologous proteins using reciprocal best hits approach.
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
        return pd.DataFrame()
    
    # Find reciprocal best hits
    rbh = pd.merge(
        best1[['qseqid', 'sseqid', 'pident', 'evalue', 'bitscore']],
        best2[['qseqid', 'sseqid', 'pident', 'evalue', 'bitscore']],
        left_on=['qseqid', 'sseqid'],
        right_on=['sseqid', 'qseqid'],
        suffixes=('_1', '_2')
    )
    
    if rbh.empty:
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

def create_ortholog_groups(apul_peve_rbh, apul_ptua_rbh, peve_ptua_rbh):
    """Create ortholog groups from pairwise reciprocal best hits."""
    
    ortholog_groups = []
    group_id = 1
    
    # Track processed proteins
    processed_apul = set()
    processed_peve = set()
    processed_ptua = set()
    
    # Three-way orthologs (present in all species)
    for _, row_ap in apul_peve_rbh.iterrows():
        apul_gene = row_ap['gene1']
        peve_gene = row_ap['gene2']
        
        if apul_gene in processed_apul or peve_gene in processed_peve:
            continue
            
        # Check if this Apul gene also has ortholog in Ptua
        apul_ptua_match = apul_ptua_rbh[apul_ptua_rbh['gene1'] == apul_gene]
        if not apul_ptua_match.empty:
            ptua_gene = apul_ptua_match.iloc[0]['gene2']
            
            if ptua_gene in processed_ptua:
                continue
            
            # Check if Peve and Ptua are also orthologs
            peve_ptua_match = peve_ptua_rbh[
                (peve_ptua_rbh['gene1'] == peve_gene) & 
                (peve_ptua_rbh['gene2'] == ptua_gene)
            ]
            
            if not peve_ptua_match.empty:
                # Three-way ortholog found
                avg_identity = (row_ap['avg_identity'] + 
                              apul_ptua_match.iloc[0]['avg_identity'] + 
                              peve_ptua_match.iloc[0]['avg_identity']) / 3
                
                ortholog_groups.append({
                    'group_id': f"OG_{group_id:05d}",
                    'apul': apul_gene,
                    'peve': peve_gene,
                    'ptua': ptua_gene,
                    'type': 'three_way',
                    'avg_identity': avg_identity
                })
                
                processed_apul.add(apul_gene)
                processed_peve.add(peve_gene)
                processed_ptua.add(ptua_gene)
                group_id += 1
    
    # Two-way orthologs (remaining pairs)
    for _, row in apul_peve_rbh.iterrows():
        if row['gene1'] not in processed_apul and row['gene2'] not in processed_peve:
            ortholog_groups.append({
                'group_id': f"OG_{group_id:05d}",
                'apul': row['gene1'],
                'peve': row['gene2'],
                'ptua': None,
                'type': 'apul_peve',
                'avg_identity': row['avg_identity']
            })
            processed_apul.add(row['gene1'])
            processed_peve.add(row['gene2'])
            group_id += 1
    
    for _, row in apul_ptua_rbh.iterrows():
        if row['gene1'] not in processed_apul and row['gene2'] not in processed_ptua:
            ortholog_groups.append({
                'group_id': f"OG_{group_id:05d}",
                'apul': row['gene1'],
                'peve': None,
                'ptua': row['gene2'],
                'type': 'apul_ptua',
                'avg_identity': row['avg_identity']
            })
            processed_apul.add(row['gene1'])
            processed_ptua.add(row['gene2'])
            group_id += 1
    
    for _, row in peve_ptua_rbh.iterrows():
        if row['gene1'] not in processed_peve and row['gene2'] not in processed_ptua:
            ortholog_groups.append({
                'group_id': f"OG_{group_id:05d}",
                'apul': None,
                'peve': row['gene1'],
                'ptua': row['gene2'],
                'type': 'peve_ptua',
                'avg_identity': row['avg_identity']
            })
            processed_peve.add(row['gene1'])
            processed_ptua.add(row['gene2'])
            group_id += 1
    
    return pd.DataFrame(ortholog_groups)

def main():
    """Main analysis function."""
    
    output_dir = "/home/runner/work/timeseries_molecular/timeseries_molecular/M-multi-species/output/11-orthology-analysis"
    
    print("=== Coral Orthology Analysis ===")
    print("Species:")
    print("  - Acropora pulchra (Apul)")
    print("  - Porites evermanni (Peve)")
    print("  - Pocillopora tuahiniensis (Ptua)")
    print()
    
    # Load all BLAST results
    print("Loading BLAST results...")
    apul_vs_peve = load_blast_results(f"{output_dir}/test_Apul_vs_Peve.blastp")
    peve_vs_apul = load_blast_results(f"{output_dir}/test_Peve_vs_Apul.blastp")
    apul_vs_ptua = load_blast_results(f"{output_dir}/test_Apul_vs_Ptua.blastp")
    ptua_vs_apul = load_blast_results(f"{output_dir}/test_Ptua_vs_Apul.blastp")
    peve_vs_ptua = load_blast_results(f"{output_dir}/test_Peve_vs_Ptua.blastp")
    ptua_vs_peve = load_blast_results(f"{output_dir}/test_Ptua_vs_Peve.blastp")
    
    print(f"BLAST hits loaded:")
    print(f"  Apul vs Peve: {len(apul_vs_peve)}")
    print(f"  Peve vs Apul: {len(peve_vs_apul)}")
    print(f"  Apul vs Ptua: {len(apul_vs_ptua)}")
    print(f"  Ptua vs Apul: {len(ptua_vs_apul)}")
    print(f"  Peve vs Ptua: {len(peve_vs_ptua)}")
    print(f"  Ptua vs Peve: {len(ptua_vs_peve)}")
    print()
    
    # Find reciprocal best hits for each pair
    print("Identifying reciprocal best hits...")
    apul_peve_rbh = find_reciprocal_best_hits(apul_vs_peve, peve_vs_apul)
    apul_ptua_rbh = find_reciprocal_best_hits(apul_vs_ptua, ptua_vs_apul)
    peve_ptua_rbh = find_reciprocal_best_hits(peve_vs_ptua, ptua_vs_peve)
    
    print(f"Reciprocal best hits:")
    print(f"  Apul-Peve: {len(apul_peve_rbh)}")
    print(f"  Apul-Ptua: {len(apul_ptua_rbh)}")
    print(f"  Peve-Ptua: {len(peve_ptua_rbh)}")
    print()
    
    # Create ortholog groups
    print("Creating ortholog groups...")
    ortholog_groups = create_ortholog_groups(apul_peve_rbh, apul_ptua_rbh, peve_ptua_rbh)
    
    # Summary statistics
    three_way = ortholog_groups[ortholog_groups['type'] == 'three_way']
    apul_peve_only = ortholog_groups[ortholog_groups['type'] == 'apul_peve']
    apul_ptua_only = ortholog_groups[ortholog_groups['type'] == 'apul_ptua']
    peve_ptua_only = ortholog_groups[ortholog_groups['type'] == 'peve_ptua']
    
    print("=== ORTHOLOGY RESULTS ===")
    print(f"Total ortholog groups: {len(ortholog_groups)}")
    print(f"Three-way orthologs: {len(three_way)}")
    print(f"Apul-Peve pairs only: {len(apul_peve_only)}")
    print(f"Apul-Ptua pairs only: {len(apul_ptua_only)}")
    print(f"Peve-Ptua pairs only: {len(peve_ptua_only)}")
    print()
    
    if len(three_way) > 0:
        print(f"Three-way orthologs (avg identity: {three_way['avg_identity'].mean():.2f}%):")
        for _, row in three_way.iterrows():
            print(f"  {row['group_id']}: {row['apul']} | {row['peve']} | {row['ptua']} ({row['avg_identity']:.1f}%)")
        print()
    
    # Save results
    ortholog_groups.to_csv(f"{output_dir}/coral_ortholog_groups_test.csv", index=False)
    apul_peve_rbh.to_csv(f"{output_dir}/apul_peve_rbh_test.csv", index=False)
    apul_ptua_rbh.to_csv(f"{output_dir}/apul_ptua_rbh_test.csv", index=False)
    peve_ptua_rbh.to_csv(f"{output_dir}/peve_ptua_rbh_test.csv", index=False)
    
    # Create summary statistics
    summary_stats = pd.DataFrame({
        'Metric': [
            'Total test proteins (Apul)',
            'Total test proteins (Peve)', 
            'Total test proteins (Ptua)',
            'Apul-Peve orthologs',
            'Apul-Ptua orthologs',
            'Peve-Ptua orthologs',
            'Three-way orthologs',
            'Two-way orthologs only',
            'Total ortholog groups'
        ],
        'Count': [
            350,  # From our test data
            203,
            1000,
            len(apul_peve_rbh),
            len(apul_ptua_rbh),
            len(peve_ptua_rbh),
            len(three_way),
            len(ortholog_groups) - len(three_way),
            len(ortholog_groups)
        ]
    })
    
    summary_stats.to_csv(f"{output_dir}/orthology_summary_test.csv", index=False)
    
    print("Results saved:")
    print(f"  - Ortholog groups: {output_dir}/coral_ortholog_groups_test.csv")
    print(f"  - Summary statistics: {output_dir}/orthology_summary_test.csv")
    print(f"  - Individual RBH files: apul_peve_rbh_test.csv, etc.")
    print()
    
    print("=== ANALYSIS NOTES ===")
    print("This is a test run on subset data (350 Apul, 203 Peve, 1000 Ptua proteins).")
    print("For complete analysis, run BLAST on full protein datasets.")
    print("Parameters used:")
    print("  - E-value threshold: 1e-5")
    print("  - Minimum identity: 30%")
    print("  - Minimum coverage: 50%")
    print("  - Method: Reciprocal best hits (RBH)")

if __name__ == "__main__":
    main()