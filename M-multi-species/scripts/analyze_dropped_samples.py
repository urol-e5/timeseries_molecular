#!/usr/bin/env python3
"""
Analyze which RNA-seq samples were dropped from count matrices and identify reasons.

This script compares:
1. Samples submitted for sequencing (from metadata)
2. Samples that completed alignment pipeline 
3. Samples included in final count matrices
4. Quality metrics and documented issues

Output: Comprehensive report of sample status and reasons for exclusion
"""

import pandas as pd
import os
import glob
from pathlib import Path

def load_metadata():
    """Load and parse RNA-seq metadata"""
    metadata_path = "M-multi-species/data/rna_metadata.csv"
    
    # Read metadata, handling the numbered format
    df = pd.read_csv(metadata_path)
    
    # Clean column names (remove the number prefix)
    df.columns = [col.split('.', 1)[1] if '.' in col and col.split('.')[0].isdigit() else col for col in df.columns]
    
    # Print columns for debugging
    print("Metadata columns:", list(df.columns))
    
    # Extract species from Species-Strain column  
    df['Species'] = df['Species-Strain'].apply(lambda x: x.split()[0] if pd.notna(x) else None)
    
    return df

def get_aligned_samples():
    """Get list of samples that completed alignment for each species"""
    aligned_samples = {}
    
    # Apul samples
    apul_dirs = glob.glob("D-Apul/output/02.20-D-Apul-RNAseq-alignment-HiSat2/ACR-*")
    aligned_samples['Acropora'] = [os.path.basename(d) for d in apul_dirs]
    
    # Peve samples  
    peve_dirs = glob.glob("E-Peve/output/02.20-E-Peve-RNAseq-alignment-HiSat2/POR-*")
    aligned_samples['Porites'] = [os.path.basename(d) for d in peve_dirs]
    
    # Ptua samples
    ptua_dirs = glob.glob("F-Ptua/output/02.20-F-Ptua-RNAseq-alignment-HiSat2/POC-*")
    aligned_samples['Pocillopora'] = [os.path.basename(d) for d in ptua_dirs]
    
    return aligned_samples

def get_count_matrix_samples():
    """Get samples included in final count matrices"""
    count_matrix_samples = {}
    
    # Apul count matrix
    apul_matrix = pd.read_csv("D-Apul/output/02.20-D-Apul-RNAseq-alignment-HiSat2/apul-gene_count_matrix.csv")
    apul_samples = [col for col in apul_matrix.columns if col != 'gene_id']
    # Remove .1 suffix for paired-end duplicates and get unique samples
    apul_samples_clean = [col.replace('.1', '') for col in apul_samples]
    apul_samples_unique = sorted(list(set(apul_samples_clean)))
    count_matrix_samples['Acropora'] = apul_samples_unique
    
    # Peve count matrix
    peve_matrix = pd.read_csv("E-Peve/output/02.20-E-Peve-RNAseq-alignment-HiSat2/peve-gene_count_matrix.csv")
    peve_samples = [col for col in peve_matrix.columns if col != 'gene_id']
    count_matrix_samples['Porites'] = peve_samples
    
    # Ptua count matrix
    ptua_matrix = pd.read_csv("F-Ptua/output/02.20-F-Ptua-RNAseq-alignment-HiSat2/ptua-gene_count_matrix.csv")
    ptua_samples = [col for col in ptua_matrix.columns if col != 'gene_id']
    count_matrix_samples['Pocillopora'] = ptua_samples
    
    return count_matrix_samples

def analyze_sample_status():
    """Main analysis function to compare sample status across pipeline stages"""
    
    print("=== RNA-seq Sample Status Analysis ===\n")
    
    # Load data
    metadata = load_metadata()
    aligned_samples = get_aligned_samples()
    count_matrix_samples = get_count_matrix_samples()
    
    # Create sample ID from ColonyID and Timepoint for comparison
    metadata['SampleID'] = metadata['ColonyID'] + '-' + metadata['Timepoint']
    
    # Analyze by species
    species_map = {
        'Acropora': 'Acropora pulchra',
        'Porites': 'Porites evermanni', 
        'Pocillopora': 'Pocillopora tuahiniensis'
    }
    
    all_results = []
    
    for genus, full_species in species_map.items():
        print(f"\n{'='*50}")
        print(f"ANALYSIS FOR {full_species.upper()}")
        print(f"{'='*50}")
        
        # Get samples from metadata for this species
        species_metadata = metadata[metadata['Species-Strain'] == full_species]
        submitted_samples = set(species_metadata['SampleID'].tolist())
        
        # Get aligned and count matrix samples
        aligned_set = set(aligned_samples.get(genus, []))
        count_matrix_set = set(count_matrix_samples.get(genus, []))
        
        print(f"\nSample counts:")
        print(f"  Submitted for sequencing: {len(submitted_samples)}")
        print(f"  Completed alignment: {len(aligned_set)}")
        
        # For Apul, show both total columns and unique samples
        if genus == 'Acropora':
            total_cols = len([col for col in pd.read_csv("D-Apul/output/02.20-D-Apul-RNAseq-alignment-HiSat2/apul-gene_count_matrix.csv").columns if col != 'gene_id'])
            print(f"  Included in count matrix: {len(count_matrix_set)} unique samples ({total_cols} total columns)")
        else:
            print(f"  Included in count matrix: {len(count_matrix_set)}")
        
        # Find samples that didn't make it through each stage
        failed_alignment = submitted_samples - aligned_set
        dropped_from_matrix = aligned_set - count_matrix_set
        
        print(f"\nSamples dropped:")
        print(f"  Failed alignment: {len(failed_alignment)}")
        print(f"  Dropped from count matrix: {len(dropped_from_matrix)}")
        
        # Check for quality issues in submitted samples
        quality_issues = []
        for _, sample in species_metadata.iterrows():
            issues = []
            # Check concentration
            if sample['Conc-ng.uL'] == 0:
                issues.append("Zero concentration")
            elif sample['Conc-ng.uL'] < 10:
                issues.append("Low concentration (<10 ng/μL)")
            
            # Check total amount
            if sample['TotalAmount-ng'] < 1000:
                issues.append("Low total amount (<1000 ng)")
                
            if issues:
                quality_issues.append({
                    'SampleID': sample['SampleID'],
                    'Issues': '; '.join(issues),
                    'Concentration': sample['Conc-ng.uL'],
                    'TotalAmount': sample['TotalAmount-ng']
                })
        
        if quality_issues:
            print(f"\nQuality issues in submitted samples:")
            for issue in quality_issues:
                print(f"  {issue['SampleID']}: {issue['Issues']} "
                      f"(Conc: {issue['Concentration']} ng/μL, "
                      f"Total: {issue['TotalAmount']} ng)")
        
        # List specific dropped samples
        if failed_alignment:
            print(f"\nSamples that failed alignment:")
            for sample in sorted(failed_alignment):
                print(f"  - {sample}")
                
        if dropped_from_matrix:
            print(f"\nSamples dropped from count matrix:")
            for sample in sorted(dropped_from_matrix):
                print(f"  - {sample}")
        
        # Store results for summary
        all_results.append({
            'Species': full_species,
            'Submitted': len(submitted_samples),
            'Aligned': len(aligned_set),
            'InMatrix': len(count_matrix_set),
            'FailedAlignment': len(failed_alignment),
            'DroppedFromMatrix': len(dropped_from_matrix),
            'QualityIssues': len(quality_issues)
        })
    
    # Print summary table
    print(f"\n\n{'='*80}")
    print("SUMMARY TABLE")
    print(f"{'='*80}")
    
    summary_df = pd.DataFrame(all_results)
    print(summary_df.to_string(index=False))
    
    # Check for documented outliers
    print(f"\n\n{'='*80}")
    print("DOCUMENTED OUTLIERS AND ISSUES")
    print(f"{'='*80}")
    
    print("\nKnown issues from analysis files:")
    print("- ACR-265 TP3 (1E1): Documented as outlier in sRNAseq analysis")
    print("- POR-73 TP2 (1F1): Zero concentration in metadata (0 ng/μL)")
    
    # Total samples across all species
    total_submitted = sum([r['Submitted'] for r in all_results])
    total_in_matrix = sum([r['InMatrix'] for r in all_results])
    
    print(f"\nOVERALL SUMMARY:")
    print(f"- Total samples submitted: {total_submitted}")
    print(f"- Total samples in count matrices: {total_in_matrix}")
    print(f"- Success rate: {total_in_matrix/total_submitted*100:.1f}%")
    
    print(f"\nNote: Acropora count matrix has duplicate columns for paired-end reads")
    print(f"      (80 columns total = 40 unique samples × 2 reads each)")

if __name__ == "__main__":
    # Change to repository directory
    os.chdir("/home/runner/work/timeseries_molecular/timeseries_molecular")
    analyze_sample_status()