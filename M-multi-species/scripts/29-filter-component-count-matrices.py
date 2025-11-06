#!/usr/bin/env python3
"""
Filter count matrices based on ortholog groups from a component.

This script:
1. Reads OG_IDs from a component file (e.g., Component_24_top100.csv)
2. Maps OG_IDs to species-specific gene IDs from ortholog_groups_annotated.csv
3. Filters count matrices for each species to include only genes from the component
4. Saves filtered matrices to the output directory
"""

import sys
import os
import csv
import pandas as pd


def read_component_og_ids(component_file):
    """
    Read OG_IDs from a component CSV file.
    
    Args:
        component_file: Path to component CSV file (e.g., Component_24_top100.csv)
        
    Returns:
        List of OG_IDs
    """
    og_ids = []
    with open(component_file, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            og_ids.append(row['OG_ID'])
    
    print(f"Read {len(og_ids)} OG_IDs from {component_file}")
    return og_ids


def map_og_to_species_genes(ortholog_file, og_ids):
    """
    Map OG_IDs to species-specific gene IDs.
    
    Args:
        ortholog_file: Path to ortholog_groups_annotated.csv
        og_ids: List of OG_IDs to map
        
    Returns:
        Dictionary with keys 'apul', 'peve', 'ptua', each containing list of gene IDs
    """
    og_set = set(og_ids)
    gene_mapping = {
        'apul': [],
        'peve': [],
        'ptua': []
    }
    
    with open(ortholog_file, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row['group_id'] in og_set:
                # Add gene IDs if they exist (not empty)
                if row['apul']:
                    gene_mapping['apul'].append(row['apul'])
                if row['peve']:
                    gene_mapping['peve'].append(row['peve'])
                if row['ptua']:
                    gene_mapping['ptua'].append(row['ptua'])
    
    print(f"\nMapped to species-specific genes:")
    print(f"  apul: {len(gene_mapping['apul'])} genes")
    print(f"  peve: {len(gene_mapping['peve'])} genes")
    print(f"  ptua: {len(gene_mapping['ptua'])} genes")
    
    return gene_mapping


def filter_count_matrix(count_matrix_file, gene_ids, output_file, gene_prefix="", strip_suffix=""):
    """
    Filter a count matrix to include only specified genes.
    
    Args:
        count_matrix_file: Path to input count matrix CSV
        gene_ids: List of gene IDs to keep
        output_file: Path to output filtered CSV
        gene_prefix: Prefix to add to gene IDs when looking them up (e.g., "gene-")
        strip_suffix: Suffix to remove from gene IDs before matching (e.g., "-T1")
    """
    # Create set of gene IDs to keep (with prefix and/or without suffix if needed)
    processed_ids = []
    for gene_id in gene_ids:
        # Strip suffix if specified
        if strip_suffix and gene_id.endswith(strip_suffix):
            gene_id = gene_id[:-len(strip_suffix)]
        # Add prefix if specified
        if gene_prefix:
            gene_id = f"{gene_prefix}{gene_id}"
        processed_ids.append(gene_id)
    
    genes_to_keep = set(processed_ids)
    
    print(f"\nFiltering {count_matrix_file}...")
    print(f"Looking for {len(genes_to_keep)} genes...")
    
    # Read the count matrix
    df = pd.read_csv(count_matrix_file)
    
    # Filter rows where gene_id is in our set
    filtered_df = df[df['gene_id'].isin(genes_to_keep)]
    
    print(f"Found {len(filtered_df)} genes in count matrix")
    
    # Save filtered matrix
    filtered_df.to_csv(output_file, index=False)
    print(f"Saved to {output_file}")
    
    return len(filtered_df)


def main():
    # Base paths
    base_dir = "/home/runner/work/timeseries_molecular/timeseries_molecular"
    
    # Input files
    component_file = os.path.join(
        base_dir,
        "M-multi-species/output/26-rank35-optimization/lambda_gene_0.2/top_genes_per_component",
        "Component_24_top100.csv"
    )
    
    ortholog_file = os.path.join(
        base_dir,
        "M-multi-species/output/12-ortho-annot",
        "ortholog_groups_annotated.csv"
    )
    
    # Count matrix files
    apul_matrix = os.path.join(
        base_dir,
        "D-Apul/output/02.20-D-Apul-RNAseq-alignment-HiSat2",
        "apul-gene_count_matrix.csv"
    )
    
    peve_matrix = os.path.join(
        base_dir,
        "E-Peve/output/02.20-E-Peve-RNAseq-alignment-HiSat2",
        "peve-gene_count_matrix.csv"
    )
    
    ptua_matrix = os.path.join(
        base_dir,
        "F-Ptua/output/02.20-F-Ptua-RNAseq-alignment-HiSat2",
        "ptua-gene_count_matrix.csv"
    )
    
    # Output directory
    output_dir = os.path.join(
        base_dir,
        "M-multi-species/output/26-rank35-optimization/lambda_gene_0.2/top_genes_per_component"
    )
    
    # Verify input files exist
    for filepath in [component_file, ortholog_file, apul_matrix, peve_matrix, ptua_matrix]:
        if not os.path.exists(filepath):
            print(f"Error: Input file not found: {filepath}")
            sys.exit(1)
    
    print("=== Generating Component 24 Count Matrices ===\n")
    
    # Step 1: Read OG_IDs from component file
    og_ids = read_component_og_ids(component_file)
    
    # Step 2: Map OG_IDs to species-specific gene IDs
    gene_mapping = map_og_to_species_genes(ortholog_file, og_ids)
    
    # Step 3: Filter count matrices for each species
    # Note: 
    #   - apul IDs in ortholog file have "-T1" suffix but count matrix doesn't
    #   - peve and ptua have "gene-" prefix in count matrix
    results = {}
    
    results['apul'] = filter_count_matrix(
        apul_matrix,
        gene_mapping['apul'],
        os.path.join(output_dir, "Component_24_apul_count_matrix.csv"),
        gene_prefix="",
        strip_suffix="-T1"
    )
    
    results['peve'] = filter_count_matrix(
        peve_matrix,
        gene_mapping['peve'],
        os.path.join(output_dir, "Component_24_peve_count_matrix.csv"),
        gene_prefix="gene-"
    )
    
    results['ptua'] = filter_count_matrix(
        ptua_matrix,
        gene_mapping['ptua'],
        os.path.join(output_dir, "Component_24_ptua_count_matrix.csv"),
        gene_prefix="gene-"
    )
    
    print("\n=== Summary ===")
    print(f"Total OG_IDs in Component_24: {len(og_ids)}")
    print(f"Genes found in count matrices:")
    print(f"  apul: {results['apul']}")
    print(f"  peve: {results['peve']}")
    print(f"  ptua: {results['ptua']}")
    print("\nFiltered count matrices saved successfully!")


if __name__ == "__main__":
    main()
