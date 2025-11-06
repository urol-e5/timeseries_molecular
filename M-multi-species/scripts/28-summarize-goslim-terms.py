#!/usr/bin/env python3
"""
Summarize GO Slim term occurrences across component gene annotation files.

This script:
1. Reads all CSV files with "annotation" in the filename from 
   M-multi-species/output/26-rank35-optimization/lambda_gene_0.2/top_genes_per_component
2. Counts occurrences of terms in the goslim_names column (semicolon-separated)
3. Creates a single CSV file with term counts
"""

import os
import sys
import pandas as pd
from pathlib import Path
from collections import Counter


def count_goslim_terms(annotation_files):
    """
    Count occurrences of GO Slim terms across all annotation files.
    
    Args:
        annotation_files: List of paths to annotation CSV files
    
    Returns:
        Counter object with term counts
    """
    all_terms = []
    
    for file_path in annotation_files:
        try:
            # Read the annotation file
            df = pd.read_csv(file_path)
            
            # Check if goslim_names column exists
            if 'goslim_names' not in df.columns:
                print(f"  Warning: 'goslim_names' column not found in {file_path.name}")
                continue
            
            # Extract all terms from goslim_names column
            for terms_str in df['goslim_names'].dropna():
                if terms_str and str(terms_str).strip():
                    # Split by semicolon and clean up each term
                    terms = [t.strip() for t in str(terms_str).split(';') if t.strip()]
                    all_terms.extend(terms)
        
        except Exception as e:
            print(f"  Error processing {file_path.name}: {e}")
            continue
    
    # Count occurrences
    term_counts = Counter(all_terms)
    return term_counts


def main():
    # Set up paths
    repo_root = Path(__file__).parent.parent
    component_dir = repo_root / 'output' / '26-rank35-optimization' / 'lambda_gene_0.2' / 'top_genes_per_component'
    output_file = component_dir / 'goslim_term_counts.csv'
    
    # Check if directory exists
    if not component_dir.exists():
        print(f"Error: Component directory not found: {component_dir}")
        sys.exit(1)
    
    # Find all annotation files
    annotation_files = sorted([f for f in component_dir.glob('*.csv') 
                              if 'annotation' in f.name.lower()])
    
    print(f"Found {len(annotation_files)} annotation files")
    
    if len(annotation_files) == 0:
        print("No annotation files found!")
        sys.exit(1)
    
    # Count GO Slim terms
    print("\nCounting GO Slim term occurrences...")
    term_counts = count_goslim_terms(annotation_files)
    
    print(f"Found {len(term_counts)} unique GO Slim terms")
    print(f"Total term occurrences: {sum(term_counts.values())}")
    
    # Create DataFrame and sort by count (descending)
    results_df = pd.DataFrame([
        {'term': term, 'count': count}
        for term, count in term_counts.items()
    ])
    results_df = results_df.sort_values('count', ascending=False).reset_index(drop=True)
    
    # Save to CSV
    results_df.to_csv(output_file, index=False)
    print(f"\nResults saved to: {output_file}")
    
    # Display top 10 terms
    print("\nTop 10 most common GO Slim terms:")
    print(results_df.head(10).to_string(index=False))
    
    print("\nDone!")


if __name__ == "__main__":
    main()
