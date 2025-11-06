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


def count_goslim_terms_per_component(annotation_files):
    """
    Count occurrences of GO Slim terms per component annotation file.
    
    Args:
        annotation_files: List of paths to annotation CSV files
    
    Returns:
        Dictionary mapping component names to Counter objects with term counts
    """
    component_term_counts = {}
    
    for file_path in annotation_files:
        try:
            # Extract component name from filename (e.g., "Component_1" from "Component_1_top100_annotation.csv")
            filename = file_path.stem  # Gets filename without extension
            component_name = filename.split('_top')[0]  # Extract "Component_X" part
            
            # Read the annotation file
            df = pd.read_csv(file_path)
            
            # Check if goslim_names column exists
            if 'goslim_names' not in df.columns:
                print(f"  Warning: 'goslim_names' column not found in {file_path.name}")
                continue
            
            # Extract all terms from goslim_names column for this component
            component_terms = []
            for terms_str in df['goslim_names'].dropna():
                terms_str = str(terms_str)
                if terms_str and terms_str.strip():
                    # Split by semicolon and clean up each term
                    terms = [t.strip() for t in terms_str.split(';') if t.strip()]
                    component_terms.extend(terms)
            
            # Count occurrences for this component
            component_term_counts[component_name] = Counter(component_terms)
        
        except Exception as e:
            print(f"  Error processing {file_path.name}: {e}")
            continue
    
    return component_term_counts


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
    
    # Count GO Slim terms per component
    print("\nCounting GO Slim term occurrences per component...")
    component_term_counts = count_goslim_terms_per_component(annotation_files)
    
    # Get all unique terms across all components
    all_terms = set()
    for term_counter in component_term_counts.values():
        all_terms.update(term_counter.keys())
    
    print(f"Found {len(all_terms)} unique GO Slim terms")
    
    # Create a matrix: rows are terms, columns are components
    # Sort component names numerically (Component_1, Component_2, ..., Component_35)
    component_names = sorted(component_term_counts.keys(), 
                            key=lambda x: int(x.split('_')[1]))
    
    # Build the results dictionary
    results = {'term': []}
    for comp in component_names:
        results[comp] = []
    results['total'] = []
    
    # Sort terms by total count (descending)
    term_totals = Counter()
    for comp_counts in component_term_counts.values():
        term_totals.update(comp_counts)
    
    sorted_terms = [term for term, count in term_totals.most_common()]
    
    # Fill in the matrix
    for term in sorted_terms:
        results['term'].append(term)
        total = 0
        for comp in component_names:
            count = component_term_counts[comp].get(term, 0)
            results[comp].append(count)
            total += count
        results['total'].append(total)
    
    # Create DataFrame
    results_df = pd.DataFrame(results)
    
    # Save to CSV
    results_df.to_csv(output_file, index=False)
    print(f"\nResults saved to: {output_file}")
    
    # Display summary statistics
    print(f"\nTotal term occurrences: {results_df['total'].sum()}")
    print(f"Number of components: {len(component_names)}")
    
    # Display top 10 terms with their totals
    print("\nTop 10 most common GO Slim terms:")
    print(results_df[['term', 'total']].head(10).to_string(index=False))
    
    print("\nDone!")



if __name__ == "__main__":
    main()
