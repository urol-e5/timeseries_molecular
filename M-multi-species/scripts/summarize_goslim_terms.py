#!/usr/bin/env python3
"""
Summarize GO slim terms across all component annotation files.

This script:
1. Reads all CSV files with "annotation" in their filename from the target directory
2. Extracts and counts terms from the 'goslim_names' column
3. Creates a summary CSV where:
   - Rows are unique GO slim terms
   - Columns are the component files (Component_1, Component_2, etc.)
   - Values are the count of occurrences of each term in each component
"""

import os
import pandas as pd
from pathlib import Path
from collections import defaultdict
import re

def extract_component_number(filename):
    """Extract component number from filename like 'Component_10_top100_annotation.csv'"""
    match = re.search(r'Component_(\d+)', filename)
    if match:
        return int(match.group(1))
    return None

def parse_goslim_terms(goslim_str):
    """
    Parse GO slim terms from a string.
    
    Terms are separated by semicolons.
    Some terms may have GO IDs in brackets which we'll remove.
    
    Args:
        goslim_str: String containing GO slim terms
    
    Returns:
        List of cleaned terms
    """
    if pd.isna(goslim_str) or not goslim_str.strip():
        return []
    
    # Split by semicolon
    terms = [t.strip() for t in str(goslim_str).split(';')]
    
    # Remove GO IDs in brackets and clean up
    cleaned_terms = []
    for term in terms:
        # Remove GO IDs like [GO:0005509]
        term = re.sub(r'\s*\[GO:\d+\]', '', term)
        term = term.strip()
        if term:
            cleaned_terms.append(term)
    
    return cleaned_terms

def main():
    # Define paths
    script_dir = Path(__file__).parent
    repo_root = script_dir.parent.parent
    input_dir = repo_root / "M-multi-species" / "output" / "26-rank35-optimization" / "lambda_gene_0.2" / "top_genes_per_component"
    output_file = input_dir / "goslim_summary.csv"
    
    # Find all annotation files
    annotation_files = sorted([
        f for f in input_dir.glob("*annotation*.csv")
        if f.name.startswith("Component_")
    ], key=lambda f: extract_component_number(f.name))
    
    print(f"Found {len(annotation_files)} annotation files")
    
    if len(annotation_files) == 0:
        print("No annotation files found!")
        return
    
    # Dictionary to store term counts: {term: {component: count}}
    term_counts = defaultdict(lambda: defaultdict(int))
    
    # Process each file
    for file_path in annotation_files:
        component_num = extract_component_number(file_path.name)
        component_name = f"Component_{component_num}"
        
        print(f"Processing {file_path.name}...")
        
        # Read the CSV file
        df = pd.read_csv(file_path)
        
        # Extract and count GO slim terms
        for goslim_str in df['goslim_names']:
            terms = parse_goslim_terms(goslim_str)
            for term in terms:
                term_counts[term][component_name] += 1
    
    # Create a DataFrame from the term counts
    # Get all unique terms (sorted alphabetically)
    all_terms = sorted(term_counts.keys())
    
    # Get all component names (sorted by component number)
    all_components = [f"Component_{extract_component_number(f.name)}" for f in annotation_files]
    
    # Build the data matrix
    data = []
    for term in all_terms:
        row = {'goslim_term': term}
        for component in all_components:
            row[component] = term_counts[term].get(component, 0)
        data.append(row)
    
    # Create DataFrame
    summary_df = pd.DataFrame(data)
    
    # Save to CSV
    summary_df.to_csv(output_file, index=False)
    
    print(f"\nSummary statistics:")
    print(f"Total unique GO slim terms: {len(all_terms)}")
    print(f"Total components: {len(all_components)}")
    print(f"Output saved to: {output_file}")
    
    # Show first few rows
    print(f"\nFirst 10 rows of summary:")
    print(summary_df.head(10).to_string())

if __name__ == "__main__":
    main()
