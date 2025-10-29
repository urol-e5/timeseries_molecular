#!/usr/bin/env python3
"""
Annotate rank CSV files with ortholog group annotations.

This script:
1. Reads all CSV files from M-multi-species/output/22-Visualizing-Rank-outs
2. Joins each with ortholog_groups_annotated.csv on OG ID
3. Creates _annotation.csv files with the merged data
4. Generates a summary markdown table with predominant GO processes and gene functions
"""

import os
import sys
import pandas as pd
from pathlib import Path
from collections import Counter
import re

def extract_predominant_terms(go_terms_list, top_n=3):
    """
    Extract the most common GO terms from a list of GO term strings.
    
    Args:
        go_terms_list: List of GO term strings (may contain semicolon-separated values)
        top_n: Number of top terms to return
    
    Returns:
        String with comma-separated top terms
    """
    all_terms = []
    for terms_str in go_terms_list:
        if pd.notna(terms_str) and terms_str.strip():
            # Split by semicolon and clean up each term
            terms = [t.strip() for t in str(terms_str).split(';') if t.strip()]
            all_terms.extend(terms)
    
    if not all_terms:
        return "N/A"
    
    # Count occurrences
    term_counts = Counter(all_terms)
    top_terms = term_counts.most_common(top_n)
    
    # Return as comma-separated string
    return ", ".join([term for term, count in top_terms])

def extract_predominant_functions(protein_names_list, top_n=3):
    """
    Extract the most common gene/protein functions.
    
    Args:
        protein_names_list: List of protein name strings
        top_n: Number of top functions to return
    
    Returns:
        String with comma-separated top functions
    """
    all_functions = []
    for name in protein_names_list:
        if pd.notna(name) and name.strip():
            # Extract the main function name (before parentheses or brackets)
            clean_name = str(name).split('(')[0].strip()
            if clean_name:
                all_functions.append(clean_name)
    
    if not all_functions:
        return "N/A"
    
    # Count occurrences
    func_counts = Counter(all_functions)
    top_funcs = func_counts.most_common(top_n)
    
    # Return as comma-separated string
    return ", ".join([func for func, count in top_funcs])

def annotate_rank_file(rank_file_path, annotation_df, output_dir):
    """
    Annotate a single rank CSV file with ortholog annotations.
    
    Args:
        rank_file_path: Path to the rank CSV file
        annotation_df: DataFrame with ortholog annotations
        output_dir: Directory to save output files
    
    Returns:
        Tuple of (output_file_path, summary_dict)
    """
    # Read the rank file
    rank_df = pd.read_csv(rank_file_path)
    
    # The first column should be the OG ID - get its name
    og_column = rank_df.columns[0]
    
    # Rename to match annotation file's group_id column
    rank_df_renamed = rank_df.rename(columns={og_column: 'group_id'})
    
    # Merge with annotations
    annotated_df = rank_df_renamed.merge(
        annotation_df, 
        on='group_id', 
        how='left'
    )
    
    # Create output filename with _annotation.csv suffix
    input_filename = os.path.basename(rank_file_path)
    output_filename = input_filename.replace('.csv', '_annotation.csv')
    output_path = os.path.join(output_dir, output_filename)
    
    # Save annotated file
    annotated_df.to_csv(output_path, index=False)
    
    # Generate summary statistics
    summary = {
        'filename': input_filename,
        'total_genes': len(annotated_df),
        'annotated_genes': annotated_df['protein_name'].notna().sum(),
        'predominant_go_process': extract_predominant_terms(annotated_df['go_bp'].dropna()),
        'predominant_function': extract_predominant_functions(annotated_df['protein_name'].dropna())
    }
    
    return output_path, summary

def main():
    # Set up paths
    repo_root = Path(__file__).parent.parent.parent
    rank_dir = repo_root / 'M-multi-species' / 'output' / '22-Visualizing-Rank-outs'
    annotation_file = repo_root / 'M-multi-species' / 'output' / '12-ortho-annot' / 'ortholog_groups_annotated.csv'
    output_dir = rank_dir  # Save annotated files in the same directory
    summary_output = repo_root / 'M-multi-species' / 'output' / '22-Visualizing-Rank-outs' / 'annotation_summary.md'
    
    # Check if paths exist
    if not rank_dir.exists():
        print(f"Error: Rank directory not found: {rank_dir}")
        sys.exit(1)
    
    if not annotation_file.exists():
        print(f"Error: Annotation file not found: {annotation_file}")
        sys.exit(1)
    
    # Load the ortholog annotations
    print(f"Loading ortholog annotations from {annotation_file}...")
    annotation_df = pd.read_csv(annotation_file)
    print(f"Loaded {len(annotation_df)} ortholog group annotations")
    
    # Find all CSV files in the rank directory (excluding annotation files)
    csv_files = [f for f in rank_dir.glob('*.csv') 
                 if not f.name.endswith('_annotation.csv') 
                 and f.name != 'annotation_summary.md']
    
    print(f"\nFound {len(csv_files)} CSV files to process")
    
    # Process each file
    summaries = []
    processed = 0
    
    for rank_file in sorted(csv_files):
        try:
            print(f"Processing {rank_file.name}...")
            output_path, summary = annotate_rank_file(rank_file, annotation_df, output_dir)
            summaries.append(summary)
            processed += 1
            
            if processed % 50 == 0:
                print(f"  Processed {processed}/{len(csv_files)} files...")
        except Exception as e:
            print(f"  Error processing {rank_file.name}: {e}")
            continue
    
    print(f"\nSuccessfully processed {processed}/{len(csv_files)} files")
    
    # Create summary markdown table
    print(f"\nGenerating summary table...")
    with open(summary_output, 'w') as f:
        f.write("# Gene List Annotation Summary\n\n")
        f.write(f"Total files processed: {len(summaries)}\n\n")
        f.write("| File | Total Genes | Annotated | Predominant GO Biological Process | Predominant Gene Function |\n")
        f.write("|------|-------------|-----------|-----------------------------------|---------------------------|\n")
        
        for summary in summaries:
            f.write(f"| {summary['filename']} | {summary['total_genes']} | {summary['annotated_genes']} | "
                   f"{summary['predominant_go_process']} | {summary['predominant_function']} |\n")
    
    print(f"Summary table saved to {summary_output}")
    print("\nDone!")

if __name__ == "__main__":
    main()
