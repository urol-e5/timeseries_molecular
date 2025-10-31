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

def generate_physiological_summary(go_bp_list, protein_names_list):
    """
    Generate a descriptive phrase summarizing the physiological processes.
    
    Analyzes GO biological processes and protein functions to create a concise
    summary of the predominant physiological themes in the gene list.
    
    Args:
        go_bp_list: List of GO biological process strings
        protein_names_list: List of protein names
    
    Returns:
        String describing the predominant physiological processes
    """
    if not go_bp_list or len(go_bp_list) == 0:
        return "Unannotated genes"
    
    # Extract all GO terms
    all_go_terms = []
    for terms_str in go_bp_list:
        if pd.notna(terms_str) and terms_str.strip():
            terms = [t.strip() for t in str(terms_str).split(';') if t.strip()]
            all_go_terms.extend(terms)
    
    if not all_go_terms:
        return "Unannotated genes"
    
    # Keywords for categorizing physiological processes
    process_keywords = {
        'transcription': ['transcription', 'gene expression', 'chromatin'],
        'cell_division': ['cell division', 'mitosis', 'cell cycle', 'cytokinesis'],
        'development': ['development', 'differentiation', 'morphogenesis', 'embryonic'],
        'metabolism': ['metabolic', 'biosynthetic', 'catabolic', 'metabolism'],
        'signaling': ['signal transduction', 'signaling pathway', 'receptor signaling'],
        'immune': ['immune', 'defense response', 'innate immune', 'inflammatory'],
        'transport': ['transport', 'localization', 'secretion', 'endocytosis'],
        'protein_processing': ['protein folding', 'proteolysis', 'protein modification', 'ubiquitin'],
        'cell_structure': ['cytoskeleton', 'cell adhesion', 'extracellular matrix'],
        'reproduction': ['spermatogenesis', 'fertilization', 'gamete', 'meiosis'],
        'cilium': ['cilium', 'flagell', 'axoneme'],
        'apoptosis': ['apoptosis', 'programmed cell death', 'cell death'],
        'stress': ['stress response', 'heat shock', 'oxidative stress']
    }
    
    # Count categories
    category_counts = Counter()
    for term in all_go_terms:
        term_lower = term.lower()
        for category, keywords in process_keywords.items():
            if any(keyword in term_lower for keyword in keywords):
                category_counts[category] += 1
    
    if not category_counts:
        # Fallback: extract key words from most common terms
        term_counts = Counter(all_go_terms)
        top_term = term_counts.most_common(1)[0][0] if term_counts else ""
        # Extract the main process name (before bracketed GO ID)
        main_process = top_term.split('[')[0].strip()
        if main_process:
            return f"Genes involved in {main_process.lower()}"
        return "Mixed physiological processes"
    
    # Generate summary based on top categories
    top_categories = category_counts.most_common(3)
    
    # Map category names to readable descriptions
    category_descriptions = {
        'transcription': 'gene expression regulation',
        'cell_division': 'cell division and proliferation',
        'development': 'development and differentiation',
        'metabolism': 'metabolic processes',
        'signaling': 'signal transduction',
        'immune': 'immune and defense responses',
        'transport': 'intracellular transport',
        'protein_processing': 'protein processing and degradation',
        'cell_structure': 'cytoskeletal organization',
        'reproduction': 'reproduction and gametogenesis',
        'cilium': 'cilium assembly and function',
        'apoptosis': 'apoptosis and cell death',
        'stress': 'stress response'
    }
    
    # Build summary
    if len(top_categories) == 1:
        cat = top_categories[0][0]
        return f"Genes primarily involved in {category_descriptions.get(cat, cat)}"
    elif len(top_categories) >= 2:
        primary = category_descriptions.get(top_categories[0][0], top_categories[0][0])
        secondary = category_descriptions.get(top_categories[1][0], top_categories[1][0])
        return f"Genes involved in {primary} and {secondary}"
    
    return "Mixed physiological processes"

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
    
    # The first column should be the OG ID - get its name and validate
    og_column = rank_df.columns[0]
    
    # Validate that the first column looks like an OG ID column
    if not (og_column.upper() in ['OG', 'GROUP_ID'] or 
            og_column.startswith('OG_') or 
            'OG' in og_column.upper()):
        # Check if first few values look like OG IDs
        sample_values = rank_df[og_column].head(3).astype(str)
        if not sample_values.str.contains('OG_', case=False).any():
            print(f"  Warning: First column '{og_column}' may not be an OG ID column")
    
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
        'predominant_function': extract_predominant_functions(annotated_df['protein_name'].dropna()),
        'physiological_summary': generate_physiological_summary(
            annotated_df['go_bp'].dropna().tolist(),
            annotated_df['protein_name'].dropna().tolist()
        )
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
                 if not f.name.endswith('_annotation.csv')]
    
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
        f.write("| File | Total Genes | Annotated | Predominant GO Biological Process | Predominant Gene Function | Physiological Summary |\n")
        f.write("|------|-------------|-----------|-----------------------------------|---------------------------|-----------------------|\n")
        
        for summary in summaries:
            f.write(f"| {summary['filename']} | {summary['total_genes']} | {summary['annotated_genes']} | "
                   f"{summary['predominant_go_process']} | {summary['predominant_function']} | "
                   f"{summary['physiological_summary']} |\n")
    
    print(f"Summary table saved to {summary_output}")
    print("\nDone!")

if __name__ == "__main__":
    main()
