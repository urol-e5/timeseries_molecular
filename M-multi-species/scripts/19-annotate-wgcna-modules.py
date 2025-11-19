#!/usr/bin/env python3
"""
Annotate WGCNA modules with functional information.

This script:
1. Joins wgcna_ortholog_module_assignments.csv with ortholog_groups_annotated.csv on OG ID
2. Creates a summary annotation report characterizing dominant GO Process, Gene function, 
   and physiological pathways for each WGCNA module
3. Outputs summary files to M-multi-species/output/18-ortholog-wgcna/
"""

import os
import sys
import pandas as pd
from pathlib import Path
from collections import Counter


def load_and_join_data(wgcna_file, annotation_file):
    """
    Load WGCNA module assignments and ortholog annotations, then join them.
    
    Args:
        wgcna_file: Path to wgcna_ortholog_module_assignments.csv
        annotation_file: Path to ortholog_groups_annotated.csv
    
    Returns:
        DataFrame with joined data
    """
    print("Loading data files...")
    
    # Load WGCNA module assignments
    wgcna_df = pd.read_csv(wgcna_file)
    print(f"  Loaded {len(wgcna_df)} WGCNA module assignments")
    
    # Load ortholog annotations
    annot_df = pd.read_csv(annotation_file)
    print(f"  Loaded {len(annot_df)} ortholog annotations")
    
    # Join on group_id
    merged_df = wgcna_df.merge(annot_df, on='group_id', how='left')
    print(f"  Merged dataset has {len(merged_df)} rows")
    
    return merged_df


def extract_terms_from_column(series):
    """
    Extract and count terms from a semicolon-separated column.
    
    Args:
        series: Pandas Series containing semicolon-separated terms
    
    Returns:
        Counter object with term counts
    """
    all_terms = []
    
    for value in series.dropna():
        value_str = str(value).strip()
        if value_str and value_str != 'nan':
            # Split by semicolon and clean up each term
            terms = [t.strip() for t in value_str.split(';') if t.strip()]
            all_terms.extend(terms)
    
    return Counter(all_terms)


def summarize_module_annotations(merged_df):
    """
    Create summary statistics for each WGCNA module.
    
    Args:
        merged_df: DataFrame with joined WGCNA and annotation data
    
    Returns:
        Dictionary with module summaries
    """
    print("\nSummarizing annotations for each module...")
    
    module_summaries = {}
    
    # Get unique modules
    modules = sorted(merged_df['wgcna_module'].unique())
    print(f"  Found {len(modules)} WGCNA modules: {modules}")
    
    for module in modules:
        print(f"\n  Processing module {module}...")
        module_data = merged_df[merged_df['wgcna_module'] == module]
        
        summary = {
            'module': module,
            'num_orthologs': len(module_data),
            'num_annotated': module_data['protein_name'].notna().sum(),
        }
        
        # GO Biological Process
        go_bp_counts = extract_terms_from_column(module_data['go_bp'])
        summary['go_bp_top10'] = go_bp_counts.most_common(10)
        summary['go_bp_total_terms'] = len(go_bp_counts)
        
        # GO Cellular Component
        go_cc_counts = extract_terms_from_column(module_data['go_cc'])
        summary['go_cc_top10'] = go_cc_counts.most_common(10)
        summary['go_cc_total_terms'] = len(go_cc_counts)
        
        # GO Molecular Function
        go_mf_counts = extract_terms_from_column(module_data['go_mf'])
        summary['go_mf_top10'] = go_mf_counts.most_common(10)
        summary['go_mf_total_terms'] = len(go_mf_counts)
        
        # GO Slim terms
        goslim_counts = extract_terms_from_column(module_data['goslim_names'])
        summary['goslim_top10'] = goslim_counts.most_common(10)
        summary['goslim_total_terms'] = len(goslim_counts)
        
        # Protein names (for gene function)
        protein_names = module_data['protein_name'].dropna()
        # Count unique protein names
        protein_counts = Counter([str(name).strip() for name in protein_names if str(name).strip()])
        summary['protein_top10'] = protein_counts.most_common(10)
        summary['protein_total_unique'] = len(protein_counts)
        
        module_summaries[module] = summary
        print(f"    - {summary['num_orthologs']} orthologs, {summary['num_annotated']} with annotations")
    
    return module_summaries


def create_detailed_report(module_summaries, output_dir):
    """
    Create a detailed text report for each module.
    
    Args:
        module_summaries: Dictionary with module summary statistics
        output_dir: Path to output directory
    """
    report_file = output_dir / 'wgcna_module_annotation_summary.txt'
    
    print(f"\nCreating detailed report: {report_file}")
    
    with open(report_file, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("WGCNA MODULE ANNOTATION SUMMARY\n")
        f.write("=" * 80 + "\n\n")
        
        for module in sorted(module_summaries.keys()):
            summary = module_summaries[module]
            
            f.write(f"\n{'=' * 80}\n")
            f.write(f"MODULE {module}\n")
            f.write(f"{'=' * 80}\n\n")
            
            f.write(f"Total orthologs: {summary['num_orthologs']}\n")
            f.write(f"Orthologs with annotations: {summary['num_annotated']}\n")
            f.write(f"Annotation coverage: {summary['num_annotated']/summary['num_orthologs']*100:.1f}%\n\n")
            
            # GO Biological Process
            f.write(f"GO BIOLOGICAL PROCESS (Top 10)\n")
            f.write(f"{'-' * 80}\n")
            if summary['go_bp_top10']:
                for term, count in summary['go_bp_top10']:
                    f.write(f"  {count:4d}  {term}\n")
            else:
                f.write("  No GO BP terms found\n")
            f.write(f"\nTotal unique GO BP terms: {summary['go_bp_total_terms']}\n\n")
            
            # GO Cellular Component
            f.write(f"GO CELLULAR COMPONENT (Top 10)\n")
            f.write(f"{'-' * 80}\n")
            if summary['go_cc_top10']:
                for term, count in summary['go_cc_top10']:
                    f.write(f"  {count:4d}  {term}\n")
            else:
                f.write("  No GO CC terms found\n")
            f.write(f"\nTotal unique GO CC terms: {summary['go_cc_total_terms']}\n\n")
            
            # GO Molecular Function
            f.write(f"GO MOLECULAR FUNCTION (Top 10)\n")
            f.write(f"{'-' * 80}\n")
            if summary['go_mf_top10']:
                for term, count in summary['go_mf_top10']:
                    f.write(f"  {count:4d}  {term}\n")
            else:
                f.write("  No GO MF terms found\n")
            f.write(f"\nTotal unique GO MF terms: {summary['go_mf_total_terms']}\n\n")
            
            # GO Slim (Physiological Pathways)
            f.write(f"GO SLIM TERMS / PHYSIOLOGICAL PATHWAYS (Top 10)\n")
            f.write(f"{'-' * 80}\n")
            if summary['goslim_top10']:
                for term, count in summary['goslim_top10']:
                    f.write(f"  {count:4d}  {term}\n")
            else:
                f.write("  No GO Slim terms found\n")
            f.write(f"\nTotal unique GO Slim terms: {summary['goslim_total_terms']}\n\n")
            
            # Gene Functions (Protein names)
            f.write(f"GENE FUNCTIONS / PROTEIN NAMES (Top 10)\n")
            f.write(f"{'-' * 80}\n")
            if summary['protein_top10']:
                for name, count in summary['protein_top10']:
                    f.write(f"  {count:4d}  {name}\n")
            else:
                f.write("  No protein names found\n")
            f.write(f"\nTotal unique protein names: {summary['protein_total_unique']}\n\n")
    
    print(f"  Report saved to: {report_file}")


def create_csv_summaries(module_summaries, output_dir):
    """
    Create CSV files with term counts for each annotation type across all modules.
    
    Args:
        module_summaries: Dictionary with module summary statistics
        output_dir: Path to output directory
    """
    print("\nCreating CSV summary files...")
    
    # Helper function to create term count matrix
    def create_term_matrix(summaries, term_key, top_n=50):
        """Create a matrix of term counts across modules."""
        # Collect all terms across all modules
        all_terms = set()
        module_term_counts = {}
        
        for module, summary in summaries.items():
            terms_list = summary[term_key]
            term_dict = dict(terms_list)
            module_term_counts[module] = term_dict
            all_terms.update(term_dict.keys())
        
        # Count total occurrences of each term
        term_totals = Counter()
        for term_dict in module_term_counts.values():
            term_totals.update(term_dict)
        
        # Get top N terms by total count
        top_terms = [term for term, _ in term_totals.most_common(top_n)]
        
        # Build matrix
        modules = sorted(summaries.keys())
        data = {'term': top_terms}
        
        for module in modules:
            data[f'module_{module}'] = [
                module_term_counts[module].get(term, 0) for term in top_terms
            ]
        
        data['total'] = [term_totals[term] for term in top_terms]
        
        return pd.DataFrame(data)
    
    # GO Biological Process
    go_bp_df = create_term_matrix(module_summaries, 'go_bp_top10')
    go_bp_file = output_dir / 'wgcna_module_go_bp_summary.csv'
    go_bp_df.to_csv(go_bp_file, index=False)
    print(f"  Saved GO BP summary: {go_bp_file}")
    
    # GO Cellular Component
    go_cc_df = create_term_matrix(module_summaries, 'go_cc_top10')
    go_cc_file = output_dir / 'wgcna_module_go_cc_summary.csv'
    go_cc_df.to_csv(go_cc_file, index=False)
    print(f"  Saved GO CC summary: {go_cc_file}")
    
    # GO Molecular Function
    go_mf_df = create_term_matrix(module_summaries, 'go_mf_top10')
    go_mf_file = output_dir / 'wgcna_module_go_mf_summary.csv'
    go_mf_df.to_csv(go_mf_file, index=False)
    print(f"  Saved GO MF summary: {go_mf_file}")
    
    # GO Slim terms
    goslim_df = create_term_matrix(module_summaries, 'goslim_top10')
    goslim_file = output_dir / 'wgcna_module_goslim_summary.csv'
    goslim_df.to_csv(goslim_file, index=False)
    print(f"  Saved GO Slim summary: {goslim_file}")
    
    # Protein names
    protein_df = create_term_matrix(module_summaries, 'protein_top10')
    protein_file = output_dir / 'wgcna_module_protein_summary.csv'
    protein_df.to_csv(protein_file, index=False)
    print(f"  Saved protein name summary: {protein_file}")
    
    # Create overview summary
    overview_data = {
        'module': [],
        'num_orthologs': [],
        'num_annotated': [],
        'annotation_coverage_%': [],
        'num_go_bp_terms': [],
        'num_go_cc_terms': [],
        'num_go_mf_terms': [],
        'num_goslim_terms': [],
        'num_unique_proteins': [],
        'top_go_bp': [],
        'top_goslim': [],
        'top_protein': []
    }
    
    for module in sorted(module_summaries.keys()):
        summary = module_summaries[module]
        overview_data['module'].append(module)
        overview_data['num_orthologs'].append(summary['num_orthologs'])
        overview_data['num_annotated'].append(summary['num_annotated'])
        overview_data['annotation_coverage_%'].append(
            round(summary['num_annotated']/summary['num_orthologs']*100, 1) if summary['num_orthologs'] > 0 else 0
        )
        overview_data['num_go_bp_terms'].append(summary['go_bp_total_terms'])
        overview_data['num_go_cc_terms'].append(summary['go_cc_total_terms'])
        overview_data['num_go_mf_terms'].append(summary['go_mf_total_terms'])
        overview_data['num_goslim_terms'].append(summary['goslim_total_terms'])
        overview_data['num_unique_proteins'].append(summary['protein_total_unique'])
        
        # Top terms
        overview_data['top_go_bp'].append(
            summary['go_bp_top10'][0][0] if summary['go_bp_top10'] else 'N/A'
        )
        overview_data['top_goslim'].append(
            summary['goslim_top10'][0][0] if summary['goslim_top10'] else 'N/A'
        )
        overview_data['top_protein'].append(
            summary['protein_top10'][0][0] if summary['protein_top10'] else 'N/A'
        )
    
    overview_df = pd.DataFrame(overview_data)
    overview_file = output_dir / 'wgcna_module_overview.csv'
    overview_df.to_csv(overview_file, index=False)
    print(f"  Saved module overview: {overview_file}")


def main():
    """Main execution function."""
    # Set up paths
    repo_root = Path(__file__).parent.parent
    wgcna_file = repo_root / 'output' / '18-ortholog-wgcna' / 'wgcna_ortholog_module_assignments.csv'
    annotation_file = repo_root / 'output' / '12-ortho-annot' / 'ortholog_groups_annotated.csv'
    output_dir = repo_root / 'output' / '18-ortholog-wgcna'
    
    print("=" * 80)
    print("WGCNA MODULE ANNOTATION")
    print("=" * 80)
    
    # Check if input files exist
    if not wgcna_file.exists():
        print(f"Error: WGCNA file not found: {wgcna_file}")
        sys.exit(1)
    
    if not annotation_file.exists():
        print(f"Error: Annotation file not found: {annotation_file}")
        sys.exit(1)
    
    # Load and join data
    merged_df = load_and_join_data(wgcna_file, annotation_file)
    
    # Summarize annotations for each module
    module_summaries = summarize_module_annotations(merged_df)
    
    # Create output files
    create_detailed_report(module_summaries, output_dir)
    create_csv_summaries(module_summaries, output_dir)
    
    print("\n" + "=" * 80)
    print("✅ WGCNA module annotation completed successfully!")
    print(f"\nOutput files saved to: {output_dir}")
    print("=" * 80)


if __name__ == "__main__":
    main()
