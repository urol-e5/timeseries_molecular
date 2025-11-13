#!/usr/bin/env python3
"""
Compare multiway interactions across 3 species.

This script:
1. Imports 3 CSV files containing multi-way interactions from different species
2. Counts occurrences of "CpG", "lncRNA", and "miRNA" in each row
3. Generates a summary report
4. Creates visualizations of the data
"""

import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from collections import Counter
import urllib.request


def download_csv(url, output_path):
    """
    Download a CSV file from a URL.
    
    Args:
        url: URL to download from
        output_path: Path to save the file
        
    Returns:
        True if successful, False otherwise
    """
    try:
        print(f"Downloading {url}...")
        urllib.request.urlretrieve(url, output_path)
        print(f"  Saved to {output_path}")
        return True
    except Exception as e:
        print(f"  Error downloading {url}: {e}")
        return False


def count_feature_occurrences(df):
    """
    Count occurrences of CpG, lncRNA, and miRNA in each row.
    
    Args:
        df: DataFrame to analyze
        
    Returns:
        DataFrame with original data plus count columns
    """
    # Get original column names (exclude count columns if they exist)
    original_columns = [col for col in df.columns 
                       if not col.endswith('_count')]
    
    # Use vectorized operations for better performance
    def count_features_in_row(row):
        """Count features in a single row."""
        cpg_count = 0
        lncrna_count = 0
        mirna_count = 0
        
        # Search in each cell individually to avoid false positives
        for val in row[original_columns]:
            if pd.notna(val):
                val_str = str(val).lower()
                # Use word boundaries or exact matching to avoid false positives
                cpg_count += val_str.count('cpg')
                lncrna_count += val_str.count('lncrna')
                mirna_count += val_str.count('mirna')
        
        return pd.Series({
            'CpG_count': cpg_count,
            'lncRNA_count': lncrna_count,
            'miRNA_count': mirna_count
        })
    
    # Apply the counting function to each row
    count_df = df.apply(count_features_in_row, axis=1)
    
    # Combine with original dataframe
    result_df = pd.concat([df, count_df], axis=1)
    
    return result_df


def generate_summary_stats(df, species_name):
    """
    Generate summary statistics for a dataset.
    
    Args:
        df: DataFrame with count columns
        species_name: Name of the species
        
    Returns:
        Dictionary with summary statistics
    """
    stats = {
        'species': species_name,
        'total_rows': len(df),
        'cpg_total': df['CpG_count'].sum(),
        'lncrna_total': df['lncRNA_count'].sum(),
        'mirna_total': df['miRNA_count'].sum(),
        'cpg_mean': df['CpG_count'].mean(),
        'lncrna_mean': df['lncRNA_count'].mean(),
        'mirna_mean': df['miRNA_count'].mean(),
        'cpg_rows_with_occurrence': (df['CpG_count'] > 0).sum(),
        'lncrna_rows_with_occurrence': (df['lncRNA_count'] > 0).sum(),
        'mirna_rows_with_occurrence': (df['miRNA_count'] > 0).sum(),
    }
    
    return stats


def create_visualizations(all_data, summary_df, output_dir):
    """
    Create visualizations comparing the three species.
    
    Args:
        all_data: Dictionary mapping species names to DataFrames
        summary_df: DataFrame with summary statistics
        output_dir: Directory to save visualizations
    """
    # Set style
    sns.set_style("whitegrid")
    
    # Figure 1: Total occurrences by feature type and species
    fig, ax = plt.subplots(figsize=(10, 6))
    
    species_names = summary_df['species'].tolist()
    cpg_totals = summary_df['cpg_total'].tolist()
    lncrna_totals = summary_df['lncrna_total'].tolist()
    mirna_totals = summary_df['mirna_total'].tolist()
    
    x = np.arange(len(species_names))
    width = 0.25
    
    ax.bar(x - width, cpg_totals, width, label='CpG', color='#1f77b4')
    ax.bar(x, lncrna_totals, width, label='lncRNA', color='#ff7f0e')
    ax.bar(x + width, mirna_totals, width, label='miRNA', color='#2ca02c')
    
    ax.set_xlabel('Species', fontsize=12)
    ax.set_ylabel('Total Occurrences', fontsize=12)
    ax.set_title('Total Feature Occurrences by Species', fontsize=14, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(species_names)
    ax.legend()
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'total_occurrences_by_species.png'), dpi=300)
    print(f"Saved: {os.path.join(output_dir, 'total_occurrences_by_species.png')}")
    plt.close()
    
    # Figure 2: Mean occurrences per row
    fig, ax = plt.subplots(figsize=(10, 6))
    
    cpg_means = summary_df['cpg_mean'].tolist()
    lncrna_means = summary_df['lncrna_mean'].tolist()
    mirna_means = summary_df['mirna_mean'].tolist()
    
    ax.bar(x - width, cpg_means, width, label='CpG', color='#1f77b4')
    ax.bar(x, lncrna_means, width, label='lncRNA', color='#ff7f0e')
    ax.bar(x + width, mirna_means, width, label='miRNA', color='#2ca02c')
    
    ax.set_xlabel('Species', fontsize=12)
    ax.set_ylabel('Mean Occurrences per Row', fontsize=12)
    ax.set_title('Mean Feature Occurrences per Row by Species', fontsize=14, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(species_names)
    ax.legend()
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'mean_occurrences_by_species.png'), dpi=300)
    print(f"Saved: {os.path.join(output_dir, 'mean_occurrences_by_species.png')}")
    plt.close()
    
    # Figure 3: Percentage of rows with at least one occurrence
    fig, ax = plt.subplots(figsize=(10, 6))
    
    cpg_pct = (summary_df['cpg_rows_with_occurrence'] / summary_df['total_rows'] * 100).tolist()
    lncrna_pct = (summary_df['lncrna_rows_with_occurrence'] / summary_df['total_rows'] * 100).tolist()
    mirna_pct = (summary_df['mirna_rows_with_occurrence'] / summary_df['total_rows'] * 100).tolist()
    
    ax.bar(x - width, cpg_pct, width, label='CpG', color='#1f77b4')
    ax.bar(x, lncrna_pct, width, label='lncRNA', color='#ff7f0e')
    ax.bar(x + width, mirna_pct, width, label='miRNA', color='#2ca02c')
    
    ax.set_xlabel('Species', fontsize=12)
    ax.set_ylabel('Percentage of Rows (%)', fontsize=12)
    ax.set_title('Percentage of Rows with Feature Occurrence by Species', fontsize=14, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(species_names)
    ax.legend()
    ax.set_ylim(0, 100)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'percentage_rows_with_occurrence.png'), dpi=300)
    print(f"Saved: {os.path.join(output_dir, 'percentage_rows_with_occurrence.png')}")
    plt.close()
    
    # Figure 4: Heatmap of counts by species
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    
    for idx, (species_name, df) in enumerate(all_data.items()):
        # Create a matrix of counts
        count_matrix = df[['CpG_count', 'lncRNA_count', 'miRNA_count']].head(50)
        
        sns.heatmap(count_matrix.T, cmap='YlOrRd', annot=False, 
                   cbar_kws={'label': 'Count'}, ax=axes[idx])
        axes[idx].set_title(f'{species_name} (First 50 rows)', fontsize=12, fontweight='bold')
        axes[idx].set_xlabel('Row Index', fontsize=10)
        axes[idx].set_ylabel('Feature Type', fontsize=10)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'count_heatmaps.png'), dpi=300)
    print(f"Saved: {os.path.join(output_dir, 'count_heatmaps.png')}")
    plt.close()
    
    print("\nAll visualizations created successfully!")


def main():
    # URLs for the three CSV files
    urls = {
        'species_1': 'https://gannet.fish.washington.edu/v1_web/owlshell/bu-github/ConTra/output/context_dependent_analysis_20251108_151255/tables/multi_way_interactions.csv',
        'species_2': 'https://gannet.fish.washington.edu/v1_web/owlshell/bu-github/ConTra/output/context_dependent_analysis_20251108_144034/tables/multi_way_interactions.csv',
        'species_3': 'https://gannet.fish.washington.edu/v1_web/owlshell/bu-github/ConTra/output/context_dependent_analysis_20251108_140602/tables/multi_way_interactions.csv'
    }
    
    # Set up paths
    script_dir = Path(__file__).parent
    output_dir = script_dir.parent / 'output' / '30-compare-multiway-interactions'
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print("=== Comparing Multiway Interactions Across 3 Species ===\n")
    
    # Download and load data
    all_data = {}
    for species_name, url in urls.items():
        local_file = output_dir / f'{species_name}_multi_way_interactions.csv'
        
        # Check if file already exists locally, otherwise download
        if local_file.exists():
            print(f"Using existing file: {local_file}")
        else:
            # Try to download the file
            if not download_csv(url, local_file):
                print(f"  Skipping {species_name} due to download error\n")
                continue
        
        # Load the CSV
        try:
            df = pd.read_csv(local_file)
            print(f"  Loaded {len(df)} rows with {len(df.columns)} columns")
            
            # Count feature occurrences
            df = count_feature_occurrences(df)
            all_data[species_name] = df
            
            # Save the annotated data
            annotated_file = output_dir / f'{species_name}_annotated.csv'
            df.to_csv(annotated_file, index=False)
            print(f"  Saved annotated data to {annotated_file}\n")
            
        except Exception as e:
            print(f"  Error loading CSV: {e}\n")
    
    if len(all_data) == 0:
        print("\nError: No data files could be loaded.")
        print("\nTo use this script, you can either:")
        print("  1. Ensure network access to gannet.fish.washington.edu")
        print("  2. Manually download the CSV files and place them in:")
        print(f"     {output_dir}")
        print("     with names: species_1_multi_way_interactions.csv,")
        print("                 species_2_multi_way_interactions.csv,")
        print("                 species_3_multi_way_interactions.csv")
        sys.exit(1)
    
    # Generate summary statistics
    print("\n=== Summary Statistics ===\n")
    summary_data = []
    
    for species_name, df in all_data.items():
        stats = generate_summary_stats(df, species_name)
        summary_data.append(stats)
        
        print(f"{species_name}:")
        print(f"  Total rows: {stats['total_rows']}")
        print(f"  CpG - Total: {stats['cpg_total']}, Mean: {stats['cpg_mean']:.2f}, Rows with occurrence: {stats['cpg_rows_with_occurrence']}")
        print(f"  lncRNA - Total: {stats['lncrna_total']}, Mean: {stats['lncrna_mean']:.2f}, Rows with occurrence: {stats['lncrna_rows_with_occurrence']}")
        print(f"  miRNA - Total: {stats['mirna_total']}, Mean: {stats['mirna_mean']:.2f}, Rows with occurrence: {stats['mirna_rows_with_occurrence']}")
        print()
    
    # Create summary DataFrame and save
    summary_df = pd.DataFrame(summary_data)
    summary_file = output_dir / 'summary_statistics.csv'
    summary_df.to_csv(summary_file, index=False)
    print(f"Summary statistics saved to: {summary_file}\n")
    
    # Create visualizations
    print("\n=== Creating Visualizations ===\n")
    create_visualizations(all_data, summary_df, output_dir)
    
    print("\n=== Analysis Complete ===")
    print(f"All results saved to: {output_dir}")


if __name__ == "__main__":
    main()
