#!/usr/bin/env python3
"""
Generate species-level DNA methylation bedgraph files from filtered WGBS CpG count data.

This script processes the existing filtered WGBS CpG count CSV files for each coral species
(D-Apul, E-Peve, F-Ptua) and creates single bedgraph files containing mean methylation
levels (0-100%) for all CpG loci.

Output format: chromosome, start, end, mean_methylation
"""

import pandas as pd
import numpy as np
import os
import sys
import re
from pathlib import Path

def parse_cpg_id(cpg_id):
    """
    Parse CpG ID to extract chromosome and position.
    
    Expected formats:
    - "CpG_chromosome_position" (e.g., "CpG_ntLink_1_22764")
    - "CpG_scaffold_position" (e.g., "CpG_Pocillopora_meandrina_HIv1___Sc0000000_13194")
    
    Returns:
        tuple: (chromosome, position) or None if parsing fails
    """
    try:
        # Remove the "CpG_" prefix
        if cpg_id.startswith('"CpG_'):
            cpg_id = cpg_id[5:-1]  # Remove quotes and CpG_ prefix
        elif cpg_id.startswith('CpG_'):
            cpg_id = cpg_id[4:]  # Remove CpG_ prefix
        
        # Split on underscores and get the last part as position
        parts = cpg_id.split('_')
        if len(parts) >= 2:
            position = int(parts[-1])
            chromosome = '_'.join(parts[:-1])
            return chromosome, position
        else:
            return None
    except (ValueError, IndexError) as e:
        print(f"Warning: Could not parse CpG ID '{cpg_id}': {e}")
        return None

def process_species_data(csv_file, species_name, output_dir):
    """
    Process a species' filtered WGBS CpG count CSV file and generate bedgraph.
    
    Args:
        csv_file (str): Path to the CSV file
        species_name (str): Species name for output file
        output_dir (str): Output directory path
    """
    print(f"Processing {species_name} data from {csv_file}")
    
    # Read the CSV file
    try:
        df = pd.read_csv(csv_file)
        print(f"Loaded data with {df.shape[0]} CpG sites and {df.shape[1]-1} samples")
    except Exception as e:
        print(f"Error reading {csv_file}: {e}")
        return
    
    # First column should be CpG IDs
    cpg_column = df.columns[0]
    sample_columns = df.columns[1:]
    
    # Parse CpG IDs and calculate mean methylation
    bedgraph_data = []
    skipped_count = 0
    
    for idx, row in df.iterrows():
        cpg_id = row[cpg_column]
        
        # Parse chromosome and position
        parsed = parse_cpg_id(cpg_id)
        if parsed is None:
            skipped_count += 1
            continue
            
        chromosome, position = parsed
        
        # Calculate mean methylation across samples (excluding NaN values)
        methylation_values = row[sample_columns].astype(float)
        mean_methylation = methylation_values.mean()
        
        # Skip if all values are NaN
        if pd.isna(mean_methylation):
            skipped_count += 1
            continue
        
        # BedGraph format: chromosome, start, end, value
        # For single nucleotide positions, end = start + 1
        bedgraph_data.append([chromosome, position, position + 1, mean_methylation])
    
    print(f"Processed {len(bedgraph_data)} CpG sites, skipped {skipped_count}")
    
    # Convert to DataFrame and sort by chromosome and position
    bedgraph_df = pd.DataFrame(bedgraph_data, columns=['chromosome', 'start', 'end', 'methylation'])
    bedgraph_df = bedgraph_df.sort_values(['chromosome', 'start'])
    
    # Write bedgraph file
    output_file = os.path.join(output_dir, f"{species_name}_methylation.bedgraph")
    with open(output_file, 'w') as f:
        # Write header
        f.write(f"track type=bedGraph name=\"{species_name} DNA Methylation\" description=\"Mean CpG methylation levels (0-100%)\"\n")
        
        # Write data
        for _, row in bedgraph_df.iterrows():
            f.write(f"{row['chromosome']}\t{int(row['start'])}\t{int(row['end'])}\t{row['methylation']:.2f}\n")
    
    print(f"Generated {output_file}")
    print(f"Bedgraph contains {len(bedgraph_df)} CpG sites")
    print(f"Methylation range: {bedgraph_df['methylation'].min():.2f}% - {bedgraph_df['methylation'].max():.2f}%")
    print()

def main():
    """Main function to process all three species."""
    
    # Define file paths for each species
    base_dir = "/home/runner/work/timeseries_molecular/timeseries_molecular"
    
    species_files = {
        "D-Apul": os.path.join(base_dir, "D-Apul", "data", "count-matrices-01", "Apul-filtered-WGBS-CpG-counts.csv"),
        "F-Ptua": os.path.join(base_dir, "F-Ptua", "output", "05-Ptua-bismark-CG", "filtered-WGBS-CpG-counts.csv"),
    }
    
    # Check for E-Peve files in multiple possible locations
    e_peve_possible_files = [
        os.path.join(base_dir, "E-Peve", "output", "03-Peve-bismark", "merged-WGBS-CpG-counts_filtered.csv"),
        os.path.join(base_dir, "E-Peve", "output", "03-Peve-bismark", "filtered-WGBS-CpG-counts.csv"),
        os.path.join(base_dir, "E-Peve", "data", "filtered-WGBS-CpG-counts.csv"),
    ]
    
    for e_peve_file in e_peve_possible_files:
        if os.path.exists(e_peve_file):
            species_files["E-Peve"] = e_peve_file
            break
    
    # Create output directory
    output_dir = os.path.join(base_dir, "methylation_bedgraphs")
    os.makedirs(output_dir, exist_ok=True)
    
    print(f"Output directory: {output_dir}")
    print("="*60)
    
    # Process each species
    for species, file_path in species_files.items():
        if os.path.exists(file_path):
            process_species_data(file_path, species, output_dir)
        else:
            print(f"Warning: {species} data file not found at {file_path}")
            print()
    
    # Check if E-Peve wasn't processed
    if "E-Peve" not in [k for k, v in species_files.items() if os.path.exists(v)]:
        print("NOTE: E-Peve data not found. You may need to run the processing script first:")
        print("cd E-Peve/output/03-Peve-bismark && bash ../../code/11.1-Peve-WGBS-merge-cov-files.sh")
        print()
    
    print("Bedgraph generation complete!")

if __name__ == "__main__":
    main()