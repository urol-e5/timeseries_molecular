#!/usr/bin/env python3
"""
Merge ortholog_groups.csv with annotation_with_goslim.tsv
Keep all records from ortholog_groups and join using FUN ID
"""

import pandas as pd
import sys
from pathlib import Path

def main():
    # File paths
    ortholog_file = Path("../11-orthology-analysis/ortholog_groups.csv")
    annotation_file = Path("run_20250831_172744/annotation_with_goslim.tsv")
    output_file = Path("ortholog_groups_annotated.csv")
    
    print(f"Reading ortholog groups from: {ortholog_file}")
    ortholog_df = pd.read_csv(ortholog_file)
    print(f"Ortholog groups shape: {ortholog_df.shape}")
    
    print(f"Reading annotations from: {annotation_file}")
    annotation_df = pd.read_csv(annotation_file, sep='\t')
    print(f"Annotations shape: {annotation_df.shape}")
    
    # The FUN ID is in the 'apul' column of ortholog_groups and 'query' column of annotations
    print("Merging on FUN ID (apul column from ortholog_groups, query column from annotations)")
    
    # Perform left join to keep all records from ortholog_groups
    merged_df = ortholog_df.merge(
        annotation_df, 
        left_on='apul', 
        right_on='query', 
        how='left'
    )
    
    print(f"Merged data shape: {merged_df.shape}")
    
    # Check how many ortholog groups got annotations
    annotated_count = merged_df['query'].notna().sum()
    print(f"Ortholog groups with annotations: {annotated_count}")
    print(f"Ortholog groups without annotations: {len(merged_df) - annotated_count}")
    
    # Save the merged file
    print(f"Saving merged data to: {output_file}")
    merged_df.to_csv(output_file, index=False)
    
    print("Done!")
    
    # Show sample of the merged data
    print("\nSample of merged data:")
    print(merged_df.head())
    
    # Show columns in the merged file
    print(f"\nColumns in merged file ({len(merged_df.columns)} total):")
    for i, col in enumerate(merged_df.columns):
        print(f"{i+1:2d}. {col}")

if __name__ == "__main__":
    main()
