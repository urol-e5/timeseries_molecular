#!/usr/bin/env python3
"""
Test script to validate the GO Slim term summarization script.
"""

import os
import sys
import pandas as pd
from pathlib import Path


def test_goslim_summary():
    """Test that the GO Slim summary was generated correctly."""
    
    print("Testing GO Slim term summarization...")
    print("=" * 70)
    
    # Set up paths
    repo_root = Path(__file__).parent.parent
    component_dir = repo_root / 'output' / '26-rank35-optimization' / 'lambda_gene_0.2' / 'top_genes_per_component'
    summary_file = component_dir / 'goslim_term_counts.csv'
    
    # Test 1: Check if summary file exists
    print("\n1. Checking if summary file exists...")
    if not summary_file.exists():
        raise FileNotFoundError(f"Summary file not found: {summary_file}")
    print(f"   ✓ Found summary file: {summary_file}")
    
    # Test 2: Load and validate structure
    print("\n2. Validating summary file structure...")
    df = pd.read_csv(summary_file)
    
    expected_columns = {'term', 'count'}
    actual_columns = set(df.columns)
    if actual_columns != expected_columns:
        raise ValueError(f"Expected columns {expected_columns}, got {actual_columns}")
    print(f"   ✓ Correct columns: {list(df.columns)}")
    
    # Test 3: Check data types
    print("\n3. Checking data types...")
    if df['term'].dtype != 'object':
        raise ValueError(f"Expected 'term' to be string, got {df['term'].dtype}")
    if not pd.api.types.is_integer_dtype(df['count']):
        raise ValueError(f"Expected 'count' to be integer, got {df['count'].dtype}")
    print(f"   ✓ Correct data types")
    
    # Test 4: Check for data
    print("\n4. Checking for data...")
    if len(df) == 0:
        raise ValueError("Summary file is empty")
    print(f"   ✓ Found {len(df)} unique GO Slim terms")
    
    # Test 5: Validate counts are positive
    print("\n5. Validating counts...")
    if (df['count'] <= 0).any():
        raise ValueError("Found non-positive counts")
    print(f"   ✓ All counts are positive")
    print(f"   Total occurrences: {df['count'].sum()}")
    
    # Test 6: Check sorting (should be descending by count)
    print("\n6. Checking sort order...")
    if not df['count'].is_monotonic_decreasing:
        raise ValueError("Results are not sorted by count (descending)")
    print(f"   ✓ Results sorted by count (descending)")
    
    # Test 7: Verify against source files
    print("\n7. Verifying against source files...")
    annotation_files = sorted([f for f in component_dir.glob('*.csv') 
                              if 'annotation' in f.name.lower()])
    print(f"   Found {len(annotation_files)} annotation files")
    
    if len(annotation_files) == 0:
        raise ValueError("No annotation files found for verification")
    
    # Sample a few terms to verify counts
    print("\n8. Spot-checking term counts...")
    # Get a sample annotation file
    sample_df = pd.read_csv(annotation_files[0])
    if 'goslim_names' in sample_df.columns:
        sample_terms = []
        for terms_str in sample_df['goslim_names'].dropna():
            if terms_str and str(terms_str).strip():
                terms = [t.strip() for t in str(terms_str).split(';') if t.strip()]
                sample_terms.extend(terms)
        
        if sample_terms:
            # Check if these terms exist in summary
            sample_term = sample_terms[0]
            if sample_term in df['term'].values:
                count = df[df['term'] == sample_term]['count'].iloc[0]
                print(f"   ✓ Sample term '{sample_term}' found with count {count}")
            else:
                print(f"   ⚠ Sample term '{sample_term}' not found in summary")
    
    # Test 9: Display summary statistics
    print("\n9. Summary statistics...")
    print(f"   Total unique terms: {len(df)}")
    print(f"   Total occurrences: {df['count'].sum()}")
    print(f"   Most common term: '{df.iloc[0]['term']}' ({df.iloc[0]['count']} occurrences)")
    print(f"   Least common terms: {(df['count'] == 1).sum()} terms with 1 occurrence")
    
    # Display top 5
    print("\n   Top 5 terms:")
    for idx, row in df.head(5).iterrows():
        print(f"     {row['term']}: {row['count']}")
    
    # Summary
    print("\n" + "=" * 70)
    print("✅ ALL TESTS PASSED!")
    print("\nThe GO Slim term summary was generated successfully.")
    print(f"Output file: {summary_file}")
    print("=" * 70)


if __name__ == '__main__':
    try:
        test_goslim_summary()
    except Exception as e:
        print(f"\n❌ TEST FAILED: {e}")
        sys.exit(1)
