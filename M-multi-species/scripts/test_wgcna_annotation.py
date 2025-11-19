#!/usr/bin/env python3
"""
Test script to validate the WGCNA module annotation output.
"""

import os
import sys
import pandas as pd
from pathlib import Path


def test_wgcna_annotation():
    """Test that the WGCNA module annotation was generated correctly."""
    
    print("Testing WGCNA module annotation...")
    print("=" * 70)
    
    # Set up paths
    repo_root = Path(__file__).parent.parent
    output_dir = repo_root / 'output' / '18-ortholog-wgcna'
    
    # Define expected output files
    expected_files = [
        'wgcna_module_annotation_summary.txt',
        'wgcna_module_overview.csv',
        'wgcna_module_go_bp_summary.csv',
        'wgcna_module_go_cc_summary.csv',
        'wgcna_module_go_mf_summary.csv',
        'wgcna_module_goslim_summary.csv',
        'wgcna_module_protein_summary.csv'
    ]
    
    # Test 1: Check if all output files exist
    print("\n1. Checking if all output files exist...")
    for filename in expected_files:
        file_path = output_dir / filename
        if not file_path.exists():
            raise FileNotFoundError(f"Output file not found: {file_path}")
        print(f"   ✓ Found: {filename}")
    
    # Test 2: Validate overview file structure
    print("\n2. Validating overview file structure...")
    overview_file = output_dir / 'wgcna_module_overview.csv'
    overview_df = pd.read_csv(overview_file)
    
    expected_columns = [
        'module', 'num_orthologs', 'num_annotated', 'annotation_coverage_%',
        'num_go_bp_terms', 'num_go_cc_terms', 'num_go_mf_terms', 
        'num_goslim_terms', 'num_unique_proteins', 'top_go_bp', 
        'top_goslim', 'top_protein'
    ]
    
    for col in expected_columns:
        if col not in overview_df.columns:
            raise ValueError(f"Missing expected column: {col}")
    
    print(f"   ✓ All expected columns present")
    print(f"   Total modules: {len(overview_df)}")
    
    # Test 3: Check number of modules
    print("\n3. Checking number of modules...")
    if len(overview_df) != 15:
        raise ValueError(f"Expected 15 modules, found {len(overview_df)}")
    print(f"   ✓ Found 15 WGCNA modules (0-14)")
    
    # Test 4: Validate module IDs
    print("\n4. Validating module IDs...")
    expected_modules = list(range(15))
    actual_modules = sorted(overview_df['module'].tolist())
    if actual_modules != expected_modules:
        raise ValueError(f"Module IDs mismatch. Expected {expected_modules}, got {actual_modules}")
    print(f"   ✓ Module IDs are correct (0-14)")
    
    # Test 5: Check data integrity
    print("\n5. Checking data integrity...")
    
    # Check that num_annotated <= num_orthologs
    if (overview_df['num_annotated'] > overview_df['num_orthologs']).any():
        raise ValueError("num_annotated is greater than num_orthologs for some modules")
    print(f"   ✓ Annotation counts are valid")
    
    # Check that annotation_coverage matches calculation
    calculated_coverage = (overview_df['num_annotated'] / overview_df['num_orthologs'] * 100).round(1)
    if not (calculated_coverage == overview_df['annotation_coverage_%']).all():
        raise ValueError("Annotation coverage percentages don't match calculations")
    print(f"   ✓ Annotation coverage percentages are correct")
    
    # Test 6: Validate CSV summary structures
    print("\n6. Validating CSV summary file structures...")
    
    csv_files = [
        'wgcna_module_go_bp_summary.csv',
        'wgcna_module_go_cc_summary.csv',
        'wgcna_module_go_mf_summary.csv',
        'wgcna_module_goslim_summary.csv',
        'wgcna_module_protein_summary.csv'
    ]
    
    for csv_file in csv_files:
        file_path = output_dir / csv_file
        df = pd.read_csv(file_path)
        
        # Check that 'term' column exists
        if 'term' not in df.columns:
            raise ValueError(f"Missing 'term' column in {csv_file}")
        
        # Check that module columns exist (module_0 through module_14)
        module_cols = [f'module_{i}' for i in range(15)]
        for col in module_cols:
            if col not in df.columns:
                raise ValueError(f"Missing {col} column in {csv_file}")
        
        # Check that 'total' column exists
        if 'total' not in df.columns:
            raise ValueError(f"Missing 'total' column in {csv_file}")
        
        # Verify that totals are sum of module columns
        calculated_total = df[module_cols].sum(axis=1)
        if not (calculated_total == df['total']).all():
            raise ValueError(f"Total column doesn't match sum of modules in {csv_file}")
        
        print(f"   ✓ {csv_file} structure is valid ({len(df)} terms)")
    
    # Test 7: Check text report
    print("\n7. Checking text report...")
    report_file = output_dir / 'wgcna_module_annotation_summary.txt'
    
    with open(report_file, 'r') as f:
        content = f.read()
    
    # Check that all modules are present in the report
    for module_id in range(15):
        if f"MODULE {module_id}" not in content:
            raise ValueError(f"MODULE {module_id} not found in text report")
    
    print(f"   ✓ All 15 modules present in text report")
    print(f"   Report size: {len(content)} characters")
    
    # Test 8: Validate data values
    print("\n8. Validating data values...")
    
    # Check that all modules have orthologs
    if (overview_df['num_orthologs'] <= 0).any():
        raise ValueError("Some modules have 0 orthologs")
    print(f"   ✓ All modules have orthologs")
    
    # Check total number of orthologs
    total_orthologs = overview_df['num_orthologs'].sum()
    print(f"   Total orthologs across all modules: {total_orthologs}")
    
    total_annotated = overview_df['num_annotated'].sum()
    print(f"   Total annotated orthologs: {total_annotated}")
    
    overall_coverage = (total_annotated / total_orthologs * 100)
    print(f"   Overall annotation coverage: {overall_coverage:.1f}%")
    
    # Test 9: Display summary statistics
    print("\n9. Summary statistics by module...")
    print("\n   Top 5 modules by number of orthologs:")
    top_modules = overview_df.nlargest(5, 'num_orthologs')[['module', 'num_orthologs', 'num_annotated', 'annotation_coverage_%']]
    print(top_modules.to_string(index=False))
    
    print("\n   Modules with highest annotation coverage:")
    top_coverage = overview_df.nlargest(5, 'annotation_coverage_%')[['module', 'num_orthologs', 'num_annotated', 'annotation_coverage_%']]
    print(top_coverage.to_string(index=False))
    
    # Summary
    print("\n" + "=" * 70)
    print("✅ ALL TESTS PASSED!")
    print("\nThe WGCNA module annotation was generated successfully.")
    print(f"Output directory: {output_dir}")
    print("=" * 70)


if __name__ == '__main__':
    try:
        test_wgcna_annotation()
    except Exception as e:
        print(f"\n❌ TEST FAILED: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
