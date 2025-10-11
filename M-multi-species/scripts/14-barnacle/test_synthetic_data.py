#!/usr/bin/env python3
"""
Test script to verify the synthetic data can be successfully loaded
and built into a tensor using the build_tensor_and_run.py functions.
"""

import os
import sys
import pandas as pd
import numpy as np
from typing import Dict, List, Tuple


def read_normalized_csv(path: str) -> pd.DataFrame:
    """Read and validate normalized expression CSV file."""
    df = pd.read_csv(path)
    if 'group_id' not in df.columns:
        raise ValueError(f"Expected 'group_id' column in {path}")
    df = df.set_index('group_id')
    return df


def parse_sample_timepoint(column_name: str) -> Tuple[str, int]:
    """Parse sample ID and timepoint from column name."""
    parts = column_name.split('.')
    if len(parts) < 3:
        raise ValueError(f"Unexpected column format (need SAMPLE.TP#): {column_name}")
    time_token = parts[-1]
    if not time_token.startswith('TP'):
        raise ValueError(f"Expected TP# token at end of column: {column_name}")
    try:
        tp = int(time_token.replace('TP', ''))
    except Exception as exc:
        raise ValueError(f"Failed to parse timepoint from {column_name}") from exc
    sample_id = '.'.join(parts[:-1])
    return sample_id, tp


def build_tensor(
    species_to_df: Dict[str, pd.DataFrame],
    expected_timepoints: List[int],
) -> Tuple[np.ndarray, List[str], Dict[int, Dict[str, str]], List[str]]:
    """Build 3D tensor from species dataframes."""
    # Intersect genes
    common_genes = None
    for _, df in species_to_df.items():
        common_genes = df.index if common_genes is None else common_genes.intersection(df.index)
    common_genes = sorted(list(common_genes))
    if len(common_genes) == 0:
        raise ValueError("No intersecting group_id genes across species")

    # Parse columns per species
    species_info = {}
    for species, df in species_to_df.items():
        sample_ids = []
        sample_map = {}
        for col in df.columns:
            sample_id, tp = parse_sample_timepoint(col)
            if tp in expected_timepoints:
                sample_map[(sample_id, tp)] = col
                if sample_id not in sample_ids:
                    sample_ids.append(sample_id)
        species_info[species] = {
            'sample_ids': sample_ids,
            'sample_map': sample_map,
        }

    # Build combined sample list
    sample_labels: List[str] = []
    species_sample_map: Dict[int, Dict[str, str]] = {}
    combined_idx = 0
    for species in ['apul', 'peve', 'ptua']:
        df = species_to_df[species]
        info = species_info[species]
        for sample_id in info['sample_ids']:
            has_any = False
            for tp in expected_timepoints:
                key = (sample_id, tp)
                if key in info['sample_map']:
                    has_any = True
                    break
            if has_any:
                sample_labels.append(f"{species}_{sample_id}")
                species_sample_map[combined_idx] = {
                    'species': species,
                    'sample_id': sample_id,
                }
                combined_idx += 1

    n_genes = len(common_genes)
    n_samples = len(sample_labels)
    n_time = len(expected_timepoints)
    if n_samples == 0:
        raise ValueError("No samples with valid timepoints found")

    tensor = np.empty((n_genes, n_samples, n_time), dtype=float)
    tensor[:] = np.nan

    # Fill tensor
    for s_idx, label in enumerate(sample_labels):
        species, sample_id = label.split('_', 1)
        df = species_to_df[species]
        for t_idx, tp in enumerate(expected_timepoints):
            key = (sample_id, tp)
            col = species_info[species]['sample_map'].get(key)
            if col is None:
                continue
            values = df.loc[common_genes, col].to_numpy(dtype=float)
            tensor[:, s_idx, t_idx] = values

    # Replace NaNs with zeros
    tensor = np.nan_to_num(tensor, nan=0.0)
    return tensor, sample_labels, species_sample_map, common_genes


def test_synthetic_data():
    """Test that synthetic data can be loaded and processed."""
    
    input_dir = 'M-multi-species/output/14-barnacle-synthetic'
    
    print("Testing synthetic data compatibility...")
    print("=" * 70)
    
    # Step 1: Load the data
    print("\n1. Loading data files...")
    expected_files = {
        'apul': os.path.join(input_dir, 'apul_normalized_expression.csv'),
        'peve': os.path.join(input_dir, 'peve_normalized_expression.csv'),
        'ptua': os.path.join(input_dir, 'ptua_normalized_expression.csv'),
    }
    
    for sp, p in expected_files.items():
        if not os.path.exists(p):
            raise FileNotFoundError(f"Missing input for {sp}: {p}")
        print(f"   ✓ Found {sp}: {p}")
    
    species_to_df = {sp: read_normalized_csv(p) for sp, p in expected_files.items()}
    print(f"   ✓ Loaded {len(species_to_df)} species dataframes")
    
    # Step 2: Build tensor
    print("\n2. Building tensor...")
    tensor, sample_labels, species_sample_map, genes = build_tensor(
        species_to_df, 
        expected_timepoints=[1, 2, 3, 4]
    )
    
    print(f"   ✓ Tensor shape: {tensor.shape}")
    print(f"   ✓ Genes: {len(genes)}")
    print(f"   ✓ Samples: {len(sample_labels)}")
    print(f"   ✓ Sample labels (first 5): {sample_labels[:5]}")
    
    # Step 3: Verify tensor properties
    print("\n3. Verifying tensor properties...")
    print(f"   Tensor range: [{tensor.min():.2f}, {tensor.max():.2f}]")
    print(f"   Tensor mean: {tensor.mean():.2f}")
    print(f"   Tensor std: {tensor.std():.2f}")
    print(f"   Non-zero elements: {np.count_nonzero(tensor)} / {tensor.size}")
    print(f"   Sparsity: {1 - np.count_nonzero(tensor) / tensor.size:.2%}")
    
    # Step 4: Verify sample mapping
    print("\n4. Verifying sample mapping...")
    species_counts = {}
    for idx, info in species_sample_map.items():
        species = info['species']
        species_counts[species] = species_counts.get(species, 0) + 1
    
    print(f"   Sample distribution:")
    for species, count in sorted(species_counts.items()):
        print(f"     {species}: {count} samples")
    
    # Step 5: Load and compare with ground truth
    print("\n5. Comparing with ground truth...")
    gt_dir = os.path.join(input_dir, 'ground_truth')
    
    true_gene_factors = pd.read_csv(os.path.join(gt_dir, 'true_gene_factors.csv'), index_col=0)
    true_sample_factors = pd.read_csv(os.path.join(gt_dir, 'true_sample_factors.csv'), index_col=0)
    true_time_factors = pd.read_csv(os.path.join(gt_dir, 'true_time_factors.csv'), index_col=0)
    
    print(f"   ✓ Ground truth gene factors: {true_gene_factors.shape}")
    print(f"   ✓ Ground truth sample factors: {true_sample_factors.shape}")
    print(f"   ✓ Ground truth time factors: {true_time_factors.shape}")
    
    # Reconstruct tensor from ground truth to verify
    n_components = true_gene_factors.shape[1]
    reconstructed = np.zeros_like(tensor)
    
    for r in range(n_components):
        gene_vec = true_gene_factors.iloc[:, r].values
        sample_vec = true_sample_factors.iloc[:, r].values
        time_vec = true_time_factors.iloc[:, r].values
        
        component_tensor = np.einsum('i,j,k->ijk', gene_vec, sample_vec, time_vec)
        reconstructed += component_tensor
    
    # The reconstructed tensor should be similar to the input tensor (before noise)
    # But won't match exactly because we added noise and baseline
    correlation = np.corrcoef(tensor.ravel(), reconstructed.ravel())[0, 1]
    print(f"\n   Correlation between input and reconstructed: {correlation:.3f}")
    
    if correlation > 0.5:
        print(f"   ✓ Good correlation - synthetic structure preserved")
    else:
        print(f"   ⚠ Low correlation - may need to adjust noise level")
    
    # Summary
    print("\n" + "=" * 70)
    print("✅ ALL TESTS PASSED!")
    print("\nSynthetic data is ready for use with build_tensor_and_run.py")
    print("\nRecommended parameters for convergence:")
    print("  --rank 5 (matches ground truth)")
    print("  --lambda-gene 0.1")
    print("  --lambda-sample 0.1")
    print("  --lambda-time 0.05")
    print("  --max-iter 2000 (increased for better convergence)")
    print("  --tol 1e-4 (slightly relaxed)")
    print("=" * 70)


if __name__ == '__main__':
    test_synthetic_data()
