#!/usr/bin/env python3
"""Smoke test for the tensor decomposition package."""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))

import logging
import numpy as np
import pandas as pd
from pathlib import Path

# Import our package modules
from new_tensor.io import load_csv, extract_sample_metadata, validate_input_data
from new_tensor.preprocess import aggregate_replicates, normalize_tensor, filter_genes
from new_tensor.tensor import build_tensor, save_tensor, save_tensor_mappings
from new_tensor.model import sparse_cp_decomposition

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def run_smoke_test():
    """Run a complete smoke test on a small subset of data."""
    logger.info("Starting smoke test...")

    # Load test data (200 genes subset)
    test_file = "test_subset.csv"
    if not Path(test_file).exists():
        logger.error(f"Test file {test_file} not found")
        return False

    # Step 1: Load and validate data
    logger.info("Step 1: Loading data...")
    df = load_csv(test_file)
    validate_input_data(df)
    logger.info(f"Loaded {df.shape[0]} genes with {df.shape[1] - 1} samples")

    # Step 2: Extract metadata
    logger.info("Step 2: Extracting metadata...")
    sample_info, species_list, timepoints_list = extract_sample_metadata(df)
    logger.info(f"Species: {species_list}")
    logger.info(f"Timepoints: {timepoints_list}")

    # Step 3: Aggregate replicates
    logger.info("Step 3: Aggregating replicates...")
    aggregated_df = aggregate_replicates(df, sample_info, method='median')
    logger.info(f"Aggregated to {aggregated_df.shape[1] - 1} species-timepoint combinations")

    # Step 4: Filter genes
    logger.info("Step 4: Filtering genes...")
    filtered_df = filter_genes(aggregated_df, min_expression=1.0, min_variance_percentile=10.0, max_genes=50)
    logger.info(f"Filtered to {filtered_df.shape[0]} genes")

    # Step 5: Normalize data
    logger.info("Step 5: Normalizing data...")
    normalized_df = normalize_tensor(filtered_df, method='zscore')
    logger.info("Data normalized")

    # Step 6: Build tensor
    logger.info("Step 6: Building tensor...")
    tensor, gene_mapping, species_mapping, timepoint_mapping = build_tensor(
        normalized_df, species_list, timepoints_list
    )
    logger.info(f"Tensor shape: {tensor.shape}")

    # Step 7: Save tensor and mappings
    logger.info("Step 7: Saving tensor and mappings...")
    output_dir = "smoke_test_output"
    save_tensor(tensor, output_dir)
    save_tensor_mappings(
        output_dir, gene_mapping, species_mapping, timepoint_mapping,
        species_list, timepoints_list
    )
    logger.info(f"Results saved to {output_dir}")

    # Step 8: Fit sparse CP decomposition
    logger.info("Step 8: Fitting CP decomposition...")
    factors, loss_history, metrics = sparse_cp_decomposition(
        tensor=tensor,
        rank=3,  # Small rank for smoke test
        lambda_A=0.01,
        lambda_B=0.01,
        lambda_C=0.01,
        non_negative=False,
        max_iter=20,  # Few iterations for speed
        verbose=True
    )

    logger.info(f"Final loss: {metrics['final_loss']:.6f}")
    logger.info(f"Explained variance: {metrics['explained_variance']:.4f}")
    logger.info(f"Sparsity A: {metrics['sparsity_A']:.4f}")

    # Step 9: Save results
    logger.info("Step 9: Saving results...")
    import json

    # Save factor matrices
    gene_factors = pd.DataFrame(factors[0])
    gene_factors.to_csv(f"{output_dir}/gene_factors.csv", index=False)

    species_factors = pd.DataFrame(factors[1])
    species_factors.to_csv(f"{output_dir}/species_factors.csv", index=False)

    time_factors = pd.DataFrame(factors[2])
    time_factors.to_csv(f"{output_dir}/time_factors.csv", index=False)

    # Save metrics and loss history
    with open(f"{output_dir}/fit_metrics.json", 'w') as f:
        json.dump(metrics, f, indent=2)

    loss_df = pd.DataFrame({'iteration': range(len(loss_history)), 'loss': loss_history})
    loss_df.to_csv(f"{output_dir}/loss_history.csv", index=False)

    logger.info("Smoke test completed successfully!")
    return True


if __name__ == "__main__":
    success = run_smoke_test()
    sys.exit(0 if success else 1)
