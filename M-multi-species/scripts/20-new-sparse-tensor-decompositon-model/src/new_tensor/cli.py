"""Command-line interface for the tensor decomposition workflow."""

import logging
import json
from pathlib import Path
from typing import Optional

import typer
import pandas as pd
import numpy as np
import yaml

from .io import load_csv, extract_sample_metadata, validate_input_data
from .preprocess import aggregate_replicates, normalize_tensor, filter_genes
from .tensor import build_tensor, save_tensor, save_tensor_mappings
from .model import sparse_cp_decomposition
from .optimize import optimize_hyperparameters, save_optimization_results

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

app = typer.Typer()


def load_config(config_path: str) -> dict:
    """Load configuration from YAML file."""
    with open(config_path, 'r') as f:
        return yaml.safe_load(f)


@app.command("build-tensor")
def build_tensor_cmd(
    input_path: str = typer.Option(..., help="Path to input CSV file"),
    output_dir: str = typer.Option("output/tensor", help="Output directory"),
    aggregation_method: str = typer.Option("median", help="Replicate aggregation method"),
    normalization: str = typer.Option("zscore", help="Normalization method"),
    min_expression: float = typer.Option(1.0, help="Minimum expression threshold"),
    min_variance_percentile: float = typer.Option(10.0, help="Variance percentile threshold"),
    config: Optional[str] = typer.Option(None, help="Path to config YAML file")
):
    """Build tensor from input CSV data."""
    logger.info("Starting tensor construction")

    # Load configuration if provided
    if config:
        config_data = load_config(config)
        # Override defaults with config values
        for key, value in config_data.get('build_tensor', {}).items():
            if key in locals():
                locals()[key] = value

    # Load and validate input data
    df = load_csv(input_path)
    validate_input_data(df)

    # Extract sample metadata
    sample_info, species_list, timepoints_list = extract_sample_metadata(df)

    # Aggregate replicates
    aggregated_df = aggregate_replicates(df, sample_info, method=aggregation_method)

    # Filter genes
    filtered_df = filter_genes(
        aggregated_df,
        min_expression=min_expression,
        min_variance_percentile=min_variance_percentile
    )

    # Normalize data
    normalized_df = normalize_tensor(filtered_df, method=normalization)

    # Build tensor
    tensor, gene_mapping, species_mapping, timepoint_mapping = build_tensor(
        normalized_df, species_list, timepoints_list
    )

    # Save results
    output_path = Path(output_dir)
    save_tensor(tensor, output_dir)
    save_tensor_mappings(
        output_dir, gene_mapping, species_mapping, timepoint_mapping,
        species_list, timepoints_list
    )

    logger.info(f"Tensor saved to {output_path}")


@app.command("optimize")
def optimize_cmd(
    tensor_path: str = typer.Option(..., help="Path to input tensor .npz file"),
    output_dir: str = typer.Option("output/optimization", help="Output directory"),
    n_trials: int = typer.Option(100, help="Number of optimization trials"),
    rank_min: int = typer.Option(2, help="Minimum rank to try"),
    rank_max: int = typer.Option(12, help="Maximum rank to try"),
    lambda_min: float = typer.Option(1e-4, help="Minimum L1 penalty"),
    lambda_max: float = typer.Option(1.0, help="Maximum L1 penalty"),
    config: Optional[str] = typer.Option(None, help="Path to config YAML file")
):
    """Run hyperparameter optimization."""
    logger.info("Starting hyperparameter optimization")

    # Load configuration if provided
    if config:
        config_data = load_config(config)
        # Override defaults with config values
        for key, value in config_data.get('optimize', {}).items():
            if key in locals():
                locals()[key] = value

    # Load tensor
    from .tensor import load_tensor
    # Handle both directory path and file path
    tensor_file = Path(tensor_path)
    if tensor_file.is_file():
        # Full file path provided
        tensor_dir = str(tensor_file.parent)
        tensor_filename = tensor_file.name
    else:
        # Directory path provided
        tensor_dir = tensor_path
        tensor_filename = 'tensor.npz'
    tensor = load_tensor(tensor_dir, tensor_filename)

    # Run optimization
    rank_range = (rank_min, rank_max)
    lambda_bounds = (lambda_min, lambda_max)

    best_params, study = optimize_hyperparameters(
        tensor=tensor,
        n_trials=n_trials,
        rank_range=rank_range,
        lambda_bounds=lambda_bounds,
        verbose=True
    )

    # Save results
    save_optimization_results(study, output_dir, best_params)

    logger.info(f"Optimization results saved to {output_dir}")


@app.command("fit")
def fit_cmd(
    tensor_path: str = typer.Option(..., help="Path to input tensor .npz file"),
    output_dir: str = typer.Option("output/fit", help="Output directory"),
    rank: int = typer.Option(..., help="Rank for decomposition"),
    lambda_a: float = typer.Option(0.1, help="L1 penalty for gene factors"),
    lambda_b: float = typer.Option(0.1, help="L1 penalty for species factors"),
    lambda_c: float = typer.Option(0.1, help="L1 penalty for time factors"),
    non_negative: bool = typer.Option(False, help="Enforce non-negativity"),
    max_iter: int = typer.Option(100, help="Maximum iterations"),
    config: Optional[str] = typer.Option(None, help="Path to config YAML file")
):
    """Fit CP decomposition with specified parameters."""
    logger.info("Starting CP decomposition fit")

    # Load configuration if provided
    if config:
        config_data = load_config(config)
        # Override defaults with config values
        for key, value in config_data.get('fit', {}).items():
            if key in locals():
                locals()[key] = value

    # Load tensor
    from .tensor import load_tensor
    # Handle both directory path and file path
    tensor_file = Path(tensor_path)
    if tensor_file.is_file():
        # Full file path provided
        tensor_dir = str(tensor_file.parent)
        tensor_filename = tensor_file.name
    else:
        # Directory path provided
        tensor_dir = tensor_path
        tensor_filename = 'tensor.npz'
    tensor = load_tensor(tensor_dir, tensor_filename)

    # Fit model
    factors, loss_history, metrics = sparse_cp_decomposition(
        tensor=tensor,
        rank=rank,
        lambda_A=lambda_a,
        lambda_B=lambda_b,
        lambda_C=lambda_c,
        non_negative=non_negative,
        max_iter=max_iter,
        verbose=True
    )

    # Save factors
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # Save factor matrices
    gene_factors = pd.DataFrame(factors[0])
    gene_factors.to_csv(output_path / 'gene_factors.csv', index=False)

    species_factors = pd.DataFrame(factors[1])
    species_factors.to_csv(output_path / 'species_factors.csv', index=False)

    time_factors = pd.DataFrame(factors[2])
    time_factors.to_csv(output_path / 'time_factors.csv', index=False)

    # Save metrics and loss history
    with open(output_path / 'fit_metrics.json', 'w') as f:
        json.dump(metrics, f, indent=2)

    loss_df = pd.DataFrame({'iteration': range(len(loss_history)), 'loss': loss_history})
    loss_df.to_csv(output_path / 'loss_history.csv', index=False)

    logger.info(f"Fit results saved to {output_path}")


@app.command("export")
def export_cmd(
    fit_dir: str = typer.Option(..., help="Directory containing fit results"),
    mappings_dir: str = typer.Option(..., help="Directory containing tensor mappings"),
    output_dir: str = typer.Option("output/export", help="Output directory"),
    plot_heatmaps: bool = typer.Option(True, help="Generate factor heatmaps")
):
    """Export results in various formats."""
    logger.info("Starting export")

    # Load fit results
    fit_path = Path(fit_dir)

    gene_factors = pd.read_csv(fit_path / 'gene_factors.csv')
    species_factors = pd.read_csv(fit_path / 'species_factors.csv')
    time_factors = pd.read_csv(fit_path / 'time_factors.csv')

    with open(fit_path / 'fit_metrics.json', 'r') as f:
        metrics = json.load(f)

    # Load mappings
    mappings_path = Path(mappings_dir)

    gene_df = pd.read_csv(mappings_path / 'genes.csv')
    gene_mapping = dict(zip(gene_df['gene_id'], gene_df['index']))

    species_df = pd.read_csv(mappings_path / 'species.csv')
    species_mapping = dict(zip(species_df['species'], species_df['index']))

    tp_df = pd.read_csv(mappings_path / 'timepoints.csv')
    timepoint_mapping = dict(zip(tp_df['timepoint'], tp_df['index']))

    # Create output directory
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # Export gene factors with gene IDs
    gene_factors_with_ids = gene_factors.copy()
    gene_factors_with_ids.insert(0, 'gene_id', list(gene_mapping.keys()))
    gene_factors_with_ids.to_csv(output_path / 'gene_factors_with_ids.csv', index=False)

    # Export species factors with species codes
    species_factors_with_codes = species_factors.copy()
    species_factors_with_codes.insert(0, 'species', list(species_mapping.keys()))
    species_factors_with_codes.to_csv(output_path / 'species_factors_with_codes.csv', index=False)

    # Export time factors with timepoint labels
    time_factors_with_labels = time_factors.copy()
    time_factors_with_labels.insert(0, 'timepoint', list(timepoint_mapping.keys()))
    time_factors_with_labels.to_csv(output_path / 'time_factors_with_labels.csv', index=False)

    # Save comprehensive results
    results = {
        'metrics': metrics,
        'gene_mapping': gene_mapping,
        'species_mapping': species_mapping,
        'timepoint_mapping': timepoint_mapping,
        'n_components': gene_factors.shape[1]
    }

    with open(output_path / 'decomposition_results.json', 'w') as f:
        json.dump(results, f, indent=2)

    # Generate plots if requested
    if plot_heatmaps:
        try:
            import matplotlib.pyplot as plt
            import seaborn as sns

            # Set up the plotting style
            plt.style.use('default')

            # Gene factors heatmap
            plt.figure(figsize=(12, 8))
            sns.heatmap(gene_factors.values, cmap='viridis', center=0)
            plt.title('Gene Factors')
            plt.xlabel('Components')
            plt.ylabel('Genes')
            plt.tight_layout()
            plt.savefig(output_path / 'gene_factors_heatmap.png', dpi=300, bbox_inches='tight')
            plt.close()

            # Species factors heatmap
            plt.figure(figsize=(10, 6))
            sns.heatmap(species_factors.values, cmap='viridis', center=0, annot=True, fmt='.2f')
            plt.title('Species Factors')
            plt.xlabel('Components')
            plt.ylabel('Species')
            plt.tight_layout()
            plt.savefig(output_path / 'species_factors_heatmap.png', dpi=300, bbox_inches='tight')
            plt.close()

            # Time factors heatmap
            plt.figure(figsize=(8, 6))
            sns.heatmap(time_factors.values, cmap='viridis', center=0, annot=True, fmt='.2f')
            plt.title('Time Factors')
            plt.xlabel('Components')
            plt.ylabel('Timepoints')
            plt.tight_layout()
            plt.savefig(output_path / 'time_factors_heatmap.png', dpi=300, bbox_inches='tight')
            plt.close()

            logger.info("Generated factor heatmaps")

        except ImportError:
            logger.warning("Matplotlib not available, skipping plots")

    logger.info(f"Export completed, results saved to {output_path}")


if __name__ == "__main__":
    app()
