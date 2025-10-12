"""Plotting utilities for tensor decomposition results."""

import logging
from typing import List, Dict, Optional, Tuple
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

logger = logging.getLogger(__name__)


def plot_factor_heatmaps(
    gene_factors: np.ndarray,
    species_factors: np.ndarray,
    time_factors: np.ndarray,
    gene_labels: Optional[List[str]] = None,
    species_labels: Optional[List[str]] = None,
    time_labels: Optional[List[str]] = None,
    output_dir: str = "output/plots",
    figsize: Tuple[int, int] = (12, 10),
    cmap: str = 'viridis'
) -> None:
    """Create heatmaps for all factor matrices.

    Args:
        gene_factors: Gene factor matrix (genes × components)
        species_factors: Species factor matrix (species × components)
        time_factors: Time factor matrix (timepoints × components)
        gene_labels: Optional gene labels for y-axis
        species_labels: Optional species labels for y-axis
        time_labels: Optional timepoint labels for y-axis
        output_dir: Output directory for plots
        figsize: Figure size (width, height)
        cmap: Colormap name
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # Set up the plotting style
    plt.style.use('default')

    # Gene factors heatmap
    fig, axes = plt.subplots(2, 2, figsize=figsize)

    # Gene factors
    ax = axes[0, 0]
    sns.heatmap(gene_factors, ax=ax, cmap=cmap, center=0)
    ax.set_title('Gene Factors', fontsize=14, fontweight='bold')
    ax.set_xlabel('Components')
    ax.set_ylabel('Genes' if gene_labels is None else '')

    # Species factors
    ax = axes[0, 1]
    sns.heatmap(species_factors, ax=ax, cmap=cmap, center=0, annot=True, fmt='.2f')
    ax.set_title('Species Factors', fontsize=14, fontweight='bold')
    ax.set_xlabel('Components')
    ax.set_ylabel('Species' if species_labels is None else '')

    # Time factors
    ax = axes[1, 0]
    sns.heatmap(time_factors, ax=ax, cmap=cmap, center=0, annot=True, fmt='.2f')
    ax.set_title('Time Factors', fontsize=14, fontweight='bold')
    ax.set_xlabel('Components')
    ax.set_ylabel('Timepoints' if time_labels is None else '')

    # Component magnitudes
    ax = axes[1, 1]
    component_magnitudes = np.sqrt(
        np.sum(gene_factors**2, axis=0) *
        np.sum(species_factors**2, axis=0) *
        np.sum(time_factors**2, axis=0)
    )
    ax.bar(range(len(component_magnitudes)), component_magnitudes)
    ax.set_title('Component Magnitudes', fontsize=14, fontweight='bold')
    ax.set_xlabel('Component')
    ax.set_ylabel('Magnitude')

    plt.tight_layout()
    plt.savefig(output_path / 'factor_heatmaps.png', dpi=300, bbox_inches='tight')
    plt.close()

    logger.info(f"Saved factor heatmaps to {output_path / 'factor_heatmaps.png'}")


def plot_loss_history(
    loss_history: List[float],
    output_dir: str = "output/plots",
    figsize: Tuple[int, int] = (10, 6)
) -> None:
    """Plot training loss history.

    Args:
        loss_history: List of loss values over iterations
        output_dir: Output directory for plots
        figsize: Figure size (width, height)
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    plt.figure(figsize=figsize)

    plt.plot(loss_history, linewidth=2)
    plt.xlabel('Iteration', fontsize=12)
    plt.ylabel('Loss', fontsize=12)
    plt.title('Training Loss History', fontsize=14, fontweight='bold')
    plt.grid(True, alpha=0.3)

    # Add some styling
    ax = plt.gca()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout()
    plt.savefig(output_path / 'loss_history.png', dpi=300, bbox_inches='tight')
    plt.close()

    logger.info(f"Saved loss history plot to {output_path / 'loss_history.png'}")


def plot_reconstruction_error(
    original_tensor: np.ndarray,
    reconstructed_tensor: np.ndarray,
    output_dir: str = "output/plots",
    figsize: Tuple[int, int] = (12, 8)
) -> None:
    """Plot reconstruction error analysis.

    Args:
        original_tensor: Original tensor
        reconstructed_tensor: Reconstructed tensor
        output_dir: Output directory for plots
        figsize: Figure size (width, height)
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # Compute errors
    error_tensor = original_tensor - reconstructed_tensor
    abs_error = np.abs(error_tensor)

    fig, axes = plt.subplots(2, 3, figsize=figsize)

    # Original tensor slice (first species, all genes and timepoints)
    im1 = axes[0, 0].imshow(original_tensor[:, 0, :].T, cmap='viridis', aspect='auto')
    axes[0, 0].set_title('Original (Species 1)', fontweight='bold')
    axes[0, 0].set_xlabel('Genes')
    axes[0, 0].set_ylabel('Timepoints')
    plt.colorbar(im1, ax=axes[0, 0])

    # Reconstructed tensor slice
    im2 = axes[0, 1].imshow(reconstructed_tensor[:, 0, :].T, cmap='viridis', aspect='auto')
    axes[0, 1].set_title('Reconstructed (Species 1)', fontweight='bold')
    axes[0, 1].set_xlabel('Genes')
    axes[0, 1].set_ylabel('Timepoints')
    plt.colorbar(im2, ax=axes[0, 1])

    # Error heatmap
    im3 = axes[0, 2].imshow(abs_error[:, 0, :].T, cmap='Reds', aspect='auto')
    axes[0, 2].set_title('Absolute Error (Species 1)', fontweight='bold')
    axes[0, 2].set_xlabel('Genes')
    axes[0, 2].set_ylabel('Timepoints')
    plt.colorbar(im3, ax=axes[0, 2])

    # Error distribution
    axes[1, 0].hist(abs_error.flatten(), bins=50, alpha=0.7, edgecolor='black')
    axes[1, 0].set_xlabel('Absolute Error')
    axes[1, 0].set_ylabel('Frequency')
    axes[1, 0].set_title('Error Distribution', fontweight='bold')
    axes[1, 0].grid(True, alpha=0.3)

    # Q-Q plot of errors
    from scipy import stats
    error_flat = abs_error.flatten()
    error_flat = error_flat[~np.isnan(error_flat)]
    if len(error_flat) > 0:
        stats.probplot(error_flat, dist="norm", plot=axes[1, 1])
        axes[1, 1].set_title('Q-Q Plot of Errors', fontweight='bold')
        axes[1, 1].grid(True, alpha=0.3)

    # Remove empty subplot
    axes[1, 2].remove()

    plt.tight_layout()
    plt.savefig(output_path / 'reconstruction_analysis.png', dpi=300, bbox_inches='tight')
    plt.close()

    logger.info(f"Saved reconstruction analysis to {output_path / 'reconstruction_analysis.png'}")


def plot_component_loadings(
    gene_factors: np.ndarray,
    species_factors: np.ndarray,
    time_factors: np.ndarray,
    gene_labels: Optional[List[str]] = None,
    species_labels: Optional[List[str]] = None,
    time_labels: Optional[List[str]] = None,
    output_dir: str = "output/plots",
    top_n: int = 10,
    figsize: Tuple[int, int] = (15, 10)
) -> None:
    """Plot top loadings for each component.

    Args:
        gene_factors: Gene factor matrix
        species_factors: Species factor matrix
        time_factors: Time factor matrix
        gene_labels: Gene labels
        species_labels: Species labels
        time_labels: Timepoint labels
        output_dir: Output directory
        top_n: Number of top loadings to show
        figsize: Figure size
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    n_components = gene_factors.shape[1]

    fig, axes = plt.subplots(2, n_components, figsize=figsize)

    for comp in range(n_components):
        # Gene loadings
        gene_loadings = gene_factors[:, comp]
        top_gene_idx = np.argsort(np.abs(gene_loadings))[-top_n:]

        axes[0, comp].barh(
            range(top_n),
            gene_loadings[top_gene_idx],
            alpha=0.7,
            edgecolor='black'
        )
        axes[0, comp].set_xlabel('Loading')
        axes[0, comp].set_title(f'Component {comp+1}\nGene Loadings', fontweight='bold')
        axes[0, comp].grid(True, alpha=0.3)

        if gene_labels is not None:
            axes[0, comp].set_yticks(range(top_n))
            axes[0, comp].set_yticklabels([gene_labels[i] for i in top_gene_idx])

        # Species loadings
        species_loadings = species_factors[:, comp]

        axes[1, comp].barh(
            range(len(species_loadings)),
            species_loadings,
            alpha=0.7,
            edgecolor='black'
        )
        axes[1, comp].set_xlabel('Loading')
        axes[1, comp].set_title(f'Species Loadings', fontweight='bold')
        axes[1, comp].grid(True, alpha=0.3)

        if species_labels is not None:
            axes[1, comp].set_yticks(range(len(species_loadings)))
            axes[1, comp].set_yticklabels(species_labels)

    plt.tight_layout()
    plt.savefig(output_path / 'component_loadings.png', dpi=300, bbox_inches='tight')
    plt.close()

    logger.info(f"Saved component loadings plot to {output_path / 'component_loadings.png'}")


def plot_optimization_history(
    study: 'optuna.Study',
    output_dir: str = "output/plots",
    figsize: Tuple[int, int] = (12, 8)
) -> None:
    """Plot optimization history and parameter importance.

    Args:
        study: Optuna study object
        output_dir: Output directory
        figsize: Figure size
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # Plot optimization history
    fig, axes = plt.subplots(2, 2, figsize=figsize)

    # Optimization history
    optuna.visualization.matplotlib.plot_optimization_history(study, ax=axes[0, 0])
    axes[0, 0].set_title('Optimization History', fontweight='bold')

    # Parameter importance
    try:
        optuna.visualization.matplotlib.plot_param_importances(study, ax=axes[0, 1])
        axes[0, 1].set_title('Parameter Importance', fontweight='bold')
    except Exception:
        logger.warning("Could not plot parameter importance")

    # Slice plot for rank
    try:
        optuna.visualization.matplotlib.plot_slice(study, params=['rank'], ax=axes[1, 0])
        axes[1, 0].set_title('Rank Slice Plot', fontweight='bold')
    except Exception:
        logger.warning("Could not plot slice for rank")

    # Contour plot for lambda parameters
    try:
        optuna.visualization.matplotlib.plot_contour(study, params=['lambda_A', 'lambda_B'], ax=axes[1, 1])
        axes[1, 1].set_title('Lambda Contour Plot', fontweight='bold')
    except Exception:
        logger.warning("Could not plot contour for lambdas")

    plt.tight_layout()
    plt.savefig(output_path / 'optimization_diagnostics.png', dpi=300, bbox_inches='tight')
    plt.close()

    logger.info(f"Saved optimization diagnostics to {output_path / 'optimization_diagnostics.png'}")
