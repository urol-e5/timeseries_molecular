#!/usr/bin/env python3
"""
Generate synthetic gene expression data for barnacle tensor decomposition.

This script creates three species datasets (apul, peve, ptua) with structured
temporal patterns designed to facilitate convergence of the sparse CP decomposition.

The synthetic data includes:
- Common genes across all species (ortholog groups)
- Multiple samples per species
- Temporal patterns with 4 timepoints
- Realistic noise levels
- Structured patterns that should converge well
"""

import argparse
import os
from typing import List, Tuple

import numpy as np
import pandas as pd


def generate_temporal_components(
    n_components: int = 5,
    n_timepoints: int = 4,
    seed: int = 42,
) -> np.ndarray:
    """
    Generate temporal loading patterns for synthetic components.
    
    These patterns are designed to be distinct and identifiable,
    which helps with convergence.
    """
    np.random.seed(seed)
    
    time_factors = np.zeros((n_timepoints, n_components))
    
    # Define patterns based on number of components
    patterns = []
    
    # Pattern 1: Linear increase
    patterns.append(np.linspace(0.2, 1.0, n_timepoints))
    
    # Pattern 2: Linear decrease
    patterns.append(np.linspace(1.0, 0.2, n_timepoints))
    
    # Pattern 3: Peak at middle timepoints
    if n_timepoints == 4:
        patterns.append(np.array([0.3, 0.8, 0.9, 0.4]))
    else:
        # General middle peak
        mid = n_timepoints // 2
        pattern = np.zeros(n_timepoints)
        pattern[:mid] = np.linspace(0.3, 0.9, mid)
        pattern[mid:] = np.linspace(0.9, 0.4, n_timepoints - mid)
        patterns.append(pattern)
    
    # Pattern 4: Early peak
    if n_timepoints == 4:
        patterns.append(np.array([1.0, 0.7, 0.4, 0.3]))
    else:
        patterns.append(np.linspace(1.0, 0.3, n_timepoints))
    
    # Pattern 5: Late peak
    if n_timepoints == 4:
        patterns.append(np.array([0.3, 0.4, 0.7, 1.0]))
    else:
        patterns.append(np.linspace(0.3, 1.0, n_timepoints))
    
    # Assign patterns to components (cycle if more components than patterns)
    for i in range(n_components):
        pattern_idx = i % len(patterns)
        time_factors[:, i] = patterns[pattern_idx]
        # Add small random variation to distinguish repeated patterns
        if i >= len(patterns):
            time_factors[:, i] += np.random.normal(0, 0.05, n_timepoints)
            time_factors[:, i] = np.maximum(time_factors[:, i], 0.1)  # Keep positive
    
    return time_factors


def generate_gene_factors(
    n_genes: int,
    n_components: int = 5,
    sparsity: float = 0.7,
    seed: int = 42,
) -> np.ndarray:
    """
    Generate gene loading matrix with controlled sparsity.
    
    Most genes are zero for most components (sparse),
    which is realistic and helps convergence.
    """
    np.random.seed(seed)
    
    gene_factors = np.zeros((n_genes, n_components))
    
    # For each component, select a subset of genes
    genes_per_component = int(n_genes * (1 - sparsity))
    
    for comp in range(n_components):
        # Select random genes for this component
        active_genes = np.random.choice(n_genes, genes_per_component, replace=False)
        # Assign positive values from gamma distribution (realistic for expression)
        gene_factors[active_genes, comp] = np.random.gamma(2, 2, genes_per_component)
    
    return gene_factors


def generate_sample_factors(
    n_samples: int,
    n_components: int = 5,
    seed: int = 42,
) -> np.ndarray:
    """
    Generate sample loading matrix.
    
    Different samples have different component weights,
    representing biological variation.
    """
    np.random.seed(seed)
    
    # Use Dirichlet distribution to ensure samples are mixtures of components
    sample_factors = np.random.dirichlet(np.ones(n_components) * 2, n_samples) * 10
    
    return sample_factors


def build_synthetic_tensor(
    n_genes: int,
    n_samples: int,
    n_timepoints: int = 4,
    n_components: int = 5,
    noise_level: float = 0.1,
    seed: int = 42,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Build a synthetic tensor from CP factors with added noise.
    
    Returns:
        tensor: (n_genes, n_samples, n_timepoints)
        gene_factors: (n_genes, n_components)
        sample_factors: (n_samples, n_components)
        time_factors: (n_timepoints, n_components)
    """
    np.random.seed(seed)
    
    # Generate factor matrices
    gene_factors = generate_gene_factors(n_genes, n_components, sparsity=0.7, seed=seed)
    sample_factors = generate_sample_factors(n_samples, n_components, seed=seed + 1)
    time_factors = generate_temporal_components(n_components, n_timepoints, seed=seed + 2)
    
    # Reconstruct tensor using CP decomposition formula
    # tensor[i,j,k] = sum_r gene[i,r] * sample[j,r] * time[k,r]
    tensor = np.zeros((n_genes, n_samples, n_timepoints))
    
    for r in range(n_components):
        # Outer product of the three factor vectors for component r
        component_tensor = np.einsum('i,j,k->ijk',
                                     gene_factors[:, r],
                                     sample_factors[:, r],
                                     time_factors[:, r])
        tensor += component_tensor
    
    # Add Gaussian noise
    noise = np.random.normal(0, noise_level * np.std(tensor), tensor.shape)
    tensor += noise
    
    # Ensure non-negative (like expression data should be)
    tensor = np.maximum(tensor, 0)
    
    # Add a baseline expression level (realistic for log-normalized data)
    tensor += np.random.gamma(1, 1, (n_genes, 1, 1))
    
    return tensor, gene_factors, sample_factors, time_factors


def create_species_dataframe(
    tensor: np.ndarray,
    species_code: str,
    n_samples_per_species: int,
    sample_offset: int,
    gene_ids: List[str],
) -> pd.DataFrame:
    """
    Create a dataframe for one species from tensor slice.
    
    Args:
        tensor: Full tensor (n_genes, total_samples, n_timepoints)
        species_code: 'apul', 'peve', or 'ptua'
        n_samples_per_species: Number of samples for this species
        sample_offset: Starting sample index for this species
        gene_ids: List of gene IDs
    
    Returns:
        DataFrame with group_id as index and columns like SAMPLE.TP#
    """
    n_genes, _, n_timepoints = tensor.shape
    
    # Sample naming conventions per species
    sample_prefixes = {
        'apul': 'ACR',  # Acropora pulchra
        'peve': 'POR',  # Porites evermanni  
        'ptua': 'POC',  # Pocillopora tuahiniensis
    }
    prefix = sample_prefixes[species_code]
    
    # Build column names
    columns = []
    data_rows = []
    
    for sample_idx in range(n_samples_per_species):
        sample_id = f"{prefix}.{100 + sample_idx + sample_offset}"
        for tp in range(1, n_timepoints + 1):
            col_name = f"{sample_id}.TP{tp}"
            columns.append(col_name)
            # Extract data from tensor
            data_rows.append(tensor[:, sample_offset + sample_idx, tp - 1])
    
    # Create dataframe
    data_matrix = np.column_stack(data_rows)
    df = pd.DataFrame(data_matrix, columns=columns, index=gene_ids)
    df.index.name = 'group_id'
    
    return df


def main():
    parser = argparse.ArgumentParser(
        description='Generate synthetic gene expression data for barnacle analysis'
    )
    parser.add_argument(
        '--output-dir',
        required=True,
        help='Directory to save synthetic CSV files'
    )
    parser.add_argument(
        '--n-genes',
        type=int,
        default=10223,
        help='Number of genes (default: 10223 to match real data)'
    )
    parser.add_argument(
        '--n-samples-per-species',
        type=int,
        default=10,
        help='Number of samples per species (default: 10)'
    )
    parser.add_argument(
        '--n-components',
        type=int,
        default=5,
        help='Number of underlying components (default: 5)'
    )
    parser.add_argument(
        '--noise-level',
        type=float,
        default=0.1,
        help='Noise level as fraction of signal std (default: 0.1)'
    )
    parser.add_argument(
        '--seed',
        type=int,
        default=42,
        help='Random seed for reproducibility (default: 42)'
    )
    
    args = parser.parse_args()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Parameters
    n_genes = args.n_genes
    n_samples_per_species = args.n_samples_per_species
    total_samples = 3 * n_samples_per_species  # Three species
    n_timepoints = 4
    n_components = args.n_components
    
    print(f"Generating synthetic data:")
    print(f"  Genes: {n_genes}")
    print(f"  Samples per species: {n_samples_per_species}")
    print(f"  Total samples: {total_samples}")
    print(f"  Timepoints: {n_timepoints}")
    print(f"  Components: {n_components}")
    print(f"  Noise level: {args.noise_level}")
    print(f"  Random seed: {args.seed}")
    
    # Generate gene IDs (ortholog groups)
    gene_ids = [f"OG_{i:05d}" for i in range(1, n_genes + 1)]
    
    # Build synthetic tensor
    tensor, gene_factors, sample_factors, time_factors = build_synthetic_tensor(
        n_genes=n_genes,
        n_samples=total_samples,
        n_timepoints=n_timepoints,
        n_components=n_components,
        noise_level=args.noise_level,
        seed=args.seed,
    )
    
    print(f"\nTensor shape: {tensor.shape}")
    print(f"Tensor range: [{tensor.min():.2f}, {tensor.max():.2f}]")
    print(f"Tensor mean: {tensor.mean():.2f}")
    
    # Create dataframes for each species
    species_configs = [
        ('apul', 0),
        ('peve', n_samples_per_species),
        ('ptua', 2 * n_samples_per_species),
    ]
    
    for species_code, offset in species_configs:
        df = create_species_dataframe(
            tensor=tensor,
            species_code=species_code,
            n_samples_per_species=n_samples_per_species,
            sample_offset=offset,
            gene_ids=gene_ids,
        )
        
        output_path = os.path.join(args.output_dir, f'{species_code}_normalized_expression.csv')
        df.to_csv(output_path)
        print(f"Saved {species_code}: {output_path}")
        print(f"  Shape: {df.shape}")
        print(f"  Columns: {df.columns[:8].tolist()} ...")
    
    # Save ground truth factors for validation
    truth_dir = os.path.join(args.output_dir, 'ground_truth')
    os.makedirs(truth_dir, exist_ok=True)
    
    pd.DataFrame(gene_factors, index=gene_ids).to_csv(
        os.path.join(truth_dir, 'true_gene_factors.csv')
    )
    pd.DataFrame(sample_factors).to_csv(
        os.path.join(truth_dir, 'true_sample_factors.csv')
    )
    pd.DataFrame(time_factors, index=[f'TP{i}' for i in range(1, n_timepoints + 1)]).to_csv(
        os.path.join(truth_dir, 'true_time_factors.csv')
    )
    
    print(f"\nGround truth factors saved to: {truth_dir}")
    print("✅ Synthetic data generation complete!")


if __name__ == '__main__':
    main()
