#!/usr/bin/env python3
"""
Alternative tensor construction: (genes × timepoints × samples)
Each sample becomes a separate slice instead of grouping by species.
"""
import argparse
import json
import os
from datetime import datetime
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from barnacle.decomposition import SparseCP


def read_normalized_csv(path: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    if 'group_id' not in df.columns:
        raise ValueError(f"Expected 'group_id' column in {path}")
    df = df.set_index('group_id')
    return df


def parse_sample_timepoint(column_name: str) -> Tuple[str, str, int]:
    tp_index = column_name.rfind('-TP')
    if tp_index == -1:
        raise ValueError(f"Expected TP# token at end of column: {column_name}")

    time_token = column_name[tp_index + 1:]
    if not time_token.startswith('TP'):
        raise ValueError(f"Expected TP# token at end of column: {column_name}")
    try:
        tp = int(time_token.replace('TP', ''))
    except Exception as exc:
        raise ValueError(f"Failed to parse timepoint from {column_name}") from exc

    sample_id = column_name[:tp_index]
    sample_prefix = sample_id.split('-')[0]
    species_map = {'ACR': 'apul', 'POR': 'peve', 'POC': 'ptua'}
    if sample_prefix not in species_map:
        raise ValueError(f"Unknown species prefix in column: {column_name}")

    species = species_map[sample_prefix]
    return species, sample_id, tp


def build_tensor_by_sample(
    df: pd.DataFrame,
    expected_timepoints: List[int],
) -> Tuple[np.ndarray, List[str], Dict[int, Dict[str, str]], List[str]]:
    """
    Build tensor with shape (genes, timepoints, samples).
    Each unique sample becomes a slice, with timepoints as rows.
    """
    common_genes = sorted(list(df.index))

    # Collect all unique samples across all species
    sample_data = {}  # {sample_id: {species: str, columns: {tp: col_name}}}
    
    for col in df.columns:
        try:
            species, sample_id, tp = parse_sample_timepoint(col)
            if tp in expected_timepoints:
                if sample_id not in sample_data:
                    sample_data[sample_id] = {
                        'species': species,
                        'columns': {}
                    }
                sample_data[sample_id]['columns'][tp] = col
        except ValueError:
            continue

    # Sort samples by species, then by sample_id
    sorted_samples = sorted(
        sample_data.keys(),
        key=lambda sid: (sample_data[sid]['species'], sid)
    )

    # Build sample labels and metadata
    sample_labels: List[str] = []
    sample_metadata: Dict[int, Dict[str, str]] = {}
    
    for idx, sample_id in enumerate(sorted_samples):
        species = sample_data[sample_id]['species']
        sample_labels.append(f"{species}_{sample_id}")
        sample_metadata[idx] = {
            'species': species,
            'sample_id': sample_id,
        }

    n_genes = len(common_genes)
    n_time = len(expected_timepoints)
    n_samples = len(sorted_samples)

    if n_samples == 0:
        raise ValueError("No samples with valid timepoints found")

    # Tensor shape: (genes, timepoints, samples)
    tensor = np.empty((n_genes, n_time, n_samples), dtype=float)
    tensor[:] = np.nan

    # Fill tensor
    for s_idx, sample_id in enumerate(sorted_samples):
        columns = sample_data[sample_id]['columns']
        for t_idx, tp in enumerate(expected_timepoints):
            col = columns.get(tp)
            if col is None:
                continue
            values = df.loc[common_genes, col].to_numpy(dtype=float)
            tensor[:, t_idx, s_idx] = values

    # Replace NaNs with zeros
    tensor = np.nan_to_num(tensor, nan=0.0)
    return tensor, sample_labels, sample_metadata, common_genes


def run_sparse_cp(
    tensor: np.ndarray,
    rank: int,
    lambda_gene: float,
    lambda_time: float,
    lambda_sample: float,
    max_iter: int,
    tol: float,
    seed: int,
):
    """
    Note: lambda order now matches tensor dimensions: [gene, time, sample]
    """
    model = SparseCP(
        rank=rank,
        lambdas=[lambda_gene, lambda_time, lambda_sample],
        nonneg_modes=[0],
        n_initializations=3,
        random_state=seed,
        n_iter_max=max_iter,
        tol=tol,
    )
    decomposition = model.fit_transform(tensor, verbose=1)
    return model, decomposition


def save_outputs(
    output_dir: str,
    tensor: np.ndarray,
    sample_labels: List[str],
    sample_metadata: Dict[int, Dict[str, str]],
    genes: List[str],
    model,
    decomposition,
):
    os.makedirs(output_dir, exist_ok=True)
    factors_dir = os.path.join(output_dir, 'barnacle_factors')
    figs_dir = os.path.join(output_dir, 'figures')
    os.makedirs(factors_dir, exist_ok=True)
    os.makedirs(figs_dir, exist_ok=True)

    np.save(os.path.join(output_dir, 'multiomics_tensor.npy'), tensor)

    # Extract factors: order is [genes, timepoints, samples]
    gene_factors = pd.DataFrame(decomposition.factors[0], index=genes)
    time_factors = pd.DataFrame(decomposition.factors[1], 
                               index=[f'TP{t}' for t in range(1, tensor.shape[1] + 1)])
    sample_factors = pd.DataFrame(decomposition.factors[2], index=sample_labels)

    gene_factors.to_csv(os.path.join(factors_dir, 'gene_factors.csv'))
    time_factors.to_csv(os.path.join(factors_dir, 'time_factors.csv'))
    sample_factors.to_csv(os.path.join(factors_dir, 'sample_factors.csv'))

    # Compute component weights
    weights_attr = getattr(decomposition, 'weights', None)
    if weights_attr is not None:
        weights = np.asarray(weights_attr).astype(float).ravel()
        if np.allclose(weights, weights[0]) if len(weights) > 0 else True:
            weights = None
    else:
        weights = None

    if weights is None:
        gene_norms = np.linalg.norm(gene_factors.values, axis=0)
        time_norms = np.linalg.norm(time_factors.values, axis=0)
        sample_norms = np.linalg.norm(sample_factors.values, axis=0)
        weights = gene_norms * time_norms * sample_norms

    pd.DataFrame({'weight': weights}).to_csv(
        os.path.join(factors_dir, 'component_weights.csv'), index=False)

    # Sample mapping
    mapping_rows = []
    for idx, label in enumerate(sample_labels):
        mapping_rows.append({
            'sample_index': idx,
            'label': label,
            'species': sample_metadata[idx]['species'],
            'sample_id': sample_metadata[idx]['sample_id'],
        })
    pd.DataFrame(mapping_rows).to_csv(
        os.path.join(factors_dir, 'sample_mapping.csv'), index=False)

    # Metadata
    raw_loss = getattr(model, 'loss_', None)
    if raw_loss is None:
        final_loss = None
    else:
        try:
            if isinstance(raw_loss, (list, tuple, np.ndarray)) and len(raw_loss) > 0:
                final_loss = float(raw_loss[-1])
            else:
                final_loss = float(raw_loss)
        except Exception:
            final_loss = None

    metadata = {
        'timestamp': datetime.utcnow().isoformat() + 'Z',
        'tensor_shape': list(map(int, tensor.shape)),
        'tensor_structure': '(genes, timepoints, samples)',
        'n_components': int(gene_factors.shape[1]),
        'model_converged': bool(getattr(model, 'converged_', False)),
        'final_loss': final_loss,
    }
    with open(os.path.join(factors_dir, 'metadata.json'), 'w') as fh:
        json.dump(metadata, fh, indent=2)

    # Figures
    plt.figure(figsize=(6, 4))
    plt.bar(np.arange(len(weights)), weights)
    plt.xlabel('Component')
    plt.ylabel('Weight')
    plt.title('Component Weights')
    plt.tight_layout()
    plt.savefig(os.path.join(figs_dir, 'component_weights.png'), dpi=200)
    plt.close()

    # Time loadings across components
    plt.figure(figsize=(7, 4))
    for k in range(time_factors.shape[1]):
        plt.plot(range(1, time_factors.shape[0] + 1), 
                time_factors.iloc[:, k], marker='o', label=f'C{k+1}')
    plt.xticks(range(1, time_factors.shape[0] + 1))
    plt.xlabel('Timepoint')
    plt.ylabel('Loading')
    plt.title('Time Loadings by Component')
    plt.legend(ncols=2, fontsize=8)
    plt.tight_layout()
    plt.savefig(os.path.join(figs_dir, 'time_loadings.png'), dpi=200)
    plt.close()

    # Sample loadings colored by species
    fig, ax = plt.subplots(figsize=(10, 6))
    species_colors = {'apul': 'C0', 'peve': 'C1', 'ptua': 'C2'}
    
    for k in range(sample_factors.shape[1]):
        for idx, label in enumerate(sample_labels):
            species = sample_metadata[idx]['species']
            ax.scatter(k, sample_factors.iloc[idx, k], 
                      c=species_colors[species], alpha=0.6, s=50)
    
    ax.set_xlabel('Component')
    ax.set_ylabel('Sample Loading')
    ax.set_title('Sample Loadings by Component and Species')
    
    # Legend
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor=species_colors[sp], label=sp) 
                      for sp in ['apul', 'peve', 'ptua']]
    ax.legend(handles=legend_elements)
    plt.tight_layout()
    plt.savefig(os.path.join(figs_dir, 'sample_loadings.png'), dpi=200)
    plt.close()


def main() -> None:
    parser = argparse.ArgumentParser(
        description='Build tensor (genes × timepoints × samples) and run Barnacle SparseCP')
    parser.add_argument('--input-file', required=True, 
                       help='Path to merged vst_counts_matrix.csv file')
    parser.add_argument('--output-dir', required=True, 
                       help='Output directory for results')
    parser.add_argument('--rank', type=int, default=5)
    parser.add_argument('--lambda-gene', type=float, default=0.1)
    parser.add_argument('--lambda-time', type=float, default=0.05)
    parser.add_argument('--lambda-sample', type=float, default=0.1)
    parser.add_argument('--max-iter', type=int, default=1000)
    parser.add_argument('--tol', type=float, default=1e-5)
    parser.add_argument('--seed', type=int, default=42)
    args = parser.parse_args()

    if not os.path.exists(args.input_file):
        raise FileNotFoundError(f"Input file not found: {args.input_file}")

    df = read_normalized_csv(args.input_file)
    tensor, sample_labels, sample_metadata, genes = build_tensor_by_sample(
        df, expected_timepoints=[1, 2, 3, 4])

    print(f"\nTensor shape: {tensor.shape}")
    print(f"  Genes: {tensor.shape[0]}")
    print(f"  Timepoints: {tensor.shape[1]}")
    print(f"  Samples: {tensor.shape[2]}")

    model, decomposition = run_sparse_cp(
        tensor=tensor,
        rank=args.rank,
        lambda_gene=args.lambda_gene,
        lambda_time=args.lambda_time,
        lambda_sample=args.lambda_sample,
        max_iter=args.max_iter,
        tol=args.tol,
        seed=args.seed,
    )

    save_outputs(args.output_dir, tensor, sample_labels, sample_metadata, 
                genes, model, decomposition)


if __name__ == '__main__':
    main()