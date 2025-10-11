#!/usr/bin/env python3
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


def parse_sample_timepoint(column_name: str) -> Tuple[str, int]:
    # Expect <SAMPLE>.<TP#>, e.g., POC.201.TP3
    parts = column_name.split('.')
    if len(parts) < 3:
        raise ValueError(f"Unexpected column format (need SAMPLE.TP#): {column_name}")
    # sample id is everything before the last token (TP#)
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

    # Build combined sample list preserving species blocks and sample order
    sample_labels: List[str] = []
    species_sample_map: Dict[int, Dict[str, str]] = {}
    all_columns: List[str] = []
    combined_idx = 0
    for species in ['apul', 'peve', 'ptua']:
        df = species_to_df[species]
        info = species_info[species]
        for sample_id in info['sample_ids']:
            cols = []
            has_any = False
            for tp in expected_timepoints:
                key = (sample_id, tp)
                if key in info['sample_map']:
                    cols.append(info['sample_map'][key])
                    has_any = True
            if has_any:
                all_columns.extend(cols)
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


def run_sparse_cp(
    tensor: np.ndarray,
    rank: int,
    lambda_gene: float,
    lambda_sample: float,
    lambda_time: float,
    max_iter: int,
    tol: float,
    seed: int,
):
    model = SparseCP(
        rank=rank,
        lambdas=[lambda_gene, lambda_sample, lambda_time],
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
    species_sample_map: Dict[int, Dict[str, str]],
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

    # Extract factors: order assumed [genes, samples, time]
    gene_factors = pd.DataFrame(decomposition.factors[0], index=genes)
    sample_factors = pd.DataFrame(decomposition.factors[1], index=sample_labels)
    time_factors = pd.DataFrame(decomposition.factors[2], index=[f'TP{t}' for t in range(1, tensor.shape[2] + 1)])

    gene_factors.to_csv(os.path.join(factors_dir, 'gene_factors.csv'))
    sample_factors.to_csv(os.path.join(factors_dir, 'sample_factors.csv'))
    time_factors.to_csv(os.path.join(factors_dir, 'time_factors.csv'))

    # Component weights (lambdas)
    # Check if decomposition provides meaningful weights (not all ones or None)
    weights_attr = getattr(decomposition, 'weights', None)
    if weights_attr is not None:
        weights = np.asarray(weights_attr).astype(float).ravel()
        # If weights are all essentially the same (like all 1.0), compute from factor norms
        if np.allclose(weights, weights[0]) if len(weights) > 0 else True:
            weights = None
    else:
        weights = None

    # Compute weights from factor matrix norms if needed
    if weights is None:
        # Get norms for each component across all three modes
        gene_norms = np.linalg.norm(gene_factors.values, axis=0)  # Frobenius norm per component
        sample_norms = np.linalg.norm(sample_factors.values, axis=0)
        time_norms = np.linalg.norm(time_factors.values, axis=0)
        # Component weight is product of norms (represents overall magnitude of each component)
        weights = gene_norms * sample_norms * time_norms

    pd.DataFrame({'weight': weights}).to_csv(os.path.join(factors_dir, 'component_weights.csv'), index=False)

    # Sample mapping
    mapping_rows = []
    for idx, label in enumerate(sample_labels):
        mapping_rows.append({
            'combined_index': idx,
            'label': label,
            'species': species_sample_map[idx]['species'],
            'sample_id': species_sample_map[idx]['sample_id'],
        })
    pd.DataFrame(mapping_rows).to_csv(os.path.join(factors_dir, 'sample_mapping.csv'), index=False)

    # Metadata
    # Handle loss, which may be a list/array over iterations
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
        'n_components': int(gene_factors.shape[1]),
        'model_converged': bool(getattr(model, 'converged_', False)),
        'final_loss': final_loss,
    }
    with open(os.path.join(factors_dir, 'metadata.json'), 'w') as fh:
        json.dump(metadata, fh, indent=2)

    # Figures: component weights and time loadings
    plt.figure(figsize=(6, 4))
    plt.bar(np.arange(len(weights)), weights)
    plt.xlabel('Component')
    plt.ylabel('Weight')
    plt.tight_layout()
    plt.savefig(os.path.join(figs_dir, 'component_weights.png'), dpi=200)
    plt.close()

    plt.figure(figsize=(7, 4))
    for k in range(time_factors.shape[1]):
        plt.plot(range(1, time_factors.shape[0] + 1), time_factors.iloc[:, k], marker='o', label=f'C{k+1}')
    plt.xticks(range(1, time_factors.shape[0] + 1))
    plt.xlabel('Timepoint')
    plt.ylabel('Loading')
    plt.legend(ncols=2, fontsize=8)
    plt.tight_layout()
    plt.savefig(os.path.join(figs_dir, 'time_loadings.png'), dpi=200)
    plt.close()


def main() -> None:
    parser = argparse.ArgumentParser(description='Build tensor and run Barnacle SparseCP')
    parser.add_argument('--input-dir', required=True, help='Directory with *normalized_expression.csv files')
    parser.add_argument('--output-dir', required=True, help='Output directory for results')
    parser.add_argument('--rank', type=int, default=5)
    parser.add_argument('--lambda-gene', type=float, default=0.1)
    parser.add_argument('--lambda-sample', type=float, default=0.1)
    parser.add_argument('--lambda-time', type=float, default=0.05)
    parser.add_argument('--max-iter', type=int, default=1000)
    parser.add_argument('--tol', type=float, default=1e-5)
    parser.add_argument('--seed', type=int, default=42)
    args = parser.parse_args()

    input_dir = args.input_dir
    output_dir = args.output_dir

    expected_files = {
        'apul': os.path.join(input_dir, 'apul_normalized_expression.csv'),
        'peve': os.path.join(input_dir, 'peve_normalized_expression.csv'),
        'ptua': os.path.join(input_dir, 'ptua_normalized_expression.csv'),
    }
    for sp, p in expected_files.items():
        if not os.path.exists(p):
            raise FileNotFoundError(f"Missing input for {sp}: {p}")

    species_to_df = {sp: read_normalized_csv(p) for sp, p in expected_files.items()}
    tensor, sample_labels, species_sample_map, genes = build_tensor(species_to_df, expected_timepoints=[1, 2, 3, 4])

    model, decomposition = run_sparse_cp(
        tensor=tensor,
        rank=args.rank,
        lambda_gene=args.lambda_gene,
        lambda_sample=args.lambda_sample,
        lambda_time=args.lambda_time,
        max_iter=args.max_iter,
        tol=args.tol,
        seed=args.seed,
    )

    save_outputs(output_dir, tensor, sample_labels, species_sample_map, genes, model, decomposition)


if __name__ == '__main__':
    main()


