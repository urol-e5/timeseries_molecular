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


# ---------------------------------------------------------
# LOADING + PARSING FUNCTIONS (unchanged)
# ---------------------------------------------------------

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
    tp = int(time_token.replace('TP', ''))

    sample_id = column_name[:tp_index]
    sample_prefix = sample_id.split('-')[0]

    species_map = {'ACR': 'apul', 'POR': 'peve', 'POC': 'ptua'}
    if sample_prefix not in species_map:
        raise ValueError(f"Unknown species prefix in column: {column_name}")

    species = species_map[sample_prefix]
    return species, sample_id, tp


def build_tensor(
    df: pd.DataFrame,
    expected_timepoints: List[int],
) -> Tuple[np.ndarray, List[str], Dict[int, Dict[str, str]], List[str]]:

    common_genes = sorted(list(df.index))

    species_info = {
        'apul': {'sample_ids': [], 'sample_map': {}},
        'peve': {'sample_ids': [], 'sample_map': {}},
        'ptua': {'sample_ids': [], 'sample_map': {}},
    }

    for col in df.columns:
        try:
            species, sample_id, tp = parse_sample_timepoint(col)
            if tp in expected_timepoints:
                key = (sample_id, tp)
                species_info[species]['sample_map'][key] = col
                if sample_id not in species_info[species]['sample_ids']:
                    species_info[species]['sample_ids'].append(sample_id)
        except ValueError:
            continue

    sample_labels: List[str] = []
    species_sample_map: Dict[int, Dict[str, str]] = {}
    all_columns: List[str] = []
    combined_idx = 0

    for species in ['apul', 'peve', 'ptua']:
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

    for s_idx, label in enumerate(sample_labels):
        species, sample_id = label.split('_', 1)
        for t_idx, tp in enumerate(expected_timepoints):
            key = (sample_id, tp)
            col = species_info[species]['sample_map'].get(key)
            if col is None:
                continue
            values = df.loc[common_genes, col].to_numpy(dtype=float)
            tensor[:, s_idx, t_idx] = values

    tensor = np.nan_to_num(tensor, nan=0.0)
    return tensor, sample_labels, species_sample_map, common_genes


# ---------------------------------------------------------
# MODEL RUNNER (unchanged)
# ---------------------------------------------------------

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
    decomposition = model.fit_transform(tensor, verbose=0)
    return model, decomposition


# ---------------------------------------------------------
# OUTPUT FUNCTIONS (unchanged)
# ---------------------------------------------------------

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

    gene_factors = pd.DataFrame(decomposition.factors[0], index=genes)
    sample_factors = pd.DataFrame(decomposition.factors[1], index=sample_labels)
    time_factors = pd.DataFrame(
        decomposition.factors[2],
        index=[f'TP{t}' for t in range(1, tensor.shape[2] + 1)]
    )

    gene_factors.to_csv(os.path.join(factors_dir, 'gene_factors.csv'))
    sample_factors.to_csv(os.path.join(factors_dir, 'sample_factors.csv'))
    time_factors.to_csv(os.path.join(factors_dir, 'time_factors.csv'))

    weights_attr = getattr(decomposition, 'weights', None)
    if weights_attr is not None:
        weights = np.asarray(weights_attr).astype(float).ravel()
        if np.allclose(weights, weights[0]):
            weights = None
    else:
        weights = None

    if weights is None:
        gene_norms = np.linalg.norm(gene_factors.values, axis=0)
        sample_norms = np.linalg.norm(sample_factors.values, axis=0)
        time_norms = np.linalg.norm(time_factors.values, axis=0)
        weights = gene_norms * sample_norms * time_norms

    pd.DataFrame({'weight': weights}).to_csv(
        os.path.join(factors_dir, 'component_weights.csv'), index=False
    )

    mapping_rows = []
    for idx, label in enumerate(sample_labels):
        mapping_rows.append({
            'combined_index': idx,
            'label': label,
            'species': species_sample_map[idx]['species'],
            'sample_id': species_sample_map[idx]['sample_id'],
        })
    pd.DataFrame(mapping_rows).to_csv(os.path.join(factors_dir, 'sample_mapping.csv'), index=False)

    raw_loss = getattr(model, 'loss_', None)
    if raw_loss is None:
        final_loss = None
    else:
        try:
            final_loss = float(raw_loss[-1]) if len(raw_loss) > 0 else float(raw_loss)
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


# ---------------------------------------------------------
# NEW SECTION: CROSS-VALIDATED RANK SEARCH
# ---------------------------------------------------------

def compute_sse(Y_true: np.ndarray, Y_pred: np.ndarray) -> float:
    return float(np.sum((Y_true - Y_pred)**2))


def reconstruct_to_shape(model, factors, shape):
    # Manually reconstruct tensor using CP decomposition formula
    # tensor[i,j,k] = sum_r gene[i,r] * sample[j,r] * time[k,r]
    n_components = factors[0].shape[1]
    reconstructed = np.zeros((factors[0].shape[0], factors[1].shape[0], factors[2].shape[0]))
    
    for r in range(n_components):
        component_tensor = np.einsum('i,j,k->ijk',
                                     factors[0][:, r],
                                     factors[1][:, r],
                                     factors[2][:, r])
        reconstructed += component_tensor
    
    return reconstructed[:, :shape[1], :shape[2]]


def split_tensor_by_species(tensor: np.ndarray, sample_labels: List[str]):
    species_blocks = {"apul": [], "peve": [], "ptua": []}

    for idx, label in enumerate(sample_labels):
        species = label.split("_")[0]
        species_blocks[species].append(idx)

    replicate_tensors = {}
    for sp, idx_list in species_blocks.items():
        idx_list = sorted(idx_list)
        if len(idx_list) > 0:
            replicate_tensors[sp] = tensor[:, idx_list, :]
    return replicate_tensors


def cross_validated_rank_search(
    tensor: np.ndarray,
    sample_labels: List[str],
    lambda_gene: float,
    lambda_sample: float,
    lambda_time: float,
    max_iter: int,
    tol: float,
    seed: int,
    R_grid: List[int],
):

    reps = split_tensor_by_species(tensor, sample_labels)
    rep_keys = list(reps.keys())

    results = {}

    for R in R_grid:
        sse_list = []

        for train_rep in rep_keys:
            Y_train = reps[train_rep]

            model = SparseCP(
                rank=R,
                lambdas=[lambda_gene, lambda_sample, lambda_time],
                nonneg_modes=[0],
                n_initializations=3,
                random_state=seed,
                n_iter_max=max_iter,
                tol=tol,
            )
            decomp = model.fit_transform(Y_train, verbose=0)

            for test_rep in rep_keys:
                if test_rep == train_rep:
                    continue

                Y_test = reps[test_rep]
                Y_pred = reconstruct_to_shape(model, decomp.factors, Y_test.shape)
                sse_list.append(compute_sse(Y_test, Y_pred))

        results[R] = float(np.mean(sse_list))
        print(f"[CV] R={R}, mean CV-SSE={results[R]:.4e}")

    best_R = min(results, key=results.get)
    print("\n=== BEST R SELECTED BY CROSS-VALIDATED SSE ===")
    print(f"R = {best_R}, mean CV-SSE = {results[best_R]:.4e}")

    return best_R, results


# ---------------------------------------------------------
# MAIN
# ---------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(description='Build tensor and run Barnacle SparseCP with CV rank selection')
    parser.add_argument('--input-file', required=True, help='Path to merged vst_counts_matrix.csv file')
    parser.add_argument('--output-dir', required=True, help='Output directory for results')
    parser.add_argument('--lambda-gene', type=float, default=0.1)
    parser.add_argument('--lambda-sample', type=float, default=0.1)
    parser.add_argument('--lambda-time', type=float, default=0.05)
    parser.add_argument('--max-iter', type=int, default=1000)
    parser.add_argument('--tol', type=float, default=1e-5)
    parser.add_argument('--seed', type=int, default=42)
    parser.add_argument('--r-grid', nargs='+', type=int, default=[5, 10, 15, 20, 25, 30])
    args = parser.parse_args()

    df = read_normalized_csv(args.input_file)
    tensor, sample_labels, species_sample_map, genes = build_tensor(df, expected_timepoints=[1, 2, 3, 4])

    # ---------------------------
    # Cross-validated rank search
    # ---------------------------
    best_R, cv_results = cross_validated_rank_search(
        tensor=tensor,
        sample_labels=sample_labels,
        lambda_gene=args.lambda_gene,
        lambda_sample=args.lambda_sample,
        lambda_time=args.lambda_time,
        max_iter=args.max_iter,
        tol=args.tol,
        seed=args.seed,
        R_grid=args.r_grid,
    )

    # ---------------------------
    # Final model using best R
    # ---------------------------
    model, decomposition = run_sparse_cp(
        tensor=tensor,
        rank=best_R,
        lambda_gene=args.lambda_gene,
        lambda_sample=args.lambda_sample,
        lambda_time=args.lambda_time,
        max_iter=args.max_iter,
        tol=args.tol,
        seed=args.seed,
    )

    save_outputs(args.output_dir, tensor, sample_labels, species_sample_map, genes, model, decomposition)


if __name__ == '__main__':
    main()
