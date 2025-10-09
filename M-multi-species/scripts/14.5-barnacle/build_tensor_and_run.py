#!/usr/bin/env python3
import argparse
import json
import os
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from barnacle.decomposition import SparseCP


# --------------------------- I/O helpers ---------------------------

def discover_species_files(input_dir: Path) -> Dict[str, Path]:
    files = list(input_dir.glob("*_normalized_expression.csv"))
    if not files:
        raise FileNotFoundError(f"No files like *_normalized_expression.csv in {input_dir}")
    out = {}
    for p in files:
        # species code = prefix before first underscore
        sp = p.name.split("_", 1)[0]
        out[sp] = p
    return out


def read_normalized_csv(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    if "group_id" not in df.columns:
        raise ValueError(f"Expected 'group_id' column in {path}")
    df = df.set_index("group_id")
    return df


def parse_sample_timepoint(column_name: str) -> Tuple[str, int]:
    parts = column_name.split(".")
    if len(parts) < 2:
        raise ValueError(f"Unexpected column format (need SAMPLE.TP#): {column_name}")
    time_token = parts[-1]
    if not time_token.startswith("TP"):
        raise ValueError(f"Expected TP# token at end of column: {column_name}")
    try:
        tp = int(time_token.replace("TP", ""))
    except Exception as exc:
        raise ValueError(f"Failed to parse timepoint from {column_name}") from exc
    sample_id = ".".join(parts[:-1])
    return sample_id, tp


# --------------------------- preprocessing ---------------------------

def scale_matrix(df: pd.DataFrame, mode: str) -> pd.DataFrame:
    """
    mode:
      - 'none': return as is
      - 'per-gene-z': z-score each row (gene) across columns (samples/timepoints)
    """
    if mode == "none":
        return df
    if mode == "per-gene-z":
        # population std (ddof=0); if std == 0, fill with 0 afterward
        means = df.mean(axis=1)
        stds = df.std(axis=1, ddof=0).replace(0, np.nan)
        z = (df.sub(means, axis=0)).div(stds, axis=0)
        return z.fillna(0.0)
    raise ValueError(f"Unknown scale mode: {mode}")


def build_tensor(
    species_to_df: Dict[str, pd.DataFrame],
    expected_timepoints: List[int],
    min_tps_per_sample: int = 1,
) -> Tuple[np.ndarray, List[str], Dict[int, Dict[str, str]], List[str], np.ndarray]:
    """
    Returns:
      tensor: (genes, samples, time)
      sample_labels: ["species_sample"]
      species_sample_map: {sample_index: {"species": sp, "sample_id": sid}}
      genes: ordered list of common group_ids
      missing_mask: same shape as tensor; 1 = observed, 0 = missing
    """
    # Intersect genes across all species
    common_genes = None
    for df in species_to_df.values():
        common_genes = df.index if common_genes is None else common_genes.intersection(df.index)
    common_genes = sorted(list(common_genes))
    if len(common_genes) == 0:
        raise ValueError("No intersecting group_id genes across species")

    # Parse columns into (sample_id, tp) and collect per species
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
        species_info[species] = {"sample_ids": sample_ids, "sample_map": sample_map}

    # Build sample axis with only samples meeting coverage threshold
    sample_labels: List[str] = []
    species_sample_map: Dict[int, Dict[str, str]] = {}
    combined_idx = 0

    for species in sorted(species_to_df.keys()):  # stable ordering
        info = species_info[species]
        for sample_id in sorted(info["sample_ids"]):
            # Count how many expected timepoints present
            present = sum((sample_id, tp) in info["sample_map"] for tp in expected_timepoints)
            if present >= min_tps_per_sample:
                sample_labels.append(f"{species}_{sample_id}")
                species_sample_map[combined_idx] = {"species": species, "sample_id": sample_id}
                combined_idx += 1

    n_genes = len(common_genes)
    n_samples = len(sample_labels)
    n_time = len(expected_timepoints)
    if n_samples == 0:
        raise ValueError("No samples met the minimum timepoint coverage.")

    tensor = np.empty((n_genes, n_samples, n_time), dtype=float)
    tensor[:] = np.nan
    missing_mask = np.zeros_like(tensor, dtype=float)

    for s_idx, label in enumerate(sample_labels):
        species, sample_id = label.split("_", 1)
        df = species_to_df[species]
        df_sub = df.loc[common_genes]
        for t_idx, tp in enumerate(expected_timepoints):
            key = (sample_id, tp)
            col = species_info[species]["sample_map"].get(key)
            if col is None:
                continue
            values = df_sub[col].to_numpy(dtype=float)
            tensor[:, s_idx, t_idx] = values
            missing_mask[:, s_idx, t_idx] = 1.0

    return tensor, sample_labels, species_sample_map, common_genes, missing_mask


def impute_missing(tensor: np.ndarray, missing_mask: np.ndarray, policy: str) -> np.ndarray:
    if policy == "zero":
        out = tensor.copy()
        out[np.isnan(out)] = 0.0
        return out

    out = tensor.copy()
    if policy == "global-mean":
        gm = np.nanmean(out)
        out[np.isnan(out)] = gm if not np.isnan(gm) else 0.0
        return out

    if policy == "gene-mean":
        # compute mean per gene across (sample,time)
        G, _, _ = out.shape
        for g in range(G):
            m = np.nanmean(out[g, :, :])
            if np.isnan(m):
                m = 0.0
            mask = np.isnan(out[g, :, :])
            out[g, :, :][mask] = m
        return out

    raise ValueError(f"Unknown missing-policy: {policy}")


# --------------------------- modeling ---------------------------

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
    np.random.seed(seed)
    model = SparseCP(
        rank=rank,
        lambdas=[lambda_gene, lambda_sample, lambda_time],
        nonneg_modes=[0],              # nonnegativity on gene mode
        n_initializations=3,
        random_state=seed,
        n_iter_max=max_iter,
        tol=tol,
    )
    # If your Barnacle supports masks/weights for missing data, pass them here instead of imputing.
    decomposition = model.fit_transform(tensor, verbose=1)
    return model, decomposition


# --------------------------- outputs ---------------------------

def save_outputs(
    output_dir: Path,
    tensor: np.ndarray,
    sample_labels: List[str],
    species_sample_map: Dict[int, Dict[str, str]],
    genes: List[str],
    model,
    decomposition,
    rank: int,
    args_namespace: argparse.Namespace,
    qc_tables: Dict[str, pd.DataFrame],
    expected_timepoints: List[int],
):
    output_dir.mkdir(parents=True, exist_ok=True)
    factors_dir = output_dir / "barnacle_factors"
    figs_dir = output_dir / "figures"
    factors_dir.mkdir(exist_ok=True, parents=True)
    figs_dir.mkdir(exist_ok=True, parents=True)

    # Save tensor
    np.save(output_dir / "multiomics_tensor.npy", tensor)

    # Factors
    gene_factors = pd.DataFrame(decomposition.factors[0], index=genes)
    sample_factors = pd.DataFrame(decomposition.factors[1], index=sample_labels)
    time_index = [f"TP{t}" for t in expected_timepoints]
    time_factors = pd.DataFrame(decomposition.factors[2], index=time_index)

    gene_factors.to_csv(factors_dir / "gene_factors.csv")
    sample_factors.to_csv(factors_dir / "sample_factors.csv")
    time_factors.to_csv(factors_dir / "time_factors.csv")

    # Component weights
    if hasattr(decomposition, "weights") and decomposition.weights is not None:
        weights = np.asarray(decomposition.weights)
    else:
        weights = np.ones(gene_factors.shape[1], dtype=float)
    pd.DataFrame({"weight": weights}).to_csv(factors_dir / "component_weights.csv", index=False)

    # Sample mapping
    mapping_rows = []
    for idx, label in enumerate(sample_labels):
        mapping_rows.append({
            "combined_index": idx,
            "label": label,
            "species": species_sample_map[idx]["species"],
            "sample_id": species_sample_map[idx]["sample_id"],
        })
    pd.DataFrame(mapping_rows).to_csv(factors_dir / "sample_mapping.csv", index=False)

    # Metadata + args
    raw_loss = getattr(model, "loss_", None)
    if raw_loss is None:
        final_loss = None
    else:
        try:
            final_loss = float(raw_loss[-1]) if isinstance(raw_loss, (list, tuple, np.ndarray)) and len(raw_loss) > 0 else float(raw_loss)
        except Exception:
            final_loss = None

    metadata = {
        "timestamp": datetime.utcnow().isoformat() + "Z",
        "tensor_shape": list(map(int, tensor.shape)),
        "n_components": int(gene_factors.shape[1]),
        "rank": int(rank),
        "model_converged": bool(getattr(model, "converged_", False)),
        "final_loss": final_loss,
        "expected_timepoints": expected_timepoints,  # add for clarity/provenance
    }
    with open(factors_dir / "metadata.json", "w") as fh:
        json.dump(metadata, fh, indent=2)

    # Save run args for provenance
    with open(output_dir / "args.json", "w") as fh:
        json.dump(vars(args_namespace), fh, indent=2)

    # QC tables
    for name, df in qc_tables.items():
        df.to_csv(output_dir / f"qc_{name}.csv", index=True)

    # Figures
    plt.figure(figsize=(6, 4))
    plt.bar(np.arange(len(weights)), weights)
    plt.xlabel("Component")
    plt.ylabel("Weight")
    plt.tight_layout()
    plt.savefig(figs_dir / "component_weights.png", dpi=200)
    plt.close()

    plt.figure(figsize=(7, 4))
    for k in range(time_factors.shape[1]):
        plt.plot(range(1, time_factors.shape[0] + 1), time_factors.iloc[:, k], marker="o", label=f"C{k+1}")
    plt.xticks(range(1, time_factors.shape[0] + 1), time_index)
    plt.xlabel("Timepoint")
    plt.ylabel("Loading")
    plt.legend(ncols=2, fontsize=8)
    plt.tight_layout()
    plt.savefig(figs_dir / "time_loadings.png", dpi=200)
    plt.close()


# --------------------------- QC helpers ---------------------------

def qc_summary(
    species_to_df: Dict[str, pd.DataFrame],
    expected_timepoints: List[int],
    sample_labels: List[str],
    species_sample_map: Dict[int, Dict[str, str]],
    tensor: np.ndarray,
    missing_mask: np.ndarray,
) -> Dict[str, pd.DataFrame]:
    # Missingness rates
    miss_rate = pd.DataFrame({
        "missing_fraction": [1.0 - np.mean(missing_mask)],
        "observed_fraction": [np.mean(missing_mask)],
    }, index=["all"])

    # Per species sample counts
    sp_counts = []
    for idx, lbl in enumerate(sample_labels):
        sp = species_sample_map[idx]["species"]
        sp_counts.append(sp)
    sp_counts = pd.Series(sp_counts).value_counts().rename_axis("species").to_frame("n_samples")

    # Coverage per species x timepoint
    records = []
    for sp, df in species_to_df.items():
        for tp in expected_timepoints:
            n = 0
            for col in df.columns:
                try:
                    _, tpp = parse_sample_timepoint(col)
                except Exception:
                    continue
                if tpp == tp:
                    n += 1
            records.append({"species": sp, "timepoint": tp, "n_columns": n})
    coverage = pd.DataFrame(records).sort_values(["species", "timepoint"])

    return {"missingness": miss_rate, "sample_counts": sp_counts, "coverage": coverage}


# --------------------------- main ---------------------------

def parse_timepoints(s: str) -> List[int]:
    try:
        return [int(x) for x in s.split(",") if x.strip() != ""]
    except Exception as e:
        raise argparse.ArgumentTypeError(f"Failed to parse --timepoints '{s}': {e}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Build tensor and run Barnacle SparseCP (improved)")
    parser.add_argument("--input-dir", required=True, help="Directory with *_normalized_expression.csv files")
    parser.add_argument("--output-dir", required=True, help="Output directory for results")

    # Modeling
    parser.add_argument("--rank", type=int, default=5)
    parser.add_argument("--lambda-gene", type=float, default=0.1)
    parser.add_argument("--lambda-sample", type=float, default=0.1)
    parser.add_argument("--lambda-time", type=float, default=0.05)
    parser.add_argument("--max-iter", type=int, default=1000)
    parser.add_argument("--tol", type=float, default=1e-5)
    parser.add_argument("--seed", type=int, default=42)

    # New features
    parser.add_argument("--timepoints", type=parse_timepoints, default=[1, 2, 3, 4],
                        help="Comma-separated list, e.g. '1,2,3,4'")
    parser.add_argument("--missing-policy", choices=["zero", "gene-mean", "global-mean"], default="gene-mean",
                        help="How to impute missing values before decomposition")
    parser.add_argument("--scale", choices=["none", "per-gene-z"], default="none",
                        help="Optional per-species scaling before stacking")
    parser.add_argument("--min-tps-per-sample", type=int, default=1,
                        help="Require a sample to have at least this many expected timepoints")
    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)

    # Discover & load species files
    species_files = discover_species_files(input_dir)
    species_to_df = {}
    for sp, p in species_files.items():
        df = read_normalized_csv(p)
        df = scale_matrix(df, args.scale)  # optional per-species scaling
        species_to_df[sp] = df

    # Build tensor
    tensor_raw, sample_labels, species_sample_map, genes, missing_mask = build_tensor(
        species_to_df, expected_timepoints=args.timepoints, min_tps_per_sample=args.min_tps_per_sample
    )

    # Impute missing if needed (Barnacle may not accept masks directly)
    tensor = impute_missing(tensor_raw, missing_mask, args.missing_policy)

    # Fit model
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

    # QC tables
    qc_tables = qc_summary(
        species_to_df=species_to_df,
        expected_timepoints=args.timepoints,
        sample_labels=sample_labels,
        species_sample_map=species_sample_map,
        tensor=tensor,
        missing_mask=missing_mask,
    )

    # Save everything
    save_outputs(
        output_dir=output_dir,
        tensor=tensor,
        sample_labels=sample_labels,
        species_sample_map=species_sample_map,
        genes=genes,
        model=model,
        decomposition=decomposition,
        rank=args.rank,
        args_namespace=args,
        qc_tables=qc_tables,
        expected_timepoints=args.timepoints,
    )


if __name__ == "__main__":
    main()