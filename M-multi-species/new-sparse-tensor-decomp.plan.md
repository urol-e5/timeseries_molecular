# New Sparse Tensor Decomposition Workflow (Python, uv)

## Scope

- Implement a Python package and CLI to load `vst_counts_matrix.csv`, aggregate replicates to a 3×4 species-time grid per gene, construct a tensor, run sparse CP decomposition with L1 penalties (optionally non-negative), and perform robust hyperparameter optimization. 
- Input: `M-multi-species/output/14-pca-orthologs/vst_counts_matrix.csv`
- Output: `M-multi-species/output/20-new-sparse-tensor-decompositon-model/`
- Code: `M-multi-species/scripts/20-new-sparse-tensor-decompositon-model/`

## Directory structure

- `M-multi-species/scripts/20-new-sparse-tensor-decompositon-model/`
- `pyproject.toml` (PEP 621; uv-compatible)
- `README.md`
- `src/new_tensor/`
 - `__init__.py`
 - `io.py` (load wide CSV, parse sample labels)
 - `preprocess.py` (filtering, normalization, replicate aggregation)
 - `tensor.py` (reshape to (genes, species, timepoints), mappings)
 - `model.py` (sparse CP with L1, nonneg option)
 - `optimize.py` (Optuna-based hyperparameter search)
 - `cli.py` (click/typer CLI: build, optimize, fit, export)
 - `plotting.py` (diagnostics, factor heatmaps)
 - `utils.py` (logging, seeding, I/O helpers)
- `config.yaml` (defaults: input path, species codes, TP mapping, agg method, search spaces)

## Data parsing assumptions

- Wide matrix with header `group_id` then columns `SPECIES-INDIVIDUAL-TP{1..4}` (e.g., `POC-40-TP3`).
- Species codes: `ACR`, `POR`, `POC`. Time points: `TP1..TP4`.
- Robust parsing via regex: `^(?P<species>[A-Z]{3})-(?P<ind>[^-]+)-TP(?P<tp>\d+)$`.
- Replicates aggregated per (species, timepoint) using configurable reducer: median (default), mean, or trimmed-mean.

## Workflow

1. Load CSV → validate header, extract species/timepoint from column names, map `group_id` to rows.
2. Aggregate replicates across individuals to 3×4 grid per gene; handle missing via robust reducer (ignore NaNs).
3. Optional normalization: per-gene z-score or min-max; optional variance filtering.
4. Build tensor `X ∈ ℝ^{G×S×T}` with `G≈10k`, `S=3`, `T=4`.
5. Fit sparse CP decomposition (rank R) with L1 on factors A (genes), B (species), C (time):

- Objective: 0.5‖X − Σ_r a_r∘b_r∘c_r‖²_F + λ_A‖A‖₁ + λ_B‖B‖₁ + λ_C‖C‖₁; optional A,B,C ≥ 0.
- Solver: ALS with proximal step (soft-threshold + nonneg projection), early stopping on relative loss.

6. Hyperparameter optimization (Optuna): search R, λ_A, λ_B, λ_C, nonneg flag, normalization, aggregator.

- CV via random masking of entries; metric: masked RMSE + sparsity regularization.
- Track trial logs, best params, and learning curves.

7. Final fit with best params on full tensor; export artifacts.

## CLI

- `build-tensor`: parse input, aggregate, normalize, save tensor `.npz` + mapping CSVs.
- `optimize`: run Optuna with config; save `study.db` or `study.json` and CV report.
- `fit`: train with selected/best R and lambdas; save factors and recon.
- `export`: write CSVs: `genes_factors.csv` (G×R), `species_factors.csv` (3×R), `time_factors.csv` (4×R); JSON: params, metrics; PNGs: heatmaps.

Example commands (from repo root):

- Install uv and deps:
- `uv venv && uv pip install -e M-multi-species/scripts/20-new-sparse-tensor-decompositon-model`
- Build tensor:
- `uv run new-tensor build-tensor --input M-multi-species/output/14-pca-orthologs/vst_counts_matrix.csv --agg median --norm zscore`
- Optimize:
- `uv run new-tensor optimize --n-trials 100 --rank-min 2 --rank-max 12`
- Fit and export:
- `uv run new-tensor fit --rank 6 --nonneg` then `uv run new-tensor export`

## Outputs (under `M-multi-species/output/20-new-sparse-tensor-decompositon-model/<run_id>/`)

- `tensor_shapes.json`, `mappings/genes.csv`, `mappings/species.csv`, `mappings/timepoints.csv`
- `study/optuna_study.json`, `cv_metrics.csv`, `best_params.json`
- `factors/genes_factors.csv`, `species_factors.csv`, `time_factors.csv`
- `reconstruction_summary.json`, `explained_variance.csv`
- `plots/factor_heatmaps.png`, `plots/recon_error.png`
- `logs/run.log`

## Performance/robustness

- Vectorized aggregation; small S×T dims keep ALS steps cheap.
- Early stopping + warm starts across trials.
- Reproducibility via fixed seeds; deterministic NumPy backend.
- Optionally enable sparse COO backend if future inputs include zeros/missing.

## Dependencies (pyproject)

- `numpy`, `pandas`, `scipy`, `tensorly>=0.9`, `optuna`, `typer[all]`, `matplotlib`, `pyyaml`
- Optional: `cupy` (GPU) behind flag; not required.

## Testing/demo

- Tiny smoke test using a 200-gene subset sampled from the provided CSV to verify end-to-end run in <1 min.

### To-dos

- [ ] Scaffold uv-compatible package and directory layout
- [ ] Implement CSV loader and sample-label parser for vst_counts_matrix.csv
- [ ] Implement replicate aggregation and normalization options
- [ ] Build (genes×species×time) tensor and save mappings
- [ ] Implement sparse CP with L1 and nonneg option
- [ ] Add Optuna search with CV masking and trial logging
- [ ] Create Typer CLI commands: build, optimize, fit, export
- [ ] Write README with uv install and usage examples
- [ ] Run smoke test on 200-gene subset and save outputs