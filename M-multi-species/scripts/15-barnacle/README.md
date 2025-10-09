# 15-barnacle

This directory replicates the 14-barnacle pipeline and runs Barnacle sparse CP decomposition at multiple ranks (3, 8, 10, 12), writing each run to its own subdirectory under `M-multi-species/output/15-barnacle/`.

## Inputs

- Normalized CSVs in `M-multi-species/output/13.00-multiomics-barnacle/`:
  - `apul_normalized_expression.csv`
  - `peve_normalized_expression.csv`
  - `ptua_normalized_expression.csv`

Each CSV must have a `group_id` column and columns named like `<SAMPLE>.<TP#>` (e.g., `POC.201.TP3`).

## Usage

Run all target ranks and generate a comparison summary:

```bash
./M-multi-species/scripts/15-barnacle/run_all.sh
```

Or run a single rank manually:

```bash
uv run python M-multi-species/scripts/15-barnacle/build_tensor_and_run.py \
  --input-dir M-multi-species/output/13.00-multiomics-barnacle \
  --output-dir M-multi-species/output/15-barnacle/rank-8 \
  --rank 8 --lambda-gene 0.1 --lambda-sample 0.1 --lambda-time 0.05 \
  --max-iter 1000 --tol 1e-5 --seed 42
```

## Outputs

Written to `M-multi-species/output/15-barnacle/`:

- `rank-3/`, `rank-8/`, `rank-10/`, `rank-12/` (each contains):
  - `multiomics_tensor.npy`
  - `barnacle_factors/`
    - `gene_factors.csv`
    - `sample_factors.csv`
    - `time_factors.csv`
    - `component_weights.csv`
    - `sample_mapping.csv`
    - `metadata.json` (includes `rank`, `final_loss`, and convergence flag)
  - `figures/`
    - `component_weights.png`
    - `time_loadings.png`
- `SUMMARY.md` (comparison across ranks)

Timepoints are assumed to be {1,2,3,4}. Missing values are filled with 0.0 before decomposition.

