# 14-barnacle

This directory contains a Python CLI to build a 3D tensor from normalized per‑species expression matrices and run Barnacle sparse CP decomposition.

## Inputs

- Normalized CSVs in `M-multi-species/output/13.00-multiomics-barnacle/`:
  - `apul_normalized_expression.csv`
  - `peve_normalized_expression.csv`
  - `ptua_normalized_expression.csv`

Each CSV must have a `group_id` column and columns named like `<SAMPLE>.<TP#>` (e.g., `POC.201.TP3`).

## Usage

```bash
uv run python M-multi-species/scripts/14.1-barnacle/build_tensor_and_run.py \
  --input-dir M-multi-species/output/13.00-multiomics-barnacle \
  --output-dir M-multi-species/output/14.1-barnacle \
  --rank 5 --lambda-gene 0.1 --lambda-sample 0.1 --lambda-time 0.05 \
  --max-iter 1000 --tol 1e-5 --seed 42
```

## Outputs

Written to `M-multi-species/output/14.1-barnacle/`:

- `multiomics_tensor.npy`
- `barnacle_factors/`
  - `gene_factors.csv`
  - `sample_factors.csv`
  - `time_factors.csv`
  - `component_weights.csv`
  - `sample_mapping.csv`
  - `metadata.json`
- `figures/`
  - `component_weights.png`
  - `time_loadings.png`

Timepoints are assumed to be {1,2,3,4}. Missing values are filled with 0.0 before decomposition.


