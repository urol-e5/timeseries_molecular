# 14.5-barnacle

This directory contains a Python CLI to build a 3D tensor from normalized per‑species expression matrices and run Barnacle sparse CP decomposition.

## Inputs

- Normalized CSVs in `M-multi-species/output/13.00-multiomics-barnacle/`:
  - `apul_normalized_expression.csv` (Acropora pulchra)
  - `peve_normalized_expression.csv` (Pocillopora verrucosa)
  - `ptua_normalized_expression.csv` (Pocillopora meandrina)

Each CSV must have:
- A `group_id` column containing ortholog group identifiers
- Sample columns named like `<SPECIES_PREFIX>.<SAMPLE_ID>.<TP#>` where:
  - `SPECIES_PREFIX` is ACR (Acropora pulchra), POR (Pocillopora verrucosa), or POC (Pocillopora meandrina)
  - `SAMPLE_ID` is the unique sample identifier (e.g., 139, 145, 150)
  - `TP#` is the timepoint (TP1, TP2, TP3, TP4)
  - Example: `ACR.139.TP3`, `POR.216.TP1`, `POC.201.TP4`

## Usage

```bash
uv run python M-multi-species/scripts/14.5-barnacle/build_tensor_and_run.py \
  --input-dir M-multi-species/output/13.00-multiomics-barnacle \
  --output-dir M-multi-species/output/14.5-barnacle \
  --rank 5 --lambda-gene 0.1 --lambda-sample 0.1 --lambda-time 0.05 \
  --max-iter 1000 --tol 1e-5 --seed 42
```

or 

```
uv run python M-multi-species/scripts/14.5-barnacle/build_tensor_and_run.py \
  --input-dir data/normalized \
  --output-dir results/sparsecp \
  --timepoints 1,2,3,4 \
  --missing-policy gene-mean \
  --scale per-gene-z \
  --min-tps-per-sample 2 \
  --rank 6 \
  --lambda-gene 0.1 \
  --lambda-sample 0.05 \
  --lambda-time 0.05 \
  --max-iter 1500 \
  --tol 1e-5 \
  --seed 42
  ```

## Outputs

Written to `M-multi-species/output/14.5-barnacle/`:

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


