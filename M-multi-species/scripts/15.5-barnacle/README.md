# 15-barnacle

This directory replicates the 14-barnacle pipeline and runs Barnacle sparse CP decomposition at multiple ranks (3, 8, 10, 12), writing each run to its own subdirectory under `M-multi-species/output/15-barnacle/`.

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

