# 14-barnacle

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
uv run python M-multi-species/scripts/14-barnacle/build_tensor_and_run.py \
  --input-dir M-multi-species/output/13.00-multiomics-barnacle \
  --output-dir M-multi-species/output/14-barnacle \
  --rank 5 --lambda-gene 0.1 --lambda-sample 0.1 --lambda-time 0.05 \
  --max-iter 1000 --tol 1e-5 --seed 42
```

## Outputs

Written to `M-multi-species/output/14-barnacle/`:

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

## Synthetic Data

For testing and parameter tuning, this directory includes tools to generate synthetic data with high convergence probability. See [SYNTHETIC_DATA.md](SYNTHETIC_DATA.md) for details.

### Quick Start with Synthetic Data

```bash
# 1. Generate synthetic data
python M-multi-species/scripts/14-barnacle/generate_synthetic_data.py \
  --output-dir M-multi-species/output/14-barnacle-synthetic \
  --n-genes 10223 \
  --n-samples-per-species 10 \
  --n-components 5 \
  --noise-level 0.1 \
  --seed 42

# 2. Test data compatibility
python M-multi-species/scripts/14-barnacle/test_synthetic_data.py

# 3. Run tensor decomposition (requires barnacle installation)
uv run python M-multi-species/scripts/14-barnacle/build_tensor_and_run.py \
  --input-dir M-multi-species/output/14-barnacle-synthetic \
  --output-dir M-multi-species/output/14-barnacle-synthetic-results \
  --rank 5 --lambda-gene 0.1 --lambda-sample 0.1 --lambda-time 0.05 \
  --max-iter 2000 --tol 1e-4 --seed 42
```

The synthetic data is designed to converge well because it has:
- Clear temporal patterns
- Controlled sparsity
- High signal-to-noise ratio
- Complete data (no missing values)

See [SYNTHETIC_DATA.md](SYNTHETIC_DATA.md) for more details on the synthetic data generation and validation.

