# 14.1-barnacle

This directory contains a Python CLI to build a 3D tensor from a merged normalized expression matrix and run Barnacle sparse CP decomposition.

## Inputs

- Merged normalized CSV in `M-multi-species/output/14-pca-orthologs/`:
  - `vst_counts_matrix.csv`

The CSV must have:
- A `group_id` column containing ortholog group identifiers
- Sample columns named like `<SPECIES_PREFIX>-<SAMPLE_ID>-<TP#>` where:
  - `SPECIES_PREFIX` is ACR (Acropora pulchra), POR (Pocillopora verrucosa), or POC (Pocillopora meandrina)
  - `SAMPLE_ID` is the unique sample identifier (e.g., 139, 145, 150)
  - `TP#` is the timepoint (TP1, TP2, TP3, TP4)
  - Example: `ACR-139-TP3`, `POR-216-TP1`, `POC-201-TP4`

This merged format combines all species into a single matrix, unlike the separate species files used by other Barnacle pipelines.

## Usage

```bash
uv run python M-multi-species/scripts/14.1-barnacle/build_tensor_and_run.py \
  --input-file M-multi-species/output/14-pca-orthologs/vst_counts_matrix.csv \
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

## Parameter Optimization and Convergence Testing

This directory includes a systematic convergence testing script (`convergence_test.py`) that helps identify optimal Barnacle parameters for your specific dataset. The script uses a two-step approach:

### Step 1: Rank Optimization
Tests different rank values (5, 10, 15, 20, 25, 30, 35) with baseline parameters to identify the rank that achieves the best convergence performance.

### Step 2: Parameter Grid Search
For the optimal rank identified in Step 1, the script tests a comprehensive grid of parameter combinations including:
- **max_iter**: [1000, 2000, 5000, 10000, 15000, 20000, 25000]
- **tol**: [1e-5, 1e-4, 1e-3]
- **lambda_gene**: Various combinations focusing on gene regularization (0.01 to 1.0)
- **lambda_sample**: Sample regularization values (0.1 to 0.5)
- **lambda_time**: Time regularization values (0.05 to 0.2)

### Usage

```bash
# Run the full convergence testing pipeline
uv run python M-multi-species/scripts/14.1-barnacle/convergence_test.py \
  --input-file M-multi-species/output/14-pca-orthologs/vst_counts_matrix.csv \
  --output-dir M-multi-species/output/14.1-barnacle-convergence-test \
  --results-file convergence_test_results.csv
```

### Results and Interpretation

The script generates:
- **Intermediate results** (saved every 10 tests)
- **Final results summary** in CSV format
- **Detailed results** in JSON format for programmatic analysis
- **Console output** showing convergence status for each parameter combination

**Success criteria:**
- **Convergence**: The algorithm reaches the specified tolerance (`tol`) within `max_iter` iterations
- **Success**: The script completes without errors
- **Optimal parameters**: Among converged runs, the combination with lowest final loss is recommended

**Parameter selection strategy:**
1. Prioritize convergence over all other criteria
2. Among converged solutions, select the one with lowest final loss
3. If no convergence is achieved, use the successful run with lowest loss
4. Fallback to rank 20 if no runs succeed

The convergence testing helps ensure your Barnacle decomposition achieves reliable, reproducible results by systematically exploring the parameter space and identifying settings that work best for your specific dataset structure.


