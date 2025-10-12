# New Sparse Tensor Decomposition Workflow

A Python package for sparse CP tensor decomposition of multi-species time series molecular data using L1 regularization and Optuna hyperparameter optimization.

## Installation

### Using uv (recommended)

```bash
# Navigate to the package directory
cd M-multi-species/scripts/20-new-sparse-tensor-decompositon-model

# Install package in editable mode
uv pip install -e .
```

### Using pip (alternative)

```bash
# Navigate to the package directory
cd M-multi-species/scripts/20-new-sparse-tensor-decompositon-model

# Install dependencies
pip install numpy pandas scipy tensorly optuna typer matplotlib pyyaml

# The package can be run directly without installation using PYTHONPATH
```

## Quick Start

### Using uv (if installed)

```bash
# Navigate to the package directory
cd M-multi-species/scripts/20-new-sparse-tensor-decompositon-model

# 1. Build Tensor
uv run new-tensor build-tensor \
    --input-path /Users/sr320/Documents/GitHub/timeseries_molecular/M-multi-species/output/14-pca-orthologs/vst_counts_matrix.csv \
    --aggregation-method median \
    --normalization zscore

# 2. Optimize Hyperparameters
uv run new-tensor optimize \
    --tensor-path output/tensor/tensor.npz \
    --n-trials 50 \
    --rank-min 2 \
    --rank-max 8

# 3. Fit Model
uv run new-tensor fit \
    --tensor-path output/tensor/tensor.npz \
    --rank 5 \
    --non-negative

# 4. Export Results
uv run new-tensor export \
    --fit-dir output/fit \
    --mappings-dir output/tensor
```

### Using direct Python execution (works without uv)

```bash
# Navigate to the package directory
cd M-multi-species/scripts/20-new-sparse-tensor-decompositon-model

# 1. Build Tensor
PYTHONPATH=src python -m new_tensor.cli build-tensor \
    --input-path /Users/sr320/Documents/GitHub/timeseries_molecular/M-multi-species/output/14-pca-orthologs/vst_counts_matrix.csv \
    --aggregation-method median \
    --normalization zscore

# 2. Optimize Hyperparameters
PYTHONPATH=src python -m new_tensor.cli optimize \
    --tensor-path output/tensor/tensor.npz \
    --n-trials 50 \
    --rank-min 2 \
    --rank-max 8

# 3. Fit Model
PYTHONPATH=src python -m new_tensor.cli fit \
    --tensor-path output/tensor/tensor.npz \
    --rank 5 \
    --non-negative

# 4. Export Results
PYTHONPATH=src python -m new_tensor.cli export \
    --fit-dir output/fit \
    --mappings-dir output/tensor
```

## Command Reference

### `build-tensor`

Build 3D tensor from CSV data with replicate aggregation and normalization.

**Options:**
- `--input-path`: Path to input CSV file (required)
- `--output-dir`: Output directory (default: output/tensor)
- `--aggregation-method`: Replicate aggregation method (default: median)
- `--normalization`: Normalization method (default: zscore)
- `--min-expression`: Minimum expression threshold (default: 1.0)
- `--min-variance-percentile`: Variance percentile threshold (default: 10.0)
- `--config`: Path to config YAML file

### `optimize`

Run hyperparameter optimization using Optuna.

**Options:**
- `--tensor-path`: Path to input tensor .npz file (required)
- `--output-dir`: Output directory (default: output/optimization)
- `--n-trials`: Number of optimization trials (default: 100)
- `--rank-min`: Minimum rank to try (default: 2)
- `--rank-max`: Maximum rank to try (default: 12)
- `--lambda-min`: Minimum L1 penalty (default: 1e-4)
- `--lambda-max`: Maximum L1 penalty (default: 1.0)
- `--config`: Path to config YAML file

### `fit`

Fit CP decomposition with specified parameters.

**Options:**
- `--tensor-path`: Path to input tensor .npz file (required)
- `--output-dir`: Output directory (default: output/fit)
- `--rank`: Rank for decomposition (required)
- `--lambda-a`: L1 penalty for gene factors (default: 0.1)
- `--lambda-b`: L1 penalty for species factors (default: 0.1)
- `--lambda-c`: L1 penalty for time factors (default: 0.1)
- `--non-negative`: Enforce non-negativity (default: false)
- `--max-iter`: Maximum iterations (default: 100)
- `--config`: Path to config YAML file

### `export`

Export results in various formats including plots.

**Options:**
- `--fit-dir`: Directory containing fit results (required)
- `--mappings-dir`: Directory containing tensor mappings (required)
- `--output-dir`: Output directory (default: output/export)
- `--plot-heatmaps`: Generate factor heatmaps (default: true)

## Configuration

Use `config.yaml` to set default parameters:

```yaml
input:
  file_path: "M-multi-species/output/14-pca-orthologs/vst_counts_matrix.csv"

preprocessing:
  aggregation_method: "median"
  normalization: "zscore"

optimization:
  n_trials: 100
  rank_min: 2
  rank_max: 12

fitting:
  rank: 6
  non_negative: false
```

## Output Structure

```
output/
├── tensor/
│   ├── tensor.npz
│   ├── genes.csv
│   ├── species.csv
│   ├── timepoints.csv
│   └── tensor_shapes.json
├── optimization/
│   ├── optuna_study.pkl
│   ├── best_params.json
│   ├── trials.csv
│   └── top_trials.csv
├── fit/
│   ├── gene_factors.csv
│   ├── species_factors.csv
│   ├── time_factors.csv
│   └── fit_metrics.json
└── export/
    ├── gene_factors_with_ids.csv
    ├── species_factors_with_codes.csv
    ├── time_factors_with_labels.csv
    ├── decomposition_results.json
    └── factor_heatmaps.png
```

## Algorithm Details

### Sparse CP Decomposition

The package implements sparse CP decomposition using Alternating Least Squares (ALS) with L1 regularization:

```
min 0.5‖X - [A,B,C]‖² + λ_A‖A‖₁ + λ_B‖B‖₁ + λ_C‖C‖₁
```

Where:
- `X ∈ ℝ^{G×S×T}`: Gene expression tensor
- `[A,B,C]`: CP decomposition with factors A (genes), B (species), C (time)
- `λ_A, λ_B, λ_C`: L1 penalty parameters for sparsity

### Hyperparameter Optimization

Uses Optuna with cross-validation via random masking to optimize:
- Rank R (number of components)
- L1 penalties λ_A, λ_B, λ_C
- Non-negativity constraints

## Alternative Execution Methods

### If uv is not available

If you don't have `uv` installed, you can run the workflow using direct Python execution:

```bash
# Navigate to the package directory
cd M-multi-species/scripts/20-new-sparse-tensor-decompositon-model

# Run commands with PYTHONPATH set
PYTHONPATH=src python -m new_tensor.cli [COMMAND] [OPTIONS]
```

### Troubleshooting

- **Module not found errors**: Make sure you're in the correct directory and using `PYTHONPATH=src`
- **Command not found**: If `uv` is not installed, use the direct Python execution method above
- **Permission errors**: Ensure you have write permissions in the output directories

## Dependencies

- **Core**: numpy, pandas, scipy
- **Tensor operations**: tensorly
- **Optimization**: optuna
- **CLI**: typer
- **Plotting**: matplotlib, seaborn
- **Configuration**: pyyaml

## License

MIT License
