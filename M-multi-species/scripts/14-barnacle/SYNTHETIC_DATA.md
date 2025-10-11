# Synthetic Data Generation for Barnacle

This directory includes a script to generate synthetic gene expression data for testing the Barnacle sparse CP decomposition.

## Purpose

The synthetic data generator creates gene expression datasets with known underlying structure, which:
- Helps validate the tensor decomposition approach
- Provides data with high probability of convergence
- Allows testing different parameter configurations
- Enables comparison between recovered and ground-truth factors

## Data Characteristics

The synthetic data mimics the real barnacle multi-species dataset:

- **Three species**: `apul` (Acropora pulchra), `peve` (Porites evermanni), `ptua` (Pocillopora tuahiniensis)
- **Common genes**: All species share the same set of ortholog groups (OG_00001, OG_00002, etc.)
- **Samples**: Configurable number of samples per species (default: 10)
- **Timepoints**: 4 timepoints (TP1-TP4) per sample
- **Format**: CSV files with `group_id` column and columns named `SAMPLE.TP#`

## Synthetic Data Structure

The synthetic data is generated from a controlled CP decomposition:

1. **Gene factors**: Sparse matrix where each component involves a subset of genes
2. **Sample factors**: Mixed component weights for each sample (biological variation)
3. **Time factors**: Distinct temporal patterns (increasing, decreasing, peaked, etc.)
4. **Noise**: Gaussian noise added to simulate measurement error
5. **Non-negativity**: All values are non-negative (like real expression data)

This structure ensures:
- Clear component separation for easier convergence
- Realistic sparsity patterns
- Interpretable temporal dynamics

## Usage

### Generate Synthetic Data

```bash
python M-multi-species/scripts/14-barnacle/generate_synthetic_data.py \
  --output-dir M-multi-species/output/14-barnacle-synthetic \
  --n-genes 10223 \
  --n-samples-per-species 10 \
  --n-components 5 \
  --noise-level 0.1 \
  --seed 42
```

### Parameters

- `--output-dir`: Directory to save synthetic CSV files (required)
- `--n-genes`: Number of genes (default: 10223, matching real data)
- `--n-samples-per-species`: Samples per species (default: 10)
- `--n-components`: Number of underlying components (default: 5)
- `--noise-level`: Noise level as fraction of signal (default: 0.1)
- `--seed`: Random seed for reproducibility (default: 42)

### Run Tensor Decomposition on Synthetic Data

After generating the data, run the standard pipeline:

```bash
python M-multi-species/scripts/14-barnacle/build_tensor_and_run.py \
  --input-dir M-multi-species/output/14-barnacle-synthetic \
  --output-dir M-multi-species/output/14-barnacle-synthetic-results \
  --rank 5 \
  --lambda-gene 0.1 \
  --lambda-sample 0.1 \
  --lambda-time 0.05 \
  --max-iter 1000 \
  --tol 1e-5 \
  --seed 42
```

## Outputs

The script generates:

1. **Three species CSV files**:
   - `apul_normalized_expression.csv`
   - `peve_normalized_expression.csv`
   - `ptua_normalized_expression.csv`

2. **Ground truth factors** (in `ground_truth/` subdirectory):
   - `true_gene_factors.csv`: True gene loadings
   - `true_sample_factors.csv`: True sample loadings
   - `true_time_factors.csv`: True temporal patterns

You can compare the recovered factors from the decomposition with these ground truth factors to validate the method.

## Example: Full Workflow

```bash
# 1. Generate synthetic data
python M-multi-species/scripts/14-barnacle/generate_synthetic_data.py \
  --output-dir M-multi-species/output/14-barnacle-synthetic \
  --n-genes 1000 \
  --n-samples-per-species 8 \
  --n-components 5 \
  --noise-level 0.05 \
  --seed 123

# 2. Run tensor decomposition
python M-multi-species/scripts/14-barnacle/build_tensor_and_run.py \
  --input-dir M-multi-species/output/14-barnacle-synthetic \
  --output-dir M-multi-species/output/14-barnacle-synthetic-results \
  --rank 5 \
  --lambda-gene 0.1 \
  --lambda-sample 0.1 \
  --lambda-time 0.05 \
  --max-iter 2000 \
  --tol 1e-4 \
  --seed 42

# 3. Check convergence in results
cat M-multi-species/output/14-barnacle-synthetic-results/barnacle_factors/metadata.json
```

## Why This Data Converges Well

The synthetic data is designed for convergence because:

1. **Clear structure**: Components have distinct temporal patterns
2. **Sparsity**: Most genes are zero for most components
3. **Non-negativity**: Data and factors are non-negative (natural for gene expression)
4. **Sufficient signal**: Signal-to-noise ratio is high enough
5. **No missing values**: Complete data for all samples and timepoints
6. **Realistic scale**: Values are in a similar range to real normalized expression data

This makes it an ideal test case for parameter tuning and validation.
