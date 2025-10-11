# Synthetic Data Generation for Barnacle - Summary

## Overview

This implementation provides a complete synthetic data generation pipeline for testing the Barnacle sparse CP decomposition on multi-species gene expression data.

## What Was Delivered

### 1. Synthetic Data Generator (`generate_synthetic_data.py`)

A Python script that creates realistic gene expression datasets with:

- **Three species**: apul, peve, ptua (matching real data)
- **Configurable parameters**:
  - Number of genes (default: 10,223 to match real data)
  - Number of samples per species (default: 10)
  - Number of underlying components (default: 5)
  - Noise level (default: 0.1)
  - Random seed for reproducibility

**Key Features**:
- Generates data from a known CP decomposition structure
- Includes distinct temporal patterns (increasing, decreasing, peaked)
- Sparse gene factors (70% sparsity by default)
- Non-negative values (realistic for expression data)
- Controlled noise for realism
- Outputs ground truth factors for validation

### 2. Test Script (`test_synthetic_data.py`)

Validates the synthetic data by:
- Loading all three species files
- Building the 3D tensor
- Verifying data properties
- Comparing with ground truth (0.989 correlation achieved)
- Confirming compatibility with `build_tensor_and_run.py`

### 3. Example Workflow (`example_workflow.sh`)

A complete end-to-end workflow script demonstrating:
1. Data generation
2. Validation
3. Tensor decomposition (requires barnacle installation)
4. Results inspection

### 4. Documentation

- **SYNTHETIC_DATA.md**: Comprehensive guide to synthetic data generation
- **Updated README.md**: Integrated synthetic data information into main README
- Clear usage examples and parameter explanations

### 5. Generated Datasets

Complete synthetic dataset in `M-multi-species/output/14-barnacle-synthetic/`:
- `apul_normalized_expression.csv` (10,223 genes × 40 columns)
- `peve_normalized_expression.csv` (10,223 genes × 40 columns)
- `ptua_normalized_expression.csv` (10,223 genes × 40 columns)
- `ground_truth/` directory with true factor matrices

## Why This Data Has High Convergence Probability

1. **Clear Structure**: Components have distinct, non-overlapping temporal patterns
2. **High Sparsity**: 70% of gene-component entries are zero (realistic and convergence-friendly)
3. **Non-negative Constraints**: All values ≥ 0 (natural for gene expression)
4. **High Signal-to-Noise Ratio**: Default noise level of 10% allows clear pattern recovery
5. **Complete Data**: No missing values (simplifies optimization)
6. **Realistic Scale**: Values in similar range to real normalized expression data (0-160)
7. **Known Ground Truth**: Can validate recovered factors against true factors

## Validation Results

The test script confirms:
- ✅ Data format matches expected structure
- ✅ All three species have 10,223 common genes
- ✅ Tensor shape: (10,223 genes, 30 samples, 4 timepoints)
- ✅ Correlation with ground truth: 0.989
- ✅ No missing data
- ✅ Compatible with build_tensor_and_run.py

## Recommended Parameters for Convergence

Based on the synthetic data structure:

```bash
--rank 5              # Matches ground truth
--lambda-gene 0.1     # Moderate gene regularization
--lambda-sample 0.1   # Moderate sample regularization  
--lambda-time 0.05    # Light time regularization
--max-iter 2000       # Increased from 1000 for better convergence
--tol 1e-4           # Slightly relaxed from 1e-5
--seed 42            # Reproducibility
```

## How to Use

### Quick Start

```bash
# Generate data
python M-multi-species/scripts/14-barnacle/generate_synthetic_data.py \
  --output-dir M-multi-species/output/14-barnacle-synthetic

# Validate
python M-multi-species/scripts/14-barnacle/test_synthetic_data.py

# Run decomposition (requires barnacle)
uv run python M-multi-species/scripts/14-barnacle/build_tensor_and_run.py \
  --input-dir M-multi-species/output/14-barnacle-synthetic \
  --output-dir M-multi-species/output/14-barnacle-synthetic-results \
  --rank 5 --lambda-gene 0.1 --lambda-sample 0.1 --lambda-time 0.05 \
  --max-iter 2000 --tol 1e-4 --seed 42
```

### Or Use the Workflow Script

```bash
bash M-multi-species/scripts/14-barnacle/example_workflow.sh
```

## Files Created

### Scripts
- `M-multi-species/scripts/14-barnacle/generate_synthetic_data.py` - Data generator
- `M-multi-species/scripts/14-barnacle/test_synthetic_data.py` - Validation script
- `M-multi-species/scripts/14-barnacle/example_workflow.sh` - Complete workflow

### Documentation
- `M-multi-species/scripts/14-barnacle/SYNTHETIC_DATA.md` - Detailed guide
- `M-multi-species/scripts/14-barnacle/README.md` - Updated with synthetic data info
- `M-multi-species/scripts/14-barnacle/SUMMARY.md` - This file

### Data
- `M-multi-species/output/14-barnacle-synthetic/apul_normalized_expression.csv`
- `M-multi-species/output/14-barnacle-synthetic/peve_normalized_expression.csv`
- `M-multi-species/output/14-barnacle-synthetic/ptua_normalized_expression.csv`
- `M-multi-species/output/14-barnacle-synthetic/ground_truth/true_gene_factors.csv`
- `M-multi-species/output/14-barnacle-synthetic/ground_truth/true_sample_factors.csv`
- `M-multi-species/output/14-barnacle-synthetic/ground_truth/true_time_factors.csv`

## Next Steps

1. **Install Barnacle**: `uv pip install git+https://github.com/blasks/barnacle.git@612b6a4`
2. **Run Decomposition**: Use the synthetic data to test convergence
3. **Compare Results**: Match recovered factors against ground truth
4. **Tune Parameters**: If needed, adjust regularization parameters
5. **Apply to Real Data**: Once parameters work well on synthetic data

## Benefits

This synthetic dataset provides:
- A reliable test case for the tensor decomposition pipeline
- Ground truth for validating the method
- A controlled environment for parameter tuning
- Confidence that the real data can converge with the right parameters
- A reproducible example for documentation and training
