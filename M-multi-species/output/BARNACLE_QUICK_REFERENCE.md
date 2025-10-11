# Barnacle Analysis Quick Reference

## Best Parameter Sets Identified

### Overall Best Performance (Lowest Loss)

| Configuration | Rank | λ_gene | λ_sample | λ_time | Max Iter | Tolerance | Final Loss | Converged | Location |
|--------------|------|--------|----------|--------|----------|-----------|------------|-----------|----------|
| **Best for Rank 5** | 5 | 0.01 | 0.1 | 0.05 | 1000+ | 1e-5 | ~840,915 | No | test_006, test_066 |
| **Best for Rank 12** | 12 | 0.1 | 0.1 | 0.05 | 1000 | 1e-5 | 587,541 | No | 15-barnacle/rank-12 |
| **Best for Rank 10** | 10 | 0.1 | 0.1 | 0.05 | 1000 | 1e-5 | 634,014 | No | 15-barnacle/rank-10 |
| **Best for Rank 8** | 8 | 0.1 | 0.1 | 0.05 | 1000 | 1e-5 | 689,401 | No | 15-barnacle/rank-8 |
| **Best for Rank 3** | 3 | 0.1 | 0.1 | 0.05 | 1000 | 1e-5 | 1,368,623 | No | 15-barnacle/rank-3 |

### Key Insights

1. **Lower λ_gene (0.01) achieves lowest loss** for Rank 5 across all convergence tests
2. **Higher ranks achieve lower loss** but require more components
3. **Iterations beyond 1000 provide minimal improvement**
4. **Relaxing tolerance does not achieve convergence**

## Quick Access to Results

### Convergence Testing Results (94 Tests)

```bash
# View all test results
cat M-multi-species/output/14.1-barnacle-convergence-tests/intermediate_results_80.csv

# View convergence summary
cat M-multi-species/output/14.1-barnacle-convergence-tests/BARNACLE_CONVERGENCE_SUMMARY.md

# Access specific test (e.g., test 6 - best performer)
ls -la M-multi-species/output/14.1-barnacle-convergence-tests/test_006/
```

### Rank Comparison Results

```bash
# View rank comparison summary
cat M-multi-species/output/15-barnacle/SUMMARY.md

# View specific rank results
ls -la M-multi-species/output/15-barnacle/rank-12/barnacle_factors/
cat M-multi-species/output/15-barnacle/rank-12/barnacle_factors/metadata.json

# View visualizations
open M-multi-species/output/15-barnacle/rank-12/figures/component_weights.png
open M-multi-species/output/15-barnacle/rank-12/figures/time_loadings.png
```

## Recommended Parameter Sets for Different Use Cases

### For Lowest Final Loss
```python
# Use Rank 12 with standard parameters
rank = 12
lambda_gene = 0.1
lambda_sample = 0.1
lambda_time = 0.05
max_iter = 1000
tol = 1e-5
```

### For Maximum Sparsity (Rank 5)
```python
# Use low gene regularization
rank = 5
lambda_gene = 0.01  # Lower lambda = less sparsity but lower loss
lambda_sample = 0.1
lambda_time = 0.05
max_iter = 1000
tol = 1e-5
```

### For Interpretability (Fewer Components)
```python
# Use Rank 3 or 5
rank = 3  # or 5
lambda_gene = 0.1
lambda_sample = 0.1
lambda_time = 0.05
max_iter = 1000
tol = 1e-5
```

## All Available Analyses

| Analysis | Script Location | Output Location | Key Results |
|----------|----------------|-----------------|-------------|
| **Initial Baseline** | `scripts/13.00-multiomics-barnacle.Rmd` | `output/13.00-multiomics-barnacle/` | Normalized data, initial factors |
| **Refined Approach** | `scripts/14-barnacle/` | `output/14-barnacle/` | Rank 5 with visualizations |
| **Single Run** | `scripts/14.1-barnacle/` | `output/14.1-barnacle/` | Rank 5 single run |
| **Convergence Tests** | `scripts/14.1-barnacle/convergence_test.py` | `output/14.1-barnacle-convergence-tests/` | 94 parameter combinations |
| **Advanced Config** | `scripts/14.5-barnacle/` | `output/14.5-barnacle/` | Advanced preprocessing |
| **Multi-Rank** | `scripts/15-barnacle/` | `output/15-barnacle/` | Ranks 3, 8, 10, 12 |
| **Additional Ranks** | `scripts/15.5-barnacle/` | `output/15.5-barnacle/` | Ranks 3, 8, 10 |
| **HPC (Raven)** | `scripts/16-barnacle-raven/` | - | Cluster implementation |
| **HPC (Klone)** | `scripts/17-barnacle-klone/` | - | Cluster implementation |

## Convergence Test Highlights

### Top 10 Best Performing Tests (Lowest Loss)

| Test ID | Loss | λ_gene | λ_sample | λ_time | Max Iter | Tolerance |
|---------|------|--------|----------|--------|----------|-----------|
| 006 | 840,915 | 0.01 | 0.1 | 0.05 | 1000 | 1e-5 |
| 016 | 840,915 | 0.01 | 0.1 | 0.05 | 2000 | 1e-5 |
| 066 | 840,915 | 0.01 | 0.1 | 0.05 | 5000 | 1e-5 |
| 005 | 845,241 | 0.05 | 0.1 | 0.05 | 1000 | 1e-5 |
| 015 | 845,241 | 0.05 | 0.1 | 0.05 | 2000 | 1e-5 |
| 065 | 845,241 | 0.05 | 0.1 | 0.05 | 5000 | 1e-5 |
| 001 | 850,175 | 0.1 | 0.1 | 0.05 | 1000 | 1e-5 |
| 011 | 850,175 | 0.1 | 0.1 | 0.05 | 2000 | 1e-5 |
| 061 | 850,175 | 0.1 | 0.1 | 0.05 | 5000 | 1e-5 |
| 009 | 850,175 | 0.1 | 0.1 | 0.1 | 1000 | 1e-5 |

**Pattern**: Lower gene regularization (λ_gene=0.01) consistently achieves lowest loss.

### Worst 5 Performing Tests (Highest Loss)

| Test ID | Loss | λ_gene | λ_sample | λ_time | Max Iter | Tolerance |
|---------|------|--------|----------|--------|----------|-----------|
| 094 | 1,006,770 | 1.0 | 0.1 | 0.05 | 10000 | 1e-3 |
| 084 | 959,313 | 1.0 | 0.1 | 0.05 | 10000 | 1e-4 |
| 074 | 911,497 | 1.0 | 0.1 | 0.05 | 5000 | 1e-4 |

**Pattern**: High gene regularization (λ_gene=1.0) combined with relaxed tolerance results in highest loss.

## Convergence Status Summary

- **Total tests**: 94
- **Converged**: 0 (0%)
- **Successful runs**: 94 (100%)
- **Loss range**: 840,915 - 1,006,770

## Visualization Files

### Available Figures by Analysis

| Analysis | Component Weights | Timepoint Loadings |
|----------|------------------|-------------------|
| **14-barnacle** | ✓ | ✓ |
| **15-barnacle/rank-3** | ✓ | ✓ |
| **15-barnacle/rank-8** | ✓ | ✓ |
| **15-barnacle/rank-10** | ✓ | ✓ |
| **15-barnacle/rank-12** | ✓ | ✓ |
| **15.5-barnacle/rank-3** | ✓ | ✓ |
| **15.5-barnacle/rank-8** | ✓ | ✓ |
| **15.5-barnacle/rank-10** | ✓ | ✓ |

All figures are PNG format located in respective `figures/` subdirectories.

## Factor Matrix Files

Each analysis produces the following factor matrices:

- `gene_factors.csv` - Gene loadings (10,223 × rank)
- `sample_factors.csv` - Sample loadings (30 × rank)
- `time_factors.csv` - Timepoint loadings (4 × rank)
- `component_weights.csv` - Component importance weights
- `sample_mapping.csv` - Sample to species mapping
- `metadata.json` - Run metadata including loss and convergence status

## Command Reference

### Rerun Best Parameter Set (Rank 5)

```bash
uv run python M-multi-species/scripts/14.1-barnacle/build_tensor_and_run.py \
  --input-file M-multi-species/output/14-pca-orthologs/vst_counts_matrix.csv \
  --output-dir M-multi-species/output/my-barnacle-run \
  --rank 5 --lambda-gene 0.01 --lambda-sample 0.1 --lambda-time 0.05 \
  --max-iter 1000 --tol 1e-5 --seed 42
```

### Rerun Best Rank (Rank 12)

```bash
uv run python M-multi-species/scripts/15-barnacle/build_tensor_and_run.py \
  --input-dir M-multi-species/output/13.00-multiomics-barnacle \
  --output-dir M-multi-species/output/my-rank12-run \
  --rank 12 --lambda-gene 0.1 --lambda-sample 0.1 --lambda-time 0.05 \
  --max-iter 1000 --tol 1e-5 --seed 42
```

### Run All Ranks Comparison

```bash
bash M-multi-species/scripts/15-barnacle/run_all.sh
```

---

**For detailed analysis and interpretation**, see: `M-multi-species/output/BARNACLE_ANALYSIS_REPORT.md`
