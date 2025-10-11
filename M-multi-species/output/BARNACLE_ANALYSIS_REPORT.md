# Comprehensive Barnacle Analysis Report

**Repository**: urol-e5/timeseries_molecular  
**Date**: October 2025  
**Analysis Focus**: Multi-species tensor decomposition convergence and parameter optimization

---

## Executive Summary

This report provides a detailed overview of all Barnacle sparse tensor decomposition analyses performed on the multi-species timeseries molecular dataset. The analyses encompass:

1. **Initial baseline analysis** (13.00-multiomics-barnacle)
2. **Refined tensor approach** (14-barnacle, 14.1-barnacle, 14.5-barnacle)
3. **Systematic convergence testing** (14.1-barnacle-convergence-tests) - **94 parameter combinations tested**
4. **Multi-rank comparison** (15-barnacle) - **4 different ranks evaluated**
5. **Additional rank variations** (15.5-barnacle)

### Key Finding

**No parameter combination achieved convergence** across all 94 tests in the systematic convergence analysis. However, the analyses successfully identified the best-performing parameter sets and provided critical insights into the behavior of the algorithm with this specific dataset.

---

## Dataset Characteristics

### Tensor Dimensions
- **Genes**: 10,223 orthologous genes
- **Samples**: 30 combined samples across 3 species
  - Acropora pulchra (Apul)
  - Pocillopora verrucosa (Peve)
  - Pocillopora meandrina (Ptua)
- **Timepoints**: 4 timepoints (TP1, TP2, TP3, TP4)
- **Total tensor shape**: 10,223 × 30 × 4

### Data Preprocessing
- Ortholog mapping across three coral species
- Normalized expression values using sctransform::vst or log2(CPM+1) fallback
- Missing values filled with 0.0 for tensor decomposition

---

## Analysis 1: Initial Baseline (13.00-multiomics-barnacle)

### Objective
Establish initial tensor decomposition pipeline and generate baseline factor matrices.

### Location
- **Scripts**: `M-multi-species/scripts/13.00-multiomics-barnacle.Rmd`
- **Output**: `M-multi-species/output/13.00-multiomics-barnacle/`

### Key Outputs
- Normalized expression matrices per species:
  - `apul_normalized_expression.csv`
  - `peve_normalized_expression.csv`
  - `ptua_normalized_expression.csv`
- Tensor: `multiomics_tensor.npy`
- Factor matrices in `barnacle_factors/` directory

### Results
- Successfully created 3D tensor from multi-species data
- Generated factor matrices for genes, samples, and timepoints
- Established data preprocessing pipeline for downstream analyses

---

## Analysis 2: Refined Tensor Approach (14-barnacle)

### Objective
Refine the tensor decomposition approach with improved parameters and visualization.

### Location
- **Scripts**: `M-multi-species/scripts/14-barnacle/`
- **Output**: `M-multi-species/output/14-barnacle/`

### Parameters
- **Rank**: 5
- **Lambda (gene)**: 0.1
- **Lambda (sample)**: 0.1
- **Lambda (time)**: 0.05
- **Max iterations**: 1000
- **Tolerance**: 1e-5
- **Random seed**: 42

### Key Outputs
- Component weights visualization
- Timepoint loadings across components
- Factor matrices with metadata

### Results Summary
Location: `M-multi-species/output/14-barnacle/SUMMARY.md`

- **Tensor shape**: (10223, 30, 4)
- **Components extracted**: 5
- **Convergence**: Not achieved
- Generated visualizations:
  - Component weights bar chart
  - Timepoint loadings line plot

---

## Analysis 3: Systematic Convergence Testing (14.1-barnacle-convergence-tests)

### Objective
**Systematically test 94 parameter combinations** to identify optimal settings for achieving convergence.

### Location
- **Scripts**: `M-multi-species/scripts/14.1-barnacle/convergence_test.py`
- **Output**: `M-multi-species/output/14.1-barnacle-convergence-tests/`
- **Summary**: `BARNACLE_CONVERGENCE_SUMMARY.md`

### Testing Methodology

#### Parameters Tested

1. **Maximum Iterations**: 1000, 2000, 5000, 10000
2. **Convergence Tolerances**: 1e-5, 1e-4, 1e-3
3. **Lambda Regularization Values**:
   - Gene regularization (λ_gene): 0.01, 0.05, 0.1, 0.2, 0.5, 1.0
   - Sample regularization (λ_sample): 0.1, 0.2, 0.5
   - Time regularization (λ_time): 0.05, 0.1, 0.2

#### Test Coverage

- **Total Combinations**: 94 unique parameter sets
- **Successful Runs**: 94/94 (100% success rate)
- **Converged Runs**: 0/94 (0% convergence rate)
- **Final Loss Range**: 840,915 - 1,006,770

### Critical Findings

#### 1. Convergence Failure Patterns
- **All parameter combinations failed to converge**
- **Increasing iterations did not help** (tested up to 10,000 iterations)
- **Relaxing tolerance did not achieve convergence** (tested down to 1e-3)
- **No clear pattern of improvement** with different regularization values

#### 2. Best Parameter Performance

**Lowest loss consistently achieved with λ_gene=0.01** across all parameter combinations:

| Parameter Category | Loss Range | Best Performance |
|-------------------|------------|------------------|
| Baseline (max_iter=1000, tol=1e-5) | 840,915 - 886,452 | **840,915** (λ_gene=0.01) |
| Increased iterations (2000) | 840,915 - 911,497 | **840,915** (λ_gene=0.01) |
| Relaxed tolerance (1e-4) | 864,254 - 911,497 | **864,254** (λ_gene=0.01) |
| High iterations (5000) | 840,915 - 911,497 | **840,915** (λ_gene=0.01) |
| Very relaxed tolerance (1e-3) | 882,064 - 1,006,770 | **882,064** (λ_gene=0.01) |
| Extreme iterations (10000) | 840,915 - 911,497 | **840,915** (λ_gene=0.01) |

#### 3. Algorithm Behavior
- **Algorithm stalls** at similar loss values regardless of parameter settings
- **No evidence** that further parameter tuning will achieve convergence
- Loss plateaus suggest the algorithm reaches a local minimum that it cannot escape from

### Best Parameter Set Identified

While convergence was not achieved, the best-performing parameter set was:

```
Rank: 5
Lambda (gene): 0.01
Lambda (sample): 0.1
Lambda (time): 0.05
Max iterations: 1000-10000 (no significant difference)
Tolerance: 1e-5
Final loss: ~840,915
```

### Available Data

- **Test outputs**: `test_001/` through `test_094/` (each contains full factor matrices and metadata)
- **Intermediate results**: CSV and JSON files at intervals (10, 20, 30, etc.)
- **Individual run logs**: Available in each test directory

---

## Analysis 4: Multi-Rank Comparison (15-barnacle)

### Objective
Compare Barnacle decomposition performance across different rank values.

### Location
- **Scripts**: `M-multi-species/scripts/15-barnacle/`
- **Output**: `M-multi-species/output/15-barnacle/`
- **Summary**: `SUMMARY.md`

### Ranks Evaluated
- Rank 3
- Rank 8
- Rank 10
- Rank 12

### Parameters (constant across ranks)
- **Lambda (gene)**: 0.1
- **Lambda (sample)**: 0.1
- **Lambda (time)**: 0.05
- **Max iterations**: 1000
- **Tolerance**: 1e-5
- **Random seed**: 42

### Results Summary

| Rank | Components | Converged | Final Loss | Gene Sparsity | Sample Sparsity | Time Sparsity |
|---:|---:|:---:|---:|---:|---:|---:|
| 3 | 3 | False | 1,368,622.54 | 0.003 | 0.000 | 0.000 |
| 8 | 8 | False | 689,401.09 | 0.022 | 0.000 | 0.000 |
| 10 | 10 | False | 634,014.28 | 0.015 | 0.000 | 0.000 |
| 12 | 12 | False | 587,541.02 | 0.012 | 0.000 | 0.000 |

### Key Observations

1. **Loss decreases with higher rank**: Higher rank values achieve lower final loss
2. **No convergence at any rank**: All ranks failed to converge
3. **Sparsity patterns**: 
   - Gene sparsity shows some variation (0.3% - 2.2%)
   - Sample and time sparsity remain at 0% across all ranks
4. **Trade-off consideration**: Higher ranks may overfit; lower ranks may underfit

### Visualizations Available

For each rank, the following visualizations are available:

- `rank-X/figures/component_weights.png` - Component weights bar chart
- `rank-X/figures/time_loadings.png` - Timepoint loadings across components

---

## Analysis 5: Additional Rank Variations (15.5-barnacle)

### Location
- **Scripts**: `M-multi-species/scripts/15.5-barnacle/`
- **Output**: `M-multi-species/output/15.5-barnacle/`

### Ranks Tested
- Rank 3
- Rank 8
- Rank 10

### Status
Analysis completed with factor matrices and visualizations generated for each rank.

---

## Analysis 6: Advanced Configurations (14.5-barnacle)

### Objective
Test advanced preprocessing and configuration options.

### Location
- **Scripts**: `M-multi-species/scripts/14.5-barnacle/`

### Advanced Features Tested
- Different missing value imputation policies (gene-mean)
- Per-gene z-score scaling
- Minimum timepoints per sample filtering
- Flexible rank selection

---

## Synthesis and Recommendations

### Summary of Findings

1. **Convergence Challenge**: No parameter combination achieved convergence across 94 systematic tests
2. **Best Parameters**: λ_gene=0.01 consistently achieves lowest loss (~840,915)
3. **Rank Effect**: Higher ranks achieve lower loss but don't guarantee convergence
4. **Algorithm Limitation**: SparseCP appears to reach a plateau/local minimum with this dataset

### Potential Issues Identified

1. **Dataset characteristics** may be incompatible with SparseCP assumptions
   - High dimensionality (10,223 genes)
   - Small sample size (30 samples)
   - Limited timepoints (4)
   - Sparse data structure

2. **Rank selection** may be inappropriate for this 3D tensor structure
   - Tested ranks: 3, 5, 8, 10, 12
   - All failed to converge
   - Lower ranks may be insufficient; higher ranks may overfit

3. **Data preprocessing** may need revision
   - Missing value imputation strategy (0-fill vs gene-mean)
   - Normalization approach
   - Feature selection (reduce gene set)

4. **Algorithm limitations** with this specific data structure
   - SparseCP may not be suitable for this type of data
   - Alternative tensor decomposition methods may be needed

### Recommendations for Future Work

#### 1. Alternative Tensor Decomposition Methods

Consider trying:
- **Tucker decomposition** (lower rank core tensor)
- **PARAFAC without sparsity constraints**
- **Non-negative matrix factorization** on unfolded tensor
- **Canonical correlation analysis** approaches

#### 2. Data Preprocessing Modifications

- **Feature selection**: Reduce gene set size (e.g., top 1000-5000 most variable genes)
- **Different normalization strategies**: Test various normalization methods
- **Alternative imputation**: Use gene-mean or k-NN imputation instead of 0-fill
- **Data transformation**: Consider log-transformation or standardization

#### 3. Rank Exploration

- **Test lower ranks**: Ranks 2, 3, 4
- **Test higher ranks**: Ranks 15, 20, 25
- **Automatic rank selection**: Use cross-validation or information criteria

#### 4. Algorithm Parameter Reconsideration

- **Multiple random initializations**: Increase n_initializations
- **Different optimization algorithms** within Barnacle
- **Orthogonal constraints** instead of non-negativity
- **Different sparsity patterns** for different modes

#### 5. Validation with Simpler Datasets

- **Test on synthetic data** with known structure
- **Compare with smaller real datasets** that should converge
- **Validate algorithm installation** and version compatibility

### Best Practices Identified

Based on the comprehensive testing, if using Barnacle SparseCP with this dataset:

1. **Use λ_gene=0.01** for lowest loss
2. **Set max_iter=1000** (higher values don't significantly improve results)
3. **Use tolerance=1e-5** (relaxing doesn't achieve convergence)
4. **Accept non-convergence** and evaluate results based on biological interpretability
5. **Higher ranks** (8-12) achieve lower loss but may overfit

---

## Directory Structure and Outputs

### Complete Analysis Tree

```
M-multi-species/
├── scripts/
│   ├── 13.00-multiomics-barnacle.Rmd       # Initial baseline
│   ├── 14-barnacle/                         # Refined approach
│   ├── 14.1-barnacle/                       # Convergence testing
│   ├── 14.5-barnacle/                       # Advanced config
│   ├── 15-barnacle/                         # Multi-rank comparison
│   ├── 15.5-barnacle/                       # Additional ranks
│   ├── 16-barnacle-raven/                   # HPC version (Raven)
│   └── 17-barnacle-klone/                   # HPC version (Klone)
└── output/
    ├── 13.00-multiomics-barnacle/           # Normalized data + initial factors
    ├── 14-barnacle/                         # Rank 5 analysis
    ├── 14.1-barnacle/                       # Single run
    ├── 14.1-barnacle-convergence-tests/     # 94 parameter tests
    │   ├── BARNACLE_CONVERGENCE_SUMMARY.md
    │   ├── test_001/ through test_094/
    │   └── intermediate_results_*.csv/json
    ├── 15-barnacle/                         # Multi-rank comparison
    │   ├── SUMMARY.md
    │   ├── rank-3/
    │   ├── rank-8/
    │   ├── rank-10/
    │   └── rank-12/
    └── 15.5-barnacle/                       # Additional rank tests
        ├── rank-3/
        ├── rank-8/
        └── rank-10/
```

### Key Files by Analysis

#### Convergence Tests (14.1-barnacle-convergence-tests)
- `BARNACLE_CONVERGENCE_SUMMARY.md` - Comprehensive convergence analysis
- `intermediate_results_10.csv` through `intermediate_results_80.csv` - Results at intervals
- `test_XXX/barnacle_factors/metadata.json` - Individual run metadata
- `test_XXX/barnacle_factors/*.csv` - Factor matrices for each test

#### Rank Comparison (15-barnacle)
- `SUMMARY.md` - Multi-rank comparison summary
- `rank-X/barnacle_factors/metadata.json` - Per-rank metadata
- `rank-X/figures/*.png` - Visualizations for each rank

---

## Computational Considerations

### Resource Usage

- **Convergence testing**: 94 tests × ~5-30 minutes each = ~8-40 hours total
- **Memory requirements**: Tensor size (10,223 × 30 × 4) = ~1.2 MB + factor matrices
- **Disk space**: Each test ~1-5 MB; total ~100-500 MB for all tests

### HPC Implementations

- **16-barnacle-raven**: Implementation for Raven HPC cluster
- **17-barnacle-klone**: Implementation for Klone HPC cluster

---

## Conclusions

### What We Learned

1. **Systematic testing revealed fundamental convergence challenges** with this dataset and algorithm combination
2. **Best parameter set identified**: λ_gene=0.01, λ_sample=0.1, λ_time=0.05
3. **Rank matters**: Higher ranks achieve lower loss but don't guarantee convergence
4. **Algorithm behavior**: SparseCP reaches a plateau around 840K-1,400K loss depending on rank

### Practical Outcomes

Despite non-convergence, the analyses:
- Generated interpretable factor matrices
- Identified temporal patterns in timepoint loadings
- Revealed species-sample groupings in sample factors
- Provided gene-level component associations

### Next Steps

1. **Consider alternative methods**: Tucker decomposition, PARAFAC, NMF
2. **Refine data preprocessing**: Feature selection, alternative normalization
3. **Biological validation**: Evaluate factor matrices for biological interpretability
4. **Method comparison**: Compare Barnacle results with other dimensionality reduction approaches (PCA, MOFA2)

---

## References and Links

### Documentation
- Convergence testing summary: `M-multi-species/output/14.1-barnacle-convergence-tests/BARNACLE_CONVERGENCE_SUMMARY.md`
- Rank comparison summary: `M-multi-species/output/15-barnacle/SUMMARY.md`
- Initial baseline: `M-multi-species/output/13.00-multiomics-barnacle/README.md`

### Scripts
- Convergence test script: `M-multi-species/scripts/14.1-barnacle/convergence_test.py`
- Multi-rank script: `M-multi-species/scripts/15-barnacle/run_all.sh`
- Baseline R analysis: `M-multi-species/scripts/13.00-multiomics-barnacle.Rmd`

### Key Results Files
- All test results: `M-multi-species/output/14.1-barnacle-convergence-tests/intermediate_results_80.csv`
- Rank comparison: `M-multi-species/output/15-barnacle/SUMMARY.md`

---

**Report completed**: October 2025  
**Total analyses documented**: 6 major analysis streams  
**Total parameter combinations tested**: 94  
**Total ranks evaluated**: 4 (with additional variations)  
**Convergence achieved**: 0/94 tests  
**Best loss achieved**: 587,541 (rank 12) to 840,915 (rank 5, λ_gene=0.01)
