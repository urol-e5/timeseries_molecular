# Sparse Tensor Decomposition Analysis Report

## Run Information

- **Timestamp**: 2025-10-15 12:28:43
- **Input File**: ../output/13.00-multiomics-stdm/vst_counts_matrix_long_format.csv
- **Output Directory**: ../output/13.00-multiomics-stdm/20251015_122334/rank_comparison/rank_09/20251015_122804
- **Method**: PARAFAC
- **Rank**: 9

## Input Data Summary

- **Shape**: (9800, 30, 4) (genes × species × timepoints)
- **Number of Genes**: 9,800
- **Number of Species**: 30
- **Number of Timepoints**: 4
- **Sparsity**: 2.50%
- **Value Range**: [0.0000, 18.4240]
- **Mean Value**: 7.2089
- **Std Deviation**: 2.3063

## Analysis Parameters

- **Decomposition Method**: PARAFAC
- **Rank**: 9
- **Sparsity Threshold**: 0.01
- **Normalization**: False
- **Log Transform**: False
- **Standardization**: True

## Results Summary

- **Reconstruction Error**: 0.401437
- **Quality Assessment**: Acceptable
- **Gene Factors Shape**: (9800, 9)
- **Species Factors Shape**: (30, 9)
- **Time Factors Shape**: (4, 9)

## Output Files

This directory contains the following files:

### 1. **gene_factors.npy**
- **Description**: Gene factor matrix (9,800 × 9)
- **Usage**: Each row represents a gene's loading on each component
- **Interpretation**: Higher absolute values indicate stronger association with that component
- **Load in Python**: `gene_factors = np.load('gene_factors.npy')`

### 2. **species_factors.npy**
- **Description**: Species factor matrix (30 × 9)
- **Usage**: Shows how each species contributes to each component
- **Interpretation**: Components with high loading in one species indicate species-specific patterns
- **Load in Python**: `species_factors = np.load('species_factors.npy')`

### 3. **time_factors.npy**
- **Description**: Time point factor matrix (4 × 9)
- **Usage**: Shows temporal patterns for each component
- **Interpretation**: Positive/negative trends indicate up/down regulation over time
- **Load in Python**: `time_factors = np.load('time_factors.npy')`

### 4. **reconstructed_tensor.npy**
- **Description**: Reconstructed tensor from decomposition (9,800 × 30 × 4)
- **Usage**: Can be compared with original data for validation
- **Load in Python**: `reconstructed = np.load('reconstructed_tensor.npy')`

### 5. **summary.json**
- **Description**: JSON file containing decomposition statistics
- **Contents**: Method, rank, reconstruction error, and factor shapes
- **Load in Python**: 
  ```python
  import json
  with open('summary.json') as f:
      summary = json.load(f)
  ```

### 6. **README.md** (this file)
- **Description**: This comprehensive report documenting the analysis

## Quality Assessment

**Reconstruction Error**: 0.401437 (Acceptable)

The fit is acceptable. Consider increasing the rank to capture more variation.

### Interpretation Guide

#### Reconstruction Error Ranges:
- **0.0 - 0.1**: Excellent fit (may be overfitting)
- **0.1 - 0.3**: Good fit
- **0.3 - 0.5**: Acceptable fit
- **> 0.5**: Poor fit (increase rank or check data quality)

## Recommendations

### Parameter Optimization

**Rank Selection**: Current rank (9) appears appropriate

**Method Comparison**:
- **PARAFAC/CP**: Best for identifying simple, interpretable patterns
- **Tucker**: Better for complex, hierarchical patterns (but more parameters)

### Next Steps

1. **Examine Factor Matrices**:
   ```python
   import numpy as np
   
   # Load factors
   gene_factors = np.load('gene_factors.npy')
   species_factors = np.load('species_factors.npy')
   time_factors = np.load('time_factors.npy')
   
   # Identify top genes for each component
   for component in range(9):
       top_genes_idx = np.argsort(np.abs(gene_factors[:, component]))[-10:][::-1]
       print(f"Component {component + 1} - Top 10 genes: {top_genes_idx}")
   ```

2. **Validate Results**:
   ```python
   # Compare original vs reconstructed
   original = np.load('path/to/original_tensor.npy')
   reconstructed = np.load('reconstructed_tensor.npy')
   
   diff = np.abs(original - reconstructed)
   print(f"Mean absolute difference: {np.mean(diff):.6f}")
   print(f"Max absolute difference: {np.max(diff):.6f}")
   ```

3. **Biological Interpretation**:
   - Look for gene modules with similar factor profiles (co-regulated genes)
   - Identify species-specific regulatory patterns
   - Examine temporal dynamics (up/down regulation over time)

## Additional Analysis Ideas

- **Gene Set Enrichment**: Use top genes from each component for pathway analysis
- **Clustering**: Cluster genes based on their factor profiles
- **Visualization**: Create heatmaps of factor matrices to identify patterns
- **Cross-validation**: Run decomposition multiple times with different parameters

## How to Load and Use Results

```python
import numpy as np
import json

# Load all results
gene_factors = np.load('gene_factors.npy')
species_factors = np.load('species_factors.npy')
time_factors = np.load('time_factors.npy')
reconstructed = np.load('reconstructed_tensor.npy')

with open('summary.json') as f:
    summary = json.load(f)

print(f"Reconstruction error: {summary['reconstruction_error']:.6f}")

# Example: Get top 20 genes for component 1
component_idx = 0
top_genes = np.argsort(np.abs(gene_factors[:, component_idx]))[-20:][::-1]
print(f"Top 20 genes for component {component_idx + 1}: {top_genes}")

# Example: Examine temporal pattern for a component
import matplotlib.pyplot as plt
plt.plot(time_factors[:, component_idx])
plt.xlabel('Time Point')
plt.ylabel('Factor Value')
plt.title(f'Temporal Pattern - Component {component_idx + 1}')
plt.show()
```

## Citation

If you use these results in a publication, please cite the workflow-stdm package.

---

*Report generated automatically by workflow-stdm*
*Timestamp: 2025-10-15 12:28:43*
