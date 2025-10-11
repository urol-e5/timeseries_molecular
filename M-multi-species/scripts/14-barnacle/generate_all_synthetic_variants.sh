#!/bin/bash
# Script to generate different variations of synthetic data for testing convergence
# with different parameter configurations

set -e

SCRIPT_DIR="M-multi-species/scripts/14-barnacle"
OUTPUT_BASE="M-multi-species/output"

echo "=========================================================================="
echo "Generating Multiple Synthetic Datasets for Convergence Testing"
echo "=========================================================================="

# Configuration 1: Default (High convergence probability)
echo ""
echo "1. Generating DEFAULT synthetic dataset..."
echo "   - 10,223 genes, 10 samples/species, 5 components"
echo "   - Noise: 0.1 (10% of signal)"
echo "--------------------------------------------------------------------------"
python3 ${SCRIPT_DIR}/generate_synthetic_data.py \
  --output-dir ${OUTPUT_BASE}/14-barnacle-synthetic-default \
  --n-genes 10223 \
  --n-samples-per-species 10 \
  --n-components 5 \
  --noise-level 0.1 \
  --seed 42

# Configuration 2: Small dataset (for quick testing)
echo ""
echo "2. Generating SMALL synthetic dataset..."
echo "   - 1,000 genes, 5 samples/species, 3 components"
echo "   - Noise: 0.05 (5% of signal)"
echo "--------------------------------------------------------------------------"
python3 ${SCRIPT_DIR}/generate_synthetic_data.py \
  --output-dir ${OUTPUT_BASE}/14-barnacle-synthetic-small \
  --n-genes 1000 \
  --n-samples-per-species 5 \
  --n-components 3 \
  --noise-level 0.05 \
  --seed 123

# Configuration 3: High noise (harder convergence)
echo ""
echo "3. Generating HIGH-NOISE synthetic dataset..."
echo "   - 10,223 genes, 10 samples/species, 5 components"
echo "   - Noise: 0.3 (30% of signal)"
echo "--------------------------------------------------------------------------"
python3 ${SCRIPT_DIR}/generate_synthetic_data.py \
  --output-dir ${OUTPUT_BASE}/14-barnacle-synthetic-noisy \
  --n-genes 10223 \
  --n-samples-per-species 10 \
  --n-components 5 \
  --noise-level 0.3 \
  --seed 456

# Configuration 4: More components (more complex)
echo ""
echo "4. Generating COMPLEX synthetic dataset..."
echo "   - 10,223 genes, 10 samples/species, 8 components"
echo "   - Noise: 0.1 (10% of signal)"
echo "--------------------------------------------------------------------------"
python3 ${SCRIPT_DIR}/generate_synthetic_data.py \
  --output-dir ${OUTPUT_BASE}/14-barnacle-synthetic-complex \
  --n-genes 10223 \
  --n-samples-per-species 10 \
  --n-components 8 \
  --noise-level 0.1 \
  --seed 789

# Configuration 5: More samples (better for convergence)
echo ""
echo "5. Generating LARGE-SAMPLE synthetic dataset..."
echo "   - 10,223 genes, 15 samples/species, 5 components"
echo "   - Noise: 0.1 (10% of signal)"
echo "--------------------------------------------------------------------------"
python3 ${SCRIPT_DIR}/generate_synthetic_data.py \
  --output-dir ${OUTPUT_BASE}/14-barnacle-synthetic-large \
  --n-genes 10223 \
  --n-samples-per-species 15 \
  --n-components 5 \
  --noise-level 0.1 \
  --seed 999

echo ""
echo "=========================================================================="
echo "Summary of Generated Datasets"
echo "=========================================================================="
echo ""
echo "Dataset Name          | Genes  | Samples/Sp | Components | Noise"
echo "----------------------|--------|------------|------------|-------"
echo "default               | 10,223 | 10         | 5          | 0.10"
echo "small                 | 1,000  | 5          | 3          | 0.05"
echo "noisy                 | 10,223 | 10         | 5          | 0.30"
echo "complex               | 10,223 | 10         | 8          | 0.10"
echo "large                 | 10,223 | 15         | 5          | 0.10"
echo ""
echo "All datasets saved to: ${OUTPUT_BASE}/14-barnacle-synthetic-*"
echo ""
echo "Recommended Testing Order:"
echo "  1. 'small' - Quick validation that pipeline works"
echo "  2. 'default' - Best convergence probability"
echo "  3. 'large' - Even better convergence with more samples"
echo "  4. 'noisy' - Test robustness to noise"
echo "  5. 'complex' - Test with more components"
echo ""
echo "To test any dataset, run:"
echo "  python3 ${SCRIPT_DIR}/test_synthetic_data.py <dataset-name>"
echo ""
echo "=========================================================================="
