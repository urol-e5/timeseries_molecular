#!/bin/bash
# Example script demonstrating the complete workflow for synthetic data generation
# and tensor decomposition with barnacle

set -e  # Exit on error

echo "=========================================================================="
echo "Barnacle Synthetic Data - Complete Workflow Example"
echo "=========================================================================="

# Configuration
SCRIPTS_DIR="M-multi-species/scripts/14-barnacle"
OUTPUT_BASE="M-multi-species/output"
SYNTHETIC_DIR="${OUTPUT_BASE}/14-barnacle-synthetic"
RESULTS_DIR="${OUTPUT_BASE}/14-barnacle-synthetic-results"

# Step 1: Generate synthetic data
echo ""
echo "Step 1: Generating synthetic data..."
echo "----------------------------------------------------------------------"
python3 ${SCRIPTS_DIR}/generate_synthetic_data.py \
  --output-dir ${SYNTHETIC_DIR} \
  --n-genes 10223 \
  --n-samples-per-species 10 \
  --n-components 5 \
  --noise-level 0.1 \
  --seed 42

# Step 2: Validate synthetic data
echo ""
echo "Step 2: Validating synthetic data..."
echo "----------------------------------------------------------------------"
python3 ${SCRIPTS_DIR}/test_synthetic_data.py

# Step 3: Run tensor decomposition (requires barnacle installation)
echo ""
echo "Step 3: Running tensor decomposition with Barnacle..."
echo "----------------------------------------------------------------------"
echo "NOTE: This requires barnacle to be installed via:"
echo "  uv pip install git+https://github.com/blasks/barnacle.git@612b6a4"
echo ""

# Uncomment the following lines after installing barnacle:
# uv run python ${SCRIPTS_DIR}/build_tensor_and_run.py \
#   --input-dir ${SYNTHETIC_DIR} \
#   --output-dir ${RESULTS_DIR} \
#   --rank 5 \
#   --lambda-gene 0.1 \
#   --lambda-sample 0.1 \
#   --lambda-time 0.05 \
#   --max-iter 2000 \
#   --tol 1e-4 \
#   --seed 42

echo ""
echo "Step 4: Inspect results (after decomposition)..."
echo "----------------------------------------------------------------------"
echo "Results will be in: ${RESULTS_DIR}"
echo ""
echo "Files to check:"
echo "  - ${RESULTS_DIR}/barnacle_factors/metadata.json (convergence status)"
echo "  - ${RESULTS_DIR}/barnacle_factors/gene_factors.csv"
echo "  - ${RESULTS_DIR}/barnacle_factors/sample_factors.csv"
echo "  - ${RESULTS_DIR}/barnacle_factors/time_factors.csv"
echo "  - ${RESULTS_DIR}/figures/component_weights.png"
echo "  - ${RESULTS_DIR}/figures/time_loadings.png"
echo ""
echo "Ground truth for comparison:"
echo "  - ${SYNTHETIC_DIR}/ground_truth/true_gene_factors.csv"
echo "  - ${SYNTHETIC_DIR}/ground_truth/true_sample_factors.csv"
echo "  - ${SYNTHETIC_DIR}/ground_truth/true_time_factors.csv"

echo ""
echo "=========================================================================="
echo "Workflow complete!"
echo "=========================================================================="
