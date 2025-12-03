# Primary Workflow

(main explanatory documentation of process; tangential details and complementary information can be found in other files in this directory)

# Expression Data

Count data from <https://urol-e5.github.io/MOSAiC/counts.html>

``` bash
# Load ortholog annotations
orthos <- read_csv("https://raw.githubusercontent.com/urol-e5/timeseries_molecular/refs/heads/main/M-multi-species/output/12-ortho-annot/ortholog_groups_annotated.csv")

# Load species-specific count matrices
apul <- read_csv("https://gannet.fish.washington.edu/gitrepos/urol-e5/timeseries_molecular/D-Apul/output/02.20-D-Apul-RNAseq-alignment-HiSat2/apul-gene_count_matrix.csv")
peve <- read_csv("https://gannet.fish.washington.edu/gitrepos/urol-e5/timeseries_molecular/E-Peve/output/02.20-E-Peve-RNAseq-alignment-HiSat2/peve-gene_count_matrix.csv")
ptua <- read_csv("https://gannet.fish.washington.edu/gitrepos/urol-e5/timeseries_molecular/F-Ptua/output/02.20-F-Ptua-RNAseq-alignment-HiSat2/ptua-gene_count_matrix.csv")
```

> Ortholog groups were selected and merged across species to create a combined count matrix containing 9,800 orthologs present at non-NA values in all 117 samples. Metadata was parsed from sample names to extract species and timepoint information.

## Variance Stabilizing Transformation (VST)

Low-count and incomplete ortholog rows were removed. A variance stabilizing transformation was applied using DESeq2::vst() to normalize count variance across the dynamic range blind to study design.

``` r
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData   = metadata,
  design    = ~ species * timepoint
)

vst_obj <- vst(dds, blind = TRUE)
vst_mat <- assay(vst_obj)
```

```         
../output/14-pca-orthologs/vst_counts_matrix.csv
```

<img src="http://gannet.fish.washington.edu/seashell/snaps/Monosnap_Image_2025-11-26_05-34-30.png" style="width:50%;"/>

# Running Barnacle

## SR version

```         
M-multi-species/scripts/14.1-barnacle/build_tensor_and_run.py
```

[Code on GitHub](https://github.com/urol-e5/timeseries_molecular/blob/main/M-multi-species/scripts/14.1-barnacle/build_tensor_and_run.py)

*example*

``` bash
uv run python M-multi-species/scripts/14.1-barnacle/build_tensor_and_run.py \
  --input-file M-multi-species/output/14-pca-orthologs/vst_counts_matrix.csv \
  --output-dir M-multi-species/output/14.1-barnacle \
  --rank 5 --lambda-gene 0.1 --lambda-sample 0.1 --lambda-time 0.05 \
  --max-iter 1000 --tol 1e-5 --seed 42
```

This script performs **tensor decomposition analysis** on multi-species time-series molecular data using the Barnacle SparseCP algorithm. Here's what it does:

1.  **Reads normalized gene expression data** from a CSV file containing multiple coral species samples across timepoints

2.  **Builds a 3D tensor** with dimensions:

    -   Genes (features)

    -   Samples (biological replicates grouped by species)

    -   Time (timepoints 1-4)

3.  **Runs Sparse CP (CANDECOMP/PARAFAC) decomposition** to identify latent patterns across all three dimensions simultaneously

4.  **Saves results** including factor matrices, visualizations, and metadata

## **Key Functions**

-   [**parse_sample_timepoint()**](vscode-file://vscode-app/Applications/Visual%20Studio%20Code.app/Contents/Resources/app/out/vs/code/electron-browser/workbench/workbench.html): Extracts species, sample ID, and timepoint from column names like "ACR-139-TP1"

    -   ACR → *A. pulchra* (apul)

    -   POR → *P. evermanni* (peve)

    -   POC → *P. tuahiniensis* (ptua)

-   [**build_tensor()**](vscode-file://vscode-app/Applications/Visual%20Studio%20Code.app/Contents/Resources/app/out/vs/code/electron-browser/workbench/workbench.html): Constructs the 3D array from the CSV, organizing samples by species blocks and filling in timepoint data

-   [**run_sparse_cp()**](vscode-file://vscode-app/Applications/Visual%20Studio%20Code.app/Contents/Resources/app/out/vs/code/electron-browser/workbench/workbench.html): Applies sparse tensor decomposition with L1 regularization penalties (lambdas) on each mode to find low-rank patterns

-   [**save_outputs()**](vscode-file://vscode-app/Applications/Visual%20Studio%20Code.app/Contents/Resources/app/out/vs/code/electron-browser/workbench/workbench.html): Exports:

    -   Factor matrices (gene/sample/time loadings per component)

    -   Component weights

    -   Figures showing component importance and temporal patterns

    -   Metadata (convergence, loss, etc.)

## Barnacle: Rank Determination

## Barnacle: Optimization

---

## SAM'S METHOD

### 1. Full Dissertation Grid Search Function

**Function:** `dissertation_grid_search_cv()`

**Location:** Lines ~1360-1750 in `13.00-multiomics-barnacle.Rmd`

**Implementation Details:**

```python
dissertation_grid_search_cv(
    tensor,
    replicate_groups,
    rank_range=[5, 10, 15, 20, 25, 30],
    lambda_values=[0.0, 0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1.0],
    n_iter_max=10000,
    random_state=42
)
```

**What it does:**

1. **Grid Search**: Tests ALL rank × lambda combinations
   - For each combination (R, λ):
     - Fits model to each CV fold (leave-one-group-out)
     - Calculates SSE on held-out data
     - Stores decomposition for FMS calculation

2. **SSE Calculation**: 
   - Each fitted model evaluated against ALL held-out groups
   - Creates cross-validated SSE scores (as dissertation describes)

3. **FMS Calculation**:
   - Pairwise Factor Match Score between successful fold models
   - Compares gene and time factors only (sample factors have dimension mismatch)
   - Creates cross-validated FMS scores (as dissertation describes)

4. **Two-Stage Selection** (exact dissertation method):
   - **Stage 1 - Rank Selection**: 
     - Filter to λ=0.0 models only
     - Select R at minimum mean CV-SSE
   - **Stage 2 - Lambda Selection**:
     - Filter to optimal rank only
     - Find maximum FMS and its standard error
     - Calculate 1SE threshold: `max_FMS - SE(max_FMS)`
     - Select maximum λ where FMS ≥ threshold

**Returns:**
- `optimal_rank`: Selected rank from Stage 1
- `optimal_lambda`: Selected lambda from Stage 2
- `grid_results_df`: Full DataFrame with all combinations and metrics

### RESULTS

- Rank 35 identified as ideal rank
  - Seems to select highest rank, since SSE appears to decrease with rank. Not useful

  <img width="3572" height="1767" alt="image" src="https://github.com/user-attachments/assets/77708c54-9323-4fef-bc0d-edd3d72725e3" />


