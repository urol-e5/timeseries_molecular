# Physiology Count Matrices Summary

This directory contains count matrices for physiology data (lipidomics and metabolomics) generated for the E5 timeseries project.

## Files Generated

### Combined Matrices (All Species)
- `combined-lipidomics-count-matrix.csv`: 371 lipids × 96 samples
- `combined-metabolomics-count-matrix.csv`: 194 metabolites × 96 samples

### Species-Specific Lipidomics Matrices
- `D-Apul-lipidomics-count-matrix.csv`: 371 lipids × 29 samples (Acropora pulchra)
- `E-Peve-lipidomics-count-matrix.csv`: 371 lipids × 38 samples (Porites evermanni)
- `F-Ptua-lipidomics-count-matrix.csv`: 371 lipids × 29 samples (Pocillopora tuahiniensis)

### Species-Specific Metabolomics Matrices
- `D-Apul-metabolomics-count-matrix.csv`: 194 metabolites × 29 samples (Acropora pulchra)
- `E-Peve-metabolomics-count-matrix.csv`: 194 metabolites × 38 samples (Porites evermanni)
- `F-Ptua-metabolomics-count-matrix.csv`: 194 metabolites × 29 samples (Pocillopora tuahiniensis)

## Format

All matrices follow the same format as CpG count matrices:
- First column (no header): Feature IDs (lipid names or metabolite names)
- Column headers: Sample names in format "ACR-139-TP1", "POR-216-TP2", etc.
- Values: Numeric concentrations
  - Lipidomics: nmol per μg protein
  - Metabolomics: Relative quantification values

## Sample Naming Convention

Samples are named as: `{COLONY}-{TIMEPOINT}`
- Examples: ACR-139-TP1, POR-216-TP2, POC-52-TP3
- Species codes:
  - ACR = Acropora pulchra (D-Apul)
  - POR = Porites evermanni (E-Peve)
  - POC = Pocillopora tuahiniensis (F-Ptua)
- Timepoints: TP1, TP2, TP3, TP4

## Data Sources

- **Lipidomics**: Processed from `M-multi-species/data/lipidomics/lipids.rds`
- **Metabolomics**: Processed from `M-multi-species/data/metabolomics/2025-04-15_Roberts-98.xlsx`

## Generation Script

Generated using: `M-multi-species/scripts/generate_count_matrices.R`

## Usage

These matrices can be used directly in analyses that expect count matrix format, similar to how CpG methylation matrices are used in the pipeline.