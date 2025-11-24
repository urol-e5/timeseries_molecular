# Visualization of BP GO Processes for Rank 35 Components

## Overview

This directory contains visualizations of Biological Process (BP) Gene Ontology terms from the top 100 genes in each of the 35 components from the rank 35 tensor decomposition analysis.

## Files Generated

### Individual Component Visualizations (35 files)
- `rank35_comp1_BP_GO_top20.png` through `rank35_comp35_BP_GO_top20.png`
  - Bar plots showing the top 20 most frequent BP GO terms for each component
  - Each visualization shows the frequency distribution of GO terms

### Summary Visualizations
- **`rank35_ALL_components_BP_GO_summary.png`**: Bar plot of the top 30 BP GO terms across all 35 components
  - Red numbers indicate in how many components each term appears
  - Sorted by total frequency across all components

- **`rank35_BP_GO_heatmap.png`**: Heatmap showing the distribution of the top 30 GO terms across all 35 components
  - Rows: GO BP terms (ordered by overall frequency)
  - Columns: Components (comp1 through comp35)
  - Color intensity: Frequency count (white = 0, dark blue = high)

### Data Tables
- **`rank35_BP_GO_summary_table.csv`**: Summary table with three columns:
  - `term`: GO BP term name
  - `total_count`: Total frequency across all components
  - `n_components`: Number of components containing the term

## Top 10 Most Common GO BP Terms

1. **Innate immune response** (80 occurrences, 27 components)
2. **Proteolysis** (79 occurrences, 25 components)
3. **Positive regulation of gene expression** (78 occurrences, 27 components)
4. **Defense response to Gram-negative bacterium** (76 occurrences, 30 components)
5. **Positive regulation of transcription by RNA polymerase II** (76 occurrences, 21 components)
6. **Adenylate cyclase-activating G protein-coupled receptor signaling pathway** (69 occurrences, 33 components)
7. **Cell surface receptor signaling pathway** (59 occurrences, 29 components)
8. **Gene expression** (59 occurrences, 25 components)
9. **Signal transduction** (59 occurrences, 21 components)
10. **Negative regulation of transcription by RNA polymerase II** (47 occurrences, 18 components)

## Key Findings

The analysis reveals several major biological themes across the rank 35 components:

1. **Immune Response**: Strong representation of innate immunity and defense against bacterial pathogens
2. **Transcriptional Regulation**: Many terms related to gene expression and transcription regulation
3. **Signal Transduction**: Significant presence of cell signaling and receptor-mediated pathways
4. **Structural/Developmental**: Terms related to cytoskeleton, cell adhesion, and basement membrane organization
5. **Proteolysis**: High frequency of protein degradation processes

## Source Data

- **Input Directory**: `/Users/sr320/Documents/GitHub/timeseries_molecular/M-multi-species/output/22-Visualizing-Rank-outs/`
- **Input Files**: 35 files matching pattern `rank_35_comp*_top100_annotation.csv`
- **Analysis Code**: `/Users/sr320/Documents/GitHub/timeseries_molecular/M-multi-species/scripts/23-visualizing-rank35/01-visualize-BP-GO-processes.Rmd`
- **HTML Report**: `/Users/sr320/Documents/GitHub/timeseries_molecular/M-multi-species/scripts/23-visualizing-rank35/01-visualize-BP-GO-processes.html`

## Methods

For each component:
1. Extracted all GO BP terms from the `go_bp` column
2. Counted frequency of each unique term
3. Selected top 20 terms for visualization
4. Created bar plots showing term frequencies

Summary analyses:
1. Combined data from all 35 components
2. Calculated total frequency and component prevalence for each term
3. Generated overall summary visualization and heatmap
4. Exported summary statistics to CSV

## Date Generated

November 2, 2025

