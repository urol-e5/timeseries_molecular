# Scripts for Visualizing Rank 35 BP GO Processes

## Overview

This directory contains the analysis script for visualizing Biological Process Gene Ontology terms from the rank 35 tensor decomposition components.

## Files

### Analysis Scripts

- **`01-visualize-BP-GO-processes.Rmd`**: R Markdown document that performs the complete analysis
  - Reads all 35 rank_35_comp*_top100_annotation.csv files
  - Parses GO BP terms from the annotation data
  - Creates individual visualizations for each component
  - Generates summary visualizations and statistics
  - Produces an interactive HTML report

- **`01-visualize-BP-GO-processes.html`**: HTML report with all visualizations and analysis

## Running the Analysis

To regenerate the analysis, run:

```bash
cd /Users/sr320/Documents/GitHub/timeseries_molecular/M-multi-species/scripts/23-visualizing-rank35
Rscript -e "rmarkdown::render('01-visualize-BP-GO-processes.Rmd')"
```

## Required R Packages

- `tidyverse`: Data manipulation and visualization
- `ggplot2`: Creating plots
- `RColorBrewer`: Color palettes
- `scales`: Formatting scales
- `rmarkdown`: Rendering R Markdown documents

## Input

- **Directory**: `/Users/sr320/Documents/GitHub/timeseries_molecular/M-multi-species/output/22-Visualizing-Rank-outs/`
- **Files**: 35 CSV files matching pattern `rank_35_comp*_top100_annotation.csv`

## Output

All outputs are saved to: `/Users/sr320/Documents/GitHub/timeseries_molecular/M-multi-species/output/23-visualizing-rank35/`

### Generated Files

1. **Individual Component Plots** (35 PNG files):
   - `rank35_comp1_BP_GO_top20.png` through `rank35_comp35_BP_GO_top20.png`
   - Dimensions: 14" x 10", 300 DPI
   - Shows top 20 BP GO terms for each component

2. **Summary Visualizations**:
   - `rank35_ALL_components_BP_GO_summary.png`: Overall top 30 terms across all components
   - `rank35_BP_GO_heatmap.png`: Heatmap of term distribution across components
   - Dimensions: 14-16" x 10-12", 300 DPI

3. **Data Tables**:
   - `rank35_BP_GO_summary_table.csv`: Summary statistics table

4. **Documentation**:
   - `README.md`: Description of outputs and findings

## Analysis Workflow

1. **File Discovery**: Identifies all rank_35 annotation files
2. **GO Term Parsing**: Extracts BP GO terms from the `go_bp` column
3. **Frequency Analysis**: Counts occurrences of each term per component
4. **Individual Visualizations**: Creates bar plots for each component's top 20 terms
5. **Cross-Component Analysis**: Aggregates data across all components
6. **Summary Visualizations**: Creates overall bar plot and heatmap
7. **Report Generation**: Compiles all results into interactive HTML

## Key Functions

- `parse_go_bp()`: Parses GO BP terms from semicolon-separated strings
- `process_file()`: Processes a single annotation file and returns term counts

## Date Created

November 2, 2025

