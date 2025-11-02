# Scripts for Visualizing Rank 05 BP GO Processes

## Overview

This directory contains analysis scripts for visualizing Biological Process Gene Ontology terms from the rank 05 tensor decomposition components, with comprehensive redundancy analysis comparing against rank 35.

## Files

### Analysis Scripts

- **`01-visualize-BP-GO-processes.Rmd`**: R Markdown document for GO term visualization
  - Reads all 5 rank_05_comp*_top100_annotation.csv files
  - Parses GO BP terms from the annotation data
  - Creates individual visualizations for each component
  - Generates summary visualizations and statistics
  - Produces an interactive HTML report

- **`02-analyze-group-id-redundancy.Rmd`**: R Markdown document for redundancy analysis
  - Analyzes how often the same group_id appears across components
  - Calculates redundancy percentages for each GO term
  - Identifies ubiquitous vs component-specific genes
  - Compares patterns with rank 35 results
  - Produces comprehensive redundancy report

- **`01-visualize-BP-GO-processes.html`**: HTML report with all GO visualizations
- **`02-analyze-group-id-redundancy.html`**: HTML report with redundancy analysis

## Running the Analysis

To regenerate the analysis, run:

```bash
cd /Users/sr320/Documents/GitHub/timeseries_molecular/M-multi-species/scripts/24-visualizing-rank05

# Run GO term visualization
Rscript -e "rmarkdown::render('01-visualize-BP-GO-processes.Rmd')"

# Run redundancy analysis
Rscript -e "rmarkdown::render('02-analyze-group-id-redundancy.Rmd')"
```

## Required R Packages

- `tidyverse`: Data manipulation and visualization
- `ggplot2`: Creating plots
- `RColorBrewer`: Color palettes
- `scales`: Formatting scales
- `rmarkdown`: Rendering R Markdown documents
- `knitr`: Report generation

## Input

- **Directory**: `/Users/sr320/Documents/GitHub/timeseries_molecular/M-multi-species/output/22-Visualizing-Rank-outs/`
- **Files**: 5 CSV files matching pattern `rank_05_comp*_top100_annotation.csv`

## Output

All outputs are saved to: `/Users/sr320/Documents/GitHub/timeseries_molecular/M-multi-species/output/24-visualizing-rank05/`

### Generated Files

1. **Individual Component Plots** (5 PNG files):
   - `rank05_comp1_BP_GO_top20.png` through `rank05_comp5_BP_GO_top20.png`
   - Dimensions: 14" x 10", 300 DPI
   - Shows top 20 BP GO terms for each component

2. **Summary Visualizations**:
   - `rank05_ALL_components_BP_GO_summary.png`: Overall top 30 terms across all components
   - `rank05_BP_GO_heatmap.png`: Heatmap of term distribution across components
   - Dimensions: 12-14" x 10-12", 300 DPI

3. **Redundancy Visualizations**:
   - `rank05_GO_redundancy_analysis.png`: Stacked bars showing unique vs redundant
   - `rank05_GO_redundancy_percentage.png`: Redundancy percentages
   - `rank05_most_frequent_group_ids.png`: Group_ids by component frequency

4. **Data Tables**:
   - `rank05_BP_GO_summary_table.csv`: Summary statistics
   - `rank05_GO_redundancy_full_analysis.csv`: Complete redundancy data
   - `rank05_redundant_group_ids_details.csv`: Group_ids in multiple components
   - `rank05_group_id_component_frequency.csv`: Group_id frequency table

5. **Documentation**:
   - `README.md`: Description of outputs and findings
   - `REDUNDANCY_ANALYSIS_SUMMARY.md`: Detailed redundancy analysis with rank 05 vs 35 comparison

## Analysis Workflow

### Script 01: GO Term Visualization
1. **File Discovery**: Identifies all rank_05 annotation files
2. **GO Term Parsing**: Extracts BP GO terms from the `go_bp` column
3. **Frequency Analysis**: Counts occurrences of each term per component
4. **Individual Visualizations**: Creates bar plots for each component's top 20 terms
5. **Cross-Component Analysis**: Aggregates data across all components
6. **Summary Visualizations**: Creates overall bar plot and heatmap
7. **Report Generation**: Compiles all results into interactive HTML

### Script 02: Redundancy Analysis
1. **Data Integration**: Loads all components with group_id tracking
2. **GO Term Expansion**: Creates long-format data linking group_ids to GO terms
3. **Redundancy Calculation**: Determines how often same group_id appears across components
4. **Component Frequency**: Identifies ubiquitous vs specific group_ids
5. **Comparison Analysis**: Contrasts patterns with rank 35
6. **Visualization**: Creates redundancy-focused plots
7. **Export**: Saves detailed results for further analysis

## Key Functions

### Script 01
- `parse_go_bp()`: Parses GO BP terms from semicolon-separated strings
- `process_file()`: Processes a single annotation file and returns term counts

### Script 02
- `expand_go_bp_terms()`: Expands GO terms while maintaining group_id tracking
- Creates redundancy metrics for each GO term

## Key Findings

### Major Discovery: Rank 05 vs Rank 35

**Rank 05 shows 30-40% LOWER redundancy than Rank 35:**
- Rank 05: 43-70% redundancy (more diverse gene sets)
- Rank 35: 85-96% redundancy (same genes repeatedly)

This indicates:
- **Lower ranks** (like 05) capture **specific biological variation** with diverse genes
- **Higher ranks** (like 35) capture **core ubiquitous functions** with repeated genes
- **Component distinctiveness** is higher in lower rank decompositions

### Biological Interpretation

**Use Rank 05 for:**
- Component-specific biological signatures
- Diverse functional modules
- Differential expression patterns
- Interpretable biological processes

**Use Rank 35 for:**
- Core conserved mechanisms
- Ubiquitous gene expression
- Comprehensive pathway coverage
- Identifying hub genes

## Date Created

November 2, 2025

