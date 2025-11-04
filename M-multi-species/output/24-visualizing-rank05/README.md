# Visualization of BP GO Processes for Rank 05 Components

## Overview

This directory contains visualizations of Biological Process (BP) Gene Ontology terms from the top 100 genes in each of the 5 components from the rank 05 tensor decomposition analysis.

## Files Generated

### Individual Component Visualizations (5 files)
- `rank05_comp1_BP_GO_top20.png` through `rank05_comp5_BP_GO_top20.png`
  - Bar plots showing the top 20 most frequent BP GO terms for each component
  - Each visualization shows the frequency distribution of GO terms

### Summary Visualizations
- **`rank05_ALL_components_BP_GO_summary.png`**: Bar plot of the top 30 BP GO terms across all 5 components
  - Red numbers indicate in how many components each term appears
  - Sorted by total frequency across all components

- **`rank05_BP_GO_heatmap.png`**: Heatmap showing the distribution of the top 30 GO terms across all 5 components
  - Rows: GO BP terms (ordered by overall frequency)
  - Columns: Components (comp1 through comp5)
  - Color intensity: Frequency count (white = 0, dark blue = high)

### Redundancy Analysis Visualizations
- **`rank05_GO_redundancy_analysis.png`**: Stacked bar plot showing unique vs redundant occurrences
- **`rank05_GO_redundancy_percentage.png`**: Redundancy percentages for each GO term
- **`rank05_most_frequent_group_ids.png`**: Top group_ids by number of components

### Data Tables
- **`rank05_BP_GO_summary_table.csv`**: Summary table with three columns:
  - `term`: GO BP term name
  - `total_count`: Total frequency across all components
  - `n_components`: Number of components containing the term

- **`rank05_GO_redundancy_full_analysis.csv`**: Complete redundancy analysis
- **`rank05_redundant_group_ids_details.csv`**: Which group_ids appear in multiple components
- **`rank05_group_id_component_frequency.csv`**: Group_ids ranked by component frequency

## Top 10 Most Common GO BP Terms

1. **Negative regulation of transcription by RNA polymerase II** (28 occurrences, 5 components)
2. **Regulation of transcription by RNA polymerase II** (26 occurrences, 5 components)
3. **Positive regulation of transcription by RNA polymerase II** (23 occurrences, 5 components)
4. **Positive regulation of DNA-templated transcription** (17 occurrences, 4 components)
5. **Positive regulation of gene expression** (17 occurrences, 5 components)
6. **Wnt signaling pathway** (14 occurrences, 5 components)
7. **Cell population proliferation** (14 occurrences, 4 components)
8. **Nervous system development** (14 occurrences, 4 components)
9. **Cell division** (13 occurrences, 4 components)
10. **Kidney development** (12 occurrences, 4 components)

## Key Findings

### Biological Themes

The rank 05 analysis reveals major biological themes:

1. **Transcriptional Regulation**: Dominant theme with regulation of transcription (both positive and negative)
2. **Development**: Strong representation of developmental processes (nervous system, kidney, Wnt signaling)
3. **Cell Division & Proliferation**: Cell cycle and population growth
4. **Differentiation**: Tissue-specific development processes

### Redundancy Analysis: CRITICAL FINDING

**Unlike rank 35, rank 05 shows MUCH LOWER redundancy (43-70% vs 85-96%)**

| GO Term | Total Occurrences | Unique group_ids | Redundancy % |
|---------|-------------------|------------------|--------------|
| **Negative regulation of transcription by RNA polymerase II** | 28 | 12 | **57.1%** |
| **Regulation of transcription by RNA polymerase II** | 26 | 12 | **53.8%** |
| **Positive regulation of transcription by RNA polymerase II** | 23 | 13 | **43.5%** |
| **Kidney development** | 12 | 4 | 66.7% |
| **Cellular response to leukemia inhibitory factor** | 10 | 3 | 70.0% |

### What This Means: Rank 05 vs Rank 35 Comparison

**Rank 05 characteristics:**
- **More diverse gene sets**: 8-13 unique group_ids per top term (vs 2-11 in rank 35)
- **Lower redundancy**: Average 43-70% (vs 85-96% in rank 35)
- **Fewer ubiquitous genes**: Only 1 group_id (OG_04439) appears in all 5 components
- **More component-specific**: Group_ids tend to appear in 4-5 components max

**Biological Interpretation:**
- Rank 05 captures **more specific biological variation** with diverse gene sets
- Rank 35 captures **core housekeeping/ubiquitous functions** with the same genes repeatedly
- Lower rank decomposition = more distinct component signatures
- Higher rank decomposition = more shared/redundant features

### Most Frequent Group IDs

Only 1 group_id appears in all components:
- **OG_04439**: Present in all 5 components (cilium movement involved in cell motility)

Several group_ids appear in 4 components:
- **OG_00887**: 4 components (6 GO terms)
- **OG_01210**: 4 components (12 GO terms)
- **OG_01475**: 4 components (5 GO terms)
- **OG_01489**: 4 components (18 GO terms)
- **OG_02599**: 4 components (64 GO terms!)
- And 8 more...

## Source Data

- **Input Directory**: `/Users/sr320/Documents/GitHub/timeseries_molecular/M-multi-species/output/22-Visualizing-Rank-outs/`
- **Input Files**: 5 files matching pattern `rank_05_comp*_top100_annotation.csv`
- **Analysis Code**: `/Users/sr320/Documents/GitHub/timeseries_molecular/M-multi-species/scripts/24-visualizing-rank05/`
  - `01-visualize-BP-GO-processes.Rmd`
  - `02-analyze-group-id-redundancy.Rmd`
- **HTML Reports**: 
  - `/Users/sr320/Documents/GitHub/timeseries_molecular/M-multi-species/scripts/24-visualizing-rank05/01-visualize-BP-GO-processes.html`
  - `/Users/sr320/Documents/GitHub/timeseries_molecular/M-multi-species/scripts/24-visualizing-rank05/02-analyze-group-id-redundancy.html`

## Methods

For each component:
1. Extracted all GO BP terms from the `go_bp` column
2. Counted frequency of each unique term
3. Selected top 20 terms for visualization
4. Created bar plots showing term frequencies

Summary analyses:
1. Combined data from all 5 components
2. Calculated total frequency and component prevalence for each term
3. Generated overall summary visualization and heatmap
4. Exported summary statistics to CSV

Redundancy analysis:
1. Tracked group_id occurrences across components for each GO term
2. Calculated redundancy as: (total occurrences - unique group_ids) / total occurrences
3. Identified which specific group_ids appear in multiple components
4. Compared results with rank 35 analysis

## Date Generated

November 2, 2025

