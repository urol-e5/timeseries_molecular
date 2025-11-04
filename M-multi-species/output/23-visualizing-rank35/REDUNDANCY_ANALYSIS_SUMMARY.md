# Redundancy Analysis: Group IDs Across Rank 35 Components

## Executive Summary

**The vast majority (85-95%) of GO term occurrences across components are due to the SAME group_id appearing in multiple components**, not unique genes with the same functional annotation.

## Key Findings

### Overall Redundancy Statistics

Looking at the top 30 most frequent GO BP terms:

| GO Term | Total Occurrences | Unique group_ids | Redundant Occurrences | Redundancy % |
|---------|-------------------|------------------|----------------------|--------------|
| **Innate immune response** | 88 | 5 | 83 | **94.3%** |
| **Proteolysis** | 94 | 11 | 83 | **88.3%** |
| **Defense response to Gram-negative bacterium** | 81 | 5 | 76 | **93.8%** |
| **Positive regulation of gene expression** | 86 | 10 | 76 | **88.4%** |
| **Positive regulation of transcription by RNA polymerase II** | 91 | 10 | 81 | **89.0%** |
| **Adenylate cyclase-activating GPCR signaling** | 71 | 4 | 67 | **94.4%** |
| **Cell surface receptor signaling pathway** | 66 | 3 | 63 | **95.5%** |
| **Signal transduction** | 79 | 6 | 73 | **92.4%** |
| **Gene expression** | 67 | 6 | 61 | **91.0%** |
| **Positive regulation of ROS biosynthesis** | 58 | 2 | 56 | **96.6%** |

### What This Means

1. **High Redundancy**: For most GO terms, 85-96% of their appearances across components come from the same small set of group_ids (orthogroups) appearing repeatedly.

2. **Few Unique Genes**: The top GO terms are driven by very few unique genes:
   - "Innate immune response": Only **5 unique group_ids** account for 88 occurrences
   - "Cell surface receptor signaling": Only **3 unique group_ids** account for 66 occurrences
   - "Positive regulation of ROS biosynthesis": Only **2 unique group_ids** account for 58 occurrences

3. **Highly Prevalent Group IDs**: Some group_ids appear in nearly ALL components:
   - **OG_00264**: Present in all 35 components
   - **OG_01620**: Present in all 35 components  
   - **OG_02537**: Present in all 35 components (with 22 different GO terms!)
   - **OG_01636**: Present in 34 components
   - **OG_01919**: Present in 31 components
   - **OG_02927**: Present in 29 components

## Specific Examples

### OG_02537 - The Most Ubiquitous Group
- **Present in**: All 35 components
- **Associated GO terms**: 22 different biological processes including:
  - Adenylate cyclase-activating G protein-coupled receptor signaling pathway
  - Cell adhesion
  - Cell surface receptor signaling pathway
  - Apoptotic cell clearance
  - Signal transduction
  - And 17 more...

### OG_01645 - Structural Gene Across 26 Components
- **Present in**: 26 components (comp1, comp2, comp3, comp4, comp6, comp7, comp8, comp9, comp11, comp13, comp14, comp15, comp16, comp18, comp20, comp21, comp22, comp24, comp25, comp26, comp27, comp28, comp30, comp31, comp32, comp34)
- **Associated with**:
  - Actin cytoskeleton organization
  - Basement membrane organization
  - Cell division
  - Cell adhesion

### OG_02927 - Signaling Gene Across 29 Components
- **Present in**: 29 components
- **Associated with**:
  - Adenylate cyclase-activating GPCR signaling pathway
  - Cell surface receptor signaling pathway
  - Circadian behavior
  - G protein-coupled receptor signaling pathway
  - Signal transduction

## Implications

### For Interpretation:
1. **The high cross-component frequencies of GO terms are NOT due to diverse sets of genes** performing similar functions independently in each component.

2. **Instead, a small core set of genes (group_ids)** is consistently selected as "top genes" across most or all of the 35 components in the tensor decomposition.

3. **This suggests these genes have high loadings/weights** in the tensor decomposition and are core features driving variation across multiple components.

### Biological Interpretation:
- The repeated genes likely represent:
  - **Core housekeeping functions** (e.g., proteolysis, gene expression)
  - **Fundamental signaling pathways** (e.g., GPCR signaling)
  - **Conserved stress responses** (e.g., innate immunity, ROS production)
  
- These genes may be:
  - Highly expressed across all samples/conditions
  - Highly variable across the dataset
  - Central hub genes in regulatory networks

### Statistical Consideration:
- When reporting "80 occurrences across 27 components" for innate immune response, this is really **5 genes appearing 17.6 times on average**, not 80 different gene annotations.
  
- The effective diversity of biological functions is much lower than the raw occurrence counts suggest.

## Data Files Generated

1. **`rank35_GO_redundancy_full_analysis.csv`**: Complete redundancy analysis for all GO terms
2. **`rank35_redundant_group_ids_details.csv`**: Which specific group_ids appear in which components
3. **`rank35_group_id_component_frequency.csv`**: Group_ids ranked by how many components they appear in
4. **`rank35_GO_redundancy_analysis.png`**: Visualization showing unique vs redundant occurrences
5. **`rank35_GO_redundancy_percentage.png`**: Bar chart of redundancy percentages
6. **`rank35_most_frequent_group_ids.png`**: Top group_ids by component frequency

## Interactive HTML Report

Full interactive report: `/Users/sr320/Documents/GitHub/timeseries_molecular/M-multi-species/scripts/23-visualizing-rank35/02-analyze-group-id-redundancy.html`

---

**Date**: November 2, 2025

