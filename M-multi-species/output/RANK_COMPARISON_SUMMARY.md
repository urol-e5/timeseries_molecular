# Rank 05 vs Rank 35: Comprehensive Comparison

**Analysis Date**: November 2, 2025

## Executive Summary

A comprehensive analysis of GO BP term distributions and gene redundancy across tensor decomposition ranks reveals **fundamentally different patterns** between rank 05 and rank 35, with critical implications for biological interpretation.

## Key Finding: Redundancy Dramatically Differs by Rank

| Metric | Rank 05 (5 components) | Rank 35 (35 components) | Difference |
|--------|------------------------|-------------------------|------------|
| **Average Redundancy** | 43-70% | 85-96% | **-30 to -40%** |
| **Unique genes per GO term** | 8-13 | 2-11 | **More diverse** |
| **Group_ids in ALL components** | 1 | 3 | **Fewer ubiquitous** |
| **Biological specificity** | High | Low | **More distinct** |

## What This Means

### Rank 05: Specific Biological Variation
- **Lower redundancy** → More unique genes per function
- **Component-specific signatures** → Distinct biological processes
- **Diverse gene sets** → 8-13 genes per top GO term
- **Better for**: Identifying differential processes, component-specific features

### Rank 35: Core Ubiquitous Functions
- **Higher redundancy** → Same genes appear repeatedly
- **Shared signatures** → Core housekeeping functions
- **Concentrated gene sets** → 2-11 genes per top GO term
- **Better for**: Identifying conserved mechanisms, hub genes

## Biological Process Themes

### Rank 05 Top Processes
1. Negative regulation of transcription by RNA polymerase II (28 occurrences, 12 unique genes)
2. Regulation of transcription by RNA polymerase II (26 occurrences, 12 unique genes)
3. Positive regulation of transcription by RNA polymerase II (23 occurrences, 13 unique genes)
4. Positive regulation of gene expression (17 occurrences, 8 unique genes)
5. Wnt signaling pathway (14 occurrences, 6 unique genes)

**Theme**: Developmental regulation and transcriptional control

### Rank 35 Top Processes
1. Innate immune response (88 occurrences, 5 unique genes)
2. Proteolysis (94 occurrences, 11 unique genes)
3. Defense response to Gram-negative bacterium (81 occurrences, 5 unique genes)
4. Positive regulation of gene expression (86 occurrences, 10 unique genes)
5. Signal transduction (79 occurrences, 6 unique genes)

**Theme**: Stress responses and core cellular functions

## Example: Same GO Term, Different Patterns

**"Positive regulation of transcription by RNA polymerase II"**

| Metric | Rank 05 | Rank 35 |
|--------|---------|---------|
| Total occurrences | 23 | 91 |
| Unique group_ids | **13** | 10 |
| Redundancy % | **43.5%** | 89.0% |
| Avg times per gene | 1.8 | 9.1 |

→ **Rank 05** shows **13 diverse transcriptional regulators** each appearing ~2 times
→ **Rank 35** shows **10 core regulators** each appearing ~9 times

## Ubiquitous Genes

### Rank 05: Only 1 gene in ALL components
- **OG_04439**: Cilium movement (5/5 components)

### Rank 35: 3 genes in ALL components
- **OG_00264**: Flagellated sperm motility (35/35 components)
- **OG_01620**: Retinoid/retinol metabolism (35/35 components)
- **OG_02537**: GPCR signaling, 22 GO terms (35/35 components)

## Recommendations

### For Biological Discovery
**→ Prefer Rank 05 (or similarly low ranks)**
- More interpretable biological modules
- Component-specific insights
- Less redundancy in features
- Clearer functional separation

### For Core Gene Identification
**→ Use Rank 35 (or similarly high ranks)**
- Identifies ubiquitous hub genes
- Finds conserved mechanisms
- Comprehensive pathway coverage
- Stable across conditions

### For Method Validation
**→ Compare multiple ranks**
- Check consistency of biological themes
- Validate component assignments
- Assess redundancy patterns
- Ensure biological interpretability

## Statistical Interpretation

### Cross-Component GO Term Frequencies

When you see "28 occurrences across 5 components":
- **In Rank 05**: Likely represents ~12 unique genes (meaningful diversity)
- **In Rank 35**: Likely represents ~5 unique genes (high redundancy)

→ **The same frequency count means very different things depending on rank**

### Effective Biological Diversity

| Rank | Apparent Diversity | Actual Diversity (unique genes) | Effective Ratio |
|------|-------------------|--------------------------------|-----------------|
| Rank 05 | High | High | **1:0.5** (2 occurrences per gene) |
| Rank 35 | High | Low | **1:0.1** (10 occurrences per gene) |

## Files and Documentation

### Rank 05 Analysis
- **Output**: `/Users/sr320/Documents/GitHub/timeseries_molecular/M-multi-species/output/24-visualizing-rank05/`
- **Code**: `/Users/sr320/Documents/GitHub/timeseries_molecular/M-multi-species/scripts/24-visualizing-rank05/`
- **Reports**: 
  - `01-visualize-BP-GO-processes.html`
  - `02-analyze-group-id-redundancy.html`

### Rank 35 Analysis
- **Output**: `/Users/sr320/Documents/GitHub/timeseries_molecular/M-multi-species/output/23-visualizing-rank35/`
- **Code**: `/Users/sr320/Documents/GitHub/timeseries_molecular/M-multi-species/scripts/23-visualizing-rank35/`
- **Reports**:
  - `01-visualize-BP-GO-processes.html`
  - `02-analyze-group-id-redundancy.html`

## Visualization Files

Both analyses include:
- Individual component GO term bar plots
- Cross-component summary plots
- Heatmaps of term distributions
- Redundancy analysis plots
- Group_id frequency plots

Total visualizations generated: **43 PNG files** (38 for rank 35 + 5 for rank 05) + **6 summary plots**

## Conclusion

**The choice of rank fundamentally affects biological interpretation:**

- **Lower ranks** → Component-specific, diverse gene sets, interpretable modules
- **Higher ranks** → Ubiquitous genes, core functions, redundant features

**For most biological questions, rank 05 provides more meaningful insights** due to its:
1. Greater functional diversity
2. Lower redundancy
3. More distinct component signatures
4. More interpretable biological processes

**However, rank 35 is valuable for:**
1. Identifying core conserved genes
2. Finding ubiquitous regulatory mechanisms
3. Comprehensive functional coverage

---

**Generated**: November 2, 2025
**Analyst**: Automated GO Term and Redundancy Analysis Pipeline
