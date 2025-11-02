# Redundancy Analysis: Group IDs Across Rank 05 Components

## Executive Summary

**Rank 05 shows DRAMATICALLY LOWER redundancy (43-70%) compared to Rank 35 (85-96%)**, indicating that lower rank decompositions capture more diverse and specific biological variation, while higher ranks capture core ubiquitous functions.

## Key Findings

### Overall Redundancy Statistics

Looking at the top GO BP terms in Rank 05:

| GO Term | Total Occurrences | Unique group_ids | Redundant Occurrences | Redundancy % |
|---------|-------------------|------------------|----------------------|--------------|
| **Negative regulation of transcription by RNA polymerase II** | 28 | 12 | 16 | **57.1%** |
| **Regulation of transcription by RNA polymerase II** | 26 | 12 | 14 | **53.8%** |
| **Positive regulation of transcription by RNA polymerase II** | 23 | 13 | 10 | **43.5%** |
| **Positive regulation of DNA-templated transcription** | 19 | 10 | 9 | **47.4%** |
| **Positive regulation of gene expression** | 17 | 8 | 9 | **52.9%** |
| **Nervous system development** | 16 | 7 | 9 | **56.3%** |
| **Cell division** | 15 | 7 | 8 | **53.3%** |
| **Cell population proliferation** | 15 | 6 | 9 | **60.0%** |
| **Wnt signaling pathway** | 14 | 6 | 8 | **57.1%** |
| **Kidney development** | 12 | 4 | 8 | **66.7%** |

### What This Means

1. **Much Lower Redundancy**: For most GO terms, 43-70% of appearances are redundant (vs 85-96% in rank 35)

2. **More Unique Genes**: Top GO terms driven by larger, more diverse gene sets:
   - "Positive regulation of transcription by RNA polymerase II": **13 unique group_ids** for 23 occurrences (vs 10 for 91 in rank 35)
   - "Negative regulation of transcription by RNA polymerase II": **12 unique group_ids** for 28 occurrences (vs 9 for 66 in rank 35)
   - Average: **8-13 unique genes per term** (vs 2-11 in rank 35)

3. **Few Ubiquitous Group IDs**: Only 1 group_id appears in ALL components:
   - **OG_04439**: Present in all 5 components (vs 3 group_ids in all 35 components for rank 35)

## Comparison: Rank 05 vs Rank 35

### Redundancy Comparison

| Metric | Rank 05 | Rank 35 | Difference |
|--------|---------|---------|------------|
| **Average Redundancy (top 30 terms)** | 43-70% | 85-96% | **-30 to -40%** |
| **Unique group_ids per term** | 8-13 | 2-11 | **More diverse** |
| **Group_ids in ALL components** | 1 (of 5) | 3 (of 35) | **Fewer ubiquitous** |
| **Max component appearance** | 5/5 (100%) | 35/35 (100%) | **Similar saturation** |

### Biological Interpretation

#### Rank 05 (Lower Rank):
- **Captures specific biological variation**
- **Component-specific gene signatures**
- **More diverse functional modules**
- **Better for:**
  - Identifying distinct biological processes
  - Finding component-specific features
  - Differential expression patterns

#### Rank 35 (Higher Rank):
- **Captures core housekeeping functions**
- **Ubiquitous gene expression patterns**
- **Shared functional modules**
- **Better for:**
  - Identifying conserved mechanisms
  - Finding core regulatory genes
  - Broadly expressed pathways

### Functional Differences

**Rank 05 Top Processes:**
- Transcriptional regulation (diverse regulators)
- Developmental processes (kidney, nervous system)
- Wnt signaling pathway
- Cell division and proliferation

**Rank 35 Top Processes:**
- Innate immune response
- Proteolysis
- Defense responses
- GPCR signaling
- General stress responses

## Specific Examples

### Rank 05: Diverse Gene Sets

**"Positive regulation of transcription by RNA polymerase II"**
- **13 unique group_ids** for 23 occurrences
- Average of 1.8 times per gene
- Represents **diverse transcriptional regulators**

### Rank 35: Same Genes Repeatedly

**"Positive regulation of transcription by RNA polymerase II"**
- **10 unique group_ids** for 91 occurrences
- Average of 9.1 times per gene
- Represents **core ubiquitous regulators**

### Most Ubiquitous Genes in Rank 05

Only **1 group_id appears in all 5 components**:
- **OG_04439**: Cilium movement involved in cell motility

**13 group_ids appear in 4 of 5 components**:
- OG_00887 (6 GO terms)
- OG_01210 (12 GO terms)
- OG_01475 (5 GO terms)
- OG_01489 (18 GO terms)
- OG_02599 (64 GO terms!)
- OG_03061, OG_04129, OG_04170, OG_04335, OG_04737, OG_04999, OG_05110, OG_05208

### OG_02599 - The Most Functionally Diverse

This group_id in rank 05 shows incredible functional diversity:
- **Present in**: 4 of 5 components
- **Associated with**: **64 different GO terms** including:
  - Basement membrane organization
  - Bone morphogenesis
  - Anterior neuropore closure
  - And 61 more distinct processes

## Implications

### For Rank Selection:
1. **Use Rank 05 (or lower ranks)** when you want:
   - Component-specific signatures
   - Diverse biological processes
   - Distinct functional modules
   - Less redundancy

2. **Use Rank 35 (or higher ranks)** when you want:
   - Core conserved functions
   - Ubiquitous gene expression
   - Broadly shared mechanisms
   - Comprehensive coverage

### For Interpretation:
1. **Lower redundancy in rank 05 means** the cross-component frequencies reflect genuinely diverse gene sets performing related functions
   
2. **Higher redundancy in rank 35 means** the cross-component frequencies largely reflect the same small set of core genes

3. **Component distinctiveness**: Rank 05 components are more biologically distinct from each other

### Statistical Consideration:
- In rank 05: "28 occurrences across 5 components" = **12 unique genes** (more meaningful diversity)
- In rank 35: "88 occurrences across 35 components" = **5 unique genes** (highly redundant)

## Data Files Generated

1. **`rank05_GO_redundancy_full_analysis.csv`**: Complete redundancy analysis for all GO terms
2. **`rank05_redundant_group_ids_details.csv`**: Which specific group_ids appear in which components
3. **`rank05_group_id_component_frequency.csv`**: Group_ids ranked by how many components they appear in
4. **`rank05_GO_redundancy_analysis.png`**: Visualization showing unique vs redundant occurrences
5. **`rank05_GO_redundancy_percentage.png`**: Bar chart of redundancy percentages
6. **`rank05_most_frequent_group_ids.png`**: Top group_ids by component frequency

## Interactive HTML Report

Full interactive report: `/Users/sr320/Documents/GitHub/timeseries_molecular/M-multi-species/scripts/24-visualizing-rank05/02-analyze-group-id-redundancy.html`

## Recommendation

**For biological interpretation, rank 05 provides more meaningful insights** with:
- Greater functional diversity
- More component-specific signals
- Less redundancy in feature representation
- More interpretable biological modules

---

**Date**: November 2, 2025

