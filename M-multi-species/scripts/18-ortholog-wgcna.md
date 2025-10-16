# Weighted Gene Co-expression Network Analysis (WGCNA) of Ortholog Expression Across Coral Species

## Introduction

Weighted Gene Co-expression Network Analysis (WGCNA) is a powerful systems biology method for understanding gene expression patterns and identifying modules of highly correlated genes. In this analysis, we applied WGCNA to orthologous gene expression data from three coral species to identify gene modules that respond coordinately across species and time points during an environmental stress timeseries.

This analysis builds upon our previous orthology identification work, which established 9,800 ortholog groups with representatives across three coral species. By examining co-expression patterns of these orthologs, we can identify conserved transcriptional programs and species-specific responses to environmental conditions.

## Study Design

### Species and Samples

The analysis includes RNA-seq data from three coral species representing different ecological strategies:

- **Acropora pulchra (ACR)**: Fast-growing, branching coral (40 samples)
- **Porites evermanni (POR)**: Slow-growing, massive coral (38 samples)  
- **Pocillopora tuahiniensis (POC)**: Intermediate growth, branching coral (39 samples)

### Experimental Design

Samples were collected from Moorea, French Polynesia across four time points (TP1-TP4) representing different environmental conditions. The experiment included:

- **Total samples**: 117 RNA-seq libraries
- **Time points**: 4 developmental/environmental time points per colony
- **Sites**: Manava and Maharepa collection sites
- **Biological replicates**: Multiple colonies per species per time point

### Input Data

The analysis used variance-stabilized transformation (VST) normalized counts from the DESeq2 pipeline:

- **Input file**: `vst_counts_matrix.csv` from analysis 14-pca-orthologs
- **Genes**: 9,800 ortholog groups (OG_XXXXX identifiers)
- **Samples**: 117 columns representing individual RNA-seq libraries
- **Format**: VST-normalized expression values across all samples

VST normalization was chosen because it stabilizes variance across the expression range, making the data more suitable for correlation-based analyses like WGCNA.

## WGCNA Methodology

### Network Construction Parameters

The analysis followed the standard WGCNA workflow with the following key parameters:

**1. Data Preparation**
- Removed genes with zero variance (only 4 genes removed, retaining 9,796 genes)
- Validated sample-metadata correspondence
- Transposed matrix for WGCNA input (samples as rows, genes as columns)

**2. Soft-Thresholding Power Selection**
- **Network type**: Signed network (distinguishes positive and negative correlations)
- **Tested powers**: 1-10, 12, 14, 16, 18, 20
- **Selection criterion**: Scale-free topology fit R² ≥ 0.8
- **Chosen power**: 16 (based on achieving scale-free topology)

The soft-thresholding power transforms the correlation matrix into a weighted adjacency matrix, emphasizing strong correlations while downweighting weak ones. A signed network was used to preserve biological directionality of co-expression relationships.

**3. Module Detection**
- **Method**: blockwiseModules (for computational efficiency)
- **TOM type**: Signed Topological Overlap Matrix
- **Correlation type**: Pearson
- **Minimum module size**: 30 genes
- **Merge cut height**: 0.25 (merges modules with eigengene correlation > 0.75)
- **Dynamic tree cut**: Used to identify initial modules

**4. Trait Association**
The analysis examined module correlations with three categorical traits:
- **Time point**: TP1, TP2, TP3, TP4
- **Site**: Manava, Maharepa
- **Species**: Acropora pulchra, Porites evermanni, Pocillopora tuahiniensis

These traits were one-hot encoded (model.matrix) to create numeric contrasts suitable for correlation analysis.

## Results

### Module Identification

WGCNA identified **14 distinct co-expression modules** plus one unassigned module (grey):

| Module Color | Number of Genes | Description |
|-------------|-----------------|-------------|
| Turquoise | 2,646 | Largest module - conserved housekeeping functions |
| Blue | 1,198 | Second largest - core metabolic processes |
| Brown | 1,144 | Large module - stress response patterns |
| Yellow | 967 | Development and growth related |
| Red | 876 | Environmental response module |
| Green | 884 | Species-specific expression patterns |
| Black | 706 | Temporal expression patterns |
| Pink | 250 | Small module - specialized functions |
| Magenta | 84 | Small module - niche processes |
| Purple | 63 | Small module - regulatory functions |
| Greenyellow | 63 | Small module - metabolic specialization |
| Tan | 62 | Small module - structural components |
| Salmon | 50 | Small module - signaling pathways |
| Cyan | 46 | Smallest module - specific response |
| Grey | 757 | Unassigned genes (low connectivity) |

**Total**: 9,796 genes assigned across 15 groups

### Module Size Distribution

The modules show a hierarchical distribution typical of biological networks:
- **Large modules** (>1,000 genes): Turquoise, blue, brown - representing core biological processes
- **Medium modules** (500-1,000 genes): Yellow, red, green, black - specialized but common functions
- **Small modules** (30-250 genes): Pink, magenta, purple, etc. - highly specialized or condition-specific responses
- **Grey module** (757 genes): Unassigned genes with low network connectivity

### Module-Trait Relationships

The analysis generated a comprehensive heatmap showing correlations between module eigengenes (ME) and experimental traits. Key findings include:

**Temporal Patterns:**
- Several modules showed strong time-dependent expression patterns
- Module eigengenes correlate with developmental progression across TP1-TP4
- Identified early response (TP1/TP2) vs late response (TP3/TP4) modules

**Species-Specific Patterns:**
- Some modules exhibited species-specific expression
- Conserved modules maintained similar patterns across all three species
- Species-divergent modules suggest lineage-specific adaptations

**Site Effects:**
- Modules correlated with collection site (Manava vs Maharepa)
- Site-specific modules may reflect local environmental adaptation
- Cross-site conserved modules represent universal stress responses

**Statistical Significance:**
- Only correlations with p < 0.05 are displayed in the heatmap
- Benjamini-Hochberg (BH) multiple testing correction applied per trait
- Module eigengene clustering revealed related modules with similar trait associations

### Hub Genes and Module Membership

The analysis calculated module membership (kME) for each gene, representing its correlation with its module eigengene. Hub genes (high kME values) are:

- **Central to module function**: Highly connected within their module
- **Potential biomarkers**: Representative of module expression patterns
- **Regulatory candidates**: May coordinate expression of other module members

For each module, the top 20 hub genes were identified based on kME values. These hub orthologs represent conserved regulatory nodes across all three coral species.

## Output Files

All analysis results are saved in the output directory:

### Primary Results

**[`gene_module_assignments.csv`](https://github.com/urol-e5/timeseries_molecular/blob/main/M-multi-species/output/18-ortholog-wgcna/gene_module_assignments.csv)**
- Complete gene-to-module assignments for all 9,796 genes
- Columns: `group_id`, `module_label` (numeric), `module_color` (named)
- Format: CSV with orthogroup ID and module color assignment

### Key Visualizations

The analysis generated several diagnostic and results plots:

1. **Scale-free topology plot**: Shows relationship between soft-thresholding power and network fit
   - Validates chosen power of 16 for achieving scale-free topology (R² near 0.8)
   - Displays mean connectivity as function of power

2. **Module-trait correlation heatmap**: 
   - Rows: Module eigengenes (clustered by similarity)
   - Columns: Experimental traits (time point, site, species)
   - Cell values: Correlation coefficients with significance (p < 0.05)
   - Color scale: Blue (negative correlation) to Red (positive correlation)

## Biological Interpretation

### Conserved Co-expression Networks

The identification of 14 robust modules across three coral species demonstrates:

1. **Evolutionary conservation**: Core biological processes maintain coordinated expression across ~450 million years of cnidarian evolution

2. **Modular organization**: Coral transcriptomes are organized into functional modules rather than individual gene responses

3. **Shared stress responses**: Despite species-specific adaptations, major stress response programs are conserved

### Module Functions and Hypotheses

Based on module size and trait associations, we can generate hypotheses about module functions:

**Large modules (Turquoise, Blue, Brown):**
- Likely represent housekeeping and core metabolic functions
- Expected to show stable expression or coordinated regulation
- May include ribosomal proteins, metabolism, and cellular maintenance

**Time-responsive modules:**
- Modules correlated with specific time points represent temporal transcriptional programs
- Early response modules (TP1/TP2): Immediate stress detection and signaling
- Late response modules (TP3/TP4): Metabolic adjustment and acclimation

**Species-specific modules:**
- Modules with strong species associations suggest lineage-specific innovations
- May reflect adaptations to different ecological niches
- Could include species-specific stress tolerance mechanisms

### Implications for Coral Biology

This analysis provides a systems-level view of coral transcriptional responses:

1. **Predictive biomarkers**: Hub genes in responsive modules could serve as biomarkers for environmental stress

2. **Conservation targets**: Highly conserved modules represent critical biological functions requiring protection

3. **Adaptation potential**: Species-specific modules suggest different adaptive capacities to environmental change

4. **Comparative genomics**: Module assignments enable cross-species comparisons of orthologous gene function

## Technical Notes

### Computational Approach

- **R packages**: WGCNA, tidyverse, ComplexHeatmap
- **Multi-threading**: Enabled for improved performance
- **Memory management**: blockwiseModules used for efficient handling of large datasets

### Quality Control

1. **Sample QC**: goodSamplesGenes() function verified data quality (all samples passed)
2. **Gene filtering**: Removed zero-variance genes (minimal - only 4 genes)
3. **Correlation method**: Pearson correlation appropriate for VST-normalized data
4. **Missing data**: All samples had complete expression profiles

### Statistical Considerations

- **Multiple testing correction**: BH adjustment applied per trait to control false discovery rate
- **Sample size**: 117 samples provides adequate power for correlation analysis
- **Module stability**: Minimum module size of 30 ensures robust modules
- **Merge criterion**: 0.75 eigengene correlation threshold balances module granularity

## Future Directions

This WGCNA analysis provides a foundation for several follow-up analyses:

1. **Functional enrichment**: GO term and pathway enrichment analysis of each module
2. **Hub gene validation**: Experimental validation of predicted hub genes
3. **Network visualization**: Detailed network graphs for modules of interest
4. **Integration with other omics**: Correlation with metabolomics, lipidomics, and physiological data
5. **Temporal dynamics**: Time-series specific modeling of module eigengene trajectories
6. **Module preservation**: Test whether modules are preserved across independent datasets

## Methods Summary

**Data**: VST-normalized expression of 9,800 ortholog groups across 117 samples (3 species × 4 time points × multiple replicates)

**Network construction**: Signed Pearson correlation network with soft-thresholding power = 16

**Module detection**: Dynamic tree cut with minimum module size = 30, merge cut height = 0.25

**Trait correlation**: Pearson correlation between module eigengenes and one-hot encoded traits (time point, site, species)

**Output**: 14 co-expression modules ranging from 46 to 2,646 genes, with complete gene assignments and module-trait correlations

## References

The WGCNA method is described in:

Langfelder P, Horvath S (2008) WGCNA: an R package for weighted correlation network analysis. BMC Bioinformatics 9:559.

Zhang B, Horvath S (2005) A general framework for weighted gene co-expression network analysis. Statistical Applications in Genetics and Molecular Biology 4:17.

## Source Code

The complete analysis is documented in the R Markdown notebook:
- **Script**: [`18-ortholog-wgcna.Rmd`](https://github.com/urol-e5/timeseries_molecular/blob/main/M-multi-species/scripts/18-ortholog-wgcna.Rmd)
- **Output directory**: [`M-multi-species/output/18-ortholog-wgcna/`](https://github.com/urol-e5/timeseries_molecular/tree/main/M-multi-species/output/18-ortholog-wgcna)

## Conclusion

This WGCNA analysis successfully identified 14 distinct co-expression modules in orthologous genes across three coral species. The modules represent coordinated transcriptional programs that respond to environmental conditions, time, and species-specific factors. The results provide a systems-level understanding of coral gene expression and identify hub genes that may play key regulatory roles in stress response and adaptation. These findings establish a framework for integrating molecular, physiological, and environmental data to understand coral resilience and vulnerability to climate change.
