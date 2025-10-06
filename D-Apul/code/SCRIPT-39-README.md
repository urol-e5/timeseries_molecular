# Gene-Body Methylation and Expression Modeling Analysis

## Overview

This document describes the newly created script `39-Apul-genebody-meth-expr-modeling.Rmd` that models the relationship between DNA methylation in distinct genic regions and gene expression.

## Script Location

- **Code**: `D-Apul/code/39-Apul-genebody-meth-expr-modeling.Rmd`
- **Output Directory**: `D-Apul/output/39-Apul-genebody-meth-expr-modeling/`

## Purpose

The script addresses the issue requirement to:
> Develop R code to model how DNA methylation in distinct genic regions (promoter, 5′ end, exons, introns, and 3′ end) influences gene expression.

## Key Features

### 1. Data Integration
- **Gene Expression**: Loads and normalizes RNA-seq count data using DESeq2 variance stabilizing transformation
- **DNA Methylation**: Processes filtered CpG methylation percentages from WGBS data
- **Genomic Annotations**: Integrates genome GFF annotations and UTR annotations

### 2. Genic Region Classification

The script annotates CpGs to five distinct genic regions:

- **Promoter**: 1kb upstream of gene transcription start site
- **5' UTR**: 5' untranslated region (up to 1kb)
- **Exons**: Protein-coding exonic sequences
- **Introns**: Intronic regions within gene boundaries
- **3' UTR**: 3' untranslated region (up to 1kb)

### 3. Methylation Calculation

For each gene and region type, calculates:
- Average methylation across all CpG sites in that region
- Handles missing data appropriately
- Creates region-specific methylation profiles

### 4. Statistical Modeling

Implements three complementary modeling approaches:

#### Linear Regression
- Simple multiple regression model
- Provides interpretable coefficients for each region
- Shows direction and magnitude of methylation effects

#### Elastic Net (Regularized Regression)
- Combines L1 (LASSO) and L2 (Ridge) regularization
- Performs feature selection to identify most important regions
- Cross-validation to optimize regularization parameter
- Handles multicollinearity between regional methylation values

#### Random Forest
- Non-parametric ensemble learning method
- Captures non-linear relationships
- Provides feature importance rankings
- Robust to overfitting

### 5. Visualizations

Creates comprehensive visualizations including:

1. **Model Coefficient Plots**: Shows effect size of each region on expression
2. **Scatter Plots**: Methylation vs. expression for each region with regression lines
3. **Correlation Heatmap**: All pairwise correlations between regions and expression
4. **Feature Importance Plot**: Random forest variable importance rankings
5. **Distribution Boxplots**: Methylation level distributions by region
6. **Model Comparison**: Performance metrics across all three models

### 6. Output Files

Saves multiple result files:

- `gene-methylation-expression-features.csv`: Complete feature matrix
- `cpg-genic-annotations.csv`: CpG-to-region mapping
- `model-performance-comparison.csv`: Model metrics (R², RMSE)
- `random-forest-feature-importance.csv`: Feature rankings
- `elastic-net-coefficients.csv`: Selected features and coefficients
- `methylation-by-region-statistics.csv`: Summary statistics

## Dependencies

Required R packages:
- `tidyverse`: Data manipulation and visualization
- `DESeq2`: RNA-seq normalization
- `rtracklayer`: GFF/GTF file parsing
- `GenomicRanges`: Genomic interval operations
- `glmnet`: Elastic net regression
- `randomForest`: Random forest modeling
- `pheatmap`: Heatmap visualization
- `caret`: Model training utilities

## Usage

### Running the Complete Analysis

```r
# In R or RStudio
rmarkdown::render("D-Apul/code/39-Apul-genebody-meth-expr-modeling.Rmd")
```

### Interactive Exploration

Open the Rmd file in RStudio and run chunks sequentially to:
- Inspect intermediate results
- Modify parameters (e.g., promoter size, filtering thresholds)
- Generate additional visualizations

## Scientific Interpretation

### Expected Patterns

Based on literature, we might expect:

1. **Promoter Methylation**: Often negatively correlated with expression (gene silencing)
2. **Gene Body Methylation**: Can show positive correlation with highly expressed genes
3. **5' UTR**: May affect translation efficiency
4. **3' UTR**: Can influence mRNA stability and expression
5. **Exon vs Intron**: Different methylation patterns reflecting functional constraints

### Model Interpretation

- **Linear Model**: Provides baseline understanding of regional effects
- **Elastic Net**: Identifies which regions are most predictive when controlling for multicollinearity
- **Random Forest**: Captures complex, potentially non-linear relationships

### Biological Insights

The analysis will reveal:
- Which genic regions show strongest methylation-expression relationships
- Whether effects are consistent across genes or context-dependent
- Relative importance of different regulatory regions

## Extensions and Future Work

The script provides a foundation for:

1. **Time-Point Analysis**: Examine how methylation-expression relationships change across developmental/stress timepoints
2. **Gene-Set Enrichment**: Focus on specific functional gene categories
3. **Multi-Omic Integration**: Combine with miRNA and lncRNA predictors
4. **Species Comparison**: Apply to P. evermanni and P. tuahiniensis data

## Notes

- Script follows repository naming convention: `39-Apul-genebody-meth-expr-modeling.Rmd`
- Uses existing processed data (no raw data reprocessing needed)
- Comprehensive documentation and comments throughout
- Reproducible with seed set for random processes
- Performance metrics allow model comparison

## Contact

For questions or issues with this script, please open an issue in the GitHub repository.
