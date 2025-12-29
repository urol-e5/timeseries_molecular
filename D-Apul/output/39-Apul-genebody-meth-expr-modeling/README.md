# Gene-Body Methylation and Expression Modeling Results

This directory contains outputs from the analysis of DNA methylation patterns in genic regions and their relationship to gene expression.

## Script

- Source: `../../code/39-Apul-genebody-meth-expr-modeling.Rmd`

## Expected Output Files

After running the script, this directory will contain:

### Data Files

- `gene-methylation-expression-features.csv`: Complete feature matrix with gene expression and regional methylation values for all genes and samples
- `cpg-genic-annotations.csv`: Annotations mapping each CpG site to genomic regions (promoter, 5'UTR, 3'UTR, exon, intron)
- `methylation-by-region-statistics.csv`: Summary statistics (mean, median, SD, min, max) for methylation levels in each genomic region

### Model Results

- `model-performance-comparison.csv`: Performance metrics (R-squared, RMSE) comparing linear regression, elastic net, and random forest models
- `random-forest-feature-importance.csv`: Feature importance scores from random forest model
- `elastic-net-coefficients.csv`: Non-zero coefficients from elastic net model showing which regions most strongly predict expression

## Analysis Overview

The analysis:

1. Loads gene expression and CpG methylation data
2. Annotates CpG sites to genic regions using genome annotations
3. Calculates average methylation for each gene in each region type
4. Builds predictive models relating regional methylation to gene expression
5. Evaluates model performance and identifies most important genomic regions

## Genomic Regions Analyzed

- **Promoter**: 1kb upstream of gene start
- **5' UTR**: 5' untranslated region (up to 1kb)
- **Exons**: Protein-coding exonic regions
- **Introns**: Intronic regions within genes
- **3' UTR**: 3' untranslated region (up to 1kb)

## Usage

To run the analysis:

```r
# In R or RStudio
rmarkdown::render("../../code/39-Apul-genebody-meth-expr-modeling.Rmd")
```

Or to run interactively, open the Rmd file in RStudio and run chunks sequentially.
