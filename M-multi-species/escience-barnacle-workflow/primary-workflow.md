# Primary Workflow

(main explanatory documentation of process; tangential details and complementary information can be found in other files in this directory)

## Expression Data

Count data from <https://urol-e5.github.io/MOSAiC/counts.html>

was transformed.

``` bash        
# Load ortholog annotations
orthos <- read_csv("https://raw.githubusercontent.com/urol-e5/timeseries_molecular/refs/heads/main/M-multi-species/output/12-ortho-annot/ortholog_groups_annotated.csv")

# Load species-specific count matrices
apul <- read_csv("https://gannet.fish.washington.edu/gitrepos/urol-e5/timeseries_molecular/D-Apul/output/02.20-D-Apul-RNAseq-alignment-HiSat2/apul-gene_count_matrix.csv")
peve <- read_csv("https://gannet.fish.washington.edu/gitrepos/urol-e5/timeseries_molecular/E-Peve/output/02.20-E-Peve-RNAseq-alignment-HiSat2/peve-gene_count_matrix.csv")
ptua <- read_csv("https://gannet.fish.washington.edu/gitrepos/urol-e5/timeseries_molecular/F-Ptua/output/02.20-F-Ptua-RNAseq-alignment-HiSat2/ptua-gene_count_matrix.csv")
```

> Ortholog groups were selected and merged across species to create a combined count matrix containing 9,800 orthologs present at non-NA values in all 117 samples.
> Metadata was parsed from sample names to extract species and timepoint information.

### Variance Stabilizing Transformation (VST)

Low-count and incomplete ortholog rows were removed. A variance stabilizing transformation was applied using DESeq2::vst() to normalize count variance across the dynamic range blind to study design.

``` r
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData   = metadata,
  design    = ~ species * timepoint
)

vst_obj <- vst(dds, blind = TRUE)
vst_mat <- assay(vst_obj)
```

```
../output/14-pca-orthologs/vst_counts_matrix.csv
```

<img src="http://gannet.fish.washington.edu/seashell/snaps/Monosnap_Image_2025-11-26_05-34-30.png" style="width:50%;">




## Barnacle: Rank Determination

## Barnacle: Optimization
