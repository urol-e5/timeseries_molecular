# Top 100 Genes Annotation

This directory contains annotated top 100 genes from BARNACLE tensor factorization analysis across different rank values.

## Directory Structure

Each `rank_*` subdirectory contains:

- **`top_100_genes_annotated.csv`** - Full annotation data for the top 100 genes by importance
  - Includes gene IDs, importance scores, species mappings (apul, peve, ptua)
  - Protein names, organisms, and Gene Ontology (GO) terms
  - GO categories: biological process (go_bp), cellular component (go_cc), molecular function (go_mf)

- **`annotation_summary.txt`** - Human-readable summary report
  - Total genes and annotation coverage statistics
  - Ortholog type distribution
  - Top organisms represented in annotations
  - Top 10 genes by importance with descriptions

## Data Source

The annotations are derived from:
- **Top genes**: `top_100_genes.csv` files from barnacle-sjw branch (rank_* directories)
- **Annotations**: `M-multi-species/output/12-ortho-annot/ortholog_groups_annotated.csv`

## Rank Directories

The following rank values were analyzed:
- rank_05 (43/100 genes annotated)
- rank_08 (48/100 genes annotated)
- rank_10 (34/100 genes annotated)
- rank_12 (54/100 genes annotated)
- rank_15 (47/100 genes annotated)
- rank_20 (37/100 genes annotated)
- rank_25 (38/100 genes annotated)
- rank_35 (32/100 genes annotated)
- rank_45 (29/100 genes annotated)
- rank_55 (29/100 genes annotated)
- rank_65 (30/100 genes annotated)
- rank_75 (23/100 genes annotated)

## Usage

To view annotation details for a specific rank:
```bash
# View summary
cat M-multi-species/output/13.00-multiomics-barnacle/rank_12/annotation_summary.txt

# View full annotations
head M-multi-species/output/13.00-multiomics-barnacle/rank_12/top_100_genes_annotated.csv
```

## Notes

- All genes in the top 100 lists are three-way orthologs (present in all three species)
- Not all genes have functional annotations available in the database
- Annotation coverage varies by rank, ranging from 23% to 54%
- Most annotations are from well-studied model organisms (Human, Mouse, Zebrafish, Rat)
