`timeseries_molecular/M-multi-species/output/13.00-multiomics-barnacle`

## Files

### Ortholog Expression Data

- **`apul_ortholog_expression.csv`** - Acropora pulchra ortholog gene expression data
- **`peve_ortholog_expression.csv`** - Pocillopora verrucosa ortholog gene expression data  
- **`ptua_ortholog_expression.csv`** - Pocillopora meandrina ortholog gene expression data

Each CSV file contains:
- `gene_id` column with species-specific gene identifiers
- Expression count columns for all RNA-seq samples
- Only genes present in complete three-way ortholog groups (genes with orthologs in all three species)
- Only genes with expression data available in all three species

