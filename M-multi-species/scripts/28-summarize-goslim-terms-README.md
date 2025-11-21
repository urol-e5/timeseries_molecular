# GO Slim Term Summarization

## Script: 28-summarize-goslim-terms.py

This script summarizes GO Slim term occurrences across component gene annotation files.

### Purpose

Counts and summarizes the occurrences of GO Slim terms from the `goslim_names` column across all annotation files in the top genes per component directory.

### Input

- All CSV files containing "annotation" in their filename from:
  - `M-multi-species/output/26-rank35-optimization/lambda_gene_0.2/top_genes_per_component/`

### Output

Creates a single CSV file:
- `M-multi-species/output/26-rank35-optimization/lambda_gene_0.2/top_genes_per_component/goslim_term_counts.csv`

The output file contains:
- `term`: The GO Slim term name
- `Component_1` through `Component_35`: Count of the term in each component
- `total`: Total number of occurrences across all components

Results are sorted by total count in descending order (most common terms first).

### Usage

```bash
cd /path/to/timeseries_molecular
python3 M-multi-species/scripts/28-summarize-goslim-terms.py
```

### Testing

A test script is provided to validate the output:

```bash
python3 M-multi-species/scripts/test_goslim_summary.py
```

The test validates:
- Output file exists
- Correct file structure (columns: term, Component_1...Component_35, total)
- Data types are correct
- All counts are positive
- Totals match sum of component columns
- Results are sorted by total count
- Terms exist in source annotation files

### Example Output

```
term,Component_1,Component_2,Component_3,...,Component_35,total
organelle,39,27,34,...,21,966
catalytic activity,24,11,17,...,6,514
nucleus,18,20,17,...,7,452
cytosol,22,13,17,...,12,449
anatomical structure development,23,16,13,...,6,439
...
```

### Notes

- GO Slim terms in the `goslim_names` column are semicolon-separated
- The script handles missing values (NaN) appropriately
- Each term is counted once per occurrence in each row
- The script processes all annotation files found in the target directory
