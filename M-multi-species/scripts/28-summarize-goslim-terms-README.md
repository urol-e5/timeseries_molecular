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

The output file contains two columns:
- `term`: The GO Slim term
- `count`: Number of occurrences across all annotation files

Results are sorted by count in descending order (most common terms first).

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
- Correct file structure (columns: term, count)
- Data types are correct
- All counts are positive
- Results are sorted by count
- Terms exist in source annotation files

### Example Output

```
term,count
organelle,966
catalytic activity,514
nucleus,452
cytosol,449
anatomical structure development,439
...
```

### Notes

- GO Slim terms in the `goslim_names` column are semicolon-separated
- The script handles missing values (NaN) appropriately
- Each term is counted once per occurrence in each row
- The script processes all annotation files found in the target directory
