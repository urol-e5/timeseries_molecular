# Gene List Annotation Script

## Overview

The `23-annotate-rank-lists.py` script annotates rank CSV files with ortholog group annotations from the comprehensive ortholog database.

## What it does

1. **Reads all CSV files** from `M-multi-species/output/22-Visualizing-Rank-outs/`
2. **Joins each file** with `M-multi-species/output/12-ortho-annot/ortholog_groups_annotated.csv` based on OG ID (ortholog group identifier)
3. **Creates annotated files** with suffix `_annotation.csv` containing:
   - Original data (OG ID and values)
   - Ortholog information (apul, peve, ptua gene IDs)
   - Protein annotations (protein name, organism, UniProt IDs)
   - GO annotations (GO biological process, cellular component, molecular function)
   - GOSlim annotations
4. **Generates a summary table** (`annotation_summary.md`) showing:
   - Total genes in each file
   - Number of annotated genes
   - Predominant GO Biological Processes (top 3)
   - Predominant Gene Functions (top 3)

## Usage

```bash
python3 M-multi-species/scripts/23-annotate-rank-lists.py
```

The script will:
- Process all 524 CSV files in the rank-outs directory
- Create 524 corresponding `_annotation.csv` files
- Generate an `annotation_summary.md` markdown table

## Output Files

### Annotation Files
Each rank CSV file (e.g., `rank_55_comp49_top100.csv`) gets a corresponding annotation file (`rank_55_comp49_top100_annotation.csv`) with the following columns:

1. `group_id` - Ortholog group ID (e.g., OG_00531)
2. Original data columns (varies by file)
3. `apul` - A. pulchra gene ID
4. `peve` - P. evermanni gene ID  
5. `ptua` - P. tuahiniensis gene ID
6. `type` - Ortholog type (e.g., three_way)
7. `avg_identity` - Average identity percentage
8. `query` - Query gene ID
9. `accession` - UniProt accession
10. `id` - UniProt ID
11. `reviewed` - UniProt review status
12. `protein_name` - Protein name and description
13. `organism` - Source organism
14. `pident`, `length`, `evalue`, `bitscore` - BLAST statistics
15. `title` - Full protein title
16. `go_ids` - GO term IDs
17. `go_bp` - GO Biological Process descriptions
18. `go_cc` - GO Cellular Component descriptions
19. `go_mf` - GO Molecular Function descriptions
20. `goslim_ids` - GOSlim IDs
21. `goslim_names` - GOSlim term names

### Summary Table
The `annotation_summary.md` file contains a markdown table with one row per input file showing:
- File name
- Total number of genes
- Number of genes with annotations
- Top 3 most common GO Biological Processes
- Top 3 most common Gene Functions

## Example

For a file like `rank_55_comp49_top100.csv` containing:
```csv
OG,Component_49
OG_00531,0.879
OG_07144,0.811
...
```

The output `rank_55_comp49_top100_annotation.csv` will contain:
```csv
group_id,Component_49,apul,peve,ptua,type,avg_identity,...,protein_name,...,go_bp,...
OG_00531,0.879,FUN_001545-T1,Peve_00033291,...,three_way,49.81,...,Protein name,...,GO processes,...
OG_07144,0.811,FUN_032474-T1,Peve_00036450,...,three_way,42.55,...,Protein name,...,GO processes,...
...
```

## Dependencies

- Python 3.x
- pandas

## Notes

- The script automatically skips files that end with `_annotation.csv` to avoid re-processing
- Genes without annotations will have empty values in the annotation columns
- The script uses a left join, so all original genes are retained even if they lack annotations
