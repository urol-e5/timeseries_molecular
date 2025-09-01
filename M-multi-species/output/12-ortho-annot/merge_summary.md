# Ortholog Groups Annotation Merge Summary

## Overview
Successfully merged the ortholog groups data with functional annotations to create a fully annotated ortholog groups file.

## Files Merged
- **Source ortholog groups**: `../11-orthology-analysis/ortholog_groups.csv`
- **Source annotations**: `run_20250831_172744/annotation_with_goslim.tsv`
- **Output file**: `ortholog_groups_annotated.csv`

## Merge Details
- **Join type**: Left join (kept all records from ortholog_groups)
- **Join key**: FUN ID (`apul` column from ortholog_groups ↔ `query` column from annotations)
- **Total ortholog groups**: 18,326
- **Ortholog groups with annotations**: 10,348 (56.5%)
- **Ortholog groups without annotations**: 8,024 (43.5%)

## Output File Structure
The merged file contains 23 columns:

### Original Ortholog Groups Columns (6)
1. `group_id` - Ortholog group identifier
2. `apul` - FUN ID (Fungia scutaria)
3. `peve` - Peve ID (Pocillopora verrucosa)
4. `ptua` - Ptua ID (Pocillopora meandrina)
5. `type` - Ortholog type (e.g., "three_way")
6. `avg_identity` - Average sequence identity

### Annotation Columns (17)
7. `query` - Query ID (matches apul)
8. `accession` - UniProt accession
9. `id` - UniProt ID
10. `reviewed` - Review status
11. `protein_name` - Protein name
12. `organism` - Source organism
13. `pident` - Percent identity
14. `length` - Alignment length
15. `evalue` - E-value
16. `bitscore` - Bit score
17. `title` - Full title
18. `go_ids` - GO term IDs
19. `go_bp` - Biological process GO terms
20. `go_cc` - Cellular component GO terms
21. `go_mf` - Molecular function GO terms
22. `goslim_ids` - GO Slim term IDs
23. `goslim_names` - GO Slim term names

## Usage
The merged file can be used for:
- Functional analysis of ortholog groups
- GO enrichment analysis
- Comparative genomics studies
- Pathway analysis

## File Size
- **Output file**: 9.7 MB
- **Format**: CSV with comma separation
