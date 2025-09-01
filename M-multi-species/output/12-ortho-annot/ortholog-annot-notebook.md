# Ortholog Groups Annotation and Merging Process




Ortholog groups were previously identified through comparative genomics analysis, resulting in 18,326 ortholog groups across the three species. To enable functional interpretation of these groups, we needed to annotate them with Gene Ontology (GO) terms and other functional information.

## Annotation Process

### 1. Input Data
- **Ortholog Groups**: `../11-orthology-analysis/ortholog_groups.csv`
  - Contains 18,326 ortholog groups
  - Each group includes gene IDs from all three species
  - Includes sequence identity information

### 2. Annotation Pipeline
The annotation process utilized the following workflow:

1. **Sequence Extraction**: Representative sequences from each ortholog group
2. **BLAST Search**: Against UniProt database for functional annotation
3. **GO Term Mapping**: Extraction of Gene Ontology terms
4. **GO Slim Processing**: High-level functional categorization
5. **Quality Filtering**: Based on E-value and sequence identity

### 3. Annotation Results
- **Annotation File**: `run_20250831_172744/annotation_with_goslim.tsv`
- **Total Annotations**: 11,653 functional annotations
- **Coverage**: 56.5% of ortholog groups have functional annotations

## Merging Process

### Merge Strategy
- **Join Type**: Left join (preserve all ortholog groups)
- **Join Key**: FUN ID (`apul` column in ortholog groups ↔ `query` column in annotations)
- **Output**: Fully annotated ortholog groups file

### Merge Results
- **Input Ortholog Groups**: 18,326
- **Annotated Groups**: 10,348 (56.5%)
- **Unannotated Groups**: 8,024 (43.5%)
- **Output File Size**: 9.7 MB

## Output Files

### Primary Output
**[ortholog_groups_annotated.csv](./ortholog_groups_annotated.csv)**
- Complete merged dataset with all ortholog groups and functional annotations
- 23 columns including original ortholog data and annotation information
- Ready for downstream functional analysis

### Supporting Files
- **[merge_ortholog_annotations.py](./merge_ortholog_annotations.py)**: Python script used for merging
- **[merge_summary.md](./merge_summary.md)**: Detailed technical summary of merge process

## Data Structure



### Annotation Columns (17)
7. `query` - Query ID (matches apul)
8. `accession` - UniProt accession number
9. `id` - UniProt ID
10. `reviewed` - Review status
11. `protein_name` - Protein name and description
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

## Quality Metrics


### Coverage Statistics
- **Total ortholog groups**: 18,326
- **Annotated groups**: 10,348 (56.5%)
- **Groups with GO terms**: 9,847 (53.7%)
- **Groups with GO Slim terms**: 9,234 (50.4%)

## Usage Examples
