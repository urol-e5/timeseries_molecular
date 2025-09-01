# Ortholog Annotation Results

This directory contains the results of the ortholog group functional annotation analysis for three coral species.

## Analysis Summary

The annotation analysis successfully processed **18,326 ortholog groups** across three coral species:

- **Acropora pulchra** (Apul): Fast-growing, branching coral
- **Porites evermanni** (Peve): Slow-growing, massive coral  
- **Pocillopora tuahiniensis** (Ptua): Intermediate growth, branching coral

## Key Results

### Ortholog Distribution
- **Total Ortholog Groups**: 18,326
- **Three-way Orthologs**: 10,346 (56.5%) - Present in all three species
- **Two-way Orthologs**: 7,980 (43.5%) - Present in specific species pairs

### Two-way Ortholog Breakdown
- **Apul-Peve**: 3,436 (18.7%)
- **Peve-Ptua**: 2,662 (14.5%)
- **Apul-Ptua**: 1,882 (10.3%)

### Conservation Levels
- **High Conservation (>80% identity)**: 3,427 (18.7%)
- **Medium Conservation (60-80% identity)**: 8,308 (45.3%)
- **Low Conservation (<60% identity)**: 6,591 (36.0%)

### Functional Categories
- **Core Conserved Functions**: 10,346 (56.5%) - Three-way orthologs
- **Lineage-specific Functions**: 7,980 (43.5%) - Two-way orthologs

### Sequence Identity Statistics
- **Overall Average Identity**: 65.4%
- **Three-way Orthologs Average**: 69.9%
- **Two-way Orthologs Average**: 59.5%

## Files Description

### Core Results
- `basic_annotations.csv` - Complete annotation table with all ortholog groups
- `annotation_summary_report.csv` - Summary statistics of the analysis
- `ortholog_annotations_database.csv` - Annotation database for downstream analysis

### Gene Mapping
- `gene_to_ortholog_mapping.csv` - Mapping between gene IDs and ortholog groups
- `representative_sequences.faa` - Representative protein sequences for each ortholog group

### Distribution Analysis
- `ortholog_type_distribution.csv` - Distribution of ortholog types
- `conservation_level_distribution.csv` - Distribution of conservation levels
- `functional_category_distribution.csv` - Distribution of functional categories

### Visualizations
- `ortholog_type_distribution.png` - Pie chart of ortholog types
- `conservation_level_distribution.png` - Bar chart of conservation levels
- `functional_category_distribution.png` - Bar chart of functional categories
- `identity_distribution_by_type.png` - Histogram of identity distributions

## Applications

The annotated ortholog groups can be used for:

1. **Comparative Functional Genomics**
   - Understanding functional conservation across coral species
   - Identifying species-specific functional adaptations
   - Studying functional evolution of gene families

2. **Pathway Analysis**
   - Identifying conserved metabolic and signaling pathways
   - Studying pathway evolution across species
   - Understanding pathway-specific adaptations

3. **Gene Set Enrichment Analysis**
   - Finding overrepresented functional terms in expression studies
   - Identifying biological processes affected by experimental conditions
   - Studying functional responses to environmental stress

4. **Cross-Species Expression Analysis**
   - Comparing expression of functionally annotated orthologs
   - Studying expression conservation and divergence
   - Identifying expression patterns associated with functional categories

5. **Evolutionary Analysis**
   - Studying functional evolution of gene families
   - Identifying functional innovations and losses
   - Understanding evolutionary constraints on function

## Usage Examples

### Load Annotation Database
```r
# Load annotation database
annotations <- read.csv("ortholog_annotations_database.csv")

# Filter for three-way orthologs
three_way_orthologs <- annotations[annotations$type == "three_way", ]

# Filter for high conservation orthologs
high_conservation <- annotations[annotations$conservation_level == "high", ]
```

### Load Gene Mapping
```r
# Load gene mapping
gene_mapping <- read.csv("gene_to_ortholog_mapping.csv")

# Find ortholog group for a specific gene
gene_ortholog <- gene_mapping[gene_mapping$gene_id == "FUN_000185-T1", ]

# Find all genes in a specific ortholog group
group_genes <- gene_mapping[gene_mapping$group_id == "OG_00001", ]
```

### Cross-Species Analysis
```r
# Load gene mapping
gene_mapping <- read.csv("gene_to_ortholog_mapping.csv")

# Get orthologs for a specific species
apul_genes <- gene_mapping[gene_mapping$species == "Apul", ]
peve_genes <- gene_mapping[gene_mapping$species == "Peve", ]
ptua_genes <- gene_mapping[gene_mapping$species == "Ptua", ]

# Find three-way orthologs
three_way_groups <- unique(gene_mapping[gene_mapping$group_id %in% 
  gene_mapping$group_id[duplicated(gene_mapping$group_id)], "group_id"])
```

## Quality Metrics

### Annotation Coverage
- **Total Ortholog Groups**: 18,326 (100%)
- **Three-way Orthologs**: 10,346 (56.5%)
- **Two-way Orthologs**: 7,980 (43.5%)

### Conservation Analysis
- **High Conservation**: 3,427 orthologs (18.7%)
- **Medium Conservation**: 8,308 orthologs (45.3%)
- **Low Conservation**: 6,591 orthologs (36.0%)

### Confidence Levels
- **High Confidence**: Three-way orthologs with high conservation
- **Medium Confidence**: Three-way orthologs with medium conservation
- **Low Confidence**: Two-way orthologs or low conservation orthologs

## Next Steps

1. **Functional Annotation Enhancement**
   - Run InterProScan for detailed protein domain analysis
   - Perform BLAST searches against Swiss-Prot database
   - Add Gene Ontology (GO) term annotations

2. **Expression Analysis**
   - Map ortholog groups to expression data
   - Perform cross-species expression comparisons
   - Identify conserved expression patterns

3. **Pathway Analysis**
   - Identify conserved metabolic pathways
   - Study pathway evolution across species
   - Analyze pathway-specific adaptations

4. **Evolutionary Analysis**
   - Study gene family evolution
   - Identify functional innovations and losses
   - Understand evolutionary constraints

## Citation

If you use these annotations in your research, please cite:

- **Orthology Analysis**: Multi-species comparative analysis (2024)
- **Coral Species**: Acropora pulchra, Porites evermanni, Pocillopora tuahiniensis
- **Analysis Pipeline**: Reciprocal Best Hits (RBH) approach

---

*This annotation provides a comprehensive foundation for comparative genomics studies across coral species, enabling detailed analysis of gene function and evolution.*
