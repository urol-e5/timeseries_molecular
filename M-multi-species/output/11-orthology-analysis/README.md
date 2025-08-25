# Orthologous Protein Identification for Three Coral Species

## Overview

This analysis identifies orthologous proteins across three coral species using a reciprocal best hits (RBH) approach:

- **Acropora pulchra** (D-Apul): Fast-growing, branching coral (~40,521 proteins)
- **Porites evermanni** (E-Peve): Slow-growing, massive coral (~40,389 proteins)  
- **Pocillopora tuahiniensis** (F-Ptua): Intermediate growth, branching coral (~31,840 proteins)

## Methodology

### Reciprocal Best Hits (RBH) Approach

1. **All-vs-all protein BLAST comparisons** between each species pair
2. **Best hit identification** for each query protein (highest bitscore)
3. **Reciprocal validation** - protein A in species 1 is orthologous to protein B in species 2 if:
   - A's best hit in species 2 is B, AND
   - B's best hit in species 1 is A
4. **Quality filters**:
   - E-value ≤ 1e-5
   - Sequence identity ≥ 30%
   - Query or subject coverage ≥ 50%

### Ortholog Group Creation

- **Three-way orthologs**: Proteins present in all three species with reciprocal relationships
- **Two-way orthologs**: Proteins found in only two species

## Files and Scripts

### Input Data
```
D-Apul/data/Apulchra-genome.pep.faa                    # Acropora pulchra proteins
E-Peve/data/Porites_evermanni_v1.annot.pep.fa          # Porites evermanni proteins
F-Ptua/data/Pocillopora_meandrina_HIv1.genes.pep.faa   # Pocillopora tuahiniensis proteins
```

### Analysis Scripts
```
M-multi-species/scripts/11-orthology-analysis.Rmd      # Main R Markdown analysis (full pipeline)
M-multi-species/scripts/comprehensive_orthology.py     # Python orthology identification
M-multi-species/scripts/identify_orthologs.py          # Simple pairwise ortholog finder
```

### Output Files
```
M-multi-species/output/11-orthology-analysis/
├── coral_ortholog_groups_test.csv          # Ortholog groups table
├── orthology_summary_test.csv              # Summary statistics  
├── apul_peve_rbh_test.csv                  # Apul-Peve reciprocal best hits
├── apul_ptua_rbh_test.csv                  # Apul-Ptua reciprocal best hits
├── peve_ptua_rbh_test.csv                  # Peve-Ptua reciprocal best hits
└── test_*.blastp                           # BLAST results files
```

## Results Summary (Test Data)

**Note**: These results are from a test run using subset data for demonstration:
- Apul: 350 proteins (first 350 from full dataset)
- Peve: 203 proteins (first 203 from full dataset)  
- Ptua: 1,000 proteins (first 1,000 from full dataset)

### Ortholog Counts
- **Total ortholog groups**: 42
- **Three-way orthologs**: 0 (none found in test subset)
- **Apul-Peve pairs**: 3
- **Apul-Ptua pairs**: 12  
- **Peve-Ptua pairs**: 27

### Example Ortholog Groups
```
OG_00001: FUN_000191-T1 (Apul) | Peve_00000096 (Peve) | - | 32.3% identity
OG_00005: FUN_000200-T1 (Apul) | - | Pocillopora_meandrina_HIv1___RNAseq.g28904.t1 (Ptua) | 82.3% identity
OG_00017: - | Peve_00000005 (Peve) | Pocillopora_meandrina_HIv1___RNAseq.g1185.t2 (Ptua) | 80.6% identity
```

## Running the Analysis

### Complete Analysis (Full Datasets)

⚠️ **Warning**: Full analysis requires significant computational resources (hours to days) due to large protein datasets.

```bash
# Navigate to multi-species directory
cd M-multi-species/scripts

# Run full R Markdown analysis (recommended)
R -e "rmarkdown::render('11-orthology-analysis.Rmd')"

# Or run Python scripts for specific comparisons
python3 comprehensive_orthology.py
```

### Test Analysis (Subset Data)

```bash
# The test analysis is already complete and demonstrates the approach
# Results are in M-multi-species/output/11-orthology-analysis/
```

## Biological Significance

### Applications

1. **Comparative Genomics**
   - Identify conserved genes across coral species
   - Study gene family evolution
   - Transfer functional annotations between species

2. **Expression Analysis**
   - Compare expression of orthologous genes across species
   - Identify species-specific expression patterns
   - Study conservation of regulatory networks

3. **Functional Annotation**
   - Transfer annotations from well-studied to poorly-annotated species
   - Predict protein function based on orthologous relationships
   - Identify essential vs. dispensable genes

4. **Evolutionary Analysis**
   - Reconstruct gene family evolution
   - Identify gene duplications and losses
   - Study molecular evolution rates

### Interpretation

- **High identity orthologs** (>70%): Likely essential, conserved functions
- **Medium identity orthologs** (40-70%): Conserved function with species adaptation
- **Low identity orthologs** (30-40%): Conserved domain/function, significant divergence
- **Species-specific genes**: No orthologs found, potential lineage-specific innovations

## Computational Requirements

### For Full Analysis
- **CPU**: Multi-core recommended (8+ cores)
- **Memory**: 16+ GB RAM
- **Storage**: 50+ GB for intermediate files
- **Time**: 6-24 hours depending on resources

### Dependencies
- BLAST+ (ncbi-blast+)
- R with tidyverse, here packages
- Python3 with pandas
- Optional: ggplot2, VennDiagram for visualizations

## Next Steps

1. **Run full analysis** on complete protein datasets
2. **Functional enrichment** analysis of ortholog groups
3. **Expression correlation** analysis using RNA-seq data
4. **Phylogenetic analysis** of ortholog groups
5. **Integration with other multi-omics** data (methylation, miRNA)

## Quality Assessment

The reciprocal best hits approach is conservative and tends to identify high-confidence orthologs. However:

- **False negatives**: May miss orthologs with complex evolutionary relationships
- **Paralogs**: May incorrectly identify paralogs as orthologs
- **Gene duplications**: Requires careful interpretation in species with recent duplications

For more comprehensive orthology, consider supplementing with:
- Synteny analysis
- Tree-based orthology prediction
- Domain-based approaches

## Authors and Contributions

- **Orthology pipeline development**: Automated BLAST-based approach
- **Method validation**: Reciprocal best hits with quality filters
- **Documentation**: Comprehensive analysis workflow
- **Integration**: Multi-species comparative framework