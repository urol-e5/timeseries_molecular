# WGCNA Module Annotation

## Overview

This script annotates WGCNA modules with functional information by joining module assignments with ortholog group annotations.

## Script

`19-annotate-wgcna-modules.py`

## Purpose

The script performs the following tasks:
1. Joins `wgcna_ortholog_module_assignments.csv` with `ortholog_groups_annotated.csv` on OG ID (group_id)
2. Summarizes dominant GO processes, gene functions, and physiological pathways for each WGCNA module
3. Generates comprehensive annotation reports in the output directory

## Input Files

- **WGCNA Module Assignments**: `M-multi-species/output/18-ortholog-wgcna/wgcna_ortholog_module_assignments.csv`
  - Contains ortholog group IDs and their assigned WGCNA modules
  - Format: `group_id,wgcna_module` (e.g., `OG_00705,7`)

- **Ortholog Annotations**: `M-multi-species/output/12-ortho-annot/ortholog_groups_annotated.csv`
  - Contains functional annotations for ortholog groups
  - Includes GO terms (BP, CC, MF), GO Slim terms, and protein names

## Output Files

All output files are saved to: `M-multi-species/output/18-ortholog-wgcna/`

### 1. Detailed Text Report
- **File**: `wgcna_module_annotation_summary.txt`
- **Content**: Comprehensive summary for each WGCNA module including:
  - Total orthologs and annotation coverage
  - Top 10 GO Biological Process terms
  - Top 10 GO Cellular Component terms
  - Top 10 GO Molecular Function terms
  - Top 10 GO Slim terms (physiological pathways)
  - Top 10 protein names (gene functions)

### 2. Module Overview
- **File**: `wgcna_module_overview.csv`
- **Content**: Summary statistics for all modules in tabular format
- **Columns**:
  - `module`: Module ID (0-14)
  - `num_orthologs`: Total orthologs in module
  - `num_annotated`: Orthologs with functional annotations
  - `annotation_coverage_%`: Percentage of orthologs with annotations
  - `num_go_bp_terms`: Number of unique GO BP terms
  - `num_go_cc_terms`: Number of unique GO CC terms
  - `num_go_mf_terms`: Number of unique GO MF terms
  - `num_goslim_terms`: Number of unique GO Slim terms
  - `num_unique_proteins`: Number of unique protein names
  - `top_go_bp`: Most common GO BP term
  - `top_goslim`: Most common GO Slim term
  - `top_protein`: Most common protein name

### 3. GO Term Summary Files
Each file contains term counts across all modules:

- **`wgcna_module_go_bp_summary.csv`**: GO Biological Process terms
- **`wgcna_module_go_cc_summary.csv`**: GO Cellular Component terms
- **`wgcna_module_go_mf_summary.csv`**: GO Molecular Function terms
- **`wgcna_module_goslim_summary.csv`**: GO Slim terms (physiological pathways)
- **`wgcna_module_protein_summary.csv`**: Protein names (gene functions)

**Format**: Each file has columns:
- `term`: The GO term or protein name
- `module_0` through `module_14`: Count of term occurrences in each module
- `total`: Total occurrences across all modules

Terms are sorted by total count (descending), showing the top 50 most common terms.

## Usage

### Run the annotation script:
```bash
cd M-multi-species/scripts
python3 19-annotate-wgcna-modules.py
```

### Test the output:
```bash
cd M-multi-species/scripts
python3 test_wgcna_annotation.py
```

## Requirements

- Python 3.10+
- pandas

Install dependencies:
```bash
pip3 install pandas
```

## Module Statistics

Based on current data:
- **Total WGCNA Modules**: 15 (numbered 0-14)
- **Total Orthologs**: 9,827 across all modules
- **Annotated Orthologs**: 3,801 (38.7% coverage)
- **Largest Module**: Module 1 with 2,649 orthologs
- **Highest Annotation Coverage**: Module 5 with 47.9% coverage

## Dominant Functional Themes

The script identifies dominant functional themes for each module:

1. **GO Biological Process**: Cellular and molecular processes (e.g., transcription regulation, cell division, proteolysis)
2. **GO Cellular Component**: Subcellular localization (e.g., nucleus, cytoplasm, plasma membrane)
3. **GO Molecular Function**: Molecular activities (e.g., ATP binding, protein binding, catalytic activity)
4. **GO Slim Terms**: High-level physiological pathways (e.g., organelle function, catalytic activity, development)
5. **Protein Names**: Specific gene functions based on protein annotations

## Notes

- Modules are numbered 0-14 (15 modules total)
- Not all orthologs have functional annotations (~39% annotation coverage overall)
- Terms are counted with multiplicity - if an ortholog has multiple GO terms, each is counted
- GO Slim terms provide a higher-level summary of physiological pathways compared to detailed GO terms
