# Molecular Methods Documentation

This directory contains comprehensive methods documentation for the timeseries molecular multi-omic analyses.

## Contents

### molecular-methods-section.md
Complete methods section describing the workflows for generating genomic count matrices from sample extraction to quality-controlled count matrices. This document covers:

1. **Sample Collection and Preparation**
   - Experimental design across three coral species
   - Nucleic acid extraction protocols

2. **RNA-seq: Gene Expression Quantification**
   - Quality control and read trimming
   - Genome indexing and alignment (HISAT2)
   - Transcript assembly and quantification (StringTie)
   - Count matrix generation

3. **Small RNA-seq: miRNA Discovery and Quantification**
   - sRNA-specific quality control and trimming
   - miRNA discovery and annotation (ShortStack)
   - Count matrix formatting and standardization

4. **Long Non-coding RNA (lncRNA) Discovery and Quantification**
   - Transcript assembly and classification (GFFcompare)
   - Coding potential assessment (CPC2)
   - lncRNA quantification (featureCounts)
   - Quality filtering approaches

5. **DNA Methylation: Whole Genome Bisulfite Sequencing (WGBS)**
   - Bisulfite genome preparation (Bismark)
   - Read alignment and deduplication
   - Methylation calling and extraction
   - Coverage filtering and count matrix generation

6. **Quality Control and Filtering Summary**
   - General QC principles
   - Data type-specific filtering criteria

7. **Count Matrix Generation and Normalization Philosophy**
   - Raw count matrix standards
   - Matrix structure conventions
   - Downstream normalization approaches

8. **Data Integration and Quality Assurance**
   - Sample naming standardization
   - Metadata association
   - Reproducibility documentation

## Usage

This methods section is designed for inclusion in manuscripts, supplementary materials, or protocols describing the timeseries molecular analyses. The comprehensive documentation ensures reproducibility and transparency of the analytical workflows.

## Related Files

Key analysis scripts referenced in the methods:
- `M-multi-species/scripts/10-format-miRNA-counts.R` - miRNA count matrix formatting
- `M-multi-species/scripts/11-filter_count_matrix.py` - Count matrix quality filtering
- Species-specific analysis pipelines in `D-Apul/code/`, `E-Peve/code/`, and `F-Ptua/code/` directories

## Contact

For questions about specific methods or analysis details, please refer to the main repository README or contact the appropriate analysis contributor listed in the main documentation.
