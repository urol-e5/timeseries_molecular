# Timeseries Molecular Analysis Repository

This repository contains comprehensive multi-omic time series analyses of three coral species examining molecular responses and regulatory networks across developmental/stress timepoints.

## Repository Organization

The repository is organized by species with shared multi-species analyses:

```
timeseries_molecular/
├── D-Apul/           # Acropora pulchra analyses
│   ├── code/         # Analysis scripts (R Markdown)
│   ├── data/         # Raw and processed data
│   ├── output/       # Analysis outputs
│   └── README.md     # Species-specific documentation
├── E-Peve/           # Porites evermanni analyses  
│   ├── code/         # Analysis scripts
│   ├── data/         # Raw and processed data
│   ├── output/       # Analysis outputs
│   └── README.md     # Species-specific documentation
├── F-Ptua/           # Pocillopora tuahiniensis analyses
│   ├── code/         # Analysis scripts
│   ├── data/         # Raw and processed data  
│   ├── output/       # Analysis outputs
│   └── README.md     # Species-specific documentation
├── M-multi-species/  # Cross-species comparative analyses
│   ├── data/         # Shared metadata and multi-species datasets
│   ├── scripts/      # Multi-species analysis scripts
│   └── output/       # Comparative analysis outputs
└── README.md         # This file
```

## File Naming Convention and Interaction Guidelines

### Code File Naming
All analysis scripts follow this standardized format:
```
XX.YY-<species_prefix>-<analysis_description>.Rmd
```

Where:
- `XX.YY`: Hierarchical numbering (XX = major analysis type, YY = sub-analysis)
- `<species_prefix>`: Species designation (D-Apul, E-Peve, F-Ptua)  
- `<analysis_description>`: Brief description of analysis functionality
- `.Rmd`: R Markdown format (other formats like .qmd accepted)

**Examples:**
- `00.00-D-Apul-RNAseq-reads-FastQC-MultiQC.Rmd` - Initial RNA-seq quality control
- `01.00-D-Apul-RNAseq-trimming-fastp-FastQC-MultiQC.Rmd` - Read trimming and QC
- `22.6-D-Apul-multiomic-machine-learning-updatedWGBS.Rmd` - Advanced multi-omic ML

### Output Organization
- Each code file generates outputs in corresponding `../output/` directories
- Output directories match script names (minus file extension)
- Intermediate files, plots, and processed data stored in output directories

### How Users Should Interact With Files

1. **Start with species-specific READMEs** in each species directory for detailed workflows
2. **Follow numerical order** of scripts within each species for logical analysis progression
3. **Check dependencies** - later numbered scripts often depend on outputs from earlier ones
4. **Use metadata files** in `M-multi-species/data/` for sample information and experimental design
5. **Contribute new analyses** by following naming conventions and documenting in appropriate READMEs

## Biological Hypotheses and Research Goals

This research addresses several key biological hypotheses about coral molecular responses:

### Primary Hypotheses
1. **Temporal Molecular Dynamics**: Coral gene expression, miRNA regulation, and DNA methylation patterns change systematically across developmental/stress timepoints
2. **Multi-omic Integration**: Transcriptomic, epigenomic, and metabolomic data are interconnected and can predict physiological responses
3. **Regulatory Networks**: miRNA-mRNA interactions and DNA methylation create complex regulatory networks that control coral responses
4. **Species-Specific Responses**: Different coral species exhibit distinct molecular strategies for responding to environmental conditions
5. **Predictive Modeling**: Machine learning approaches can identify molecular signatures predictive of coral physiological state

### Research Components

#### Experimental Design
- **Species**: Three coral species with different life strategies
  - *Acropora pulchra* (fast-growing, branching)
  - *Porites evermanni* (slow-growing, massive) 
  - *Pocillopora tuahiniensis* (intermediate growth, branching)
- **Timepoints**: Four seasonal sampling timepoints in 2020
  - **TP1**: January 2020 (Austral Summer)
  - **TP2**: March 2020 (Austral Summer/Autumn)
  - **TP3**: September 2020 (Austral Winter/Spring)
  - **TP4**: November 2020 (Austral Spring)
- **Location**: Mo'orea, French Polynesia (three lagoon sites on North shore)
- **Colonies**: Multiple tagged colonies per species tracked across all timepoints
- **Data Types**: RNA-seq, sRNA-seq/miRNA, WGBS (methylation), metabolomics, lipidomics, physiological measurements

#### Analysis Approaches
1. **Quality Control & Processing**: Read QC, trimming, alignment, quantification
2. **Differential Expression**: Time-series gene and miRNA expression analysis  
3. **Target Prediction**: miRNA-mRNA interaction prediction and validation
4. **Co-expression Networks**: WGCNA and correlation network analysis
5. **Epigenetic Analysis**: DNA methylation patterns and CpG island annotation
6. **Machine Learning**: Predictive modeling of phenotype from molecular data
7. **Multi-omic Integration**: Cross-platform correlation and network analysis
8. **Functional Annotation**: GO enrichment and pathway analysis

## Schematic Overview

```
┌─────────────────────────────────────────────────────────────────┐
│                    TIMESERIES MOLECULAR PIPELINE                │
└─────────────────────────────────────────────────────────────────┘

                           ┌─────────────┐
                           │   SAMPLES   │
                           │   TP1-TP4   │
                           │ 3 Species   │
                           │Multi-Colony │
                           └─────┬───────┘
                                 │
                   ┌─────────────┼─────────────┐
                   │             │             │
              ┌────▼───┐     ┌───▼───┐    ┌───▼────┐
              │D-Apul  │     │E-Peve │    │F-Ptua  │
              │A.pulchra│     │P.everm│    │P.tuahin│
              └────┬───┘     └───┬───┘    └───┬────┘
                   │             │            │
           ┌───────┼─────────────┼────────────┼───────┐
           │       │             │            │       │
      ┌────▼──┐ ┌──▼──┐    ┌────▼──┐    ┌───▼───┐ ┌──▼──┐
      │RNA-seq│ │sRNA │    │ WGBS  │    │Metabol│ │Lipid│
      │       │ │miRNA│    │Methyl │    │ omics │ │omics│
      └───┬───┘ └──┬──┘    └───┬───┘    └───┬───┘ └──┬──┘
          │        │           │            │        │
          └────────┼───────────┼────────────┼────────┘
                   │           │            │
               ┌───▼───────────▼────────────▼───┐
               │     INTEGRATED ANALYSES       │
               │                               │
               │ ○ Differential Expression     │
               │ ○ Co-expression Networks      │
               │ ○ miRNA Target Prediction     │
               │ ○ Machine Learning Models     │
               │ ○ Multi-omic Correlation      │
               │ ○ Functional Annotation       │
               └───────────────┬───────────────┘
                               │
                    ┌──────────▼──────────┐
                    │   BIOLOGICAL        │
                    │   INSIGHTS          │
                    │                     │
                    │ • Regulatory Networks│
                    │ • Temporal Dynamics │
                    │ • Species Differences│
                    │ • Predictive Models │
                    └─────────────────────┘
```

## Contributors and Their Contributions

### Core Contributors

**Sam White** - Lead Data Engineer
- RNA-seq data processing and quality control
- Genome indexing and read alignment (HISAT2)
- WGBS data processing and bisulfite genome preparation
- Established computational pipelines and file organization standards

**Kathleen Durkin** - miRNA and Regulatory Networks Specialist  
- sRNA-seq/miRNA discovery and expression analysis
- miRNA target prediction (miRanda, RNAhybrid)
- miRNA-mRNA co-expression and correlation networks
- Machine learning approaches for multi-omic integration
- Gene set enrichment and functional annotation

**Steven Roberts** - Epigenomics Lead
- WGBS methylation analysis and Bismark processing
- CpG island annotation and methylation pattern analysis
- SNP calling from bisulfite sequencing data
- Multi-omic machine learning model development
- Protein annotation and comparative analyses

**Ariana S. Huffmyer** - Metabolomics and Physiological Integration
- Metabolomics data processing and analysis
- Lipidomics data analysis and integration
- Physiological data integration and PC analysis
- Cross-platform correlation analysis
- Multi-species comparative approaches

**Jill Ashey** - Network Analysis Specialist
- WGCNA analysis for metabolomics data
- Co-expression network construction and analysis

**R. Cunning** - Multi-omic Integration
- Cross-platform data integration approaches
- Metabolomics-lipidomics-genomics correlation analysis

### Analysis Categories by Contributor

| Analysis Type | Primary Contributors | Key Scripts |
|---------------|---------------------|-------------|
| Data QC & Processing | Sam White | 00.XX, 01.XX, 02.XX series |
| Gene Expression | Sam White, Kathleen Durkin | 03.XX series |
| miRNA Analysis | Kathleen Durkin | 04.XX, 06.XX, 07.XX series |
| Co-expression Networks | Kathleen Durkin, Jill Ashey | 12.XX, 16.XX, 18.XX series |
| Methylation Analysis | Steven Roberts | 14.XX, 15.XX, 19.XX, 26.XX, 28.XX series |
| Machine Learning | Kathleen Durkin, Steven Roberts | 20.XX, 22.XX, 29.XX series |
| Metabolomics | Ariana Huffmyer, R. Cunning | M-multi-species scripts |
| Functional Annotation | Kathleen Durkin, Steven Roberts | 21.XX, 23.XX, 25.XX, 27.XX series |

## Key Datasets and Metadata

- **Sampling Dates**: [`SAMPLING_DATES.md`](SAMPLING_DATES.md) - Specific dates for timepoints TP1-TP4
- **RNA Metadata**: [`M-multi-species/data/rna_metadata.csv`](M-multi-species/data/rna_metadata.csv)
- **Physiological Data**: [E5 Timeseries Repository](https://github.com/urol-e5/timeseries/raw/refs/heads/master/time_series_analysis/Output/master_timeseries.csv)
- **Sample Information**: Detailed colony IDs, timepoints, and experimental design in metadata files

## Getting Started

1. **Clone the repository**
2. **Review species-specific READMEs** for detailed analysis workflows
3. **Check metadata files** to understand experimental design
4. **Follow numerical script order** within each species directory
5. **Ensure computational dependencies** are installed (R, Bioconductor packages, alignment tools)

For questions or contributions, please refer to individual script authors or open an issue in the repository.    
