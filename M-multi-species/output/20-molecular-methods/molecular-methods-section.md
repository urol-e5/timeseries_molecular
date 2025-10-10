# Molecular Methods

## Overview

This document describes the comprehensive molecular methods used to develop genomic count matrices from sample extraction through quality-controlled count matrices for multi-omic analyses of three coral species (*Acropora pulchra*, *Porites evermanni*, and *Pocillopora tuahiniensis*). The workflow encompasses RNA-seq, small RNA-seq (miRNA), long non-coding RNA (lncRNA), and DNA methylation (WGBS) data generation and processing.

## Sample Collection and Preparation

### Experimental Design
- **Species**: Three coral species representing different life history strategies
  - *Acropora pulchra* (fast-growing branching coral)
  - *Porites evermanni* (slow-growing massive coral)
  - *Pocillopora tuahiniensis* (intermediate growth branching coral)
- **Timepoints**: Four sampling timepoints (TP1-TP4) across the experimental period
- **Biological Replication**: Multiple colonies per species per timepoint
- **Sample Types**: Flash-frozen tissue samples for RNA extraction and WGBS library preparation

### Nucleic Acid Extraction
Samples were processed for multiple molecular assays from the same biological material:
- **Total RNA extraction**: For RNA-seq and small RNA-seq library preparation
- **Genomic DNA extraction**: For whole genome bisulfite sequencing (WGBS)

## RNA-seq: Gene Expression Quantification

### Quality Control of Raw Reads
Initial quality assessment of raw RNA-seq reads was performed using:
- **FastQC** (v0.12.1): Individual sample quality metrics including per-base quality scores, sequence length distribution, GC content, adapter contamination, and overrepresented sequences
- **MultiQC**: Aggregate quality reports across all samples for comparative assessment

### Read Trimming and Quality Filtering
Sequencing reads were trimmed using **fastp** (Chen et al., 2023):
- **Adapter removal**: Automatic detection and removal of adapter sequences
- **Quality trimming**: 15bp trimmed from each read end based on initial quality assessment
- **Poly-G trimming**: Removal of poly-G tails common in two-color chemistry sequencing platforms
- **Length filtering**: Retention of reads meeting minimum length requirements
- **Output**: Trimmed paired-end FastQ files (`*fastp-trim.fq.gz`)

Post-trimming quality verification was performed using FastQC and MultiQC to confirm improvement in read quality metrics.

### Genome Indexing
Reference genome indexing was performed using **HISAT2** (Kim et al., 2019):
- **Input**: Species-specific reference genome FASTA files
- **Index generation**: Created HISAT2-formatted genome indexes optimized for splice-aware alignment
- **Species-specific indexes**: Separate indexes generated for each coral species

### Read Alignment
Trimmed paired-end reads were aligned to species-specific reference genomes using **HISAT2**:
- **Alignment parameters**: 
  - Splice-aware alignment for accurate intron-exon boundary detection
  - Paired-end mode maintaining mate pair information
  - Multi-threading for computational efficiency
- **Output**: Sequence Alignment Map (SAM) files containing alignment coordinates and quality scores

### SAM to BAM Conversion and Sorting
Alignment files were converted and processed using **SAMtools** (v1.12):
- **SAM to BAM conversion**: `samtools view` for binary compression and space efficiency
- **Sorting**: `samtools sort` to organize alignments by genomic coordinates
- **Indexing**: `samtools index` to create BAI index files enabling rapid random access
- **Output**: Sorted and indexed BAM files (`*.sorted.bam` and `*.sorted.bam.bai`)

### Transcript Assembly and Quantification
**StringTie** (v2.2.1) (Pertea et al., 2015, 2016) was used for transcript assembly and quantification:

#### Individual Sample Assembly
- **Input**: Sorted BAM files and reference genome annotation (GFF/GTF)
- **Parameters**: 
  - Reference-guided assembly mode
  - Multi-threading for efficiency
- **Output**: Sample-specific GTF files containing assembled transcripts

#### Transcript Merging
- **StringTie merge**: Combined individual sample GTF files into a unified transcript catalog
- **Reference integration**: Merged novel transcripts with known reference annotations
- **Output**: Consensus GTF file (`stringtie_merged.gtf`)

#### Count Matrix Generation
StringTie **prepDE.py** script generated count matrices for differential expression analysis:
- **Gene-level counts**: `apul-gene_count_matrix.csv` (or species-specific equivalents)
- **Transcript-level counts**: `apul-transcript_count_matrix.csv`
- **Matrix format**: Rows = genes/transcripts, Columns = samples, Values = read counts
- **Compatibility**: Formatted for use with DESeq2 and other R-based differential expression tools

## Small RNA-seq: miRNA Discovery and Quantification

### Quality Control and Trimming
Small RNA-seq data processing followed a specialized workflow optimized for short reads:

#### Initial Quality Assessment
- **FastQC**: Assessment of raw sRNA-seq reads (typically 15-30bp)
- **MultiQC**: Aggregated quality metrics across samples

#### Adapter Trimming
**fastp** was used with parameters optimized for small RNAs:
- **Adapter removal**: Aggressive adapter trimming for short inserts
- **Poly-G removal**: Trimming of poly-G tails
- **Length filtering**: Retention of 18-31bp reads typical of mature miRNAs
- **Quality threshold**: Minimum quality score filtering
- **Output**: Trimmed FastQ files (`*fastp-adapters-polyG-31bp-merged.fq.gz`)

### miRNA Discovery and Annotation
**ShortStack** (v4.1.0) (Axtell, 2013; Shahid & Axtell, 2014; Johnson et al., 2016) performed comprehensive small RNA analysis:

#### Reference Database
- **miRBase**: Customized cnidarian miRNA database (v22.1) curated by Jill Ashley
- **Inclusion criteria**: Published cnidarian miRNAs from multiple species
- **Database file**: `cnidarian-mirbase-mature-v22.1.fasta`

#### ShortStack Analysis
- **Alignment**: Short read alignment to species-specific reference genomes
- **Locus identification**: Detection of small RNA-producing loci
- **miRNA annotation**: Comparison against known miRNA database
- **Novel miRNA prediction**: De novo identification of miRNA candidates
- **Hairpin structure analysis**: Secondary structure prediction for miRNA precursors
- **Strand bias assessment**: Evaluation of 5p/3p arm expression patterns

#### Count Matrix Generation
ShortStack generated raw count matrices:
- **Format**: Tab-delimited text file with genomic coordinates and read counts per sample
- **Features**: miRNA loci, genomic coordinates, strand information, MIRNA classification column

### miRNA Count Matrix Formatting
Custom R script (`10-format-miRNA-counts.R`) processed raw ShortStack output:

#### Filtering Steps
1. **miRNA selection**: Retained only rows with "Y" in MIRNA column (confirmed miRNAs)
2. **Column removal**: Removed coordinate and classification columns (columns 1 and 3)
3. **Sample name extraction**: Extracted sample identifiers from column headers

#### Sample Naming Standardization
- **Metadata integration**: Matched sequencing sample names to biological metadata
- **Naming convention**: Reformatted to `AzentaSampleName_ColonyID_Timepoint`
- **Validation**: Verified accurate mapping between original and renamed samples

#### Output
- **Formatted count matrix**: Tab-delimited file (`*_miRNA_counts_formatted.txt`)
- **Row names**: miRNA identifiers
- **Column names**: Standardized sample identifiers
- **Values**: Raw read counts per miRNA per sample

## Long Non-coding RNA (lncRNA) Discovery and Quantification

lncRNA discovery employed a multi-step computational pipeline combining alignment, assembly, classification, and coding potential assessment.

### Read Alignment
Trimmed RNA-seq reads (same data as gene expression analysis) were aligned using **HISAT2**:
- **Strategy**: Splice-aware alignment to detect novel transcripts
- **Output**: Coordinate-sorted BAM files

### Transcript Assembly
**StringTie** assembled transcripts from alignments:
- **Mode**: Reference-guided assembly using genome annotation
- **Individual assemblies**: Generated per-sample GTF files
- **Merged assembly**: Combined all samples into unified transcript catalog

### Transcript Classification
**GFFcompare** (v0.12.6) classified assembled transcripts relative to reference annotation:

#### Classification Categories
Novel transcript identification focused on specific class codes:
- **Class code "u"**: Unknown intergenic transcripts (potential lncRNAs)
- **Class code "x"**: Exonic overlap with reference on opposite strand
- **Class code "o"**: Generic exonic overlap with reference
- **Class code "i"**: Intronic (fully contained within reference intron)

#### Length Filtering
- **Minimum length**: 200 nucleotides (standard lncRNA definition)
- **Rationale**: Distinguishes lncRNAs from small regulatory RNAs and degradation products

### Coding Potential Assessment
**Coding Potential Calculator 2 (CPC2)** (v1.0.1) predicted coding potential:

#### FASTA Extraction
**BEDtools getfasta** extracted sequences for candidate lncRNA transcripts:
- **Input**: Genome FASTA and filtered GTF coordinates
- **Parameters**: `-name -split` for proper handling of multi-exonic transcripts
- **Output**: Candidate lncRNA sequences in FASTA format

#### CPC2 Classification
- **Algorithm**: Machine learning-based coding potential prediction
- **Features**: ORF length, ORF coverage, Fickett score, isoelectric point
- **Classification**: Binary coding/noncoding prediction
- **Threshold**: Default CPC2 thresholds for noncoding classification

#### Noncoding Transcript Selection
- **Filtering**: Retained transcripts classified as "noncoding" by CPC2
- **Output**: List of confirmed lncRNA transcript IDs

### lncRNA Quantification
**featureCounts** (Subread v2.0.5) quantified lncRNA expression:

#### Count Generation Parameters
- **Annotation**: Custom lncRNA GTF file from discovery pipeline
- **Feature type**: `-t lncRNA` (custom feature designation)
- **Gene ID**: `-g gene_id` attribute for grouping
- **Paired-end mode**: `-p` flag for proper fragment counting
- **Multi-overlap handling**: `-O --fraction` for fractional counting of overlapping features
- **Multi-threading**: `-T 42` for parallel processing

#### Count Matrix Processing
Generated count matrices were processed through multiple cleaning steps:

1. **Column name cleaning** (using `awk`):
   - Removed full file paths from sample names
   - Stripped `.sorted.bam` suffix
   - Retained only base sample identifiers

2. **Quality filtering** (using custom Python script `11-filter_count_matrix.py`):
   - **Filtering criteria**: Removed lncRNAs with <10 reads in >50% of samples
   - **Rationale**: Eliminates low-confidence lncRNAs with insufficient coverage
   - **Output**: `lncRNA_counts.clean.filtered.txt`

#### Final Output
- **Clean count matrix**: Tab-delimited text file
- **Row format**: lncRNA gene IDs
- **Column format**: Sample identifiers
- **Values**: Raw read counts per lncRNA per sample
- **Quality**: Filtered for reliable quantification across samples

## DNA Methylation: Whole Genome Bisulfite Sequencing (WGBS)

### Quality Control of Bisulfite-Converted Reads
WGBS data quality assessment followed similar protocols to RNA-seq:
- **FastQC**: Evaluated raw bisulfite-converted sequencing reads
  - Note: C-to-T conversion artifacts expected in bisulfite-treated sequences
  - Base composition skew expected due to cytosine deamination
- **MultiQC**: Aggregated quality metrics across all WGBS samples

### Read Trimming
**fastp** trimmed bisulfite-converted reads:
- **Adapter removal**: Detection and removal of methylation adapter sequences
- **Quality filtering**: Trimmed low-quality bases from read ends
- **Output**: Trimmed paired-end FastQ files (`*fastp-trim.fq.gz`)

### Bisulfite Genome Preparation
**Bismark** (v0.24.0) prepared bisulfite-converted reference genomes:
- **Genome conversion**: Generated C-to-T and G-to-A converted reference genomes
- **Indexing**: Created Bowtie2 indexes (v2.4.4) for both converted strands
- **Parameters**: 
  - Multi-threading (`--parallel 28`)
  - Verbose output for tracking conversion progress
- **Output**: Bisulfite-converted genome indexes for alignment

### Bisulfite Read Alignment
**Bismark** aligned trimmed reads to bisulfite-converted genomes:

#### Alignment Strategy
- **Paired-end mode**: Maintained mate pair information for improved mapping accuracy
- **Aligner**: Bowtie2 backend for efficient alignment
- **Multi-threading**: Parallel processing (`-p 8`)
- **Strand-specific mapping**: Separate alignments to OT (original top), CTOT, OB, and CTOB strands

#### Score Threshold Optimization
Multiple alignment stringency parameters were tested:
- **Score_min values tested**:
  - `L,0,-0.4`
  - `L,0,-0.6`
  - `L,0,-0.8`
  - `L,0,-1.0`
  - `L,-1,-0.6`
- **Selection criteria**: Mapping efficiency and alignment quality metrics
- **Output**: BAM files with bisulfite alignment information (`*_bismark_bt2_pe.bam`)

### Deduplication
**Bismark deduplicate_bismark** removed PCR duplicates:
- **Strategy**: Identified and removed reads with identical mapping positions
- **Paired-end deduplication**: Considered both mate positions for duplicate determination
- **Output**: Deduplicated BAM files (`*deduplicated.bam`)

### BAM File Processing
**SAMtools** sorted and indexed deduplicated alignments:
- **Sorting**: Organized by genomic coordinates (`samtools sort --threads 24`)
- **Indexing**: Created BAI index files for rapid access
- **Output**: Sorted, indexed BAM files ready for methylation extraction

### Methylation Calling
**Bismark methylation_extractor** quantified methylation levels:

#### Extraction Parameters
- **Context**: Extracted methylation in CpG, CHG, and CHH contexts
- **Strand-specific**: Maintained original strand information
- **Coverage files**: Generated coverage files for each cytosine position
- **Bedgraph output**: Created genome-wide methylation tracks

#### Methylation Metrics
For each cytosine position, extracted:
- **Chromosome/Position**: Genomic coordinates
- **Strand**: DNA strand orientation (+/-)
- **Methylated count**: Number of reads supporting methylation
- **Unmethylated count**: Number of reads supporting no methylation
- **Context**: CpG, CHG, or CHH
- **Methylation percentage**: Calculated as (methylated / (methylated + unmethylated)) × 100

### Coverage Filtering and Count Matrix Generation
Methylation data were filtered and formatted for downstream analysis:

#### Coverage Thresholds
- **Minimum coverage**: Sites required minimum read depth for reliable quantification
- **Sample representation**: Sites retained if present in sufficient number of samples
- **Example threshold**: ≥20 reads coverage in sufficient samples (`merged-WGBS-CpG-counts_filtered_n20.csv`)

#### Count Matrix Format
- **Rows**: CpG sites (genomic coordinates)
- **Columns**: Samples
- **Values**: Methylation percentages or raw methylated/unmethylated counts
- **Output**: CSV format for compatibility with R-based analysis tools

## Quality Control and Filtering Summary

### General QC Principles Applied Across All Data Types

#### Pre-processing QC
1. **FastQC/MultiQC**: Comprehensive quality metrics visualization
2. **Adapter identification**: Detection of technical sequences for removal
3. **Quality score assessment**: Evaluation of per-base and per-read quality
4. **Length distribution**: Verification of expected fragment sizes

#### Post-processing QC
1. **Alignment statistics**: Mapping rates, duplicate rates, properly paired reads
2. **Coverage assessment**: Sequencing depth across features/genomic regions
3. **Count distribution**: Examination of count matrix distributions
4. **Sample correlation**: Verification of biological replicate similarity

### Data Type-Specific Filtering

#### RNA-seq and lncRNA
- **Low expression filtering**: Removal of genes/transcripts with insufficient reads
- **Sample coverage**: Required expression in minimum percentage of samples
- **Count threshold**: Typically ≥10 reads in ≥50% of samples

#### miRNA
- **miRNA classification**: Retained only confirmed miRNAs (MIRNA = "Y")
- **Database matching**: Filtered based on miRBase homology
- **Read count minimums**: Applied to remove spurious low-count features

#### DNA Methylation
- **Coverage thresholds**: Minimum read depth per CpG site (e.g., 10-20× coverage)
- **Sample representation**: Required coverage across minimum number of samples
- **Context-specific filtering**: Separate filters for CpG, CHG, CHH contexts

## Count Matrix Generation and Normalization Philosophy

### Raw Count Matrices
All pipelines generated raw, unnormalized count matrices preserving:
- **Integer counts**: Whole number read/fragment counts per feature
- **No pre-normalization**: Raw counts retained for statistical modeling
- **Complete feature set**: Initial matrices include all detected features before filtering

### Matrix Structure Standards
Consistent format across all data types:
- **Rows**: Genomic features (genes, transcripts, lncRNAs, miRNAs, CpG sites)
- **Columns**: Samples with standardized naming (SampleID_ColonyID_Timepoint)
- **First column**: Feature identifiers (gene IDs, miRNA names, genomic coordinates)
- **File formats**: Tab-delimited text (.txt) or comma-separated values (.csv)

### Downstream Normalization
Raw count matrices were prepared for various normalization approaches:
- **DESeq2**: Size factor normalization accounting for library size and composition
- **TMM**: Trimmed mean of M-values normalization (edgeR)
- **TPM/FPKM**: Transcript per million or fragments per kilobase per million (length normalization)
- **Methylation percentages**: Intrinsic normalization via (methylated/total) calculation

## Data Integration and Quality Assurance

### Sample Naming Standardization
Consistent sample identifiers were applied across all data types:
- **Format**: `SampleID_ColonyID_Timepoint`
- **Example**: `A001_Colony1_TP1`
- **Metadata integration**: Sample names linked to experimental metadata for downstream analysis
- **Cross-platform consistency**: Identical naming enables multi-omic integration

### Metadata Association
Comprehensive metadata tracked:
- **Sample information**: Species, colony ID, timepoint, extraction date
- **Sequencing information**: Run ID, lane, barcode, read length
- **Quality metrics**: Reads passed filter, alignment rate, library size
- **File provenance**: Input files, software versions, parameter settings

### Reproducibility and Documentation
All analyses maintained:
- **Script preservation**: Analysis code in version-controlled R Markdown/Quarto documents
- **Parameter documentation**: Complete parameter sets recorded in code chunks
- **File checksums**: MD5 checksums for all input and output files
- **Software versions**: Documented versions of all software tools used

## Software and Tool Versions

### Quality Control
- FastQC v0.12.1
- MultiQC (latest via conda/mamba)

### Trimming and Preprocessing
- fastp (Chen et al., 2023)

### Alignment and Assembly
- HISAT2 v2.2.1 (Kim et al., 2019)
- SAMtools v1.12
- StringTie v2.2.1 (Pertea et al., 2015, 2016)
- GFFcompare v0.12.6

### Feature Quantification
- featureCounts (Subread v2.0.5)

### Small RNA Analysis
- ShortStack v4.1.0 (Axtell, 2013; Shahid & Axtell, 2014; Johnson et al., 2016)

### lncRNA Discovery
- BEDtools v2 
- CPC2 v1.0.1 (Coding Potential Calculator 2)

### Bisulfite Sequencing
- Bismark v0.24.0
- Bowtie2 v2.4.4

### Data Processing
- R v4.x with packages:
  - dplyr, readr, stringr, tibble (tidyverse ecosystem)
- Python v3 with standard libraries

## References

Axtell, M. J. (2013). ShortStack: comprehensive annotation and quantification of small RNA genes. *RNA*, 19(6), 740-751.

Chen, S., Zhou, Y., Chen, Y., & Gu, J. (2018). fastp: an ultra-fast all-in-one FASTQ preprocessor. *Bioinformatics*, 34(17), i884-i890.

Johnson, N. R., Yeoh, J. M., Coruh, C., & Axtell, M. J. (2016). Improved placement of multi-mapping small RNAs. *G3: Genes, Genomes, Genetics*, 6(7), 2103-2111.

Kim, D., Paggi, J. M., Park, C., Bennett, C., & Salzberg, S. L. (2019). Graph-based genome alignment and genotyping with HISAT2 and HISAT-genotype. *Nature Biotechnology*, 37(8), 907-915.

Pertea, M., Pertea, G. M., Antonescu, C. M., Chang, T. C., Mendell, J. T., & Salzberg, S. L. (2015). StringTie enables improved reconstruction of a transcriptome from RNA-seq reads. *Nature Biotechnology*, 33(3), 290-295.

Pertea, M., Kim, D., Pertea, G. M., Leek, J. T., & Salzberg, S. L. (2016). Transcript-level expression analysis of RNA-seq experiments with HISAT, StringTie and Ballgown. *Nature Protocols*, 11(9), 1650-1667.

Shahid, S., & Axtell, M. J. (2014). Identification and annotation of small RNA genes using ShortStack. *Methods*, 67(1), 20-27.
