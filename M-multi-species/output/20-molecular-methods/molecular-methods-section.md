# Molecular Methods

## Sample Collection and Preparation

We collected samples from three coral species representing different life history strategies: *Acropora pulchra* (fast-growing branching coral), *Porites evermanni* (slow-growing massive coral), and *Pocillopora tuahiniensis* (intermediate growth branching coral). Samples were collected at four timepoints (TP1-TP4) across the experimental period with multiple colonies per species per timepoint to ensure biological replication. All samples were flash-frozen immediately upon collection for subsequent RNA extraction and whole genome bisulfite sequencing (WGBS) library preparation. From each biological sample, we extracted both total RNA for RNA-seq and small RNA-seq library preparation, as well as genomic DNA for WGBS analysis, enabling comprehensive multi-omic profiling from the same biological material

## RNA-seq Gene Expression Quantification

We performed initial quality assessment of raw RNA-seq reads using FastQC (v0.12.1) to evaluate per-base quality scores, sequence length distribution, GC content, adapter contamination, and overrepresented sequences. MultiQC was used to generate aggregate quality reports across all samples for comparative assessment. Sequencing reads were then trimmed using fastp (Chen et al., 2023) with automatic detection and removal of adapter sequences. Based on initial quality assessment, we trimmed 15bp from each read end and removed poly-G tails common in two-color chemistry sequencing platforms. Reads were filtered to retain only those meeting minimum length and quality requirements, producing trimmed paired-end FastQ files. Post-trimming quality verification was performed using FastQC and MultiQC to confirm improvement in read quality metrics.

For read alignment, we first generated species-specific reference genome indexes using HISAT2 (Kim et al., 2019), creating HISAT2-formatted genome indexes optimized for splice-aware alignment from species-specific reference genome FASTA files. Trimmed paired-end reads were then aligned to these reference genomes using HISAT2 with splice-aware alignment for accurate intron-exon boundary detection, maintaining mate pair information in paired-end mode while using multi-threading for computational efficiency. The resulting Sequence Alignment Map (SAM) files containing alignment coordinates and quality scores were converted to binary BAM format using SAMtools (v1.12) for space efficiency. BAM files were then sorted by genomic coordinates and indexed to create BAI index files enabling rapid random access.

We used StringTie (v2.2.1) (Pertea et al., 2015, 2016) for transcript assembly and quantification. For each sample, sorted BAM files and reference genome annotations (GFF/GTF) were used as input for reference-guided assembly, generating sample-specific GTF files containing assembled transcripts. Individual sample GTF files were then merged using StringTie merge to create a unified transcript catalog, integrating novel transcripts with known reference annotations in a consensus GTF file. Finally, we generated count matrices for differential expression analysis using the StringTie prepDE.py script, producing both gene-level and transcript-level count matrices in CSV format with rows representing genes or transcripts, columns representing samples, and values representing raw read counts. These matrices were formatted for compatibility with DESeq2 and other R-based differential expression analysis tools.

## Small RNA-seq miRNA Discovery and Quantification

Small RNA-seq data processing followed a specialized workflow optimized for short reads. We assessed raw sRNA-seq reads (typically 15-30bp) using FastQC, with MultiQC providing aggregated quality metrics across samples. Adapter trimming was performed using fastp with parameters optimized for small RNAs, including aggressive adapter trimming for short inserts, removal of poly-G tails, retention of 18-31bp reads typical of mature miRNAs, and minimum quality score filtering.

For miRNA discovery and annotation, we used ShortStack (v4.1.0) (Axtell, 2013; Shahid & Axtell, 2014; Johnson et al., 2016) with a customized cnidarian miRNA database (miRBase v22.1) curated by Jill Ashley, which includes published cnidarian miRNAs from multiple species. ShortStack performed short read alignment to species-specific reference genomes, detected small RNA-producing loci, compared them against the known miRNA database for annotation, and predicted novel miRNA candidates de novo. The analysis included secondary structure prediction for miRNA precursor hairpins and evaluation of 5p/3p arm expression patterns to assess strand bias.

ShortStack generated raw count matrices as tab-delimited text files containing genomic coordinates and read counts per sample, with features including miRNA loci, genomic coordinates, strand information, and a MIRNA classification column. We processed these raw outputs using a custom R script to filter for confirmed miRNAs (retaining only rows with "Y" in the MIRNA column), remove coordinate and classification columns, and extract sample identifiers from column headers. Sample names were standardized by matching sequencing sample identifiers to biological metadata and reformatting to the convention AzentaSampleName_ColonyID_Timepoint. We validated accurate mapping between original and renamed samples before generating the final formatted count matrix with miRNA identifiers as row names, standardized sample identifiers as column names, and raw read counts as values.

## Long Non-coding RNA Discovery and Quantification

Long non-coding RNA (lncRNA) discovery employed a multi-step computational pipeline combining alignment, assembly, classification, and coding potential assessment. Using the same trimmed RNA-seq reads as the gene expression analysis, we performed splice-aware alignment to reference genomes using HISAT2 to detect novel transcripts, producing coordinate-sorted BAM files. We then used StringTie for reference-guided transcript assembly using genome annotations, generating per-sample GTF files which were subsequently combined into a unified transcript catalog.

To classify assembled transcripts relative to reference annotations, we used GFFcompare (v0.12.6), focusing on novel transcripts with specific class codes: "u" (unknown intergenic transcripts representing potential lncRNAs), "x" (exonic overlap with reference on opposite strand), "o" (generic exonic overlap with reference), and "i" (intronic, fully contained within reference intron). We applied a minimum length filter of 200 nucleotides, consistent with the standard lncRNA definition, to distinguish lncRNAs from small regulatory RNAs and degradation products.

To assess coding potential of candidate lncRNAs, we first extracted sequences for candidate transcripts using BEDtools getfasta with the genome FASTA and filtered GTF coordinates, using the -name and -split parameters for proper handling of multi-exonic transcripts. We then used Coding Potential Calculator 2 (CPC2 v1.0.1), which employs machine learning-based prediction using features including ORF length, ORF coverage, Fickett score, and isoelectric point to provide binary coding/noncoding classification. Transcripts classified as "noncoding" by CPC2 using default thresholds were retained as confirmed lncRNAs.

We quantified lncRNA expression using featureCounts (Subread v2.0.5) with our custom lncRNA GTF file from the discovery pipeline. Count generation used the lncRNA feature type designation with gene_id attribute for grouping, paired-end mode for proper fragment counting, and fractional counting for overlapping features (-O --fraction flags), with multi-threading for parallel processing. Generated count matrices were processed through multiple cleaning steps: column names were cleaned using awk to remove full file paths and .sorted.bam suffixes, retaining only base sample identifiers. We then applied quality filtering using a custom Python script (11-filter_count_matrix.py) to remove lncRNAs with fewer than 10 reads in more than 50% of samples, eliminating low-confidence lncRNAs with insufficient coverage. The final clean count matrix was formatted as a tab-delimited text file with lncRNA gene IDs as rows, sample identifiers as columns, and raw read counts per lncRNA per sample as values.

## Whole Genome Bisulfite Sequencing DNA Methylation Analysis

We performed quality assessment of raw WGBS data using FastQC and MultiQC, noting that base composition skew and apparent C-to-T conversion artifacts are expected in bisulfite-treated sequences due to cytosine deamination. Reads were trimmed using fastp with detection and removal of methylation adapter sequences and trimming of low-quality bases from read ends. 

For alignment, we first prepared bisulfite-converted reference genomes using Bismark (v0.24.0), which generated both C-to-T and G-to-A converted reference genomes and created Bowtie2 (v2.4.4) indexes for both converted strands using multi-threading (--parallel 28) with verbose output for tracking conversion progress. Trimmed reads were aligned to the bisulfite-converted genomes using Bismark in paired-end mode with the Bowtie2 backend, maintaining mate pair information for improved mapping accuracy while using parallel processing (-p 8) and performing strand-specific mapping to separate alignments to OT (original top), CTOT, OB, and CTOB strands. We tested multiple alignment stringency parameters (score_min values of L,0,-0.4; L,0,-0.6; L,0,-0.8; L,0,-1.0; and L,-1,-0.6) to optimize mapping efficiency and alignment quality metrics.

Following alignment, we removed PCR duplicates using Bismark deduplicate_bismark, which identifies and removes reads with identical mapping positions while considering both mate positions for duplicate determination in paired-end mode. Deduplicated BAM files were then sorted by genomic coordinates and indexed using SAMtools (--threads 24) to create BAI index files for rapid access.

We quantified methylation levels using Bismark methylation_extractor, which extracted methylation in CpG, CHG, and CHH contexts while maintaining original strand information and generating coverage files for each cytosine position along with genome-wide methylation tracks in bedgraph format. For each cytosine position, we extracted chromosome and position coordinates, strand orientation, counts of reads supporting methylation and unmethylation, sequence context (CpG, CHG, or CHH), and calculated methylation percentage as (methylated reads / (methylated + unmethylated reads)) × 100.

Methylation data were filtered and formatted for downstream analysis by applying coverage thresholds requiring minimum read depth for reliable quantification at each site and minimum representation across samples. For example, sites were required to have at least 20 reads coverage in a sufficient number of samples. The final count matrices were formatted as CSV files with rows representing CpG sites (genomic coordinates), columns representing samples, and values representing methylation percentages or raw methylated/unmethylated counts, ensuring compatibility with R-based analysis tools.

## Quality Control and Data Processing Standards

Across all molecular data types, we applied consistent quality control principles throughout the analytical workflow. Pre-processing quality control included comprehensive quality metrics visualization using FastQC and MultiQC, detection of technical sequences for removal, evaluation of per-base and per-read quality scores, and verification of expected fragment size distributions. Post-processing quality control assessed alignment statistics including mapping rates, duplicate rates, and properly paired read percentages, evaluated sequencing depth across features and genomic regions, examined count matrix distributions, and verified biological replicate similarity through sample correlation analysis.

Data type-specific filtering criteria were applied to ensure data quality. For RNA-seq and lncRNA data, we removed genes and transcripts with insufficient expression, requiring a minimum of 10 reads in at least 50% of samples to ensure reliable quantification. For miRNA data, we retained only confirmed miRNAs based on ShortStack classification (MIRNA = "Y") and miRBase homology, with additional read count minimums applied to remove spurious low-count features. For DNA methylation data, we required minimum read depth per CpG site (typically 10-20× coverage), adequate coverage across a minimum number of samples, and applied separate filters for CpG, CHG, and CHH sequence contexts.

All analytical pipelines generated raw, unnormalized count matrices that preserved integer read or fragment counts per feature without pre-normalization, retaining complete feature sets before filtering to enable appropriate statistical modeling. Count matrices followed a consistent structure across all data types: rows represented genomic features (genes, transcripts, lncRNAs, miRNAs, or CpG sites), columns represented samples with standardized naming following the format SampleID_ColonyID_Timepoint, the first column contained feature identifiers (gene IDs, miRNA names, or genomic coordinates), and files were saved in tab-delimited text (.txt) or comma-separated values (.csv) formats. These raw count matrices were prepared for various downstream normalization approaches including DESeq2 size factor normalization accounting for library size and composition, trimmed mean of M-values (TMM) normalization in edgeR, transcripts per million (TPM) or fragments per kilobase per million (FPKM) for length normalization, and methylation percentages which provide intrinsic normalization via (methylated/total) calculation.

To enable multi-omic integration, we applied consistent sample identifiers across all data types using the format SampleID_ColonyID_Timepoint (e.g., A001_Colony1_TP1). Sample names were systematically linked to experimental metadata tracking sample information (species, colony ID, timepoint, extraction date), sequencing information (run ID, lane, barcode, read length), quality metrics (reads passed filter, alignment rate, library size), and file provenance (input files, software versions, parameter settings). All analyses maintained reproducibility through preservation of analysis code in version-controlled R Markdown and Quarto documents, complete documentation of parameter sets in code chunks, generation of MD5 checksums for all input and output files, and thorough documentation of software versions for all tools used.

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
