`urol-e5/timeseries_molecular/D-Apul/code`

Code directory for _A.pulchra_.

All file names should adhere to this format:

00.00-<species_designation><code_topic>.Rmd

- `00.00`: Numbers to left of decimal should be inremented for each new type of analysis, which is unrelated to other existing analyses. Numbers to the right of the decimal should be incremented for analyses related to existing analyses.

	- E.g. `01.00-exon-counts.Rmd`
	- E.g. `01.01-exon-counts-per-gene.Rmd`
	
- `<species_designation>`: Should match top level of repo. E.g. `D-Apul`.

- `<code_topic>`: Provide a bried description of code functionality. E.g. `exon_counts`.

- `.Rmd`: Most code is anticipated to be run using R Markdown. However, other formats are fine (e.g. Jupyter Notebooks). 

All outputs from a script should be directed to a corresponding directory with the same name (excluding the suffix) in the `../ouput` directory.

---


- [`00.00-D-Apul-RNAseq-reads-FastQC-MultiQC.Rmd`](./00.00-D-Apul-RNAseq-reads-FastQC-MultiQC.Rmd): Quality check of raw RNA-seq reads using FastQC and MultiQC.

[`01.00-D-Apul-RNAseq-trimming-fastp-FastQC-MultiQC.Rmd`](./01.00-D-Apul-RNAseq-trimming-fastp-FastQC-MultiQC.Rmd): Quality trimming and adapter removal of RNA-seq reads using [fastp](https://github.com/OpenGene/fastp), followed by quality checks with [FastQC](https://github.com/s-andrews/FastQC) and [MultiQC](https://github.com/MultiQC/MultiQC).

[`02.00-D-Apul-RNAseq-gff-to-gtf.Rmd`](./02.00-D-Apul-RNAseq-gff-to-gtf.Rmd): Generate a GTF file from the _A.pulchra_ genome GFF, as the GTF is needed for downstream analyses.

[`02.10-D-Apul-RNAseq-genome-index-HiSat2.Rmd`](./02.10-D-Apul-RNAseq-genome-index-HiSat2.Rmd): Generate a genome index file with exons and splice sites, using [HISAT2](https://daehwankimlab.github.io/hisat2/manual/) for subsequent RNA-seq reads alignments with [HISAT2](https://daehwankimlab.github.io/hisat2/manual/).