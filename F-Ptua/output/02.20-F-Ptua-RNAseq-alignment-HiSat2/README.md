`urol-e5/timeseries_molecular/F-Ptua/output/02.20-F-Ptua-RNAseq-alignment-HiSat2`

Output directory from [`02.20-F-Ptua-RNAseq-alignment-HiSat2.Rmd`](../../code/02.20-F-Ptua-RNAseq-alignment-HiSat2.Rmd).

NOTE: Due to large file sizes of some output files (e.g. BAMs), they are not actually present in this repo.

They are accessible here:

- https://gannet.fish.washington.edu/gitrepos/urol-e5/timeseries_molecular/F-Ptua/output/02.20-F-Ptua-RNAseq-alignment-HiSat2

---
  
  - `checksums.md5`: MD5 checksum for all files in this directory. Excludes subdirectories.
  
  - `ptua-gene_count_matrix.csv`: Gene count matrix for use in [DESeq2](https://github.com/thelovelab/DESeq2).
  
  - `ptua-transcript_count_matrix.csv`: Transcript count matrix for use in [DESeq2](https://github.com/thelovelab/DESeq2). 
  
  - `prepDE-sample_list.txt`: Sample file list provided as input to [StringTie](https://ccb.jhu.edu/software/stringtie/index.shtml?t=manual) for [DESeq2](https://github.com/thelovelab/DESeq2) count matrix generation. Also serves as documentation of which files were used for this step. 
  
  - `Pocillopora_meandrina_HIv1.assembly.stringtie.gtf`: Canonical [StringTie](https://ccb.jhu.edu/software/stringtie/index.shtml?t=manual) GTF file compiled from all individual sample GTFs. 
  
  - `sorted-bams-merged.bam`: Merged (and sorted) BAM consisting of all individual sample BAMs. 
  
  - `sorted-bams-merged.bam.bai`: BAM index file. Useful for visualizing assemblies in IGV. 
  
  - `sorted_bams.list`: List file needed for merging of BAMS with samtools. Also serves as documentation of which files were used for this step. 
  
  - `multiqc_report.html`: [MultiQC](https://github.com/MultiQC/MultiQC) report aggregating all individual [HISAT2](https://daehwankimlab.github.io/hisat2/manual/) alignment stats and samtools flagstats. 
  
  
  - `gtf_list.txt`: List file needed for merging of GTF files with [StringTie](https://ccb.jhu.edu/software/stringtie/index.shtml?t=manual). Also serves as documentation of which files were used for this step. 


## Subdirectories

Each subdirectory is labelled based on sample name and each contains individual [HISAT2](https://daehwankimlab.github.io/hisat2/manual/) alignment and [StringTie](https://ccb.jhu.edu/software/stringtie/index.shtml?t=manual) output files:
  
  - `<sample_name>_checksums.md5`: MD5 checksums for all files in the directory. 
  
  - `*.ctab`: Data tables formatted for import into https://github.com/alyssafrazee/ballgown. 
  
  - `<sample_name>.cov_refs.gtf`: [StringTie](https://ccb.jhu.edu/software/stringtie/index.shtml?t=manual) genome reference sequnce coverage GTF. 
  
  - `<sample_name>.gtf`: [StringTie](https://ccb.jhu.edu/software/stringtie/index.shtml?t=manual) GTF. 
  
  - `<sample_name>.sorted.bam`: [HISAT2](https://daehwankimlab.github.io/hisat2/manual/) assembly BAM. 
  
  - `<sample_name>.sorted.bam.bai`: BAM index file. Useful for visualizing assemblies in [IGV](https://igv.org/). 
  
  - `<sample_name>-hisat2_output.flagstat`: samtools flagstat output file. 
  
  - `<sample_name>_hisat2.stats`: [HISAT2](https://daehwankimlab.github.io/hisat2/manual/) assembly stats. 
  
  - `input_fastqs_checksums.md5`: MD5 checksums of files used as input for assembly. Primarily serves as documentation to track/verify which files were actually used.
