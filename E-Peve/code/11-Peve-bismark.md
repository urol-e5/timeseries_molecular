### Background 

This code will align trimmed _P.evermanni_ EM-seq data to the _P.evermanni_ genome using Bismark and create loci counts matrices using the [nf-core Methylseq pipeline](https://nf-co.re/methylseq/3.0.0/).


**Input data:**

- [Trimmed fastqs](https://gannet.fish.washington.edu/gitrepos/urol-e5/timeseries_molecular/E-Peve/output/01.00-E-Peve-WGBS-trimming-fastp-FastQC-MultiQC)
- [genome](https://gannet.fish.washington.edu/seashell/snaps/Porites_evermanni_v1.fa)


**Files needed to run nextflow pipeline:**

- [config file](https://gannet.fish.washington.edu/metacarcinus/E5/Pevermanni/20250619_methylseq/uw_hyak_srlab.config)
- [samplesheet](https://gannet.fish.washington.edu/metacarcinus/E5/Pevermanni/20250619_methylseq/samplesheet.csv)

### Methods
**Copy genome files to klone**

```
# show path
pwd
/gscratch/srlab/strigg/GENOMES

# copy genome
wget https://gannet.fish.washington.edu/seashell/snaps/Porites_evermanni_v1.fa

```


**Copy EM-seq data to Klone**

```
# open screen session
screen -S methylseq

#make and change into dir

mkdir /gscratch/scrubbed/strigg/analyses/20250529_methylseq/data

cd /gscratch/scrubbed/strigg/analyses/20250529_methylseq/data

# copy data 
rsync --progress --verbose --archive shellytrigg@gannet.fish.washington.edu:/volume2/web/gitrepos/urol-e5/timeseries_molecular/E-Peve/output/01.00-E-Peve-WGBS-trimming-fastp-FastQC-MultiQC .

```

**Run nf-core Methylseq pipeline on Klone**

```
# open screen session 
screen -S methylseq

# start interactive node
salloc -A srlab -p cpu-g2-mem2x -N 1 -c 1 --mem=16GB --time=72:00:00

#activate conda environment
mamba activate nextflow

nextflow run nf-core/methylseq \
-c /gscratch/srlab/strigg/bin/uw_hyak_srlab.config \
--input /gscratch/scrubbed/strigg/analyses/20250619_methylseq/samplesheet.csv \
--outdir /gscratch/scrubbed/strigg/analyses/20250619_methylseq \
--fasta /gscratch/srlab/strigg/GENOMES/Porites_evermanni_v1.fa \
--em_seq \
-resume \
-with-report nf_report.html \
-with-trace \
-with-timeline nf_timeline.html \
--skip_trimming \
--nomeseq 

```

### Results

- Multiqc report: [https://gannet.fish.washington.edu/metacarcinus/E5/Pevermanni/20250619_methylseq/multiqc/bismark/multiqc_report.html](https://gannet.fish.washington.edu/metacarcinus/E5/Pevermanni/20250619_methylseq/multiqc/bismark/multiqc_report.html)
- Bismark summary report: [https://gannet.fish.washington.edu/metacarcinus/E5/Pevermanni/20250619_methylseq/bismark/summary/bismark_summary_report.html](https://gannet.fish.washington.edu/metacarcinus/E5/Pevermanni/20250619_methylseq/bismark/summary/bismark_summary_report.html)
- Pipeline report: [https://gannet.fish.washington.edu/metacarcinus/E5/Pevermanni/20250619_methylseq/nf_report.html](https://gannet.fish.washington.edu/metacarcinus/E5/Pevermanni/20250619_methylseq/nf_report.html)
- Pipeline timeline: [https://gannet.fish.washington.edu/metacarcinus/E5/Pevermanni/20250619_methylseq/nf_timeline.html](https://gannet.fish.washington.edu/metacarcinus/E5/Pevermanni/20250619_methylseq/nf_timeline.html)
- Counts matrices: [https://gannet.fish.washington.edu/metacarcinus/E5/Pevermanni/20250619_methylseq/bismark/methylation_calls/methylation_coverage/](https://gannet.fish.washington.edu/metacarcinus/E5/Pevermanni/20250619_methylseq/bismark/methylation_calls/methylation_coverage/)
	- <sample_name>.fastp-trim_bismark_bt2_pe.deduplicated.bismark.cov.gz
- Deduplicated sorted bam files: [https://gannet.fish.washington.edu/metacarcinus/E5/Pevermanni/20250619_methylseq/bismark/deduplicated/](https://gannet.fish.washington.edu/metacarcinus/E5/Pevermanni/20250619_methylseq/bismark/deduplicated/) 
	- <sample_name>.deduplicated.sorted.bam 
- Other bismark output: [https://gannet.fish.washington.edu/metacarcinus/E5/Pevermanni/20250619_methylseq/bismark/](https://gannet.fish.washington.edu/metacarcinus/E5/Pevermanni/20250619_methylseq/bismark/)
