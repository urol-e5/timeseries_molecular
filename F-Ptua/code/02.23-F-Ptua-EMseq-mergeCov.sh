#!/bin/bash
## Job Name
#SBATCH --job-name=20250821_mergeCov
## Allocation Definition
#SBATCH --account=coenv
#SBATCH --partition=cpu-g2
## Resources
## Nodes
#SBATCH --nodes=1
## Walltime (days-hours:minutes:seconds format)
#SBATCH --time=1-00:00:00
## Memory per node
#SBATCH --mem=300G
##turn on e-mail notification
#SBATCH --mail-type=ALL
#SBATCH --mail-user=strigg@uw.edu
## Specify the working directory for this job 
#SBATCH --chdir=/gscratch/scrubbed/strigg/analyses/20250821_meth_Ptua

#ran this rsync command first to copy over the merged .cov.gz files to the outdir
rsync --archive --verbose --progress shellytrigg@gannet.fish.washington.edu:/volume2/web/metacarcinus/E5/Ptuahiniensis/20250422_methylseq/bismark/coverage2cytosine/coverage/ .

%%bash

set -ex


# make bed file from cov file keeping only CpGs w. 10x cov
for f in *.CpG_report.merged_CpG_evidence.cov.gz
do
  STEM=$(basename "${f}") # Get the entire filename including the long suffix
  STEM="${STEM%*.CpG_report.merged_CpG_evidence.cov.gz}" # Remove the suffix using parameter expansion
  zcat "${f}" | awk -F $'\t' 'BEGIN {OFS = FS} {if ($5+$6 >= 10) {print $1, $2, $3, $4}}' \
  > "${STEM}"_10x.bedgraph
done

# create modified tables with two columns; one is the CpG ID which is merged chrom and start site; one is the %meth 
for file in *10x.bedgraph; do
    awk -F"\t" -v fname="${file%_10x*}" 'BEGIN {print "CpG\t" fname}{print "CpG_"_$1"_"$2"\t"$4}' "$file" > "${file%.bedgraph}_processed.txt"
done

python /gscratch/srlab/strigg/scripts/merge_processed_txt.py
