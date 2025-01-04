#!/bin/bash

# This script is designed to be called by a SLURM script which
# runs this script across an array of HPC nodes.

### IMPORTANT ###

# INPUT FILES
repo_dir="/gscratch/srlab/sam/gitrepos/urol-e5/timeseries_molecular"
trimmed_fastqs_dir="${repo_dir}/D-Apul/output/01.20-D-Apul-WGBS-trimming-fastp-FastQC-MultiQC"
bisulfite_genome_dir="${repo_dir}/D-Apul/data"

# OUTPUT FILES
output_dir_top="${repo_dir}/D-Apul/output/02.20-D-Apul-WGBS-alignment-SLURM_array-bismark"

# PARAMETERS
bowtie2_min_score="L,0,-0.6"

# CPU threads
# Bismark already spawns multiple instances and additional threads are multiplicative."
bismark_threads=5

###################################################################################


## SET ARRAY TASKS ##
cd "${output_dir_top}"

# Get the FastQ file pair for this task
# the `p` sets the line number to process
# which corresponds to the array task ID
pair=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${output_dir_top}/fastq_pairs.txt")

echo "Contents of pair:"
echo "${pair}"
echo ""

R1_file=$(echo $pair | awk '{print $1}')
R2_file=$(echo $pair | awk '{print $2}')

# Get just the sample name (excludes the _R[12]_001*)
sample_name=$(echo "$R1_file" | awk -F"_" '{print $1}')

# Check if R1_file and R2_file are not empty
if [ -z "$R1_file" ] || [ -z "$R2_file" ]; then
  echo "Error: R1_file or R2_file is empty. Exiting."
  exit 1
fi

# Check if sample_name is not empty
if [ -z "$sample_name" ]; then
  echo "Error: sample_name is empty. Exiting."
  exit 1
fi

echo "Contents of sample_name: ${sample_name}"
echo ""


## RUN BISMARK ALIGNMENTS ##
bismark \
--genome ${bisulfite_genome_dir} \
--score_min "${bowtie2_min_score}" \
--parallel "${bismark_threads}" \
--non_directional \
--gzip \
-p "${bismark_threads}" \
-1 ${trimmed_fastqs_dir}/${R1_file} \
-2 ${trimmed_fastqs_dir}/${R2_file} \
--output_dir "${output_dir_top}" \
2> "${output_dir_top}"/${sample_name}-${SLURM_ARRAY_TASK_ID}-bismark_summary.txt
