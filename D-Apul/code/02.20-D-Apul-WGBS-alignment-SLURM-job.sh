#!/bin/bash
#SBATCH --job-name=bismark_job_array
#SBATCH --account=coenv
#SBATCH --partition=cpu-g2
#SBATCH --output=bismark_job_%A_%a.out
#SBATCH --error=bismark_job_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=15
#SBATCH --mem=100G
#SBATCH --time=72:00:00
##turn on e-mail notification
#SBATCH --mail-type=ALL
#SBATCH --mail-user=samwhite@uw.edu
## Specify the working directory for this job
#SBATCH --chdir=/gscratch/srlab/sam/gitrepos/urol-e5/timeseries_molecular/D-Apul/output/02.20-D-Apul-WGBS-alignment-SLURM_array-bismark/

# Execute Roberts Lab bioinformatics container
# Binds home directory
# Binds /gscratch directory
# Directory bindings allow outputs to be written to the hard drive.

# Executes Bismark alignment using 02.01-bismark-bowtie2-alignment-SLURM-array.sh script.

# To execute this SLURM script as an array, start the script with the following command:

# sbatch --array=0-$(($$(wc -l < fastq_pairs.txt) - 1)) 02.02-bismark-SLURM-job.sh

# IMPORTANT: Requires fastq_pairs.txt to exist prior to submission!
apptainer exec \
--home "$PWD" \
--bind /mmfs1/home/ \
--bind /gscratch \
/gscratch/srlab/sr320/srlab-bioinformatics-container-586bf21.sif \
/gscratch/srlab/sam/gitrepos/urol-e5/timeseries_molecular/D-Apul/code/02.20-D-Apul-WGBS-alignment-SLURM_array-bismark.sh