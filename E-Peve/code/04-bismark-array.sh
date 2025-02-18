#!/bin/sh

#SBATCH --job-name=bismark_array          # Job name
#SBATCH --output=%x_%A_%a.slurm.out             # Standard output and error log
#SBATCH --error=%x_%A_%a.slurm.err              # Error log
#SBATCH --account=srlab
#SBATCH --partition=ckpt #update this line - use hyakalloc to find partitions you can use
#SBATCH --time=07-02:00:00
#SBATCH --ntasks=1                        # Run a single task
#SBATCH --cpus-per-task=24                # Number of CPU cores per task
#SBATCH --array=0-35                      # Array range (adjust based on the number of samples)
#SBATCH --mem=100G                         # Memory per node
#SBATCH --chdir=/gscratch/srlab/sr320/github/timeseries_molecular/E-Peve/output/04-Peve-bismark-array



# Execute Roberts Lab bioinformatics container
# Binds home directory
# Binds /gscratch directory
# Directory bindings allow outputs to be written to the hard drive.
apptainer exec \
--home "$PWD" \
--bind /mmfs1/home/ \
--bind /mmfs1/gscratch/ \
/gscratch/srlab/sr320/srlab-bioinformatics-container-586bf21.sif \
../../code/04.2.sh
