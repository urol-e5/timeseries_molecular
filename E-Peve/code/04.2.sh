#!/bin/bash
# Set directories and files
reads_dir="../../data/03-Peve-bismark/"
genome_folder="../../data/"
output_dir="."
checkpoint_file="completed_samples.log"

#

# Run Bismark for this sample
# bismark \
#     --genome /gscratch/srlab/sr320/github/timeseries_molecular/E-Peve/data \
#     -p 8 \
#     -u 10000 \
#     -score_min L,0,-0.6 \
#     -1 /gscratch/srlab/sr320/github/timeseries_molecular/E-Peve/data/03-Peve-bismark/POR-69-TP1_R1_001.fastp-trim.fq.gz \
#     -2 /gscratch/srlab/sr320/github/timeseries_molecular/E-Peve/data/03-Peve-bismark/POR-69-TP1_R2_001.fastp-trim.fq.gz \
#     -o /gscratch/srlab/sr320/github/timeseries_molecular/E-Peve/output/04-Peve-bismark-array
    
bismark \
    --genome ../../data \
    --multicore 8 \
    -u 10000 \
    --score_min L,0,-0.6 \
    -1 ../../data/03-Peve-bismark/POR-69-TP1_R1_001.fastp-trim.fq.gz \
    -2 /mmfs1/gscratch/srlab/sr320/github/timeseries_molecular/E-Peve/data/03-Peve-bismark/POR-69-TP1_R2_001.fastp-trim.fq.gz \
    --output_dir /gscratch/srlab/sr320/github/timeseries_molecular/E-Peve/output/04-Peve-bismark-array