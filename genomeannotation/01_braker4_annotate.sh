#!/usr/bin/env bash
#1: Structural annotation with BRAKER4 (ETP mode: GeneMark-ETP + AUGUSTUS)
GENOME="~/CricketBaseCode/genomeannotation/data/genome.fa.masked" # Path to the soft-masked genome FASTA file
#Arthropoda from OrthoDBv12 (https://bioinf.uni-greifswald.de/bioinf/partitioned_odb12/Arthropoda.fa.gz)
PROTEINS="~/CricketBaseCode/genomeannotation/data/Arthropoda.fa" #gunzip first 

# (_1.fastq.gz, _2.fastq.gz) pairs 
#7 paired-end RNA-seq libraries (SRR14026720 ... SRR14026726) from different developmental stages 
# 9 DRR (DRR358356 ... DRR358364) from different tissue 
CSV="~/CricketBaseCode/genomeannotation/samples.csv"
CONFIG="~/CricketBaseCode/genomeannotation/config.ini" 

source ~/miniconda3/etc/profile.d/conda.sh
conda activate braker4

cd ~/CricketBaseCode/genomeannotation
export BRAKER4_CONFIG=$PWD/config.ini
export SINGULARITYENV_PREPEND_PATH=/opt/conda/bin

snakemake \
    --snakefile ~/tools/BRAKER4/Snakefile \
    --cores 16 --jobs 200 \
    --printshellcmds \
    --rerun-incomplete \
    --latency-wait 120 \
    --restart-times 3 \
    --use-singularity \
    --singularity-prefix ~/CricketBaseCode/genomeannotation/.singularity_cache \
    --singularity-args "-B /lustre12 -B /home --env PREPEND_PATH=/opt/conda/bin" \
    --executor slurm \
    --default-resources slurm_partition=l40s slurm_account=l40s mem_mb=90000 \
    --keep-going








