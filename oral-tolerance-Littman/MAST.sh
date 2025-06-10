#!/bin/bash
#SBATCH -J MAST
#SBATCH -p cpu
#SBATCH -N 1
#SBATCH -c 16
#SBATCH --mem=128G
#SBATCH -t 2-00:00:00
#SBATCH -o logs/MAST.stdout.txt
#SBATCH -e logs/MAST.stderr.txt
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=parkt1@mskcc.org

source ~/.bash_profile

mkdir -p logs

mamba activate multiome

echo $(date "+%m/%d/%Y-%H.%M.%S")
Rscript MAST.R
echo $(date "+%m/%d/%Y-%H.%M.%S")

sstat -j $SLURM_JOB_ID.batch --format=JobID,MaxRss%8,MinCPU,MaxDiskWrite,nodelist%8