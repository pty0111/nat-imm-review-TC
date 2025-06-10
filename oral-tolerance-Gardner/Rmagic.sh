#!/bin/bash
#SBATCH -J Rmagic
#SBATCH -p cpu
#SBATCH -N 1
#SBATCH -c 16
#SBATCH --mem=256G
#SBATCH -t 1-00:00:00
#SBATCH -o logs/Rmagic.stdout.txt
#SBATCH -e logs/Rmagic.stderr.txt
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=parkt1@mskcc.org

source ~/.bash_profile

mkdir -p logs

mamba activate multiome

echo "Running now"
echo $(date "+%m/%d/%Y-%H.%M.%S")
Rscript Rmagic.R
echo "done"
echo $(date "+%m/%d/%Y-%H.%M.%S")
