#!/bin/bash
#SBATCH -p shared
#SBATCH --mem=48G
#SBATCH -c 2
#SBATCH --job-name=fig_sup4
#SBATCH -o logs/fig_sup4.txt
#SBATCH -e logs/fig_sup4.txt
#SBATCH -t 1:00:00
# Run from figures/fig_sup4/ :  sbatch predictor_distributions.sh
mkdir -p logs
set +u
ml load conda_R/4.5.x
set -u
Rscript predictor_distributions.R
