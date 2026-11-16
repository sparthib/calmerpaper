#!/bin/bash
#SBATCH -p shared
#SBATCH --mem=8G
#SBATCH -c 1
#SBATCH --job-name=fig5_comb
#SBATCH -o logs/fig5_combined.txt
#SBATCH -e logs/fig5_combined.txt
#SBATCH -t 0:20:00
# Run from figures/fig5_partial_effects/ :  sbatch forest_plot_combined.sh
mkdir -p logs
set +u
ml load conda_R/4.5.x
set -u
Rscript forest_plot_combined.R
