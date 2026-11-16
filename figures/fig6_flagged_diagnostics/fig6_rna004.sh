#!/bin/bash
#SBATCH -p shared
#SBATCH --mem=32G
#SBATCH -c 2
#SBATCH --job-name=fig6_rna004
#SBATCH -o logs/fig6_rna004.txt
#SBATCH -e logs/fig6_rna004.txt
#SBATCH -t 1:00:00
# Run from figures/fig6_flagged_diagnostics/ AFTER compute_diagnostics_rna004.sh finishes:
#   sbatch fig6_rna004.sh
mkdir -p logs
set +u
ml load conda_R/4.5.x
set -u
Rscript flagged_diagnostics_figure_rna004.R
