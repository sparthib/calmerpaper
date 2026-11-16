#!/bin/bash
#SBATCH -p shared
#SBATCH --mem=90G
#SBATCH -c 4
#SBATCH --job-name=flagdiag_rna004
#SBATCH -o logs/flagdiag_rna004.txt
#SBATCH -e logs/flagdiag_rna004.txt
#SBATCH -t 4:00:00
# Run from figures/fig6_flagged_diagnostics/ :  sbatch compute_diagnostics_rna004.sh
mkdir -p logs
set +u
ml load conda_R/4.5.x
set -u
Rscript compute_diagnostics_rna004.R
