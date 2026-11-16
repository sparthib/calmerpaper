#!/bin/bash

#SBATCH -p shared
#SBATCH --mem=64G
#SBATCH -c 4
#SBATCH --job-name=fig2_rna004
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o logs/fig2_rna004.txt
#SBATCH -e logs/fig2_rna004.txt
#SBATCH -t 4:00:00

echo "**** Job starts ****"
date +"%Y-%m-%d %T"

mkdir -p logs

# conda_R activate scripts reference unset env vars; relax -u while loading.
set +u
ml load conda_R/4.5.x
set -u

Rscript interactions_rna004.R

echo "**** Job ends ****"
date +"%Y-%m-%d %T"
