#!/bin/bash

#SBATCH -p shared
#SBATCH --mem=64G
#SBATCH -c 4
#SBATCH --job-name=mapq_vs_tx_pcr
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o logs/mapq_vs_tx_pcr.txt
#SBATCH -e logs/mapq_vs_tx_pcr.txt
#SBATCH -t 4:00:00

echo "**** Job starts ****"
date +"%Y-%m-%d %T"

mkdir -p logs

# conda_R activate scripts reference unset env vars; relax -u while loading.
set +u
ml load conda_R/4.5.x
set -u

Rscript mapq_vs_tx_features_pcr_kits.R

echo "**** Job ends ****"
date +"%Y-%m-%d %T"
