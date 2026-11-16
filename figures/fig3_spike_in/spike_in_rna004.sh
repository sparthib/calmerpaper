#!/bin/bash

#SBATCH -p shared
#SBATCH --mem=64G
#SBATCH -c 4
#SBATCH --job-name=fig3_rna004
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/fig3_rna004.%j.txt
#SBATCH -e logs/fig3_rna004.%j.txt
#SBATCH -t 4-00:00:00

mkdir -p logs

echo "**** Job starts ****"
date +"%Y-%m-%d %T"
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"

module load conda_R/4.5.x

# Pooled Sequins Mix A figure; loads all MAPQ-filter oarfish quant internally.
Rscript spike_in_rna004.R

echo "**** Job ends ****"
date +"%Y-%m-%d %T"
