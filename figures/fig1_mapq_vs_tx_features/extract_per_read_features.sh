#!/bin/bash

#SBATCH -p shared
#SBATCH --mem=50G
#SBATCH -c 4
#SBATCH --job-name=per_read_feat
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH --output=logs/per_read_feat.%a.txt
#SBATCH --error=logs/per_read_feat.%a.txt
#SBATCH --array=1-113
#SBATCH -t 6:00:00

mkdir -p logs

echo "**** Job starts ****"
date +"%Y-%m-%d %T"
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

ml load conda_R/4.5.x

Rscript extract_per_read_features.R

echo "**** Job ends ****"
date +"%Y-%m-%d %T"
