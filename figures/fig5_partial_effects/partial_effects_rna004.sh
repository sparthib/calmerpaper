#!/bin/bash
#SBATCH --job-name=fig5p_rna004
#SBATCH --partition=shared
#SBATCH --mem=40G
#SBATCH --cpus-per-task=4
#SBATCH --time=04:00:00
#SBATCH --output=logs/fig5_partial_effects.%j.txt
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=sparthi1@jhu.edu

# NB: no `set -u` — the conda_R/4.5.x activation script references unbound vars
# (e.g. ADDR2LINE) and would abort module load under nounset.
module load conda_R/4.5.x
Rscript partial_effects_rna004.R "${1:-SRR30947490}"
