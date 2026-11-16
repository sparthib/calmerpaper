#!/bin/bash
#SBATCH --partition=shared
#SBATCH --cpus-per-task=4
#SBATCH --mem=40G
#SBATCH --time=04:00:00
#SBATCH --job-name=fig5f_rna004
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH -o /dcs04/hicks/data/sparthib/lrs_coverage_project/lrs_cov/figures/fig5_partial_effects/forest_plot.%j.txt
set -eo pipefail
module load conda_R/4.5.x
cd /dcs04/hicks/data/sparthib/lrs_coverage_project/lrs_cov/figures/fig5_partial_effects
Rscript forest_plot_rna004.R
echo WRAPPER_DONE

