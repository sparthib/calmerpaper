#!/bin/bash
#SBATCH --partition=shared
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=0:30:00
#SBATCH --job-name=fig6_figure
#SBATCH -o /dcs04/hicks/data/sparthib/lrs_coverage_project/lrs_cov/figures/fig6_flagged_diagnostics/fig6_figure.%j.txt
set -eo pipefail
module load conda_R/4.5.x
cd /dcs04/hicks/data/sparthib/lrs_coverage_project/lrs_cov/figures/fig6_flagged_diagnostics
Rscript flagged_diagnostics_figure.R
echo WRAPPER_DONE

