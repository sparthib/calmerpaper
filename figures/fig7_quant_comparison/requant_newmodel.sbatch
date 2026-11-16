#!/bin/bash
#SBATCH --partition=shared
#SBATCH --cpus-per-task=4
#SBATCH --mem=90G
#SBATCH --time=8:00:00
#SBATCH --job-name=requant_newmodel
#SBATCH -o /dcs04/hicks/data/sparthib/lrs_coverage_project/lrs_cov/temp/requant_newmodel_out.txt
set -eo pipefail
module load conda_R/4.5.x
Rscript /dcs04/hicks/data/sparthib/lrs_coverage_project/lrs_cov/temp/requant_newmodel.R
echo WRAPPER_DONE
