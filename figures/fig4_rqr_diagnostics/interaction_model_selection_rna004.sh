#!/bin/bash
#SBATCH --job-name=fig4_rna004
#SBATCH --partition=shared
#SBATCH --mem=48G
#SBATCH --cpus-per-task=4
#SBATCH --time=1-00:00:00
#SBATCH --output=logs/fig4_rqr_selection.%j.txt

echo "**** Job starts ****"; date
# conda_R activate scripts reference unset env vars; relax -u while loading.
set +u
module load conda_R/4.5.x
set -u

SAMPLE="${1:-SRR30947490}"
echo "==== R: $(which Rscript) | sample: $SAMPLE ===="
Rscript interaction_model_selection_rna004.R "$SAMPLE"
RC=$?
echo "==== Rscript exit: $RC ===="
echo "**** Job ends ****"; date
exit $RC
