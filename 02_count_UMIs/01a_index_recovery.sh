#!/bin/bash
################################################################################
################################################################################
# File name: 01a_index_recovery.sh
# Author: M. Rotival
################################################################################
################################################################################
# call: sbatch 01a_index_recovery.sh

################################################################################
# SLURM options
#SBATCH --job-name=MPRA_index
#SBATCH -o /pasteur/helix/projects/evo_immuno_pop/MPRA/logs/%x_logs/01a_getIndex_%A_%a.log
#SBATCH -e /pasteur/helix/projects/evo_immuno_pop/MPRA/logs/%x_logs/01a_getIndex_%A_%a.log
#SBATCH -c 1
#SBATCH --mem=230G
#SBATCH --qos=geh
#SBATCH --partition=gehbigmem,geh,depgg
#SBATCH --array=1-88
#SBATCH --mail-user=mrotival@pasteur.fr
#SBATCH --mail-type=END
################################################################################

ID=${SLURM_ARRAY_TASK_ID} 

MPRA_DIR="/pasteur/helix/projects/evo_immuno_pop/MPRA"
ANALYSIS_DIR="MPRA_count_exp6_analysis3"
SCRIPT_DIR="${MPRA_DIR}/scripts/${ANALYSIS_DIR}/"

module purge
module load R/4.1.0
Rscript ${SCRIPT_DIR}/01a_index_recovery.R --sample $ID

#
#
