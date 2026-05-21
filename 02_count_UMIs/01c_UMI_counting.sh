#!/bin/bash
################################################################################
################################################################################
# File name: 01c_UMI_counting.sh
# Author: M. Rotival
################################################################################
################################################################################

# call : sbatch 01c_UMI_counting.sh
################################################################################
# SLURM options
#SBATCH --job-name=MPRA_count
#SBATCH -o /pasteur/helix/projects/evo_immuno_pop/MPRA/logs/%x_logs/01c_getOR_%A_%a.log
#SBATCH -e /pasteur/helix/projects/evo_immuno_pop/MPRA/logs/%x_logs/01c_getOR_%A_%a.log
#SBATCH -c 1
#SBATCH --mem=50G
#SBATCH --qos=geh
#SBATCH --partition=gehbigmem,geh,depgg
#SBATCH --array=1-396
#SBATCH --mail-user=mrotival@pasteur.fr
#SBATCH --mail-type=END
################################################################################


COUNTER=${SLURM_ARRAY_TASK_ID}

MAX_SAMPLE=44
SUB_VALUES=('100%' '50%' '10%' '5%' '1%' '0.5%' '0.1%' '0.05%' '0.01%')

NTEST=$(($MAX_SAMPLE*$((${#SUB_VALUES[@]})))) #  396

# for COUNTER in `seq 1 $NTEST`; do
SUB=$(echo $(($COUNTER-1))/$MAX_SAMPLE | bc);
ID=$(echo "$COUNTER-$MAX_SAMPLE*$SUB" | bc); 
echo "$COUNTER" - ID= "$ID", SUB = "$SUB":"${SUB_VALUES[$SUB]}";
# done;

EVO_IMMUNO_POP_ZEUS="/pasteur/helix/projects/evo_immuno_pop"
SCRIPT_DIR="${EVO_IMMUNO_POP_ZEUS}/MPRA/scripts/MPRA_count_exp6_analysis3/"

module purge
module load R/4.1.0

if [[ ${SUB_VALUES[$SUB]} == "100%" ]]; then
  Rscript ${SCRIPT_DIR}/01c_UMI_counting.R --sample "$ID"
else
  Rscript ${SCRIPT_DIR}/01c_UMI_counting.R --sample "$ID" --sub_sample "${SUB_VALUES[$SUB]}"
fi

#
#
