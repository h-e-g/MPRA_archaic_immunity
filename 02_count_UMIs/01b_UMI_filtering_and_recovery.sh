#!/bin/bash
################################################################################
################################################################################
# File name: 01b_UMI_filtering_and_recovery.sh
# Author: M. Rotival
################################################################################
################################################################################
# call : sbatch 01b_UMI_filtering_and_recovery.sh
################################################################################
# SLURM options
#SBATCH --job-name=MPRA_UMI
#SBATCH -o /pasteur/helix/projects/evo_immuno_pop/MPRA/logs/%x_logs/01b_getUMI_%A_%a.log
#SBATCH -e /pasteur/helix/projects/evo_immuno_pop/MPRA/logs/%x_logs/01b_getUMI_%A_%a.log
#SBATCH -c 1
#SBATCH --mem=230G
#SBATCH --qos=geh
#SBATCH --partition=gehbigmem,geh,depgg
#SBATCH --array=1-792
#SBATCH --mail-user=mrotival@pasteur.fr
#SBATCH --mail-type=END
################################################################################

COUNTER=${SLURM_ARRAY_TASK_ID}

MAX_SAMPLE=88
SUB_VALUES=('100%' '50%' '10%' '5%' '1%' '0.5%' '0.1%' '0.05%' '0.01%')

NTEST=$(($MAX_SAMPLE*$((${#SUB_VALUES[@]})))) # 576

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
  Rscript ${SCRIPT_DIR}/01b_UMI_filtering_and_recovery.R --sample "$ID"
else
  Rscript ${SCRIPT_DIR}/01b_UMI_filtering_and_recovery.R --sample "$ID" --sub_sample "${SUB_VALUES[$SUB]}"
fi
#
#

# MAX_SAMPLE=64
# SUB_VALUES=('100%' '50%' '10%' '5%' '1%' '0.5%' '0.1%' '0.05%' '0.01%')

# NTEST=$(($MAX_SAMPLE*$((${#SUB_VALUES[@]})))) 
# FAILED=('0.1%__4' '0.1%__2' '0.1%__3' '0.1%__1' '0.1%__9' '0.1%__10' '0.1%__7' '0.1%__8' '0.1%__13' '0.1%__14' '0.1%__5' '0.1%__6' '0.1%__11' '0.1%__12' '0.5%__2' '0.5%__1' '0.5%__10' '0.5%__7' '0.5%__8' '0.5%__25' '0.5%__26' '0.5%__6' '0.5%__32' '0.5%__36' '0.5%__12' '0.5%__30' '10%__2' '10%__1' '10%__10' '10%__8' '10%__12' '1%__10' '1%__12' '50%__2' '50%__1' '50%__10' '50%__8' '50%__6' '50%__12' '5%__2' '5%__1' '5%__12')


# for COUNTER in `seq 1 $NTEST`; do
#   SUB=`echo $(($COUNTER-1))/$MAX_SAMPLE | bc`;
#   ID=`echo "$COUNTER-$MAX_SAMPLE*$SUB" | bc`; 
#   #echo $COUNTER - ID= $ID, SUB = $SUB:${SUB_VALUES[$SUB]};
#   # echo "${SUB_VALUES[$SUB]}__$ID"
#   for i in ${FAILED[@]}; do 
#     # echo $i 
#     if [[ $i == "${SUB_VALUES[$SUB]}__$ID" ]]; then echo -n "$COUNTER,"; fi
#   done;
# done;
