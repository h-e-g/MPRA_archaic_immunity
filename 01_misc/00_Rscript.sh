#!/bin/bash
################################################################################
################################################################################
# File name: 00_Rscript.sh
# Author: M. Rotival
################################################################################
################################################################################

#SBATCH --mail-user=mrotival@pasteur.fr
#SBATCH --mail-type=END
#SBATCH -o /pasteur/helix/projects/evo_immuno_pop/MPRA/logs/0_logfile_Rscript_%A_%a.log


module purge # remove modules that may have been loaded by mistake
module load R/4.1.0
module load samtools # needed to call bcftools from within R (see ${SCRIPT_DIR}/../querySNPs.R-)
module load SLiM/4.0.1
module load gdal # needed for map plots (06a)
module load geos # needed for map plots (06a)
module load proj # needed for map plots (06a)
SCRIPT_DIR="/pasteur/helix/projects/evo_immuno_pop/MPRA/scripts/"
SCRIPT=$1
shift;
echo $*
Rscript ${SCRIPT_DIR}/${SCRIPT} $* ${SLURM_ARRAY_TASK_ID}
exit
