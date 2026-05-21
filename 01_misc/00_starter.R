################################################################################
################################################################################
# File name: 00_starter.R
# Author: M. Rotival
################################################################################
################################################################################

if(!grepl('Volumes', MPRA_DIR)) {
.libPaths(c("/pasteur/helix/projects/evo_immuno_pop/single_cell/resources/R_libs/4.1.0/", "/pasteur/appa/homes/mrotival/R/x86_64-pc-linux-gnu-library/4.0"))
}
shhh <- suppressPackageStartupMessages

cat('loading libraries')
shhh(library(dplyr))
try(shhh(library(ShortRead)))
shhh(library(data.table))
shhh(library(Biostrings))
shhh(library(ggplot2))
shhh(library(seqinr))
shhh(library(tictoc))
shhh(library(stringr))
shhh(library(ggrastr))
shhh(library(GenomicRanges))
shhh(library(liftOver))
shhh(library(tidyr))
shhh(library(DescTools))
shhh(library(scales))
cat('libraries loaded\n')


# BAMfile=${MPRA_DIR}/outs/alignedCRS_sorted.bam
#
# FASTQ_1=${MPRA_DIR}/HN00152006/MPRA-BC_1.fastq
# FASTQ_2=${MPRA_DIR}/HN00152006/MPRA-BC_2.fastq
# UMI=${MPRA_DIR}/HN00152006/MPRA-BC_UMI.fastq

IGSR <- "/pasteur/helix/projects/IGSR"
#MPRA_DIR <- "/pasteur/helix/projects/evo_immuno_pop/MPRA"
ANALYSIS_DIR <- "MPRA_count_exp6_analysis3"
FIGURE_DIR <- sprintf("%s/figures/%s", MPRA_DIR, ANALYSIS_DIR)

MPRA_DIR_GASPARD <- sprintf("%s/outs/RESULTS_gakerner", MPRA_DIR)
OUTS_DIR_GASPARD <- sprintf("%s/MPRAflow/outs", MPRA_DIR_GASPARD)

# OUT_THREELANES_ONEIDENTIFIER_raw=sprintf('%s/run2_3lanes_OneIdentifier/run2_3lanes_OneIdentifier_coords_to_barcodes.tsv',OUTS_DIR_GASPARD)
# OUT_THREELANES_ONEIDENTIFIER_filtered=sprintf('%s/run2_3lanes_OneIdentifier/run2_3lanes_OneIdentifier_filtered_coords_to_barcodes.tsv',OUTS_DIR_GASPARD)
# OUT_THREELANES_ONEIDENTIFIER_filtered_strict=sprintf('%s/run2_3lanes_OneIdentifier/run2_3lanes_OneIdentifier_filtered_coords_to_barcodes_strict.tsv',OUTS_DIR_GASPARD)
OUT_LIB_0_1_2_RAW <- sprintf("%s/MPRAflow/outs/LIB012_coords_to_barcodes_raw.tsv.gz", MPRA_DIR)
OUT_LIB_0_1_2_FILTERED <- sprintf("%s/MPRAflow/outs/LIB012_coords_to_barcodes_filtered.tsv.gz", MPRA_DIR)
OUT_LIB_0_1_2_FILTERED_STRICT <- sprintf("%s/MPRAflow/outs/LIB012_coords_to_barcodes_filtered_strict.tsv.gz", MPRA_DIR)

Indexes = fread(sprintf("%s/data/Index_lists/Index_list_exp1-6quater.txt", MPRA_DIR))
Indexes[,celline:=gsub('-ACE2','',celltype)]
Indexes[,COND_ID:=paste(celline,condition,sep='_')]
Barcode_DIR <- c("LIB0" = "01_barcode_lib_1", "LIB1" = "01_barcode_lib_2", "LIB2" = "01_barcode_lib_2", "LIB1+2" = "01_barcode_lib_2")
CountData_DIR <- c("exp1" = "02_MPRA_count_exp1_B7075-1", "exp2" = "02_MPRA_count_exp2_HN00196048", "exp3" = "02_MPRA_count_exp3_data",  "exp4" = "02_MPRA_count_exp4_data", "exp5" = "02_MPRA_count_exp5_data", "exp6" = "02_MPRA_count_exp6_data")

LIB_summary <- unique(Indexes[, .(library = Sample_ID, SID, celltype, condition, material, replicate, SETUP_ID, barcode_lib,celline,COND_ID,experiment)])
N_sample=LIB_summary[,.N]/2

LIB_summary <- LIB_summary[order(experiment,barcode_lib,celltype,factor(condition,c('NS','IFNA2b','DEX','SARS','IAV','TNFa')),replicate,material),]

# LIB_summary[, SETUP_ID := factor(SETUP_ID, names(color_setup))]
# LIB_summary <- LIB_summary[order(SETUP_ID, material)]

condition_summary=unique(LIB_summary[celltype!='Calu3' & celltype!='A549', .(celltype, celline, condition, COND_ID,analysis_name=paste0(celltype,'_',condition,'_all'))])[order(celline,factor(condition,c('NS','IFNA2b','DEX','SARS','IAV','TNFa')))]
condition_summary_reps=unique(LIB_summary[celltype!='Calu3' & celltype!='A549', .(celltype, celline, condition, COND_ID,analysis_name=SETUP_ID, experiment,barcode_lib)])[order(celline,factor(condition,c('NS','IFNA2b','DEX','SARS','IAV','TNFa')))]
technical_covariates=fread(sprintf('%s/data/00_technical_covariates_MPRA.txt',MPRA_DIR))

source(sprintf("%s/scripts/misc_plots.R", MPRA_DIR))
