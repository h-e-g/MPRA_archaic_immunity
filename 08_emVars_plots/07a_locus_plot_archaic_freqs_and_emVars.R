# check /Volumes/evo_immuno_pop/users/mrotival/Misc/plot_archaic_eQTLS for examples scripts for plot archaixc eQTLs and colocalization
# TODO :check for allele flipping between 1KG and SNP_info

# calls
# Rscript 06e_plot_archaic_freqs.R --chrom 6 --start 135.6e6 --end 140.1e6 --population "STU|PJL|GBR|IBS|CHB|JPT|AGTA|PAPUAN" --keyword IFNGR1
# Rscript 06e_plot_archaic_freqs.R --chrom 3 --start 45.8e6 --end 46.8e6 --population "STU|PJL|GBR|IBS|CHB|JPT|AGTA|PAPUAN" --keyword LZTFL1
# Rscript 06e_plot_archaic_freqs.R --chrom 4 --start 38.8e6 --end 38.9e6 --population "STU|PJL|GBR|IBS|CHB|JPT|AGTA|PAPUAN" --keyword TLR1
# Rscript 06e_plot_archaic_freqs.R --chrom 12 --start 113.1e6 --end 113.7e6 --population "STU|PJL|GBR|IBS|CHB|JPT|AGTA|PAPUAN" --keyword OAS1
# Rscript 06e_plot_archaic_freqs.R --chrom 18 --start 59.0e6 --end 62.0e6 --population "STU|PJL|GBR|IBS|CHB|JPT|AGTA|PAPUAN" --keyword BCL2
MPRA_DIR <- "/pasteur/helix/projects/evo_immuno_pop/MPRA"
ANALYSIS_DIR <- "MPRA_count_exp6_analysis3"
SCRIPT_DIR <- sprintf("%s/scripts/%s/", MPRA_DIR, ANALYSIS_DIR)


source(sprintf("%s/scripts/%s/00_starter.R", MPRA_DIR, ANALYSIS_DIR))

minBC <- 10

RUN_ID <- "RUN3_Z2_nBC10"
CRITERIA_ACTIVE_OUT <- "FDR5_scrambled5pct_FC0.2_GCnorm"
CRITERIA_EMVARS_OUT <- "EmVar_FDR5_FC.2"
# CRITERIA_DIFF_OUT <- "Diff_FDR5_FC.2"
# CRITERIA_EMVARS_DIFF <- "EmVarDiff_lfdr20_FC.2"
# STIM_ONLY <- FALSE


myCHROM <- 3
mySTART <- 45.7e6
myEND <- 46.7e6
KEYWORD <- "LZTFL1"
# myPOPs <- c('AGTA','PJL')
myPOPs <- c("STU", "PJL", "GBR", "IBS", "CHB", "JPT", "AGTA", "PAPUAN")

cmd <- commandArgs(trailingOnly = TRUE)
print(cmd)

for (i in seq_along(cmd)) {
  if (cmd[i] == "--run_dir" || cmd[i] == "-r") {
    RUN_ID <- cmd[i + 1]
  }
  if (cmd[i] == "--criteria_active" || cmd[i] == "-a") {
    CRITERIA_ACTIVE_OUT <- cmd[i + 1]
    CRITERIA_ACTIVE <- gsub("noActivityFilter_", "", CRITERIA_ACTIVE_OUT)
  }
  # if (cmd[i] == "--criteria_diff" || cmd[i] == "-d") {
  #   CRITERIA_DIFF_OUT <- cmd[i + 1]
  #   CRITERIA_DIFF <- gsub('noDiffFilter_','',CRITERIA_DIFF_OUT)
  # }
  if (cmd[i] == "--criteria_emvar" || cmd[i] == "-e") {
    CRITERIA_EMVARS_OUT <- cmd[i + 1]
    CRITERIA_EMVARS <- gsub("noEmVarFilter_", "", CRITERIA_EMVARS_OUT)
  }
  # if (cmd[i] == "--criteria_emvar_diff" || cmd[i] == "-ed") {
  #   CRITERIA_EMVARS_DIFF <- cmd[i + 1]
  # }
  # if(cmd[i] == "--stim_only" || cmd[i] == "-s") {
  #   STIM_ONLY <- as.logical(cmd[i + 1])
  # }
  if (cmd[i] == "--chrom" || cmd[i] == "-c") {
    myCHROM <- as.numeric(cmd[i + 1])
  } # ID of the CHROM to plot
  if (cmd[i] == "--start") {
    mySTART <- as.numeric(cmd[i + 1])
  } # start of the window to plot (bp)
  if (cmd[i] == "--end") {
    myEND <- as.numeric(cmd[i + 1])
  } # end of the window to plot (bp)
  if (cmd[i] == "--population" || cmd[i] == "--pop" || cmd[i] == "-p") {
    myPOPs <- str_split(cmd[i + 1], "\\|", simplify = TRUE)
  } # end of the window to plot (bp)
  if (cmd[i] == "--keyword" || cmd[i] == "-k") {
    KEYWORD <- cmd[i + 1]
  } # end of the window to plot (bp)
}
# STIM_ONLY_OUT <- ifelse(STIM_ONLY, "stimOnly", "all")



tic("loading oligo & SNP annotations")
# load annotation of oligos (beforz: CRS)
oligo_source <- fread(sprintf("%s/data/%s/00_oligo_annot_v2.txt", MPRA_DIR, ANALYSIS_DIR))
# load annotations of SNPs
SNP_annot_v4 <- fread(sprintf("%s/data/%s/00_SNP_annot_v4.txt", MPRA_DIR, ANALYSIS_DIR))
# load annotation per POP
selected_annot <- fread(sprintf("%s/data/%s/selected_snp_annotation.tsv.gz", MPRA_DIR, ANALYSIS_DIR))
# load clean annotation of introgression
selected_annot_wide <- fread(sprintf("%s/data/%s/selected_snp_annotation_wide.tsv.gz", MPRA_DIR, ANALYSIS_DIR))
SNP_annot_v5 <- merge(SNP_annot_v4[, -"Introgression_scenario"], selected_annot_wide, by = c("ID", "posID"))
toc()


##### active parameters
source(sprintf("%s/scripts/%s/03_00_parameter_activity.R", MPRA_DIR, ANALYSIS_DIR))
##### diff parameters
# source(sprintf("%s/scripts/%s/04_00_parameter_diff_activity.R", MPRA_DIR, ANALYSIS_DIR))
##### emVars parameters
source(sprintf("%s/scripts/%s/05_00_parameter_emVars.R", MPRA_DIR, ANALYSIS_DIR))
##### emVar diff parameters
# source(sprintf("%s/scripts/%s/06_00_parameter_emVars_Diff.R", MPRA_DIR, ANALYSIS_DIR))

IN_DIR <- sprintf("%s/data/%s/03a_aggMPRA_analyse/%s/", MPRA_DIR, ANALYSIS_DIR, RUN_ID)
ACTIVE_DIR <- sprintf("%s/data/%s/03b_computeFDR/%s/Activity/%s/", MPRA_DIR, ANALYSIS_DIR, RUN_ID, CRITERIA_ACTIVE)
# DIFF_DIR <- sprintf("%s/data/%s/03b_computeFDR/%s/Diff/%s/%s/", MPRA_DIR, ANALYSIS_DIR, RUN_ID, CRITERIA_DIFF, CRITERIA_ACTIVE_OUT)
EMVAR_DIR <- sprintf("%s/data/%s/03b_computeFDR/%s/emVars/%s/%s", MPRA_DIR, ANALYSIS_DIR, RUN_ID, CRITERIA_EMVARS, CRITERIA_ACTIVE_OUT)
# EMVARDIFF_DIR <- sprintf("%s/data/%s/03b_computeFDR/%s/EmVar_Diff/%s/%s/%s/%s/%s", MPRA_DIR, ANALYSIS_DIR, RUN_ID, CRITERIA_EMVARS_DIFF, CRITERIA_EMVARS_OUT, CRITERIA_DIFF_OUT, CRITERIA_ACTIVE_OUT, STIM_ONLY_OUT)


library(janitor)
library(data.table)
library(rtracklayer)
library(GenomicRanges)
library(stringr)
library(ggplot2)
library(VariantAnnotation)
library(dplyr)
library(ggrepel)
# library(dplyr)
# library(data.table)

FIGURE_DIR <- sprintf("%s/figures/%s/%s/06_locusView/%s/%s/", MPRA_DIR, ANALYSIS_DIR, RUN_ID, CRITERIA_EMVARS_OUT, CRITERIA_ACTIVE_OUT)
# FIGURE_DIR_LOCUS <- sprintf("%s/06_locusView/locus_aSNP_freqs", FIGURE_DIR)
dir.create(FIGURE_DIR, showWarnings = FALSE, recursive = TRUE)

# instruction_file <- sprintf("%s/scripts/%s/02a_MPRA_analyse_instructions_full_sorted.txt", MPRA_DIR, ANALYSIS_DIR)
# instruction_DT <- fread(instruction_file)
# VERSION_IN <- 3
# sub_sample <- 2000
# ANALYSIS_TYPE <- "emVars"
# results_signif_1cond <- fread(sprintf("%s/data/%s/MPRA_analyses/%s/allConditions_alphas_and_pvalues_v%s_max%sBC_annot.txt.gz", MPRA_DIR, ANALYSIS_DIR, ANALYSIS_TYPE, VERSION_IN, sub_sample))
# results_signif_int <- fread(sprintf("%s/data/%s/MPRA_analyses/%s/allInteractions_alphas_and_pvalues_v%s_max%sBC_annot.txt.gz", MPRA_DIR, ANALYSIS_DIR, ANALYSIS_TYPE, VERSION_IN, sub_sample))
# results_signif_1cond[, log2FC_der_vs_anc := ifelse(allele.1 == ANCESTRAL, logFC_2vs1, -logFC_2vs1)]
# results_signif_1cond[, log2FC_archaic_vs_modern := ifelse(!grepl("ancestral reintrogressed", SNP_type), log2FC_der_vs_anc, -log2FC_der_vs_anc)]
# results_signif_1cond_all <- fread(sprintf("%s/data/%s/MPRA_analyses/%s/", MPRA_DIR, ANALYSIS_DIR, ANALYSIS_TYPE, VERSION_IN, sub_sample))
# results_signif_1cond_all <- merge(results_signif_1cond, condition_summary, by = "analysis_name")


# FDR_DIR <- sprintf("%s/data/%s/03b_computeFDR/%s/", MPRA_DIR, ANALYSIS_DIR, RUN_ID)


###########################################################################################
############# @@ TODO load the correct version of the data #################################
###########################################################################################


# read oligo activity
oligo_activity_file <- sprintf("%s/oligo_activity__all__%s.tsv.gz", ACTIVE_DIR, CRITERIA_ACTIVE)
oligo_activity_obs <- fread(oligo_activity_file)

# read oligo targets
oligo_target <- fread(sprintf("%s/data/%s/00_oligo_targets_v4.txt", MPRA_DIR, ANALYSIS_DIR), sep = "\t")
oligo_target_collapsed <- fread(sprintf("%s/data/%s/00_oligo_targets_collapsed_v4.txt", MPRA_DIR, ANALYSIS_DIR), sep = "\t")
crs_targets <- unique(oligo_target_collapsed[, -"oligo"])

# read emVars
emVARs_annot_obs <- fread(file = sprintf("%s/all_emVars_annotated_celltype__%s.tsv", EMVAR_DIR, CRITERIA_ACTIVE_OUT))

# myPosID <- "10:122299806"
setnames(emVARs_annot_obs, "ANALYSIS_NAME", "analysis_name")

# response_emVars_obs <- fread(sprintf("%s/all_emVars_diff_obs_response.tsv", EMVARDIFF_DIR))

# load annotation of SNPs
SNP_annot <- SNP_annot_v5
# SNP_info = fread(sprintf("%s/single_cell/project/pop_eQTL/Z_OWEY/00_private_notshared/SNP_info_basics.tsv.gz", EIP))
# SNP_info_covid = fread(sprintf("%s/single_cell/project/pop_eQTL/Z_OWEY/00_private_notshared/SNP_info_covid.tsv.gz", EIP))
# SNP_info = merge(SNP_info, SNP_info_covid, by = c("posID", "rsID"))

# load annotation of archaic
Archaic_annot_clean <- fread(sprintf("%s/data/vcfs/all_common_introgressed/ALL_common_introgressed_allChr_v3.txt.gz", MPRA_DIR), sep = "\t")
# setnames(Archaic_annot_clean, "posID", "posID_hg19")
# myPOPs=c('STU','PJL','GBR','IBS','CHB','JPT','AGTA','PAPUAN')

Archaic_annot_locus <- Archaic_annot_clean[CHROM == myCHROM & POS >= mySTART & POS <= myEND, ]
Archaic_annot_locus[,Introgressed_from := allele_match]
#Archaic_annot_locus[,Introgressed_from := Introgression_source_final]
names(color_archaic_MPRA)[names(color_archaic_MPRA)=='Archaic (unknown)']='no match'
#names(color_archaic_MPRA)[4]='Undeterminded'

#Archaic_annot_locus[Introgressed_from == "Archaic (shared)", Introgressed_from := "Vindija/Denisova"]
Archaic_annot_locus[, label := "."]

rsID_file <- sprintf("%s/data/vcfs/common_introgressed_or_MPRA_tested_499086snps.vcf.gz", MPRA_DIR)
# rsID_file <- sprintf("%s/data/ALL_common_introgressed_allChr.vcf.gz", MPRA_DIR)
rsID_DT <- fread(rsID_file, skip = 116)
setnames(rsID_DT, "#CHROM", "CHROM")
rsID_DT <- rsID_DT[order(CHROM, POS, REF, ALT)]
rsID_DT[, INFO := NULL]
rsID_DT[, rsID := ID]
rsID_DT[, ID := paste(CHROM, POS, REF, ALT, "b37", sep = ":")]
Archaic_annot_locus <- merge(Archaic_annot_locus, rsID_DT[, .(ID, rsID)], by = "ID")

if (KEYWORD == "IFNGR1") {
  target_SNP <- Archaic_annot_locus[posID_hg19 == "6:137540719" & POP == "STU", ]
  target_SNP[, POP := NULL]
  target_SNP[, NAF := NULL]
  target_SNP_POP <- melt(target_SNP, id.vars = c("CHROM", "POS", "REF", "ALT", "posID_hg19"), measure.vars = c("YRI", "SGDP_African", "GBR", "IBS", "SGDP_WestEurasian", "CHB", "JPT", "SGDP_Eastasian", "SGDP_Agta", "SGDP_Papuan"), variable.name = "POP", value.name = "NAF")
  target_SNP_POP <- merge(target_SNP_POP, target_SNP, by = c("CHROM", "POS", "REF", "ALT", "posID_hg19"))
  Archaic_annot_locus <- rbind(Archaic_annot_locus, target_SNP_POP)
  Archaic_annot_locus[posID_hg19 == "6:137540719", label := rsID]
}

if (KEYWORD == "BCL2") {
  target_SNP <- Archaic_annot_locus[posID_hg19 == "18:60902328" & POP == "SGDP_Papuan", ]
  target_SNP[, POP := NULL]
  target_SNP[, NAF := NULL]
  target_SNP_POP <- melt(target_SNP, id.vars = c("CHROM", "POS", "REF", "ALT", "posID_hg19"), measure.vars = c("YRI", "SGDP_African", "GBR", "IBS", "SGDP_WestEurasian", "CHB", "JPT", "SGDP_Eastasian", "SGDP_Agta", "SGDP_Papuan"), variable.name = "POP", value.name = "NAF")
  target_SNP_POP <- merge(target_SNP_POP, target_SNP, by = c("CHROM", "POS", "REF", "ALT", "posID_hg19"))
  Archaic_annot_locus <- rbind(Archaic_annot_locus, target_SNP_POP)
  Archaic_annot_locus[posID_hg19 == "18:60902328", label := rsID]
}

if (KEYWORD == "LZTFL1") {
  target_SNP <- Archaic_annot_locus[posID_hg19 == "3:45859651" & POP == "PJL", ]
  target_SNP[, POP := NULL]
  target_SNP[, NAF := NULL]
  target_SNP_POP <- melt(target_SNP, id.vars = c("CHROM", "POS", "REF", "ALT", "posID_hg19"), measure.vars = c("YRI", "SGDP_African", "GBR", "IBS", "SGDP_WestEurasian", "CHB", "JPT", "SGDP_Eastasian", "SGDP_Agta", "SGDP_Papuan"), variable.name = "POP", value.name = "NAF")
  target_SNP_POP <- merge(target_SNP_POP, target_SNP, by = c("CHROM", "POS", "REF", "ALT", "posID_hg19"))
  Archaic_annot_locus <- rbind(Archaic_annot_locus, target_SNP_POP)
  Archaic_annot_locus[posID_hg19 == "3:45859651", label := rsID]
}

if (KEYWORD == "OAS1") {
  target_SNP <- Archaic_annot_locus[posID_hg19 == "12:113364366" & POP == "PJL", ]
  target_SNP[, POP := NULL]
  target_SNP[, NAF := NULL]
  target_SNP_POP <- melt(target_SNP, id.vars = c("CHROM", "POS", "REF", "ALT", "posID_hg19"), measure.vars = c("YRI", "SGDP_African", "GBR", "IBS", "SGDP_WestEurasian", "CHB", "JPT", "SGDP_Eastasian", "SGDP_Agta", "SGDP_Papuan"), variable.name = "POP", value.name = "NAF")
  target_SNP_POP <- merge(target_SNP_POP, target_SNP, by = c("CHROM", "POS", "REF", "ALT", "posID_hg19"))
  Archaic_annot_locus <- rbind(Archaic_annot_locus, target_SNP_POP)
  Archaic_annot_locus[posID_hg19 == "12:113364366", label := rsID]
  Archaic_annot_locus[,length(unique(crsID))] # 173
}

if (KEYWORD == "TLR1") {
  target_SNP <- Archaic_annot_locus[posID_hg19 == "4:38806019" & POP == "CHB", ]
  target_SNP[, POP := NULL]
  target_SNP[, NAF := NULL]
  target_SNP_POP <- melt(target_SNP, id.vars = c("CHROM", "POS", "REF", "ALT", "posID_hg19"), measure.vars = c("YRI", "SGDP_African", "GBR", "IBS", "SGDP_WestEurasian", "CHB", "JPT", "SGDP_Eastasian", "SGDP_Agta", "SGDP_Papuan"), variable.name = "POP", value.name = "NAF")
  target_SNP_POP <- merge(target_SNP_POP, target_SNP, by = c("CHROM", "POS", "REF", "ALT", "posID_hg19"))
  Archaic_annot_locus <- rbind(Archaic_annot_locus, target_SNP_POP)
  Archaic_annot_locus[posID_hg19 == "4:38806019", label := rsID]
}

if (KEYWORD == "STAT5B") {
  target_SNP <- Archaic_annot_locus[posID_hg19 == "17:40400760" & POP == "SGDP_Agta", ]
  target_SNP[, POP := NULL]
  target_SNP[, NAF := NULL]
  target_SNP_POP <- melt(target_SNP, id.vars = c("CHROM", "POS", "REF", "ALT", "posID_hg19"), measure.vars = c("YRI", "SGDP_African", "GBR", "IBS", "SGDP_WestEurasian", "CHB", "JPT", "SGDP_Eastasian", "SGDP_Agta", "SGDP_Papuan"), variable.name = "POP", value.name = "NAF")
  target_SNP_POP <- merge(target_SNP_POP, target_SNP, by = c("CHROM", "POS", "REF", "ALT", "posID_hg19"))
  Archaic_annot_locus <- rbind(Archaic_annot_locus, target_SNP_POP)
  Archaic_annot_locus[posID_hg19 == "17:40400760", label := rsID]
}



Archaic_annot_locus[, POP := toupper(gsub("SGDP_", "", POP))]

if (!is.null(myPOPs)) {
  Archaic_annot_locus <- Archaic_annot_locus[POP %chin% toupper(myPOPs), ]
}



# liftover to b38
library(liftOver)
chain <- import.chain(sprintf("%s//Chain_Files_For_Liftover/hg19ToHg38/hg19ToHg38.over.modified.chain", IGSR))
Archaic_annot_locus_GR37 <- makeGRangesFromDataFrame(Archaic_annot_locus, seqnames = "CHROM", start.field = "POS", end.field = "POS")
seqlevelsStyle(Archaic_annot_locus_GR37) <- "UCSC"
Archaic_annot_locus_GR38 <- liftOver(Archaic_annot_locus_GR37, chain)
Archaic_annot_locus[, POS_b38 := as.numeric(start(Archaic_annot_locus_GR38))]
mySTART_b38 <- mySTART - Archaic_annot_locus[, min(POS)] + Archaic_annot_locus[, min(POS_b38)]
myEND_b38 <- myEND - Archaic_annot_locus[, max(POS)] + Archaic_annot_locus[, max(POS_b38)]

#linecolor='lightgrey'
#linecolor='#00000000' # transparent
 linecolor=''


p_asnp <- ggplot()
if (KEYWORD == "LZTFL1") {
  p_asnp <- p_asnp + geom_vline(xintercept = Archaic_annot_locus[rsID == "rs17713054", POS_b38] / 1e6, lty = 2, col = "lightgrey")
  p_asnp <- p_asnp + scale_size_manual("target SNP", values = c("." = 1, "rs17713054" = 3)) + ylim(c(0,0.6))
}else{
	  p_asnp <- p_asnp + ylim(c(0,max(Archaic_annot_locus[NAF > 0,NAF])*1.05))
}
if(linecolor==''){
	#p_asnp <- p_asnp + geom_linerange(data = Archaic_annot_locus[NAF > 0, ], aes(x = POS_b38 / 1e6, ymax = NAF, ymin = 0, color = Introgressed_from), alpha = 0.5) 
	p_asnp <- p_asnp + geom_point(data = Archaic_annot_locus[NAF > 0, ], aes(x = POS_b38 / 1e6, y = NAF, color = Introgressed_from, fill = Introgressed_from, size = label), alpha = 0.8, shape = 21) 
	#p_asnp <- p_asnp + scale_color_manual("Population", values = color_populations_MPRA) 
	p_asnp <- p_asnp + scale_color_manual("Archaic", values = color_archaic_MPRA) 
}else{
	p_asnp <- p_asnp + geom_linerange(data = Archaic_annot_locus[NAF > 0, ], aes(x = POS_b38 / 1e6, ymax = NAF, ymin = 0), alpha = 0.5, color = linecolor) 
	p_asnp <- p_asnp + geom_point(data = Archaic_annot_locus[NAF > 0, ], aes(x = POS_b38 / 1e6, y = NAF, fill = Introgressed_from, size = label), alpha = 0.8, shape = 21, col=linecolor) 
}

p_asnp <- p_asnp + facet_grid(rows = vars(POP)) 
p_asnp <- p_asnp + scale_fill_manual("Archaic", values = color_archaic_MPRA) 
p_asnp <- p_asnp + xlab(paste0("Chromosome ", myCHROM, " (Mb)")) 
p_asnp <- p_asnp + ylab("Archaic allele frequency")
p_asnp <- p_asnp + xlim(c(mySTART_b38, myEND_b38) / 1e6) + theme_plot()
# if (KEYWORD == "IFNGR1") {
#   p_asnp <- p_asnp + geom_vline(xintercept = 137540719 / 1e6, lty = 2, col = "lightgrey")
#   p_asnp <- p_asnp + scale_size_manual("target SNP", values = c("." = 1, "rs121913171" = 3))
# }
# p_asnp <- p_asnp + geom_label_repel(data = Archaic_annot_locus[label != ".", ], aes(x = POS / 1e6, y = NAF, label = label, color = POP), box.padding = 0.8, max.overlaps = Inf, size = 2, segment.size = 0.2, show.legend = FALSE, fill = "white")
p_asnp <- p_asnp + theme(legend.box = "vertical") 
pdf(sprintf("%s/locus_%s-%s-%s_%s.pdf", FIGURE_DIR, myCHROM, round(mySTART / 1e6, 1), round(myEND / 1e6, 1), KEYWORD), height = 0.2 * 6.7 * (length(myPOPs) + .5), width = 7.2 * .55)
print(p_asnp)
dev.off()

# pdf(sprintf("%s/locus_%s-%s-%s_%s.pdf", FIGURE_DIR, myCHROM, round(mySTART / 1e6, 1), round(myEND / 1e6, 1), KEYWORD), height = 0.2 * 6.7 *  (length(myPOPs) + .5), width = 7.2 * .55  )
# print(p_asnp+ylim(c(0,0.6)))
# dev.off()

# VERSION_IN <- 4
# sub_sample <- 2000
# ANALYSIS_TYPE <- "emVars"
# results_signif_1cond <- fread(sprintf("%s/data/%s/MPRA_analyses/%s/allConditions_alphas_and_pvalues_v%s_max%sBC_annot.txt.gz", MPRA_DIR, ANALYSIS_DIR, ANALYSIS_TYPE, VERSION_IN, sub_sample))
# results_signif_int <- fread(sprintf("%s/data/%s/MPRA_analyses/%s/allInteractions_alphas_and_pvalues_v%s_max%sBC_annot.txt.gz", MPRA_DIR, ANALYSIS_DIR, ANALYSIS_TYPE, VERSION_IN, sub_sample))
# results_signif_1cond=merge(results_signif_1cond,SNP_annot[,.(posID,rsID,Introgressed_from,SNP_type,ANCESTRAL,DERIVED,INTROGRESSED)],by="posID")
# results_signif_1cond[, log2FC_der_vs_anc := ifelse(allele.1 == ANCESTRAL, logFC_2vs1, -logFC_2vs1)]
# results_signif_1cond[, log2FC_archaic_vs_modern := ifelse(!grepl("ancestral reintrogressed", SNP_type), log2FC_der_vs_anc, -log2FC_der_vs_anc)]
# results_signif_1cond_all <- merge(results_signif_1cond, condition_summary, by = "analysis_name")

# load annotation of CRSs

plot_region <- function(figdata, chr, start, end, minBC = 30, cellines, add_nonsig = FALSE, build = "b37") {
  color_direction <- c("up" = reddish, "ns" = "grey", "down" = blueish)

  region_DT <- figdata[CHROM == chr & POS_b38 > start & POS_b38 < end & (nBCs_g1a1 + nBCs_g1a2 + nBCs_g2a1 + nBCs_g1a2) > minBC & celline %in% cellines, ]

  p <- ggplot()
  p <- p + geom_point(data = region_DT[is_emVar_CRITERIA == FALSE, ], aes(
    x = POS_b38 / 1e6, y = log2FC_archaic_vs_modern
  ), alpha = 0.3, col = "lightgrey", size = 0.5)
  p <- p + scale_x_continuous(limit = c(start / 1e6, end / 1e6))
  if (add_nonsig) {
    signif_dots <- region_DT[posID %in% region_DT[is_emVar_CRITERIA == TRUE, posID]]
  } else {
    signif_dots <- region_DT[is_emVar_CRITERIA == TRUE, ]
  }
  p <- p + geom_pointrange(data = signif_dots, aes(
    x = POS_b38 / 1e6, y = log2FC_archaic_vs_modern, alpha = ifelse(is_emVar_CRITERIA, "sig", "ns"),
    ymin = log2FC_archaic_vs_modern - 1.96 * abs(log2FC.se),
    ymax = log2FC_archaic_vs_modern + 1.96 * abs(log2FC.se),
    color = COND_ID
  ), size = 0.2)
  p <- p + scale_color_manual(values = c(color_setup_simplified_norep))
  p <- p + scale_alpha_manual(values = c("sig" = 0.8, "ns" = 0.2))
  p <- p + theme_plot() + geom_hline(yintercept = 0, linetype = 2)
  p
}

cellines <- "A549|K562|HepG2"
cellines <- unlist(str_split(cellines, "\\|"))

fig_data <- emVARs_annot_obs
fig_data <- merge(fig_data, Archaic_annot_locus[, .(ID, POS_b38)], by = "ID")
region_DT <- fig_data[CHROM == myCHROM & POS_b38 > mySTART_b38 & POS_b38 < myEND_b38 & (nBCs_g1a1 + nBCs_g1a2 + nBCs_g2a1 + nBCs_g1a2) > minBC & celline %in% cellines, ]
p_emVars <- plot_region(fig_data, myCHROM, mySTART_b38, myEND_b38, minBC = 100, cellines = cellines, add_nonsig = FALSE, build = "b38") + theme(legend.box = "vertical") + guides(col = guide_legend(ncol = 4, byrow = TRUE))
pdf(sprintf("%s/locus_%s-%s-%s_%s_emVars_b38.pdf", FIGURE_DIR, myCHROM, round(mySTART_b38 / 1e6, 1), round(myEND_b38 / 1e6, 1), KEYWORD), height = 0.45 * 6.7, width = 7.2 * .55)
print(p_emVars)
dev.off()

p_emVars <- p_emVars + facet_grid(rows = vars(celline))
pdf(sprintf("%s/locus_%s-%s-%s_%s_emVars_split_b38.pdf", FIGURE_DIR, myCHROM, round(mySTART_b38 / 1e6, 1), round(myEND_b38 / 1e6, 1), KEYWORD), height = 0.65 * 6.7, width = 7.2 * .55)
print(p_emVars + ylab("log2FC archaic VS modern"))
dev.off()


OPENTARGETS_DIR <- sprintf("%s/open_targets/opentargets/platform/25.03/", IGSR)
list_COLOC_parquet_files <- dir(sprintf("%s/coloc", OPENTARGETS_DIR), pattern = "parquet")
list_CS_parquet_files <- dir(sprintf("%s/credibleSets", OPENTARGETS_DIR), pattern = "parquet")

library(data.table)
library(ggplot2)
library(arrow)
library(tictoc)
# coloc_results <- list()
# for (FILE in list_COLOC_parquet_files) {
#   coloc_results[[FILE]] <- as.data.table(read_parquet(sprintf("%s/coloc/%s", OPENTARGETS_DIR, FILE)))
# }
# coloc_results <- rbindlist(coloc_results)



# there seems to be a problem with CS of LZTFL1 (we thus used a manually downloaded version instead of theparquet files)
# CS_results <- list()
# for (FILE in list_CS_parquet_files) {
#   tic(FILE)
#   CS_results[[FILE]] <- as.data.table(read_parquet(sprintf("%s/credibleSets/%s", OPENTARGETS_DIR, FILE)))
#   toc()
# }
# CS_results <- rbindlist(CS_results)

if (KEYWORD == "LZTFL1" | KEYWORD == "OAS1") {
  ymax=ifelse(KEYWORD == "LZTFL1", 100, 20)
  #CS_covid <- rbindlist(CS_results[studyId == "GCST90134600" & chromosome == 3 & position > mySTART_b38 & position < myEND_b38, locus], idcol = "CS")
  CS_covid <- rbind(fread(sprintf("%s/data/covid_GWAS/credibleSets/LZTFL1_COVID_d9c38588fe4a5da237d148e00d37953b-credibleset.tsv", MPRA_DIR)),
                  fread(sprintf("%s/data/covid_GWAS/credibleSets/OAS1-3_COVID_6284eb71622bc629c4d47a0f49d55337-credibleset.tsv", MPRA_DIR)))
  
  covid_GWAS <- fread(sprintf("%s/data/covid_GWAS/GCST90134600.h.tsv.gz", MPRA_DIR))
  covid_GWAS[, variant_id:=gsub(":","_",SNP) ]
  covid_GWAS <- covid_GWAS[chromosome == myCHROM & base_pair_location > mySTART_b38 & base_pair_location < myEND_b38, ]
  covid_GWAS[, CS := ifelse(variant_id%in%CS_covid[,variant],1,0)]
  # covid_GWAS[, CS := 0]
  # covid_GWAS[match(CS_covid$variantId, gsub(":", "_", covid_GWAS$SNP)), CS := CS_covid[, CS]]
  # covid_GWAS_GR38=makeGRangesFromDataFrame(covid_GWAS,seqname='chromosome',start.field='base_pair_location',end.filed='base_pair_location')

  p_gwas <- ggplot(covid_GWAS[order(CS)], aes(x = base_pair_location / 1e6, y = -log10(p_value), col = as.character(CS))) +
    geom_point(size = 0.5, alpha = 0.8)
  p_gwas <- p_gwas + xlab(paste0("Chromosome ", myCHROM, " (Mb)")) + ylab("COVID-19 association (-log10(P))") 
  # p_gwas <- p_gwas + scale_color_manual(values = c("0" = "#000000", "1" = "#FF0000", "2" = "#0033AA", "3" = "#AAAA33", "4" = "#008800"))
  p_gwas <- p_gwas + scale_color_manual(values = c("0" = "#000000", "1" = "#FF0000"))
  p_gwas <- p_gwas + xlim(c(mySTART_b38, myEND_b38) / 1e6) + ylim(c(0, ymax)) + guides(color = "none") + theme_plot()

  pdf(sprintf("%s/locus_%s-%s-%s_%s_covid_gwas_b38.pdf", FIGURE_DIR, myCHROM, round(mySTART_b38 / 1e6, 1), round(myEND_b38 / 1e6, 1), KEYWORD), height = 0.2 * 6.7, width = 7.2 * .55)
  print(p_gwas)
  dev.off()


fig_data[,length(unique(crsID))]
# 1:     A549_IAV    40
# 2:    A549_SARS    48
# 3: HepG2_IFNA2b    14
# 4:     K562_DEX     9
# 5:  K562_IFNA2b     8

    SupTable6ab <- unique(fig_data[order(-abs(log2FC_archaic_vs_modern)),.(rsID,is_emVar_CRITERIA,COND_ID,log2FC_archaic_vs_modern,pval_emp,crsID,Introgression_source_top,max_intro_allele_freq, INTROGRESSED.allele, ID, cre_class_CRITERIA, POP_introgressed, is_COVID=ifelse(variantId_hg38%in%CS_covid[,variant],'yes',''))])
    SupTable6ab[,length(unique(crsID))]
    SupTable6ab[is_COVID=='yes',length(unique(crsID))]
    SupTable6ab[is_emVar_CRITERIA==TRUE,length(unique(crsID))]
    SupTable6ab[is_COVID=='yes' & is_emVar_CRITERIA==TRUE,length(unique(crsID))]
    SupTable6ab[is_COVID=='yes' & is_emVar_CRITERIA==TRUE,1]
    SupTable6ab[is_COVID=='yes' & pval_emp<0.05 & COND_ID=='A549_SARS' & crsID %in% SupTable6ab[is_emVar_CRITERIA==TRUE,crsID],1]
if(KEYWORD=='OAS1'){
  fwrite(SupTable6ab, file = sprintf("%s/SupTable6a_OAS1_emVars.tsv", FIGURE_DIR), sep = "\t")
  }
if(KEYWORD=='LZTFL1'){
  fwrite(SupTable6ab, file = sprintf("%s/SupTable6b_LZTFL1_emVars.tsv", FIGURE_DIR), sep = "\t")
  }


fig_data[pval_emp<.05,length(unique(crsID))]

fig_data[pval_emp<.05,length(unique(crsID)),by=COND_ID]


}
# 1:     A549_IAV    40
# 2:    A549_SARS    48
# 3: HepG2_IFNA2b    14
# 4:     K562_DEX     9
# 5:  K562_IFNA2b     8



################################################################################################################################
################################ make regional plot ############################################################################
################################################################################################################################

# plot_region <- function(chr, start, end) {
#   color_direction <- c("up" = reddish, "ns" = "grey", "down" = blueish)
#   region_DT <- emVARs_all_active[chromosome == chr & position > start & position < end & grepl("_all", ANALYSIS_NAME), ]
#   p <- ggplot(region_DT, aes(
#     x = position, y = abs(logFC_2vs1),
#     ymin = abs(logFC_2vs1) - 1.95 * abs(logFC.se),
#     ymax = abs(logFC_2vs1) + 1.95 * abs(logFC.se),
#     alpha = ifelse(FDR < 0.05, "sig", "ns"),
#     color = case_when(
#       FDR < 0.05 & logFC_2vs1 > 0 ~ "up",
#       FDR < 0.05 & logFC_2vs1 < 0 ~ "down",
#       FDR > 0.05 ~ "ns"
#     )
#   ))
#   p <- p + geom_pointrange(size = .5)
#   p <- p + facet_grid(rows = vars(ANALYSIS_NAME))
#   p <- p + scale_color_manual(values = color_direction)
#   p <- p + scale_alpha_manual(values = c("sig" = 1, "ns" = 0.2))
#   p <- p + theme_plot() + geom_hline(yintercept = 0, linetype = 2)
#   print(p)
# }


# pdf(sprintf("%s/9_region_plot_chr3_46Mb.pdf", FIGURE_DIR), height = N_sample)
# plot_region(3, 45e6, 47e6)
# dev.off()

# pdf(sprintf("%s/9_region_plot_chr3_46Mb.pdf", FIGURE_DIR), height = 4)
# plot_region(3, 45e6, 47e6)
# dev.off()



# pdf(sprintf("%s/9_region_plot_chr12_113Mb.pdf", FIGURE_DIR), height = 4)
# plot_region(12, 112e6, 114e6)
# dev.off()


# myCHROM <- 6
# mySTART <- 135.6e6
# myEND <- 140.1e6
# KEYWORD="IFNGR1"

# myCHROM <- 18
# mySTART <- 59.0e6
# myEND <- 62.0e6
# KEYWORD="BCL2"

# myCHROM <- 12
# mySTART <- 113.3e6
# myEND <- 113.6e6
# KEYWORD="OAS1"


# myCHROM <- 4
# mySTART <- 38.8e6
# myEND <- 38.9e6
# KEYWORD="TLR1"
