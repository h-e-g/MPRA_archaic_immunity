MPRA_DIR <- "/pasteur/helix/projects/evo_immuno_pop/MPRA"
ANALYSIS_DIR <- "MPRA_count_exp6_analysis3"

source(sprintf("%s/scripts/%s/00_starter.R", MPRA_DIR, ANALYSIS_DIR))
source(sprintf("%s/scripts/misc_plots.R", MPRA_DIR))
source(sprintf("%s/scripts/%s/00_load_associations.R", MPRA_DIR, ANALYSIS_DIR))

RUN_ID <- "RUN3_Z2_nBC10"
CRITERIA_ACTIVE_OUT <- "FDR5_scrambled5pct_FC0.2_GCnorm"
CRITERIA_EMVARS_OUT <- "EmVar_FDR5_FC.2"
EMVAR_LIST <- 'significant'

cmd <- commandArgs(trailingOnly = TRUE)
print(cmd)

for (i in seq_along(cmd)) {
  if (cmd[i] == "--run_dir" || cmd[i] == "-r") {
    RUN_ID <- cmd[i + 1]
  }
  if (cmd[i] == "--criteria_active" || cmd[i] == "-a") {
    CRITERIA_ACTIVE_OUT <- cmd[i + 1]
  }
  if (cmd[i] == "--criteria_emvar" || cmd[i] == "-e") {
    CRITERIA_EMVARS_OUT <- cmd[i + 1]
  }
  if (cmd[i] == "--emVar_list" || cmd[i] == "-l") {
    EMVAR_LIST <- cmd[i + 1] # either 'suggestive' or 'significant'
  }
}

if(!EMVAR_LIST %in% c('suggestive','significant')) {
  cat("WARNING: --emVar_list must be either 'suggestive' or 'significant', setting to 'significant' and continuing\n")
  EMVAR_LIST <- 'significant'
}

CRITERIA_ACTIVE <- gsub("noActivityFilter_", "", CRITERIA_ACTIVE_OUT)
CRITERIA_EMVARS <- gsub("noEmVarFilter_", "", CRITERIA_EMVARS_OUT)


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
# load annotations of SNPs
SNP_annot <- SNP_annot_v5
SNP_annot[, NON_INTROGRESSED := ifelse(INTROGRESSED.allele == ANCESTRAL, DERIVED, ANCESTRAL)]


##### active parameters
source(sprintf("%s/scripts/%s/03_00_parameter_activity.R", MPRA_DIR, ANALYSIS_DIR))
##### emVars parameters
source(sprintf("%s/scripts/%s/05_00_parameter_emVars.R", MPRA_DIR, ANALYSIS_DIR))

IN_DIR <- sprintf("%s/data/%s/03a_aggMPRA_analyse/%s/", MPRA_DIR, ANALYSIS_DIR, RUN_ID)
ACTIVE_DIR <- sprintf("%s/data/%s/03b_computeFDR/%s/Activity/%s/", MPRA_DIR, ANALYSIS_DIR, RUN_ID, CRITERIA_ACTIVE)
EMVAR_DIR <- sprintf("%s/data/%s/03b_computeFDR/%s/emVars/%s/%s", MPRA_DIR, ANALYSIS_DIR, RUN_ID, CRITERIA_EMVARS, CRITERIA_ACTIVE_OUT)

library(ggplot2)

FIGURE_DIR <- sprintf("%s/figures/%s/%s/04_emVars/%s/%s/", MPRA_DIR, ANALYSIS_DIR, RUN_ID, CRITERIA_EMVARS, CRITERIA_ACTIVE_OUT)

dir.create(sprintf("%s/Frequency_vs_emVars", FIGURE_DIR), showWarnings = FALSE, recursive = TRUE)

##### load oligo activity
oligo_activity_file <- sprintf("%s/oligo_activity__all__%s.tsv.gz", ACTIVE_DIR, CRITERIA_ACTIVE)
oligo_activity_class <- fread(oligo_activity_file)

oligo_activity_file <- sprintf("%s/all_oligos_annotated__%s.tsv.gz", ACTIVE_DIR, CRITERIA_ACTIVE)
oligo_activity_obs <- fread(oligo_activity_file)

emVARs_obs <- fread(file = sprintf("%s/all_emVars_annotated_celltype__%s.tsv", EMVAR_DIR, CRITERIA_ACTIVE_OUT))
emVARs_obs[, posID := crsID]
# emVARs_obs <- emVARs_obs[is_emVar_CRITERIA==TRUE,]

crs_targets <- fread(sprintf("%s/data/%s/00_oligo_targets_v4.txt", MPRA_DIR, ANALYSIS_DIR))

SNP_annot_emVars <- fread(file = sprintf("%s/SNP_annot_emVars__%s.tsv", EMVAR_DIR, CRITERIA_ACTIVE_OUT))
SNP_annot_emVars <- merge(SNP_annot_emVars, unique(crs_targets[method == "NearestTSS" & celltype == "any", .(posID = crsID, DistTSS_kb = -Score1, NearestTSS = TargetGene)]), by = c("posID"), allow.cartesian = TRUE, all.x = TRUE)


###### load frequency data
# from 00f
minSNP=0
minHap_length=0
introgressed_SNPs_annot <- list()
for (CHR in 1:22) {
  cat("CHR", CHR, "\n")
  #introgressed_SNPs_annot[[CHR]] <- try(fread(sprintf("%s/data/vcfs/annotate_introgressed/full_annotation/haplotypes_aSNP_with_introgressed_and_frequency_%s.tsv.gz", MPRA_DIR, CHR)))
  introgressed_SNPs_annot[[CHR]] <- try(fread(sprintf("%s/data/vcfs/annotate_introgressed/full_annotation/haplotypes_aSNP_with_introgressed_and_frequency_%s_minSNP_%s_minHapLength_%s.tsv.gz", MPRA_DIR, CHR, minSNP, minHap_length)))
}
# introgressed_SNPs_annot[[6]]=NULL
introgressed_SNPs_annot <- rbindlist(introgressed_SNPs_annot)

INTRO.freqs <- introgressed_SNPs_annot[, .(ID, POP, freq_introgressed_allele = target_SNP_freq.fixed, freq_introgressed_haplotype = target_SNP_on_archaic_hap_freq.fixed, GUESSED_INTROGRESSED.fixed, GUESSED_INTROGRESSED.nonfixed)]
worldwide_freqs <- INTRO.freqs[, .(allele_freqWW = mean(freq_introgressed_allele, na.rm = T), haplo_freqWW = mean(freq_introgressed_haplotype, na.rm = T)), by = ID]
cuts_haplo <- c(0, 0.03, 0.05, 0.07, 0.1, 0.15, 1)
cuts_allele <- c(0, 0.03, 0.05, 0.07, 0.1, 0.15, 0.2, 1)
worldwide_freqs[, haplotype_freq_bin := cut(haplo_freqWW, cuts_haplo, include.lowest = TRUE)]
worldwide_freqs[, allele_freq_bin := cut(allele_freqWW, cuts_allele, include.lowest = TRUE)]
worldwide_freqs[, haplotype_freq_10bin := cut(haplo_freqWW, quantile(haplo_freqWW, seq(0, 1, l = 10), na.rm = T), include.lowest = TRUE)]
worldwide_freqs[, allele_freq_10bin := cut(allele_freqWW, quantile(allele_freqWW, seq(0, 1, l = 10), na.rm = T), include.lowest = TRUE)]

###### add frequency data to emVars
emVARs_obs_freqs <- merge(emVARs_obs, worldwide_freqs, by = c("ID"), allow.cartesian = TRUE, all.x = TRUE)


########################################
###### Pct EmVar by frequency ##########
########################################


minSNP <- 5

Pct_emVar_byFreq <- get_Pct_emVars(emVARs_obs_freqs, split_by = c("haplotype_freq_10bin"), total_by = NULL, emVarlist = EMVAR_LIST)
Pct_emVar_byFreq <- make_it_long(Pct_emVar_byFreq, split_by = c("haplotype_freq_10bin"))
Pct_emVar_byFreq[N_test > minSNP, FDR := p.adjust(Pvalue, "fdr"), by = measure]
Pct_emVar_byFreq[N_test > minSNP][order(Pvalue), head(.SD, 5), by = measure]
fwrite(Pct_emVar_byFreq, sprintf("%s/Frequency_vs_emVars/Pct_emVar_%s_byFreq_haplo.tsv", FIGURE_DIR, EMVAR_LIST), sep = "\t", quote = FALSE)


Pct_emVar_byFreq <- get_Pct_emVars(emVARs_obs_freqs, split_by = c("allele_freq_10bin"), total_by = NULL, emVarlist = EMVAR_LIST)
Pct_emVar_byFreq <- make_it_long(Pct_emVar_byFreq, split_by = c("allele_freq_10bin"))
Pct_emVar_byFreq[N_test > minSNP, FDR := p.adjust(Pvalue, "fdr"), by = measure]
Pct_emVar_byFreq[N_test > minSNP][order(Pvalue), head(.SD, 5), by = measure]
fwrite(Pct_emVar_byFreq, sprintf("%s/Frequency_vs_emVars/Pct_emVar_%s_byFreq_allele.tsv", FIGURE_DIR,  EMVAR_LIST), sep = "\t", quote = FALSE)

for (i_MEASURE in seq_len(Measure_table[, .N])) {
  MEASURE <- Measure_table[i_MEASURE, measure]
  MEASURE_LABEL <- Measure_table[i_MEASURE, measure_label]
  MEASURE_NULL <- Measure_table[i_MEASURE, null]

  FigData <- Pct_emVar_byFreq[measure == MEASURE & N_test > minSNP, ]
  p <- ggplot(FigData)
  p <- p + geom_hline(yintercept = MEASURE_NULL, col = "lightgrey", linetype = "dashed")
  p <- p + geom_pointrange(aes(x = allele_freq_10bin, y = Pct * 100, ymin = Pct_lo * 100, ymax = Pct_hi * 100, alpha = FDR < 0.05), size = .3)
  # p <- p + facet_grid(rows=vars(measure)) + ylab("Percentage of alleles") + xlab("Population")
  if (is.na(MEASURE_NULL)) {
    MEASURE_NULL <- as.numeric(Pct_emVar_byFreq[measure == MEASURE, .(Pct = median(Pct))])
  }
  p <- p + geom_hline(yintercept = MEASURE_NULL * 100, linetype = 2)
  p <- p + theme_plot(rotate.x = 90, fontsize = 9) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + guides(col = "none")
  p <- p + scale_alpha_manual(values = alpha_TRUEFALSE)
  # p <- p + scale_color_manual(values=color_populations_MPRA)
  p <- p + ylab(MEASURE_LABEL) + xlab("Worldwide allele frequency bin (archaic allele)")
  pdf(sprintf("%s/Frequency_vs_emVars/01a_Pct_%s_%s_by_Freq.pdf", FIGURE_DIR, MEASURE, EMVAR_LIST), height = 2.3, width = 1.5)
  print(p)
  dev.off()
}



#####################################################
###### corrlation effect size vs frequency ##########
#####################################################

correlation_results <- list()
# cor allele vs absolute FC (shrinked)
correlation_results[["allele_shrinked"]] <- emVARs_obs_freqs[is_tested_cre == TRUE, .(
  r = cor(allele_freqWW, abs(log2FC_archaic_vs_modern) * (1 - lfdr), method = "p", use = "p"),
  P_pearson = cor.test(allele_freqWW, abs(log2FC_archaic_vs_modern) * (1 - lfdr), method = "p")$p.value,
  rho = cor(allele_freqWW, abs(log2FC_archaic_vs_modern) * (1 - lfdr), method = "s", use = "p"),
  P_spearman = cor.test(allele_freqWW, abs(log2FC_archaic_vs_modern) * (1 - lfdr), method = "s")$p.value,
  type = "shrinked",
  Freq = "allele_freq_WorldWide",
  stat = "correlation_absoluteFC_vs_freq"
), by = COND_ID]

# cor allele vs  absolute FC (not shrinked)
correlation_results[["allele_notshrinked"]] <- emVARs_obs_freqs[is_tested_cre == TRUE, .(
  r = cor(allele_freqWW, abs(log2FC_archaic_vs_modern), method = "p", use = "p"),
  P_pearson = cor.test(allele_freqWW, abs(log2FC_archaic_vs_modern), method = "p")$p.value,
  rho = cor(allele_freqWW, abs(log2FC_archaic_vs_modern), method = "s", use = "p"),
  P_spearman = cor.test(allele_freqWW, abs(log2FC_archaic_vs_modern), method = "s")$p.value,
  type = "not shrinked",
  Freq = "allele_freq_WorldWide",
  stat = "correlation_absoluteFC_vs_freq"
), by = COND_ID]

# cor haplo vs absolute FC (shrinked)
correlation_results[["haplo_shrinked"]] <- emVARs_obs_freqs[is_tested_cre == TRUE, .(
  r = cor(haplo_freqWW, abs(log2FC_archaic_vs_modern) * (1 - lfdr), method = "p", use = "p"),
  P_pearson = cor.test(haplo_freqWW, abs(log2FC_archaic_vs_modern) * (1 - lfdr), method = "p")$p.value,
  rho = cor(haplo_freqWW, abs(log2FC_archaic_vs_modern) * (1 - lfdr), method = "s", use = "p"),
  P_spearman = cor.test(haplo_freqWW, abs(log2FC_archaic_vs_modern) * (1 - lfdr), method = "s")$p.value,
  type = "shrinked",
  Freq = "haplo_freq_WorldWide",
  stat = "correlation_absoluteFC_vs_freq"
), by = COND_ID]

# cor haplo vs  absolute FC (not shrinked)
correlation_results[["haplo_notshrinked"]] <- emVARs_obs_freqs[is_tested_cre == TRUE, .(
  r = cor(haplo_freqWW, abs(log2FC_archaic_vs_modern), method = "p", use = "p"),
  P_pearson = cor.test(haplo_freqWW, abs(log2FC_archaic_vs_modern), method = "p")$p.value,
  rho = cor(haplo_freqWW, abs(log2FC_archaic_vs_modern), method = "s", use = "p"),
  P_spearman = cor.test(haplo_freqWW, abs(log2FC_archaic_vs_modern), method = "s")$p.value,
  type = "shrinked",
  Freq = "haplo_freq_WorldWide",
  stat = "correlation_absoluteFC_vs_freq"
), by = COND_ID]

# aggregate & print
correlation_results <- rbindlist(correlation_results)
correlation_results[, FDR_pearson := p.adjust(P_pearson, "fdr"), by = .(Freq, type, stat)]
correlation_results[, FDR_spearman := p.adjust(P_pearson, "fdr"), by = .(Freq, type, stat)]

fwrite(correlation_results, sprintf("%s/Frequency_vs_emVars/correlation_Freq_vs_AbsoluteFC.tsv", FIGURE_DIR), sep = "\t", quote = FALSE)


diff_freq_results <- list()
# diff allele freq emVar vs emVar
diff_freq_results[["allele"]] <- emVARs_obs_freqs[is_tested_cre == TRUE, .(
  median_freq_emVar = median(allele_freqWW[is_emVar_CRITERIA], na.rm = T),
  median_freq_non_emVar = median(allele_freqWW[!is_emVar_CRITERIA], na.rm = T),
  wilcoxP = wilcox.test(allele_freqWW[is_emVar_CRITERIA], allele_freqWW[!is_emVar_CRITERIA])$p.value,
  Freq = "allele_freq_WorldWide",
  stat = "diff_emVar_vs_non_emVar"
),
by = COND_ID
]


diff_freq_results[["haplo"]] <- emVARs_obs_freqs[is_tested_cre == TRUE, .(
  median_freq_emVar = median(haplo_freqWW[is_emVar_CRITERIA], na.rm = T),
  median_freq_non_emVar = median(haplo_freqWW[!is_emVar_CRITERIA], na.rm = T),
  wilcoxP = wilcox.test(haplo_freqWW[is_emVar_CRITERIA], haplo_freqWW[!is_emVar_CRITERIA])$p.value,
  Freq = "haplo_freq_WorldWide",
  stat = "diff_emVar_vs_non_emVar"
),
by = COND_ID
]

diff_freq_results <- rbindlist(diff_freq_results)

fwrite(diff_freq_results, sprintf("%s/Frequency_vs_emVars/diff_freq_emVar_vs_non_emVar.tsv", FIGURE_DIR), sep = "\t", quote = FALSE)







# # TODO : compute the frequency of expression associated emVars in each bin of frequency, for each tissue.
# ActiveFreq <- emVARs_annot_freqs[!ctrl & is_introgressed & !excluded,.(nSNP_tested=length(unique(crsID)),n_emVAR=length(unique(crsID[FDR<.05 & abs(log2FC_archaic_vs_modern)>.5]))),keyby=.(celltype,ANALYSIS_NAME,POP,cut(freq_introgressed_haplotype,c(0,0.05,0.2,1),include.lowest=TRUE))][!is.na(cut),]
# ActiveFreq[,Pct:=n_emVAR/nSNP_tested*100]
# ActiveFreq[,c('Pct_lo','Pct_hi'):=as.list(binom.test(n_emVAR,nSNP_tested)$conf.int*100),keyby=.(celltype,ANALYSIS_NAME,POP,cut)]


# ActiveFreq <- emVARs_annot_freqs[!ctrl & is_introgressed & !excluded,.(nSNP_tested=length(unique(crsID)),n_emVAR=length(unique(crsID[FDR<.05]))),keyby=.(celltype,ANALYSIS_NAME,cut(freq_introgressed_haplotype,c(0,0.05,0.2,1),include.lowest=TRUE))][!is.na(cut),]
# ActiveFreq[,Pct:=n_emVAR/nSNP_tested*100]
# ActiveFreq[,c('Pct_lo','Pct_hi'):=as.list(binom.test(n_emVAR,nSNP_tested)$conf.int*100),keyby=.(celltype,ANALYSIS_NAME,cut)]


# ActiveFreq <- emVARs_annot_freqs[!ctrl & is_introgressed & !excluded,.(nSNP_tested=length(unique(crsID)),n_emVAR=length(unique(crsID[FDR<.05 & abs(log2FC_archaic_vs_modern)>.5]))),keyby=.(POP,cut(freq_introgressed_haplotype,c(0,0.05,0.2,1),include.lowest=TRUE))][!is.na(cut),]
# ActiveFreq[,Pct:=n_emVAR/nSNP_tested*100]
# ActiveFreq[,c('Pct_lo','Pct_hi'):=as.list(binom.test(n_emVAR,nSNP_tested)$conf.int*100),keyby=.(POP,cut)]


# ActiveFreq <- emVARs_annot_freqs[!ctrl & is_introgressed & !excluded,.(nSNP_tested=length(unique(crsID)),n_emVAR=length(unique(crsID[FDR<.05 & abs(log2FC_archaic_vs_modern)>.5]))),keyby=.(celltype, cut(mean(freq_introgressed_haplotype),c(0,0.05,0.2,1),include.lowest=TRUE))][!is.na(cut),]
# ActiveFreq[,Pct:=n_emVAR/nSNP_tested*100]
# ActiveFreq[,c('Pct_lo','Pct_hi'):=as.list(binom.test(n_emVAR,nSNP_tested)$conf.int*100),keyby=.(celltype,cut)]


cat("All done\n")
q("no")
