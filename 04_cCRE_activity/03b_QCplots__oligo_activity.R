MPRA_DIR <- "/pasteur/helix/projects/evo_immuno_pop/MPRA"
ANALYSIS_DIR <- "MPRA_count_exp6_analysis3"

source(sprintf("%s/scripts/%s/00_starter.R", MPRA_DIR, ANALYSIS_DIR))
source(sprintf("%s/scripts/%s/00_load_associations.R", MPRA_DIR, ANALYSIS_DIR))



RUN_ID <- "RUN3_Z2_nBC10"
CRITERIA_ACTIVE <- "FDR5_scrambled5pct_FC0.2_GCnorm"

cmd <- commandArgs(trailingOnly = TRUE)
print(cmd)

for (i in seq_along(cmd)) {
  if (cmd[i] == "--run_dir" || cmd[i] == "-r") {
    RUN_ID <- cmd[i + 1]
  }
  if (cmd[i] == "--criteria" || cmd[i] == "-r") {
    CRITERIA_ACTIVE <- cmd[i + 1]
  }
}

source(sprintf("%s/scripts/%s/03_00_parameter_activity.R", MPRA_DIR, ANALYSIS_DIR))

IN_DIR <- sprintf("%s/data/%s/03a_aggMPRA_analyse/%s/", MPRA_DIR, ANALYSIS_DIR, RUN_ID)
ACTIVE_DIR <- sprintf("%s/data/%s/03b_computeFDR/%s/Activity/%s/", MPRA_DIR, ANALYSIS_DIR, RUN_ID, CRITERIA_ACTIVE)
FIGURE_DIR <- sprintf("%s/figures/%s/%s", MPRA_DIR, ANALYSIS_DIR, RUN_ID)
dir.create(FIGURE_DIR, showWarnings = FALSE, recursive = TRUE)

library(ggtext)

################################################################################
########    QC plots                                                     #######
################################################################################

cat("QC plots , oligo activity ")

ACTIVITY_DIR <- sprintf("%s/02_oligo_activity/%s", FIGURE_DIR, CRITERIA_ACTIVE)
dir.create(ACTIVITY_DIR, showWarnings = FALSE, recursive = TRUE)

# load oligo annotations
source(sprintf("%s/scripts/%s/02z__define_includedCRS.R", MPRA_DIR, ANALYSIS_DIR))

# oligo_source <- fread(sprintf("%s/data/%s/00_oligo_annot_v3_withTSS.txt", MPRA_DIR, ANALYSIS_DIR))
# #oligo_source[, allele.num := paste0("A", 1:.N, sep = ""), by = .(CRS, shift, strand)]
# oligo_type <- oligo_source[, .(type = paste(unique(type), collapse = "\\")), by = oligo]
# # list oligos associated with >1 source
# dup_oligo <- oligo_source[duplicated(oligo), unique(oligo)]



###############################################################################
###########            reading oligo activity                 #################
###############################################################################



#### TODO check these: numbers generation
# oligo_activity_obs <- fread(sprintf("%s/all_oligos_annotated.tsv.gz", ACTIVE_DIR))

# oligo_activity_obs[, SETUP_ID := factor(SETUP_ID, names(color_setup))]
# print(oligo_activity_obs[, .(N_active = sum(FDR < .01), Pct_active = mean(FDR < .01), N_up = sum(FDR < .05 & alpha.GC_norm > 1), Pct_up = mean(FDR < .01 & alpha.GC_norm > 1), N_down = sum(FDR < .05 & alpha.GC_norm < 1), Pct_down = mean(FDR < .05 & alpha.GC_norm < 1)), by = ANALYSIS_NAME])
# print(oligo_activity_obs[, .(N_strongly_active = sum(FDR_scrambled < .05), Pct_strongly_active = mean(FDR_scrambled < .05), N_strongly_up = sum(FDR_scrambled < .05 & alpha.GC_norm > 1), Pct_strongly_up = mean(FDR_scrambled < .05 & alpha.GC_norm > 1), N_strongly_down = sum(FDR_scrambled < .05 & alpha.GC_norm < 1), Pct_down = mean(FDR_scrambled < .05 & alpha.GC_norm < 1)), by = ANALYSIS_NAME])

# print(oligo_activity_obs[grep("all", ANALYSIS_NAME), .(N_active = sum(FDR < .01 & abs(log2(alpha.GC_norm)) > 0.5), Pct_active = mean(FDR < .01 & abs(log2(alpha.GC_norm)) > .5), N_up = sum(FDR < .01 & log2(alpha.GC_norm) > 0.5), Pct_up = mean(FDR < .01 & log2(alpha.GC_norm) > 0.5), N_down = sum(FDR < .01 & -log2(alpha.GC_norm) > 0.5), Pct_down = mean(FDR < .01 & -log2(alpha.GC_norm) > 0.5)), by = ANALYSIS_NAME])

# oligo_activity_obs[, .(active_CRS_any = any(FDR_scrambled < .05)), by = .(crsID, type_simple)][, .(sum(active_CRS_any), mean(active_CRS_any)), by = type_simple]
# oligo_activity_obs[, .(active_CRS_any = any(FDR < .01 & abs(log2(alpha.GC_norm)) > 0.5)), by = .(crsID, type_simple)][, .(sum(active_CRS_any), mean(active_CRS_any)), by = type_simple]
# oligo_activity_obs[, length(unique(crsID)), by = type]

# oligo_activity_obs[, .(active_CRS_lfdr = any(lfdr < .2 & abs(log2(alpha.GC_norm)) > 0.5)), by = .(crsID, type_simple)][, .(sum(active_CRS_any), mean(active_CRS_any)), by = type_simple]


oligo_activity_perm <- fread(sprintf("%s/all_oligos_annotated_perm__%s.tsv.gz", ACTIVE_DIR, CRITERIA_ACTIVE))
oligo_activity_perm_ctc <- oligo_activity_perm[power == 0 & boot == 0 & ANALYSIS_SUBTYPE == "celltype_cond", ]

oligo_activity_obs <- fread(sprintf("%s/all_oligos_annotated__%s.tsv.gz", ACTIVE_DIR, CRITERIA_ACTIVE))
oligo_activity_obs_ctc <- oligo_activity_obs[power == 0 & boot == 0 & ANALYSIS_SUBTYPE == "celltype_cond", ]
oligo_activity_obs_ctc_archaic <- oligo_activity_obs_ctc[oligo %in% tested_and_ctrl_oligos_final[type == "tested", oligo], ]

oligo_activity_obs_sample <- oligo_activity_obs[power == 0 & boot == 0 & ANALYSIS_SUBTYPE == "sample", ]


# setnames(condition_summary,'analysis_name','ANALYSIS_NAME')
# oligo_activity_perm=merge(oligo_activity_perm,by="ANALYSIS_NAME",allow.cartesian=TRUE)

p <- ggplot(oligo_activity_perm_ctc[condition == "NS", ], aes(x = pmin(2, alpha_CRITERIA), fill = celline, alpha = oligo_active_CRITERIA)) +
  geom_histogram()
p <- p + facet_grid(rows = vars(celline)) + scale_fill_manual(values = color_celline) + scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.5))
p <- p + theme_plot(lpos = "right") + xlab("normalized, GC-corrected activity (capped at 2)")

pdf(sprintf("%s/02a__permCRS_by_celltype__%s.pdf", ACTIVITY_DIR, CRITERIA_ACTIVE), height = 2.5, width = 3)
print(p)
dev.off()

# p <- ggplot(oligo_activity_perm[condition == "NS", ], aes(x = pmin(2, alpha.GC), fill = celline, alpha = FDR < 0.05)) +
#   geom_histogram()
# p <- p + facet_grid(rows = vars(celline)) + scale_fill_manual(values = color_celline) + scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.5))
# p <- p + theme_plot(lpos = "right") + xlab("GC-corrected activity (capped at 2)")

# pdf(sprintf("%s/02b__permCRS_alpha_GCadj_by_celltype__%s.pdf", ACTIVITY_DIR), height = 2.5, width = 3)
# print(p)
# dev.off()

# p <- ggplot(oligo_activity_perm[condition == "NS", ], aes(x = pmin(2, alpha.GC_stabVar), fill = celline, alpha = FDR < 0.05)) +
#   geom_histogram()
# p <- p + facet_grid(rows = vars(celline)) + scale_fill_manual(values = color_celline) + scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.5))
# p <- p + theme_plot(lpos = "right") + xlab("normalized, GC-corrected activity\n(Stabilized Variance, capped at 2)")

# pdf(sprintf("%s/02d__permCRS_alpha_GCadj_StabVar__by_celltype.pdf", ACTIVITY_DIR), height = 2.5, width = 3)
# print(p)
# dev.off()

#####################################################################################
######## distribution of alpha.GC_norm by oligo class ###############################
#####################################################################################

fig_data <- oligo_activity_obs_ctc_archaic[1:.N]
# fig_data[,alpha.GC:=alpha.GC_stabVar]

p <- ggplot(fig_data[condition == "NS"], aes(x = pmax(-2, pmin(2, log2(alpha_CRITERIA))), fill = factor(paste(celline, ifelse(oligo_strong_CRITERIA, "strong", "weak")), names(color_celline_2levels)), alpha = ifelse(oligo_active_CRITERIA, "sig", "ns")))
p <- p + geom_histogram(col = "white", binwidth = 0.1, linewidth = .1) + facet_grid(rows = vars(celline))
p <- p + theme_plot(fontsize = 12) + scale_fill_manual(values = rev(color_celline_2levels)) + scale_alpha_manual(values = c("ns" = 0.5, "sig" = 1))
p <- p + xlab("log2 (CRS activity)\n (RNA/DNA ratio)") + guides(fill = guide_legend(ncol = 2, byrow = FALSE), alpha = FALSE)
pdf(sprintf("%s/2b_alpha_by_celltype__NSonly__%s.pdf", ACTIVITY_DIR, CRITERIA_ACTIVE), height = 4, width = 3)
print(p)
dev.off()


p <- ggplot(fig_data, aes(x = pmax(-2, pmin(2, log2(alpha_CRITERIA))), fill = factor(paste(celline, ifelse(oligo_strong_CRITERIA, "strong", "weak")), names(color_celline_2levels)), alpha = ifelse(oligo_active_CRITERIA, "sig", "ns")))
p <- p + geom_histogram(col = "white", binwidth = 0.1, linewidth = .1) + facet_grid(rows = vars(celline))
p <- p + theme_plot(fontsize = 12) + scale_fill_manual(values = color_celline_2levels) + scale_alpha_manual(values = c("ns" = 0.5, "sig" = 1))
p <- p + xlab("log2 (CRS activity)\n (RNA/DNA ratio)") + guides(fill = guide_legend(ncol = 2, byrow = FALSE), alpha = FALSE)
pdf(sprintf("%s/2b_alpha_by_celltype__all_conditions__%s.pdf", ACTIVITY_DIR, CRITERIA_ACTIVE), height = 4, width = 3)
print(p)
dev.off()

p <- ggplot(oligo_activity_perm[condition == "NS", ], aes(x = pmax(-10, pmin(10, Zscore_CRITERIA)), fill = celline, alpha = oligo_active_CRITERIA)) +
  geom_histogram()
p <- p + facet_grid(rows = vars(celline)) + scale_fill_manual(values = color_celline) + scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.5))
p <- p + theme_plot(lpos = "right") + xlab("GC-corrected activity Z scores")

pdf(sprintf("%s/02c__permCRS_Zscores_by_celltype.pdf", ACTIVITY_DIR), height = 2.5, width = 3)
print(p)
dev.off()

FigData <- list(observed = oligo_activity_obs_ctc[condition == "NS", ], permuted = oligo_activity_perm_ctc[condition == "NS", ])
FigData <- rbindlist(FigData, idcol = "data_source")

p <- ggplot(FigData, aes(x = pmax(-10, pmin(10, Zscore_CRITERIA)), fill = celline, alpha = data_source)) +
  geom_histogram(position = "identity")
p <- p + facet_grid(rows = vars(celline)) + scale_fill_manual(values = color_celline) + scale_alpha_manual(values = c("observed" = .4, "permuted" = .7))
p <- p + theme_plot(lpos = "right") + xlab("GC-corrected activity Z scores")

pdf(sprintf("%s/02d__Zscores_by_celltype_obs_vs_perm__%s.pdf", ACTIVITY_DIR, CRITERIA_ACTIVE), height = 2.5, width = 3)
print(p)
dev.off()


FigData <- list(observed = oligo_activity_obs_ctc[condition == "NS", ], permuted = oligo_activity_perm_ctc[condition == "NS", ])
FigData <- rbindlist(FigData, idcol = "data_source")

p <- ggplot(FigData, aes(x = 2*pnorm(abs(Zscore_CRITERIA),lower.tail = FALSE), fill = celline)) + geom_histogram(position = "identity")
p <- p + facet_grid(rows = vars(celline),cols=vars(data_source)) + scale_fill_manual(values = color_celline) + scale_alpha_manual(values = c("observed" = .4, "permuted" = .7))
p <- p + theme_plot(lpos = "right") + xlab("activity p-values")

pdf(sprintf("%s/02e__pvalues_by_celltype_obs_vs_perm__%s.pdf", ACTIVITY_DIR, CRITERIA_ACTIVE), height = 2.5, width = 4)
print(p)
dev.off()


##############################################################################################################################
######################## Oligo activity : POSTIVE & NEGATIVE CONTROLs ########################################################
##############################################################################################################################


fig_data <- oligo_activity_obs_ctc[condition == "NS", ]
fig_data[, class := factor(type_simple, labels = c("Scrambled", "Promoter", "Archaic"), levels = c("scrambled", "Promoter", "other"))]

##########################################
##############  Figure 1D   ##############
##########################################

# p <- ggplot(fig_data, aes(x = class, y = log2(alpha.GC), fill = celline, alpha = class))
# p <- p + geom_violin(scale = "width") + geom_boxplot(fill = "white", alpha = .5, notch = TRUE, outlier.size = 0.7)
# p <- p + scale_alpha_manual(values = c("Archaic" = 0.5, "Promoter" = 1, "Scrambled" = 0.2)) + scale_fill_manual(values = color_celline)
# p <- p + theme_plot(rotate.x = 90, fontsize = 12) + facet_grid(cols = vars(celline))
# p <- p + ylab("log2 (CRS activity)\nRNA/DNA ratio") + xlab("CRS type") + guides(fill = FALSE, alpha = FALSE)

# pdf(sprintf("%s/03b__Positive_negative_Ctrl_alpha_GCadj_by_celltype.pdf", ACTIVITY_DIR), height = 3, width = 3)
# print(p)
# dev.off()

p <- ggplot(fig_data, aes(x = class, y = log2(alpha_CRITERIA), fill = celline, alpha = class))
p <- p + geom_violin(scale = "width") + geom_boxplot(fill = "white", alpha = .5, notch = TRUE, outlier.size = 0.7)
p <- p + scale_alpha_manual(values = c("Archaic" = 0.5, "Promoter" = 1, "Scrambled" = 0.2)) + scale_fill_manual(values = color_celline)
p <- p + theme_plot(rotate.x = 90, fontsize = 12) + facet_grid(cols = vars(celline))
p <- p + ylab("log2 (CRS activity)\nRNA/DNA ratio") + xlab("CRS type") + guides(fill = FALSE, alpha = FALSE)


pdf(sprintf("%s/03a__Positive_negative_Ctrl_alpha_by_celltype__%s.pdf", ACTIVITY_DIR, CRITERIA_ACTIVE), height = 3, width = 3)
print(p)
dev.off()


fig_data <- oligo_activity_obs[condition == "NS", ]
fig_data[, class := factor(type_simple, labels = c("Scrambled", "Promoter", "Archaic"), levels = c("scrambled", "Promoter", "other"))]

fig_data[, .(
  Prom_pvalue = wilcox.test(alpha_CRITERIA[class == "Promoter"], alpha_CRITERIA[class == "Scrambled"])$p.value,
  Archaic_pvalue = wilcox.test(alpha_CRITERIA[class == "Scrambled"], alpha_CRITERIA[class == "Archaic"])$p.value
), by = .(celline)]
#   celline   Prom_pvalue Archaic_pvalue
# 1:    A549  1.287812e-77   1.284009e-01
# 2:   HepG2 5.328201e-103   4.038343e-06
# 3:    K562 1.086746e-101   4.979533e-04
fig_data[, .(
  Prom = median(alpha_CRITERIA[class == "Promoter"]),
  Scrambled = median(alpha_CRITERIA[class == "Scrambled"]),
  Archaic = median(alpha_CRITERIA[class == "Archaic"])
),
by = .(celline)
]
#    celline     Prom Scrambled   Archaic
# 1:    A549 2.579454  1.020174 0.9961586
# 2:   HepG2 2.246831  1.063436 0.9955061
# 3:    K562 5.021059  1.034965 0.9963511

##########################################
##############  END Fig 1D  ##############
##########################################

# p <- ggplot(fig_data, aes(x = class, y = log2(alpha.GC_stabVar), fill = celline, alpha = class))
# p <- p + geom_violin(scale = "width") + geom_boxplot(fill = "white", alpha = .5, notch = TRUE, outlier.size = 0.7)
# p <- p + scale_alpha_manual(values = c("Archaic" = 0.5, "Promoter" = 1, "Scrambled" = 0.2)) + scale_fill_manual(values = color_celline)
# p <- p + theme_plot(rotate.x = 90) + facet_grid(cols = vars(celline))
# p <- p + ylab("log2 (CRS activity)\nRNA/DNA ratio") + xlab("CRS type") + guides(fill = FALSE, alpha = FALSE)


# pdf(sprintf("%s/03c__Positive_negative_Ctrl_alpha_GCadj_StabVar_by_celltype.pdf", ACTIVITY_DIR), height = 3, width = 3)
# print(p)
# dev.off()

# fig_data <- oligo_activity_obs[ANALYSIS_SUBTYPE=='celltype', ][(type == "archaic" | type_simple %in% c("Promoter", "scrambled")) & nBC_LIB0 > 30 & nBC_LIB2 > 30, ]
# p <- ggplot(fig_data, aes(x = class, y = log2(alpha), fill = celline, alpha = class))
# p <- p + geom_violin(scale = "width") + geom_boxplot(fill = "white", alpha = .5, notch = TRUE, outlier.size = 0.7)
# p <- p + scale_alpha_manual(values = c("Archaic" = 0.5, "Promoter" = 1, "Scrambled" = 0.2)) + scale_fill_manual(values = color_celline)
# p <- p + theme_plot(rotate.x = 90) + facet_grid(cols = vars(celline))
# p <- p + ylab("log2 (CRS activity)\nRNA/DNA ratio") + xlab("CRS type") + guides(fill = FALSE, alpha = FALSE)


# pdf(sprintf("%s/03d__Positive_negative_Ctrl_alpha_by_celltype.pdf", ACTIVITY_DIR), height = 3, width = 3)
# print(p)
# dev.off()

##############################################################################################################################
##############################################################################################################################
##############################################################################################################################

# oligo_activity_obs <- merge(oligo_activity_obs[!oligo %in% dup_oligo, ], oligo_source[!duplicated(oligo) & allele.num %in% c("A1", "A2"), .(oligo, allele.num, shift, strand)], by = c("oligo"))
# oligo_activity_obs_wide <- dcast(oligo_activity_obs, CRS + SETUP_ID + shift + strand ~ allele.num, value.var = list("OR_GCadj", "nBC_DNA", "nBC_RNA"))
# oligo_activity_obs_wide[is.na(nBC_DNA_A1),nBC_DNA_A1:=0]
# oligo_activity_obs_wide[is.na(nBC_DNA_A2),nBC_DNA_A2:=0]
# oligo_activity_obs_wide[is.na(nBC_RNA_A1),nBC_RNA_A1:=0]
# oligo_activity_obs_wide[is.na(nBC_RNA_A2),nBC_RNA_A2:=0]


# oligo_activity_obs <- fread(sprintf("%s/all_oligos_annotated.tsv.gz", ACTIVE_DIR))
# oligo_activity_obs <- oligo_activity_obs[power==0 & boot==0 & perm=='perm_0_0_0',]

nBC_list <- c(0, 1, 2, 3, 5, 10, 20, 30, 50, 100, 200, 300, 500)

# oligo_activity_obs <- merge(oligo_activity_obs, unique(LIB_summary[,.(SETUP_ID,barcode_lib,experiment)]), by = "SETUP_ID")
# # oligo_activity_obs <- merge(oligo_activity_obs, unique(Associations_Filtered[,.(oligo, total_count_BC, barcode_lib=barcode_library)]) , by = c("oligo", "barcode_lib"),allow.cartesian=TRUE)
# # oligo_activity_obs[,.(Cor=mean(cor(dcast(.SD,oligo+total_count_BC~replicate,value.var="OR_GCadj")[,-(1:2)],use='p'))),by=.(celltype,condition,barcode_lib,cut(total_count_BC,breaks=nBC_list))]

###################################################
##### correlation between replicates (Matrix) #####
###################################################

oligo_activity_obs <- merge(oligo_activity_obs, Associations_Filtered[barcode_library %in% c("LIB0", "LIB2"), .N, by = .(oligo, barcode_library)][, .(meanBC = mean(N)), by = oligo], by = c("oligo"), allow.cartesian = TRUE)

oligo_activity_sample <- merge(oligo_activity_obs_sample, technical_covariates[, .(SETUP_ID, replicate, plasmid_library, lentivirus_preparation)], by.x = "ANALYSIS_NAME", by.y = "SETUP_ID")
oligo_activity_sample[, nBC_sample := case_when(
  plasmid_library == "LIB2" ~ nBC_LIB2,
  plasmid_library == "LIB0" ~ nBC_LIB0,
  TRUE ~ -1
)]
oligo_activity_sample[, meanBC := mean(nBC_sample), by = oligo]

oligo_activity_by_cond <- oligo_activity_obs[ANALYSIS_SUBTYPE == "celltype_cond", ][!duplicated(paste(posID, COND_ID, type, shift, strand, allele.num)), ]

dir.create(sprintf("%s/Activity_correlation", ACTIVITY_DIR), showWarnings = FALSE, recursive = TRUE)
for (nbC_threshold in nBC_list) {
  cat("\ncomputing correlation with over ", nbC_threshold, "barcodes")
  oligo_activity_sample_wide <- dcast(oligo_activity_sample[nBC_sample > nbC_threshold], oligo ~ ANALYSIS_NAME, value.var = ALPHA_CRITERIA)
  COR_matrix <- cor(oligo_activity_sample_wide[, -c("oligo")], method = "s", use = "p")
  COR_matrix <- COR_matrix[condition_summary_reps$analysis_name, condition_summary_reps$analysis_name]
  rownames(COR_matrix) <- gsub("-ACE2", "", rownames(COR_matrix))
  colnames(COR_matrix) <- gsub("-ACE2", "", colnames(COR_matrix))

  library(corrplot)
  pdf(sprintf("%s/Activity_correlation/01_corr_btw_sample_activity_matrix_nBCover%s__%s.pdf", ACTIVITY_DIR, nbC_threshold, CRITERIA_ACTIVE))
  # corrplot(COR_matrix,tl.col=color_setup[condition_summary_reps$analysis_name],labels=names(color_setup_simplified)[match(condition_summary_reps$analysis_name,names(color_setup))])
  corrplot(COR_matrix, tl.col = color_setup[condition_summary_reps$analysis_name])
  dev.off()

  basal_oligo_activity <- apply(oligo_activity_sample_wide[, -c("oligo")], 1, mean, na.rm = T)
  COR_matrix <- cor(oligo_activity_sample_wide[, -c("oligo")] - basal_oligo_activity, method = "s", use = "p")
  COR_matrix <- COR_matrix[condition_summary_reps$analysis_name, condition_summary_reps$analysis_name]
  pdf(sprintf("%s/Activity_correlation/02_corr_btw_sample_activity_matrix_nBCover%s_deviation_from_average__%s.pdf", ACTIVITY_DIR, nbC_threshold, CRITERIA_ACTIVE))
  # corrplot(COR_matrix,tl.col=color_setup[condition_summary_reps$analysis_name],labels=names(color_setup_simplified)[match(condition_summary_reps$analysis_name,names(color_setup))])
  corrplot(COR_matrix, tl.col = color_setup[condition_summary_reps$analysis_name])
  dev.off()
}



##########################################
##############  Figure 1C   ##############
##########################################


######### spearman
nbC_threshold <- 10
oligo_activity_sample_wide <- dcast(oligo_activity_sample[nBC_sample > nbC_threshold], oligo ~ ANALYSIS_NAME, value.var = ALPHA_CRITERIA)
COR_matrix <- cor(oligo_activity_sample_wide[, -c("oligo")], method = "s", use = "p")
COR_matrix <- COR_matrix[condition_summary_reps$analysis_name, condition_summary_reps$analysis_name]
COR_DT <- as.data.table(melt(COR_matrix))
COR_DT <- merge(COR_DT, condition_summary_reps[, .(analysis_name, celline, condition)], by.x = "Var1", by.y = "analysis_name")
COR_DT <- merge(COR_DT, condition_summary_reps[, .(analysis_name, celline, condition)], by.x = "Var2", by.y = "analysis_name", suffix = c(".1", ".2"))
COR_REPS <- COR_DT[Var1 != Var2 & condition.1 == condition.2 & celline.1 == celline.2, ]
COR_REPS <- COR_REPS[!duplicated(paste(sort(c(Var1, Var2)), collapse = "_")), ]

COR_REPS[, ID := paste(sort(c(Var1, Var2)), collapse = "_"), by = .(Var1, Var2)]
COR_REPS <- COR_REPS[!duplicated(ID)]

# COR_REPS[,range(value)]
# [1] 0.6238572 0.9002973

pdf(sprintf("%s/Activity_correlation/03a_corr_btw_reps_spearman_10BC_per_cellline__%s.pdf", ACTIVITY_DIR, CRITERIA_ACTIVE), height = 3.3, 2.4)
p <- ggplot(COR_REPS, aes(x = factor(condition.1, c("NS", "IAV", "SARS", "IFNA2b", "DEX", "TNFa")), y = value, fill = paste(celline.1, condition.1, sep = "_")))
p <- p + geom_violin(scale = "width", col = NA) + geom_point(col = "black", alpha = .5, size = .5) + ylim(c(0, 1))
p <- p + xlab("") + ylab("spearman's rho") + scale_fill_manual(values = color_setup_simplified_norep)
p <- p + facet_grid(~celline.1, scale = "free_x", space = "free_x") + theme_plot(fontsize = 12, rotate.x = 90, lpos = "none")
print(p)
dev.off()


######### pearson
nbC_threshold <- 10
oligo_activity_sample_wide <- dcast(oligo_activity_sample[nBC_sample > nbC_threshold], oligo ~ ANALYSIS_NAME, value.var = ALPHA_CRITERIA)
COR_matrix <- cor(oligo_activity_sample_wide[, -c("oligo")], method = "p", use = "p")
COR_matrix <- COR_matrix[condition_summary_reps$analysis_name, condition_summary_reps$analysis_name]
COR_DT <- as.data.table(melt(COR_matrix))
COR_DT <- merge(COR_DT, condition_summary_reps[, .(analysis_name, celline, condition)], by.x = "Var1", by.y = "analysis_name")
COR_DT <- merge(COR_DT, condition_summary_reps[, .(analysis_name, celline, condition)], by.x = "Var2", by.y = "analysis_name", suffix = c(".1", ".2"))
COR_REPS <- COR_DT[Var1 != Var2 & condition.1 == condition.2 & celline.1 == celline.2, ]
COR_REPS <- COR_REPS[!duplicated(paste(sort(c(Var1, Var2)), collapse = "_")), ]

COR_REPS[, ID := paste(sort(c(Var1, Var2)), collapse = "_"), by = .(Var1, Var2)]
COR_REPS <- COR_REPS[!duplicated(ID)]

# COR_REPS[,range(value)]
# [1] 0.7872376 0.9833233

pdf(sprintf("%s/Activity_correlation/03a_corr_btw_reps_pearson_10BC_per_cellline__%s.pdf", ACTIVITY_DIR, CRITERIA_ACTIVE), height = 3.3, 2.4)
p <- ggplot(COR_REPS, aes(x = factor(condition.1, c("NS", "SARS", "IAV", "IFNA2b", "DEX", "TNFa")), y = value, fill = paste(celline.1, condition.1, sep = "_")))
p <- p + geom_violin(scale = "width", col = NA) + geom_point(col = "black", alpha = .5, size = .5) + ylim(c(0, 1))
p <- p + xlab("") + ylab("correlation betwen replicates\n(Pearson's r)") + scale_fill_manual(values = color_setup_simplified_norep)
p <- p + facet_grid(~celline.1, scale = "free_x", space = "free_x") + theme_plot(fontsize = 12, rotate.x = 90, lpos = "none")
p <- p + geom_hline(col = "grey", linetype = 2, yintercept = c(0, 1))
print(p)
dev.off()
##########################################
##############  END Fig1C   ##############
##########################################


#################################################################################################################################
################################################ correlation between replicates  ################################################
#################################################################################################################################

# oligo_activity_sample <- merge(oligo_activity_obs[ANALYSIS_SUBTYPE == "sample", ], technical_covariates[, .(SETUP_ID, plasmid_library, lentivirus_preparation, replicate)], by.x = "ANALYSIS_NAME", by.y = "SETUP_ID")
# oligo_activity_sample[, nBC_sample := case_when(
#   plasmid_library == "LIB2" ~ nBC_LIB2,
#   plasmid_library == "LIB0" ~ nBC_LIB0,
#   TRUE ~ -1
# )]
# oligo_activity_sample[,meanBC:=mean(nBC_sample),by=oligo]

pairwise_cor <- function(rep_matrix, method = "p") {
  COR <- cor(rep_matrix, method = method, use = "p")
  median(COR[upper.tri(COR)])
}

cor_reps <- list()
for (nbC_threshold in nBC_list) {
  cor_reps[[paste("nBC>", nbC_threshold)]] <- oligo_activity_sample[meanBC > nbC_threshold, .(.N, minBC = nbC_threshold, Cor = pairwise_cor(dcast(.SD, oligo ~ replicate, value.var = ALPHA_CRITERIA)[, -c("oligo", "meanBC")])), keyby = .(celltype, condition)]
}
cor_reps <- rbindlist(cor_reps)


########################################################################
##### Fig S3: correlation between replicates (by bin of BC number)  ####
########################################################################

cor_reps <- oligo_activity_sample[meanBC > 10, .(.N,
  meanBC = mean(meanBC),
  Cor = pairwise_cor(
    dcast(.SD, oligo + meanBC ~ replicate, value.var = ALPHA_CRITERIA)[, -c("oligo", "meanBC")],
    method = "s"
  )
), keyby = .(celline,
  condition,
  BC_bin = cut(meanBC, breaks = c(nBC_list, Inf), include.lowest = TRUE)
)]

########### spearman

pdf(sprintf("%s/01_correlation_btw_replicates_vs_barcodes_spearman__%s.pdf", ACTIVITY_DIR, CRITERIA_ACTIVE), height = 3.7, width = 4)
p <- ggplot(cor_reps, aes(x = as.numeric(BC_bin), y = Cor, col = factor(paste(celline, condition, sep = "_"), names(color_setup_simplified_norep)))) +
  ylab("Correlation between replicates (median)") +
  xlab("Number of barcodes")
p <- p + geom_point() + scale_color_manual(values = color_setup_simplified_norep) + theme_plot(rotate.x = 90, fontsize = 12, lpos = "right")
p <- p + geom_line() + scale_x_continuous(breaks = cor_reps[, unique(as.numeric(BC_bin))], labels = cor_reps[, unique(BC_bin)])
p <- p + guides(colour = guide_legend(byrow = FALSE, ncol = 1))
print(p)
dev.off()

########### pearson

cor_reps <- oligo_activity_sample[meanBC > 10, .(.N,
  meanBC = mean(meanBC),
  Cor = pairwise_cor(
    dcast(.SD, oligo + meanBC ~ replicate, value.var = ALPHA_CRITERIA)[, -c("oligo", "meanBC")],
    method = "p"
  )
), keyby = .(celline,
  condition,
  BC_bin = cut(meanBC, breaks = c(nBC_list, Inf), include.lowest = TRUE)
)]

###########

pdf(sprintf("%s/01_correlation_btw_replicates_vs_barcodes_pearson__%s.pdf", ACTIVITY_DIR, CRITERIA_ACTIVE), height = 3.7, width = 4)
p <- ggplot(cor_reps, aes(x = as.numeric(BC_bin), y = Cor, col = factor(paste(celline, condition, sep = "_"), names(color_setup_simplified_norep)))) +
  ylab("Correlation between replicates (median)") +
  xlab("Number of barcodes")
p <- p + geom_point() + scale_color_manual(values = color_setup_simplified_norep) + theme_plot(rotate.x = 90, fontsize = 12, lpos = "right")
p <- p + geom_line() + scale_x_continuous(breaks = cor_reps[, unique(as.numeric(BC_bin))], labels = cor_reps[, unique(BC_bin)])
p <- p + guides(colour = guide_legend(byrow = FALSE, ncol = 1))
print(p)
dev.off()
# dcast(oligo_activity_obs[celltype=='HepG2' & condition=='NS' & cut(meanBC,c(nBC_list,Inf))=='(5,10]',],oligo+meanBC~replicate,value.var="OR_GCadj")[,-(1:2)][,cor(R1,R2,use='p',method='s')]

############################################################################
##### END correlation between replicates (by bin of BC number), Fig S3 #####
############################################################################


#################################################################################
##### Fig SX: concordance of activity between replicates (by activityt bin)  ####
#################################################################################

get_sign_concordance=function(alpha){
 if(length(alpha)>1){
   mean(abs(apply(combn(sign(alpha),2),2,sum)/2),na.rm=T)
 }else{-999}
}


sign_contrast <- oligo_activity_obs_sample[,.(activity=mean(log2(alpha_CRITERIA)),concordance=get_sign_concordance(log2(alpha_CRITERIA)),mad=mad(log2(alpha_CRITERIA)),na.rm=T),keyby=.(COND_ID,oligo, GC)]
sign_contrast[concordance>0,]


fig_data=sign_contrast[concordance>0,.(activity=mean(activity),concordance=mean(concordance)),keyby=.(COND_ID,cut(activity,c(-Inf,-2,-1,-.5,seq(-0.3,0.3,by=0.05),.5,1,2,Inf)))]
p <- ggplot(fig_data,aes(x=activity,y=concordance,col=COND_ID)) + geom_line() + theme_plot(lpos='right')
p <- p + geom_vline(xintercept=c(-0.2,0.2),col='lightgrey',linetype=1)
p <- p + xlab('mean CRE activity (log2)')+ scale_color_manual(values=color_setup_simplified_norep) + coord_cartesian(xlim=c(-.6,.6),ylim=c(0.5,1)) + ylab('concordance in activity (positive/negative CRE) between replicates')
pdf(sprintf("%s/07d_oligo_activity_concordance_across_replicates_unbiased.pdf", ACTIVITY_DIR), height = 4, width = 4)
print(p)
dev.off()


active_oligos <- oligo_activity_obs_ctc_archaic[oligo_class_CRITERIA!='inactive',paste(oligo,COND_ID)]
fig_data_active <- sign_contrast[concordance>0 , .(activity=mean(activity),concordance=mean(concordance)),keyby=.(COND_ID,is_active=ifelse(paste(oligo,COND_ID)%in%active_oligos & abs(activity)>.2,'yes','no'),cut(activity,c(-Inf,-2,-1,-.5,seq(-0.3,0.3,by=0.05),.5,1,2,Inf)))]

p <- ggplot(fig_data_active,aes(x=activity,y=concordance,col=COND_ID)) + geom_line() + theme_plot(lpos='right') + facet_grid(cols=vars(is_active))
p <- p + geom_vline(xintercept=c(-0.2,0.2),col='lightgrey',linetype=1)
p <- p + xlab('mean CRE activity (log2)')+ scale_color_manual(values=color_setup_simplified_norep) + coord_cartesian(xlim=c(-2,2),ylim=c(0,1)) + ylab('concordance in activity (positive/negative CRE) between replicates')

pdf(sprintf("%s/07e_oligo_activity_concordance_across_alleles_unbiased_splitActiveCREs.pdf", ACTIVITY_DIR), height = 4, width = 6)
print(p)
dev.off()

fig_data=sign_contrast[concordance>0,.(GC=mean(GC),concordance=mean(concordance)),keyby=.(COND_ID,cut(GC,10))]
p <- ggplot(fig_data,aes(x=GC,y=concordance,col=COND_ID)) + geom_line() + theme_plot(lpos='right')
p <- p + xlab('mean CRE activity (log2)')+ scale_color_manual(values=color_setup_simplified_norep) + coord_cartesian(ylim=c(0.5,1)) + ylab('concordance in activity (positive/negative CRE) between replicates')
pdf(sprintf("%s/07f_oligo_activity_concordance_across_replicates_unbiased_byGC.pdf", ACTIVITY_DIR), height = 4, width = 4)
print(p)
dev.off()


fig_data=sign_contrast[concordance>0,.(activity=mean(activity),mad=mean(mad)),keyby=.(COND_ID,cut(activity,c(-Inf,-2,-1,-.5,seq(-0.3,0.3,by=0.05),.5,1,2,Inf)))]
p <- ggplot(fig_data,aes(x=activity,y=mad,col=COND_ID)) + geom_line() + theme_plot(lpos='right')
p <- p + geom_vline(xintercept=c(-0.2,0.2),col='lightgrey',linetype=1)
p <- p + xlab('mean CRE activity (log2)')+ scale_color_manual(values=color_setup_simplified_norep) + coord_cartesian(xlim=c(-.6,.6)) + ylab('Mean absolute deviation of activity between replicates')
pdf(sprintf("%s/07g_oligo_activity_mad_across_replicates_unbiased.pdf", ACTIVITY_DIR), height = 4, width = 4)
print(p)
dev.off()


fig_data=sign_contrast[concordance>0,.(GC=mean(GC),mad=mean(mad)),keyby=.(COND_ID,cut(GC,10))]
p <- ggplot(fig_data,aes(x=GC,y=mad,col=COND_ID)) + geom_line() + theme_plot(lpos='right')
#p <- p + geom_vline(xintercept=c(-0.2,0.2),col='lightgrey',linetype=1)
p <- p + xlab('mean GC content')+ scale_color_manual(values=color_setup_simplified_norep) + ylab('Mean absolute deviation of activity between replicates')
pdf(sprintf("%s/07h_oligo_activity_mad_across_replicates_unbiased_byGC.pdf", ACTIVITY_DIR), height = 4, width = 4)
print(p)
dev.off()

N_sample <- length(unique(oligo_activity_sample$ANALYSIS_NAME))

############################################################################################################
############################# correlation between DNA & RNA UMI (raw activity) #############################
############################################################################################################

# pdf(sprintf("%s/03b_UMI_per_oligo_DNA_vs_RNA__%s.pdf", ACTIVITY_DIR), height = 2.5, 2 + N_sample)
# p <- ggplot(oligo_activity_sample, aes(x = nUMI_DNA_g1, y = nUMI_RNA_g1, col = oligo_active_CRITERIA)) +
#   scale_y_continuous(trans = "log10") +
#   scale_x_continuous(trans = "log10")
# p <- p + rasterize(geom_point(alpha = .5), dpi = 200) + facet_grid(~COND_ID)
# p <- p + ylab("nUMI (RNA)") + xlab("nUMI (DNA), 1 dot = 1 oligo)")
# p <- p + theme_plot() + geom_abline(a = 0, b = 1) + geom_smooth(method = "lm", col = "black")
# print(p)
# dev.off()



# pdf(sprintf("%s/03b_UMI_per_oligo_DNA_vs_RNA__%s.pdf", ACTIVITY_DIR), height = 2.5, 2 + N_sample)
# p <- ggplot(oligo_activity_by_cond[condition == "NS", ], aes(x = nUMI_DNA_g1, y = nUMI_RNA_g1, col = oligo_active_CRITERIA))
# p <- p + scale_y_continuous(trans = "log10") + scale_x_continuous(trans = "log10")
# p <- p + rasterize(geom_point(alpha = .5), dpi = 200) + facet_grid(~COND_ID)
# p <- p + ylab("nUMI (RNA)") + xlab("nUMI (DNA), 1 dot = 1 oligo)")
# p <- p + theme_plot() + geom_abline(a = 0, b = 1) + geom_smooth(method = "lm", col = "black")
# print(p)
# dev.off()


# color_celline_2levels <- c(sapply(color_celline[c("HepG2", "A549", "K562")], mergeCols, "white"), sapply(color_celline[c("HepG2", "A549", "K562")], mergeCols, "black"))
# names(color_celline_2levels) <- paste(rep(c("HepG2", "A549", "K562"), 2), rep(c("FALSE", "TRUE"), e = 3))

# pdf(sprintf("%s/03b_UMI_per_oligo_DNA_vs_RNA_by_cellline.pdf", ACTIVITY_DIR), height = 4, 3.3)
# p <- ggplot(oligo_activity_by_cond[condition == "NS", ], aes(x = nUMI_DNA_g1, y = nUMI_RNA_g1, col = celline))
# p <- p + scale_y_continuous(trans = "log10") + scale_x_continuous(trans = "log10") + scale_alpha_manual(values = c("TRUE" = .5, "FALSE" = 0.3))
# p <- p + rasterize(geom_point(size = 0.5, aes(col = paste(celline, FDR < 0.01), alpha = FDR < .01)), dpi = 200) + facet_grid(rows = vars(celline)) + scale_color_manual(values = color_celline_2levels)
# p <- p + ylab("nUMI (RNA)") + xlab("nUMI (DNA), 1 dot = 1 oligo)")
# p <- p + theme_plot(lpos = "right") + geom_abline(a = 0, b = 1) + geom_smooth(method = "lm", col = "black")
# print(p)
# dev.off()
##########################################
##############  Figure 1B   ##############
##########################################

oligo_activity_by_cond[, nUMI_per_million_RNA := nUMI_RNA_g1 / sum(nUMI_RNA_g1) * 1e6, by = .(celline, condition)]
oligo_activity_by_cond[, nUMI_per_million_DNA := nUMI_DNA_g1 / sum(nUMI_DNA_g1) * 1e6, by = .(celline, condition)]

color_celline_2levels <- c(sapply(color_celline[c("HepG2", "A549", "K562")], mergeCols, "white"), sapply(color_celline[c("HepG2", "A549", "K562")], mergeCols, "black"))
names(color_celline_2levels) <- paste(rep(c("HepG2", "A549", "K562"), 2), rep(c("FALSE", "TRUE"), e = 3))

pdf(sprintf("%s/03b_UMI_per_oligo_DNA_vs_RNA_by_cellline_norm__%s.pdf", ACTIVITY_DIR, CRITERIA_ACTIVE), height = 3.3, 2)
p <- ggplot(oligo_activity_by_cond[condition == "NS", ], aes(x = nUMI_per_million_DNA, y = nUMI_per_million_RNA, col = celline))
p <- p + scale_y_continuous(trans = "log10") + scale_x_continuous(trans = "log10") + scale_alpha_manual(values = c("TRUE" = .5, "FALSE" = 0.3))
p <- p + rasterize(geom_point(size = 0.5, aes(col = paste(celline, oligo_active_CRITERIA), alpha = oligo_active_CRITERIA)), dpi = 200) + facet_grid(rows = vars(celline)) + scale_color_manual(values = color_celline_2levels)
p <- p + ylab("nUMI (RNA)") + xlab("nUMI (DNA),\n 1 dot = 1 oligo)")
p <- p + theme_plot(lpos = "none", fontsize = 12, rotate.x = 90) + geom_abline(a = 0, b = 1) #+ geom_smooth(method = "lm", col = "black")
print(p)
dev.off()


##########################################
##############  END Fig1B   ##############
##########################################

##########################################################################################################################
############################## END   correlation between DNA & RNA UMI (raw activity) ####################################
##########################################################################################################################

# oligo_activity_sample_wide <- dcast(oligo_activity_sample[allele.num!='',], crsID + COND_ID + shift + strand + type + replicate ~ allele.num, value.var = list("alpha.GC_norm", "nUMI_DNA_g1", "nUMI_RNA_g1"))
# oligo_activity_by_cond_wide <- dcast(oligo_activity_by_cond, crsID + COND_ID + cellline + condition+ shift + strand + type ~ allele.num, value.var = list("alpha.GC_norm"))



# pdf(sprintf("%s/03a_activity_per_oligo_allele_compared.pdf", ACTIVITY_DIR), height = 2.5, 2 + N_sample)
# p <- ggplot(oligo_activity_by_cond_wide, aes(x = A1, y = A2, col = FDR < 0.01)) +
#   scale_y_continuous(trans = "log10") +
#   scale_x_continuous(trans = "log10")
# p <- p + rasterize(geom_point(alpha = .5), dpi = 200) + facet_grid(~COND_ID)
# p <- p + ylab("activity allele 2") + xlab("activity allele 1, 1 dot = 1 oligo)")
# p <- p + theme_plot() + geom_abline(a = 0, b = 1) + geom_smooth(method = "lm", col = "black")
# print(p)
# dev.off()

##############################################################################################################################
############################################### GC content effect on activity) ###############################################
##############################################################################################################################

oligo_activity_obs_ctc[allele.label=='scrambled',cor.test(alpha,GC,method='s')[c('estimate','p.value')],by=celline]
oligo_activity_obs_ctc[allele.label=='scrambled' & condition=='NS',cor.test(alpha,GC,method='s')[c('estimate','p.value')],by=celline]



pdf(sprintf("%s/04a_GC_content_of_oligo_per_group.pdf", ACTIVITY_DIR), height = 3, width = 1 + N_sample)
p <- ggplot(oligo_activity_obs, aes(x = type, y = GC, fill = type))
p <- p + rasterize(geom_violin(scale = "width"), dpi = 200) + geom_boxplot(notch = TRUE, fill = "#FFFFFF88") + facet_grid(~COND_ID)
p <- p + ylab("GC content") + xlab("oligo type") + geom_hline(yintercept = 0, col = "grey")
p <- p + theme_plot(rotate.x = 90, lpos = "none")
print(p)
dev.off()

pdf(sprintf("%s/04b_GC_content_of_oligo_vs_alpha_NS__%s.pdf", ACTIVITY_DIR, CRITERIA_ACTIVE), height = 2.2, width = 7)
p <- ggplot(oligo_activity_by_cond[condition == "NS"], aes(x = GC, y = log2(alpha_CRITERIA), col = factor(type_simple, c("other", "Promoter", "scrambled"), labels = c("Promoter", "Scrambled", "Other"))))
p <- p + rasterize(geom_point(data = oligo_activity_by_cond[condition == "NS" & type_simple == "other", ], size = 0.3, alpha = .5, aes(col = type_simple)), dpi = 200)
p <- p + rasterize(geom_point(data = oligo_activity_by_cond[condition == "NS" & type_simple != "other", ], size = 0.3, alpha = .5, aes(col = type_simple)), dpi = 200)
p <- p + xlab("GC content") + ylab("oligo activity, log2(RNA/DNA)") + geom_hline(yintercept = 0, col = "grey")
p <- p + theme_plot(rotate.x = 90, fontsize = 11, lpos = "right") + facet_grid(~celline)
p <- p + geom_smooth(data = oligo_activity_by_cond[condition == "NS" & type_simple == "scrambled", ], method = "lm", col = "black", aes(fill = celline))
p <- p + scale_color_manual(values = c("other" = "grey", "Promoter" = "#1c81f4", "scrambled" = "#da4545"))
p <- p + guides(colour = guide_legend(override.aes = list(size = 2), ncol = 3), fill = guide_legend(ncol = 3), legend.box = "vertical") + scale_fill_manual(values = color_celline)
print(p)
dev.off()


pdf(sprintf("%s/04b_GC_content_of_oligo_vs_alpha_all__%s.pdf", ACTIVITY_DIR, CRITERIA_ACTIVE), height = 2.2, width = 7)
p <- ggplot(oligo_activity_by_cond, aes(x = GC, y = log2(alpha_CRITERIA), col = factor(type_simple, c("other", "Promoter", "scrambled"), labels = c("Promoter", "Scrambled", "Other"))))
p <- p + rasterize(geom_point(data = oligo_activity_by_cond[type_simple == "other", ], size = 0.3, alpha = .5, aes(col = type_simple)), dpi = 200)
p <- p + rasterize(geom_point(data = oligo_activity_by_cond[type_simple != "other", ], size = 0.3, alpha = .5, aes(col = type_simple)), dpi = 200)
p <- p + xlab("GC content") + ylab("oligo activity, log2(RNA/DNA)") + geom_hline(yintercept = 0, col = "grey")
p <- p + theme_plot(rotate.x = 90, fontsize = 11, lpos = "right") + facet_grid(~celline)
p <- p + geom_smooth(data = oligo_activity_by_cond[type_simple == "scrambled", ], method = "lm", col = "black", aes(fill = celline))
p <- p + scale_color_manual(values = c("other" = "grey", "Promoter" = "#1c81f4", "scrambled" = "#da4545"))
p <- p + guides(colour = guide_legend(override.aes = list(size = 2), ncol = 3), fill = guide_legend(ncol = 3), legend.box = "vertical") + scale_fill_manual(values = color_celline)
print(p)
dev.off()

# pdf(sprintf("%s/05a_UMI_peroligo_DNA_vs_RNA_FoldEnrich_vsBackGd.pdf", ACTIVITY_DIR), height = 3, width = 1 + N_sample)
# p <- ggplot(oligo_activity_obs[!oligo %in% dup_oligo, ], aes(x = type, y = log2(OR + 0.001), fill = type))
# p <- p + rasterize(geom_violin(scale = "width"), dpi = 200) + geom_boxplot(notch = TRUE, fill = "#FFFFFF88") + facet_grid(~SETUP_ID)
# p <- p + ylab("FE (RNA/DNA) /vs background") + xlab("oligo type") + geom_hline(yintercept = 0, col = "grey")
# p <- p + theme_plot(rotate.x = 90, lpos = "none")
# print(p)
# dev.off()

# pdf(sprintf("%s/05b_UMI_per_oligo_DNA_vs_RNA_FoldEnrich_vsBackGd_over100BC.pdf", ACTIVITY_DIR), height = 3, width = 1 + N_sample)
# p <- ggplot(oligo_activity_obs[!oligo %in% dup_oligo & nBC_DNA > 100 & nBC_RNA > 100, ], aes(x = type, y = log2(OR + 0.001), fill = type))
# p <- p + rasterize(geom_violin(scale = "width"), dpi = 200) + geom_boxplot(notch = TRUE, fill = "#FFFFFF88") + facet_grid(~SETUP_ID)
# p <- p + ylab("FE (RNA/DNA) /vs background") + xlab("oligo type") + geom_hline(yintercept = 0, col = "grey")
# p <- p + theme_plot(rotate.x = 90, lpos = "none")
# print(p)
# dev.off()

# pdf(sprintf("%s/05c_UMI_per_oligo_DNA_vs_RNA_FoldEnrich_vsBackGd_simplified.pdf", ACTIVITY_DIR), height = 3, width = 1 + N_sample)
# p <- ggplot(oligo_activity_obs[!oligo %in% dup_oligo, ], aes(x = type_simple, y = log2(OR + 0.001), fill = type_simple))
# p <- p + rasterize(geom_violin(scale = "width"), dpi = 200) + geom_boxplot(notch = TRUE, fill = "#FFFFFF88") + facet_grid(~SETUP_ID)
# p <- p + ylab("FE (RNA/DNA) /vs background") + xlab("oligo type") + geom_hline(yintercept = 0, col = "grey")
# p <- p + theme_plot(rotate.x = 90, lpos = "none") + ylim(c(-4, 4))
# print(p)
# dev.off()



# pdf(sprintf("%s/06a_OR_GCadj_vsBackGd.pdf", ACTIVITY_DIR), height = 3, width = 1 + N_sample)
# p <- ggplot(oligo_activity_obs[!oligo %in% dup_oligo, ], aes(x = type, y = log2(OR_GCadj + 0.001), fill = type))
# p <- p + rasterize(geom_violin(scale = "width"), dpi = 200) + geom_boxplot(notch = TRUE, fill = "#FFFFFF88") + facet_grid(~SETUP_ID)
# p <- p + ylab("FE (RNA/DNA) /vs background") + xlab("oligo type") + geom_hline(yintercept = 0, col = "grey")
# p <- p + theme_plot(rotate.x = 90, lpos = "none")
# print(p)
# dev.off()

# pdf(sprintf("%s/06b_OR_GCadj_vsBackGd_over100BC.pdf", ACTIVITY_DIR), height = 3, width = 1 + N_sample)
# p <- ggplot(oligo_activity_obs[!oligo %in% dup_oligo & nBC_DNA > 100 & nBC_RNA > 100, ], aes(x = type, y = log2(OR_GCadj + 0.001), fill = type))
# p <- p + rasterize(geom_violin(scale = "width"), dpi = 200) + geom_boxplot(notch = TRUE, fill = "#FFFFFF88") + facet_grid(~SETUP_ID)
# p <- p + ylab("FE (RNA/DNA) /vs background") + xlab("oligo type") + geom_hline(yintercept = 0, col = "grey")
# p <- p + theme_plot(rotate.x = 90, lpos = "none")
# print(p)
# dev.off()

# pdf(sprintf("%s/06c_OR_GCadj_vsBackGd_simplified.pdf", ACTIVITY_DIR), height = 3, width = 1 + N_sample)
# p <- ggplot(oligo_activity_obs[!oligo %in% dup_oligo, ], aes(x = type_simple, y = log2(OR_GCadj + 0.001), fill = type_simple))
# p <- p + rasterize(geom_violin(scale = "width"), dpi = 200) + geom_boxplot(notch = TRUE, fill = "#FFFFFF88") + facet_grid(~SETUP_ID)
# p <- p + ylab("FE (RNA/DNA) /vs background") + xlab("oligo type") + geom_hline(yintercept = 0, col = "grey")
# p <- p + theme_plot(rotate.x = 90, lpos = "none") + ylim(c(-4, 4))
# print(p)
# dev.off()

# pdf(sprintf("%s/06d_OR_GCadj_vsBackGd_simplified_over100BC.pdf", ACTIVITY_DIR), height = 3, width = 1 + N_sample)
# p <- ggplot(oligo_activity_obs[!oligo %in% dup_oligo & nBC_DNA > 100 & nBC_RNA > 100, ], aes(x = type_simple, y = log2(OR_GCadj + 0.001), fill = type_simple))
# p <- p + rasterize(geom_violin(scale = "width"), dpi = 200) + geom_boxplot(notch = TRUE, fill = "#FFFFFF88") + facet_grid(~SETUP_ID)
# p <- p + ylab("FE (RNA/DNA) /vs background") + xlab("oligo type") + geom_hline(yintercept = 0, col = "grey")
# p <- p + theme_plot(rotate.x = 90, lpos = "none") + ylim(c(-4, 4))
# print(p)
# dev.off()



exclude <- dcast(oligo_activity_obs_ctc_archaic[,.N,by=.(crsID, celline, condition, allele.label)], crsID + celline + condition ~ allele.label, value.var = "N")[INTROGRESSED != 1 | `NON-INTROGRESSED` != 1, .(crsID, celline, condition)]

oligo_activity_obs_ctc_archaic_compare_alleles <- oligo_activity_obs_ctc_archaic[!paste(crsID, celline, condition) %in% exclude[, paste(crsID, celline, condition)], ]
oligo_activity_obs_ctc_archaic_compare_alleles <- dcast(oligo_activity_obs_ctc_archaic_compare_alleles, crsID + celline + condition + COND_ID ~ allele.label, value.var = "alpha_CRITERIA")
setnames(oligo_activity_obs_ctc_archaic_compare_alleles, c("NON-INTROGRESSED"), "NON_INTROGRESSED")

oligo_activity_obs_ctc_archaic_compare_alleles[, COND_ID2 := factor(gsub("_", "-", COND_ID), levels = gsub("_", "-", names(color_setup_simplified_norep)))]
p <- ggplot(oligo_activity_obs_ctc_archaic_compare_alleles, aes(x = log2(NON_INTROGRESSED), y = log2(INTROGRESSED), col = COND_ID))
p <- p + geom_abline(slope = 1, intercept = 0, col = "grey", linetype = "dashed") + rasterize(geom_point(size = .3, alpha = .3), dpi = 300)
p <- p + scale_color_manual(values = color_setup_simplified_norep)

cor_df <- oligo_activity_obs_ctc_archaic_compare_alleles[, .(
  pearson = cor(log2(NON_INTROGRESSED), log2(INTROGRESSED), method = "p", use='p'),
  spearman = cor(log2(NON_INTROGRESSED), log2(INTROGRESSED), method = "s", use='p')
), by = COND_ID2]
cor_df[,labels:=paste0("r = ", round(pearson, 2), "\n", "rho = ", round(spearman, 2))]
p <- p + geom_text(
  data = cor_df,
  aes(x = -2, y = 4, label = labels,col=COND_ID2),
  hjust = 0, vjust = 1, size = 2
)
p <- p + facet_wrap(~COND_ID2, ncol = 4) + theme_plot(fontsize = 12, lpos = "none")

pdf(sprintf("%s/07a_oligo_activity_across_alleles.pdf", ACTIVITY_DIR), height = 4, width = 4)
print(p)
dev.off()


oligo_activity_obs_ctc_archaic_compare_alleles[,.(rho=cor(INTROGRESSED,NON_INTROGRESSED,use='p',method='s')),by=COND_ID][,range(rho)]

sign_contrast <- oligo_activity_obs_ctc_archaic_compare_alleles[,.(activity=mean((log2(INTROGRESSED)+log2(NON_INTROGRESSED))/2),concordance=mean(sign(log2(INTROGRESSED))==sign(log2(NON_INTROGRESSED)),na.rm=T)),keyby=.(COND_ID,cut((log2(INTROGRESSED)+log2(NON_INTROGRESSED)/2),c(-Inf,-2,-1,-.5,seq(-0.3,0.3,by=0.05),.5,1,2,Inf)))]
p <- ggplot(sign_contrast,aes(x=activity,y=concordance,col=COND_ID)) + geom_line() + theme_plot(lpos='right')
p <- p + xlab('mean CRE activity (log2)')+ scale_color_manual(values=color_setup_simplified_norep) + xlim(c(-0.3,0.3)) + ylab('concordance in activity (positive/negative CRE) between alleles')

pdf(sprintf("%s/07b_oligo_activity_concordance_across_alleles.pdf", ACTIVITY_DIR), height = 4, width = 4)
print(p)
dev.off()


sign_contrast <- oligo_activity_obs_ctc_archaic_compare_alleles[,.(activity=mean(log2(INTROGRESSED)),concordance=mean(sign(log2(INTROGRESSED))==sign(log2(NON_INTROGRESSED)),na.rm=T)),keyby=.(COND_ID,cut(log2(INTROGRESSED),c(-Inf,-2,-1,-.5,seq(-0.3,0.3,by=0.05),.5,1,2,Inf)))]
p <- ggplot(sign_contrast,aes(x=activity,y=concordance,col=COND_ID)) + geom_line() + theme_plot(lpos='right')
p <- p + geom_vline(xintercept=c(-0.2,0.2),col='lightgrey',linetype=1)
p <- p + xlab('mean CRE activity (log2)')+ scale_color_manual(values=color_setup_simplified_norep) + coord_cartesian(xlim=c(-.6,.6)) + ylab('concordance in activity (positive/negative CRE) between alleles')
pdf(sprintf("%s/07c_oligo_activity_concordance_across_alleles_unbiased.pdf", ACTIVITY_DIR), height = 4, width = 4)
print(p)
dev.off()	

active_cres <- oligo_activity_obs_ctc_archaic[cre_class_CRITERIA!='inactive',paste(crsID,COND_ID)]
sign_contrast_active <- oligo_activity_obs_ctc_archaic_compare_alleles[,.(activity=mean(log2(INTROGRESSED)),concordance=mean(sign(log2(INTROGRESSED))==sign(log2(NON_INTROGRESSED)),na.rm=T)),keyby=.(is_active=ifelse(paste(crsID,COND_ID)%in%active_cres,'yes','no'), COND_ID,cut(log2(INTROGRESSED),c(-Inf,-2,-1,-.5,seq(-0.3,0.3,by=0.05),.5,1,2,Inf)))]
p <- ggplot(sign_contrast_active,aes(x=activity,y=concordance,col=COND_ID)) + geom_line() + theme_plot(lpos='right') + facet_grid(cols=vars(is_active))
p <- p + geom_vline(xintercept=c(-0.2,0.2),col='lightgrey',linetype=1)
p <- p + xlab('mean CRE activity (log2)')+ scale_color_manual(values=color_setup_simplified_norep) + coord_cartesian(xlim=c(-.6,.6)) + ylab('concordance in activity (positive/negative CRE) between alleles')
pdf(sprintf("%s/07d_oligo_activity_concordance_across_alleles_unbiased_splitActiveCREs.pdf", ACTIVITY_DIR), height = 4, width = 6)
print(p)
dev.off()

# get_correlation <- function(nBC) {
#   oligo_activity_obs_wide[nBC_DNA_A1 > nBC & nBC_DNA_A2 > nBC & nBC_RNA_A1 > nBC & nBC_RNA_A2 > nBC, .(r = cor(OR_GCadj_A1, OR_GCadj_A2, use = "p"), rho = cor(OR_GCadj_A1, OR_GCadj_A2, method = "s", use = "p"), nBC), by = SETUP_ID]
# }
# correlation_CRS_btw_alleles <- rbindlist(lapply(nBC_list, get_correlation))
# correlation_CRS_btw_alleles <- merge(correlation_CRS_btw_alleles, condition_summary_reps, by.x = "SETUP_ID", by.y = "analysis_name")
# pdf(sprintf("%s/07a_correlation_oligo_activity_across_alleles_byBCthreshold.pdf", FIGURE_DIR), height = 5, width = 7)
# p <- ggplot(correlation_CRS_btw_alleles, aes(x = nBC, y = rho, col = SETUP_ID))
# p <- p + rasterize(geom_point(alpha = .5), dpi = 200) + geom_line(alpha = 0.5) + facet_grid(~ paste(experiment, barcode_lib, sep = " - "))
# p <- p + ylab("correlation of CRS activity between alleles\n(spearman's rho)") + xlab("minimum number of barcodes considered")
# p <- p + theme_plot(rotate.x = 90) + scale_color_manual(values = color_setup)
# print(p)
# dev.off()
# ## do something with ggpairs

# oligo_activity_obs[, mean_BC_oligo := mean(nBC_DNA + nBC_RNA) / 2, by = oligo]
# oligo_activity_obs[, logOR_GC_adj := log2(OR_GCadj)]
# oligo_activity_obs_wide_setups <- dcast(oligo_activity_obs[order(mean_BC_oligo)][mean_BC_oligo > 20, ], oligo + cut(mean_BC_oligo, c(20, 50, 100, 200, 500, Inf)) ~ SETUP_ID, value.var = list("logOR_GC_adj"))

# pdf(sprintf("%s/08a_oligo_activity_across_celltypes_byBCthreshold.pdf", FIGURE_DIR), width = 7, height = 7)
# p <- ggpairs(oligo_activity_obs_wide_setups[order(mean_BC_oligo)], columns = 2 + seq_len(N_sample), aes(color = mean_BC_oligo), lower = list(continuous = wrap("points", alpha = 0.3, size = 0.05)))
# p <- p + theme_plot(rotate.x = 90)
# p <- p + scale_color_viridis(discrete = TRUE, direction = -1) + scale_fill_viridis(discrete = TRUE, direction = -1)
# print(p)
# dev.off()




# ggpairs(oligo_activity_obs_wide_setups)
# p <- ggplot(UMI_perBC_annot_DNA_vs_RNA,aes(x=BC_class,y=DNA,fill=celltype))
# p <- p + theme_plot(rotate.x=90,lpos='none') + geom_violin(scale='width') + geom_boxplot(fill='#FFFFFF88',notch=T)
# p <- p + scale_y_continuous(trans='log10') + scale_fill_manual(values=color_celline) + facet_grid(~celltype)

# pdf(sprintf('%s/Number_of_UMI_perDNA_BC_vs_celltype.pdf',FIGURE_DIR),height=3,width=3)
# print(p)
# dev.off()

# UMI_perBC[material == "DNA", .(
#   r = cor(nUMI_perBC, N_assoc, use = "p", method = "p"),
#   rho = cor(nUMI_perBC, N_assoc, use = "p", method = "s")
# ), by = SETUP_ID]

# CT_compare <- dcast(UMI_perBC_annot_DNA_vs_RNA, barcode1 + known_BC1 + shannon + condition ~ SETUP_ID, value.var = c("DNA", "RNA"))

# pdf(sprintf('%s/UMI_perBC_vs_knownBC_RNA_to_DNA_ratio_v2.pdf',FIGURE_DIR),height=3,width=2)
# p <- ggplot(UMI_perBC_annot_DNA_vs_RNA[sample(seq_len(.N),1e6)],aes(x=known_BC1,y=log2((1+RNA)/(1+DNA)),fill=known_BC1))
# p <- p + rasterize(geom_violin(scale='width',bw=0.5),dpi=200) + rasterize(geom_boxplot(notch=TRUE,fill='#FFFFFF',alpha=.5),dpi=200)
# p <- p + ylab('log2 UMI ratio (RNA/DNA)')+xlab('is known (BC)')
# p <- p + theme_plot() + facet_grid(~celltype) +geom_hline(yintercept=log2(3),col='grey')
# print(p)
# dev.off()

# pdf(sprintf('%s/UMI_perBC_vs_knownBC_RNA_to_DNA_ratio_highC_v2.pdf',FIGURE_DIR),height=3,width=2)
# p <- ggplot(UMI_perBC_annot_DNA_vs_RNA[shannon>10][sample(seq_len(.N),1e6)],aes(x=known_BC1,y=log2((1+RNA)/(1+DNA)),fill=known_BC1))
# p <- p + rasterize(geom_violin(scale='width',bw=0.5),dpi=200) + rasterize(geom_boxplot(notch=TRUE,fill='#FFFFFF',alpha=.5),dpi=200)
# p <- p + ylab('log2 UMI ratio (RNA/DNA)')+xlab('is known (BC)')
# p <- p + theme_plot() + facet_grid(~celltype) +geom_hline(yintercept=log2(3),col='grey')
# print(p)
# dev.off()

## what percentage of BC-associated oligos are found in each library
# N_assoc_oligo <- Associations_Filtered[, length(unique(oligo))]
# N_assoc_BC <- Associations_Filtered[, length(unique(BC))]


# UMI_per_oligo[!is.na(oligo),.(N_oligo=length(unique(oligo)),recovery=100*length(unique(oligo))/N_assoc_oligo),keyby=library]
# >99.8%
# library N_oligo recovery
# 1: Lib3 (HepG2-DNA) 11441 99.81679
# 2: Lib4 (Calu3-DNA) 11443 99.83423
# 3: Lib2 (Calu3-RNA) 11456 99.94765
# 4: Lib1 (HepG2-RNA) 11452 99.91276

## what percentage of BC-associated oligo are found in each index
# UMI_per_oligo_index[!is.na(oligo),.(N_oligo=length(unique(oligo)),recovery=100*length(unique(oligo))/N_assoc_oligo),by=.(library,index)]
# >99.4%

# library      index N_oligo recovery
# 1: Lib3 (HepG2-DNA) AAGTTGACGA 11404 99.49398
# 2: Lib4 (Calu3-DNA) AGATCGAATA 11403 99.48526
# 3: Lib2 (Calu3-RNA) AGTAGGCGGA 11448 99.87786
# 4: Lib1 (HepG2-RNA) ATCCAAGTTG 11437 99.78189
# 5: Lib2 (Calu3-RNA) CGCGCGACTT 11449 99.88658
# 6: Lib4 (Calu3-DNA) CGTCCGGCAT 11394 99.40674
# 7: Lib2 (Calu3-RNA) GAATAGACGC 11448 99.87786
# 8: Lib4 (Calu3-DNA) TGCCTTCCTA 11437 99.78189
# 9: Lib1 (HepG2-RNA) TTGAACTAAG 11444 99.84296
# 10: Lib1 (HepG2-RNA) ACTACTTGCT 11431 99.72954
# 11: Lib3 (HepG2-DNA) GAAGCTGAAG 11398 99.44163
# 12: Lib3 (HepG2-DNA) TGGAATGACG 11414 99.58122


## what percentage of oligo-associated BC are found in each library
# UMI_perBC_lib[!is.na(oligo) & known_BC1,.(N_BC=length(unique(barcode1)),recovery=100*length(unique(barcode1))/N_assoc_BC),keyby=library]
# ~85% for DNA
# ~95% for RNA

# library    N_BC recovery
# 1: Lib3 (HepG2-DNA) 2531669 85.31161
# 2: Lib4 (Calu3-DNA) 2527482 85.17052
# 3: Lib2 (Calu3-RNA) 2832875 95.46158
# 4: Lib1 (HepG2-RNA) 2807621 94.61058
# UMI_perBC_lib[,.(N_BC=length(unique(barcode1)),recovery=100*length(unique(barcode1))/N_assoc_BC),keyby=library]
# library    N_BC recovery
# 1: Lib3 (HepG2-DNA) 5968825 201.1361
# 2: Lib1 (HepG2-RNA) 8715079 293.6788
# 3: Lib4 (Calu3-DNA) 5963932 200.9712
# 4: Lib2 (Calu3-RNA) 8702266 293.2470

## what percentage of oligo-associated BC are found in each index
# UMI_perBC_annot[!is.na(oligo),.(N_BC=length(unique(barcode1)),recovery=100*length(unique(barcode1))/N_assoc_BC),keyby=.(library,index)]
# ~57-79% for DNA
# ~85-92% for RNA

# library      index    N_BC recovery
# 1: Lib1 (HepG2-RNA) ACTACTTGCT 2544148 85.73213
# 2: Lib1 (HepG2-RNA) ATCCAAGTTG 2526745 85.14568
# 3: Lib1 (HepG2-RNA) TTGAACTAAG 2697903 90.91333
# 4: Lib2 (Calu3-RNA) AGTAGGCGGA 2726437 91.87486
# 5: Lib2 (Calu3-RNA) CGCGCGACTT 2661685 89.69286
# 6: Lib2 (Calu3-RNA) GAATAGACGC 2625000 88.45666
# 7: Lib3 (HepG2-DNA) AAGTTGACGA 2034935 68.57278
# 8: Lib3 (HepG2-DNA) GAAGCTGAAG 1933468 65.15357
# 9: Lib3 (HepG2-DNA) TGGAATGACG 2088033 70.36207
# 10: Lib4 (Calu3-DNA) AGATCGAATA 1831304 61.71087
# 11: Lib4 (Calu3-DNA) CGTCCGGCAT 1719397 57.93985
# 12: Lib4 (Calu3-DNA) TGCCTTCCTA 2368369 79.80877



# > UMI_per_oligo_DNA_vs_RNA[which(is.na(nUMI_DNA)),]
#                            oligo celltype condition nUMI_DNA nUMI_RNA nBC_DNA nBC_RNA nUMI_DNA_tot nUMI_RNA_tot
#  1:   10:6490221_A_-_0_archaic    HepG2        NS       NA       12      NA       1     21692596     66784415
#  2:    10:7526003_G_+_0_Purged    HepG2        NS       NA        3      NA       1     21692596     66784415
#  3: 12:104319850_G_+_0_archaic    HepG2        NS       NA        5      NA       1     21692596     66784415
#  4: 12:125402872_A_-_0_archaic    HepG2        NS       NA       13      NA       1     21692596     66784415
#  5:  12:16589359_A_+_0_archaic    Calu3        NS       NA        6      NA       2     25953045     61417883
#  6:  12:99055672_T_+_0_archaic    Calu3        NS       NA        1      NA       1     25953045     61417883
#  7:  16:24093936_A_+_0_archaic    Calu3        NS       NA        5      NA       1     25953045     61417883
#  8:  16:71937979_A_-_0_archaic    Calu3        NS       NA        3      NA       1     25953045     61417883
#  9:  16:84914656_A_+_0_archaic    HepG2        NS       NA        4      NA       1     21692596     66784415
# 10:   19:44252098_T_-_0_Purged    Calu3        NS       NA        1      NA       1     25953045     61417883
# 11: 1:57202008_G_-_+50_archaic    HepG2        NS       NA        4      NA       1     21692596     66784415
# 12:   1:57229256_G_-_0_archaic    Calu3        NS       NA        7      NA       1     25953045     61417883
# 13:   1:89651039_A_-_0_archaic    HepG2        NS       NA        2      NA       1     21692596     66784415
# 14:  3:112673011_T_-_0_archaic    Calu3        NS       NA        9      NA       1     25953045     61417883
# 15:  3:132031058_T_+_0_archaic    HepG2        NS       NA        4      NA       1     21692596     66784415
# 16:     3:161739090_A_-_0_MPRA    Calu3        NS       NA        5      NA       1     25953045     61417883
# 17:   3:46591625_A_-_0_archaic    HepG2        NS       NA        4      NA       1     21692596     66784415
# 18:  4:114645255_T_-_0_archaic    Calu3        NS       NA        9      NA       1     25953045     61417883
# 19:   4:38794041_C_-_0_archaic    HepG2        NS       NA       15      NA       2     21692596     66784415
# 20:   5:74674554_A_-_0_archaic    Calu3        NS       NA        4      NA       1     25953045     61417883
# 21:   5:74730820_C_-_0_archaic    HepG2        NS       NA        1      NA       1     21692596     66784415
# 22:    6:10839759_A_+_0_Purged    HepG2        NS       NA       14      NA       1     21692596     66784415
# 23:  6:133032714_C_-_0_archaic    HepG2        NS       NA        1      NA       1     21692596     66784415
# 24:   6:32359507_A_-_0_archaic    Calu3        NS       NA       14      NA       4     25953045     61417883
# 25:   6:32367046_T_-_0_archaic    Calu3        NS       NA       12      NA       1     25953045     61417883
# 26:   6:32570880_T_+_0_archaic    HepG2        NS       NA        2      NA       1     21692596     66784415
# 27:  7:105175068_C_+_0_archaic    Calu3        NS       NA       23      NA       1     25953045     61417883
#                            oligo celltype condition nUMI_DNA nUMI_RNA nBC_DNA nBC_RNA nUMI_DNA_tot nUMI_RNA_tot
# > UMI_per_oligo_DNA_vs_RNA[which(is.na(nUMI_RNA)),]
#                            oligo celltype condition nUMI_DNA nUMI_RNA nBC_DNA nBC_RNA nUMI_DNA_tot nUMI_RNA_tot
# 1:   12:99055672_T_+_0_archaic    HepG2        NS        2       NA       1      NA     21692596     66784415
# 2: 16:22310315_A_+_+50_archaic    HepG2        NS        2       NA       1      NA     21692596     66784415
# 3:     8:67005490_T_-_0_Purged    HepG2        NS        2       NA       1      NA     21692596     66784415



# oligo_activity_obs[,.(.N,Pct_up=mean(OR_GCadj>1 & FDR_scrambled<.01,na.rm=T),Pct_dwn=mean(OR_GCadj<1 & FDR_scrambled<.01,na.rm=T)),by=.(type_simple,SETUP_ID)]
# oligo_activity_obs[type=='scrambled' & OR_GCadj>2,][order(-OR_GCadj)]
#### check with
## https://meme-suite.org/meme/tools/meme
# UMI_per_oligo_DNA_vs_RNA_seq[type=='scrambled' & p.value>0.01,paste0('>',sequence)]

# hot question 1: is there any correlation between FE & GC content ?
# oligo_activity_obs[type=='scrambled',cor.test(OR,GC),by=SETUP_ID] # not really
cat("\nAll done!\n")
q("no")



