MPRA_DIR <- "/pasteur/helix/projects/evo_immuno_pop/MPRA"
ANALYSIS_DIR <- "MPRA_count_exp6_analysis3"
SCRIPT_DIR <- sprintf("%s/scripts/%s/", MPRA_DIR, ANALYSIS_DIR)

source(sprintf("%s/scripts/%s/00_starter.R", MPRA_DIR, ANALYSIS_DIR))
source(sprintf("%s/scripts/%s/00_load_associations.R", MPRA_DIR, ANALYSIS_DIR))
source(sprintf("%s/scripts/misc_plots.R", MPRA_DIR))

RUN_ID <- "RUN3_Z2_nBC10"
CRITERIA_ACTIVE <- "FDR5_scrambled5pct_FC0.2_GCnorm"
CRITERIA_EMVARS <- "EmVar_FDR5_FC.2"
CRITERIA_DIFF <- "Diff_FDR5_FC.2"
CRITERIA_EMVARS_DIFF <- "EmVarDiff_FDR5_FC.2"
FILTER_ACTIVE <- TRUE
FILTER_DIFF <- FALSE
FILTER_EMVAR <- TRUE
STIM_ONLY <- FALSE

cmd <- commandArgs(trailingOnly = TRUE)
print(cmd)

for (i in seq_along(cmd)) {
  if (cmd[i] == "--run_dir" || cmd[i] == "-r") {
    RUN_ID <- cmd[i + 1]
  }
  if (cmd[i] == "--criteria_active" || cmd[i] == "-a") {
    CRITERIA_ACTIVE <- cmd[i + 1]
  }
  if (cmd[i] == "--criteria_diff" || cmd[i] == "-d") {
    CRITERIA_DIFF <- cmd[i + 1]
  }
  if (cmd[i] == "--criteria_emvar" || cmd[i] == "-e") {
    CRITERIA_EMVARS <- cmd[i + 1]
  }
  if (cmd[i] == "--criteria_emvar_diff" || cmd[i] == "-ed") {
    CRITERIA_EMVARS_DIFF <- cmd[i + 1]
  }
  if (cmd[i] == "--filter_active" || cmd[i] == "-fa") {
    FILTER_ACTIVE <- as.logical(cmd[i + 1])
  }
  if (cmd[i] == "--filter_diff" || cmd[i] == "-fd") {
    FILTER_DIFF <- as.logical(cmd[i + 1])
  }
  if (cmd[i] == "--filter_emvar" || cmd[i] == "-fe") {
    FILTER_EMVAR <- as.logical(cmd[i + 1])
  }
  if (cmd[i] == "--stim_only" || cmd[i] == "-s") {
    STIM_ONLY <- as.logical(cmd[i + 1])
  }
}

if (!FILTER_ACTIVE) {
  CRITERIA_ACTIVE_OUT <- paste0("noActivityFilter_", CRITERIA_ACTIVE)
} else {
  CRITERIA_ACTIVE_OUT <- CRITERIA_ACTIVE
}

if (!FILTER_DIFF) {
  CRITERIA_DIFF_OUT <- paste0("noDiffFilter_", CRITERIA_DIFF)
} else {
  CRITERIA_DIFF_OUT <- CRITERIA_DIFF
}

if (!FILTER_EMVAR) {
  CRITERIA_EMVARS_OUT <- paste0("noEmVarFilter_", CRITERIA_EMVARS)
} else {
  CRITERIA_EMVARS_OUT <- CRITERIA_EMVARS
}
STIM_ONLY_OUT <- ifelse(STIM_ONLY, "stimOnly", "all")

tic("loading oligo & SNP annotations")
# load annotation of oligos (beforz: CRS)
oligo_source <- fread(sprintf("%s/data/%s/00_oligo_annot_v2.txt", MPRA_DIR, ANALYSIS_DIR))
# load annotations of SNPs
SNP_annot_v4 <- fread(sprintf("%s/data/%s/00_SNP_annot_v4.txt", MPRA_DIR, ANALYSIS_DIR))
# load annotation per POP
selected_annot <- fread(sprintf("%s/data/%s/selected_snp_annotation.tsv.gz", MPRA_DIR, ANALYSIS_DIR))
# load clean annotation of introgression
selected_annot_wide <- fread(sprintf("%s/data/%s/selected_snp_annotation_wide.tsv.gz", MPRA_DIR, ANALYSIS_DIR))
selected_annot_wide[, Introgression_source_top := Introgression_source_top_initial]
selected_annot_wide[Introgression_source_top == "both", Introgression_source_top := "Vindija/Denisova"]
selected_annot_wide[Introgression_source_top == "", Introgression_source_top := "Undetermined"]

SNP_annot_v5 <- merge(SNP_annot_v4[, -"Introgression_scenario"], selected_annot_wide, by = c("ID", "posID"))
toc()

# load oligo annotations
source(sprintf("%s/scripts/%s/02z__define_includedCRS.R", MPRA_DIR, ANALYSIS_DIR))


##### active parameters
source(sprintf("%s/scripts/%s/03_00_parameter_activity.R", MPRA_DIR, ANALYSIS_DIR))
##### diff parameters
source(sprintf("%s/scripts/%s/04_00_parameter_diff_activity.R", MPRA_DIR, ANALYSIS_DIR))
##### emVars parameters
source(sprintf("%s/scripts/%s/05_00_parameter_emVars.R", MPRA_DIR, ANALYSIS_DIR))
##### emVars diff parameters
source(sprintf("%s/scripts/%s/06_00_parameter_emVars_Diff.R", MPRA_DIR, ANALYSIS_DIR))

IN_DIR <- sprintf("%s/data/%s/03a_aggMPRA_analyse/%s/", MPRA_DIR, ANALYSIS_DIR, RUN_ID)
ACTIVE_DIR <- sprintf("%s/data/%s/03b_computeFDR/%s/Activity/%s/", MPRA_DIR, ANALYSIS_DIR, RUN_ID, CRITERIA_ACTIVE)
DIFF_DIR <- sprintf("%s/data/%s/03b_computeFDR/%s/Diff/%s/%s/", MPRA_DIR, ANALYSIS_DIR, RUN_ID, CRITERIA_DIFF, CRITERIA_ACTIVE_OUT)
EMVAR_DIR <- sprintf("%s/data/%s/03b_computeFDR/%s/emVars/%s/%s", MPRA_DIR, ANALYSIS_DIR, RUN_ID, CRITERIA_EMVARS, CRITERIA_ACTIVE_OUT)
EMVARDIFF_DIR <- sprintf("%s/data/%s/03b_computeFDR/%s/EmVar_Diff/%s/%s/%s/%s/%s", MPRA_DIR, ANALYSIS_DIR, RUN_ID, CRITERIA_EMVARS_DIFF, CRITERIA_EMVARS_OUT, CRITERIA_DIFF_OUT, CRITERIA_ACTIVE_OUT, STIM_ONLY_OUT)

dir.create(EMVARDIFF_DIR, recursive = TRUE)


FIGURE_EMVAR_DIFF_DIR <- sprintf("%s/figures/%s/%s/05_emVarDiff/%s/%s/%s/%s/%s", MPRA_DIR, ANALYSIS_DIR, RUN_ID, CRITERIA_EMVARS_DIFF, CRITERIA_EMVARS_OUT, CRITERIA_DIFF_OUT, CRITERIA_ACTIVE_OUT, STIM_ONLY_OUT)
POWER_EMVAR_DIFF_DIR <- sprintf("%s/power/", FIGURE_DIR)
dir.create(POWER_EMVAR_DIFF_DIR, recursive = TRUE)

FIGURE_EMVAR_DIR <- sprintf("%s/figures/%s/%s/04_emVars/%s/%s/", MPRA_DIR, ANALYSIS_DIR, RUN_ID, CRITERIA_EMVARS, CRITERIA_ACTIVE_OUT)
POWER_EMVAR_DIR <- sprintf("%s/power/", FIGURE_EMVAR_DIR)
dir.create(POWER_EMVAR_DIR, recursive = TRUE)

# dir.create(FIGURE_DIR, recursive = TRUE)
# dir.create(sprintf('%s/SupTables/', FIGURE_DIR, recursive = TRUE))

# emVARs_all_results <- fread(sprintf("%s/all_emVars_results.tsv.gz", IN_DIR), sep = "\t")
# emVARs_all_results <- emVARs_all_results[crsID %in% tested_and_ctrl_oligos_final[type=='tested',posID],]

# emVARs_Diff_all_results <- fread(sprintf("%s/all_emVars_Diff_results.tsv.gz", IN_DIR), sep = "\t")
# emVARs_Diff_all_results <- emVARs_Diff_all_results[crsID %in% tested_and_ctrl_oligos_final[type=='tested',posID], ]

oligo_activity <- fread(file = sprintf("%s/all_oligo_results.tsv.gz", IN_DIR), sep = "\t")

emVARs_obs <- fread(sprintf("%s/all_emVars_annotated__%s.tsv.gz", EMVAR_DIR, CRITERIA_ACTIVE_OUT), sep = "\t")
emVARs_obs_ctc <- emVARs_obs[ANALYSIS_SUBTYPE == "celltype_cond", ]
emVARs_obs_ctc <- merge(emVARs_obs_ctc, condition_summary[, -"celltype"], by.x = "ANALYSIS_NAME", by.y = "analysis_name")

FDR_TH <- emVARs_obs_ctc[FDR < .05, .(FDR_threshold = max(pval_emp)), by = .(celltype, condition)]

emVARs_power <- fread(sprintf("%s/all_emVars_power__%s.tsv.gz", EMVAR_DIR, CRITERIA_ACTIVE_OUT), sep = "\t")
emVARs_power_ctc <- emVARs_power[power == 1 & ANALYSIS_SUBTYPE == "celltype_cond", ]
emVARs_power_ctc <- merge(emVARs_power_ctc, condition_summary[, -"celltype"], by.x = "ANALYSIS_NAME", by.y = "analysis_name")


true_log2FC <- merge(emVARs_power_ctc[perm == "perm_0_0_0" & power == 1, .(perm, power, boot, crsID, oligo1, oligo2, ANALYSIS_NAME)], oligo_activity[perm == "perm_0_0_0" & power == 0, .(oligo1 = oligo, alpha1 = alpha, ANALYSIS_NAME)], by = c("oligo1", "ANALYSIS_NAME"), all.x = TRUE)
true_log2FC <- merge(true_log2FC, oligo_activity[perm == "perm_0_0_0" & power == 0, .(oligo2 = oligo, alpha2 = alpha, ANALYSIS_NAME)], by = c("oligo2", "ANALYSIS_NAME"), all.x = TRUE)
true_log2FC[, true_log2FC := log2(alpha2) - log2(alpha1)]
emVARs_power_ctc <- merge(emVARs_power_ctc, true_log2FC[, .(crsID, true_a1 = alpha1, true_a2 = alpha2, true_log2FC, ANALYSIS_NAME)], by = c("crsID", "ANALYSIS_NAME"), allow.cartesian = TRUE)

nBC_breaks <- c(0, 100, 200, 800, 1600, Inf)
nBC_breaks <- c(0, 100, 300, 500, 1000, Inf)
nUMIperBC_breaks <- c(0, 2, 4, Inf)

# celltype data
FigData <- list(H1 = emVARs_power_ctc[condition == "NS" & perm == "perm_0_0_0", ], H0 = emVARs_power_ctc[condition == "NS" & perm != "perm_0_0_0", ])
FigData <- rbindlist(FigData, idcol = "data_source")

p <- ggplot(FigData, aes(x = pval.LRT, fill = celline, alpha = data_source)) +
  geom_histogram(position = "dodge")
p <- p + facet_grid(cols = vars(celline)) + guides(fill = "none")
p <- p + scale_fill_manual(values = color_celline)
p <- p + scale_alpha_manual(values = c("H1" = .4, "H0" = .7))
p <- p + theme_plot(lpos = "top", fontsize = 12, rotate.x = 90) + xlab("Differential activity p-values (emVar)")

pdf(sprintf("%s/1a__emVar_pvalues__H1_vs_H0.pdf", POWER_EMVAR_DIR), height = 2.5, width = 4)
print(p)
dev.off()

p <- ggplot(FigData, aes(x = pval.LRT, fill = celline, alpha = data_source)) +
  geom_histogram(position = "dodge")
p <- p + facet_grid(rows = vars(celline)) + guides(fill = "none")
p <- p + scale_fill_manual(values = color_celline)
p <- p + scale_alpha_manual(values = c("H1" = .4, "H0" = .7))
p <- p + theme_plot(lpos = "top", fontsize = , rotate.x = 90) + xlab("Differential activity\n p-values (emVar)")

pdf(sprintf("%s/1a__emVar_pvalues__H1_vs_H0_VERTICAL.pdf", POWER_EMVAR_DIR), height = 3, width = 2)
print(p)
dev.off()


FigData <- list(H1 = emVARs_power_ctc[condition == "NS" & perm == "perm_0_0_0", ], H0 = emVARs_power_ctc[condition == "NS" & perm != "perm_0_0_0", ])
FigData <- rbindlist(FigData, idcol = "data_source")

p <- ggplot(FigData, aes(x = pval_emp, fill = celline, alpha = data_source)) +
  geom_histogram(position = "dodge", breaks = seq(0, 1, l = 31))
p <- p + facet_grid(cols = vars(celline)) + guides(fill = "none")
p <- p + scale_fill_manual(values = color_celline)
p <- p + scale_alpha_manual(values = c("H1" = .4, "H0" = .7))
p <- p + theme_plot(lpos = "top", fontsize = 12, rotate.x = 90) + xlab("Differential activity p-values (emVar, empirical)")

pdf(sprintf("%s/1a__emVar_pval_emp__H1_vs_H0.pdf", POWER_EMVAR_DIR), height = 2.5, width = 4)
print(p)
dev.off()

pdf(sprintf("%s/1b__emVar_pval_emp_by_nBC.pdf", POWER_EMVAR_DIR), height = 4, width = 6)
p <- ggplot(FigData, aes(x = pval_emp, fill = celline, alpha = data_source)) +
  geom_histogram(aes(y = after_stat(density)), position = "dodge", breaks = seq(0, 1, l = 21))
p <- p + facet_grid(rows = vars(celline), cols = vars(cut((nBCs_g1a1 + nBCs_g1a2) / 2, nBC_breaks / 2))) + guides(fill = "none")
p <- p + scale_fill_manual(values = color_celline)
p <- p + scale_alpha_manual(values = c("H1" = .4, "H0" = .7))
p <- p + theme_plot(lpos = "top", fontsize = 12, rotate.x = 90) + xlab("Differential activity p-values (emVar, empirical)")

# p <- p + facet_grid( ~ ifelse(perm=='perm_0_0_0','H1','H0') + cut((nBCs_g1a1+nBCs_g1a2)/2,nBC_breaks/2),scale='free',ncol=5)
# p <- p + theme_plot(lpos='right') + xlab('emVar Pvalue')
# p <- p + guides(color=guide_legend(ncol=1),linetype=guide_legend(ncol=1))+labs(color="nBC per CRS",linetype="nUMI per BC") + theme(legend.title=element_text(size=6))
print(p)
dev.off()


pdf(sprintf("%s/1b__emVar_pvalues_by_nBC.pdf", POWER_EMVAR_DIR), height = 4, width = 6)
p <- ggplot(FigData, aes(x = pval.LRT, fill = celline, alpha = data_source)) +
  geom_histogram(aes(y = after_stat(density)), position = "dodge", breaks = seq(0, 1, l = 21))
p <- p + facet_grid(rows = vars(celline), cols = vars(cut((nBCs_g1a1 + nBCs_g1a2) / 2, nBC_breaks / 2))) + guides(fill = "none")
p <- p + scale_fill_manual(values = color_celline)
p <- p + scale_alpha_manual(values = c("H1" = .4, "H0" = .7))
p <- p + theme_plot(lpos = "top", fontsize = 12, rotate.x = 90) + xlab("Differential activity p-values (emVar, empirical)")

# p <- p + facet_grid( ~ ifelse(perm=='perm_0_0_0','H1','H0') + cut((nBCs_g1a1+nBCs_g1a2)/2,nBC_breaks/2),scale='free',ncol=5)
# p <- p + theme_plot(lpos='right') + xlab('emVar Pvalue')
# p <- p + guides(color=guide_legend(ncol=1),linetype=guide_legend(ncol=1))+labs(color="nBC per CRS",linetype="nUMI per BC") + theme(legend.title=element_text(size=6))
print(p)
dev.off()

# nBC_breaks=c(0,10,50,100,150,200,300,400,600,800,1000,1500,2000,Inf)
# nUMIperBC_breaks=c(0,1,2,3,5,Inf)
# logFC_breaks=c(0,.1,.2,.3,0.4,.5,.6,.7,.8,.9,1,2,Inf)

FigData <- list(H1 = emVARs_power_ctc[perm == "perm_0_0_0", ], H0 = emVARs_power_ctc[perm != "perm_0_0_0", ])
FigData <- rbindlist(FigData, idcol = "data_source")
FigData <- merge(FigData, FDR_TH, by = c("celltype", "condition"))


# power_estimate <- FigData[boot==1,.(.N,
#                                 nBC=mean(nBCs_g1a1+nBCs_g1a2),
#                                 logFC=mean(abs(true_log2FC)),
#                                 Power=mean(FDR<.05),
#                                 Power_P05=mean(pval_emp<.05)
#                                 ),keyby=.(perm,
#                                           nBC_group=cut((nBCs_g1a1+nBCs_g1a2)/2,nBC_breaks/2),
#                                           logFC_group=cut(abs(true_log2FC),logFC_breaks),
#                                           COND_ID)]

nBC_breaks <- c(0, 100, 300, 800, Inf)
logFC_breaks <- c(0, 0.01, 0.05, .1, 0.15, .2, 0.25, .3, 0.35, 0.4, .5, .6, .8, Inf)

nUMIperBC_breaks <- c(0, 2, 4, Inf)
power_estimate_collapsed <- FigData[boot == 1, .(.N,
  nBC = mean(nBCs_g1a1 + nBCs_g1a2),
  logFC = mean(abs(true_log2FC), na.rm = T),
  Power = mean(pval_emp < FDR_threshold & abs(log2FC_archaic_vs_modern) > LOG2FC_TH_EMVAR),
  Power_active = mean((pval_emp < FDR_threshold & abs(log2FC_archaic_vs_modern) > LOG2FC_TH_EMVAR)[is_tested_cre == TRUE]),
  Power_P05 = mean(pval_emp[is_tested_cre == TRUE] < .05)
), keyby = .(perm,
  nBC_group = cut((nBCs_g1a1 + nBCs_g1a2) / 2, nBC_breaks / 2),
  logFC_group = cut(abs(true_log2FC), logFC_breaks)
)]

power_estimate_by_cond <- FigData[boot == 1, .(.N,
  nBC = mean(nBCs_g1a1 + nBCs_g1a2),
  logFC = mean(abs(true_log2FC), na.rm = T),
  Power = mean(pval_emp < FDR_threshold & abs(log2FC_archaic_vs_modern) > LOG2FC_TH_EMVAR),
  Power_active = mean((pval_emp < FDR_threshold & abs(log2FC_archaic_vs_modern) > LOG2FC_TH_EMVAR)[is_tested_cre == TRUE]),
  Power_P05 = mean(pval_emp < .05),
  Power_active_P05 = mean(pval_emp[is_tested_cre == TRUE] < .05)
), keyby = .(perm,
  COND_ID,
  logFC_group = cut(abs(true_log2FC), logFC_breaks)
)]


### by nBC
pdf(sprintf("%s/2a_power_collapse_logFC_nBC_FDR5.pdf", POWER_EMVAR_DIR), height = 4, width = 4)
p <- ggplot(power_estimate_collapsed[N > 50, ], aes(x = logFC, y = Power * 100, color = nBC_group)) +
  geom_line(alpha = .5) +
  geom_point(size = 0.4, alpha = .5)
p <- p + theme_plot(lpos = "right") + facet_grid(rows = vars(ifelse(perm == "perm_0_0_0", "H1", "H0"))) + scale_color_viridis(discrete = TRUE, option = "B", direction = -1, end = .9)
p <- p + guides(color = guide_legend(ncol = 1), linetype = guide_legend(ncol = 1)) + labs(color = "nBC per CRS") + theme(legend.title = element_text(size = 6))
print(p)
dev.off()

pdf(sprintf("%s/2b_power_collapse_logFC_nBC_FDR5_active.pdf", POWER_EMVAR_DIR), height = 4, width = 4)
p <- ggplot(power_estimate_collapsed[N > 50, ], aes(x = logFC, y = Power_active * 100, color = nBC_group)) +
  geom_line(alpha = .5) +
  geom_point(size = 0.4, alpha = .5)
p <- p + theme_plot(lpos = "right") + facet_grid(rows = vars(ifelse(perm == "perm_0_0_0", "H1", "H0"))) + scale_color_viridis(discrete = TRUE, option = "B", direction = -1, end = .9)
p <- p + guides(color = guide_legend(ncol = 1), linetype = guide_legend(ncol = 1)) + labs(color = "nBC per CRS") + theme(legend.title = element_text(size = 6))
print(p)
dev.off()

pdf(sprintf("%s/2c_power_collapse_logFC_nBC_P05.pdf", POWER_EMVAR_DIR), height = 4, width = 4)
p <- ggplot(power_estimate_collapsed[N > 50, ], aes(x = logFC, y = Power_P05 * 100, color = nBC_group)) +
  geom_line(alpha = .5) +
  geom_point(size = 0.4, alpha = .5)
p <- p + theme_plot(lpos = "right") + facet_grid(rows = vars(ifelse(perm == "perm_0_0_0", "H1", "H0"))) + scale_color_viridis(discrete = TRUE, option = "B", direction = -1, end = .9) + ylab("Power  /  type I error")
p <- p + guides(color = guide_legend(ncol = 1), linetype = guide_legend(ncol = 1)) + labs(color = "nBC per CRS", linetype = "nUMI per BC") + theme(legend.title = element_text(size = 6))
print(p)
dev.off()


### by condition

pdf(sprintf("%s/3a_power_byCOND_logFC_nBC_FDR5.pdf", POWER_EMVAR_DIR), height = 4, width = 4)
p <- ggplot(power_estimate_by_cond[N > 50, ], aes(x = logFC, y = Power * 100, color = COND_ID)) +
  geom_line(alpha = .5) +
  geom_point(size = 0.4, alpha = .5)
p <- p + theme_plot(lpos = "right") + facet_grid(rows = vars(ifelse(perm == "perm_0_0_0", "H1", "H0"))) + scale_color_manual(values = color_setup_simplified_norep) + ylab("Power  /  type I error")
p <- p + guides(color = guide_legend(ncol = 1), linetype = guide_legend(ncol = 1)) + labs(color = "nBC per CRS", linetype = "nUMI per BC") + theme(legend.title = element_text(size = 6))
print(p)
dev.off()

pdf(sprintf("%s/3b_power_byCOND_logFC_nBC_FDR5_active.pdf", POWER_EMVAR_DIR), height = 4, width = 4)
p <- ggplot(power_estimate_by_cond[N > 50, ], aes(x = logFC, y = Power_active * 100, color = COND_ID)) +
  geom_line(alpha = .5) +
  geom_point(size = 0.4, alpha = .5)
p <- p + theme_plot(lpos = "right") + facet_grid(rows = vars(ifelse(perm == "perm_0_0_0", "H1", "H0"))) + scale_color_manual(values = color_setup_simplified_norep) + ylab("Power  /  type I error")
p <- p + guides(color = guide_legend(ncol = 1), linetype = guide_legend(ncol = 1)) + labs(color = "nBC per CRS", linetype = "nUMI per BC") + theme(legend.title = element_text(size = 6))
print(p)
dev.off()

pdf(sprintf("%s/3c_power_byCOND_logFC_nBC_P05.pdf", POWER_EMVAR_DIR), height = 4, width = 4)
p <- ggplot(power_estimate_by_cond[N > 50, ], aes(x = logFC, y = Power_P05 * 100, color = COND_ID)) +
  geom_line(alpha = .5) +
  geom_point(size = 0.4, alpha = .5)
p <- p + theme_plot(lpos = "right") + facet_grid(rows = vars(ifelse(perm == "perm_0_0_0", "H1", "H0"))) + scale_color_manual(values = color_setup_simplified_norep) + ylab("Power  /  type I error")
p <- p + guides(color = guide_legend(ncol = 1), linetype = guide_legend(ncol = 1)) + labs(color = "nBC per CRS") + theme(legend.title = element_text(size = 6))
print(p)
dev.off()


Fig_data_power <- melt(power_estimate_collapsed[perm == "perm_0_0_0"], measure.vars = c("Power", "Power_P05"), variable.name = "threshold")
Fig_data_power[threshold == "Power", threshold := "5% FDR"]
Fig_data_power[threshold == "Power_P05", threshold := "P_emp < 0.05"]

pdf(sprintf("%s/4a_power_collapse_logFC_nBC_double_threshold.pdf", POWER_EMVAR_DIR), height = 2.5, width = 4)
p <- ggplot(Fig_data_power[N > 50, ], aes(x = logFC, y = value * 100, color = nBC_group, linetype = threshold))
p <- p + geom_line(alpha = .5) + geom_point(size = 0.4, alpha = .5)
p <- p + scale_color_viridis(discrete = TRUE, option = "B", direction = -1, end = .9)
p <- p + theme_plot(lpos = "right") + ylab("Power") + xlab("|log2 FC|")
p <- p + guides(color = guide_legend(ncol = 1), linetype = guide_legend(ncol = 1))
p <- p + labs(color = "nBC per CRS", linetype = "threshold") + theme(legend.title = element_text(size = 6))
print(p)
dev.off()



Fig_data_power <- melt(power_estimate_collapsed[perm == "perm_0_0_0"], measure.vars = c("Power_active", "Power_P05"), variable.name = "threshold")
Fig_data_power[threshold == "Power_active", threshold := "5% FDR"]
Fig_data_power[threshold == "Power_P05", threshold := "P_emp < 0.05"]

pdf(sprintf("%s/4b_power_collapse_logFC_nBC_double_threshold_activeOnly.pdf", POWER_EMVAR_DIR), height = 2.5, width = 4)
p <- ggplot(Fig_data_power[N > 50, ], aes(x = logFC, y = value * 100, color = nBC_group, linetype = threshold))
p <- p + geom_line(alpha = .5) + geom_point(size = 0.4, alpha = .5)
p <- p + scale_color_viridis(discrete = TRUE, option = "B", direction = -1, end = .9)
p <- p + theme_plot(lpos = "right") + ylab("Power") + xlab("|log2 FC|")
p <- p + guides(color = guide_legend(ncol = 1), linetype = guide_legend(ncol = 1))
p <- p + labs(color = "nBC per CRS", linetype = "threshold") + theme(legend.title = element_text(size = 6))
print(p)
dev.off()

# by condition
Fig_data_power <- melt(power_estimate_by_cond[perm == "perm_0_0_0"], measure.vars = c("Power", "Power_P05"), variable.name = "threshold")
Fig_data_power[threshold == "Power", threshold := "5% FDR"]
Fig_data_power[threshold == "Power_P05", threshold := "P_emp < 0.05"]

pdf(sprintf("%s/4c_power_byCOND_logFC_nBC_double_threshold.pdf", POWER_EMVAR_DIR), height = 2.5, width = 4)
p <- ggplot(Fig_data_power[N > 50, ], aes(x = logFC, y = value * 100, color = COND_ID, linetype = threshold))
p <- p + geom_line(alpha = .5) + geom_point(size = 0.4, alpha = .5)
p <- p + scale_color_manual(values = color_setup_simplified_norep)
p <- p + theme_plot(lpos = "right") + ylab("Power") + xlab("|log2 FC|")
p <- p + guides(color = guide_legend(ncol = 1), linetype = guide_legend(ncol = 1))
p <- p + labs(color = "condition", linetype = "threshold") + theme(legend.title = element_text(size = 6))
print(p)
dev.off()


Fig_data_power <- melt(power_estimate_by_cond[perm == "perm_0_0_0"], measure.vars = c("Power_active", "Power_P05"), variable.name = "threshold")
Fig_data_power[threshold == "Power_active", threshold := "5% FDR"]
Fig_data_power[threshold == "Power_P05", threshold := "P_emp < 0.05"]

pdf(sprintf("%s/4d_power_byCOND_logFC_nBC_double_threshold_activeOnly.pdf", POWER_EMVAR_DIR), height = 2.5, width = 4)
p <- ggplot(Fig_data_power[N > 50, ], aes(x = logFC, y = value * 100, color = COND_ID, linetype = threshold))
p <- p + geom_line(alpha = .5) + geom_point(size = 0.4, alpha = .5)
p <- p + scale_color_manual(values = color_setup_simplified_norep)
p <- p + theme_plot(lpos = "right") + ylab("Power") + xlab("|log2 FC|")
p <- p + guides(color = guide_legend(ncol = 1), linetype = guide_legend(ncol = 1))
p <- p + labs(color = "condition", linetype = "threshold") + theme(legend.title = element_text(size = 6))
print(p)
dev.off()


Pct_nBC <- emVARs_obs_ctc[, .N, by = .(nBC = cut((nBCs_g1a1 + nBCs_g1a2) / 2, nBC_breaks / 2))][, Pct := N / sum(N)][1:.N, ]
# make a pie chart
Pct_nBC[, label := paste0("  ", round(Pct * 100, ifelse(Pct > 0.02, 0, 1)), "%")]
# N_Context[, label := c(paste0("  ", round(Pct * 100), "%")[Pct > .03], "", "2.4%", "")]

Pct_nBC[, x_pos := 1.65]
p <- ggplot(Pct_nBC, aes(x = "", y = Pct, fill = nBC)) +
  geom_bar(stat = "identity", width = 1, col = "black")
p <- p + coord_polar("y", start = 0) + scale_fill_viridis(discrete = TRUE, option = "B", direction = -1, end = .9) + theme_void()
p <- p + labs(fill = "Average number\nof barcodes\nper allele")
p <- p + geom_text(aes(x = x_pos, label = label), position = position_stack(vjust = 0.5)) + theme(legend.position = "bottom")
pdf(sprintf("%s/05c_piechart_number_of_barcodes_per_snp.pdf", POWER_EMVAR_DIR), height = 5, width = 4)
print(p)
dev.off()


cat("All done!\n")
q("no")

# power_estimate_collapse_boot=results_power_annot[,.(.N,
#                                 nBC=mean(nBC_CRS1+nBC_CRS2),
#                                 logFC=mean(abs(true_logFC_2vs1)),
#                                 Power=mean(pval.LRT<.01)
#                                 ),keyby=.(perm,
#                                           nBC_group=cut((nBC_CRS1+nBC_CRS2)/2,nBC_breaks/2),
#                                           logFC_group=cut(abs(true_logFC_2vs1),logFC_breaks),
#                                           nUMI_RNA=factor(boot,levels=c(1,2),labels=c('current','2x increase')))]


# pdf(sprintf('%s/power_collapse_logFC_nBC_nUMIperBC_RNA.pdf',POWER_DIR),height=2.5,width=4)
# p <- ggplot(power_estimate_collapse_boot[N>100 & perm==0,],aes(x=logFC,y=Power*100,color=nBC_group,linetype=nUMI_RNA)) + geom_line(alpha=.5) + geom_point(size=0.4,alpha=.5)
# p <- p + theme_plot(lpos='right') + facet_grid(rows=vars(ifelse(perm=='perm0','H1','H0'))) + scale_color_viridis(discrete=TRUE,option = "B",direction=-1,end=.9) + ylab('Power  /  type I error')
# p <- p + guides(color=guide_legend(ncol=1),linetype=guide_legend(ncol=1))+labs(color="nBC per CRS",linetype="nUMI per BC") + theme(legend.title=element_text(size=6))
# print(p)
# dev.off()
