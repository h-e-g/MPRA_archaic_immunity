# 04b_response_Diff.R

# running: sbatch -p geh,common --mem=30G  00_Rscript.sh MPRA_count_exp6_analysisZ/04b_upset_plot_oligo_Diff.R

MPRA_DIR <- "/pasteur/helix/projects/evo_immuno_pop/MPRA"
ANALYSIS_DIR <- "MPRA_count_exp6_analysis3"

source(sprintf("%s/scripts/%s/00_starter.R", MPRA_DIR, ANALYSIS_DIR))
source(sprintf("%s/scripts/%s/00_load_associations.R", MPRA_DIR, ANALYSIS_DIR))
source(sprintf("%s/scripts/misc_plots.R", MPRA_DIR))

RUN_ID <- "RUN3_Z2_nBC10"
CRITERIA_ACTIVE <- "FDR5_scrambled5pct_FC0.2_GCnorm"
CRITERIA_DIFF <- "Diff_FDR5_FC.2"
FILTER_ACTIVE <- TRUE

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
  if (cmd[i] == "--filter_active" || cmd[i] == "-fa") {
    FILTER_ACTIVE <- as.logical(cmd[i + 1])
  }
}

if (!FILTER_ACTIVE) {
  CRITERIA_ACTIVE_OUT <- paste0("noActivityFilter_", CRITERIA_ACTIVE)
} else {
  CRITERIA_ACTIVE_OUT <- CRITERIA_ACTIVE
}


tic("loading oligo & SNP annotations")
# load annotation of oligos (beforz: CRS)
oligo_source <- fread(sprintf("%s/data/%s/00_oligo_annot_v2.txt", MPRA_DIR, ANALYSIS_DIR))

# load annotations of SNPs
SNP_annot <- fread(sprintf("%s/data/%s/00_SNP_annot_v4.txt", MPRA_DIR, ANALYSIS_DIR))
# load annotation per POP
selected_annot <- fread(sprintf("%s/data/%s/selected_snp_annotation.tsv.gz", MPRA_DIR, ANALYSIS_DIR))
toc()

# source(sprintf("%s/scripts/%s/03_00_parameter_activity.R", MPRA_DIR, ANALYSIS_DIR))
source(sprintf("%s/scripts/%s/04_00_parameter_diff_activity.R", MPRA_DIR, ANALYSIS_DIR))

IN_DIR <- sprintf("%s/data/%s/03a_aggMPRA_analyse/%s/", MPRA_DIR, ANALYSIS_DIR, RUN_ID)
ACTIVE_DIR <- sprintf("%s/data/%s/03b_computeFDR/%s/Activity/%s/", MPRA_DIR, ANALYSIS_DIR, RUN_ID, CRITERIA_ACTIVE)
DIFF_DIR <- sprintf("%s/data/%s/03b_computeFDR/%s/Diff/%s/%s/", MPRA_DIR, ANALYSIS_DIR, RUN_ID, CRITERIA_DIFF, CRITERIA_ACTIVE_OUT)

dir.create(DIFF_DIR, recursive = TRUE)
FIGURE_DIR <- sprintf("%s/figures/%s/%s/03_oligo_diff/%s/%s", MPRA_DIR, ANALYSIS_DIR, RUN_ID, CRITERIA_DIFF, CRITERIA_ACTIVE_OUT)
dir.create(FIGURE_DIR, recursive = TRUE)

dir.create(sprintf("%s/SupTables/", FIGURE_DIR), recursive = TRUE)

# load oligo annotations
source(sprintf("%s/scripts/%s/02z__define_includedCRS.R", MPRA_DIR, ANALYSIS_DIR))

# library(ggsignif)

#################################################################
####### load diffential activity results         ################
#################################################################

##### observed data #####
oligo_activity_Diff_obs <- fread(sprintf("%s/all_oligos_diff_annotated__%s.tsv.gz", DIFF_DIR, CRITERIA_DIFF))

# response data
oligo_response_annot_obs <- oligo_activity_Diff_obs[ANALYSIS_SUBTYPE == "response", ]
oligo_response_annot_obs[, COND_ID := gsub("-ACE2", "", gsub("response_(.*)", "\\1", ANALYSIS_NAME))]
oligo_response_annot_obs <- merge(oligo_response_annot_obs, condition_summary, by = "COND_ID")
oligo_response_annot_obs <- merge(oligo_response_annot_obs, SNP_annot[, .(posID, POP_adaptive, Introgressed_from = Adaptive_from, Introgression_scenario_v2)], by = "posID")
oligo_response_annot_obs_archaic <- oligo_response_annot_obs[oligo %chin% tested_and_ctrl_oligos_final_annot[type == "tested", oligo], ]

# celltype comparison data (NS)
oligo_celltype_comp_NS_annot_obs <- oligo_activity_Diff_obs[grepl("celltype_comp", ANALYSIS_SUBTYPE) & grepl("_NS$", ANALYSIS_NAME), ]
oligo_celltype_comp_NS_annot_obs_archaic <- oligo_celltype_comp_NS_annot_obs[oligo %chin% tested_and_ctrl_oligos_final_annot[type == "tested", oligo], ]

# celltype comparison data (all)
oligo_celltype_comp_all_annot_obs <- oligo_activity_Diff_obs[grepl("celltype_comp", ANALYSIS_SUBTYPE) & grepl("_all$", ANALYSIS_NAME), ]
oligo_celltype_comp_all_annot_obs_archaic <- oligo_celltype_comp_NS_annot_obs[oligo %chin% tested_and_ctrl_oligos_final_annot[type == "tested", oligo], ]

##### permuted data #####
oligo_activity_Diff_perm <- fread(sprintf("%s/all_oligos_diff_annotated_perm__%s.tsv.gz", DIFF_DIR, CRITERIA_DIFF))

# response data
oligo_response_annot_perm <- oligo_activity_Diff_perm[ANALYSIS_SUBTYPE == "response", ]
oligo_response_annot_perm[, COND_ID := gsub("-ACE2", "", gsub("response_(.*)", "\\1", ANALYSIS_NAME))]
oligo_response_annot_perm <- merge(oligo_response_annot_perm, condition_summary, by = "COND_ID")
oligo_response_annot_perm <- merge(oligo_response_annot_perm, SNP_annot[, .(posID, POP_adaptive, Introgressed_from = Adaptive_from, Introgression_scenario_v2)], by = "posID")
oligo_response_annot_perm_archaic <- oligo_response_annot_perm[oligo %chin% tested_and_ctrl_oligos_final_annot[type == "tested", oligo], ]

# celltype comparison data (NS)
oligo_celltype_comp_NS_annot_perm <- oligo_activity_Diff_perm[grepl("celltype_comp", ANALYSIS_SUBTYPE) & grepl("_NS$", ANALYSIS_NAME), ]
oligo_celltype_comp_NS_annot_perm_archaic <- oligo_celltype_comp_NS_annot_perm[oligo %chin% tested_and_ctrl_oligos_final_annot[type == "tested", oligo], ]

# celltype comparison data (all)
oligo_celltype_comp_all_annot_perm <- oligo_activity_Diff_perm[grepl("celltype_comp", ANALYSIS_SUBTYPE) & grepl("_all$", ANALYSIS_NAME), ]
oligo_celltype_comp_all_annot_perm_archaic <- oligo_celltype_comp_NS_annot_perm[oligo %chin% tested_and_ctrl_oligos_final_annot[type == "tested", oligo], ]

# target genes
oligo_target <- fread(sprintf("%s/data/%s/00_oligo_targets_v4.txt", MPRA_DIR, ANALYSIS_DIR), sep = "\t")
oligo_target_collapsed <- fread(sprintf("%s/data/%s/00_oligo_targets_collapsed_v4.txt", MPRA_DIR, ANALYSIS_DIR), sep = "\t")
crs_targets <- unique(oligo_target_collapsed[, -"oligo"])
NearestGene <- unique(oligo_target[method == "NearestGene" & gene_type != "", .(posID = crsID, crsID, TargetGene)])

oligo_activity_obs <- fread(sprintf("%s/all_oligos_annotated__%s.tsv.gz", ACTIVE_DIR, CRITERIA_ACTIVE))
oligo_activity_obs_NS=oligo_activity_obs[condition=='NS' & ANALYSIS_SUBTYPE=='celltype_cond',]


#########################################################
####### Pvalue distribution              ################
#########################################################

# response data
FigData <- list(observed = oligo_response_annot_obs_archaic, permuted = oligo_response_annot_perm_archaic)
FigData <- rbindlist(FigData, idcol = "data_source")

p <- ggplot(FigData, aes(x = pval.LRT, fill = COND_ID, alpha = data_source)) +
  geom_histogram(position = "dodge")
p <- p + facet_wrap(~factor(COND_ID,levels=names(color_setup_simplified_norep)), ncol = 3)
p <- p + scale_fill_manual(values = color_setup_simplified_norep) + scale_alpha_manual(values = c("observed" = .4, "permuted" = .7))
p <- p + theme_plot(lpos = "top", fontsize = 12, rotate.x = 90) + xlab("Differential activity p-values") + guides(fill = "none")

pdf(sprintf("%s/1a__response_pvalues_by_condition_observed_vs_permuted.pdf", FIGURE_DIR), height = 3.5, width = 4)
print(p)
dev.off()


# response data (empirical pvalues)
FigData <- list(observed = oligo_response_annot_obs_archaic, permuted = oligo_response_annot_perm_archaic)
FigData <- rbindlist(FigData, idcol = "data_source")

p <- ggplot(FigData, aes(x = pval_emp, fill = COND_ID, alpha = data_source)) +
  geom_histogram(position = "dodge", breaks = seq(0, 1, l = 31))
p <- p + facet_wrap(~~factor(COND_ID,levels=names(color_setup_simplified_norep)), ncol = 3)
p <- p + scale_fill_manual(values = color_setup_simplified_norep) + scale_alpha_manual(values = c("observed" = .4, "permuted" = .7))
p <- p + theme_plot(lpos = "top", fontsize = 12, rotate.x = 90) + xlab("Differential activity p-values (empirical)") + guides(fill = "none")

pdf(sprintf("%s/1a__response_pval_emp_by_condition__observed_vs_permuted.pdf", FIGURE_DIR), height = 3.5, width = 4)
print(p)
dev.off()

# celltype data
FigData <- list(observed = oligo_celltype_comp_NS_annot_obs_archaic, permuted = oligo_celltype_comp_NS_annot_perm_archaic)
FigData <- rbindlist(FigData, idcol = "data_source")
FigData[, comparison := gsub("-ACE2", "", paste(gsub("_NS_all", "", group1_labels), gsub("_NS_all", "", group2_labels), sep = "-"))]

p <- ggplot(FigData, aes(x = pval.LRT, fill = comparison, alpha = data_source)) +
  geom_histogram(position = "dodge")
p <- p + facet_grid(cols = vars(comparison)) + guides(fill = "none")
p <- p + scale_alpha_manual(values = c("observed" = .4, "permuted" = .7))
p <- p + theme_plot(lpos = "top", fontsize = 12, rotate.x = 90) + xlab("Differential activity p-values")

pdf(sprintf("%s/1b__cellline_comparison_pvalues__observed_vs_permuted.pdf", FIGURE_DIR), height = 2.5, width = 4)
print(p)
dev.off()


# celltype data (empirical pvalues)
FigData <- list(observed = oligo_celltype_comp_NS_annot_obs_archaic, permuted = oligo_celltype_comp_NS_annot_perm_archaic)
FigData <- rbindlist(FigData, idcol = "data_source")
FigData[, comparison := gsub("-ACE2", "", paste(gsub("_NS_all", "", group1_labels), gsub("_NS_all", "", group2_labels), sep = "-"))]

p <- ggplot(FigData, aes(x = pval_emp, fill = comparison, alpha = data_source))
p <- p + geom_histogram(position = "dodge", breaks = seq(0, 1, l = 31))
p <- p + facet_grid(cols = vars(comparison)) + guides(fill = "none")
p <- p + scale_alpha_manual(values = c("observed" = .4, "permuted" = .7))
p <- p + theme_plot(lpos = "top", fontsize = 12, rotate.x = 90) + xlab("Differential activity p-values (empirical)")

pdf(sprintf("%s/1b__cellline_comparison_pval_emp__observed_vs_permuted.pdf", FIGURE_DIR), height = 2.5, width = 4)
print(p)
dev.off()

# FigData <- list(observed = oligo_celltype_comp_all_annot_obs_archaic, permuted = oligo_celltype_comp_all_annot_perm_archaic)
# FigData <- rbindlist(FigData, idcol = "data_source")

# p <- ggplot(FigData, aes(x = pval.LRT, fill = celline, alpha=data_source)) + geom_histogram(stat='Identity', position = "dodge")
# p <- p + facet_grid(rows = vars(celline))
# p <- p + scale_fill_manual(values = color_celline) + scale_alpha_manual(values = c("observed" = .4, "permuted" = .7))
# p <- p + theme_plot(lpos = "right") + xlab("Differential activity p-values")

# pdf(sprintf("%s/1b__cellline_compatison_pvalues__observed_vs_permuted.pdf", FIGURE_DIR), height = 2.5, width = 3)
# print(p)
# dev.off()


###############################################################
################# count responsive CRSs #######################
###############################################################

### count responsive cre per POP
oligo_response_annot_perPOP <- merge(oligo_response_annot_obs_archaic, selected_annot, by = "posID", allow.cartesian = TRUE)
oligo_response_annot_perPOP <- merge(oligo_response_annot_obs_archaic, selected_annot[is_adaptive == TRUE, ], by = "posID", allow.cartesian = TRUE)

Nsig <- oligo_response_annot_perPOP[, .(N = length(unique(oligo)), N_sig = length(unique(oligo[oligo_diff_CRITERIA]))), keyby = .(Source_introgression, condition)]
Nsig[, Pct_sig := N_sig / N * 100, ]
Nsig[, Pct_sig_pop := N_sig / sum(N_sig) * 100, by = .(condition)]
Nsig[, Nsig_cond := sum(N_sig), by = .(condition)]
Nsig[, Ntest_cond := sum(N), by = .(condition)]
Nsig[, as.list(unlist(fisher.test(matrix(c(Ntest_cond - N, Nsig_cond - N_sig, N, N_sig), 2, byrow = TRUE))[c("estimate", "conf.int", "p.value")])), keyby = .(Source_introgression, condition)]

Nsig <- oligo_response_annot_perPOP[, .(N = length(unique(oligo)), N_sig = length(unique(oligo[oligo_diff_CRITERIA]))), keyby = .(allele_match, condition)]
Nsig[, Pct_sig := N_sig / N * 100, ]
Nsig[, Pct_sig_pop := N_sig / sum(N_sig) * 100, by = .(condition)]
Nsig[, Nsig_cond := sum(N_sig), by = .(condition)]
Nsig[, Ntest_cond := sum(N), by = .(condition)]
Nsig[, as.list(unlist(fisher.test(matrix(c(Ntest_cond - N, Nsig_cond - N_sig, N, N_sig), 2, byrow = TRUE))[c("estimate", "conf.int", "p.value")])), keyby = .(allele_match, condition)]

Nsig <- oligo_response_annot_perPOP[, .(N = length(unique(posID)), N_sig = length(unique(posID[oligo_diff_CRITERIA]))), keyby = .(POP, condition)]
Nsig[, Pct_sig := N_sig / N * 100, ]
Nsig[, Pct_sig_pop := N_sig / sum(N_sig) * 100, by = .(condition)]
Nsig[, Nsig_cond := sum(N_sig), by = .(condition)]
Nsig[, Ntest_cond := sum(N), by = .(condition)]
Nsig[, as.list(unlist(fisher.test(matrix(c(Ntest_cond - N, Nsig_cond - N_sig, N, N_sig), 2, byrow = TRUE))[c("estimate", "conf.int", "p.value")])), keyby = .(POP, condition)]

Nsig <- oligo_response_annot_perPOP[, .(N = length(unique(posID)), N_sig = length(unique(posID[oligo_diff_CRITERIA]))), keyby = .(Introgression_scenario_short, condition)]
Nsig[, Pct_sig := N_sig / N * 100, ]
Nsig[, Pct_sig_pop := N_sig / sum(N_sig) * 100, by = .(condition)]
Nsig[, Nsig_cond := sum(N_sig), by = .(condition)]
Nsig[, Ntest_cond := sum(N), by = .(condition)]
Nsig[, as.list(unlist(fisher.test(matrix(c(Ntest_cond - N, Nsig_cond - N_sig, N, N_sig), 2, byrow = TRUE))[c("estimate", "conf.int", "p.value")])), keyby = .(Introgression_scenario_short, condition)]


##### count frequency of responsive oligo
# Nsig <- oligo_response_annot_obs[, .(N=length(unique(oligo)), N_sig = length(unique(oligo[oligo_diff_CRITERIA]))), keyby = .(Introgressed_from, condition)]
Nsig <- oligo_response_annot_obs[, .(.N, N_sig = sum(oligo_diff_CRITERIA)), keyby = .(Introgressed_from, condition)]
Nsig[, Pct_sig := N_sig / N * 100, ]
Nsig[, Pct_sig_pop := N_sig / sum(N_sig) * 100, by = .(condition)]
Nsig[, Nsig_cond := sum(N_sig), by = .(condition)]
Nsig[, Ntest_cond := sum(N), by = .(condition)]
Nsig[, as.list(unlist(fisher.test(matrix(c(Ntest_cond - N, Nsig_cond - N_sig, N, N_sig), 2, byrow = TRUE))[c("estimate", "conf.int", "p.value")])), keyby = .(Introgressed_from, condition)]

Nsig[, Pct_sig_pop := N_sig / sum(N_sig) * 100, by = .(condition)]
dcast(Nsig, Introgressed_from ~ condition, value.var = "N_sig")
dcast(Nsig, Introgressed_from ~ condition, value.var = "Pct_sig")
Nsig <- oligo_response_annot_obs[, .(.N, N_sig = sum(oligo_diff_CRITERIA)), keyby = .(POP_adaptive, condition)]
Nsig[, Pct_sig := N_sig / N * 100, ]
Nsig[, Pct_sig_pop := N_sig / sum(N_sig) * 100, by = .(condition)]
Nsig[, Nsig_cond := sum(N_sig), by = .(condition)]
Nsig[, Ntest_cond := sum(N), by = .(condition)]
Nsig[, as.list(unlist(fisher.test(matrix(c(Ntest_cond - N, Nsig_cond - N_sig, N, N_sig), 2, byrow = TRUE))[c("estimate", "conf.int", "p.value")])), keyby = .(POP_adaptive, condition)]
dcast(Nsig, POP_adaptive ~ condition, value.var = "N_sig")
dcast(Nsig, POP_adaptive ~ condition, value.var = "Pct_sig")
dcast(Nsig, POP_adaptive ~ condition, value.var = "Pct_sig_pop")

# create a clean table of response to stim
SupTable2e <- oligo_response_annot_obs_archaic[, .(oligo, crsID, INTROGRESSED.allele, POP_adaptive, Introgressed_from, Introgression_scenario_v2, allele.label, celline, condition, log2FC_2vs1, log2FC.se, pval.LRT, pval_emp, lfdr, FDR, oligo_diff_CRITERIA)]
fwrite(SupTable2e, file = sprintf("%s/SupTables/SupTable2e__%s.tsv.gz", FIGURE_DIR, CRITERIA_DIFF), sep = "\t")
SupTable2e <- fread( sprintf("%s/SupTables/SupTable2e__%s.tsv.gz", FIGURE_DIR, CRITERIA_DIFF))
SupTable2e <- SupTable2e[crsID%in% SupTable2e[oligo_diff_CRITERIA==TRUE,crsID],]

SupTable1b_4161_final=fread(sprintf("%s/figures/%s/Z_figures/SupTable_1/SupTable1b_4161_testedSNPs.tsv", MPRA_DIR, ANALYSIS_DIR))
SupTable1b_4161_final[,prioritized_gene:=paste(prioritized_distance_10kb, prioritized_myeloid_enhancer, prioritized_lymphoid_enhancer, prioritized_lung_enhancer, prioritized_covid19,sep=':')]
SupTable1b_4161_final[,prioritized_gene:=paste(setdiff(unlist(str_split(prioritized_gene,':')),c('_','-')),collapse=':'),by=posID]
SupTable1b_4161_final[prioritized_gene=='',prioritized_gene:=nearestGene]

cols=c('log2FC_2vs1','log2FC.se','pval_emp','FDR')
my_levels<- names(color_setup_simplified_norep)[!grepl('NS',names(color_setup_simplified_norep))]
new_col_order <- CJ(my_levels, cols)
oo=order(factor(new_col_order$my_levels, levels = my_levels),factor(new_col_order$cols, levels = cols))
new_col_order <- new_col_order[oo,paste(cols, my_levels, sep = "_")]

SupTable2e[,score:=max(abs(log2FC_2vs1)-2*log2FC.se),by=crsID]


SupTable2e_wide=dcast(SupTable2e,score+crsID+oligo+allele.label~celline+factor(condition,condition_order),value.var=list('log2FC_2vs1','log2FC.se','pval_emp','FDR'))
setcolorder(SupTable2e_wide, c(setdiff(names(SupTable2e_wide), new_col_order), new_col_order))

SupTable2e_wide=merge(SupTable2e_wide,dcast(oligo_activity_obs_NS,crsID+oligo~celline,value.var=c('alpha_CRITERIA','oligo_class_CRITERIA')),by=c('crsID','oligo'))
SupTable2e_wide=merge(SupTable2e_wide,SupTable1b_4161_final[,.(crsID=posID, nearestGene, prioritized_gene)],by='crsID')
SupTable2e_wide=SupTable2e_wide[order(-score,allele.label),]
fwrite(SupTable2e_wide, file = sprintf("%s/SupTables/SupTable2e__%s.tsv", FIGURE_DIR, CRITERIA_DIFF), sep = "\t")

pi1_cond <- SupTable2e[, .(pi1 = 1 - mean(pval_emp > .05) / .95, pi1.se = sd(replicate(1000, 1 - mean(pval_emp[sample(seq_len(.N), replace = T)] > .05) / .95))), keyby = .(celline, condition)]
pi1_cond

fwrite(pi1_cond, file = sprintf("%s/pi1_cond__%s.tsv", FIGURE_DIR, CRITERIA_DIFF), sep = "\t")

pi1_byPop <- SupTable2e[, .(pi1 = 1 - mean(pval_emp > .05) / .95, pi1.se = sd(replicate(1000, 1 - mean(pval_emp[sample(seq_len(.N), replace = T)] > .05) / .95))), keyby = .(celline, condition, POP_adaptive)]
pi1_byPop[, Z := (pi1 - median(pi1)) / pi1.se, by = .(celline, condition)]
pi1_byPop[, pval := pnorm(abs(Z), lower.tail = F), by = .(celline, condition)]
pi1_byPop[, FDR := p.adjust(pval, method = "fdr")]

fwrite(pi1_byPop, file = sprintf("%s/pi1_byPop__%s.tsv", FIGURE_DIR, CRITERIA_DIFF), sep = "\t")
#################################################################
####### pairwise comparison between NS and stim #################
#################################################################

getOutlierP <- function(x, y) {
  XY <- cbind(x, y)
  D2 <- apply(scale(XY, center = TRUE, scale = FALSE), 1, function(xy) {
    matrix(xy, nrow = 1) %*% solve(var(XY, na.rm = TRUE)) %*% matrix(xy, ncol = 1)
  })
  pchisq(D2, df = 2, lower = FALSE)
}

#### load activity stats
oligo_activity_obs <- fread(sprintf("%s/all_oligos_annotated__%s.tsv.gz", ACTIVE_DIR, CRITERIA_ACTIVE))
oligo_activity_obs_ctc <- oligo_activity_obs[power == 0 & boot == 0 & ANALYSIS_SUBTYPE == "celltype_cond", ]


group_labels <- oligo_response_annot_obs_archaic[, .N, by = .(ANALYSIS_NAME, group1_labels, group2_labels, COND_ID, celline, condition)]
# i_response <- 5
# group1_label <- 'HepG2_NS_all'
# group2_label <- 'HepG2_IFNA2b_all'
compare_all <- list()
dir.create(sprintf("%s/scatter_response/", FIGURE_DIR), recursive = TRUE)
for (i_response in seq_len(nrow(group_labels))) {
  group1_label <- group_labels[i_response, group1_labels]
  group2_label <- group_labels[i_response, group2_labels]
  myCOND_ID <- group_labels[i_response, COND_ID]
  CELLINE <- group_labels[i_response, celline]
  CONDITION <- group_labels[i_response, condition]
  compare <- oligo_activity_obs_ctc[ANALYSIS_NAME %in% c(group1_label, group2_label), .(oligo, crsID, ANALYSIS_NAME, alpha_CRITERIA, log_alpha_se_CRITERIA, oligo_class_CRITERIA, celline, condition, COND_ID)]
  # option 1
  # compare <- dcast(compare, oligo + crsID + celline ~ relevel(as.factor(condition), "NS"), value.var = "alpha_CRITERIA")
  # compare <- merge(compare, oligo_response_annot_obs_archaic[ANALYSIS_NAME == group_labels[i_response, ANALYSIS_NAME], .(oligo, crsID, oligo_diff_CRITERIA, meanBC = nBC_LIB0 + nBC_LIB2)], by = c("oligo", "crsID"))
  # setnames(compare, group_labels[i_response, condition], "response")
  # option 2
  compare <- dcast(compare, oligo + crsID + celline ~ relevel(as.factor(condition), "NS"), value.var = list("alpha_CRITERIA", "log_alpha_se_CRITERIA"))
  compare <- merge(compare, oligo_response_annot_obs_archaic[ANALYSIS_NAME == group_labels[i_response, ANALYSIS_NAME], .(oligo, crsID, oligo_diff_CRITERIA, meanBC = nBC_LIB0 + nBC_LIB2)], by = c("oligo", "crsID"))
  setnames(compare, paste("alpha_CRITERIA", group_labels[i_response, condition], sep = "_"), "response")
  setnames(compare, paste("log_alpha_se_CRITERIA", group_labels[i_response, condition], sep = "_"), "log_alpha_se_CRITERIA_response")
  setnames(compare, "alpha_CRITERIA_NS", "NS")
  compare[, T_diff := (log(response) - log(NS)) / sqrt(log_alpha_se_CRITERIA_NS**2 + log_alpha_se_CRITERIA_response**2)]
  compare[, T_diff_norm := T_diff / IQR(T_diff) * 1.35]
  compare[, P_diff_norm := pnorm(abs(T_diff_norm), low = FALSE) * 2]
  compare[, FDR_diff_norm := p.adjust(P_diff_norm, method = "fdr")]
  compare_all[[group_labels[i_response, COND_ID]]] <- compare
}
compare_all <- rbindlist(compare_all, idcol = "COND_ID")

for (i_response in seq_len(nrow(group_labels))) {
  group1_label <- group_labels[i_response, group1_labels]
  group2_label <- group_labels[i_response, group2_labels]
  myCOND_ID <- group_labels[i_response, COND_ID]
  CELLINE <- group_labels[i_response, celline]
  CONDITION <- group_labels[i_response, condition]

  compare <- compare_all[COND_ID == myCOND_ID, ]
  compare <- merge(compare, NearestGene, by = "crsID", all.x = TRUE)
  compare[, P_outlier := getOutlierP(log2(NS), log2(response))]

  p <- ggplot(compare[order(oligo_diff_CRITERIA),], aes(x = log2(NS), y = log2(response), alpha = oligo_diff_CRITERIA, col = oligo_diff_CRITERIA, size = oligo_diff_CRITERIA)) 
  p <- p + geom_hline(yintercept = 0, col = "lightgrey") + geom_hline(yintercept = 0, col = "lightgrey")
  p <- p + geom_abline(aes(intercept = 0, slope = 1), col = "lightgrey", linetype = "dashed")
	p <- p + rasterize(geom_point(), dpi = 200)
  p <- p + xlab(paste("log oligo activity (", CELLINE, "- NS )"))
  p <- p + ylab(paste("log oligo activity (", CELLINE, "-", CONDITION, ")"))
  p <- p + scale_alpha_manual(values = c("FALSE" = 0.3, "TRUE" = 1))
  color_code <- c("FALSE" = "lightgrey", "TRUE" = unname(color_setup_simplified_norep[myCOND_ID]))
  p <- p + scale_color_manual(values = color_code)
  p <- p + scale_size_manual(values = c("FALSE" = 0.5, "TRUE" = 1))
  p <- p + theme_plot(fontsize = 12, lpos = "none")
  pdf(sprintf("%s/scatter_response/02a_activity_%s.pdf", FIGURE_DIR, myCOND_ID), height = 3, width = 3)
  print(p)
  dev.off()

  p0 <- p + geom_text_repel(data = compare[!is.na(TargetGene) & TargetGene != "" & oligo_diff_CRITERIA == TRUE & P_outlier < .01], aes(log2(NS), y = log2(response), label = TargetGene), col = color_code["TRUE"], size=2, show.legend = FALSE)
  p0 <- p0 + guides(text = FALSE)

  pdf(sprintf("%s/scatter_response/02b_activity_%s_labels.pdf", FIGURE_DIR, myCOND_ID), height = 3, width = 3)
  print(p0)
  dev.off()

  p0 <- p + geom_text_repel(data = compare[!is.na(TargetGene) & TargetGene != "" & oligo_diff_CRITERIA == TRUE & P_outlier < .01], aes(log2(NS), y = log2(response), label = TargetGene), col = color_code["TRUE"], size=1, show.legend = FALSE)
  p0 <- p0 + guides(text = FALSE)

  pdf(sprintf("%s/scatter_response/02b_activity_%s_labels_smaller.pdf", FIGURE_DIR, myCOND_ID), height = 3, width = 3)
  print(p0)
  dev.off()

  p <- ggplot(compare, aes(x = log2(NS), y = log2(response), alpha = oligo_diff_CRITERIA, size = oligo_diff_CRITERIA, col = log10(meanBC))) +
    rasterize(geom_point(), dpi = 200)
  p <- p + xlab(paste("log oligo activity (", CELLINE, "- NS )"))
  p <- p + ylab(paste("log oligo activity (", CELLINE, "-", CONDITION, ")"))
  p <- p + scale_alpha_manual(values = c("FALSE" = 0.1, "TRUE" = 1))
  p <- p + scale_size_manual(values = c("FALSE" = 0.5, "TRUE" = 1))
  p <- p + theme_plot(fontsize = 12, lpos = "right")
  p <- p + scale_color_viridis(option = "D", direction = 1)
  pdf(sprintf("%s/scatter_response/02c_activity_%s_BCnumber.pdf", FIGURE_DIR, myCOND_ID), height = 3, width = 4)
  print(p)
  dev.off()
}


# annotate oligos based on activity i the NS state
compare <- merge(oligo_activity_obs_ctc[condition == "NS", .(oligo, crsID, oligo_class_CRITERIA, cre_class_CRITERIA, ANALYSIS_NAME)], oligo_response_annot_obs_archaic[, .(oligo, crsID, INTROGRESSED.allele, POP_adaptive, Introgressed_from, Introgression_scenario_v2, allele.label, celline, condition, log2FC_2vs1, log2FC.se, pval.LRT, pval_emp, lfdr, FDR, oligo_diff_CRITERIA, ANALYSIS_NAME = group1_labels)], by = c("oligo", "crsID", "ANALYSIS_NAME"))


THRESHOLD <- 0
EffectStim <- compare[oligo_diff_CRITERIA == TRUE, .(
  N = length(unique(crsID)),
  N_inactive_up = length(unique(crsID[log2FC_2vs1 > THRESHOLD & grepl("inactive", oligo_class_CRITERIA)])),
  N_inactive_down = length(unique(crsID[log2FC_2vs1 < -THRESHOLD & grepl("inactive", oligo_class_CRITERIA)])),
  N_enhancer_up = length(unique(crsID[log2FC_2vs1 > THRESHOLD & grepl("enhancer", oligo_class_CRITERIA)])),
  N_enhancer_down = length(unique(crsID[log2FC_2vs1 < -THRESHOLD & grepl("enhancer", oligo_class_CRITERIA)])),
  N_silencer_up = length(unique(crsID[log2FC_2vs1 < -THRESHOLD & grepl("silencer", oligo_class_CRITERIA)])),
  N_silencer_down = length(unique(crsID[log2FC_2vs1 > THRESHOLD & grepl("silencer", oligo_class_CRITERIA)]))
), keyby = .(celline, condition)]
# EffectStim
#    celline condition     N N_inactive_up N_inactive_down N_enhancer_up N_enhancer_down N_silencer_up N_silencer_down
#     <char>    <char> <int>         <int>           <int>         <int>           <int>         <int>           <int>
# 1:    A549       IAV   180             3               0             0              65             0             112
# 2:    A549      SARS   421             7               3             2             167             0             246
# 3:    A549      TNFa    97            27               5            24              27             5              13
# 4:   HepG2       DEX     5             2               0             3               0             0               0
# 5:   HepG2    IFNA2b    97            14               8            34               8            31               6
# 6:   HepG2      TNFa    19             2               0             8               8             1               1
# 7:    K562       DEX    17             3               2             5               1             1               5
# 8:    K562    IFNA2b    23             5               2             9               1             5               2
# 9:    K562      TNFa    20            10               0             4               0             0               7

EffectStim[, .(Pct_enhancer_induced = (N_inactive_up + N_enhancer_up) / N, Pct_enhancer_diff = (N_enhancer_up + N_enhancer_down + N_inactive_up) / N, Pct_previously_inactive = (N_inactive_up + N_inactive_down) / N), keyby = .(celline, condition)]
# [1] 0.5257732 0.5263158 0.7000000

EffectStim[, (N_inactive_up + N_inactive_down) / N, keyby = .(celline, condition)]
# Key: <celline, condition>
#    celline condition         V1
#     <char>    <char>      <num>
# 1:    A549       IAV 0.04575163
# 2:    A549      SARS 0.05081301
# 3:    A549      TNFa 0.23809524
# 4:   HepG2       DEX 0.21428571
# 5:   HepG2    IFNA2b 0.28915663
# 6:   HepG2      TNFa 0.07692308
# 7:    K562       DEX 0.40000000
# 8:    K562    IFNA2b 0.40625000
# 9:    K562      TNFa 0.52941176

EffectStim[condition %in% c("IAV", "SARS"), (N_enhancer_down + N_silencer_down) / N]
#  0.9833333 0.9809976
fwrite(EffectStim, file = sprintf("%s/03a_direction_stim_effect__%s__source_data.tsv", FIGURE_DIR, CRITERIA_DIFF), sep = "\t")

compare[oligo_diff_CRITERIA == TRUE, .(
  N = length(unique(crsID)),
  N_inactive_up = length(unique(crsID[log2FC_2vs1 > THRESHOLD & grepl("inactive", oligo_class_CRITERIA)])),
  N_inactive_down = length(unique(crsID[log2FC_2vs1 < -THRESHOLD & grepl("inactive", oligo_class_CRITERIA)])),
  N_enhancer_up = length(unique(crsID[log2FC_2vs1 > THRESHOLD & grepl("enhancer", oligo_class_CRITERIA)])),
  N_enhancer_down = length(unique(crsID[log2FC_2vs1 < -THRESHOLD & grepl("enhancer", oligo_class_CRITERIA)])),
  N_silencer_up = length(unique(crsID[log2FC_2vs1 < -THRESHOLD & grepl("silencer", oligo_class_CRITERIA)])),
  N_silencer_down = length(unique(crsID[log2FC_2vs1 > THRESHOLD & grepl("silencer", oligo_class_CRITERIA)]))
)]
#      N N_inactive_up N_inactive_down N_enhancer_up N_enhancer_down N_silencer_up N_silencer_down
#    <int>         <int>           <int>         <int>           <int>         <int>           <int>
# 1:   605            63              20            66             209            41             284

compare[oligo_diff_CRITERIA == TRUE & !condition %in% c("IAV", "SARS"), .(N = length(unique(crsID))), by = .(cre_class_CRITERIA)][, Pct := N / sum(N)][1:.N, ]
#    cre_class_CRITERIA     N        Pct
#                <char> <int>      <num>
# 1:           silencer    62 0.25726141
# 2:           inactive    51 0.21161826
# 3:           enhancer    76 0.31535270
# 4:    strong enhancer    44 0.18257261
# 5:    strong silencer     8 0.03319502




names_stimulation <- names(color_setup_simplified_norep)[!grepl('-NS',names(color_setup_simplified_norep))]
names_stimulation <- gsub('_','-',names_stimulation)
names_stimulation <- names_stimulation[!grepl('-NS',names_stimulation)]

FigData <- melt(EffectStim, id.vars = c("celline", "condition"), measure.vars = c("N_enhancer_up", "N_enhancer_down", "N_silencer_up", "N_silencer_down", "N_inactive_up", "N_inactive_down"))
FigData[, element := gsub("N_(.*)_(up|down)", "\\1", variable)]
FigData[, direction := ifelse(grepl("up", variable), "Increased", "Decreased")]
p <- ggplot(FigData, aes(x = element, y = value, fill = paste(celline, condition, sep = "_"), alpha = direction))
p <- p + geom_bar(stat = "Identity", width = 0.8) + xlab("Activity in non-stimulated cells") + ylab("Number of differentially active CRE")
p <- p + scale_fill_manual(values = color_setup_simplified_norep)
p <- p + scale_alpha_manual(values = c("Increased" = 1, "Decreased" = 0.5))
p <- p + facet_wrap(~ factor(paste(celline, condition, sep = "-"),names_stimulation), scales = "free_y")
p <- p + theme_plot(fontsize = 10, lpos = "top", rotate.x = 90) + guides(fill = "none")
p <- p + guides(alpha = guide_legend(title = "Change in activity\nupon stimulation", title.position = "top", title.hjust = 0.5))
p <- p + theme(legend.title = element_text(size = 10, family = "sans")) # scale_y_sqrt(breaks=c(0,10,50,150,300))
pdf(sprintf("%s/scatter_response/03a_direction_stim_effect.pdf", FIGURE_DIR), height = 5, width = 2.5)
print(p)
dev.off()

FigData[,value:=as.double(value)]
FigData[value==0.01,value:=0.03]
FigData[value==0.03,value:=0.00]

p <- ggplot(FigData, aes(x = element, y = value, fill = paste(celline, condition, sep = "_"), alpha = direction))
p <- p + geom_bar(stat = "Identity", width = 0.8, position='dodge') + xlab("Activity in non-stimulated cells") + ylab("Number of differentially active CRE")
p <- p + scale_fill_manual(values = color_setup_simplified_norep)
p <- p + scale_alpha_manual(values = c("Increased" = 1, "Decreased" = 0.5))
p <- p + facet_wrap(~ factor(paste(celline, condition, sep = "-"),names_stimulation), scales = "free_y")
p <- p + theme_plot(fontsize = 10, lpos = "top", rotate.x = 90) + guides(fill = "none")
p <- p + guides(alpha = guide_legend(title = "Change in activity\nupon stimulation", title.position = "top", title.hjust = 0.5))
p <- p + theme(legend.title = element_text(size = 10, family = "sans")) # scale_y_sqrt(breaks=c(0,10,50,150,300))
p <- p + geom_hline(yintercept=0,col='lightgrey',linetype=3, linewidth=0.1)

pdf(sprintf("%s/scatter_response/03a_direction_stim_effect_v2.pdf", FIGURE_DIR), height = 5, width = 2.5)
print(p)
dev.off()

#################################################################
####### plot of CRS responsiveness VS TF binding ################
#################################################################

TFscore_annot <- fread(sprintf("%s/data/%s/00_TFscore/0_TFscore_all.tsv.gz", MPRA_DIR, ANALYSIS_DIR))
# TFBM <- fread(sprintf("%s/data/%s/00_TFBM/0_TFBM_all.tsv.gz", MPRA_DIR, ANALYSIS_DIR))


response_vs_score <- merge(compare, TFscore_annot, by = "oligo", allow.cartesian = TRUE)
# response_vs_score_test <- response_vs_score[, .(
#   rho = cor(log2FC_2vs1, score, method = "spearman"),
#   pval_spearman = cor.test(log2FC_2vs1, score, method = "spearman")$p.value,
#   pval_wilcox = wilcox.test(score[oligo_active_CRITERIA], score[!oligo_active_CRITERIA])$p.value
# ),
# by = .(TF = names, matrix_id, celline, condition)
# ]

# dir.create(sprintf("%s/",FIGURE_DIR))
# fwrite(response_vs_score_test, file = sprintf("%s/oligo_response_vs_TFscore.tsv.gz", FIGURE_DIR), sep = "\t")

wilcoxP <- function(x, y, ...) {
  as.numeric(try(wilcox.test(x, y, ...)$p.value))
}
medianDiff <- function(x, y, ...) {
  as.numeric(try(median(x, ...) - median(y, ...)))
}

response_vs_score_test_v2 <- response_vs_score[, .(
  rho = cor(log2FC_2vs1, score, method = "spearman"),
  r = cor(log2FC_2vs1, score),
  pval_spearman = cor.test(log2FC_2vs1, score, method = "spearman")$p.value,
  pval_pearson = cor.test(log2FC_2vs1, score, method = "pearson")$p.value,
  pval_wilcox = wilcoxP(score[oligo_diff_CRITERIA], score[!oligo_diff_CRITERIA]),
  pval_wilcox_up = wilcoxP(score[oligo_diff_CRITERIA & log2FC_2vs1 > 0], score[!oligo_diff_CRITERIA | log2FC_2vs1 < 0]),
  pval_wilcox_down = wilcoxP(score[oligo_diff_CRITERIA & log2FC_2vs1 < 0], score[!oligo_diff_CRITERIA | log2FC_2vs1 > 0]),
  delta_median = medianDiff(score[oligo_diff_CRITERIA], score[!oligo_diff_CRITERIA]),
  delta_median_up = medianDiff(score[oligo_diff_CRITERIA & log2FC_2vs1 > 0], score[!oligo_diff_CRITERIA | log2FC_2vs1 < 0]),
  delta_median_down = medianDiff(score[oligo_diff_CRITERIA & log2FC_2vs1 < 0], score[!oligo_diff_CRITERIA | log2FC_2vs1 > 0]),
  pval_wilcox_Pemp = wilcoxP(score[pval_emp <= 0.05], score[pval_emp > 0.05]),
  pval_wilcox_Pemp_up = wilcoxP(score[pval_emp <= 0.05 & log2FC_2vs1 > 0], score[pval_emp > 0.05 | log2FC_2vs1 < 0]),
  pval_wilcox_Pemp_down = wilcoxP(score[pval_emp <= 0.05 & log2FC_2vs1 < 0], score[pval_emp > 0.05 | log2FC_2vs1 > 0]),
  delta_median_Pemp = medianDiff(score[pval_emp <= 0.05], score[pval_emp > 0.05]),
  delta_median_Pemp_up = medianDiff(score[pval_emp <= 0.05 & log2FC_2vs1 > 0], score[pval_emp > 0.05 | log2FC_2vs1 < 0]),
  delta_median_Pemp_down = medianDiff(score[pval_emp <= 0.05 & log2FC_2vs1 < 0], score[pval_emp > 0.05 | log2FC_2vs1 > 0])
),
by = .(TF = names, matrix_id, celline, condition)
]

fwrite(response_vs_score_test_v2, file = sprintf("%s/SupTables/SupTable2f_oligo_response_vs_TFscore_quantitative_detail.tsv.gz", FIGURE_DIR), sep = "\t")
response_vs_score_test_v2= fread(sprintf("%s/SupTables/SupTable2f_oligo_response_vs_TFscore_quantitative_detail.tsv.gz", FIGURE_DIR))
# response_vs_score_test_v2[p.adjust(pval_wilcox,'fdr')<0.05,][order(celline,condition,pval_wilcox),.N,by=.(celline,condition)]
response_vs_score_test_short= response_vs_score_test_v2[order(factor(celline,levels=c('A549','HepG2','K562')),factor(condition,levels=c('IAV','SARS','IFNa','DEX','TNFa'))),.(celline,condition,matrix_id,TF,delta_median,pval_wilcox,FDR=p.adjust(pval_wilcox,'fdr'))]

fwrite(response_vs_score_test_short[FDR<0.05,][order(celline,condition,pval_wilcox)],file=sprintf("%s/SupTables/SupTable2f_oligo_response_vs_TFscore_quantitative_FDR5.tsv", FIGURE_DIR),sep='\t')

# only report where P wilcoxon passes FDR<5%%

response_vs_score_test_v2[matrix_id %in% c("MA0105.1", "MA0653.1", "MA0113.3"), ]

response_vs_score_test_v2[delta_median_up > 0, ][order(pval_wilcox_up), head(.SD, 5), keyby = .(celline, condition)]
response_vs_score_test_v2[delta_median_down > 0, ][order(pval_wilcox_down), head(.SD, 5), keyby = .(celline, condition)]

response_vs_score[, score_scaled := scale(score), by = matrix_id]

p <- ggplot(response_vs_score[matrix_id %in% c("MA0105.1", "MA0653.1"), ], aes(x = ifelse(oligo_diff_CRITERIA, "yes", "no"), y = score_scaled, fill = paste(celline, condition, sep = "_")))
p <- p + geom_violin(scale = "width", alpha = .5, linewidth = 0) 
p <- p + facet_grid(names ~ celline + factor(condition, levels = condition_order))
p <- p + theme_plot(rotate.x = 90, lpos = "none", fontsize = 12) 
p <- p + scale_fill_manual(values = color_setup_simplified_norep) + xlab("responsive CRS")
# p <- p + geom_boxplot(notch=TRUE,alpha=.5,fill='white')
p <- p + geom_boxplot(alpha = .5, fill = "black", width = .2, outlier.shape = NA) + stat_summary(fun.y = median, geom = "point", shape = "-", col = "white")
p <- p + theme(panel.border = element_blank()) + ylab("TF binding affiinity") + theme(strip.text.x = element_text(size = 8))
# p <- p + geom_signif(
# 		comparisons = list(c("yes", "no")),
# 		map_signif_level = c('*'=0.05,'**'=0.01,'***'=0.001),
# # 		textsize = 3 )
pdf(sprintf("%s/02c_boxplot_responseCRS_TFBS.pdf", FIGURE_DIR), height = 3, width = 4)
print(p)
dev.off()

cat("All done !")
q("no")

# oligo_response[ FDR < .05 & abs(log2FC_2vs1)> logFC_TH & excluded==FALSE & ctrl==FALSE, .(posID, celline, condition, type, POP_adaptive, Introgressed_from)]
# oligo_response_any <- oligo_response[ FDR < .05 & abs(log2FC_2vs1)> logFC_TH & excluded==FALSE & ctrl==FALSE, .(posID, celline, condition, type, POP_adaptive, Introgressed_from)]

# get motif for TRIM62 emVAR
# oligo_source[crsID=='1:33616310',]

# get motif for TNFSF13B emVAR
# oligo_source[crsID=='13:108995449',]


# To explore and to ADD?: there is a trend toward a weak enrichment of Denisova-introgressed variants in TNFa responsive CREs. it seems to be driven by Agta adaptive alleles, but …
# •	So far I see this based on allele matching between vindija and Denisova. Need to confirm with calls made on introgression source of the full haplotype rather than allele match.
# •	So far, I assign variants to a single population. I should instead use the same assignment as for Figure 1.
# •	This is done at the oligo x cell type level (not CRE), so there is a risk of overconfidence: if multiple oligos from the same CRE respond to TNFa similarly, or the same oligo responds similarly across cell types, we count it as independent observations, when they are not truly independant.
# It would seem here that the signal was driven by some Agta/Denisova variants sharing same response to TNF across cell lines. Maybe the best is to see if we can detect this signal at the level of response emVars so we can discuss what introgression actually did to immune response.
