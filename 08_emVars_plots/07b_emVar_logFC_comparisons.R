MPRA_DIR <- "/pasteur/helix/projects/evo_immuno_pop/MPRA"
ANALYSIS_DIR <- "MPRA_count_exp6_analysis3"
SCRIPT_DIR <- sprintf("%s/scripts/%s/", MPRA_DIR, ANALYSIS_DIR)

source(sprintf("%s/scripts/%s/00_starter.R", MPRA_DIR, ANALYSIS_DIR))

RUN_ID <- "RUN3_Z2_nBC10"
CRITERIA_ACTIVE_OUT <- "FDR5_scrambled5pct_FC0.2_GCnorm"
CRITERIA_DIFF_OUT <- "noDiffFilter_Diff_FDR5_FC.2"
CRITERIA_EMVARS_OUT <- "EmVar_FDR5_FC.2"
CRITERIA_EMVARS_DIFF <- "EmVarDiff_FDR5_FC.2"
STIM_ONLY <- FALSE

cmd <- commandArgs(trailingOnly = TRUE)
print(cmd)

for (i in seq_along(cmd)) {
  if (cmd[i] == "--run_dir" || cmd[i] == "-r") {
    RUN_ID <- cmd[i + 1]
  }
  if (cmd[i] == "--criteria_active" || cmd[i] == "-a") {
    CRITERIA_ACTIVE_OUT <- cmd[i + 1]
  }
  if (cmd[i] == "--criteria_diff" || cmd[i] == "-d") {
    CRITERIA_DIFF_OUT <- cmd[i + 1]
  }
  if (cmd[i] == "--criteria_emvar" || cmd[i] == "-e") {
    CRITERIA_EMVARS_OUT <- cmd[i + 1]
  }
  if (cmd[i] == "--criteria_emvar_diff" || cmd[i] == "-ed") {
    CRITERIA_EMVARS_DIFF <- cmd[i + 1]
  }
  if (cmd[i] == "--stim_only" || cmd[i] == "-s") {
    STIM_ONLY <- as.logical(cmd[i + 1])
  }
}

CRITERIA_ACTIVE <- gsub("noActivityFilter_", "", CRITERIA_ACTIVE_OUT)
CRITERIA_DIFF <- gsub("noDiffFilter_", "", CRITERIA_DIFF_OUT)
CRITERIA_EMVARS <- gsub("noEmVarFilter_", "", CRITERIA_EMVARS_OUT)
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
SNP_annot_v5 <- merge(SNP_annot_v4[, -"Introgression_scenario"], selected_annot_wide, by = c("ID", "posID"))
toc()


##### active parameters
source(sprintf("%s/scripts/%s/03_00_parameter_activity.R", MPRA_DIR, ANALYSIS_DIR))
##### diff parameters
source(sprintf("%s/scripts/%s/04_00_parameter_diff_activity.R", MPRA_DIR, ANALYSIS_DIR))
##### emVars parameters
source(sprintf("%s/scripts/%s/05_00_parameter_emVars.R", MPRA_DIR, ANALYSIS_DIR))
##### diff parameters
source(sprintf("%s/scripts/%s/06_00_parameter_emVars_Diff.R", MPRA_DIR, ANALYSIS_DIR))

IN_DIR <- sprintf("%s/data/%s/03a_aggMPRA_analyse/%s/", MPRA_DIR, ANALYSIS_DIR, RUN_ID)
ACTIVE_DIR <- sprintf("%s/data/%s/03b_computeFDR/%s/Activity/%s/", MPRA_DIR, ANALYSIS_DIR, RUN_ID, CRITERIA_ACTIVE)
DIFF_DIR <- sprintf("%s/data/%s/03b_computeFDR/%s/Diff/%s/%s/", MPRA_DIR, ANALYSIS_DIR, RUN_ID, CRITERIA_DIFF, CRITERIA_ACTIVE_OUT)
EMVAR_DIR <- sprintf("%s/data/%s/03b_computeFDR/%s/emVars/%s/%s", MPRA_DIR, ANALYSIS_DIR, RUN_ID, CRITERIA_EMVARS, CRITERIA_ACTIVE_OUT)
EMVARDIFF_DIR <- sprintf("%s/data/%s/03b_computeFDR/%s/EmVar_Diff/%s/%s/%s/%s/%s", MPRA_DIR, ANALYSIS_DIR, RUN_ID, CRITERIA_EMVARS_DIFF, CRITERIA_EMVARS_OUT, CRITERIA_DIFF_OUT, CRITERIA_ACTIVE_OUT, STIM_ONLY_OUT)


IN_DIR <- sprintf("%s/data/%s/03a_aggMPRA_analyse/%s/", MPRA_DIR, ANALYSIS_DIR, RUN_ID)
FIGURE_DIR <- sprintf("%s/figures/%s/%s/04_emVars/%s/%s", MPRA_DIR, ANALYSIS_DIR, RUN_ID, CRITERIA_EMVARS, CRITERIA_ACTIVE_OUT)

dir.create(FIGURE_DIR, recursive = TRUE)
dir.create(sprintf("%s/comparison_effect_size/", FIGURE_DIR), recursive = TRUE)

dir.create(sprintf("%s/SupTables/", FIGURE_DIR), recursive = TRUE)


emVARs_all_results <- fread(sprintf("%s/all_emVars_results.tsv.gz", IN_DIR), sep = "\t")

oligo_activity_file <- sprintf("%s/oligo_activity__all__%s.tsv.gz", ACTIVE_DIR, CRITERIA_ACTIVE)
oligo_activity <- fread(oligo_activity_file)

oligo_target <- fread(sprintf("%s/data/%s/00_oligo_targets_v4.txt", MPRA_DIR, ANALYSIS_DIR), sep = "\t")
oligo_target_collapsed <- fread(sprintf("%s/data/%s/00_oligo_targets_collapsed_v4.txt", MPRA_DIR, ANALYSIS_DIR), sep = "\t")
crs_targets <- unique(oligo_target_collapsed[, -"oligo"])
NearestGene <- unique(oligo_target[method == "NearestGene" & gene_type != "", .(posID = crsID, TargetGene)])

emVARs_annot_obs <- fread(file = sprintf("%s/all_emVars_annotated_celltype__%s.tsv", EMVAR_DIR, CRITERIA_ACTIVE_OUT))

SNP_annot_emVars <- fread(file = sprintf("%s/SNP_annot_emVars__%s.tsv", EMVAR_DIR, CRITERIA_ACTIVE_OUT))
any_emVars <- SNP_annot_emVars[nCelltype>0,posID]
celline_specific_emVars <- SNP_annot_emVars[nCelltype==1,posID]


logFC_emVars <- dcast(emVARs_annot_obs[is_introgressed & posID %chin% posID[is_emVar_CRITERIA == TRUE], ], posID ~ ANALYSIS_NAME, value.var = "log2FC_archaic_vs_modern", fun.aggregate = mean)
FDR_emVars <- dcast(emVARs_annot_obs[is_introgressed & posID %chin% posID[is_emVar_CRITERIA == TRUE], ], posID ~ ANALYSIS_NAME, value.var = "FDR", fun.aggregate = mean)
Pval_emVars <- dcast(emVARs_annot_obs[is_introgressed & posID %chin% posID[is_emVar_CRITERIA == TRUE], ], posID ~ ANALYSIS_NAME, value.var = "pval_emp", fun.aggregate = mean)

COR_matrix_shared <- cor(merge(SNP_annot_emVars[nContext == 3, .(posID, nContext)], logFC_emVars, by = "posID")[, -(1:2)], use = "p")
COR_matrix_all <- cor(logFC_emVars[, -1], use = "p")
COR_matrix <- COR_matrix_all[condition_summary$analysis_name, condition_summary$analysis_name]
COR_matrix[lower.tri(COR_matrix)] <- COR_matrix_shared[condition_summary$analysis_name, condition_summary$analysis_name][lower.tri(COR_matrix)]
rownames(COR_matrix) <- condition_summary$COND_ID
colnames(COR_matrix) <- condition_summary$COND_ID

pdf(sprintf("%s/comparison_effect_size/01_corr_logFC_matrix_split_emVars_vs_shared_emVars.pdf", FIGURE_DIR))
# corrplot(COR_matrix,tl.col=color_setup[condition_summary_reps$analysis_name],labels=names(color_setup_simplified)[match(condition_summary_reps$analysis_name,names(color_setup))])
corrplot(COR_matrix, tl.col = color_setup_simplified_norep[condition_summary$COND_ID])
dev.off()

COR_matrix_shared <- cor(merge(SNP_annot_emVars[nContext == 3, .(posID, nContext)], logFC_emVars, by = "posID")[, -(1:2)], use = "p", method = "s")
COR_matrix_all <- cor(logFC_emVars[, -1], use = "p", method = "s")
COR_matrix <- COR_matrix_all[condition_summary$analysis_name, condition_summary$analysis_name]
COR_matrix[lower.tri(COR_matrix)] <- COR_matrix_shared[condition_summary$analysis_name, condition_summary$analysis_name][lower.tri(COR_matrix)]
rownames(COR_matrix) <- condition_summary$COND_ID
colnames(COR_matrix) <- condition_summary$COND_ID

range(COR_matrix_all[upper.tri(COR_matrix)])
# 0.4034166 0.7487846
range(COR_matrix_shared[upper.tri(COR_matrix)])
# 0.3973363 0.8562329
all_pairs <- combn(logFC_emVars[, -1], 2, simplify = FALSE)
sapply(all_pairs, function(x) cor.test(x[[1]], x[[2]], method = "s")$p.value)


pdf(sprintf("%s/comparison_effect_size/01_corr_logFC_matrix_split_emVars_vs_shared_emVars_spearman.pdf", FIGURE_DIR))
# corrplot(COR_matrix,tl.col=color_setup[condition_summary_reps$analysis_name],labels=names(color_setup_simplified)[match(condition_summary_reps$analysis_name,names(color_setup))])
corrplot(COR_matrix, tl.col = color_setup_simplified_norep[condition_summary$COND_ID])
dev.off()


# emVARs_annot_obs <- merge(emVARs_annot_obs, condition_summary[, .(COND_ID, ANALYSIS_NAME = analysis_name)], by = "ANALYSIS_NAME")


########################################################################################################################
############################## COMPARISON BETWEEN EMVARs : HEPG2 vs A549 ######################################################
########################################################################################################################

logFC_all <- dcast(emVARs_annot_obs[is_introgressed & !ctrl & !excluded & ANALYSIS_SUBTYPE == "celltype_cond" & COND_ID %in% c("HepG2_NS", "A549_NS"), ], posID + oligo1 + oligo2 ~ COND_ID, value.var = list("log2FC_archaic_vs_modern", "is_emVar_CRITERIA", "FDR", 'pval_emp')) # , fun.aggregate = list(mean, any, mean)
logFC_all <- merge(logFC_all, NearestGene, by = "posID", all.x=TRUE)

emVARs_Diff_obs_celltype <- fread(sprintf("%s/all_emVars_diff_obs_celltype.tsv", EMVARDIFF_DIR))
n_emVars=emVARs_Diff_obs_celltype[analysis_subtype=='celltype_comp' & posID %in% any_emVars,length(unique(posID))]

NS_detected_emVars <- emVARs_annot_obs[posID %in% any_emVars &pval_emp < 0.05 & condition=='NS',unique(posID)]
STIM_detected_emVars <- emVARs_annot_obs[posID %in% any_emVars & pval_emp < 0.05 & condition!='NS',unique(posID)]
length(union(NS_detected_emVars, STIM_detected_emVars))
length(setdiff(any_emVars, STIM_detected_emVars))
length(setdiff(any_emVars, NS_detected_emVars))

Table_S4a=emVARs_annot_obs[posID %in% any_emVars & pval_emp < 0.05,][order(factor(COND_ID, levels = names(color_setup_simplified_norep))),][,.( 
    condition_A549=paste(condition[celline=='A549'], collapse = ","),
    condition_HepG2=paste(condition[celline=='HepG2'], collapse = ","),
    condition_K562=paste(condition[celline=='K562'], collapse = ","),
    is_detected_NS_any=ifelse(any(condition=='NS'),'yes',''),
    is_detected_STIM_any=ifelse(any(condition!='NS'),'yes',''),
    is_detected_NS_A549=ifelse(any(condition=='NS' & celline=='A549'),'yes',''),
    is_detected_STIM_A549=ifelse(any(condition!='NS' & celline=='A549'),'yes',''),
    is_detected_NS_HepG2=ifelse(any(condition=='NS' & celline=='HepG2'),'yes',''),
    is_detected_STIM_HepG2=ifelse(any(condition!='NS' & celline=='HepG2'),'yes',''),
    is_detected_NS_K562=ifelse(any(condition=='NS' & celline=='K562'),'yes',''),
    is_detected_STIM_K562=ifelse(any(condition!='NS' & celline=='K562'),'yes','')),by=.(posID,rsID)]
Table_S4a[, emVar_type_A549 := case_when(is_detected_NS_A549=='yes' & is_detected_STIM_A549=='yes' ~ 'shared',
                                    is_detected_NS_A549=='yes' & is_detected_STIM_A549=='' ~ 'NS_specific',
                                    is_detected_NS_A549=='' & is_detected_STIM_A549=='yes' ~ 'STIM_specific')]
Table_S4a[, emVar_type_HepG2 := case_when(is_detected_NS_HepG2=='yes' & is_detected_STIM_HepG2=='yes' ~ 'shared',
                                    is_detected_NS_HepG2=='yes' & is_detected_STIM_HepG2=='' ~ 'NS_specific',
                                    is_detected_NS_HepG2=='' & is_detected_STIM_HepG2=='yes' ~ 'STIM_specific')]
Table_S4a[, emVar_type_K562 := case_when(is_detected_NS_K562=='yes' & is_detected_STIM_K562=='yes' ~ 'shared',
                                    is_detected_NS_K562=='yes' & is_detected_STIM_K562=='' ~ 'NS_specific',
                                    is_detected_NS_K562=='' & is_detected_STIM_K562=='yes' ~ 'STIM_specific')]
Table_S4a[, emVar_type_any := case_when(is_detected_NS_any=='yes' & is_detected_STIM_any=='yes' ~ 'shared',
                                    is_detected_NS_any=='yes' & is_detected_STIM_any=='' ~ 'NS_specific',
                                    is_detected_NS_any=='' & is_detected_STIM_any=='yes' ~ 'STIM_specific')]

SupTable1b_4161_final=fread(sprintf("%s/figures/%s/Z_figures/SupTable_1/SupTable1b_4161_testedSNPs.tsv", MPRA_DIR, ANALYSIS_DIR))
SupTable1b_4161_final[,prioritized_gene:=paste(prioritized_distance_10kb, prioritized_myeloid_enhancer, prioritized_lymphoid_enhancer, prioritized_lung_enhancer, prioritized_covid19,sep=':')]
SupTable1b_4161_final[,prioritized_gene:=paste(setdiff(unlist(str_split(prioritized_gene,':')),c('_','-')),collapse=':'),by=posID]
SupTable1b_4161_final[prioritized_gene=='',prioritized_gene:=nearestGene]

Table_S4a= merge(Table_S4a,SupTable1b_4161_final[,.(posID, rsID, ANCESTRAL, DERIVED, INTROGRESSED, CHROM,POS_b37,POP_adaptive,maxFreq, Source_introgression,nearestGene,prioritized_gene)], by=c('posID','rsID'))
Table_S4a=Table_S4a[order(CHROM,POS_b37),.(rsID, posID, INTROGRESSED_OTHER=paste(INTROGRESSED,ifelse(INTROGRESSED==DERIVED,ANCESTRAL,DERIVED),sep='/'),white1='',
																condition_A549, is_detected_NS_A549, is_detected_STIM_A549, emVar_type_A549, white2='',
																condition_HepG2, is_detected_NS_HepG2, is_detected_STIM_HepG2, emVar_type_HepG2, white3='',
																condition_K562, is_detected_NS_K562, is_detected_STIM_K562, emVar_type_K562, white4='',
																is_detected_NS_any, is_detected_STIM_any, emVar_type_any)]
write.table(Table_S4a, file = sprintf("%s/SupTables/SupTable4a_condition_specific_emVars.tsv", FIGURE_DIR), sep = "\t", quote = FALSE, row.names = FALSE)




n_diff_CT=emVARs_Diff_obs_celltype[analysis_subtype=='celltype_comp' & posID %in% any_emVars & is_emVar_Diff_CRITERIA==TRUE,length(unique(posID))]
n_diff_CT/n_emVars
# 0.4092888

n_diff_CT_suggestive=emVARs_Diff_obs_celltype[analysis_subtype=='celltype_comp' & posID %in% any_emVars & pval_emp<0.05,length(unique(posID))]
n_diff_CT_suggestive/n_emVars
#0.7024673

n_Cellline_specific=emVARs_Diff_obs_celltype[analysis_subtype=='celltype_comp' & posID %in% celline_specific_emVars,length(unique(posID))]
n_Cellline_specific
#537
n_Cellline_specific_diff_CT=emVARs_Diff_obs_celltype[analysis_subtype=='celltype_comp' & posID %in% celline_specific_emVars & is_emVar_Diff_CRITERIA==TRUE,length(unique(posID))]
n_Cellline_specific_diff_CT
# 196
n_Cellline_specific_diff_CT/n_Cellline_specific
# 0.3649907

emVARs_Diff_obs_HepG2_vs_A549 <- emVARs_Diff_obs_celltype[group2_labels == "A549-ACE2_NS_all" & group1_labels == "HepG2_NS_all", ]
logFC_all <- merge(logFC_all, emVARs_Diff_obs_HepG2_vs_A549[, .(posID, is_emVar_Diff_CRITERIA, Delta_log2FC_archaic_vs_modern, pval_interaction=pval_emp)], by = "posID")
# original
logFC_all[, sharing_HepG2_A549 := case_when(
  (is_emVar_CRITERIA_HepG2_NS == TRUE | is_emVar_CRITERIA_A549_NS == TRUE) & pval_emp_HepG2_NS < 0.05 & pval_emp_A549_NS < 0.05 & sign(log2FC_archaic_vs_modern_HepG2_NS) == sign(log2FC_archaic_vs_modern_A549_NS) & pval_interaction>.05 ~ "shared, identical effect",
  (is_emVar_CRITERIA_HepG2_NS == TRUE | is_emVar_CRITERIA_A549_NS == TRUE) & pval_emp_HepG2_NS < 0.05 & pval_emp_A549_NS < 0.05 & sign(log2FC_archaic_vs_modern_HepG2_NS) == sign(log2FC_archaic_vs_modern_A549_NS) & pval_interaction<.05 & sign(Delta_log2FC_archaic_vs_modern) != sign(log2FC_archaic_vs_modern_HepG2_NS) ~ "shared, stronger in HepG2",
  (is_emVar_CRITERIA_HepG2_NS == TRUE | is_emVar_CRITERIA_A549_NS == TRUE) & pval_emp_HepG2_NS < 0.05 & pval_emp_A549_NS < 0.05 & sign(log2FC_archaic_vs_modern_HepG2_NS) == sign(log2FC_archaic_vs_modern_A549_NS) & pval_interaction<.05 & sign(Delta_log2FC_archaic_vs_modern) == sign(log2FC_archaic_vs_modern_HepG2_NS) ~ "shared, stronger in A549",
  (is_emVar_CRITERIA_HepG2_NS == TRUE & is_emVar_CRITERIA_A549_NS == TRUE) & pval_emp_HepG2_NS < 0.05 & pval_emp_A549_NS < 0.05 & sign(log2FC_archaic_vs_modern_HepG2_NS) != sign(log2FC_archaic_vs_modern_A549_NS) & pval_interaction<.05 ~ "shared, opposite effect",
  (is_emVar_CRITERIA_HepG2_NS == TRUE & pval_emp_A549_NS > 0.05 & pval_interaction<.05) ~ "HepG2 specific",
  (is_emVar_CRITERIA_HepG2_NS == TRUE & pval_emp_A549_NS > 0.05 & pval_interaction>.05) ~ "shared, identical effect",
  (is_emVar_CRITERIA_A549_NS == TRUE & pval_emp_HepG2_NS > 0.05  & pval_interaction<.05) ~ "A549 specific",
  (is_emVar_CRITERIA_A549_NS == TRUE & pval_emp_HepG2_NS > 0.05 & pval_interaction>.05) ~ "shared, identical effect",
  is_emVar_CRITERIA_HepG2_NS == FALSE & is_emVar_CRITERIA_A549_NS == FALSE ~ "non significant"
)]

# # modified
# logFC_all[, sharing_HepG2_A549 := case_when(
#   (is_emVar_CRITERIA_HepG2_NS == TRUE | is_emVar_CRITERIA_A549_NS == TRUE) & pval_emp_HepG2_NS < 0.05 & pval_emp_A549_NS < 0.05 & sign(log2FC_archaic_vs_modern_HepG2_NS) == sign(log2FC_archaic_vs_modern_A549_NS) & pval_interaction>.05 ~ "shared, identical effect",
#   (is_emVar_CRITERIA_HepG2_NS == TRUE | is_emVar_CRITERIA_A549_NS == TRUE) & pval_emp_HepG2_NS < 0.05 & pval_emp_A549_NS < 0.05 & sign(log2FC_archaic_vs_modern_HepG2_NS) == sign(log2FC_archaic_vs_modern_A549_NS) & pval_interaction<.05 & sign(Delta_log2FC_archaic_vs_modern) != sign(log2FC_archaic_vs_modern_HepG2_NS) ~ "shared, stronger in HepG2",
#   (is_emVar_CRITERIA_HepG2_NS == TRUE | is_emVar_CRITERIA_A549_NS == TRUE) & pval_emp_HepG2_NS < 0.05 & pval_emp_A549_NS < 0.05 & sign(log2FC_archaic_vs_modern_HepG2_NS) == sign(log2FC_archaic_vs_modern_A549_NS) & pval_interaction<.05 & sign(Delta_log2FC_archaic_vs_modern) == sign(log2FC_archaic_vs_modern_HepG2_NS) ~ "shared, stronger in A549",
#   (is_emVar_CRITERIA_HepG2_NS == TRUE & is_emVar_CRITERIA_A549_NS == TRUE) & pval_emp_HepG2_NS < 0.05 & pval_emp_A549_NS < 0.05 & sign(log2FC_archaic_vs_modern_HepG2_NS) != sign(log2FC_archaic_vs_modern_A549_NS) & pval_interaction<.05 ~ "shared, opposite effect",
#   (is_emVar_CRITERIA_HepG2_NS == TRUE & pval_emp_A549_NS > 0.05) ~ "HepG2 specific",
#   (is_emVar_CRITERIA_A549_NS == TRUE & pval_emp_HepG2_NS > 0.05 ) ~ "A549 specific",
#   is_emVar_CRITERIA_HepG2_NS == FALSE & is_emVar_CRITERIA_A549_NS == FALSE ~ "non significant"
# )]



color_code=c( "shared, identical effect" = "#d6c60080",
  "shared, undetected in A549" = mergeCols("#d6c600",unname(color_celline_2levels_light["HepG2 weak"])),
  "shared, undetected in HepG2" = mergeCols("#d6c600",unname(color_celline_2levels_light["A549 weak"])),
  "shared, opposite effect" = "#DD4466",
  "shared, stronger in HepG2" = unname(color_celline_2levels_light["HepG2 weak"]),
  "shared, stronger in A549" = unname(color_celline_2levels_light["A549 weak"]),
  "HepG2 specific" = unname(color_celline_2levels_light["HepG2 strong"]),
  "A549 specific" = unname(color_celline_2levels_light["A549 strong"])
)



color_code_dark=c( "shared, identical effect" = mergeCols("#d6c600",'black',0.25),
  "shared, undetected in A549" = mergeCols("#d6c600",unname(color_celline_2levels["HepG2 weak"])),
  "shared, undetected in HepG2" = mergeCols("#d6c600",unname(color_celline_2levels["A549 weak"])),
  "shared, opposite effect" = "#DD4466",
  "shared, stronger in HepG2" = unname(color_celline_2levels["HepG2 weak"]),
  "shared, stronger in A549" = unname(color_celline_2levels["A549 weak"]),
  "HepG2 specific" = mergeCols(color_celline_2levels["HepG2 weak"],'black',0.25),
  "A549 specific" = mergeCols(color_celline_2levels["A549 weak"],'black',0.25)
)

pdf(sprintf("%s/comparison_effect_size/02_logFC_HepG2_vs_A549_NS.pdf", FIGURE_DIR), height = 2.2, width = 4)
p <- ggplot(logFC_all[!is.na(sharing_HepG2_A549) & sharing_HepG2_A549 != "non significant"], aes(x = log2FC_archaic_vs_modern_HepG2_NS, y = log2FC_archaic_vs_modern_A549_NS))
p <- p + rasterize(geom_point(aes(col = factor(sharing_HepG2_A549,names(color_code)))), dpi = 200)
p <- p + xlab("log2 FC (HepG2)") + ylab("log2 FC (A549)")
p <- p + geom_hline(yintercept = 0, col = "grey") + geom_vline(xintercept = 0, col = "grey")
p <- p + theme_plot(rotate.x = 90, lpos = "right", fontsize = 11)
p <- p + geom_abline(intercept = 0, slope = 1, col = "lightgrey", linetype = 2)
p <- p + scale_color_manual(values = color_code)
p <- p + guides(col = guide_legend(ncol=1, byrow = FALSE)) + theme(legend.box.spacing = unit(.1, "mm"), legend.spacing.y = unit(.1, 'cm'))
print(p)
dev.off()


	getOutlierP <- function(x, y) {
		XY <- cbind(x, y)
		D2 <- apply(scale(XY, center = TRUE, scale = FALSE), 1, function(xy) {
			matrix(xy, nrow = 1) %*% solve(var(XY, na.rm = TRUE)) %*% matrix(xy, ncol = 1)
		})
		pchisq(D2, df = 2, lower = FALSE)
	}
	library(ggnewscale)
	logFC_all[, P_outlier := getOutlierP(log2FC_archaic_vs_modern_HepG2_NS, log2FC_archaic_vs_modern_A549_NS)]
	p <- p + new_scale_colour() + scale_color_manual(values = color_code_dark)
	p <- p + geom_text_repel(data = logFC_all[!is.na(sharing_HepG2_A549) & sharing_HepG2_A549 != "non significant" & P_outlier < .01 | sharing_HepG2_A549 == "shared, opposite effect"], aes(x = log2FC_archaic_vs_modern_HepG2_NS, y = log2FC_archaic_vs_modern_A549_NS,  col = sharing_HepG2_A549, label = TargetGene), size=2, show.legend = FALSE) #,
	p <- p + guides(text = FALSE)
	pdf(sprintf("%s/comparison_effect_size/02b_logFC_HepG2_vs_A549_withnames.pdf", FIGURE_DIR),  height = 3.2, width = 5)
	print(p)
	dev.off()


########################################################################################################################
############################## COMPARISON BETWEEN EMVARs : HEPG2 vs K562 ###############################################
########################################################################################################################


logFC_all <- dcast(emVARs_annot_obs[is_introgressed & !ctrl & !excluded & ANALYSIS_SUBTYPE == "celltype_cond" & COND_ID %in% c("HepG2_NS", "K562_NS"), ], posID + oligo1 + oligo2 ~ COND_ID, value.var = list("log2FC_archaic_vs_modern", "is_emVar_CRITERIA", "FDR", 'pval_emp')) # , fun.aggregate = list(mean, any, mean)
logFC_all <- merge(logFC_all, NearestGene, by = "posID", all.x=TRUE)

emVARs_Diff_obs_celltype <- fread(sprintf("%s/all_emVars_diff_obs_celltype.tsv", EMVARDIFF_DIR))
emVARs_Diff_obs_HepG2_vs_K562 <- emVARs_Diff_obs_celltype[group2_labels == "K562_NS_all" & group1_labels == "HepG2_NS_all", ]
logFC_all <- merge(logFC_all, emVARs_Diff_obs_HepG2_vs_K562[, .(posID, is_emVar_Diff_CRITERIA, Delta_log2FC_archaic_vs_modern, pval_interaction=pval_emp)], by = "posID")
logFC_all[, sharing_HepG2_K562 := case_when(
  (is_emVar_CRITERIA_HepG2_NS == TRUE | is_emVar_CRITERIA_K562_NS == TRUE) & pval_emp_HepG2_NS < 0.05 & pval_emp_K562_NS < 0.05 & sign(log2FC_archaic_vs_modern_HepG2_NS) == sign(log2FC_archaic_vs_modern_K562_NS) & pval_interaction>.05 ~ "shared, identical effect",
  (is_emVar_CRITERIA_HepG2_NS == TRUE | is_emVar_CRITERIA_K562_NS == TRUE) & pval_emp_HepG2_NS < 0.05 & pval_emp_K562_NS < 0.05 & sign(log2FC_archaic_vs_modern_HepG2_NS) == sign(log2FC_archaic_vs_modern_K562_NS) & pval_interaction<.05 & sign(Delta_log2FC_archaic_vs_modern) != sign(log2FC_archaic_vs_modern_HepG2_NS) ~ "shared, stronger in HepG2",
  (is_emVar_CRITERIA_HepG2_NS == TRUE | is_emVar_CRITERIA_K562_NS == TRUE) & pval_emp_HepG2_NS < 0.05 & pval_emp_K562_NS < 0.05 & sign(log2FC_archaic_vs_modern_HepG2_NS) == sign(log2FC_archaic_vs_modern_K562_NS) & pval_interaction<.05 & sign(Delta_log2FC_archaic_vs_modern) == sign(log2FC_archaic_vs_modern_HepG2_NS) ~ "shared, stronger in K562",
  (is_emVar_CRITERIA_HepG2_NS == TRUE & is_emVar_CRITERIA_K562_NS == TRUE) & pval_emp_HepG2_NS < 0.05 & pval_emp_K562_NS < 0.05 & sign(log2FC_archaic_vs_modern_HepG2_NS) != sign(log2FC_archaic_vs_modern_K562_NS) & pval_interaction<.05 ~ "shared, opposite effect",
  (is_emVar_CRITERIA_HepG2_NS == TRUE & pval_emp_K562_NS > 0.05 & pval_interaction<.05) ~ "HepG2 specific",
  (is_emVar_CRITERIA_HepG2_NS == TRUE & pval_emp_K562_NS > 0.05 & pval_interaction>.05) ~ "shared, identical effect",
  (is_emVar_CRITERIA_K562_NS == TRUE & pval_emp_HepG2_NS > 0.05  & pval_interaction<.05) ~ "K562 specific",
  (is_emVar_CRITERIA_K562_NS == TRUE & pval_emp_HepG2_NS > 0.05 & pval_interaction>.05) ~ "shared, identical effect",
  is_emVar_CRITERIA_HepG2_NS == FALSE & is_emVar_CRITERIA_K562_NS == FALSE ~ "non significant"
)]


color_code=c( "shared, identical effect" = "#d6c600",
  "shared, undetected in K562" = mergeCols("#d6c600",unname(color_celline_2levels_light["HepG2 weak"])),
  "shared, undetected in HepG2" = mergeCols("#d6c600",unname(color_celline_2levels_light["K562 weak"])),
  "shared, opposite effect" = "#DD4466",
  "shared, stronger in HepG2" = unname(color_celline_2levels_light["HepG2 weak"]),
  "shared, stronger in K562" = unname(color_celline_2levels_light["K562 weak"]),
  "HepG2 specific" = unname(color_celline_2levels_light["HepG2 strong"]),
  "K562 specific" = unname(color_celline_2levels_light["K562 strong"])
)

pdf(sprintf("%s/comparison_effect_size/02_logFC_HepG2_vs_K562_NS.pdf", FIGURE_DIR), height = 2.2, width = 4)
p <- ggplot(logFC_all[!is.na(sharing_HepG2_K562) & sharing_HepG2_K562 != "non significant"], aes(x = log2FC_archaic_vs_modern_HepG2_NS, y = log2FC_archaic_vs_modern_K562_NS))
p <- p + rasterize(geom_point(aes(col = factor(sharing_HepG2_K562,names(color_code))),alpha = .5), dpi = 200)
p <- p + xlab("log2 FC (HepG2)") + ylab("log2 FC (K562)")
p <- p + geom_hline(yintercept = 0, col = "grey") + geom_vline(xintercept = 0, col = "grey")
p <- p + theme_plot(rotate.x = 90, lpos = "right", fontsize = 11)
p <- p + geom_abline(intercept = 0, slope = 1, col = "lightgrey", linetype = 2)
p <- p + scale_color_manual(values = color_code)
p <- p +guides(col = guide_legend(ncol=1, byrow = FALSE)) + theme(legend.box.spacing = unit(.1, "mm"), legend.spacing.y = unit(.1, 'cm'))
print(p)
dev.off()

getOutlierP <- function(x, y) {
  XY <- cbind(x, y)
  D2 <- apply(scale(XY, center = TRUE, scale = FALSE), 1, function(xy) {
    matrix(xy, nrow = 1) %*% solve(var(XY, na.rm = TRUE)) %*% matrix(xy, ncol = 1)
  })
  pchisq(D2, df = 2, lower = FALSE)
}
logFC_all[, P_outlier := getOutlierP(log2FC_archaic_vs_modern_HepG2_NS, log2FC_archaic_vs_modern_K562_NS)]

p <- p + geom_text_repel(data = logFC_all[!is.na(sharing_HepG2_K562) & sharing_HepG2_K562 != "non significant" & P_outlier < .01 | sharing_HepG2_K562 == "shared, opposite effect"], aes(x = log2FC_archaic_vs_modern_HepG2_NS, y = log2FC_archaic_vs_modern_K562_NS,  col = sharing_HepG2_K562, label = TargetGene), size=2, show.legend = FALSE) #,
p <- p + guides(text = FALSE)
pdf(sprintf("%s/comparison_effect_size/02b_logFC_HepG2_vs_K562_withnames.pdf", FIGURE_DIR),  height = 3.2, width = 5)
print(p)
dev.off()



########################################################################################################################
############################## COMPARISON BETWEEN EMVARs : A549 VS K562 ######################################################
########################################################################################################################

logFC_all <- dcast(emVARs_annot_obs[is_introgressed & !ctrl & !excluded & ANALYSIS_SUBTYPE == "celltype_cond" & COND_ID %in% c("K562_NS", "A549_NS"), ], posID + oligo1 + oligo2 ~ COND_ID, value.var = list("log2FC_archaic_vs_modern", "is_emVar_CRITERIA", "FDR", 'pval_emp')) # , fun.aggregate = list(mean, any, mean)
logFC_all <- merge(logFC_all, NearestGene, by = "posID" , all.x=TRUE)

emVARs_Diff_obs_celltype <- fread(sprintf("%s/all_emVars_diff_obs_celltype.tsv", EMVARDIFF_DIR))
emVARs_Diff_obs_A549_vs_K562 <- emVARs_Diff_obs_celltype[group2_labels == "A549-ACE2_NS_all" & group1_labels == "K562_NS_all", ]
logFC_all <- merge(logFC_all, emVARs_Diff_obs_A549_vs_K562[, .(posID, is_emVar_Diff_CRITERIA, Delta_log2FC_archaic_vs_modern, pval_interaction=pval_emp)], by = "posID")
logFC_all[, sharing_A549_K562 := case_when(
  (is_emVar_CRITERIA_A549_NS == TRUE | is_emVar_CRITERIA_K562_NS == TRUE) & pval_emp_A549_NS < 0.05 & pval_emp_K562_NS < 0.05 & sign(log2FC_archaic_vs_modern_A549_NS) == sign(log2FC_archaic_vs_modern_K562_NS) & pval_interaction>.05 ~ "shared, identical effect",
  (is_emVar_CRITERIA_A549_NS == TRUE | is_emVar_CRITERIA_K562_NS == TRUE) & pval_emp_A549_NS < 0.05 & pval_emp_K562_NS < 0.05 & sign(log2FC_archaic_vs_modern_A549_NS) == sign(log2FC_archaic_vs_modern_K562_NS) & pval_interaction<.05 & sign(Delta_log2FC_archaic_vs_modern) != sign(log2FC_archaic_vs_modern_A549_NS) ~ "shared, stronger in A549",
  (is_emVar_CRITERIA_A549_NS == TRUE | is_emVar_CRITERIA_K562_NS == TRUE) & pval_emp_A549_NS < 0.05 & pval_emp_K562_NS < 0.05 & sign(log2FC_archaic_vs_modern_A549_NS) == sign(log2FC_archaic_vs_modern_K562_NS) & pval_interaction<.05 & sign(Delta_log2FC_archaic_vs_modern) == sign(log2FC_archaic_vs_modern_A549_NS) ~ "shared, stronger in K562",
  (is_emVar_CRITERIA_A549_NS == TRUE & is_emVar_CRITERIA_K562_NS == TRUE) & pval_emp_A549_NS < 0.05 & pval_emp_K562_NS < 0.05 & sign(log2FC_archaic_vs_modern_A549_NS) != sign(log2FC_archaic_vs_modern_K562_NS) & pval_interaction<.05 ~ "shared, opposite effect",
  (is_emVar_CRITERIA_A549_NS == TRUE & pval_emp_K562_NS > 0.05 & pval_interaction<.05) ~ "A549 specific",
  (is_emVar_CRITERIA_A549_NS == TRUE & pval_emp_K562_NS > 0.05 & pval_interaction>.05) ~ "shared, identical effect",
  (is_emVar_CRITERIA_K562_NS == TRUE & pval_emp_A549_NS > 0.05  & pval_interaction<.05) ~ "K562 specific",
  (is_emVar_CRITERIA_K562_NS == TRUE & pval_emp_A549_NS > 0.05 & pval_interaction>.05) ~ "shared, identical effect",
  is_emVar_CRITERIA_A549_NS == FALSE & is_emVar_CRITERIA_K562_NS == FALSE ~ "non significant"
)]


color_code=c( "shared, identical effect" = "#d6c600",
  "shared, undetected in K562" = mergeCols("#d6c600",unname(color_celline_2levels_light["A549 weak"])),
  "shared, undetected in A549" = mergeCols("#d6c600",unname(color_celline_2levels_light["K562 weak"])),
  "shared, opposite effect" = "#DD4466",
  "shared, stronger in A549" = unname(color_celline_2levels_light["A549 weak"]),
  "shared, stronger in K562" = unname(color_celline_2levels_light["K562 weak"]),
  "A549 specific" = unname(color_celline_2levels_light["A549 strong"]),
  "K562 specific" = unname(color_celline_2levels_light["K562 strong"])
)

pdf(sprintf("%s/comparison_effect_size/02_logFC_A549_vs_K562_NS.pdf", FIGURE_DIR), height = 2.2, width = 4)
p <- ggplot(logFC_all[!is.na(sharing_A549_K562) & sharing_A549_K562 != "non significant"], aes(x = log2FC_archaic_vs_modern_A549_NS, y = log2FC_archaic_vs_modern_K562_NS))
p <- p + rasterize(geom_point(aes(col = factor(sharing_A549_K562,names(color_code))),alpha = .5), dpi = 200)
p <- p + xlab("log2 FC (A549)") + ylab("log2 FC (K562)")
p <- p + geom_hline(yintercept = 0, col = "grey") + geom_vline(xintercept = 0, col = "grey")
p <- p + theme_plot(rotate.x = 90, lpos = "right", fontsize = 11)
p <- p + geom_abline(intercept = 0, slope = 1, col = "lightgrey", linetype = 2)
p <- p + scale_color_manual(values = color_code)
p <- p +guides(col = guide_legend(ncol=1, byrow = FALSE)) + theme(legend.box.spacing = unit(.1, "mm"), legend.spacing.y = unit(.1, 'cm'))
print(p)
dev.off()

getOutlierP <- function(x, y) {
  XY <- cbind(x, y)
  D2 <- apply(scale(XY, center = TRUE, scale = FALSE), 1, function(xy) {
    matrix(xy, nrow = 1) %*% solve(var(XY, na.rm = TRUE)) %*% matrix(xy, ncol = 1)
  })
  pchisq(D2, df = 2, lower = FALSE)
}
logFC_all[, P_outlier := getOutlierP(log2FC_archaic_vs_modern_A549_NS, log2FC_archaic_vs_modern_K562_NS)]

p <- p + geom_text_repel(data = logFC_all[!is.na(sharing_A549_K562) & sharing_A549_K562 != "non significant" & P_outlier < .01 | sharing_A549_K562 == "shared, opposite effect"], aes(x = log2FC_archaic_vs_modern_A549_NS, y = log2FC_archaic_vs_modern_K562_NS,  col = sharing_A549_K562, label = TargetGene), size=2, show.legend = FALSE) #,
p <- p + guides(text = FALSE)
pdf(sprintf("%s/comparison_effect_size/02b_logFC_A549_vs_K562_withnames.pdf", FIGURE_DIR),  height = 3.2, width = 5)
print(p)
dev.off()




##### OLD but used as example for developement
######### UMI_perBC_DNA_vs_RNA
# DataPlot <- merge(RES_CRS, UMI_perCRS_DNA_vs_RNA_wide, by = "CRS")
# setnames(DataPlot, c("celltypeCalu3", "celltypeHepG2"), c("alpha_Calu3", "alpha_HepG2"))
# DataPlot[, FDR_celltype := p.adjust(Pvalue_celltypes, "fdr")]
# DataPlot[, class := case_when(
#   FDR_scrambled_Calu3 > 0.05 & FDR_scrambled_HepG2 < 0.05 & FDR_celltype < 0.05 ~ "HepG2 specific",
#   FDR_scrambled_Calu3 < 0.05 & FDR_scrambled_HepG2 > 0.05 & FDR_celltype < 0.05 ~ "Calu3 specific",
#   FDR_scrambled_Calu3 < 0.05 & FDR_scrambled_HepG2 < 0.05 & FDR_celltype > 0.05 ~ "shared, identical effect",
#   FDR_scrambled_Calu3 < 0.05 & FDR_scrambled_HepG2 < 0.05 & FDR_celltype < 0.05 & abs(OR_GCadj_Calu3) > abs(OR_GCadj_HepG2) & sign(log2(OR_GCadj_Calu3)) == sign(log2(OR_GCadj_HepG2)) ~ "shared, stronger in Calu3",
#   FDR_scrambled_Calu3 < 0.05 & FDR_scrambled_HepG2 < 0.05 & FDR_celltype < 0.05 & abs(OR_GCadj_HepG2) > abs(OR_GCadj_Calu3) & sign(log2(OR_GCadj_Calu3)) == sign(log2(OR_GCadj_HepG2)) ~ "shared, stronger in HepG2",
#   FDR_scrambled_Calu3 < 0.05 & FDR_scrambled_HepG2 < 0.05 & FDR_celltype < 0.05 & sign(log2(OR_GCadj_Calu3)) != sign(log2(OR_GCadj_HepG2)) ~ "shared, opposite effect",
#   FDR_scrambled_Calu3 > 0.05 & FDR_scrambled_HepG2 < 0.05 & FDR_celltype > 0.05 ~ "detected only in HepG2, likely shared",
#   FDR_scrambled_Calu3 < 0.05 & FDR_scrambled_HepG2 > 0.05 & FDR_celltype > 0.05 ~ "detected only in Calu3, likely shared",
#   FDR_scrambled_Calu3 > 0.05 & FDR_scrambled_HepG2 > 0.05 ~ "non significant"
# )]
# DataPlot[!is.na(alpha_Calu3), logalpha_GCadj_Calu3 := lm(log2(alpha_Calu3) ~ GC)$res]
# DataPlot[!is.na(alpha_HepG2), logalpha_GCadj_HepG2 := lm(log2(alpha_HepG2) ~ GC)$res]

# DataPlot[log2(OR_GCadj_HepG2) > 3 | abs(log2(OR_GCadj_Calu3)) > .5 & nBC_DNA_Calu3 > 50 & nBC_DNA_HepG2 > 50, ]

#############################################################################
#############################################################################
#############################################################################

cat("\nAll done\n")
q("no")
