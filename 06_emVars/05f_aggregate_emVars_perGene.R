MPRA_DIR <- "/pasteur/helix/projects/evo_immuno_pop/MPRA"
ANALYSIS_DIR <- "MPRA_count_exp6_analysis3"

source(sprintf("%s/scripts/%s/00_starter.R", MPRA_DIR, ANALYSIS_DIR))
source(sprintf("%s/scripts/%s/00_load_associations.R", MPRA_DIR, ANALYSIS_DIR))
source(sprintf("%s/scripts/misc_plots.R", MPRA_DIR))

RUN_ID <- "RUN3_Z2_nBC10"
CRITERIA_ACTIVE <- "lfdr20_scrambled1pct_GCnorm"
FILTER_ACTIVE <- TRUE
CRITERIA_EMVARS <- "EmVar_lfdr20_FC.2"

cmd <- commandArgs(trailingOnly = TRUE)
print(cmd)


for (i in seq_along(cmd)) {
  if (cmd[i] == "--run_dir" || cmd[i] == "-r") {
    RUN_ID <- cmd[i + 1]
  }
  if (cmd[i] == "--criteria_active" || cmd[i] == "-a") {
    CRITERIA_ACTIVE <- cmd[i + 1]
  }
  if (cmd[i] == "--criteria_emvar" || cmd[i] == "-e") {
    CRITERIA_EMVARS <- cmd[i + 1]
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
SNP_annot_v4 <- fread(sprintf("%s/data/%s/00_SNP_annot_v4.txt", MPRA_DIR, ANALYSIS_DIR))

# load annotation per POP
selected_annot <- fread(sprintf("%s/data/%s/selected_snp_annotation.tsv.gz", MPRA_DIR, ANALYSIS_DIR))

# load clean annotation of introgression
selected_annot_wide <- fread(sprintf("%s/data/%s/selected_snp_annotation_wide.tsv.gz", MPRA_DIR, ANALYSIS_DIR))
toc()



##### active parameters
source(sprintf("%s/scripts/%s/03_00_parameter_activity.R", MPRA_DIR, ANALYSIS_DIR))
##### emVars parameters
source(sprintf("%s/scripts/%s/05_00_parameter_emVars.R", MPRA_DIR, ANALYSIS_DIR))

IN_DIR <- sprintf("%s/data/%s/03a_aggMPRA_analyse/%s/", MPRA_DIR, ANALYSIS_DIR, RUN_ID)
ACTIVE_DIR <- sprintf("%s/data/%s/03b_computeFDR/%s/Activity/%s/", MPRA_DIR, ANALYSIS_DIR, RUN_ID, CRITERIA_ACTIVE)
EMVAR_DIR <- sprintf("%s/data/%s/03b_computeFDR/%s/emVars/%s/%s", MPRA_DIR, ANALYSIS_DIR, RUN_ID, CRITERIA_EMVARS, CRITERIA_ACTIVE_OUT)

dir.create(EMVAR_DIR, recursive = TRUE)
FIGURE_DIR <- sprintf("%s/figures/%s/%s/04_emVars/%s/%s", MPRA_DIR, ANALYSIS_DIR, RUN_ID, CRITERIA_EMVARS, CRITERIA_ACTIVE_OUT)

dir.create(FIGURE_DIR, recursive = TRUE)
dir.create(sprintf("%s/SupTables/", FIGURE_DIR), recursive = TRUE)


# load tested oligos
source(sprintf("%s/scripts/%s/02z__define_includedCRS.R", MPRA_DIR, ANALYSIS_DIR))

# emVARs_all_results <- fread( sprintf("%s/all_emVars_results.tsv.gz", IN_DIR), sep = "\t")

##########################################################################################################
################### TODO: see what is missing to run the code ############################################
##########################################################################################################

# TODO:
# clean the code to simplify,
# consider all TSS < 1Mb away & weight by distance
# add contacts from ABC or rE2G (with an evidence weight)
# compare with eQTL effect sizes
# combine both metrics with data-driven weights

oligo_activity_file <- sprintf("%s/oligo_activity__all.tsv.gz", ACTIVE_DIR)
oligo_activity <- fread(oligo_activity_file)
oligo_activity_obs <- fread(sprintf("%s/all_oligos_annotated.tsv.gz", ACTIVE_DIR)) # where does this come from	?

oligo_target <- fread(sprintf("%s/data/%s/00_oligo_targets_v4.txt", MPRA_DIR, ANALYSIS_DIR), sep = "\t")
oligo_target_collapsed <- fread(sprintf("%s/data/%s/00_oligo_targets_collapsed_v4.txt", MPRA_DIR, ANALYSIS_DIR), sep = "\t")

CS_results_finemap_eQTL_tested <- fread(sprintf("%s/data/%s/CredibleSets/CredibleSets_eQTLs.tsv.gz", MPRA_DIR, ANALYSIS_DIR))
CS_results_finemap_eQTL_tested <- merge(SNP_annot_v4[, .(ID, variantId_hg38, crsID, rsID, REF, ALT, INTROGRESSED.allele, allele2_is_REF)], CS_results_finemap_eQTL_tested, all.x = TRUE, by.x = "variantId_hg38", by.y = "variantId")


emVARs_obs_ctc <- fread(file = sprintf("%s/all_emVars_annotated_celltype__%s.tsv", EMVAR_DIR, CRITERIA_ACTIVE_OUT), sep = "\t")

########################################################################################################
######################################## eQTL catalog  only ############################################
########################################################################################################

eQTL_list <- CS_results_finemap_eQTL_tested[studyType == "eqtl" | studyType == "sceqtl", ][order(variantId_hg38, eQTL_gene, eQTL_tissue, pValueExponent)]
eQTL_list[, INTROGRESSED_BETA := ifelse(INTROGRESSED.allele == ALT, beta, ifelse(INTROGRESSED.allele == REF, -beta, NA))]
eQTL_list <- eQTL_list[!is.na(INTROGRESSED_BETA), ]

# emVARs_obs_ctc[,.(.N,Pct_blood=mean(!is.na(blood_eQTL) & blood_eQTL),Pct_lung=mean(!is.na(lung_eQTL) & lung_eQTL),Pct_liver=mean(!is.na(liver_eQTL) & liver_eQTL)),keyby=.(celline,is_emVar_CRITERIA)]
# emVARs_obs_ctc[!is.na(lung_eQTL),.(.N,Pct_K562=mean(is_emVar_CRITERIA[celline=='K562']),Pct_HEPG2=mean(is_emVar_CRITERIA[celline=='HepG2']),Pct_A549=mean(is_emVar_CRITERIA[celline=='A549'])),keyby=.(blood_eQTL,liver_eQTL,lung_eQTL)]

eQTL_vs_emVar <- merge(eQTL_list, emVARs_obs_ctc[, .(ID, is_emVar_CRITERIA, log2FC_archaic_vs_modern, celline, condition, COND_ID, Nearest)], by = "ID", allow.cartesian = TRUE, all.y = TRUE)
eQTL_vs_emVar[, mean(sign(INTROGRESSED_BETA) == sign(log2FC_archaic_vs_modern)), by = COND_ID]
eQTL_vs_emVar[order(Nearest, -abs(log2FC_archaic_vs_modern))][!duplicated(paste(Nearest, COND_ID))][, cor(INTROGRESSED_BETA, log2FC_archaic_vs_modern), by = COND_ID]

has_eQTL <- eQTL_vs_emVar[is_emVar_CRITERIA == TRUE, .(
  has_eQTL = any(!is.na(eQTL_study)),
  has_same_direction_eQTL = any(!is.na(eQTL_study) & sign(INTROGRESSED_BETA) == sign(log2FC_archaic_vs_modern)),
  has_opposite_direction_eQTL = any(!is.na(eQTL_study) & sign(INTROGRESSED_BETA) != sign(log2FC_archaic_vs_modern))
), by = .(ID, COND_ID)]
has_eQTL_all <- eQTL_vs_emVar[is_emVar_CRITERIA == TRUE, .(
  has_eQTL = any(!is.na(eQTL_study)),
  has_same_direction_eQTL = any(!is.na(eQTL_study) & sign(INTROGRESSED_BETA) == sign(log2FC_archaic_vs_modern)),
  has_opposite_direction_eQTL = any(!is.na(eQTL_study) & sign(INTROGRESSED_BETA) != sign(log2FC_archaic_vs_modern))
), by = .(ID)]

cat("Percentage emVars with an eQTL (any eQTL catalog):\n")


cat("alltogether\n")
has_eQTL_all[, mean(has_eQTL)] #

cat("by tissue\n")
has_eQTL[, mean(has_eQTL), by = COND_ID]

cat("detail by condition\n")
has_eQTL_all[, .(
  has_same_direction_eQTL = sum(has_same_direction_eQTL),
  has_opposite_direction_eQTL = sum(has_opposite_direction_eQTL),
  has_same_direction_eQTL_only = sum(has_same_direction_eQTL & !has_opposite_direction_eQTL),
  has_opposite_direction_eQTL_only = sum(has_opposite_direction_eQTL & !has_same_direction_eQTL),
  has_same_and_opposite_direction_eQTL = sum(has_opposite_direction_eQTL & has_same_direction_eQTL),
  has_any_eQTL = sum(has_eQTL)
)]


########################################################################################################
######################################## GTEx only #####################################################
########################################################################################################

eQTL_list <- CS_results_finemap_eQTL_tested[eQTL_study == "gtex" & eQTL_type == "ge", ][order(variantId_hg38, eQTL_gene, eQTL_tissue, pValueExponent)]
eQTL_list[, INTROGRESSED_BETA := ifelse(INTROGRESSED.allele == ALT, beta, ifelse(INTROGRESSED.allele == REF, -beta, NA))]
eQTL_list <- eQTL_list[!is.na(INTROGRESSED_BETA), ]


cat("Percentage emVars with an eQTL (any GTEx only):\n")
cat("alltogether\n")
has_eQTL_all[, mean(has_eQTL)] #

cat("by tissue\n")
has_eQTL[, mean(has_eQTL), by = COND_ID]

cat("detail by condition\n")
has_eQTL_all[, .(
  has_same_direction_eQTL = sum(has_same_direction_eQTL),
  has_opposite_direction_eQTL = sum(has_opposite_direction_eQTL),
  has_same_direction_eQTL_only = sum(has_same_direction_eQTL & !has_opposite_direction_eQTL),
  has_opposite_direction_eQTL_only = sum(has_opposite_direction_eQTL & !has_same_direction_eQTL),
  has_same_and_opposite_direction_eQTL = sum(has_opposite_direction_eQTL & has_same_direction_eQTL),
  has_any_eQTL = sum(has_eQTL)
)]


# best_eQTL_list <- eQTL_list[,head(.SD,1),by=.(variantId_hg38)]
ctr_emVARs_eQTL <- emVARs_all_results[grepl("eQTL", oligo1) & power == 0 & pval.LRT < .01 & ANALYSIS_SUBTYPE == "celltype_cond" & perm == "perm_0_0_0", ]
# ctr_emVARs_eQTL=merge(ctr_emVARs_eQTL,SNP_annot_v4[,.(ID,variantId_hg38,crsID, rsID, REF, ALT,INTROGRESSED.allele)],by='crsID',allow.cartesian=TRUE)
eQTL_vs_emVar <- merge(eQTL_list, ctr_emVARs_eQTL, by = "crsID", allow.cartesian = TRUE)

ctr_emVARs_eQTL_2 <- ctr_emVARs_eQTL[, .(Pos_emVar = sum(logFC_2vs1 > 0), Neg_emVar = sum(logFC_2vs1 < 0)), keyby = .(crsID, oligo1, oligo2)]
eQTL_list_2 <- eQTL_list[crsID %in% ctr_emVARs_eQTL_2$crsID, .(Pos_eQTL = sum(beta > 0), Neg_eQTL = sum(beta < 0)), keyby = .(crsID, REF, ALT)]

eQTL_vs_emVar[, .(Pos_emVar = sum((logFC_2vs1 * ifelse(allele2_is_REF, -1, 1)) > 0), Neg_emVar = sum((logFC_2vs1 * ifelse(allele2_is_REF, -1, 1)) < 0), Pos_eQTL = sum(beta > 0), Neg_eQTL = sum(beta < 0)), keyby = .(crsID, REF, ALT)] <- merge(eQTL_list, ctr_emVARs_eQTL, by = "crsID", allow.cartesian = TRUE)

Pct_concordant_GTex <- eQTL_vs_emVar[is_emVar_CRITERIA == TRUE, .(nSNP = .N, nCOncordant_effect = sum(sign(INTROGRESSED_BETA) == sign(log2FC_archaic_vs_modern)), Pct_concondant = mean(sign(INTROGRESSED_BETA) == sign(log2FC_archaic_vs_modern))), by = COND_ID]
Pct_concordant_GTex[, P := binom.test(nCOncordant_effect, nSNP)$p.value, by = COND_ID]

Pct_concordant_GTex


cat("All done !\n")
q("no")


########################################################################################################
######################################## tests and aggregation per gene (deprecated) ###################
########################################################################################################



# emVARs_obs_ctc_full <- fread(sprintf("%s/all_emVars_annotated_full_celltype__%s.tsv", EMVAR_DIR, CRITERIA_ACTIVE_OUT))

# caQTL=list();for(file in dir('/Volumes/evo_immuno_pop/MPRA/data/caQTL_ENCODE',full=T)){caQTL[[file]]=fread(file)}


caQTL <- fread(sprintf("%s/data/caQTL_ENCODE/caQTL_TableS12_PMID34038741.csv", MPRA_DIR), skip = 3)
caQTL[lead_var_ID %in% SNP_annot$rsID, ]

eQTLgen <- fread("/pasteur/helix/projects/evo_immuno_pop/single_cell/resources/references/eQTLgen/2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz")
eQTLgen <- eQTLgen[SNP %in% SNP_annot$rsID, ]
eQTL_vs_emVar <- merge(eQTLgen, emVARs_obs_ctc[, .(rsID, REF, ALT, INTROGRESSED.allele, is_emVar_CRITERIA, log2FC_archaic_vs_modern, celline, condition, COND_ID, Nearest)], by.x = "SNP", by.y = "rsID", allow.cartesian = TRUE)
eQTL_vs_emVar <- eQTL_vs_emVar[INTROGRESSED.allele != "", ]
eQTL_vs_emVar[, Zscore := ifelse(AssessedAllele == INTROGRESSED.allele, Zscore, -Zscore)]

crs_targets <- unique(oligo_target_collapsed[, -"oligo"])
oligo_target_plus <- oligo_target[grepl("\\+", method)]
oligo_target_plus[, method := gsub("\\+", "", method)]
crs_target_v2 <- rbind(oligo_target, oligo_target_plus)
crs_target_v2 <- unique(crs_target_v2[, -"oligo"])
crs_target_v2[, .(N = length(unique(crsID))), by = .(TargetGene, method, celltype)]
crs_target_v2[, .(N = length(unique(crsID))), by = .(TargetGene, method, celltype)][, .N, keyby = .(method, celltype, pmin(N, 10))]
crs_target_v2[celltype == "any", method := paste0(method, "_any")]
dispatch <- data.table(celltype = "any", new_celltype = c("HepG2", "K562", "A549"))
crs_target_any <- merge(crs_target_v2[celltype == "any", ], dispatch, by = "celltype", allow.cartesian = TRUE, all.x = TRUE)
crs_target_any[, celltype := new_celltype]
crs_target_any[, new_celltype := NULL]
crs_target_v2 <- rbind(crs_target_v2[celltype != "any", ], crs_target_any)
crs_target_v2 <- unique(crs_target_v2, by = c("crsID", "method", "celltype", "TargetGene"))

Fig_data <- crs_target_v2[, .N, by = .(TargetGene, method, celltype)][, .N, keyby = .(method, celltype, N_tested_SNP = pmin(N, 15))]
p <- ggplot(Fig_data, aes(x = N_tested_SNP, y = N, fill = method)) +
  geom_bar(stat = "Identity", position = "dodge", border = "white", stroke = .1) +
  theme_plot(rotate.x = 90)
p <- p + facet_grid(cols = vars(celltype))
pdf(sprintf("%s/04_emVars/3_emVARs_agg_by_target/00_number_of_targets_perGene.pdf", FIGURE_DIR), height = 3, width = 7)
print(p)
dev.off()

p <- ggplot(Fig_data, aes(x = N_tested_SNP, y = N, fill = celltype)) +
  geom_bar(stat = "Identity", position = "dodge", border = "white", stroke = .1) +
  theme_plot(rotate.x = 90)
p <- p + facet_wrap(~ factor(method, c("NearestGene_any", "NearestTSS_any", "NearestTSS", "ABC_any", "ABC+_any", "HiC", "rE2G", "rE2G+"))) + scale_fill_manual(values = color_celline) + scale_y_continuous(transform = "sqrt", breaks = c(2000, 1000, 500, 200, 100, 50, 20, 0))
pdf(sprintf("%s/04_emVars/3_emVARs_agg_by_target/00_number_of_targets_perGene_split_method.pdf", FIGURE_DIR), height = 5, width = 5)
print(p)
dev.off()

emVARs_obs_ctc_targets <- merge(emVARs_obs_ctc, crs_target_v2, by = c("crsID", "celltype"), allow.cartesian = TRUE, all.x = TRUE)

# emVARs_obs_ctc_targets[,activity_archaic:=ifelse(INTROGRESSED.allele==allele.1,a1,a2)]
# emVARs_obs_ctc_targets[,activity_modern:=ifelse(INTROGRESSED.allele==allele.1,a2,a1)]

# #emVARs_annot_active_ok_targets <- merge(emVARs_annot_active_ok_targets,merge(oligo_activity_obs[,.(oligo_arcahaic=oligo,ANALYSIS_NAME,alpha_2=alpha.GC_norm)]),by='oligo',allow.cartesian=TRUE,all.x=TRUE)
# emVARs_obs_ctc_targets[,mean_activity:=(log2(a1)+log2(a2))/2]
# emVARs_obs_ctc_targets[,Score1_norm:=ifelse(Score1_type=='negative distance in kb',exp(Score1/1000),Score1)/sum(ifelse(Score1_type=='negative distance in kb',exp(Score1/1000),Score1)),by=.(TargetGene, method,ANALYSIS_NAME)]

compute_haplotype_effects <- function(alpha_dist, alpha_ABC) {
  alpha_dist <- 1
  alpha_ABC <- 1
  TargetGene_effect <- emVARs_obs_ctc_targets[method %in% c("ABC_any", "NearestGene_any"), .(beta_haplotype = sum(log2FC_archaic_vs_modern * ifelse(Score1_type == "negative distance in kb", exp(alpha_dist * Score1 / 1000), alpha_ABC * Score1) * (1 - lfdr))),
    by = .(TargetGene, COND_ID = gsub("-ACE2", "", gsub("_all", "", ANALYSIS_NAME)))
  ]
  merge(TargetGene_effect, best_eQTL_list)
}



emVARs_obs_ctc_targets <- emVARs_obs_ctc_targets[, .(
  archaic_sum = sum(activity_archaic),
  modern_sum = sum(activity_modern),
  haplo_sum = sum(log2FC_archaic_vs_modern),
  haplo_sum.se = sqrt(sum(log2FC.se^2)),
  haplo_count_pos = sum(log2FC_archaic_vs_modern > 0),
  haplo_count_total = .N,
  haplo_count_pos_FDR5 = sum(log2FC_archaic_vs_modern[FDR < 0.05] > 0),
  haplo_count_total_FDR5 = sum(FDR < 0.05),
  haplo_count_crs = length(unique(crsID)),
  haplo_sum_weighted = sum(log2FC_archaic_vs_modern * ifelse(Score1_type == "negative distance in kb", exp(Score1 / 1000), Score1)),
  haplo_sum_weighted.se = sqrt(sum((ifelse(Score1_type == "negative distance in kb", exp(Score1 / 1000), Score1) * log2FC.se)^2)),
  haplo_sum_weightedLFDR = sum(log2FC_archaic_vs_modern * ifelse(Score1_type == "negative distance in kb", exp(Score1 / 1000), Score1) * (1 - lfdr)),
  haplo_sum_weightedLFDR.se = sqrt(sum((ifelse(Score1_type == "negative distance in kb", exp(Score1 / 1000), Score1) * (1 - lfdr) * log2FC.se)^2)),
  haplo_sum_weightedLFDR_norm = sum(log2FC_archaic_vs_modern * Score1_norm * (1 - lfdr)),
  haplo_sum_weightedLFDR_norm.se = sqrt(sum((Score1_norm * (1 - lfdr) * log2FC.se)^2)),
  haplo_sum_FDR5 = sum(log2FC_archaic_vs_modern * (FDR < .05)),
  haplo_sum_FDR5.se = sqrt(sum((FDR < .05) * log2FC.se)^2),
  haplo_sum_LFDR = sum(log2FC_archaic_vs_modern * (1 - lfdr)),
  haplo_sum_LFDR.se = sqrt(sum((1 - lfdr) * log2FC.se)^2),
  haplo_sum_weighted_activity = sum(log2FC_archaic_vs_modern * abs(mean_activity)),
  haplo_sum_weighted_activity.se = sqrt(sum(abs(mean_activity) * log2FC.se)^2),
  haplo_sum_weighted_activity_FDR5 = sum(log2FC_archaic_vs_modern * abs(mean_activity) * (FDR < .05)),
  haplo_sum_weighted_activity_FDR5.se = sqrt(sum(abs(mean_activity) * log2FC.se * (FDR < .05))^2),
  haplo_sum_weighted_activity_LFDR = sum(log2FC_archaic_vs_modern * abs(mean_activity) * (1 - lfdr)),
  haplo_sum_weighted_activity_LFDR.se = sqrt(sum(abs(mean_activity) * log2FC.se * (1 - lfdr))^2)
), by = .(TargetGene, method, COND_ID = gsub("-ACE2", "", gsub("_all", "", ANALYSIS_NAME)))]

emVARs_agg_by_target[, haplo_sum_Z := haplo_sum / haplo_sum.se]
emVARs_agg_by_target[, haplo_sum_weighted_Z := haplo_sum_weighted / haplo_sum_weighted.se]
emVARs_agg_by_target[, haplo_sum_FDR_Z := haplo_sum_FDR / haplo_sum_FDR.se]
emVARs_agg_by_target[, haplo_sum_weighted_activity_Z := haplo_sum_weighted_activity / haplo_sum_weighted_activity.se]
emVARs_agg_by_target[, haplo_sum_weighted_activity_FDR_Z := haplo_sum_weighted_activity_FDR / haplo_sum_weighted_activity_FDR.se]
emVARs_agg_by_target[, haplo_sum_weightedLFDR_Z := haplo_sum_weightedLFDR / haplo_sum_weightedLFDR.se]

emVARs_agg_by_target[, haplo_count_P := binom.test(haplo_count_pos, haplo_count_total)$p.value, by = .(TargetGene, method, COND_ID)]
emVARs_agg_by_target[, haplo_count_Z := qnorm(haplo_count_P / 2) * sign(haplo_count_pos - haplo_count_total / 2), by = .(TargetGene, method, COND_ID)]
emVARs_agg_by_target[, haplo_count_FDR5_P := 1]
emVARs_agg_by_target[, haplo_count_FDR5_Z := 0]
emVARs_agg_by_target[haplo_count_total_FDR5 > 0, haplo_count_FDR5_P := binom.test(haplo_count_pos_FDR5, haplo_count_total_FDR5)$p.value, by = .(TargetGene, method, COND_ID)]
emVARs_agg_by_target[haplo_count_total_FDR5 > 0, haplo_count_FDR5_Z := qnorm(haplo_count_FDR5_P / 2) * sign(haplo_count_pos_FDR5 - haplo_count_total_FDR5 / 2), by = .(TargetGene, method, COND_ID)]


emVARs_agg_by_target[, haplo_sum_Z0 := qnorm(rank(haplo_sum_Z, ties = "random") / (.N + 1), lower = TRUE), by = .(method, COND_ID)]
emVARs_agg_by_target[, haplo_sum_weighted_Z0 := qnorm(rank(haplo_sum_weighted_Z, ties = "random") / (.N + 1), lower = TRUE), by = .(method, COND_ID)]
emVARs_agg_by_target[, haplo_count_Z0 := qnorm(rank(haplo_count_Z, ties = "random") / (.N + 1), lower = TRUE), by = .(method, COND_ID)]
emVARs_agg_by_target[, haplo_count_FDR5_Z0 := qnorm(rank(haplo_count_FDR5_Z, ties = "random") / (.N + 1), lower = TRUE), by = .(method, COND_ID)]
emVARs_agg_by_target[, haplo_sum_weightedLFDR_Z0 := qnorm(rank(haplo_sum_weightedLFDR_Z, ties = "random") / (.N + 1), lower = TRUE), by = .(method, COND_ID)]


dir.create(sprintf("%s/04_emVars/3_emVARs_agg_by_target", FIGURE_DIR), showWarnings = FALSE, recursive = TRUE)
p <- ggplot(emVARs_agg_by_target[COND_ID %chin% names(color_setup_simplified_norep), ], aes(x = haplo_sum_Z0, y = haplo_sum_Z, col = COND_ID)) +
  geom_point(size = .05) +
  geom_abline(intercept = 0, slope = 1, col = "lightgrey") +
  facet_grid(cols = vars(method))
p <- p + theme_plot() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + scale_color_manual(values = color_setup_simplified_norep)
pdf(sprintf("%s/04_emVars/3_emVARs_agg_by_target/QQplot__haplo_sum_Z.pdf", FIGURE_DIR), height = 3, width = 7)
print(p)
dev.off()

p <- ggplot(emVARs_agg_by_target[COND_ID %chin% names(color_setup_simplified_norep), ], aes(x = haplo_sum_weighted_Z0, y = haplo_sum_weighted_Z, col = COND_ID)) +
  geom_point(size = .05) +
  geom_abline(intercept = 0, slope = 1, col = "lightgrey") +
  facet_grid(cols = vars(method))
p <- p + theme_plot() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + scale_color_manual(values = color_setup_simplified_norep)
pdf(sprintf("%s/04_emVars/3_emVARs_agg_by_target/QQplot__haplo_sum_weighted_Z.pdf", FIGURE_DIR), height = 3, width = 7)
print(p)
dev.off()


p <- ggplot(emVARs_agg_by_target[COND_ID %chin% names(color_setup_simplified_norep), ], aes(x = haplo_count_Z0, y = haplo_count_Z, col = COND_ID)) +
  geom_point(size = .05) +
  geom_abline(intercept = 0, slope = 1, col = "lightgrey") +
  facet_grid(cols = vars(method), rows = vars())
p <- p + theme_plot() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + scale_color_manual(values = color_setup_simplified_norep)
pdf(sprintf("%s/04_emVars/3_emVARs_agg_by_target/QQplot___haplo_count_Z.pdf", FIGURE_DIR), height = 3, width = 7)
print(p)
dev.off()

p <- ggplot(emVARs_agg_by_target[COND_ID %chin% names(color_setup_simplified_norep) & haplo_count_total > 5, ], aes(x = haplo_sum_Z0, y = haplo_sum_Z, col = COND_ID)) +
  geom_point(size = .05) +
  geom_abline(intercept = 0, slope = 1, col = "lightgrey") +
  facet_grid(cols = vars(method))
p <- p + theme_plot() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + scale_color_manual(values = color_setup_simplified_norep)
pdf(sprintf("%s/04_emVars/3_emVARs_agg_by_target/QQplot_Nsnp_over_5__haplo_sum_Z.pdf", FIGURE_DIR), height = 3, width = 7)
print(p)
dev.off()

p <- ggplot(emVARs_agg_by_target[COND_ID %chin% names(color_setup_simplified_norep) & haplo_count_total > 5, ], aes(x = haplo_sum_weighted_Z0, y = haplo_sum_weighted_Z, col = COND_ID)) +
  geom_point(size = .05) +
  geom_abline(intercept = 0, slope = 1, col = "lightgrey") +
  facet_grid(cols = vars(method))
p <- p + theme_plot() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + scale_color_manual(values = color_setup_simplified_norep)
pdf(sprintf("%s/04_emVars/3_emVARs_agg_by_target/QQplot_Nsnp_over_5__haplo_sum_weighted_Z.pdf", FIGURE_DIR), height = 3, width = 7)
print(p)
dev.off()


p <- ggplot(emVARs_agg_by_target[COND_ID %chin% names(color_setup_simplified_norep) & haplo_count_total > 5, ], aes(x = haplo_sum_weightedLFDR_Z0, y = haplo_sum_weightedLFDR_Z, col = COND_ID)) +
  geom_point(size = .05) +
  geom_abline(intercept = 0, slope = 1, col = "lightgrey") +
  facet_grid(cols = vars(method))
p <- p + theme_plot() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + scale_color_manual(values = color_setup_simplified_norep)
pdf(sprintf("%s/04_emVars/3_emVARs_agg_by_target/QQplot_Nsnp_over_5__haplo_sum_weighted_LFDR_Z.pdf", FIGURE_DIR), height = 3, width = 7)
print(p)
dev.off()
p <- ggplot(emVARs_agg_by_target[COND_ID %chin% names(color_setup_simplified_norep), ], aes(x = haplo_sum_weightedLFDR_Z0, y = haplo_sum_weightedLFDR_Z, col = COND_ID)) +
  geom_point(size = .05) +
  geom_abline(intercept = 0, slope = 1, col = "lightgrey") +
  facet_grid(cols = vars(method))
p <- p + theme_plot() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + scale_color_manual(values = color_setup_simplified_norep)
pdf(sprintf("%s/04_emVars/3_emVARs_agg_by_target/QQplot__haplo_sum_weighted_LFDR_Z.pdf", FIGURE_DIR), height = 3, width = 7)
print(p)
dev.off()

p <- ggplot(emVARs_agg_by_target[COND_ID %chin% names(color_setup_simplified_norep) & haplo_count_total > 5, ], aes(x = haplo_count_Z0, y = haplo_count_Z, col = COND_ID)) +
  geom_point(size = .05) +
  geom_abline(intercept = 0, slope = 1, col = "lightgrey") +
  facet_grid(cols = vars(method))
p <- p + theme_plot() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + scale_color_manual(values = color_setup_simplified_norep)
pdf(sprintf("%s/04_emVars/3_emVARs_agg_by_target/QQplot_Nsnp_over_5__haplo_count_Z.pdf", FIGURE_DIR), height = 3, width = 7)
print(p)
dev.off()


p <- ggplot(emVARs_agg_by_target[COND_ID %chin% names(color_setup_simplified_norep) & haplo_count_total > 5, ], aes(x = haplo_sum_weighted_Z0, y = haplo_sum_weighted_Z, col = COND_ID)) +
  geom_point(size = .05) +
  geom_abline(intercept = 0, slope = 1, col = "lightgrey") +
  facet_grid(cols = vars(method))
p <- p + theme_plot() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + scale_color_manual(values = color_setup_simplified_norep)
pdf(sprintf("%s/04_emVars/3_emVARs_agg_by_target/QQplot_Nsnp_over_5__haplo_sum_weighted_Z.pdf", FIGURE_DIR), height = 3, width = 7)
print(p)
dev.off()


#################### haplo_sum_weighted_activity_FDR
locus_genes <- c("OAS1", "OAS2", "OAS3")

p <- ggplot(emVARs_agg_by_target[COND_ID %chin% names(color_setup_simplified_norep) & TargetGene %chin% locus_genes, ]) +
  geom_hline(yintercept = 0, col = "lightgrey")
p <- p + geom_pointrange(aes(
  x = factor(COND_ID, names(color_setup_simplified_norep)),
  y = haplo_sum_weighted_activity_FDR,
  ymin = haplo_sum_weighted_activity_FDR - 1.96 * haplo_sum_weighted_activity_FDR.se,
  ymax = haplo_sum_weighted_activity_FDR + 1.96 * haplo_sum_weighted_activity_FDR.se,
  col = COND_ID
), size = .1) + facet_grid(rows = vars(method), cols = vars(TargetGene))
p <- p + theme_plot(rotate.x = 90) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + scale_color_manual(values = color_setup_simplified_norep) + guides(col = "none")
pdf(sprintf("%s/04_emVars/3_emVARs_agg_by_target/Effect_size__OAS_genes__haplo_sum_weighted_activity_FDR.pdf", FIGURE_DIR), height = 7, width = 7)
print(p)
dev.off()

locus_genes <- c("TLR1", "TLR6", "TLR10")


p <- ggplot(emVARs_agg_by_target[COND_ID %chin% names(color_setup_simplified_norep) & TargetGene %chin% locus_genes, ]) +
  geom_hline(yintercept = 0, col = "lightgrey")
p <- p + geom_pointrange(aes(
  x = factor(COND_ID, names(color_setup_simplified_norep)),
  y = haplo_sum_weighted_activity_FDR,
  ymin = haplo_sum_weighted_activity_FDR - 1.96 * haplo_sum_weighted_activity_FDR.se,
  ymax = haplo_sum_weighted_activity_FDR + 1.96 * haplo_sum_weighted_activity_FDR.se,
  col = COND_ID
), size = .1) + facet_grid(rows = vars(method), cols = vars(TargetGene))
p <- p + theme_plot(rotate.x = 90) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + scale_color_manual(values = color_setup_simplified_norep) + guides(col = "none")
pdf(sprintf("%s/04_emVars/3_emVARs_agg_by_target/Effect_size__TLR1_genes__haplo_sum_weighted_activity_FDR.pdf", FIGURE_DIR), height = 7, width = 7)
print(p)
dev.off()

#################### haplo_sum_LFDR

locus_genes <- c("OAS1", "OAS2", "OAS3")
p <- ggplot(emVARs_agg_by_target[COND_ID %chin% names(color_setup_simplified_norep) & TargetGene %chin% locus_genes, ]) +
  geom_hline(yintercept = 0, col = "lightgrey")
p <- p + geom_pointrange(aes(x = factor(COND_ID, names(color_setup_simplified_norep)), y = haplo_sum_LFDR, ymin = haplo_sum_LFDR - 1.96 * haplo_sum_LFDR.se, ymax = haplo_sum_LFDR + 1.96 * haplo_sum_LFDR.se, col = COND_ID), size = .1)
p <- p + facet_grid(rows = vars(method), cols = vars(TargetGene))
p <- p + theme_plot(rotate.x = 90) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + scale_color_manual(values = color_setup_simplified_norep) + guides(col = "none")
pdf(sprintf("%s/04_emVars/3_emVARs_agg_by_target/Effect_size__OAS_genes__haplo_sum_LFDR.pdf", FIGURE_DIR), height = 7, width = 7)
print(p)
dev.off()

locus_genes <- c("TLR1", "TLR6", "TLR10")
p <- ggplot(emVARs_agg_by_target[COND_ID %chin% names(color_setup_simplified_norep) & TargetGene %chin% locus_genes, ]) +
  geom_hline(yintercept = 0, col = "lightgrey")
p <- p + geom_pointrange(aes(x = factor(COND_ID, names(color_setup_simplified_norep)), y = haplo_sum_LFDR, ymin = haplo_sum_LFDR - 1.96 * haplo_sum_LFDR.se, ymax = haplo_sum_LFDR + 1.96 * haplo_sum_LFDR.se, col = COND_ID), size = .1)
p <- p + facet_grid(rows = vars(method), cols = vars(TargetGene))
p <- p + theme_plot(rotate.x = 90) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + scale_color_manual(values = color_setup_simplified_norep) + guides(col = "none")
pdf(sprintf("%s/04_emVars/3_emVARs_agg_by_target/Effect_size__TLR1_genes__haplo_sum_LFDR.pdf", FIGURE_DIR), height = 7, width = 7)
print(p)
dev.off()



#################### haplo_sum_weightedLFDR


locus_genes <- c("OAS1", "OAS2", "OAS3")
p <- ggplot(emVARs_agg_by_target[COND_ID %chin% names(color_setup_simplified_norep) & TargetGene %chin% locus_genes, ]) +
  geom_hline(yintercept = 0, col = "lightgrey")
p <- p + geom_pointrange(aes(x = factor(COND_ID, names(color_setup_simplified_norep)), y = haplo_sum_weightedLFDR, ymin = haplo_sum_weightedLFDR - 1.96 * haplo_sum_weightedLFDR.se, ymax = haplo_sum_LFDR + 1.96 * haplo_sum_weightedLFDR.se, col = COND_ID), size = .1)
p <- p + facet_grid(rows = vars(method), cols = vars(TargetGene))
p <- p + theme_plot(rotate.x = 90) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + scale_color_manual(values = color_setup_simplified_norep) + guides(col = "none")
pdf(sprintf("%s/04_emVars/3_emVARs_agg_by_target/Effect_size__OAS_genes__haplo_sum_LFDR_weighted.pdf", FIGURE_DIR), height = 7, width = 7)
print(p)
dev.off()


locus_genes <- c("OAS1", "OAS2", "OAS3")
p <- ggplot(emVARs_agg_by_target[COND_ID %chin% names(color_setup_simplified_norep) & TargetGene %chin% locus_genes, ]) +
  geom_hline(yintercept = 0, col = "lightgrey")
p <- p + geom_pointrange(aes(x = factor(COND_ID, names(color_setup_simplified_norep)), y = haplo_sum_weightedLFDR, ymin = haplo_sum_weightedLFDR - 1.96 * haplo_sum_weightedLFDR.se, ymax = haplo_sum_LFDR + 1.96 * haplo_sum_weightedLFDR.se, col = COND_ID), size = .1)
p <- p + facet_grid(rows = vars(method), cols = vars(TargetGene)) + coord_cartesian(ylim = c(-2, 2))
p <- p + theme_plot(rotate.x = 90) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + scale_color_manual(values = color_setup_simplified_norep) + guides(col = "none")
pdf(sprintf("%s/04_emVars/3_emVARs_agg_by_target/Effect_size__OAS_genes__haplo_sum_LFDR_weighted_2.pdf", FIGURE_DIR), height = 7, width = 7)
print(p)
dev.off()

locus_genes <- c("TLR1", "TLR6", "TLR10")
p <- ggplot(emVARs_agg_by_target[COND_ID %chin% names(color_setup_simplified_norep) & TargetGene %chin% locus_genes, ]) +
  geom_hline(yintercept = 0, col = "lightgrey")
p <- p + geom_pointrange(aes(x = factor(COND_ID, names(color_setup_simplified_norep)), y = haplo_sum_weightedLFDR, ymin = haplo_sum_weightedLFDR - 1.96 * haplo_sum_weightedLFDR.se, ymax = haplo_sum_weightedLFDR + 1.96 * haplo_sum_weightedLFDR.se, col = COND_ID), size = .1)
p <- p + facet_grid(rows = vars(method), cols = vars(TargetGene))
p <- p + theme_plot(rotate.x = 90) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + scale_color_manual(values = color_setup_simplified_norep) + guides(col = "none")
pdf(sprintf("%s/04_emVars/3_emVARs_agg_by_target/Effect_size__TLR1_genes__haplo_sum_LFDR_weighted.pdf", FIGURE_DIR), height = 7, width = 7)
print(p)
dev.off()


locus_genes <- c("PAN2", "STAT2", "IL23A", "APOF")
p <- ggplot(emVARs_agg_by_target[COND_ID %chin% names(color_setup_simplified_norep) & TargetGene %chin% locus_genes, ]) +
  geom_hline(yintercept = 0, col = "lightgrey")
p <- p + geom_pointrange(aes(x = factor(COND_ID, names(color_setup_simplified_norep)), y = haplo_sum_weightedLFDR, ymin = haplo_sum_weightedLFDR - 1.96 * haplo_sum_weightedLFDR.se, ymax = haplo_sum_weightedLFDR + 1.96 * haplo_sum_weightedLFDR.se, col = COND_ID), size = .1)
p <- p + facet_grid(rows = vars(method), cols = vars(TargetGene))
p <- p + theme_plot(rotate.x = 90) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + scale_color_manual(values = color_setup_simplified_norep) + guides(col = "none")
pdf(sprintf("%s/04_emVars/3_emVARs_agg_by_target/Effect_size__STAT2_genes__haplo_sum_LFDR_weighted.pdf", FIGURE_DIR), height = 7, width = 7)
print(p)
dev.off()

p <- ggplot(emVARs_agg_by_target[COND_ID %chin% names(color_setup_simplified_norep) & TargetGene %chin% locus_genes & method != "HiC", ]) +
  geom_hline(yintercept = 0, col = "lightgrey")
p <- p + geom_pointrange(aes(x = factor(COND_ID, names(color_setup_simplified_norep)), y = haplo_sum_weightedLFDR, ymin = haplo_sum_weightedLFDR - 1.96 * haplo_sum_weightedLFDR.se, ymax = haplo_sum_weightedLFDR + 1.96 * haplo_sum_weightedLFDR.se, col = COND_ID), size = .1)
p <- p + facet_grid(rows = vars(method), cols = vars(TargetGene))
p <- p + theme_plot(rotate.x = 90) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + scale_color_manual(values = color_setup_simplified_norep) + guides(col = "none")
pdf(sprintf("%s/04_emVars/3_emVARs_agg_by_target/Effect_size__STAT2_genes__haplo_sum_LFDR_weighted_noHiC.pdf", FIGURE_DIR), height = 7, width = 7)
print(p)
dev.off()

locus_genes <- c("STAT3", "STAT5A", "STAT5B")
p <- ggplot(emVARs_agg_by_target[COND_ID %chin% names(color_setup_simplified_norep) & TargetGene %chin% locus_genes & method != "HiC", ]) +
  geom_hline(yintercept = 0, col = "lightgrey")
p <- p + geom_pointrange(aes(x = factor(COND_ID, names(color_setup_simplified_norep)), y = haplo_sum_weightedLFDR, ymin = haplo_sum_weightedLFDR - 1.96 * haplo_sum_weightedLFDR.se, ymax = haplo_sum_weightedLFDR + 1.96 * haplo_sum_weightedLFDR.se, col = COND_ID), size = .1)
p <- p + facet_grid(rows = vars(method), cols = vars(TargetGene))
p <- p + theme_plot(rotate.x = 90) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + scale_color_manual(values = color_setup_simplified_norep) + guides(col = "none")
pdf(sprintf("%s/04_emVars/3_emVARs_agg_by_target/Effect_size__STAT5A_5B_3_genes__haplo_sum_LFDR_weighted_noHiC.pdf", FIGURE_DIR), height = 7, width = 7)
print(p)
dev.off()


locus_genes <- c("SLC6A20", "LZTFL1", "CCR5", "CCR1", "CCR2", "CCR5", "CXCR6", "LTF", "FYCO1", "CCR9", "XCL1", "LRRC2")
p <- ggplot(emVARs_agg_by_target[COND_ID %chin% names(color_setup_simplified_norep) & TargetGene %chin% locus_genes & method != "HiC", ]) +
  geom_hline(yintercept = 0, col = "lightgrey")
p <- p + geom_pointrange(aes(x = factor(COND_ID, names(color_setup_simplified_norep)), y = haplo_sum_weightedLFDR, ymin = haplo_sum_weightedLFDR - 1.96 * haplo_sum_weightedLFDR.se, ymax = haplo_sum_weightedLFDR + 1.96 * haplo_sum_weightedLFDR.se, col = COND_ID), size = .1)
p <- p + facet_grid(rows = vars(method), cols = vars(TargetGene))
p <- p + theme_plot(rotate.x = 90) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + scale_color_manual(values = color_setup_simplified_norep) + guides(col = "none")
pdf(sprintf("%s/04_emVars/3_emVARs_agg_by_target/Effect_size__chr3_46Mb__haplo_sum_LFDR_weighted_noHiC.pdf", FIGURE_DIR), height = 7, width = 7)
print(p)
dev.off()


locus_genes <- c("TNFAIP3", "MAP3K5", "IFNGR1", "IL20RA", "IL22RA", "OLIG3")
p <- ggplot(emVARs_agg_by_target[COND_ID %chin% names(color_setup_simplified_norep) & TargetGene %chin% locus_genes & method != "HiC", ]) +
  geom_hline(yintercept = 0, col = "lightgrey")
p <- p + geom_pointrange(aes(x = factor(COND_ID, names(color_setup_simplified_norep)), y = haplo_sum_weightedLFDR, ymin = haplo_sum_weightedLFDR - 1.96 * haplo_sum_weightedLFDR.se, ymax = haplo_sum_weightedLFDR + 1.96 * haplo_sum_weightedLFDR.se, col = COND_ID), size = .1)
p <- p + facet_grid(rows = vars(method), cols = vars(TargetGene))
p <- p + theme_plot(rotate.x = 90) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + scale_color_manual(values = color_setup_simplified_norep) + guides(col = "none")
pdf(sprintf("%s/04_emVars/3_emVARs_agg_by_target/Effect_size__IFNGR1_IL20RA__haplo_sum_LFDR_weighted_noHiC.pdf", FIGURE_DIR), height = 7, width = 7)
print(p)
dev.off()


locus_genes <- c("TANK", "DPP4", "IFIH1", "PSMD14")
p <- ggplot(emVARs_agg_by_target[COND_ID %chin% names(color_setup_simplified_norep) & TargetGene %chin% locus_genes & method != "HiC", ]) +
  geom_hline(yintercept = 0, col = "lightgrey")
p <- p + geom_pointrange(aes(x = factor(COND_ID, names(color_setup_simplified_norep)), y = haplo_sum_weightedLFDR, ymin = haplo_sum_weightedLFDR - 1.96 * haplo_sum_weightedLFDR.se, ymax = haplo_sum_weightedLFDR + 1.96 * haplo_sum_weightedLFDR.se, col = COND_ID), size = .1)
p <- p + facet_grid(rows = vars(method), cols = vars(TargetGene))
p <- p + theme_plot(rotate.x = 90) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + scale_color_manual(values = color_setup_simplified_norep) + guides(col = "none")
pdf(sprintf("%s/04_emVars/3_emVARs_agg_by_target/Effect_size__TANK_DPP4_IFIH1_PSMD14__haplo_sum_LFDR_weighted_noHiC.pdf", FIGURE_DIR), height = 7, width = 7)
print(p)
dev.off()

#################### haplo_sum_weightedLFDR_norm


locus_genes <- c("OAS1", "OAS2", "OAS3")
p <- ggplot(emVARs_agg_by_target[COND_ID %chin% names(color_setup_simplified_norep) & TargetGene %chin% locus_genes, ]) +
  geom_hline(yintercept = 0, col = "lightgrey")
p <- p + geom_pointrange(aes(x = factor(COND_ID, names(color_setup_simplified_norep)), y = haplo_sum_weightedLFDR_norm, ymin = haplo_sum_weightedLFDR_norm - 1.96 * haplo_sum_weightedLFDR_norm.se, ymax = haplo_sum_weightedLFDR_norm + 1.96 * haplo_sum_weightedLFDR_norm.se, col = COND_ID), size = .1)
p <- p + facet_grid(rows = vars(method), cols = vars(TargetGene))
p <- p + theme_plot(rotate.x = 90) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + scale_color_manual(values = color_setup_simplified_norep) + guides(col = "none")
pdf(sprintf("%s/04_emVars/3_emVARs_agg_by_target/Effect_size__OAS_genes__haplo_sum_LFDR_weighted_norm.pdf", FIGURE_DIR), height = 7, width = 7)
print(p)
dev.off()

pdf(sprintf("%s/04_emVars/3_emVARs_agg_by_target/Effect_size__TLR1_genes__haplo_sum_LFDR_weighted_norm2.pdf", FIGURE_DIR), height = 7, width = 7)
p <- p + facet_grid(rows = vars(method), cols = vars(TargetGene), scales = "free_y")
print(p)
dev.off()

locus_genes <- c("TLR1", "TLR6", "TLR10")
p <- ggplot(emVARs_agg_by_target[COND_ID %chin% names(color_setup_simplified_norep) & TargetGene %chin% locus_genes, ]) +
  geom_hline(yintercept = 0, col = "lightgrey")
p <- p + geom_pointrange(aes(x = factor(COND_ID, names(color_setup_simplified_norep)), y = haplo_sum_weightedLFDR_norm, ymin = haplo_sum_weightedLFDR_norm - 1.96 * haplo_sum_weightedLFDR_norm.se, ymax = haplo_sum_weightedLFDR_norm + 1.96 * haplo_sum_weightedLFDR_norm.se, col = COND_ID), size = .1)
p <- p + facet_grid(rows = vars(method), cols = vars(TargetGene))
p <- p + theme_plot(rotate.x = 90) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + scale_color_manual(values = color_setup_simplified_norep) + guides(col = "none")
pdf(sprintf("%s/04_emVars/3_emVARs_agg_by_target/Effect_size__TLR1_genes__haplo_sum_LFDR_weighted_norm.pdf", FIGURE_DIR), height = 7, width = 7)
print(p)
dev.off()
pdf(sprintf("%s/04_emVars/3_emVARs_agg_by_target/Effect_size__TLR1_genes__haplo_sum_LFDR_weighted_norm2.pdf", FIGURE_DIR), height = 7, width = 7)
p <- p + facet_grid(rows = vars(method), cols = vars(TargetGene), scales = "free_y")
print(p)
dev.off()

locus_genes <- c("PAN2", "STAT2", "IL23A", "APOF")
p <- ggplot(emVARs_agg_by_target[COND_ID %chin% names(color_setup_simplified_norep) & TargetGene %chin% locus_genes, ]) +
  geom_hline(yintercept = 0, col = "lightgrey")
p <- p + geom_pointrange(aes(x = factor(COND_ID, names(color_setup_simplified_norep)), y = haplo_sum_weightedLFDR_norm, ymin = haplo_sum_weightedLFDR_norm - 1.96 * haplo_sum_weightedLFDR_norm.se, ymax = haplo_sum_weightedLFDR_norm + 1.96 * haplo_sum_weightedLFDR_norm.se, col = COND_ID), size = .1)
p <- p + facet_grid(rows = vars(method), cols = vars(TargetGene))
p <- p + theme_plot(rotate.x = 90) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + scale_color_manual(values = color_setup_simplified_norep) + guides(col = "none")
pdf(sprintf("%s/04_emVars/3_emVARs_agg_by_target/Effect_size__STAT2_genes__haplo_sum_LFDR_weighted_norm.pdf", FIGURE_DIR), height = 7, width = 7)
print(p)
dev.off()

pdf(sprintf("%s/04_emVars/3_emVARs_agg_by_target/Effect_size__STAT2_genes__haplo_sum_LFDR_weighted_norm2.pdf", FIGURE_DIR), height = 7, width = 7)
p <- p + facet_grid(rows = vars(method), cols = vars(TargetGene), scales = "free_y")
print(p)
dev.off()

locus_genes <- c("STAT3", "STAT5A", "STAT5B")
p <- ggplot(emVARs_agg_by_target[COND_ID %chin% names(color_setup_simplified_norep) & TargetGene %chin% locus_genes, ]) +
  geom_hline(yintercept = 0, col = "lightgrey")
p <- p + geom_pointrange(aes(x = factor(COND_ID, names(color_setup_simplified_norep)), y = haplo_sum_weightedLFDR_norm, ymin = haplo_sum_weightedLFDR_norm - 1.96 * haplo_sum_weightedLFDR_norm.se, ymax = haplo_sum_weightedLFDR_norm + 1.96 * haplo_sum_weightedLFDR_norm.se, col = COND_ID), size = .1)
p <- p + facet_grid(rows = vars(method), cols = vars(TargetGene))
p <- p + theme_plot(rotate.x = 90) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + scale_color_manual(values = color_setup_simplified_norep) + guides(col = "none")
pdf(sprintf("%s/04_emVars/3_emVARs_agg_by_target/Effect_size__STAT5A_5B_3_genes__haplo_sum_LFDR_weighted_norm.pdf", FIGURE_DIR), height = 7, width = 7)
print(p)
dev.off()
pdf(sprintf("%s/04_emVars/3_emVARs_agg_by_target/Effect_size__STAT5A_5B_3_genes__haplo_sum_LFDR_weighted_norm2.pdf", FIGURE_DIR), height = 7, width = 7)
p <- p + facet_grid(rows = vars(method), cols = vars(TargetGene), scales = "free_y")
print(p)
dev.off()


locus_genes <- c("SLC6A20", "LZTFL1", "CCR5", "CCR1", "CCR2", "CCR5", "CXCR6", "LTF", "FYCO1", "CCR9", "XCL1", "LRRC2")
p <- ggplot(emVARs_agg_by_target[COND_ID %chin% names(color_setup_simplified_norep) & TargetGene %chin% locus_genes, ]) +
  geom_hline(yintercept = 0, col = "lightgrey")
p <- p + geom_pointrange(aes(x = factor(COND_ID, names(color_setup_simplified_norep)), y = haplo_sum_weightedLFDR_norm, ymin = haplo_sum_weightedLFDR_norm - 1.96 * haplo_sum_weightedLFDR_norm.se, ymax = haplo_sum_weightedLFDR_norm + 1.96 * haplo_sum_weightedLFDR_norm.se, col = COND_ID), size = .1)
p <- p + facet_grid(rows = vars(method), cols = vars(TargetGene))
p <- p + theme_plot(rotate.x = 90) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + scale_color_manual(values = color_setup_simplified_norep) + guides(col = "none")
pdf(sprintf("%s/04_emVars/3_emVARs_agg_by_target/Effect_size__chr3_46Mb__haplo_sum_LFDR_weighted_norm.pdf", FIGURE_DIR), height = 7, width = 7)
print(p)
dev.off()

pdf(sprintf("%s/04_emVars/3_emVARs_agg_by_target/Effect_size__chr3_46Mb__haplo_sum_LFDR_weighted_norm2.pdf", FIGURE_DIR), height = 7, width = 7)
p <- p + facet_grid(rows = vars(method), cols = vars(TargetGene), scales = "free_y")
print(p)
dev.off()

locus_genes <- c("TNFAIP3", "MAP3K5", "IFNGR1", "IL20RA", "IL22RA", "OLIG3")
p <- ggplot(emVARs_agg_by_target[COND_ID %chin% names(color_setup_simplified_norep) & TargetGene %chin% locus_genes, ]) +
  geom_hline(yintercept = 0, col = "lightgrey")
p <- p + geom_pointrange(aes(x = factor(COND_ID, names(color_setup_simplified_norep)), y = haplo_sum_weightedLFDR_norm, ymin = haplo_sum_weightedLFDR_norm - 1.96 * haplo_sum_weightedLFDR_norm.se, ymax = haplo_sum_weightedLFDR_norm + 1.96 * haplo_sum_weightedLFDR_norm.se, col = COND_ID), size = .1)
p <- p + facet_grid(rows = vars(method), cols = vars(TargetGene))
p <- p + theme_plot(rotate.x = 90) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + scale_color_manual(values = color_setup_simplified_norep) + guides(col = "none")
pdf(sprintf("%s/04_emVars/3_emVARs_agg_by_target/Effect_size__IFNGR1_IL20RA__haplo_sum_LFDR_weighted_norm.pdf", FIGURE_DIR), height = 7, width = 7)
print(p)
dev.off()

pdf(sprintf("%s/04_emVars/3_emVARs_agg_by_target/Effect_size__IFNGR1_IL20RA__haplo_sum_LFDR_weighted_norm2.pdf", FIGURE_DIR), height = 7, width = 7)
p <- p + facet_grid(rows = vars(method), cols = vars(TargetGene), scales = "free_y")
print(p)
dev.off()

locus_genes <- c("TANK", "DPP4", "IFIH1", "PSMD14")
p <- ggplot(emVARs_agg_by_target[COND_ID %chin% names(color_setup_simplified_norep) & TargetGene %chin% locus_genes, ]) +
  geom_hline(yintercept = 0, col = "lightgrey")
p <- p + geom_pointrange(aes(x = factor(COND_ID, names(color_setup_simplified_norep)), y = haplo_sum_weightedLFDR_norm, ymin = haplo_sum_weightedLFDR_norm - 1.96 * haplo_sum_weightedLFDR_norm.se, ymax = haplo_sum_weightedLFDR_norm + 1.96 * haplo_sum_weightedLFDR_norm.se, col = COND_ID), size = .1)
p <- p + facet_grid(rows = vars(method), cols = vars(TargetGene))
p <- p + theme_plot(rotate.x = 90) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + scale_color_manual(values = color_setup_simplified_norep) + guides(col = "none")
pdf(sprintf("%s/04_emVars/3_emVARs_agg_by_target/Effect_size__TANK_DPP4_IFIH1_PSMD14__haplo_sum_LFDR_weighted_norm.pdf", FIGURE_DIR), height = 7, width = 7)
print(p)
dev.off()

pdf(sprintf("%s/04_emVars/3_emVARs_agg_by_target/Effect_size__TANK_DPP4_IFIH1_PSMD14__haplo_sum_LFDR_weighted_norm2.pdf", FIGURE_DIR), height = 7, width = 7)
p <- p + facet_grid(rows = vars(method), cols = vars(TargetGene), scales = "free_y")
print(p)
dev.off()











## ARCHAIC ALLELES ARE ASSOCIATED WITH LOWER ENHANCERS ACTIVITY AND STRONGER REPRESSOR ACTIVITY
emVARs_annot_active_ok[FDR < .01 & abs(log2FC_archaic_vs_modern) > .5, .(.N, mean(ifelse(grepl("enhancer", oligo_class_loose), log2FC_archaic_vs_modern > 0, log2FC_archaic_vs_modern < 0))), by = .(grepl("enhancer", oligo_class_loose), ANALYSIS_NAME)][grepl("all", ANALYSIS_NAME)][grepl == FALSE, ]

emVARs_annot_active_ok[FDR < .01 & abs(log2FC_archaic_vs_modern) > .5, .(.N, mean(ifelse(grepl("enhancer", oligo_class_loose), log2FC_archaic_vs_modern > 0, log2FC_archaic_vs_modern < 0))), by = .(grepl("enhancer", oligo_class_loose), ANALYSIS_NAME)][grepl("all", ANALYSIS_NAME)][grepl == FALSE, ]

emVARs_annot_active_ok[FDR < .01 & abs(log2FC_archaic_vs_modern) > .5, .(.N, mean(ifelse(grepl("enhancer", oligo_class_loose), log2FC_archaic_vs_modern > 0, log2FC_archaic_vs_modern < 0))), keyby = .(grepl("enhancer", oligo_class_loose), POP_adaptive)]

emVARs_annot_active_ok[grepl("all", ANALYSIS_NAME) & FDR < .01 & abs(log2FC_archaic_vs_modern) > .5 & !ctrl & is_introgressed, .(.N, mean(ifelse(grepl("enhancer", oligo_class_loose), log2FC_archaic_vs_modern > 0, log2FC_archaic_vs_modern < 0))), keyby = .(Introgression_scenario)]


# ARCHAIC DERIVED ALLELES ARE ASSOCIATED WITH LOWER ENHANCERS ACTIVITY
FigData <- emVARs_annot_active_ok[grepl("all", ANALYSIS_NAME) & FDR < .05 & abs(log2FC_archaic_vs_modern) > .2 & !ctrl & is_introgressed, ]
FigData <- FigData[Introgression_scenario_v2 != "WTF", .(
  N_emVars = .N,
  N_higher_activy_archaic = sum(ifelse(grepl("enhancer", oligo_class_loose), log2FC_archaic_vs_modern > 0, log2FC_archaic_vs_modern < 0)),
  N_higher_activy_modern = sum(ifelse(grepl("enhancer", oligo_class_loose), log2FC_archaic_vs_modern < 0, log2FC_archaic_vs_modern > 0))
), keyby = .(crs_type = ifelse(grepl("enhancer", oligo_class_loose), "enhancer", "silencer"), Introgression_scenario_v2)]
FigData[, pct_higher_activy_archaic_inf := binom.test(N_higher_activy_archaic, N_emVars)$conf.int[1], keyby = .(crs_type, Introgression_scenario_v2)]
FigData[, pct_higher_activy_archaic_sup := binom.test(N_higher_activy_archaic, N_emVars)$conf.int[2], keyby = .(crs_type, Introgression_scenario_v2)]
FigData[, binom_P := binom.test(N_higher_activy_archaic, N_emVars)$p.value, keyby = .(crs_type, Introgression_scenario_v2)]

p <- ggplot(FigData) +
  geom_hline(yintercept = 50, col = "lightgrey", linetype = "dashed")
p <- p + geom_pointrange(aes(x = Introgression_scenario_v2, y = N_higher_activy_archaic / N_emVars * 100, ymin = pct_higher_activy_archaic_inf * 100, ymax = pct_higher_activy_archaic_sup * 100, col = Introgression_scenario_v2), size = .5)
p <- p + facet_grid(cols = vars(crs_type)) + ylab("Percentage of emVars that have\nhigher activity in archaics\n compared to modern humans")
p <- p + theme_plot(rotate.x = 90) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + guides(col = "none") + xlab("")

pdf(sprintf("%s/Agg_emVars_by_type_introgression_v2.pdf", FIGURE_DIR), height = 3.5, width = 3)
print(p)
dev.off()

cat("All done !\n")
q("no")
#############################################################################
#############################################################################
#############################################################################


# A note on the analysis startegy for power related analyses :
# For each CRS, we are going to :
# 1. report the number of barcode in the library and
# 2. estimate the power to detect activity at the reported level of activity
#    2a. this can be computed for the Winner's curse adjusted beta values
#        (use FIQT: F DR I nverse Q uantile T ransformation).
#         doi.org/10.1093/bioinformatics/btw303
#    2b. perhaps use local FDR from fdrtool rather than FDR ?

# 3. test whether probability of having an emVar varies across tissue/conditionb, adjiusting for difference in power.
