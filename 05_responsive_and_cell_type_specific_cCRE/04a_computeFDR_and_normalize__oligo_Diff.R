# 04a_computeFDR_and_normalize_oligo_Diff_activity.R

# running: sbatch -p geh,common --mem=30G  00_Rscript.sh MPRA_count_exp6_analysisZ/03b_aggMPRAnalyse_oligo_activity.R

MPRA_DIR <- "/pasteur/helix/projects/evo_immuno_pop/MPRA"
ANALYSIS_DIR <- "MPRA_count_exp6_analysis3"
SCRIPT_DIR <- sprintf("%s/scripts/%s/", MPRA_DIR, ANALYSIS_DIR)

source(sprintf("%s/scripts/%s/00_starter.R", MPRA_DIR, ANALYSIS_DIR))
source(sprintf("%s/scripts/%s/00_load_associations.R", MPRA_DIR, ANALYSIS_DIR))
# source(sprintf("%s/scripts/misc_plots.R", MPRA_DIR))

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

FIGURE_DIR <- sprintf("%s/figures/%s/%s", MPRA_DIR, ANALYSIS_DIR, RUN_ID)
dir.create(FIGURE_DIR, recursive = TRUE)

tic("loading oligo & SNP annotations")
# load annotation of oligos (beforz: CRS)
oligo_source <- fread(sprintf("%s/data/%s/00_oligo_annot_v2.txt", MPRA_DIR, ANALYSIS_DIR))

# load annotations of SNPs
SNP_annot <- fread(sprintf("%s/data/%s/00_SNP_annot_v4.txt", MPRA_DIR, ANALYSIS_DIR))
toc()

# source(sprintf("%s/scripts/%s/03_00_parameter_activity.R", MPRA_DIR, ANALYSIS_DIR))
source(sprintf("%s/scripts/%s/04_00_parameter_diff_activity.R", MPRA_DIR, ANALYSIS_DIR))

IN_DIR <- sprintf("%s/data/%s/03a_aggMPRA_analyse/%s/", MPRA_DIR, ANALYSIS_DIR, RUN_ID)
DIFF_DIR <- sprintf("%s/data/%s/03b_computeFDR/%s/Diff/%s/%s", MPRA_DIR, ANALYSIS_DIR, RUN_ID, CRITERIA_DIFF, CRITERIA_ACTIVE_OUT)
ACTIVE_DIR <- sprintf("%s/data/%s/03b_computeFDR/%s/Activity/%s/", MPRA_DIR, ANALYSIS_DIR, RUN_ID, CRITERIA_ACTIVE)
dir.create(DIFF_DIR, recursive = TRUE)

# load oligo annotations
source(sprintf("%s/scripts/%s/02z__define_includedCRS.R", MPRA_DIR, ANALYSIS_DIR))

####################################################################
################# compute FDR, CRS difference ######################
####################################################################

comparison_labels <- fread(sprintf("%s/04a_comparison_labels.txt", SCRIPT_DIR))

oligo_activity_Diff <- fread(file = sprintf("%s/all_oligo_diff_results.tsv.gz", IN_DIR))
oligo_activity_Diff <- oligo_activity_Diff[oligo %in% tested_and_ctrl_oligos_final_annot$oligo, ]
oligo_activity_Diff[, log2FC_2vs1 := logFC_2vs1 / log(2)]
oligo_activity_Diff[, log2FC.se := logFC.se / log(2)]
oligo_activity_Diff[, .(.N, sum(df.test > 0)), keyby = .(power, boot, ANALYSIS_NAME, Z_th, perm)]

oligo_activity_Diff_annot <- merge(oligo_activity_Diff[df.test > 0, ], oligo_source[excluded == FALSE, ], by = "oligo", allow.cartesian = TRUE)
oligo_activity_Diff_annot <- merge(oligo_activity_Diff_annot, comparison_labels, by.x = "ANALYSIS_NAME", by.y = "analysis_name")
# setnames(oligo_activity_Diff,'CRS','oligo')
if (FILTER_ACTIVE) {
  # subset to oligos that are active in at least one condition
  oligo_activity_obs <- fread(sprintf("%s/all_oligos_annotated__%s.tsv.gz", ACTIVE_DIR, CRITERIA_ACTIVE))
  #  oligo_activity_obs_ctc <- oligo_activity_obs[power == 0 & boot == 0, ]
  all_active_oligos <- oligo_activity_obs[oligo_active_CRITERIA == TRUE, .(oligo, crsID, ANALYSIS_NAME, ANALYSIS_SUBTYPE)]
  oligo_activity_Diff_annot[, is_active_oligo_group1 := paste(crsID, group1_labels) %chin% all_active_oligos[, paste(crsID, ANALYSIS_NAME)]]
  oligo_activity_Diff_annot[, is_active_oligo_group2 := paste(crsID, group2_labels) %chin% all_active_oligos[, paste(crsID, ANALYSIS_NAME)]]
  oligo_activity_Diff_annot[, is_tested_oligo := is_active_oligo_group1 | is_active_oligo_group2]
} else {
  oligo_activity_Diff_annot[, is_tested_oligo := TRUE]
}
cat(
  "\n", oligo_activity_Diff_annot[, length(unique(oligo))], "oligo tested from",
  oligo_activity_Diff_annot[, length(unique(crsID))], "CRE, across ",
  oligo_activity_Diff_annot[, length(unique(ANALYSIS_NAME))], "comparison, for a total of",
  oligo_activity_Diff_annot[, .N], "tests (",
  oligo_activity_Diff_annot[, length(unique(paste(crsID, ANALYSIS_NAME)))], "cre x comparison)\n"
)

# oligo_activity_Diff[, log2FC_der_vs_anc := ifelse(allele.1 == ANCESTRAL, logFC_2vs1 / log(2), -logFC_2vs1 / log(2))]
# oligo_activity_Diff[, log2FC_archaic_vs_modern := ifelse(!grepl("ancestral reintrogressed", SNP_type), log2FC_der_vs_anc, -log2FC_der_vs_anc)]

oligo_activity_Diff_annot <- oligo_activity_Diff_annot[order(pval.LRT), .SD, by = .(Z_th, power, boot, ANALYSIS_NAME)]

# estimate global FDR based on permutations
oligo_activity_Diff_annot[, N_FP := (.1 + cumsum(perm != "perm_0_0_0" & df.test > 0 & is_tested_oligo)), by = .(Z_th, power, boot, ANALYSIS_NAME)]
oligo_activity_Diff_annot[, N_POS := (.1 + cumsum(perm == "perm_0_0_0" & df.test > 0 & is_tested_oligo)), by = .(Z_th, power, boot, ANALYSIS_NAME)]
oligo_activity_Diff_annot[, pval_emp := N_FP / (.1 + sum(perm != "perm_0_0_0" & df.test > 0 & is_tested_oligo)), by = .(Z_th, power, boot, ANALYSIS_NAME)]
oligo_activity_Diff_annot[is_tested_oligo == TRUE, FDR := pmin(1, rev(cummin(rev(N_FP / N_POS)))), by = .(Z_th, power, boot, ANALYSIS_NAME)]
oligo_activity_Diff_annot[is_tested_oligo == FALSE, FDR := 1]

# estimate local FDR based on permutations
oligo_activity_Diff_annot[df.test > 0 & is_tested_oligo == TRUE, bin_logP := cut(-log10(pval.LRT), breaks = c(seq(0, 25, by = 0.25), Inf), include.lowest = TRUE)]
lfdr_estim <- oligo_activity_Diff_annot[df.test > 0, .(n0 = sum(perm != "perm_0_0_0" & is_tested_oligo), n1 = sum(perm == "perm_0_0_0" & is_tested_oligo)), keyby = .(Z_th, power, boot, ANALYSIS_NAME, bin_logP)]
lfdr_estim[, lfdr := pmin(1, n0 / n1), keyby = .(Z_th, power, boot, ANALYSIS_NAME, bin_logP)]
lfdr_estim[, monotonic_bins := cumsum(!duplicated(rev(cummax(rev(lfdr))))), keyby = .(Z_th, power, boot, ANALYSIS_NAME)]
lfdr_estim[, lfdr := pmin(1, sum(n0) / sum(n1)), keyby = .(Z_th, power, boot, ANALYSIS_NAME, monotonic_bins)]

oligo_activity_Diff_annot <- merge(oligo_activity_Diff_annot, lfdr_estim[, .(Z_th, power, boot, ANALYSIS_NAME, bin_logP, lfdr)], by = c("Z_th", "power", "boot", "ANALYSIS_NAME", "bin_logP"), all.x = TRUE)
oligo_activity_Diff_annot[is_tested_oligo == FALSE, lfdr := 1]

# define differential oligos:
oligo_activity_Diff_annot[, oligo_diff_CRITERIA := lfdr < LFDR_TH_DIFF & FDR < FDR_TH_DIFF & abs(log2FC_2vs1) > LOG2FC_TH_DIFF]
oligo_activity_Diff_annot[, cre_diff_CRITERIA := any(oligo_diff_CRITERIA), by = .(Z_th, power, boot, ANALYSIS_NAME, crsID)]

# emVARs_Diff_annot[FDR < .05 & power == 0 & perm == "perm_0_0_0" & grepl("celltype", ANALYSIS_SUBTYPE), .N, by = ANALYSIS_NAME]
oligo_activity_Diff_obs <- oligo_activity_Diff_annot[power == 0 & perm == "perm_0_0_0", ]
oligo_activity_Diff_perm <- oligo_activity_Diff_annot[power == 0 & perm == "perm_0_1_0", ]

fwrite(oligo_activity_Diff_obs, file = sprintf("%s/all_oligos_diff_annotated__%s.tsv.gz", DIFF_DIR, CRITERIA_DIFF), sep = "\t")
fwrite(oligo_activity_Diff_perm, file = sprintf("%s/all_oligos_diff_annotated_perm__%s.tsv.gz", DIFF_DIR, CRITERIA_DIFF), sep = "\t")
oligo_activity_Diff_obs <- fread(sprintf("%s/all_oligos_diff_annotated__%s.tsv.gz", DIFF_DIR, CRITERIA_DIFF))
oligo_response_annot_obs <- oligo_activity_Diff_obs[grepl("response", ANALYSIS_NAME), ]
oligo_response_annot_obs_archaic <- oligo_response_annot_obs[oligo %chin% tested_and_ctrl_oligos_final_annot[type == "tested", oligo], ]

oligo_celltype_comp_NS_annot_obs <- oligo_activity_Diff_obs[grepl("celltype_comp", ANALYSIS_SUBTYPE) & grepl("_NS$", ANALYSIS_NAME), ]
oligo_celltype_comp_NS_annot_obs_archaic <- oligo_celltype_comp_NS_annot_obs[oligo %chin% tested_and_ctrl_oligos_final_annot[type == "tested", oligo], ]

oligo_celltype_comp_all_annot_obs <- oligo_activity_Diff_obs[grepl("celltype_comp", ANALYSIS_SUBTYPE) & grepl("_all$", ANALYSIS_NAME), ]
oligo_celltype_comp_all_annot_obs_archaic <- oligo_celltype_comp_NS_annot_obs[oligo %chin% tested_and_ctrl_oligos_final_annot[type == "tested", oligo], ]

# oligo_activity_Diff_annot[, .N, by = .(ANALYSIS_NAME, FDR < .05, FDR_analytic < 0.05)]
# oligo_activity_Diff_annot[, .(.N, N_perm = sum(FDR < .05), N_up = sum(FDR < .05 & log2FC_2vs1 > 0), N_down = sum(FDR < .05 & log2FC_2vs1 < 0)), keyby = .(power, ANALYSIS_NAME)]
# oligo_activity_Diff_annot[, .(.N, N_perm = sum(FDR < .05 & abs(log2FC_2vs1) > 0.2), N_up = sum(FDR < .05 & log2FC_2vs1 > 0.2), N_down = sum(FDR < .05 & log2FC_2vs1 < -0.2)), keyby = .(power, ANALYSIS_NAME)]
# oligo_activity_Diff_annot[, .(.N, N_perm = sum(lfdr < .2 & abs(log2FC_2vs1) > 0.2), N_up = sum(lfdr < .2 & log2FC_2vs1 > 0.2), N_down = sum(lfdr < .2 & log2FC_2vs1 < -0.2)), keyby = .(power, ANALYSIS_NAME)]

# response
oligo_response_annot_obs[, .(
  N_oligo_tested = .N,
  N_oligo_diff = sum(oligo_diff_CRITERIA),
  N_oligo_up = sum(oligo_diff_CRITERIA & log2FC_2vs1 > 0),
  N_oligo_down = sum(oligo_diff_CRITERIA & log2FC_2vs1 < 0),
  N_cre_tested = length(unique(crsID)),
  N_cre_diff = length(unique(crsID[oligo_diff_CRITERIA])),
  N_cre_up = length(unique(crsID[oligo_diff_CRITERIA & log2FC_2vs1 > 0])),
  N_cre_down = length(unique(crsID[oligo_diff_CRITERIA & log2FC_2vs1 < 0]))
), keyby = .(power, ANALYSIS_NAME)]

oligo_response_annot_obs_archaic[, .(
  N_oligo_tested = .N,
  N_oligo_diff = sum(oligo_diff_CRITERIA),
  N_oligo_up = sum(oligo_diff_CRITERIA & log2FC_2vs1 > 0),
  N_oligo_down = sum(oligo_diff_CRITERIA & log2FC_2vs1 < 0),
  N_cre_tested = length(unique(crsID)),
  N_cre_diff = length(unique(crsID[oligo_diff_CRITERIA])),
  N_cre_up = length(unique(crsID[oligo_diff_CRITERIA & log2FC_2vs1 > 0])),
  N_cre_down = length(unique(crsID[oligo_diff_CRITERIA & log2FC_2vs1 < 0]))
), keyby = .(power, ANALYSIS_NAME)]

# response
oligo_celltype_comp_NS_annot_obs[, .(
  N_oligo_tested = .N,
  N_oligo_diff = sum(oligo_diff_CRITERIA),
  N_oligo_up = sum(oligo_diff_CRITERIA & log2FC_2vs1 > 0),
  N_oligo_down = sum(oligo_diff_CRITERIA & log2FC_2vs1 < 0),
  N_cre_tested = length(unique(crsID)),
  N_cre_diff = length(unique(crsID[oligo_diff_CRITERIA])),
  N_cre_up = length(unique(crsID[oligo_diff_CRITERIA & log2FC_2vs1 > 0])),
  N_cre_down = length(unique(crsID[oligo_diff_CRITERIA & log2FC_2vs1 < 0]))
), keyby = .(power, ANALYSIS_NAME)]

oligo_celltype_comp_NS_annot_obs_archaic[, .(
  N_oligo_tested = .N,
  N_oligo_diff = sum(oligo_diff_CRITERIA),
  N_oligo_up = sum(oligo_diff_CRITERIA & log2FC_2vs1 > 0),
  N_oligo_down = sum(oligo_diff_CRITERIA & log2FC_2vs1 < 0),
  N_cre_tested = length(unique(crsID)),
  N_cre_diff = length(unique(crsID[oligo_diff_CRITERIA])),
  N_cre_up = length(unique(crsID[oligo_diff_CRITERIA & log2FC_2vs1 > 0])),
  N_cre_down = length(unique(crsID[oligo_diff_CRITERIA & log2FC_2vs1 < 0]))
), keyby = .(power, ANALYSIS_NAME)]


### reading outputs:

# MPRA_DIR <- "/pasteur/helix/projects/evo_immuno_pop/MPRA"
# ANALYSIS_DIR <- "MPRA_count_exp6_analysis2"
# RUN_ID <- "RUN2_Z2_nBC10"
# OUT_DIR <- sprintf("%s/data/%s/03b_computeFDR/%s/", MPRA_DIR, ANALYSIS_DIR, RUN_ID)
# oligo_Diff_activity_obs <- fread(sprintf("%s/all_oligos_diff_annotated.tsv.gz", OUT_DIR))


# oligo_activity_Diff_annot[, .(.N, N_perm = sum(FDR < .05)), keyby = .(power, ANALYSIS_NAME)]

# fwrite( oligo_activity_DT, file = sprintf("%s/data/%s/MPRA_analyses/%s/allConditions_alphas_and_pvalues_v%s.1_max%sBC_annot.txt.gz", MPRA_DIR, ANALYSIS_DIR, ' oligo_activity', VERSION_OUT, sub_sample), sep = "\t")

# TODO : rely on permutation-based FDR when available
# trend : unclear, but more significant stuff in permutation-based FDR for cell type comparisons and HepG2-IFN & K562-NS replicates, and less for response_A549_SARS response

# CRS_DiffActivity_DT[analysis_name == "HepG2_IFNA2b_replicates" & FDR < .001 & FDR_analytic < 0.0001 & (nBC_g1 + nBC_g2) > 100, ][order(-activity_group1 - activity_group2)][1:3, ]
# x <- CRS_DiffActivity_DT[analysis_name == "HepG2_IFNA2b_replicates", ]

# CRS_DiffActivity_DT[, ROADMAP_current_celline := case_when(
#   grepl("response_A549", analysis_name) ~ E114,
#   grepl("A549_.*_replicates", analysis_name) ~ E114,
#   grepl("response_HepG2", analysis_name) ~ E118,
#   grepl("HepG2_.*_replicates", analysis_name) ~ E118,
#   grepl("response_K562", analysis_name) ~ E123,
#   grepl("K562_.*_replicates", analysis_name) ~ E123,
#   TRUE ~ NA
# )]
# CRS_DiffActivity_DT[, ROADMAP_ctrl_brain := E071]
# CRS_DiffActivity_DT[, celline := case_when(
#   grepl("response_A549", analysis_name) ~ "A549",
#   grepl("A549_.*_replicates", analysis_name) ~ "A549",
#   grepl("response_HepG2", analysis_name) ~ "HepG2",
#   grepl("HepG2_.*_replicates", analysis_name) ~ "HepG2",
#   grepl("response_K562", analysis_name) ~ "K562",
#   grepl("K562_.*_replicates", analysis_name) ~ "K562",
#   TRUE ~ NA
# )]

# fwrite(CRS_DiffActivity_DT, file = sprintf("%s/data/%s/MPRA_analyses/%s/allInteractions_alphas_and_pvalues_v%s.1_max%sBC_annot.txt.gz", MPRA_DIR, ANALYSIS_DIR, "CRS_activity", VERSION_IN, sub_sample), sep = "\t")
# CRS_DiffActivity_DT[,mean(FDR<.05),by=ROADMAP_current_celline]
cat("All done!\n")
q("no")
