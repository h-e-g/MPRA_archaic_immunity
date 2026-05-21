# running: sbatch -p geh,common --mem=30G  00_Rscript.sh MPRA_count_exp6_analysisZ/03b_aggMPRAnalyse_oligo_activity.R

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
dir.create(ACTIVE_DIR, recursive = TRUE)
FIGURE_DIR <- sprintf("%s/figures/%s/%s", MPRA_DIR, ANALYSIS_DIR, RUN_ID)

source(sprintf("%s/scripts/%s/02z__define_includedCRS.R", MPRA_DIR, ANALYSIS_DIR))

####################################################################
################# compute FDR, oligo activity ######################
####################################################################
cat('compute FDR, oligo activity ')

oligo_activity <- fread(file = sprintf("%s/all_oligo_results.tsv.gz", IN_DIR), sep = "\t")
oligo_activity <- oligo_activity[oligo %in% tested_and_ctrl_oligos_final$oligo, ]
# setnames(oligo_activity,'CRS','oligo')
# setnames(oligo_activity,'posID','crsID')

############# compute activity pvalue #############
oligo_activity_annot <- merge(oligo_activity[, -"GC"], oligo_source[excluded == FALSE, -c("sequence")], by = "oligo", allow.cartesian = TRUE)

if (GC_CORRECT == "none") {
  # do not adjust on GC
  oligo_activity_annot[, alpha.GC := alpha, by = .(ANALYSIS_NAME, perm, power, boot, Z_th)]
}
if (GC_CORRECT == "observed") {
  # adjust on GC based on all oligos
  oligo_activity_annot[, alpha.GC := exp(lm(log(alpha) ~ GC)$residuals), by = .(ANALYSIS_NAME, perm, power, boot, Z_th)]
}
if (GC_CORRECT == "scrambled") {
  # adjust on GC based on scrambled trend
  oligo_activity_annot[, alpha.GC := exp(log(alpha) - scale(GC, T, F) * lm(log(alpha)[type == "scrambled"] ~ scale(GC[type == "scrambled"], T, F))$coef[2]), by = .(ANALYSIS_NAME, perm, power, boot, Z_th)]
}
# NOT NEEDED: oligo_activity_annot[,summary(lm(I(log(alpha.GC)-median(log(alpha.GC)))**2~nBC_g1+I(1/nBC_g1)+I(nBC_g1**2))),by=.(ANALYSIS_NAME,perm,power,boot)]


# 1# compare_GC effects
# compare_GC=oligo_activity_annot[, .(scrambled_estimated_GC_effect = lm(log(alpha)[type == "scrambled"] ~ scale(GC[type == "scrambled"], T, F))$coef[2],
#                         all_estimated_GC_effect = lm(log(alpha) ~ scale(GC, T, F))$coef[2]), by = .(ANALYSIS_NAME, perm, power, boot, Z_th)][perm=='perm_0_0_0' & power==0]
# compare_GC[grepl('all',ANALYSIS_NAME),]

# normalize activity to have same mean absolute deviation on observed data

if (NORMALIZATION == "none") {
  oligo_activity_annot[, alpha.se_norm := alpha.se, by = .(ANALYSIS_NAME, perm, power, boot, Z_th)]
  oligo_activity_annot[, alpha.GC_norm := alpha.GC, by = .(ANALYSIS_NAME, perm, power, boot, Z_th)]
}


if (NORMALIZATION == "center") {
  oligo_activity_annot[, alpha.se_norm := alpha.se, by = .(ANALYSIS_NAME, perm, power, boot, Z_th)]
  oligo_activity_annot[, alpha.GC_norm := exp((log(alpha.GC) - median(log(alpha.GC)))), by = .(ANALYSIS_NAME, perm, power, boot, Z_th)]
}

if (NORMALIZATION == "scale_all") {
  oligo_activity_annot[, VAR_analysis := mad(log(alpha.GC)[perm == "perm_0_0_0"]), by = .(ANALYSIS_NAME, power, boot, Z_th)] # variance of activity related to the number of barcodes per oligo
  oligo_activity_annot[, VAR_inflation_factor := VAR_analysis / mean(VAR_analysis)]
  oligo_activity_annot[, alpha.se_norm := alpha.se / VAR_inflation_factor, by = .(ANALYSIS_NAME, perm, power, boot, Z_th)]
  oligo_activity_annot[, alpha.GC_norm := exp((log(alpha.GC) - median(log(alpha.GC))) / VAR_inflation_factor), by = .(ANALYSIS_NAME, perm, power, boot, Z_th)]
}

if (NORMALIZATION == "scale_perm") {
  oligo_activity_annot[, VAR_analysis := mad(log(alpha.GC)[perm == "perm_1_0_0"]), by = .(ANALYSIS_NAME, power, boot, Z_th)] # variance of activity related to the number of barcodes per oligo
  oligo_activity_annot[, VAR_inflation_factor := VAR_analysis / mean(VAR_analysis)]
  oligo_activity_annot[, alpha.se_norm := alpha.se / VAR_inflation_factor, by = .(ANALYSIS_NAME, perm, power, boot, Z_th)]
  oligo_activity_annot[, alpha.GC_norm := exp((log(alpha.GC) - median(log(alpha.GC))) / VAR_inflation_factor), by = .(ANALYSIS_NAME, perm, power, boot, Z_th)]
}

if (NORMALIZATION == "scale_nBC_perm") {
  # estimate relation ship between variance an Number of barcodes
  Var_estim <- oligo_activity_annot[perm == "perm_1_0_0" & boot == 0 & power == 0, ]
  Var_estim[, dev := (log(alpha.GC) - median(log(alpha.GC)))**2, by = ANALYSIS_NAME]
  Var_param <- Var_estim[, as.list(lm(dev ~ I(1 / (nBC_g1)))$coeff), by = .(ANALYSIS_NAME)]
  oligo_activity_annot <- merge(oligo_activity_annot, Var_param, by = "ANALYSIS_NAME", all.x = TRUE)
  setnames(oligo_activity_annot, colnames(Var_param)[-1], c("VAR_cste", "VAR_nBC"))
  oligo_activity_annot[, VAR_estimated_nBC := sqrt(VAR_cste + VAR_nBC / nBC_g1), by = .(ANALYSIS_NAME)]
  oligo_activity_annot[, VAR_inflation_nBC := VAR_estimated_nBC / mean(VAR_estimated_nBC)]
  # adjust log activity to have same mean absolute deviation regardless of barcode number.
  oligo_activity_annot[, alpha.GC_norm := exp((log(alpha.GC) - median(log(alpha.GC))) / VAR_inflation_nBC), by = .(ANALYSIS_NAME, perm, power, boot, Z_th)]
  oligo_activity_annot[, alpha.se_norm := alpha.se / VAR_inflation_nBC, by = .(ANALYSIS_NAME, perm, power, boot, Z_th)]
}
### compute SE based pvalues : test log(alpha.GC) = median(log(alpha.GC))
oligo_activity_annot[, pval.SE_norm := 2 * pnorm(abs(log(alpha.GC_norm)), 0, alpha.se_norm, lower.tail = FALSE), by = .(ANALYSIS_NAME, perm, power, boot, Z_th)]
# based on scrambled
oligo_activity_annot[, P_scrambled := 2 * pnorm(abs(log2(alpha.GC_norm)),
  median(log2(alpha.GC_norm[type == "scrambled"]), na.rm = TRUE),
  mad(log2(alpha.GC_norm[type == "scrambled"]), na.rm = TRUE),
  low = FALSE
), by = .(ANALYSIS_NAME, perm, power, boot, Z_th)]

oligo_activity_annot[, bin_logP := cut(-log10(pval.SE_norm), breaks = c(seq(0, 25, by = 0.25), Inf), include.lowest = TRUE)]

# estimate local FDR based on permutations
lfdr_estim <- oligo_activity_annot[, .(n0 = sum(perm != "perm_0_0_0"), n1 = sum(perm == "perm_0_0_0")), keyby = .(Z_th, power, boot, ANALYSIS_NAME, bin_logP)]
lfdr_estim[, lfdr := pmin(1, n0 / n1), keyby = .(Z_th, power, boot, ANALYSIS_NAME, bin_logP)]
lfdr_estim[, monotonic_bins := cumsum(!duplicated(rev(cummax(rev(lfdr))))), keyby = .(Z_th, power, boot, ANALYSIS_NAME)]
lfdr_estim[, lfdr := pmin(1, sum(n0) / sum(n1)), keyby = .(Z_th, power, boot, ANALYSIS_NAME, monotonic_bins)]
oligo_activity_annot <- merge(oligo_activity_annot, lfdr_estim[, .(Z_th, power, boot, ANALYSIS_NAME, bin_logP, lfdr)], by = c("Z_th", "power", "boot", "ANALYSIS_NAME", "bin_logP"), all.x = TRUE)

# compute FDR ordering by pvalues
oligo_activity_annot <- oligo_activity_annot[order(pval.SE_norm), .SD, by = .(Z_th, power, boot, ANALYSIS_NAME)]
oligo_activity_annot[, N_FP := (.1 + cumsum(perm != "perm_0_0_0")), by = .(Z_th, power, boot, ANALYSIS_NAME)]
oligo_activity_annot[, N_POS := (.1 + cumsum(perm == "perm_0_0_0")), by = .(Z_th, power, boot, ANALYSIS_NAME)]
oligo_activity_annot[, pval_emp := N_FP/(.1 + sum(perm == "perm_0_0_0")), by = .(Z_th, power, boot, ANALYSIS_NAME)]
oligo_activity_annot[, FDR := pmin(1, rev(cummin(rev(N_FP / N_POS)))), by = .(Z_th, power, boot, ANALYSIS_NAME)]

# compute FDR on scrambled
oligo_activity_annot[, FDR_scrambled := p.adjust(P_scrambled, "fdr"), by = .(Z_th, perm, power, boot, ANALYSIS_NAME)]

### check number of positives
oligo_activity_annot[FDR < .05 & power == 0 & perm == "perm_0_0_0" & grepl("celltype", ANALYSIS_SUBTYPE), .N, by = ANALYSIS_NAME]
### compute analytic FDR for comparison and assessment of alpha.se
# oligo_activity_annot[, FDR_analytic := p.adjust(pval.SE, "fdr"), by = .(Z_th, perm, power, boot, ANALYSIS_NAME)]
# oligo_activity_annot[, FDR_norm := p.adjust(pval.SE_norm, "fdr"), by = .(Z_th, perm, power, boot, ANALYSIS_NAME)]

### compare number of positives for each analysis
# oligo_activity_annot[, mean(FDR < .05), by = .(Z_th, perm, power, boot, ANALYSIS_NAME)][perm == "perm_0_0_0" & power == 0 & boot == 0 & grepl("all", ANALYSIS_NAME)]
# oligo_activity_annot[, mean(FDR_norm < .05), by = .(Z_th, perm, power, boot, ANALYSIS_NAME)][perm == "perm_0_0_0" & power == 0 & boot == 0 & grepl("all", ANALYSIS_NAME)]
# oligo_activity_annot[, mean(FDR_scrambled < .05), by = .(Z_th, perm, power, boot, ANALYSIS_NAME)][perm == "perm_0_0_0" & power == 0 & boot == 0 & grepl("all", ANALYSIS_NAME)]


# add details on experiment, library, and COND_ID
celline <- c("A549", "HepG2", "K562")
celltype_summary <- data.table(celltype = ifelse(celline == "A549", "A549-ACE2", celline), celline, condition = "all", COND_ID = paste0(celline, "_all"), analysis_name = paste0(ifelse(celline == "A549", "A549-ACE2", celline), "_all"))
analysis_summaries <- rbind(celltype_summary, condition_summary, condition_summary_reps[, mget(colnames(condition_summary))])
setnames(analysis_summaries, "analysis_name", "ANALYSIS_NAME")
oligo_activity_DT <- merge(oligo_activity_annot, analysis_summaries, by = "ANALYSIS_NAME", all.x = TRUE, allow.cartesian = TRUE)
# oligo_activity_DT=merge(oligo_activity_DT, oligo_source[,.(oligo,E071,type_simple)], by="oligo" , all.x=TRUE, allow.cartesian=TRUE)

oligo_activity_DT[, ROADMAP_current_celline := case_when(
  celline == "A549" ~ E114,
  celline == "HepG2" ~ E118,
  celline == "K562" ~ E123,
  TRUE ~ NA
)]
oligo_activity_DT[, ROADMAP_ctrl_brain := E071]

# # define class of oligo
# oligo_activity_DT[, oligo_class := case_when(
#   FDR > .01 ~ "inactive",
#   FDR <= .01 & alpha.GC_norm > 1 & P_scrambled > .01 ~ "enhancer",
#   FDR <= .01 & alpha.GC_norm > 1 & P_scrambled <= .01 ~ "strong enhancer",
#   FDR <= .01 & alpha.GC_norm < 1 & P_scrambled > .01 ~ "silencer",
#   FDR <= .01 & alpha.GC_norm < 1 & P_scrambled <= .01 ~ "strong silencer",
#   TRUE ~ NA
# )]
# oligo_activity_DT[, oligo_class_strict := case_when(
#   FDR > .01 ~ "inactive",
#   FDR <= .01 & alpha.GC_norm > 1 & FDR_scrambled > .01 ~ "enhancer",
#   FDR <= .01 & alpha.GC_norm > 1 & FDR_scrambled <= .01 ~ "strong enhancer",
#   FDR <= .01 & alpha.GC_norm < 1 & FDR_scrambled > .01 ~ "silencer",
#   FDR <= .01 & alpha.GC_norm < 1 & FDR_scrambled <= .01 ~ "strong silencer",
#   TRUE ~ NA
# )]
# oligo_activity_DT[, oligo_class_loose := case_when(
#   FDR > .05 ~ "inactive",
#   FDR <= .01 & alpha.GC_norm > 1 & P_scrambled > .01 ~ "enhancer",
#   FDR <= .01 & alpha.GC_norm > 1 & P_scrambled <= .01 ~ "strong enhancer",
#   FDR <= .05 & alpha.GC_norm > 1 & FDR > .01 ~ "weak enhancer",
#   FDR <= .05 & alpha.GC_norm < 1 & FDR > .01 ~ "weak silencer",
#   FDR <= .01 & alpha.GC_norm < 1 & P_scrambled > .01 ~ "silencer",
#   FDR <= .01 & alpha.GC_norm < 1 & P_scrambled <= .01 ~ "strong silencer",
#   TRUE ~ NA
# )]

# oligo_activity_DT[, oligo_class_localFDR := case_when(
#   lfdr > .2 ~ "inactive",
#   lfdr <= .2 & alpha.GC_norm > 1 & P_scrambled > .01 ~ "enhancer",
#   lfdr <= .2 & alpha.GC_norm > 1 & P_scrambled <= .01 ~ "strong enhancer",
#   lfdr <= .2 & alpha.GC_norm < 1 & P_scrambled > .01 ~ "silencer",
#   lfdr <= .2 & alpha.GC_norm < 1 & P_scrambled <= .01 ~ "strong silencer",
#   TRUE ~ NA
# )]

############

oligo_activity_DT[, alpha_CRITERIA := get(ALPHA_CRITERIA)]
oligo_activity_DT[, log_alpha_se_CRITERIA := get(ALPHA_SE_CRITERIA)]

oligo_activity_DT[, Zscore_CRITERIA := log(get(ALPHA_CRITERIA)) / get(ALPHA_SE_CRITERIA)]

oligo_activity_DT[, oligo_active_CRITERIA := lfdr <= LFDR_TH & FDR < FDR_TH & abs(log2(alpha_CRITERIA)) > LOG2FC_TH]
oligo_activity_DT[, oligo_strong_CRITERIA := oligo_active_CRITERIA & P_scrambled <= SCRAMBLED_TH & abs(log2(alpha_CRITERIA)) > LOG2FC_TH_STRONG]

oligo_activity_DT[, oligo_class_CRITERIA := case_when(
  !oligo_active_CRITERIA ~ "inactive",
  oligo_active_CRITERIA & alpha_CRITERIA > 1 & !oligo_strong_CRITERIA ~ "enhancer",
  oligo_active_CRITERIA & alpha_CRITERIA > 1 & oligo_strong_CRITERIA ~ "strong enhancer",
  oligo_active_CRITERIA & alpha_CRITERIA < 1 & !oligo_strong_CRITERIA ~ "silencer",
  oligo_active_CRITERIA & alpha_CRITERIA < 1 & oligo_strong_CRITERIA ~ "strong silencer",
  TRUE ~ NA
)]

oligo_activity_DT[,  cre_class_CRITERIA := case_when(
    any(oligo_class_CRITERIA == "strong enhancer") ~ "strong enhancer",
    any(oligo_class_CRITERIA == "strong silencer") ~ "strong silencer",
    any(oligo_class_CRITERIA == "enhancer") ~ "enhancer",
    any(oligo_class_CRITERIA == "silencer") ~ "silencer",
    all(oligo_class_CRITERIA == "inactive") ~ "inactive",
    TRUE ~ "WTF"), by = .(type, crsID, ANALYSIS_NAME, ANALYSIS_SUBTYPE, COND_ID,perm,power)]

oligo_activity_obs <- oligo_activity_DT[power == 0 & perm == "perm_0_0_0", ]
oligo_activity_perm <- oligo_activity_DT[power == 0 & perm == "perm_1_0_0", ]

fwrite(oligo_activity_obs, file = sprintf("%s/all_oligos_annotated__%s.tsv.gz", ACTIVE_DIR, CRITERIA_ACTIVE), sep = "\t")
fwrite(oligo_activity_perm, file = sprintf("%s/all_oligos_annotated_perm__%s.tsv.gz", ACTIVE_DIR, CRITERIA_ACTIVE), sep = "\t")

# oligo_activity <- oligo_activity_obs[, .(
#   oligo_class = case_when(
#     any(oligo_class == "strong enhancer") ~ "strong enhancer",
#     any(oligo_class == "strong silencer") ~ "strong silencer",
#     any(oligo_class == "enhancer") ~ "enhancer",
#     any(oligo_class == "silencer") ~ "silencer",
#     all(oligo_class == "inactive") ~ "inactive",
#     TRUE ~ "WTF"
#   ),
#   oligo_class_localFDR = case_when(
#     any(oligo_class_localFDR == "strong enhancer") ~ "strong enhancer",
#     any(oligo_class_localFDR == "strong silencer") ~ "strong silencer",
#     any(oligo_class_localFDR == "enhancer") ~ "enhancer",
#     any(oligo_class_localFDR == "silencer") ~ "silencer",
#     all(oligo_class_localFDR == "inactive") ~ "inactive",
#     TRUE ~ "WTF"
#   ),
#   oligo_class_loose = case_when(
#     any(oligo_class_loose == "strong enhancer") ~ "strong enhancer",
#     any(oligo_class_loose == "strong silencer") ~ "strong silencer",
#     any(oligo_class_loose == "enhancer") ~ "enhancer",
#     any(oligo_class_loose == "silencer") ~ "silencer",
#     any(oligo_class_loose == "weak enhancer") ~ "weak enhancer",
#     any(oligo_class_loose == "weak silencer") ~ "weak silencer",
#     all(oligo_class_loose == "inactive") ~ "inactive",
#     TRUE ~ "WTF"
#   )
# ), by = .(type, crsID, ANALYSIS_NAME, ANALYSIS_SUBTYPE, COND_ID)]


oligo_activity <- oligo_activity_obs[, .(
  oligo_class_CRITERIA = case_when(
    any(oligo_class_CRITERIA == "strong enhancer") ~ "strong enhancer",
    any(oligo_class_CRITERIA == "strong silencer") ~ "strong silencer",
    any(oligo_class_CRITERIA == "enhancer") ~ "enhancer",
    any(oligo_class_CRITERIA == "silencer") ~ "silencer",
    all(oligo_class_CRITERIA == "inactive") ~ "inactive",
    TRUE ~ "WTF"
  )
), by = .(type, crsID, ANALYSIS_NAME, ANALYSIS_SUBTYPE, COND_ID)]

fwrite(oligo_activity, file = sprintf("%s/oligo_activity__all__%s.tsv.gz", ACTIVE_DIR, CRITERIA_ACTIVE), sep = "\t")

fwrite(oligo_activity[ANALYSIS_SUBTYPE == "celltype_cond"], file = sprintf("%s/oligo_activity__celltype_cond__%s.tsv.gz", ACTIVE_DIR, CRITERIA_ACTIVE), sep = "\t")


#  oligo_class_by_COND <- oligo_activity[grepl("all", ANALYSIS_NAME), .N, by = .(COND_ID, oligo_class_loose)]
#  oligo_class_by_COND[, Pct := N / sum(N), by = .(COND_ID)]


cat("All done !")
q("no")
