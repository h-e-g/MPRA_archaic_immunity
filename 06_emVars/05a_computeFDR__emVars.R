MPRA_DIR <- "/pasteur/helix/projects/evo_immuno_pop/MPRA"
ANALYSIS_DIR <- "MPRA_count_exp6_analysis3"

source(sprintf("%s/scripts/%s/00_starter.R", MPRA_DIR, ANALYSIS_DIR))
source(sprintf("%s/scripts/%s/00_load_associations.R", MPRA_DIR, ANALYSIS_DIR))
source(sprintf("%s/scripts/misc_plots.R", MPRA_DIR))

RUN_ID <- "RUN3_Z2_nBC10"
CRITERIA_ACTIVE <- "FDR5_scrambled5pct_FC0.2_GCnorm"
FILTER_ACTIVE <- TRUE
CRITERIA_EMVARS <- "EmVar_FDR5_FC.2"

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
   CRITERIA_ACTIVE_OUT <- paste0('noActivityFilter_',CRITERIA_ACTIVE)
 }else{
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
selected_annot_wide[,Introgression_source_top:=Introgression_source_top_initial]
selected_annot_wide[Introgression_source_top=='both',Introgression_source_top:='Vindija/Denisova']
selected_annot_wide[Introgression_source_top=='',Introgression_source_top:='Undetermined']

SNP_annot_v5 <- merge(SNP_annot_v4[,-"Introgression_scenario"],selected_annot_wide,by=c('ID','posID'))
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


###############################################################
################# count emVars FDR ############################	
###############################################################
		
emVARs_all_results <- fread(sprintf("%s/all_emVars_results.tsv.gz", IN_DIR), sep = "\t")
emVARs_all_results <- emVARs_all_results[crsID %in% tested_and_ctrl_oligos_final[type=='tested',posID], ]
# oligo_activity_file <- sprintf("%s/oligo_activity__all.tsv.gz", OUT_DIR)
# oligo_activity <- fread(oligo_activity_file)

oligo_target <- fread(sprintf("%s/data/%s/00_oligo_targets_v4.txt", MPRA_DIR, ANALYSIS_DIR), sep = "\t")
oligo_target_collapsed <- fread(sprintf("%s/data/%s/00_oligo_targets_collapsed_v4.txt", MPRA_DIR, ANALYSIS_DIR), sep = "\t")
crs_targets <- unique(oligo_target_collapsed[, -"oligo"])

####################################################################
######################## compute FDR, emVars #######################
####################################################################

emVARs_annot <- merge(emVARs_all_results[df.test > 0, ], SNP_annot_v5, by = c("crsID"), allow.cartesian = TRUE)
emVARs_annot[, celltype := str_split(ANALYSIS_NAME, "_|-", simplify = T)[, 1]]
emVARs_annot <- merge(emVARs_annot, crs_targets, by = c("crsID", "celltype"), all.x = TRUE)

oligo_activity_obs <- fread(sprintf("%s/all_oligos_annotated__%s.tsv.gz", ACTIVE_DIR, CRITERIA_ACTIVE))

if (FILTER_ACTIVE==TRUE) {
  # subset to cres that are active in at least one condition
  all_active_cres <- oligo_activity_obs[oligo_active_CRITERIA == TRUE, .(oligo, crsID, ANALYSIS_NAME, ANALYSIS_SUBTYPE)]
  emVARs_annot[, is_tested_cre := paste(crsID, ANALYSIS_NAME) %chin% all_active_cres[, paste(crsID, ANALYSIS_NAME)]]
} else {
  emVARs_annot[, is_tested_cre := TRUE]
}

emVARs_annot <- merge(emVARs_annot, unique(oligo_activity_obs[, .(cre_class_CRITERIA, posID = crsID, ANALYSIS_NAME)])
  , by = c("posID", "ANALYSIS_NAME"), all.x = TRUE, suffix = c(".emVAR", ".oligo"))
# emVARs_annot <- emVARs_annot[is_tested_cre == TRUE, ]

if (any(grepl(".x$", colnames(emVARs_annot)))) {
  stop("merge issue")
}

emVARs_annot[, log2FC_der_vs_anc := ifelse(allele.1 == ANCESTRAL, logFC_2vs1 / log(2), -logFC_2vs1 / log(2))]
emVARs_annot[, log2FC_archaic_vs_modern := ifelse(allele.2 == INTROGRESSED.allele, logFC_2vs1 / log(2), -logFC_2vs1 / log(2))]
emVARs_annot[, log2FC.se := logFC.se / log(2)]

emVARs_annot[, bin_logP := cut(-log10(pval.LRT), breaks = c(seq(0, 5, by = 0.25), Inf))]

lfdr_estim <- emVARs_annot[, .(nPerm = sum(perm != "perm_0_0_0" & df.test > 0 & is_tested_cre), nObs = sum(perm == "perm_0_0_0" & df.test > 0 & is_tested_cre)), keyby = .(Z_th, power, boot, ANALYSIS_NAME, bin_logP)]
lfdr_estim[, lfdr := pmin(1, nPerm / nObs), keyby = .(Z_th, power, boot, ANALYSIS_NAME, bin_logP)]
lfdr_estim[, monotonic_bins := cumsum(!duplicated(rev(cummax(rev(lfdr))))), keyby = .(Z_th, power, boot, ANALYSIS_NAME)]
lfdr_estim[, lfdr := pmin(1, sum(nPerm) / sum(nObs)), keyby = .(Z_th, power, boot, ANALYSIS_NAME, monotonic_bins)]

emVARs_annot <- merge(emVARs_annot, lfdr_estim[, .(Z_th, power, boot, ANALYSIS_NAME, bin_logP, lfdr)], by = c("Z_th", "power", "boot", "ANALYSIS_NAME", "bin_logP"), all.x = TRUE)
emVARs_annot[is_tested_cre==FALSE, lfdr := 1]

emVARs_annot <- emVARs_annot[order(pval.LRT), .SD, by = .(Z_th, power, boot, ANALYSIS_NAME)]
emVARs_annot[, N_FP := (.1 + cumsum(perm != "perm_0_0_0" & df.test > 0 & is_tested_cre)), by = .(Z_th, power, boot, ANALYSIS_NAME)]
emVARs_annot[, N_POS := (.1 + cumsum(perm == "perm_0_0_0" & df.test > 0 & is_tested_cre)), by = .(Z_th, power, boot, ANALYSIS_NAME)]
emVARs_annot[, pval_emp := N_FP/(.1 + sum(perm != "perm_0_0_0" & df.test > 0 & is_tested_cre)), by = .(Z_th, power, boot, ANALYSIS_NAME)]
emVARs_annot[is_tested_cre==TRUE, FDR := pmin(1, rev(cummin(rev(N_FP / N_POS)))), by = .(Z_th, power, boot, ANALYSIS_NAME)]
emVARs_annot[is_tested_cre==FALSE, FDR := 1]
emVARs_annot[FDR < .05 & power == 0 & perm == "perm_0_0_0" & grepl("celltype", ANALYSIS_SUBTYPE), .N, by = ANALYSIS_NAME]

emVARs_annot[, is_emVar_CRITERIA := FDR < FDR_TH_EMVAR & lfdr < LFDR_TH_EMVAR & abs(log2FC_archaic_vs_modern) > LOG2FC_TH_EMVAR]

emVARs_obs <- emVARs_annot[power == 0 & perm == "perm_0_0_0", ]
emVARs_perm <- emVARs_annot[power == 0 & perm == "perm_0_0_1", ]
emVARs_power <- emVARs_annot[power == 1, ]

fwrite(emVARs_obs, file = sprintf("%s/all_emVars_annotated__%s.tsv.gz", EMVAR_DIR, CRITERIA_ACTIVE_OUT), sep = "\t")
fwrite(emVARs_perm, file = sprintf("%s/all_emVars_annotated_perm__%s.tsv.gz", EMVAR_DIR, CRITERIA_ACTIVE_OUT), sep = "\t")
fwrite(emVARs_power, file = sprintf("%s/all_emVars_power__%s.tsv.gz", EMVAR_DIR, CRITERIA_ACTIVE_OUT), sep = "\t")

emVARs_obs_ctc <- emVARs_obs[ANALYSIS_SUBTYPE == "celltype_cond"]
emVARs_obs_ctc <- merge(emVARs_obs_ctc, condition_summary[, -"celltype"], by.x = "ANALYSIS_NAME", by.y = "analysis_name")

emVARs_obs_sample <- emVARs_obs[ANALYSIS_SUBTYPE == "sample"]
emVARs_obs_sample <- merge(emVARs_obs_sample, condition_summary_reps[, -"celltype"], by.x = "ANALYSIS_NAME", by.y = "analysis_name")

fwrite(emVARs_obs_ctc, file = sprintf("%s/all_emVars_annotated_celltype__%s.tsv", EMVAR_DIR, CRITERIA_ACTIVE_OUT), sep = "\t")
fwrite(emVARs_obs_sample, file = sprintf("%s/all_emVars_by_sample_annotated__%s.tsv", EMVAR_DIR, CRITERIA_ACTIVE_OUT), sep = "\t")

# emVARs_annot_active <- emVARs_annot[oligo_class_loose != "inactive", ]
# emVARs_annot_active[, N_FP := (.1 + cumsum(perm != "perm_0_0_0" & df.test > 0)), by = .(Z_th, power, boot, ANALYSIS_NAME)]
# emVARs_annot_active[, N_POS := (.1 + cumsum(perm == "perm_0_0_0" & df.test > 0)), by = .(Z_th, power, boot, ANALYSIS_NAME)]
# emVARs_annot_active[, FDR := pmin(1, rev(cummin(rev(N_FP / N_POS)))), by = .(Z_th, power, boot, ANALYSIS_NAME)]
# emVARs_annot_active[FDR < .05 & power == 0 & perm == "perm_0_0_0" & grepl("celltype", ANALYSIS_SUBTYPE) & abs(log2FC_archaic_vs_modern) > 0.2, .N, by = ANALYSIS_NAME]

# fwrite(emVARs_annot_active[power == 0 & perm == "perm_0_0_0", ], file = sprintf("%s/all_emVars_activeCRS_annotated.tsv.gz", OUT_DIR), sep = "\t")
# fwrite(emVARs_annot_active[power == 0 & perm == "perm_0_0_0" & grepl("celltype", ANALYSIS_SUBTYPE), ], file = sprintf("%s/all_emVars_activeCRS_annotated_celltype_v2.tsv", OUT_DIR), sep = "\t")
# fwrite(emVARs_annot_active[power == 0 & perm == "perm_0_0_0" & grepl("sample", ANALYSIS_SUBTYPE), ], file = sprintf("%s/all_emVars_by_sample_activeCRS_annotated.tsv", OUT_DIR), sep = "\t")

# emVARs_annot_active_final <- merge(emVARs_annot_active, SNP_annot_v4[, .(ID, gwas, strongest_effect_in, effect_in, has_eQTL, has_gwas, blood_eQTL, lung_eQTL, liver_eQTL)], by = "ID", allow.cartesian = TRUE, all.x = TRUE)
# emVARs_annot_active_final <- merge(emVARs_annot_active, crs_targets, by = c("crsID", "celltype"), allow.cartesian = TRUE, all.x = TRUE)
# fwrite(emVARs_annot_active_final[power == 0 & perm == "perm_0_0_0" & grepl("celltype", ANALYSIS_SUBTYPE), ], file = sprintf("%s/all_emVars_activeCRS_annotated_celltype_v2.tsv", OUT_DIR), sep = "\t")
# emVARs_annot_active_final <- fread(sprintf("%s/all_emVars_activeCRS_annotated_celltype_v2.tsv", OUT_DIR))

# emVARs_annot_active_ok <- emVARs_annot[oligo_class_loose != "inactive" & excluded == FALSE & (!is_introgressed | INTROGRESSED.allele != ""), ]
# emVARs_annot_active_ok[, N_FP := (.1 + cumsum(perm != "perm_0_0_0" & df.test > 0)), by = .(Z_th, power, boot, ANALYSIS_NAME)]
# emVARs_annot_active_ok[, N_POS := (.1 + cumsum(perm == "perm_0_0_0" & df.test > 0)), by = .(Z_th, power, boot, ANALYSIS_NAME)]
# emVARs_annot_active_ok[, FDR := pmin(1, rev(cummin(rev(N_FP / N_POS)))), by = .(Z_th, power, boot, ANALYSIS_NAME)]
# emVARs_annot_active_ok[FDR < .05 & power == 0 & perm == "perm_0_0_0" & grepl("celltype", ANALYSIS_SUBTYPE) & abs(log2FC_archaic_vs_modern) > 0.2, .N, by = ANALYSIS_NAME]

# fwrite(emVARs_annot_active_ok[power == 0 & perm == "perm_0_0_0", ], file = sprintf("%s/all_emVars_activeCRS_ok_annotated.tsv.gz", OUT_DIR), sep = "\t")
# fwrite(emVARs_annot_active_ok[power == 0 & perm == "perm_0_0_0" & grepl("celltype", ANALYSIS_SUBTYPE), ], file = sprintf("%s/all_emVars_activeCRS_ok_annotated_celltype.tsv", OUT_DIR), sep = "\t")
# fwrite(emVARs_annot_active_ok[power == 0 & perm == "perm_0_0_0" & grepl("sample", ANALYSIS_SUBTYPE), ], file = sprintf("%s/all_emVars_by_sample_activeCRS_ok_annotated.tsv", OUT_DIR), sep = "\t")

# # tested SNPx cond, SNP, introgression_locus
# emVARs_annot[power == 0 & perm == "perm_0_0_0" & grepl("celltype_cond", ANALYSIS_SUBTYPE) & !ctrl & !excluded, .(.(N_snp_x_cond = .N, N_snp = length(unique(crsID)), N_loci = length(unique(introgression_locus))))]
# # tested SNPx cond, SNP, introgression_locus (excluding ambuguous introgressed alleles)
# emVARs_annot[power == 0 & perm == "perm_0_0_0" & grepl("celltype_cond", ANALYSIS_SUBTYPE) & !ctrl & !excluded & (!is_introgressed | INTROGRESSED.allele != ""), .(N_snp_x_cond = .N, N_snp = length(unique(crsID)), N_loci = length(unique(introgression_locus)))]
# # active crs x cond, crs, introgression_locus (excluding ambuguous introgressed alleles)
# emVARs_annot[power == 0 & perm == "perm_0_0_0" & grepl("celltype_cond", ANALYSIS_SUBTYPE) & !ctrl & !excluded & (!is_introgressed | INTROGRESSED.allele != "") & oligo_class_loose != "inactive", .(N_snp_x_cond = .N, N_snp = length(unique(crsID)), N_loci = length(unique(introgression_locus)))]
# emVARs_annot_active_ok[power == 0 & perm == "perm_0_0_0" & grepl("celltype_cond", ANALYSIS_SUBTYPE) & !ctrl & !excluded & (!is_introgressed | INTROGRESSED.allele != "") & oligo_class_loose != "inactive", .(N_snp_x_cond = .N, N_snp = length(unique(crsID)), N_loci = length(unique(introgression_locus)))]
# # emVars SNP x cond, SNP, introgression_locus (excluding ambuguous introgressed alleles)
# emVARs_annot_active_ok[FDR < .05 & abs(log2FC_archaic_vs_modern) > 0.2 & power == 0 & perm == "perm_0_0_0" & grepl("celltype_cond", ANALYSIS_SUBTYPE) & !ctrl & !excluded & (!is_introgressed | INTROGRESSED.allele != ""), .(N_snp_x_cond = .N, N_snp = length(unique(crsID)), N_loci = length(unique(introgression_locus)))]

#emVARs_annot_full=merge(emVARs_obs,SNP_annot_v4[,.(ID,gwas,strongest_effect_in,effect_in,has_eQTL,has_gwas,blood_eQTL,lung_eQTL,liver_eQTL)],by='ID',allow.cartesian=TRUE,all.x=TRUE)
#emVARs_annot_full=merge(emVARs_annot_full,crs_targets,by=c('crsID','celltype'),allow.cartesian=TRUE,all.x=TRUE)


emVARs_obs_ctc_full=merge(emVARs_obs_ctc[,-c('Source_introgression','allele_match')], selected_annot,by=c('ID','posID'),allow.cartesian=TRUE,all.x=TRUE)

if (any(grepl(".x$", colnames(emVARs_obs_ctc_full)))) {
  stop("merge issue")
}

fwrite(emVARs_obs_ctc_full, file = sprintf("%s/all_emVars_annotated_full_celltype__%s.tsv", EMVAR_DIR, CRITERIA_ACTIVE_OUT), sep = "\t")
# emVARs_annot_full <- fread(sprintf("%s/all_emVars_annotated_full_celltype__%s.tsv", EMVAR_DIR, CRITERIA_ACTIVE_OUT))


# emVar_FC <- dcast(emVARs_annot[power == 0 & perm == "perm_0_0_0", ], posID ~ ANALYSIS_NAME, value.var = "log2FC_archaic_vs_modern", fun.aggregate = mean, na.rm = TRUE)
SNP_annot_emVars <- emVARs_obs_ctc[, .(
  nCelltype = length(unique(celltype[is_emVar_CRITERIA])),
  nCelltype_P01 = length(unique(celltype[pval_emp<.01])),
  nCelltype_P05 = length(unique(celltype[pval_emp<.05])),
  nCelltype_NS = length(unique(celltype[is_emVar_CRITERIA & condition=='NS'])),
  nCelltype_P01_NS = length(unique(celltype[pval_emp<.01 & condition=='NS'])),
  nCelltype_P05_NS = length(unique(celltype[pval_emp<.05 & condition=='NS'])),
	nContext = length(unique(COND_ID[is_emVar_CRITERIA])),
  nContext_P01 = length(unique(COND_ID[pval_emp<.01])),
  nContext_P05 = length(unique(COND_ID[pval_emp<.05])),
	nContext_exp = sum(1 - lfdr),
  best_FC_FDR5 = max(c(abs(log2FC_archaic_vs_modern)[is_emVar_CRITERIA], 0)),
  best_expFC = max(c(abs(log2FC_archaic_vs_modern) * (1 - lfdr))),
  best_context_specificty_BoW = sd(log2FC_archaic_vs_modern, na.rm = TRUE) / mean(log2FC.se, na.rm = TRUE)
), by = .(ID, posID)]

fwrite(SNP_annot_emVars, file = sprintf("%s/SNP_annot_emVars__%s.tsv", EMVAR_DIR, CRITERIA_ACTIVE_OUT), sep = "\t")

cat("\nAll done\n")
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
