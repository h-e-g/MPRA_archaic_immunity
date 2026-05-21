# running: sbatch -p geh,common --mem=30G  00_Rscript.sh MPRA_count_exp6_analysis3/02c_aggMPRAnalyse.R

# TODO :  add nBC filters
# DONE :  split analysis of CRS and emVars

MPRA_DIR <- "/pasteur/helix/projects/evo_immuno_pop/MPRA"
ANALYSIS_DIR <- "MPRA_count_exp6_analysis3"

source(sprintf("%s/scripts/%s/00_starter.R", MPRA_DIR, ANALYSIS_DIR))
source(sprintf("%s/scripts/%s/00_load_associations.R", MPRA_DIR, ANALYSIS_DIR))


library(MPRAnalyze)

FIGURE_DIR <- sprintf("%s/figures/%s", MPRA_DIR, ANALYSIS_DIR)
VERSION_IN <- 1
VERSION_OUT <- VERSION_IN
sub_sample <- 2000
RUN_ID <- "RUN2_Z2_nBC10"
# minBC per oligo
minBC <- 10
instruction_file <- sprintf("%s/scripts/%s/02a_MPRA_analyse_instructions_full_sorted.txt", MPRA_DIR, ANALYSIS_DIR)


cmd <- commandArgs(trailingOnly = TRUE)
print(cmd)

for (i in seq_along(cmd)) {
  if (cmd[i] == "--version_in" || cmd[i] == "-vi") {
    VERSION_IN <- cmd[i + 1]
  } # version of input files
  if (cmd[i] == "--version_out" || cmd[i] == "-vo") {
    VERSION_OUT <- cmd[i + 1]
  } # version of output files
  if (cmd[i] == "--sub_sample" || cmd[i] == "-n") {
    sub_sample <- as.numeric(cmd[i + 1])
  } # max number of barcodes to use per oligo
  if (cmd[i] == "--min_nBC" || cmd[i] == "-n") {
    min_nBC <- as.numeric(cmd[i + 1])
  } # min number of barcodes to use per oligo
  if (cmd[i] == "--instruction_file" || cmd[i] == "-i") {
    instruction_file <- cmd[i + 1]
  } # instruction_file containing instructions for MPRAnalyze
  if (cmd[i] == "--run_dir" || cmd[i] == "-r") {
    RUN_ID <- cmd[i + 1]
  }
}

tic("loading oligo & SNP annotations")
# load annotation of oligos
oligo_source <- fread(sprintf("%s/data/%s/00_oligo_annot_v2.txt", MPRA_DIR, ANALYSIS_DIR))
# define order of priority and assign each oligo to a single annotation
oligo_source <- oligo_source[order(factor(type, c(
  "archaic", "COVID", "Expecto", "MPRA", "eQTLs", "eQTL_repeat", "Purged", "randomSNP",
  "Promoter_Lung", "Promoter_Tcell", "Promoter_Mono_NS", "Promoter_Mono_STIM", "scrambled"
))), ]
oligo_source <- oligo_source[!duplicated(oligo), ]
# load annotations of SNPs
SNP_annot <- fread(sprintf("%s/data/%s/00_SNP_annot_v1.txt", MPRA_DIR, ANALYSIS_DIR))
toc()


#### aggregate MPRAnalyze
instruction_DT <- fread(instruction_file, header = TRUE, fill = TRUE)
if (all(is.na(instruction_DT[, group2_samples]))) {
  instruction_DT[, group2_samples := NULL]
  instruction_DT[, group2_samples := ""]
}
oligo_activity <- list()
emVARs <- list()
oligo_activity_Diff <- list()
emVARs_Diff <- list()

INPUT_DIR <- sprintf("%s/data/%s/02b_runMPRAnalyze/%s/", MPRA_DIR, ANALYSIS_DIR, RUN_ID)
for (instruction_line in seq_len(instruction_DT[, .N])) {
  ANALYSIS_NAME <- instruction_DT[instruction_line, analysis_name]
  ANALYSIS_SUBTYPE <- instruction_DT[instruction_line, analysis_subtype]
  ANALYSIS_TYPE <- instruction_DT[instruction_line, analysis_type]
  GROUP2_SAMPLES <- instruction_DT[instruction_line, group2_samples]

  IN_DIR <- sprintf("%s/%s/%s/%s/", INPUT_DIR, ANALYSIS_TYPE, ANALYSIS_SUBTYPE, ANALYSIS_NAME)
  IN_FILES <- dir(IN_DIR, recursive = TRUE, pattern = "alphas_and_pvalues", full.names = FALSE)
  results_MPRA <- list()
  for (IN_FILE in IN_FILES) {
    PROV <- try(fread(file = sprintf("%s/%s", IN_DIR, IN_FILE)))
    FILE_CHAR <- as.character(str_split(IN_FILE, "/", simplify = T))
    reg_expr_filename <- "alphas_and_pvalues_v(.*)_(SNPs|oligos)_(chunk[0-9]+)_max([0-9]+)BC_Z([0-9]+).txt"
    VERSION_IN <- gsub(reg_expr_filename, "\\1", FILE_CHAR[length(FILE_CHAR)])
    DATA_TYPE <- gsub(reg_expr_filename, "\\2", FILE_CHAR[length(FILE_CHAR)])
    Chunk_nb <- gsub(reg_expr_filename, "\\3", FILE_CHAR[length(FILE_CHAR)])
    sub_sample <- gsub(reg_expr_filename, "\\4", FILE_CHAR[length(FILE_CHAR)])
    if (nrow(PROV) > 0) {
      PROV[, perm := FILE_CHAR[1]]
      if (grepl("power", IN_FILE)) {
        PROV[, power := gsub("power(.*)_boot(.*)", "\\1", FILE_CHAR[2])]
        PROV[, boot := gsub("power(.*)_boot(.*)", "\\2", FILE_CHAR[2])]
      } else {
        PROV[, power := "0"]
        PROV[, boot := "0"]
      }
      PROV[, Z_th := as.numeric(gsub(reg_expr_filename, "\\5", FILE_CHAR[length(FILE_CHAR)]))]
    }
    results_MPRA[[IN_FILE]] <- PROV
  }
  results_MPRA <- rbindlist(results_MPRA)
  if (ANALYSIS_TYPE == "oligo_activity") {
    if (GROUP2_SAMPLES == "") {
      oligo_activity[[ANALYSIS_NAME]] <- results_MPRA
    } else {
      oligo_activity_Diff[[ANALYSIS_NAME]] <- results_MPRA
    }
  } else {
    if (GROUP2_SAMPLES == "") {
      emVARs[[ANALYSIS_NAME]] <- results_MPRA
    } else {
      emVARs_Diff[[ANALYSIS_NAME]] <- results_MPRA
    }
  }
  cat(ANALYSIS_TYPE, ANALYSIS_NAME, ":", results_MPRA[, .N], "\n")
}
oligo_activity <- rbindlist(oligo_activity)
oligo_activity_Diff <- rbindlist(oligo_activity_Diff)
emVARs <- rbindlist(emVARs)
emVARs_Diff <- rbindlist(emVARs_Diff)
setnames(emVARs, "logFC", "logFC_2vs1")
setnames(emVARs_Diff, "delta_logFC", "delta_logFC_2vs1")
setnames(oligo_activity_Diff, "logFC", "logFC_2vs1")

OUT_DIR <- sprintf("%s/data/%s/03a_aggMPRA_analyse/%s/", MPRA_DIR, ANALYSIS_DIR, RUN_ID)
dir.create(OUT_DIR, recursive = TRUE)

fwrite(emVARs, file = sprintf("%s/all_emVars_results.tsv.gz", OUT_DIR), sep = "\t")
fwrite(emVARs_Diff, file = sprintf("%s/all_emVars_Diff_results.tsv.gz", OUT_DIR), sep = "\t")
fwrite(oligo_activity, file = sprintf("%s/all_oligo_results.tsv.gz", OUT_DIR), sep = "\t")
fwrite(oligo_activity_Diff, file = sprintf("%s/all_oligo_diff_results.tsv.gz", OUT_DIR), sep = "\t")

cat('All done\n')
q('no')


























emVARs_annot_active[type %in% c("archaic", "COVID") & power == 0 & perm == "perm_0_0_0" & FDR < 0.05 & grepl("celltype_cond", ANALYSIS_SUBTYPE), length(unique(posID))]
emVARs_annot_active[FDR < 0.01 & type %in% c("archaic", "COVID") & power == 0 & perm == "perm_0_0_0", length(unique(posID))]


####################################################################
################# compute FDR,  response emVars ####################
####################################################################

emVARs_Diff[, .(.N, sum(df.test_int > 0)), keyby = .(power, boot, ANALYSIS_NAME, Z_th, perm)]
emVARs_Diff_annot <- merge(emVARs_Diff[df.test_int > 0, ], SNP_annot, by = "posID", allow.cartesian = TRUE)
emVARs_Diff_annot[, Delta_log2FC_der_vs_anc := ifelse(allele.1 == ANCESTRAL, delta_logFC_2vs1 / log(2), -delta_logFC_2vs1 / log(2))]
emVARs_Diff_annot[, Delta_log2FC_archaic_vs_modern := ifelse(!grepl("ancestral reintrogressed", SNP_type), Delta_log2FC_der_vs_anc, -Delta_log2FC_der_vs_anc)]
emVARs_Diff_annot[, delta_log2FC.se := delta_logFC.se / log(2)]

emVARs_Diff_annot <- merge(emVARs[FDR <= 0.05, .(Z_th, power, boot, ANALYSIS_NAME, posID)], emVARs_Diff_annot, by = c("Z_th", "power", "boot", "ANALYSIS_NAME", "posID"), all = TRUE)
emVARs_Diff_annot <- emVARs_Diff_annot[order(pval_int), .SD, by = .(Z_th, power, boot, ANALYSIS_NAME)]
emVARs_Diff_annot[, N_FP := (.1 + cumsum(perm != "perm_0_0_0" & df.test_int > 0)), by = .(Z_th, power, boot, ANALYSIS_NAME)]
emVARs_Diff_annot[, N_POS := (.1 + cumsum(perm == "perm_0_0_0" & df.test_int > 0)), by = .(Z_th, power, boot, ANALYSIS_NAME)]
emVARs_Diff_annot[, FDR := pmin(1, rev(cummin(rev(N_FP / N_POS)))), by = .(Z_th, power, boot, ANALYSIS_NAME)]
# emVARs_Diff_annot[FDR < .05 & power == 0 & perm == "perm_0_0_0" & grepl("celltype", ANALYSIS_SUBTYPE), .N, by = ANALYSIS_NAME]
fwrite(emVARs_Diff_annot[power == 0 & perm == "perm_0_0_0", ], file = sprintf("%s/all_emVars_diff_annotated.tsv.gz", OUT_DIR), sep = "\t")
fwrite(emVARs_Diff_annot[power == 0 & perm == "perm_0_0_0" & grepl("response", ANALYSIS_SUBTYPE), ], file = sprintf("%s/all_emVars_diff_annotated.tsv", OUT_DIR), sep = "\t")

####################################################################
################# compute FDR,  response emVars ####################
#################### significant emVars only #######################
####################################################################

# extract significant emVars at 5% FDR
emVars_signif <- emVARs_annot[power == 0 & perm == "perm_0_0_0" & grepl("celltype_cond", ANALYSIS_SUBTYPE) & FDR < 0.05, .(posID, ANALYSIS_NAME, oligo1, oligo2, GC, pval.LRT, log2FC_archaic_vs_modern, log2FC.se,
  FDR_emVAR = FDR, oligo_class = oligo_class_loose, type.oligo, Lung, Mono, Nearest, Tcell, chromosome, position, allele.1, allele.2, strand, allele2_is_REF, rsID, Introgressed_from, SNP_type,
  ANCESTRAL, DERIVED, INTROGRESSED, REF, ALT, Vindija.der, Chagyrskaya.der, Altai.der, Denisova.der, NAF_V, NAF_D, MaxPosterior_D_VY, MaxPosterior_V_DY,
  YRI.der, SGDP_African.der, GBR.der, IBS.der, SGDP_WestEurasian.der, CHB.der, JPT.der, SGDP_Eastasian.der, PJL.der, STU.der, SGDP_Agta.der, SGDP_Papuan.der
)]
emVars_signif <- merge(emVars_signif, condition_summary, by.x = "ANALYSIS_NAME", by.y = "analysis_name")

# compute FDR for response emVars, focusing on significant emVars
response_emVars <- emVARs_Diff_annot[power == 0 & grepl("response", ANALYSIS_SUBTYPE), ]
response_emVars[, analysis_name := gsub("response_(.*)", "\\1_all", ANALYSIS_NAME)]
response_emVars <- merge(response_emVars, condition_summary, by = "analysis_name")
response_emVars <- merge(emVars_signif, response_emVars, by = c("posID", "COND_ID"), suffix = c(".emVar", ".response"))
response_emVars <- response_emVars[paste(posID, COND_ID) %chin% response_emVars[perm == "perm_0_0_0", paste(posID, COND_ID)], ]
response_emVars <- response_emVars[order(pval_int), .SD, by = ANALYSIS_NAME.response]
response_emVars[, N_FP := (.1 + cumsum(perm != "perm_0_0_0" & df.test_int > 0)), by = .(ANALYSIS_NAME.response)]
response_emVars[, N_POS := (.1 + cumsum(perm == "perm_0_0_0" & df.test_int > 0)), by = .(ANALYSIS_NAME.response)]
response_emVars[, FDR := pmin(1, rev(cummin(rev(N_FP / N_POS)))), by = .(ANALYSIS_NAME.response)]
fwrite(response_emVars[power == 0 & perm == "perm_0_0_0", ], file = sprintf("%s/all_response_emVars.tsv", OUT_DIR), sep = "\t")

emVars_signif[!grepl("weak|inactive", oligo_class) & type.oligo %in% c("archaic", "COVID"), length(unique(posID))]
####################################################################
################# compute FDR,  response emVars ####################
############ active oligo with significant emVars only ###############
####################################################################

# extract significant emVars at 5% FDR
emVars_signif_active <- emVARs_annot_active[power == 0 & perm == "perm_0_0_0" & grepl("celltype_cond", ANALYSIS_SUBTYPE) & FDR < 0.05, .(posID, ANALYSIS_NAME, oligo1, oligo2, GC, pval.LRT, log2FC_archaic_vs_modern, log2FC.se,
  FDR_emVAR = FDR, oligo_class_loose,
  type.emVAR, type.oligo
  # ,Lung, Mono, Nearest, Tcell, chromosome, position, allele.1, allele.2, strand, allele2_is_REF, rsID, Introgressed_from, SNP_type,
  # ANCESTRAL, DERIVED, INTROGRESSED, REF, ALT, Vindija.der, Chagyrskaya.der, Altai.der, Denisova.der, NAF_V, NAF_D, MaxPosterior_D_VY, MaxPosterior_V_DY,
  # YRI.der, SGDP_African.der, GBR.der, IBS.der, SGDP_WestEurasian.der, CHB.der, JPT.der, SGDP_Eastasian.der, PJL.der, STU.der, SGDP_Agta.der, SGDP_Papuan.der
)]
emVars_signif_active <- merge(emVars_signif_active, condition_summary, by.x = "ANALYSIS_NAME", by.y = "analysis_name")

# compute FDR for response emVars, focusing on significant emVars
response_emVars_active <- emVARs_Diff_annot[power == 0 & grepl("response", ANALYSIS_SUBTYPE), ]
response_emVars_active[, analysis_name := gsub("response_(.*)", "\\1_all", ANALYSIS_NAME)]
response_emVars_active <- merge(response_emVars_active, condition_summary, by = "analysis_name")
response_emVars_active <- merge(emVars_signif_active, response_emVars_active, by = c("posID", "COND_ID"), suffix = c(".emVar", ".response"))
response_emVars_active <- response_emVars_active[paste(posID, COND_ID) %chin% response_emVars[perm == "perm_0_0_0", paste(posID, COND_ID)], ]
response_emVars_active <- response_emVars_active[order(pval_int), .SD, by = ANALYSIS_NAME.response]
response_emVars_active[, N_FP := (.1 + cumsum(perm != "perm_0_0_0" & df.test_int > 0)), by = .(ANALYSIS_NAME.response)]
response_emVars_active[, N_POS := (.1 + cumsum(perm == "perm_0_0_0" & df.test_int > 0)), by = .(ANALYSIS_NAME.response)]
response_emVars_active[, FDR := pmin(1, rev(cummin(rev(N_FP / N_POS)))), by = .(ANALYSIS_NAME.response)]
fwrite(response_emVars_active[power == 0 & perm == "perm_0_0_0", ], file = sprintf("%s/all_response_emVars_active.tsv", OUT_DIR), sep = "\t")

response_emVars_active_signif <- response_emVars_active[FDR < .05 & type.oligo %in% c("archaic", "COVID") & power == 0 & perm == "perm_0_0_0", ][order(posID), .(posID, rsID, COND_ID, oligo_class_loose, type.oligo, pval.LRT, log2FC_archaic_vs_modern, log2FC.se, FDR_emVAR, Nearest, Lung, Mono, Tcell, SNP_type, Introgressed_from, INTROGRESSED, Delta_log2FC_archaic_vs_modern, delta_log2FC.se, pval_int, FDR, POP_adaptive)]
fwrite(response_emVars_active_signif, file = sprintf("%s/all_response_emVars_active_signif.tsv", OUT_DIR), sep = "\t")

# extract significant emVars at 5% FDR
emVars_signif_active <- emVARs_annot_active[power == 0 & perm == "perm_0_0_0" & grepl("celltype_cond", ANALYSIS_SUBTYPE) & FDR < 0.05, .(posID, ANALYSIS_NAME, oligo1, oligo2, GC, pval.LRT, log2FC_archaic_vs_modern, log2FC.se,
  FDR_emVAR = FDR, oligo_class_loose, oligo_class,
  type.emVAR, type.oligo, Lung, Mono, Nearest, Tcell, chromosome, position, allele.1, allele.2, strand, allele2_is_REF, rsID, Introgressed_from, SNP_type,
  ANCESTRAL, DERIVED, INTROGRESSED, REF, ALT, Vindija.der, Chagyrskaya.der, Altai.der, Denisova.der, NAF_V, NAF_D, MaxPosterior_D_VY, MaxPosterior_V_DY,
  YRI.der, SGDP_African.der, GBR.der, IBS.der, SGDP_WestEurasian.der, CHB.der, JPT.der, SGDP_Eastasian.der, PJL.der, STU.der, SGDP_Agta.der, SGDP_Papuan.der
)]
emVars_signif_active <- merge(emVars_signif_active, condition_summary, by.x = "ANALYSIS_NAME", by.y = "analysis_name")
fwrite(emVars_signif_active, file = sprintf("%s/all_emVars_active_signif.tsv", OUT_DIR), sep = "\t")


emVARS_annot_obs <- emVARs_annot[power == 0 & perm == "perm_0_0_0", ]
emVARS_annot_obs <- merge(emVARS_annot_obs[grepl("celltype_cond", ANALYSIS_SUBTYPE), ], condition_summary, by.x = "ANALYSIS_NAME", by.y = "analysis_name")
Pct_emVars_byActivity <- emVARS_annot_obs[!is.na(oligo_class_loose), .(.N, N_emVARs = sum(FDR < .05)), keyby = factor(CRS_class_loose, c("strong enhancer", "enhancer", "weak enhancer", "inactive", "weak silencer", "silencer", "strong silencer"))]
Pct_emVars_byActivity[, pct_emVars := N_emVARs / N]
Pct_emVars_byActivity[, pct_emVars_lo := binom.test(N_emVARs, N)$conf.int[1], by = factor]
Pct_emVars_byActivity[, pct_emVars_hi := binom.test(N_emVARs, N)$conf.int[2], by = factor]
N_emVARs_byCond <- emVARS_annot_obs[, .(N_emVARs = sum(FDR < .05)), by = COND_ID]

# emVARS_annot_obs[grepl("celltype_cond", ANALYSIS_SUBTYPE) & !is.na(log2FC_archaic_vs_modern),.(.N,sum(FDR<0.05 & abs(log2FC_archaic_vs_modern)>.3),mean(FDR<.05 & abs(log2FC_archaic_vs_modern)>.3)),by=CRS_class_loose]

emVARS_annot_obs <- emVARs_annot_active[power == 0 & perm == "perm_0_0_0", ]
emVARS_annot_obs <- merge(emVARS_annot_obs[grepl("celltype_cond", ANALYSIS_SUBTYPE), ], condition_summary, by.x = "ANALYSIS_NAME", by.y = "analysis_name")
Pct_emVars_active_byActivity <- emVARS_annot_obs[!is.na(CRS_class_loose), .(N_emVARs = sum(FDR < .05)), keyby = factor(CRS_class_loose, oligo_activity_levels)]
Pct_emVars_active_byActivity <- merge(Pct_emVars_byActivity[, .(factor, N, N_emVARs)], Pct_emVars_active_byActivity, by = "factor", all.x = TRUE, suffixes = c("", "_active"))
Pct_emVars_active_byActivity[factor == "inactive", N_emVARs_active := 0]
Pct_emVars_active_byActivity <- melt(Pct_emVars_active_byActivity, id.vars = c("factor", "N"))
Pct_emVars_active_byActivity[, pct_emVars := value / N]
Pct_emVars_active_byActivity[, pct_emVars_lo := binom.test(value, N)$conf.int[1], by = .(factor, variable)]
Pct_emVars_active_byActivity[, pct_emVars_hi := binom.test(value, N)$conf.int[2], by = .(factor, variable)]

N_emVARs_active_byCond <- emVARS_annot_obs[, .(N_emVARs = sum(FDR < .05)), by = COND_ID]

dir.create(sprintf("%s/05_emVARs/", FIGURE_DIR), recursive = TRUE)
col_strategy <- c("N_emVARs" = "#00000088", "N_emVARs_active" = "#DD555588")
p <- ggplot(Pct_emVars_active_byActivity, aes(x = factor, y = pct_emVars, col = variable, ymin = pct_emVars_lo, ymax = pct_emVars_hi))
p <- p + geom_pointrange(size = .2) + theme_plot(rotate.x = 90) + xlab("")
p <- p + scale_color_manual(values = col_strategy)
pdf(sprintf("%s/05_emVARs/01a_active_emVars_by_activity_pct.pdf", FIGURE_DIR), height = 3.5, width = 2.3)
print(p)
dev.off()

p <- ggplot(Pct_emVars_active_byActivity, aes(x = factor, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_plot(rotate.x = 90) +
  xlab("") +
  ylab("N_emVars")
p <- p + scale_fill_manual(values = col_strategy)
pdf(sprintf("%s/05_emVARs/01b_active_emVars_by_activity_counts.pdf", FIGURE_DIR), height = 3.5, width = 2.3)
print(p)
dev.off()

N_emVARs_byCond <- rbindlist(list(all = N_emVARs_byCond, activeOnly = N_emVARs_active_byCond), idcol = "type")
p <- ggplot(N_emVARs_byCond, aes(x = COND_ID, y = N_emVARs, fill = COND_ID, col = type, alpha = type))
p <- p + geom_bar(stat = "identity", position = "dodge", width = .7, stroke = 1.5) + theme_plot(rotate.x = 90) + xlab("")
p <- p + scale_fill_manual(values = color_setup_simplified_norep)
p <- p + scale_alpha_manual(values = c(all = 1, activeOnly = .7))
p <- p + scale_color_manual(values = c(all = "#00000088", activeOnly = "#DD555588"))
p <- p + guides(fill = "none")
pdf(sprintf("%s/05_emVARs/01c_N_emVars_by_condition.pdf", FIGURE_DIR), height = 3, width = 3.5)
print(p)
dev.off()

####################################################################
################# compute FDR, CRS difference ######################
####################################################################

oligo_activity_Diff[, log2FC_der_vs_anc := ifelse(allele.1 == ANCESTRAL, logFC_2vs1 / log(2), -logFC_2vs1 / log(2))]
oligo_activity_Diff[, log2FC_archaic_vs_modern := ifelse(!grepl("ancestral reintrogressed", SNP_type), log2FC_der_vs_anc, -log2FC_der_vs_anc)]
oligo_activity_Diff[, log2FC.se := logFC.se / log(2)]

oligo_activity_Diff[, .(.N, sum(df.test > 0)), keyby = .(power, boot, ANALYSIS_NAME, Z_th, perm)]
oligo_activity_Diff_annot <- merge(oligo_activity_Diff[df.test > 0, ], oligo_source, by = "CRS", allow.cartesian = TRUE)
oligo_activity_Diff_annot <- oligo_activity_Diff_annot[order(pval.LRT), .SD, by = .(Z_th, power, boot, ANALYSIS_NAME)]
oligo_activity_Diff_annot[, N_FP := (.1 + cumsum(perm != "perm_0_0_0" & df.test > 0)), by = .(Z_th, power, boot, ANALYSIS_NAME)]
oligo_activity_Diff_annot[, N_POS := (.1 + cumsum(perm == "perm_0_0_0" & df.test > 0)), by = .(Z_th, power, boot, ANALYSIS_NAME)]
oligo_activity_Diff_annot[, FDR := pmin(1, rev(cummin(rev(N_FP / N_POS)))), by = .(Z_th, power, boot, ANALYSIS_NAME)]

# emVARs_Diff_annot[FDR < .05 & power == 0 & perm == "perm_0_0_0" & grepl("celltype", ANALYSIS_SUBTYPE), .N, by = ANALYSIS_NAME]
fwrite(oligo_activity_Diff_annot[power == 0 & perm == "perm_0_0_0", ], file = sprintf("%s/all_oligos_diff_annotated.tsv.gz", OUT_DIR), sep = "\t")


# fwrite(oligo_activity_DT, file = sprintf("%s/data/%s/MPRA_analyses/%s/allConditions_alphas_and_pvalues_v%s.1_max%sBC_annot.txt.gz", MPRA_DIR, ANALYSIS_DIR, 'oligo_activity', VERSION_OUT, sub_sample), sep = "\t")

# TODO : rely on permutation-based FDR when available
CRS_DiffActivity_DT <- CRS_DiffActivity_DT[order(analysis_name, pval.LRT), ]
CRS_DiffActivity_DT <- CRS_DiffActivity_DT[, FDR := pmin(1, rev(cummin(rev(N_FP / N_POS)))), by = analysis_name]

CRS_DiffActivity_DT[df.test > 0, FDR_analytic := p.adjust(pval.LRT, "fdr"), by = analysis_name]
CRS_DiffActivity_DT[, .N, by = .(analysis_name, FDR < .05, FDR_analytic < 0.05)]
CRS_DiffActivity_DT[, .(.N, N_perm = sum(FDR < .01), N_analytic = sum(FDR_analytic < .01), N_both = sum(FDR < .01 & FDR_analytic < .01)), by = .(analysis_name)]
# trend : unclear, but more significant stuff in permutation-based FDR for cell type comparisons and HepG2-IFN & K562-NS replicates, and less for response_A549_SARS response

CRS_DiffActivity_DT[analysis_name == "HepG2_IFNA2b_replicates" & FDR < .001 & FDR_analytic < 0.0001 & (nBC_g1 + nBC_g2) > 100, ][order(-activity_group1 - activity_group2)][1:3, ]
x <- CRS_DiffActivity_DT[analysis_name == "HepG2_IFNA2b_replicates", ]
# myCRS= CRS_DiffActivity_DT[analysis_name=='HepG2_IFNA2b_replicates' & FDR<.001 & FDR_analytic<0.0001  & (nBC_g1+nBC_g2)>100,][order(-activity_group1-activity_group2)][1,CRS]
# instruction_line=8


# CRS_DiffActivity_DT=merge(CRS_DiffActivity_DT[grepl('all',analysis_name)], condition_summary, by="analysis_name" , all.x=TRUE, allow.cartesian=TRUE)
# CRS_DiffActivity_DT=merge(CRS_DiffActivity_DT, oligo_source[,.(CRS,E071)], by="CRS" , all.x=TRUE, allow.cartesian=TRUE)

CRS_DiffActivity_DT[, ROADMAP_current_celline := case_when(
  grepl("response_A549", analysis_name) ~ E114,
  grepl("A549_.*_replicates", analysis_name) ~ E114,
  grepl("response_HepG2", analysis_name) ~ E118,
  grepl("HepG2_.*_replicates", analysis_name) ~ E118,
  grepl("response_K562", analysis_name) ~ E123,
  grepl("K562_.*_replicates", analysis_name) ~ E123,
  TRUE ~ NA
)]
CRS_DiffActivity_DT[, ROADMAP_ctrl_brain := E071]
CRS_DiffActivity_DT[, celline := case_when(
  grepl("response_A549", analysis_name) ~ "A549",
  grepl("A549_.*_replicates", analysis_name) ~ "A549",
  grepl("response_HepG2", analysis_name) ~ "HepG2",
  grepl("HepG2_.*_replicates", analysis_name) ~ "HepG2",
  grepl("response_K562", analysis_name) ~ "K562",
  grepl("K562_.*_replicates", analysis_name) ~ "K562",
  TRUE ~ NA
)]

fwrite(CRS_DiffActivity_DT, file = sprintf("%s/data/%s/MPRA_analyses/%s/allInteractions_alphas_and_pvalues_v%s.1_max%sBC_annot.txt.gz", MPRA_DIR, ANALYSIS_DIR, "oligo_activity", VERSION_IN, sub_sample), sep = "\t")
# CRS_DiffActivity_DT[,mean(FDR<.05),by=ROADMAP_current_celline]


# emVARs_Diff_annot[power == 0 & perm == "perm_0_0_0" & grepl("replicate", ANALYSIS_SUBTYPE), .(.N, sum(FDR < .1)), keyby = .(ANALYSIS_NAME)]
# oligo_activity_Diff[power == 0 & perm == "perm_0_0_0" & grepl("replicate_all", ANALYSIS_SUBTYPE), .(.N, sum(pval.LRT < .001)), keyby = .(ANALYSIS_NAME)]
# oligo_activity_Diff[power == 0 & perm == "perm_0_1_0" & grepl("replicate_all", ANALYSIS_SUBTYPE), .(.N, sum(pval.LRT < .001)), keyby = .(ANALYSIS_NAME)]

cat("All done!\n")
q("no")
# #### aggregate power results
# instruction_DT <- fread(instruction_file, header = TRUE, fill = TRUE)
# for (instruction_line in instruction_DT[,which(analysis_type=="emVars")]) {
#   ANALYSIS_NAME <- instruction_DT[instruction_line, analysis_name]
#   ANALYSIS_TYPE <- instruction_DT[instruction_line, analysis_type]
#   if (grepl("EMVARS", toupper(ANALYSIS_TYPE))) {
#     OUTS_DIR <- sprintf("%s/outs/%s/%s/%s", MPRA_DIR, ANALYSIS_DIR, ANALYSIS_TYPE, ANALYSIS_NAME)
#     results_MPRA <- list()
#     for (perm in 0:1) {
#       for (boot in c(1,1.2,1.5,2)) {
#         IN_DIR <- sprintf("%s/outs/%s/%s/%s/perm%s/power1_boot%s/", MPRA_DIR, ANALYSIS_DIR, ANALYSIS_TYPE, ANALYSIS_NAME, perm, boot)
#         results_MPRA[[paste0("perm", perm,'_boot',boot)]] <- list()
#         for (crs_num in seq_len(instruction_DT[instruction_line, N_test])) {
#           results_MPRA[[paste0("perm", perm,'_boot',boot)]][[crs_num]] <- try(fread(file = sprintf("%s/alphas_and_pvalues_v%s_SNP%s_max%sBC.txt", IN_DIR, VERSION_IN, crs_num, sub_sample), sep = "\t"), silent = TRUE)
#           if (class(results_MPRA[[paste0("perm", perm,'_boot',boot)]][[crs_num]]) == "try-error") {
#             results_MPRA[[paste0("perm", perm,'_boot',boot)]][[crs_num]] <- NULL
#           }
#         }
#         results_MPRA[[paste0("perm", perm,'_boot',boot)]] <- rbindlist(results_MPRA[[paste0("perm", perm,'_boot',boot)]], idcol = "CRS_num", fill = TRUE)
#       }
#     }
#     results_MPRA <- rbindlist(results_MPRA, idcol = "perm_boot")
#     results_MPRA[,perm:=gsub('perm(.*)_boot(.*)', '\\1', perm_boot)]
#     results_MPRA[,boot:=gsub('perm(.*)_boot(.*)', '\\2', perm_boot)]
#     results_MPRA[,perm_boot:=NULL]
#     results_MPRA[, analysis_name := ANALYSIS_NAME]
#     results_MPRA[, analysis_type := ANALYSIS_TYPE]
#     fwrite(results_MPRA, file = sprintf("%s/alphas_and_pvalues_POWER_v%s_allSNPs_max%sBC.txt.gz", OUTS_DIR, VERSION_OUT, sub_sample))
#     cat(ANALYSIS_TYPE, ANALYSIS_NAME, ":", results_MPRA[, .N], "\n")
#   }
# }



#### aggregate stats
for (instruction_line in seq_len(instruction_DT[, .N])) {
  ANALYSIS_NAME <- instruction_DT[instruction_line, analysis_name]
  ANALYSIS_TYPE <- instruction_DT[instruction_line, analysis_type]

  OUTS_DIR <- sprintf("%s/outs/%s/%s/%s/", MPRA_DIR, ANALYSIS_DIR, ANALYSIS_TYPE, ANALYSIS_NAME)
  if (grepl("OLIGO", to_upper(ANALYSIS_TYPE))) {
    OUTS_DIR <- sprintf("%s/outs/%s/%s/%s/", MPRA_DIR, ANALYSIS_DIR, ANALYSIS_TYPE, ANALYSIS_NAME)
    CRS_stats <- list()
    for (perm in 0:1) {
      IN_DIR <- sprintf("%s/outs/%s/%s/%s/perm%s", MPRA_DIR, ANALYSIS_DIR, ANALYSIS_TYPE, ANALYSIS_NAME, perm)
      CRS_stats[[paste0("perm", perm)]] <- list()
      for (crs_num in seq_len(instruction_DT[instruction_line, N_test])) {
        CRS_stats[[paste0("perm", perm)]][[crs_num]] <- try(fread(file = sprintf("%s/stats_v%s_CRS%s_max%sBC.txt", IN_DIR, VERSION_IN, crs_num, sub_sample), sep = "\t"), silent = TRUE)
        if (class(CRS_stats[[paste0("perm", perm)]][[crs_num]]) == "try-error") {
          CRS_stats[[paste0("perm", perm)]][[crs_num]] <- NULL
        }
      }
      CRS_stats[[paste0("perm", perm)]] <- rbindlist(CRS_stats[[paste0("perm", perm)]], idcol = "CRS_num", fill = TRUE)
    }
    CRS_stats <- rbindlist(CRS_stats, idcol = "perm")
    CRS_stats[, analysis_name := ANALYSIS_NAME]
    CRS_stats[, analysis_type := ANALYSIS_TYPE]
    fwrite(CRS_stats, file = sprintf("%s/stats_v%s_allCRS_max%sBC.txt.gz", OUTS_DIR, VERSION_OUT, sub_sample))
    cat(ANALYSIS_TYPE, ANALYSIS_NAME, ":", CRS_stats[, .N], "\n")
  }
  if (grepl("EMVARS", toupper(ANALYSIS_TYPE))) {
    OUTS_DIR <- sprintf("%s/outs/%s/%s/%s", MPRA_DIR, ANALYSIS_DIR, ANALYSIS_TYPE, ANALYSIS_NAME)
    SNP_stats <- list()
    for (perm in 0:1) {
      IN_DIR <- sprintf("%s/outs/%s/%s/%s/perm%s", MPRA_DIR, ANALYSIS_DIR, ANALYSIS_TYPE, ANALYSIS_NAME, perm)
      SNP_stats[[paste0("perm", perm)]] <- list()
      for (crs_num in seq_len(instruction_DT[instruction_line, N_test])) {
        SNP_stats[[paste0("perm", perm)]][[crs_num]] <- try(fread(file = sprintf("%s/stats_v%s_SNP%s_max%sBC.txt", IN_DIR, VERSION_IN, crs_num, sub_sample), sep = "\t"), silent = TRUE)
        if (class(SNP_stats[[paste0("perm", perm)]][[crs_num]]) == "try-error") {
          SNP_stats[[paste0("perm", perm)]][[crs_num]] <- NULL
        }
      }
      SNP_stats[[paste0("perm", perm)]] <- rbindlist(SNP_stats[[paste0("perm", perm)]], idcol = "CRS_num", fill = TRUE)
    }
    SNP_stats <- rbindlist(SNP_stats, idcol = "perm")
    SNP_stats[, analysis_name := ANALYSIS_NAME]
    SNP_stats[, analysis_type := ANALYSIS_TYPE]
    fwrite(SNP_stats, file = sprintf("%s/stats_v%s_allSNPs_max%sBC.txt.gz", OUTS_DIR, VERSION_OUT, sub_sample))
    cat(ANALYSIS_TYPE, ANALYSIS_NAME, ":", SNP_stats[, .N], "\n")
  }
}
cat("All done\n")
q("no")
