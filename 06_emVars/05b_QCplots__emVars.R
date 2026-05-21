# running: sbatch -p geh,common --mem=30G  00_Rscript.sh MPRA_count_exp6_analysisZ/03b_aggMPRAnalyse_oligo_activity.R

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
  if (cmd[i] == "--filter_active" || cmd[i] == "-f") {
    FILTER_ACTIVE <- as.logical(cmd[i + 1])
  }
  if (cmd[i] == "--criteria_emvar" || cmd[i] == "-e") {
    CRITERIA_EMVARS <- cmd[i + 1]
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
SNP_annot <- SNP_annot_v4
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
#FIGURE_DIR <- sprintf("%s/figures/%s/04_emVars/%s/%s", MPRA_DIR, ANALYSIS_DIR, CRITERIA_EMVARS, CRITERIA_ACTIVE_OUT)

dir.create(FIGURE_DIR, recursive = TRUE)
dir.create(sprintf("%s/SupTables/", FIGURE_DIR), recursive = TRUE)
cat("\nprinting output in:", FIGURE_DIR)


# load tested oligos
source(sprintf("%s/scripts/%s/02z__define_includedCRS.R", MPRA_DIR, ANALYSIS_DIR))

library(ggrepel)

########## load emVar data  ##########

emVARs_annot_obs <- fread(file = sprintf("%s/all_emVars_annotated__%s.tsv.gz", EMVAR_DIR, CRITERIA_ACTIVE_OUT), sep = "\t")
emVARs_annot_obs_ctc <- emVARs_annot_obs[ANALYSIS_SUBTYPE == "celltype_cond"]
emVARs_annot_obs_ctc <- merge(emVARs_annot_obs_ctc, condition_summary[, -"celltype"], by.x = "ANALYSIS_NAME", by.y = "analysis_name")

emVARs_annot_perm <- fread(file = sprintf("%s/all_emVars_annotated_perm__%s.tsv.gz", EMVAR_DIR, CRITERIA_ACTIVE_OUT), sep = "\t")
emVARs_annot_perm_ctc <- emVARs_annot_perm[ANALYSIS_SUBTYPE == "celltype_cond"]
emVARs_annot_perm_ctc <- merge(emVARs_annot_perm_ctc, condition_summary[, -"celltype"], by.x = "ANALYSIS_NAME", by.y = "analysis_name")

emVARs_annot_obs_sample <- emVARs_annot_obs[ANALYSIS_SUBTYPE == "sample"]
emVARs_annot_obs_sample <- merge(emVARs_annot_obs_sample, condition_summary_reps[, -"celltype"], by.x = "ANALYSIS_NAME", by.y = "analysis_name")

emVARs_obs_ctc_full <- fread(sprintf("%s/all_emVars_annotated_full_celltype__%s.tsv", EMVAR_DIR, CRITERIA_ACTIVE_OUT))

N_emVARs_byCond <- emVARs_annot_obs_ctc[, .(N_emVARs_detected = sum(is_emVar_CRITERIA), N_emVARs_expected = sum(1 - lfdr)), by = COND_ID]
N_emVARs_byCond

oligo_target_collapsed <- fread(sprintf("%s/data/%s/00_oligo_targets_collapsed_v4.txt", MPRA_DIR, ANALYSIS_DIR), sep = "\t")
oligo_target <- fread(sprintf("%s/data/%s/00_oligo_targets_v4.txt", MPRA_DIR, ANALYSIS_DIR))
crs_targets <- unique(oligo_target_collapsed[, -"oligo"])
NearestGene <- unique(oligo_target[method == "NearestGene" & gene_type != "", .(posID = crsID, TargetGene)])
locus_genes <- merge(SNP_annot[,.(posID, introgression_locus )],NearestGene, by="posID")
locus_genes <- locus_genes[introgression_locus!="",.(locus_genes=paste(unique(TargetGene),collapse='/')),by=introgression_locus]


#          COND_ID N_emVARs_detected N_emVARs_expected
#  1:     A549_IAV               110               255
#  2:      A549_NS               132               381
#  3:    A549_SARS                81               285
#  4:    A549_TNFa               137               279
#  5:    HepG2_DEX               186               373
#  6: HepG2_IFNA2b               149               461
#  7:     HepG2_NS               159               462
#  8:   HepG2_TNFa               103               286
#  9:     K562_DEX                62               145
# 10:  K562_IFNA2b                65               231
# 11:      K562_NS                46               228
# 12:    K562_TNFa                61               219

# IN_DIR <- sprintf("%s/data/%s/03a_aggMPRA_analyse/%s/", MPRA_DIR, ANALYSIS_DIR, RUN_ID)
# OUT_DIR <- sprintf("%s/data/%s/03b_computeFDR/%s/", MPRA_DIR, ANALYSIS_DIR, RUN_ID)
# dir.create(OUT_DIR, recursive = TRUE)

# emVARs_all_results <- fread(sprintf("%s/all_emVars_results.tsv.gz", IN_DIR), sep = "\t")
# emVARs_all_results[, crsID := posID]

# CRS_activity_file <- sprintf("%s/CRS_activity__celltype_cond.tsv.gz", OUT_DIR)
# CRS_activity <- fread(CRS_activity_file)

# emVARs_active_obs <- fread(sprintf("%s/all_emVars_activeCRS_annotated.tsv.gz", IN_DIR), sep = "\t")

#######################################################################################################################
#######################################################################################################################
########################################## QC PLOTS, EMVARS ###########################################################
#######################################################################################################################
#######################################################################################################################


####################################################################################################
########################### pvalue distribution emvars ############################################
####################################################################################################
# celltype data
FigData <- list(observed = emVARs_annot_obs_ctc[condition == "NS", ], permuted = emVARs_annot_perm_ctc[condition == "NS", ])
FigData <- rbindlist(FigData, idcol = "data_source")

p <- ggplot(FigData, aes(x = pval.LRT, fill = celline, alpha = data_source)) +
  geom_histogram(position = "dodge")
p <- p + facet_grid(cols = vars(celline)) + guides(fill = "none")
p <- p + scale_fill_manual(values = color_celline)
p <- p + scale_alpha_manual(values = c("observed" = .4, "permuted" = .7))
p <- p + theme_plot(lpos = "top", fontsize = 12, rotate.x = 90) + xlab("Differential activity p-values (emVar)")

pdf(sprintf("%s/1a__emVar_pvalues__observed_vs_permuted.pdf", FIGURE_DIR), height = 2.5, width = 4)
print(p)
dev.off()


# celltype data
FigData <- list(observed = emVARs_annot_obs_ctc[condition == "NS", ], permuted = emVARs_annot_perm_ctc[condition == "NS", ])
FigData <- rbindlist(FigData, idcol = "data_source")

p <- ggplot(FigData, aes(x = pval_emp, fill = celline, alpha = data_source)) +
  geom_histogram(position = "dodge")
p <- p + facet_grid(cols = vars(celline)) + guides(fill = "none")
p <- p + scale_fill_manual(values = color_celline)
p <- p + scale_alpha_manual(values = c("observed" = .4, "permuted" = .7))
p <- p + theme_plot(lpos = "top", fontsize = 12, rotate.x = 90) + xlab("Differential activity p-values (emVar)")

pdf(sprintf("%s/1e__emVar_empirical_pvalues__observed_vs_permuted.pdf", FIGURE_DIR), height = 2.5, width = 4)
print(p)
dev.off()

FigData <- list(observed = emVARs_annot_obs_ctc[condition == "NS", ], permuted = emVARs_annot_perm_ctc[condition == "NS", ])
FigData <- rbindlist(FigData, idcol = "data_source")

p <- ggplot(FigData, aes(x = pval_emp, fill = celline, alpha = data_source)) +
  geom_histogram(position = "dodge", breaks=seq(0, 1, l=31))
p <- p + facet_grid(cols = vars(celline)) + guides(fill = "none")
p <- p + scale_fill_manual(values = color_celline)
p <- p + scale_alpha_manual(values = c("observed" = .4, "permuted" = .7))
p <- p + theme_plot(lpos = "top", fontsize = 12, rotate.x = 90) + xlab("Differential activity p-values (emVar, empirical)")

pdf(sprintf("%s/1a__emVar_pval_emp__observed_vs_permuted.pdf", FIGURE_DIR), height = 2.5, width = 4)
print(p)
dev.off()


p <- ggplot(FigData, aes(x = pval.LRT, fill = celline, alpha = data_source)) +
  geom_histogram()
p <- p + facet_grid(rows = vars(celline), cols = vars(data_source)) + guides(fill = "none", alpha = "none")
p <- p + scale_fill_manual(values = color_celline)
p <- p + scale_alpha_manual(values = c("observed" = .4, "permuted" = .7))
p <- p + theme_plot(lpos = "top", fontsize = 12, rotate.x = 90) + xlab("Differential activity p-values (emVar)")

pdf(sprintf("%s/1b__emVar_pvalues__observed_vs_permuted_split.pdf", FIGURE_DIR), height = 3, width = 3.5)
print(p)
dev.off()


p <- ggplot(FigData, aes(x = pmax(-4, pmin(4, -qnorm(pval.LRT) * sign(log2FC_archaic_vs_modern))), fill = celline, alpha = data_source)) +
  geom_histogram(position = "dodge")
p <- p + facet_grid(cols = vars(celline)) + guides(fill = "none")
p <- p + scale_fill_manual(values = color_celline)
p <- p + scale_alpha_manual(values = c("observed" = .4, "permuted" = .7))
p <- p + theme_plot(lpos = "top", fontsize = 12, rotate.x = 90) + xlab("Differential activity Z score (emVar)")

pdf(sprintf("%s/1c__emVar_Zscore__observed_vs_permuted.pdf", FIGURE_DIR), height = 2.5, width = 4)
print(p)
dev.off()

p <- ggplot(FigData, aes(x = pmax(-4, pmin(4, -qnorm(pval.LRT) * sign(log2FC_archaic_vs_modern))), fill = celline, alpha = data_source)) +
  geom_histogram()
p <- p + facet_grid(rows = vars(celline), cols = vars(data_source)) + guides(fill = "none", alpha = "none")
p <- p + scale_fill_manual(values = color_celline)
p <- p + scale_alpha_manual(values = c("observed" = .4, "permuted" = .7))
p <- p + theme_plot(lpos = "top", fontsize = 12, rotate.x = 90) + xlab("Differential activity Z score (emVar)")

pdf(sprintf("%s/1d__emVar_Zscore__observed_vs_permuted_split.pdf", FIGURE_DIR), height = 3, width = 3.5)
print(p)
dev.off()


# celltype data
# FigData <- list(observed = emVARs_annot_obs_ctc[condition == "NS", ], permuted = emVARs_annot_perm_ctc[condition == "NS", ])
# FigData <- rbindlist(FigData, idcol = "data_source")

# p <- ggplot(FigData, aes(x = pval_emp, fill = celline, alpha = data_source)) +
#   geom_histogram(position = "dodge")
# p <- p + facet_grid(cols = vars(celline)) + guides(fill = "none")
# p <- p + scale_fill_manual(values = color_celline)
# p <- p + scale_alpha_manual(values = c("observed" = .4, "permuted" = .7))
# p <- p + theme_plot(lpos = "top", fontsize = 12, rotate.x = 90) + xlab("Differential activity p-values (emVar)")

# pdf(sprintf("%s/1e__emVar_empirical_pvalues__observed_vs_permuted.pdf", FIGURE_DIR), height = 2.5, width = 4)
# print(p)
# dev.off()


####################################################################################################
########################### volcano plot emvars ####################################################
####################################################################################################

crs_target <- fread(sprintf("%s/data/%s/00_oligo_targets_v4.txt", MPRA_DIR, ANALYSIS_DIR))
# celltype data
# FigData <- merge(emVARs_annot_obs_ctc,unique(crs_target[method=='NearestTSS' & celltype!='any',.(crsID,TargetGene,gene_type,celline=celltype)]),by=c('celline','crsID'))
FigData <- merge(emVARs_annot_obs_ctc, unique(crs_target[method == "NearestTSS" & celltype == "any", .(crsID, TargetGene, gene_type)]), by = c("crsID"))

nPlot <- 10
# FigData_annot <- FigData[is_emVar_CRITERIA==TRUE & gene_type!="",][order(rank(-abs(log2FC_archaic_vs_modern))+rank(FDR)),head(.SD,nPlot),by=COND_ID]
FigData_annot <- FigData[is_emVar_CRITERIA == TRUE, ][order(pmin(rank(-abs(log2FC_archaic_vs_modern)), rank(pval.LRT))), head(.SD, nPlot), by = COND_ID]

  p <- ggplot(FigData[condition == "NS", ], aes(x = log2FC_archaic_vs_modern, y = -log10(pval.LRT)))
  p <- p + rasterize(geom_point(size = 0.1, alpha = 0.4, aes(col = COND_ID)), dpi = 200)
  p <- p + facet_grid(cols = vars(celline)) + guides(fill = "none")
  p <- p + scale_color_manual(values = color_setup_simplified_norep)
  p <- p + theme_plot(lpos = "top", fontsize = 12, rotate.x = 90)
  p <- p + xlab("Log2FC (Introgressed vs non-introgressed)") + ylab("-log10(p-value)")
  p <- p + geom_text_repel(data = FigData_annot[condition == "NS", ], aes(x = log2FC_archaic_vs_modern, y = -log10(pval.LRT), label = TargetGene), size = 2, segment.size = 0.1)
  p <- p + guides(color = "none", fill = "none", alpha = "none")

pdf(sprintf("%s/2a__emVar_volcano_NS_conditions.pdf", FIGURE_DIR), height = 2.5, width = 4)
print(p)
dev.off()
p <- ggplot(FigData[condition == "NS", ], aes(x = log2FC_archaic_vs_modern, y = -log10(pval.LRT)))
  p <- p + rasterize(geom_point(size = 0.1, alpha = 0.4, aes(col = COND_ID)), dpi = 200)
  p <- p + facet_grid(cols = vars(celline)) + guides(fill = "none")
  p <- p + scale_color_manual(values = color_setup_simplified_norep)
  p <- p + theme_plot(lpos = "top", fontsize = 12, rotate.x = 90)
  p <- p + xlab("Log2FC (Introgressed vs non-introgressed)") + ylab("-log10(p-value)")
  p <- p + geom_text_repel(data = FigData_annot[condition == "NS", ], aes(x = log2FC_archaic_vs_modern, y = -log10(pval.LRT), label = TargetGene), size = 2, segment.size = 0.1)
  p <- p + guides(color = "none", fill = "none", alpha = "none")

pdf(sprintf("%s/2b__emVar_volcano_all_conditions.pdf", FIGURE_DIR), height = 5, width = 4)
print(p)
dev.off()

FigData_no_dup <- FigData[order(!is_emVar_CRITERIA, -abs(log2FC_archaic_vs_modern)),head(.SD,1),by=.(ID,celline)]
FigData_annot_no_dup <- FigData_no_dup[is_emVar_CRITERIA == TRUE, ][order(pmin(rank(-abs(log2FC_archaic_vs_modern)), rank(pval.LRT))),][,head(.SD, nPlot), by = celline]

p <- ggplot(FigData_no_dup, aes(x = log2FC_archaic_vs_modern, y = -log10(pval.LRT)))
p <- p + rasterize(geom_point(size = 0.1, alpha = 0.4, aes(col = COND_ID)), dpi = 200)
p <- p + facet_grid(cols = vars(celline)) + guides(fill = "none")
p <- p + scale_color_manual(values = color_setup_simplified_norep)
p <- p + theme_plot(lpos = "top", fontsize = 12, rotate.x = 90)
p <- p + xlab("Log2FC (Introgressed vs non-introgressed)") + ylab("-log10(p-value)")
p <- p + geom_text_repel(data = FigData_annot_no_dup, aes(x = log2FC_archaic_vs_modern, y = -log10(pval.LRT), label = TargetGene), size = 2, segment.size = 0.1)
p <- p + guides(color = "none", fill = "none", alpha = "none")

pdf(sprintf("%s/2c__emVar_volcano_all_conditions_1volcano_per_cellline.pdf", FIGURE_DIR), height = 2.5, width = 4)
print(p)
dev.off()

FigData_no_dup[is_emVar_CRITERIA == TRUE & abs(log2FC_archaic_vs_modern) > 0.5, .N, by=POP_adaptive_top]
########################################################################################################################
########################################## Pct emVars across replicates ################################################
########################################################################################################################

nEmVars <- emVARs_annot_obs[, .(nEmVars = sum(is_emVar_CRITERIA), .N), keyby = .(ANALYSIS_SUBTYPE, ANALYSIS_NAME)]
nEmVars[, Pct := binom.test(nEmVars, N)$estimate * 100, keyby = .(ANALYSIS_SUBTYPE, ANALYSIS_NAME)]
nEmVars[, Pct_lo := binom.test(nEmVars, N)$conf.int[1] * 100, keyby = .(ANALYSIS_SUBTYPE, ANALYSIS_NAME)]
nEmVars[, Pct_hi := binom.test(nEmVars, N)$conf.int[2] * 100, keyby = .(ANALYSIS_SUBTYPE, ANALYSIS_NAME)]
nEmVars <- nEmVars[ANALYSIS_NAME != "A549_NS_all" & ANALYSIS_SUBTYPE != "celltype", ]
nEmVars[, ANALYSIS_NAME := gsub("_all", "", ANALYSIS_NAME)]

p <- ggplot(nEmVars, aes(x = ANALYSIS_NAME, y = Pct, ymin = Pct_lo, ymax = Pct_hi, col = ANALYSIS_NAME)) +
  geom_pointrange(size = .1)
p <- p + scale_color_manual(values = c(color_celltype, color_setup, color_setup_norep)) + theme_plot(rotate.x = 90, fontsize = 12)
p <- p + ylab("Pct significant emVars") + facet_grid(cols = vars(ANALYSIS_SUBTYPE), scale = "free_x", space = "free_x")
p <- p + guides(col = "none") + ylim(c(0, 35))
pdf(sprintf("%s/3a_Pct_EmVars_per_Replicate.pdf", FIGURE_DIR), height = 4, width = 7)
print(p)
dev.off()

########################################################################################################################
######################################## plot Pct emVars vs CRS activity ###############################################
########################################################################################################################

# Sup Table 3a
####### extract significant emVars at 5% FDR
# emVars_signif_active <- emVARs_annot_full[is_emVar_CRITERIA == TRUE, .(posID, ANALYSIS_NAME, oligo1, oligo2, GC, pval.LRT, log2FC_archaic_vs_modern, log2FC.se,
#   FDR_emVAR = FDR, cre_class_CRITERIA, Lung, Mono, Nearest, Tcell, chromosome, position, allele.1, allele.2, strand, allele2_is_REF, rsID, Adaptive_from, Introgression_scenario=Introgression_scenario_v2,
#   ANCESTRAL, DERIVED, INTROGRESSED = INTROGRESSED.allele, REF, ALT, Vindija.der, Chagyrskaya.der, Altai.der, Denisova.der, NAF_V, NAF_D, MaxPosterior_D_VY, MaxPosterior_V_DY,
#   YRI.der, SGDP_African.der, GBR.der, IBS.der, SGDP_WestEurasian.der, CHB.der, JPT.der, SGDP_Eastasian.der, PJL.der, STU.der, SGDP_Agta.der, SGDP_Papuan.der
# )]
# emVars_signif_active <- merge(emVars_signif_active, condition_summary, by.x = "ANALYSIS_NAME", by.y = "analysis_name")
# fwrite(emVars_signif_active, file = sprintf("%s/all_emVars_active_signif.tsv", EMVAR_DIR), sep = "\t")

# TableS3A <- emVars_signif_active[, .(rsID , chromosome, position,   ANCESTRAL, DERIVED, INTROGRESSED = INTROGRESSED.allele, CRE_class=cre_class_CRITERIA, celline, condition, log2FC_archaic_vs_modern, log2FC.se, pval.LRT,  lfdr_emVar=lfdr, FDR_emVAR = FDR)]
# fwrite(TableS3A, file = sprintf("%s/SupTable/TableS3A_emVars.tsv", EMVAR_DIR), sep = "\t")




# subset to cres that are active in at least one condition
oligo_activity_obs <- fread(sprintf("%s/all_oligos_annotated__%s.tsv.gz", ACTIVE_DIR, CRITERIA_ACTIVE))
emVARs_obs_ctc <- merge(emVARs_annot_obs_ctc,
  oligo_activity_obs[ANALYSIS_SUBTYPE == "celltype_cond", .(log2_CRE_activity_CRITERIA = mean(log2(alpha_CRITERIA))), by = .(posID = crsID, ANALYSIS_NAME)],
  by = c("posID", "ANALYSIS_NAME"),
  all.x = TRUE,
  suffix = c(".emVAR", ".oligo")
)


Pct_emVars_byActivity <- emVARs_obs_ctc[cre_class_CRITERIA != '', .(.N, N_emVARs = sum(is_emVar_CRITERIA)), keyby = factor(cre_class_CRITERIA, c("strong silencer", "silencer", "inactive",  "enhancer", "strong enhancer"), labels = rev(c("strong enhancer", "weak enhancer", "inactive", "weak silencer", "strong silencer")))]
Pct_emVars_byActivity[, pct_emVars := N_emVARs / N * 100]
Pct_emVars_byActivity[, pct_emVars_lo := binom.test(N_emVARs, N)$conf.int[1] * 100, by = factor]
Pct_emVars_byActivity[, pct_emVars_hi := binom.test(N_emVARs, N)$conf.int[2] * 100, by = factor]

p <- ggplot(Pct_emVars_byActivity, aes(x = factor, y = pct_emVars, ymin = pct_emVars_lo, ymax = pct_emVars_hi, col = factor))
p <- p + geom_pointrange(size = .1)
p <- p + scale_color_manual(values = c(CRE_class_colors)) + theme_plot(rotate.x = 90, fontsize = 12)
p <- p + ylab("Pct significant emVars") + xlab("CRE activity") # + facet_grid(cols = vars(ANALYSIS_SUBTYPE), scale = "free_x", space = "free_x")
p <- p + guides(col = "none") + ylim(c(0, 35))

pdf(sprintf("%s/4b_Pct_EmVars_per_CRE_class.pdf", FIGURE_DIR), height = 4, width = 1.7)
print(p)
dev.off()


CUTS <- c(-Inf, -1, -0.5, -0.3, -0.2, 0.2, 0.3, 0.5, 1, Inf)
Pct_emVars_byActivity <- emVARs_obs_ctc[, .(.N, N_emVARs = sum(is_emVar_CRITERIA)), keyby = cut(log2_CRE_activity_CRITERIA, CUTS)]
Pct_emVars_byActivity[, pct_emVars := N_emVARs / N * 100]
Pct_emVars_byActivity[, pct_emVars_lo := binom.test(N_emVARs, N)$conf.int[1] * 100, by = cut]
Pct_emVars_byActivity[, pct_emVars_hi := binom.test(N_emVARs, N)$conf.int[2] * 100, by = cut]
Pct_emVars_byActivity <- Pct_emVars_byActivity[!is.na(cut)]
cut_colors <- setNames(colorRampPalette(rev(CRE_class_colors))(length(CUTS) - 1), levels(Pct_emVars_byActivity[, cut]))

p <- ggplot(Pct_emVars_byActivity, aes(x = cut, y = pct_emVars, ymin = pct_emVars_lo, ymax = pct_emVars_hi, col = cut))
p <- p + geom_pointrange(size = .1)
p <- p + theme_plot(rotate.x = 90, fontsize = 12) + scale_color_manual(values = c(cut_colors))
p <- p + ylab("Pct significant emVars") # + facet_grid(cols = vars(ANALYSIS_SUBTYPE), scale = "free_x", space = "free_x")
p <- p + guides(col = "none") + ylim(c(0, 35)) + xlab("mean CRE activity\n(log2 scale)")

pdf(sprintf("%s/4d_Pct_EmVars_per_Activity_level.pdf", FIGURE_DIR), height = 4, width = 2.5)
print(p)
dev.off()


########################################################################################################################
########################### comparison of positive and negative ctrl emvars ############################################
########################################################################################################################

# TODO



# emVars_annot_obs[grepl("celltype_cond", ANALYSIS_SUBTYPE) & !is.na(log2FC_archaic_vs_modern),.(.N,sum(FDR<0.05 & abs(log2FC_archaic_vs_modern)>.3),mean(FDR<.05 & abs(log2FC_archaic_vs_modern)>.3)),by=oligo_class_loose]
# emVars_annot_obs <- emVARs_annot_active[power == 0 & perm == "perm_0_0_0", ]
# emVars_annot_obs <- merge(emVars_annot_obs[grepl("celltype_cond", ANALYSIS_SUBTYPE), ], condition_summary, by.x = "ANALYSIS_NAME", by.y = "analysis_name")
# Pct_emVars_active_byActivity <- emVars_annot_obs[!is.na(oligo_class_loose), .(N_emVARs = sum(FDR < .05)), keyby = factor(oligo_class_loose, oligo_activity_levels)]
# Pct_emVars_active_byActivity <- merge(Pct_emVars_byActivity[, .(factor, N, N_emVARs)], Pct_emVars_active_byActivity, by = "factor", all.x = TRUE, suffixes = c("", "_active"))
# Pct_emVars_active_byActivity[factor == "inactive", N_emVARs_active := 0]
# Pct_emVars_active_byActivity <- melt(Pct_emVars_active_byActivity, id.vars = c("factor", "N"))
# Pct_emVars_active_byActivity[, pct_emVars := value / N]
# Pct_emVars_active_byActivity[, pct_emVars_lo := binom.test(value, N)$conf.int[1], by = .(factor, variable)]
# Pct_emVars_active_byActivity[, pct_emVars_hi := binom.test(value, N)$conf.int[2], by = .(factor, variable)]

# N_emVARs_active_byCond <- emVars_annot_obs[, .(N_emVARs = sum(FDR < .05)), by = COND_ID]

#############################################################################
################ log2FC distributions by celltype ###########################
#############################################################################

figdata <- emVARs_obs_ctc[condition == "NS"]
p <- ggplot(figdata, aes(x = pmax(-1, pmin(1, log2FC_archaic_vs_modern)), fill = celline, alpha = ifelse(is_emVar_CRITERIA, "sig", "ns")))
p <- p + geom_histogram(col = "white", binwidth = 0.05, linewidth = .1) + facet_grid(rows = vars(celline))
p <- p + theme_plot(fontsize = 12) + scale_fill_manual(values = rev(color_celline)) + scale_alpha_manual(values = c("ns" = 0.5, "sig" = 1))
p <- p + xlab("log2FC introgressed/modern)") + guides(fill = guide_legend(ncol = 3, byrow = FALSE), alpha = FALSE)

pdf(sprintf("%s/5a_log2FC_by_celltype__NSonly__%s.pdf", FIGURE_DIR, CRITERIA_EMVARS), height = 4, width = 3)
print(p)
dev.off()


figdata <- emVARs_obs_ctc
p <- ggplot(figdata, aes(x = pmax(-1, pmin(1, log2FC_archaic_vs_modern)), fill = factor(COND_ID, names(color_setup_simplified_norep)), alpha = ifelse(is_emVar_CRITERIA, "sig", "ns")))
p <- p + geom_histogram(col = "white", binwidth = 0.05, linewidth = .1) + facet_wrap(~ factor(COND_ID, names(color_setup_simplified_norep)), nrow = 3)
p <- p + theme_plot(fontsize = 12) + scale_fill_manual(values = color_setup_simplified_norep) + scale_alpha_manual(values = c("ns" = 0.5, "sig" = 1))
p <- p + xlab("log2FC (introgressed/modern)") + guides(fill = guide_legend(ncol = 4, byrow = TRUE), alpha = FALSE)
pdf(sprintf("%s/5b_log2FC_by_celltype__all__%s.pdf", FIGURE_DIR, CRITERIA_EMVARS), height = 4, width = 6)
print(p)
dev.off()


figdata <- emVARs_obs_ctc
p <- ggplot(figdata, aes(x = pmax(-.5, pmin(.5, log2FC_archaic_vs_modern)), fill = factor(COND_ID, names(color_setup_simplified_norep)), alpha = ifelse(pval_emp<0.05, "sig", "ns")))
p <- p + geom_histogram(col = "white", binwidth = 0.05, linewidth = .1) + facet_wrap(~ factor(COND_ID, names(color_setup_simplified_norep)), nrow = 3)
p <- p + theme_plot(fontsize = 10) + scale_fill_manual(values = color_setup_simplified_norep) + scale_alpha_manual(values = c("ns" = 0.5, "sig" = 1))
p <- p + xlab("log2FC (introgressed/modern)") + guides(fill = FALSE)
pdf(sprintf("%s/5c_log2FC_by_celltype__all__%s.pdf", FIGURE_DIR, CRITERIA_EMVARS), height = 4, width = 4)
print(p)
dev.off()

#############################################################################
################ comparison with Siraj's results ###########################
#############################################################################

Siraj <- fread(sprintf("%s/data/Siraj_2024/TableS3_emVars.csv", MPRA_DIR), skip = 1, header = T)
Siraj_clean <- Siraj[, .(ID = paste(gsub("chr", "", variant), "b37", sep = ":"), active_A549, log2FC_A549, emVar_A549, log2Skew_A549, active_HEPG2, log2FC_HEPG2, emVar_HEPG2, log2Skew_HEPG2, active_K562, log2FC_K562, emVar_K562, log2Skew_K562, mean_Plasmid_K562 = mean_Plasmid_alt_K562 + mean_Plasmid_alt_K562, mean_Plasmid_HEPG2 = mean_Plasmid_alt_HEPG2 + mean_Plasmid_alt_HEPG2, mean_Plasmid_A549 = mean_Plasmid_alt_A549 + mean_Plasmid_alt_A549)]
Siraj_clean_long <- melt(Siraj_clean, id.vars = "ID")
Siraj_clean_long[, celltype := gsub("(.*)_(A549|HEPG2|K562)", "\\2", variable)]
Siraj_clean_long[, stat := gsub("(.*)_(A549|HEPG2|K562)", "\\1", variable)]
Siraj_clean <- dcast(Siraj_clean_long, ID + celltype ~ stat, value.var = "value")
Siraj_clean[celltype == "HEPG2", celltype := "HepG2"]
Siraj_clean[, log2FC := as.numeric(log2FC)]
Siraj_clean[, log2Skew := as.numeric(log2Skew)]
# merge(Siraj_clean,Archaic_annot,by='ID')[,length(unique(ID)),by=.(POP)]


##### direct comparison : emVars log2FC, every CRE
figdata <- merge(emVARs_obs_ctc[condition == "NS"], Siraj_clean, by.x = c("ID", "celline"), by.y = c("ID", "celltype"))
p <- ggplot(figdata, aes(x = log2Skew, y = logFC_2vs1 / log(2), col = paste(emVar, is_emVar_CRITERIA))) +
  geom_point()
p <- p + facet_grid(rows = vars(celline)) + theme_plot(fontsize = 12) + xlab("log2 FC (Siraj et al.)") + ylab("log2 FC (this study)")

pdf(sprintf("%s/6c_log2FC_by_celltype__comparison_Siraj.pdf", FIGURE_DIR, CRITERIA_EMVARS), height = 4, width = 3)
print(p)
dev.off()



cat("\n===================================================================================================================\n")
cat("\n=== emVars only comparison : emVars log2FC, every CRE with an emvar in any condition in our data or Siraj et al ===\n")
cat("\n===================================================================================================================\n")

figdata <- merge(emVARs_obs_ctc, Siraj_clean, by.x = c("ID", "celline"), by.y = c("ID", "celltype"))
cre_select <- figdata[(emVar == "VRAI" | is_emVar_CRITERIA == TRUE), crsID]
cre_select_Siraj <- figdata[emVar == "VRAI", crsID]
cre_select_Li <- figdata[is_emVar_CRITERIA == TRUE, crsID]
cre_select_both <- figdata[(emVar == "VRAI" & is_emVar_CRITERIA == TRUE), crsID]


N_shared_tested <- figdata[, length(unique(crsID))]
N_emVar_any <- length(unique(cre_select))
N_emVar_both <- length(unique(cre_select_both))
N_emVar_Li <- length(unique(cre_select_Li))
N_emVar_Siraj <- length(unique(cre_select_Siraj))

N_emVar_both <- length(unique(cre_select_both))


cat(sprintf("\n\n Among %s shared tested variants between our setting and Siraj et al.,
			%s were detected as emVars in at least one study,
			with %s emVars being detected in both studies\n\n ", N_shared_tested, N_emVar_any, N_emVar_both))

Tab <- figdata[condition == "NS", .(Siraj = any(emVar == "VRAI", na.rm = T), Li = any(is_emVar_CRITERIA == TRUE)), by = crsID][, table(Siraj, Li)]
cat("\nNumber of variants called as emVars:\n")
Tab # Number of variants called as emVars
fisher.test(Tab)

cat("\n Number of emVars (variants x cellines, NS only) :\n")
Tab <- figdata[!is.na(emVar) & condition == "NS", table(emVar == "VRAI", is_emVar_CRITERIA)]
Tab
fisher.test(Tab)

cat("\n Number of emVars (variants x cellines, nay condition) :\n")
Tab <- figdata[!is.na(emVar), .(Siraj = any(emVar == "VRAI", na.rm = T), Li = any(is_emVar_CRITERIA == TRUE)), by = .(crsID, celline)][, table(Siraj, Li)]
Tab
fisher.test(Tab)


figdata[condition == "NS", .(Siraj = any(emVar == "VRAI", na.rm = T), Li = any(is_emVar_CRITERIA == TRUE), both = any(is_emVar_CRITERIA == TRUE & emVar == "VRAI")), by = "crsID"][Siraj & Li]

figdata <- figdata[crsID %in% cre_select, ]
figdata[, study := case_when(
  crsID %in% setdiff(cre_select_Siraj, cre_select_Li) ~ "Siraj et al",
  crsID %in% setdiff(cre_select_Li, cre_select_Siraj) ~ "This study",
  crsID %in% intersect(cre_select_Siraj, cre_select_Li) ~ "Both",
  TRUE ~ "WTF"
)]

p <- ggplot(figdata[condition == "NS", ], aes(x = log2Skew, y = logFC_2vs1))
p <- p + geom_hline(col = "lightgrey", yintercept = 0) + geom_vline(col = "lightgrey", xintercept = 0) + geom_abline(col = "lightgrey", intercept = 0, slope = 1, linetype = 2)
p <- p + geom_smooth(method = "lm", data = figdata[condition == "NS", ], aes(x = log2Skew, y = logFC_2vs1 / log(2)), col = "black", alpha = 0.3, size = .4, shape = 16)
p <- p + geom_point(aes(, col = COND_ID, alpha = study, shape = study))
p <- p + scale_color_manual(values = color_setup_simplified_norep) + scale_alpha_manual(values = c("Siraj et al" = 0.4, "This study" = 0.4, "Both" = 1))
p <- p + theme_plot(fontsize = 12) + xlab("log2 FC (Siraj et al.)") + ylab("log2 FC (this study)")
p <- p + guides(color = guide_legend(ncol = 1, byrow = TRUE), alpha = guide_legend(ncol = 1, byrow = TRUE))
p <- p + scale_shape_manual(values = c("Siraj et al" = 18, "This study" = 17, "Both" = 16))
pdf(sprintf("%s/6d_log2FC_by_celltype__comparison_Siraj__emVars_only.pdf", FIGURE_DIR, CRITERIA_EMVARS), height = 4, width = 3)
print(p)
dev.off()

cat("\n\nSpearman correlation of estimated allelic Log2FC between Siraj and this study (all conditions) \n")
figdata[, cor.test(log2Skew, logFC_2vs1, method = "s")]

cat("\n\nSpearman correlation of estimated allelic Log2FC between Siraj and this study (NS only) \n")

figdata[condition == "NS", cor.test(log2Skew, logFC_2vs1, method = "s")]

cat("\n\nPercentage of emVars with concordant effect size between Siraj and this study (all conditions) \n")

figdata[, binom.test(sum(sign(logFC_2vs1) == sign(log2Skew), na.rm = T), sum(!is.na(log2Skew)))]

cat("\n\nPercentage of emVars with concordant effect size between Siraj and this study (NS only) \n")

figdata[condition == "NS", binom.test(sum(sign(logFC_2vs1) == sign(log2Skew), na.rm = T), sum(!is.na(log2Skew)))]


cat("\n===================================================================================================================\n")
cat("\n===   emVars only comparison : unstimulated only, in the cell line where it's significant (in either study)     ===\n")
cat("\n===================================================================================================================\n")

##### emVars only comparison : emVars log2FC,
figdata <- merge(emVARs_obs_ctc, Siraj_clean, by.x = c("ID", "celline"), by.y = c("ID", "celltype"))
figdata <- figdata[condition == "NS" & (emVar == "VRAI" | is_emVar_CRITERIA == TRUE), ]
# figdata[,cor.test(log2Skew,logFC_2vs1,method='s')]
figdata[, cor.test(log2Skew, logFC_2vs1, method = "s")]
figdata[, study := case_when(
  emVar == "VRAI" & !is_emVar_CRITERIA ~ "Siraj et al",
  emVar == "FAUX" & is_emVar_CRITERIA ~ "This study",
  emVar == "VRAI" & is_emVar_CRITERIA ~ "Both",
  TRUE ~ "WTF"
)]

p <- ggplot(figdata[condition == "NS", ], aes(x = log2Skew, y = logFC_2vs1 / log(2)))
p <- p + geom_hline(col = "lightgrey", yintercept = 0) + geom_vline(col = "lightgrey", xintercept = 0) + geom_abline(col = "lightgrey", intercept = 0, slope = 1, linetype = 2)
p <- p + geom_smooth(method = "lm", data = figdata[condition == "NS", ], aes(x = log2Skew, y = logFC_2vs1), col = "black", alpha = 0.3, size = .4, shape = 16)
p <- p + geom_point(aes(, col = COND_ID, alpha = study, shape = study))
p <- p + scale_color_manual(values = color_setup_simplified_norep) + scale_alpha_manual(values = c("Siraj et al" = 0.4, "This study" = 0.4, "Both" = 1))
p <- p + theme_plot(fontsize = 12) + xlab("log2 FC (Siraj et al.)") + ylab("log2 FC (this study)")
p <- p + guides(color = guide_legend(ncol = 1, byrow = TRUE), alpha = guide_legend(ncol = 1, byrow = TRUE))
p <- p + scale_shape_manual(values = c("Siraj et al" = 18, "This study" = 17, "Both" = 16)) + ylim(c(-1, 1)) + xlim(c(-1, 1))

pdf(sprintf("%s/6e_log2FC_by_celltype__comparison_Siraj__emVars_NS_only.pdf", FIGURE_DIR, CRITERIA_EMVARS), height = 4, width = 3)
print(p)
dev.off()

cat("\n\nPercentage of emVars with concordant effect size between Siraj and this study (all conditions) \n")

figdata[, binom.test(sum(sign(logFC_2vs1) == sign(log2Skew), na.rm = T), sum(!is.na(log2Skew)))]
# number of successes = 15, number of trials = 16, p-value = 0.0005188
# alternative hypothesis: true probability of success is not equal to 0.5
# 95 percent confidence interval:
#  0.6976793 0.9984189
# sample estimates:
# probability of success
#                 0.9375


##### activity comparison
# p <- ggplot(figdata, aes(x=cre_class_CRITERIA,y=log2FC, fill=celline)) + geom_violin() + scale_fill_manual(values=color_celline)
# p <- p + facet_grid(rows = vars(celline)) + theme_plot(fontsize=12) + ylab('log2 FC Siraj') + xlab('CRE class (this study)')
# p <- p + geom_boxplot(width=0.2,fill='white',alpha=.5)
# pdf(sprintf("%s/2f_activity_by_celltype__comparison_Siraj.pdf", FIGURE_DIR, CRITERIA_EMVARS), height = 4, width = 3)
# print(p)
# dev.off()

cat("\n===================================================================================================================\n")
cat("\n===   emVars only comparison : unstimulated only, in the cell line where it's significant (in either study)     ===\n")
cat("\n===================================================================================================================\n")

figdata <- merge(emVARs_obs_ctc, Siraj_clean[!is.na(emVar), ], by.x = c("ID", "celline"), by.y = c("ID", "celltype"))
figdata <- figdata[condition == "NS", ]
figdata[,length(unique(crsID))]
# 160
Comparison_counts=figdata[,length(unique(crsID)),keyby=.(emVar,is_emVar_CRITERIA)]
Odds_ratio=Comparison_counts[,fisher.test(matrix(V1,2,2))]

Comparison_counts
#     emVar is_emVar_CRITERIA    V1
# 1:   FAUX             FALSE   159
# 2:   FAUX              TRUE    11
# 3:   VRAI             FALSE    14
# 4:   VRAI              TRUE     3

Odds_ratio

# Fisher's Exact Test for Count Data
# p-value = 0.121
# 95 percent confidence interval:
#   0.4933926 13.6606327
# odds ratio 
#     3.0705 

figdata <- figdata[(emVar == "VRAI" | is_emVar_CRITERIA == TRUE), ]
figdata[,length(unique(crsID))]
# 28

figdata[, log2FC_2vs1 := logFC_2vs1 / log(2)]
figdata[, study := case_when(
  emVar == "VRAI" & !is_emVar_CRITERIA ~ "Siraj et al.",
  emVar == "FAUX" & is_emVar_CRITERIA ~ "This study",
  emVar == "VRAI" & is_emVar_CRITERIA ~ "Both",
  TRUE ~ "WTF"
)]
figdata[,length(unique(crsID)),keyby=.(study)]
#           study    V1
# 1:         Both     3
# 2: Siraj et al.    14
# 3:   This study    11

figdata <- melt(figdata[study != "Both", .(celline, condition, ID, posID, study, log2FC_2vs1, log2Skew)], id.vars = c("study", "celline", "condition", "ID", "posID"))
figdata[, status := ifelse((study == "Siraj et al." & variable == "log2Skew") | (study == "This study" & variable == "log2FC_2vs1"), "Discovery", "Replication")]
figdata[, wilcox.test(abs(value) ~ status)$p.value, by = study]
#           study           Wilcox P-value
# 1: Siraj et al. 9.025415e-08
# 2:   This study 3.358339e-01

p <- ggplot(figdata, aes(x = status, y = abs(value), fill = status))
p <- p + geom_violin(scale = "width") + geom_boxplot(fill = "black", width = 0.2)
p <- p + scale_fill_manual(values = c("Discovery" = grey(0.8), "Replication" = grey(0.4)))
p <- p + theme_plot(fontsize = 12, lpos = "none", rotate.x = 90) + xlab("") + ylab("estimated |log2FC|\n(emVars)")
p <- p + facet_grid(cols = vars(study))

pdf(sprintf("%s/6f_log2FC_comparison_Siraj__Winners_curse_NS_only.pdf", FIGURE_DIR, CRITERIA_EMVARS), height = 2.5, width = 2)
print(p)
dev.off()


#############################################################################
#############################################################################
#############################################################################

# SupTable3ab <- emVARs_annot_obs_ctc[,.(celline,condition,crsID, ID,variantId_hg38,rsID,CHROM,POS_b37,
# INTROGRESSED=INTROGRESSED.allele, `NON-INTROGRESSED`=ifelse(INTROGRESSED.allele==REF,ALT,REF),
# introgression_locus,
# activity_allele1=a1, activity_allele2=a2, 
# cre_class_CRITERIA, is_tested_cre, log2FC_archaic_vs_modern, log2FC.se, pval.LRT, pval_emp, lfdr, FDR, is_emVar_CRITERIA,
# nObs_allele1 = nBCs_g1a1, nObs_allele2 = nBCs_g1a2,
# eQTL_in=effect_in,
# gwas_trait=gwas,
# POP_introgressed, Source_introgression, allele_match, POP_adaptive_inital_def, POP_adaptive_top, Introgression_scenario_top, Introgression_source_top, max_intro_allele_freq, max_intro_haplo_freq, 
# target_any_proximity,target_any_contact,target_celltype_proximity,target_celltype_contact, Nearest, Lung, Mono, Tcell)]


SupTable3ab <- emVARs_annot_obs_ctc[,.(variantId_hg38, rsID, crsID, CHROM,POS_b37,
INTROGRESSED_AND_OTHER=ifelse(INTROGRESSED.allele==REF,paste0(REF,'/',ALT),paste0(ALT,'/',REF)),celline,condition, cre_class_CRITERIA,log2FC_archaic_vs_modern,log2FC.se,pval_emp,FDR,is_emVar=ifelse(is_emVar_CRITERIA,'x',''))]

emVars_any=SupTable3ab[is_emVar=='x',unique(crsID)]
SupTable1b_4161_final=fread(sprintf("%s/figures/%s/Z_figures/SupTable_1/SupTable1b_4161_testedSNPs.tsv", MPRA_DIR, ANALYSIS_DIR))
SupTable1b_4161_final[,prioritized_gene:=paste(prioritized_distance_10kb, prioritized_myeloid_enhancer, prioritized_lymphoid_enhancer, prioritized_lung_enhancer, prioritized_covid19,sep=':')]
SupTable1b_4161_final[,prioritized_gene:=paste(setdiff(unlist(str_split(prioritized_gene,':')),c('_','-')),collapse=':'),by=posID]
SupTable1b_4161_final[prioritized_gene=='',prioritized_gene:=nearestGene]
SupTable3ab <- merge(SupTable3ab,SupTable1b_4161_final[,.(crsID=posID, POP_adaptive,maxFreq, Source_introgression,nearestGene,prioritized_gene)],by='crsID')

cols=c('cre_class_CRITERIA','log2FC_archaic_vs_modern','log2FC.se','pval_emp','FDR')
SupTable3ab_wide=dcast(SupTable3ab,CHROM+POS_b37+rsID+crsID+INTROGRESSED_AND_OTHER+POP_adaptive+maxFreq+Source_introgression+nearestGene+prioritized_gene~celline+condition,value.var=cols)
SupTable3ab_wide= SupTable3ab_wide[order(CHROM,POS_b37),]

new_col_order <- CJ(condition=condition_summary$COND_ID, vars=cols)
new_col_order <- new_col_order[order(factor(condition, levels = condition_summary$COND_ID),factor(vars, levels = unique(cols))),paste(vars,condition, sep = "_")]
setcolorder(SupTable3ab_wide, c(setdiff(names(SupTable3ab_wide), new_col_order), new_col_order))
SupTable3ab_wide[,is_emVar_any:=ifelse(crsID%in%emVars_any,'x','')]

if(!FILTER_ACTIVE){
  fwrite(SupTable3ab_wide , file=sprintf("%s/SupTables/SupTable3a_emVars.tsv", FIGURE_DIR),sep='\t')
}else{
  fwrite(SupTable3ab_wide , file=sprintf("%s/SupTables/SupTable3b_emVars.tsv", FIGURE_DIR),sep='\t')
}

strong_emVars_any=emVARs_annot_obs_ctc[is_emVar_CRITERIA==TRUE & abs(log2FC_archaic_vs_modern)>0.5,.(strong_emVar_in=paste(gsub('_','-',COND_ID),collapse=' / ')),by=crsID]
SupTable3c <- merge(SupTable3ab_wide,strong_emVars_any,by='crsID',all.x=TRUE)
SupTable3c <- SupTable3c[order(CHROM,POS_b37)]

fwrite(SupTable3c , file=sprintf("%s/SupTables/SupTable3c_strong_emVars.tsv", FIGURE_DIR),sep='\t')


SupTable5e_emVarts_short <- SupTable3ab[is_emVar=='x',.(where_emVar=paste(celline,condition,sep='-',collapse='/')),by=.(crsID,rsID,variantId_hg38,CHROM,POS_b37,POP_adaptive,maxFreq,Source_introgression,nearestGene,prioritized_gene,INTROGRESSED_AND_OTHER)]
CS_results_finemap_gwas_tested <- fread(sprintf("%s/data/%s/CredibleSets/CredibleSets_gwas.tsv.gz", MPRA_DIR, ANALYSIS_DIR))

SupTable5e <- merge(SupTable5e_emVarts_short,CS_results_finemap_gwas_tested[,.(finemappingMethod,variantId,full_name)],by.x='variantId_hg38',by.y='variantId',allow.cartesian=TRUE)
fwrite(SupTable5e[order(CHROM,POS_b37),] , file=sprintf("%s/SupTables/SupTable5e_gwas_emVars.tsv", FIGURE_DIR),sep='\t')





# #nIntegration_allele1 = nUMI_DNA_g1a1, nIntegration_allele2 = nUMI_DNA_g1a2,
# # oligo_num,ANALYSIS_SUBTYPE,ANALYSIS_TYPE,

# # nUMI_DNA_g1a1,nUMI_DNA_g1a2,nUMI_DNA_g2a1,nUMI_DNA_g2a2,
# # nUMI_RNA_g1a1,nUMI_RNA_g1a2,nUMI_RNA_g2a1,nUMI_RNA_g2a2,
# # nUMI_DNA_perBC_g1a1,nUMI_DNA_perBC_g1a2,nUMI_DNA_perBC_g2a1,nUMI_DNA_perBC_g2a2,
# # nUMI_RNA_perBC_g1a1,nUMI_RNA_perBC_g1a2,nUMI_RNA_perBC_g2a1,nUMI_RNA_perBC_g2a2,
# selection_criteria, selection_criteria.simple, Introgressed_in, Introgressed_from, POP_adaptive, Adaptive_from, 
# MaxPosterior_D_VY,MaxPosterior_V_DY,hap_length_kb_aSNP_A,hap_aSNP_nb_aSNP_A, linked_aSNP, cor_linked_aSNP, 
# size_introgression_locus, LD_block_AFR, 
# Vindija.der,Chagyrskaya.der,Altai.der,Denisova.der,
# YRI.der,SGDP_African.der,GBR.der,IBS.der,SGDP_WestEurasian.der,CHB.der,JPT.der,SGDP_Eastasian.der,PJL.der,STU.der,SGDP_Agta.der,SGDP_Papuan.der,
# allele.1,allele.2,allele2_is_REF,
# fixed_DER_AFR,fixed_ANC_AFR, 
# has_eQTL,has_gwas,
# strongest_effect_in,effect_in,
# blood_eQTL,lung_eQTL,liver_eQTL,
# Introgression_scenario_v2, selection_criteria, selection_criteria.simple,Introgressed_in, Introgressed_from, POP_adaptive, Adaptive_from, 
# POP_introgressed, Source_introgression, allele_match, 
# Introgression_scenario,POP_adaptive_inital_def,
# POP_introgressed, Source_introgression, POP_adaptive_top, Introgression_scenario_top, Introgression_source_top , max_intro_allele_freq,max_intro_haplo_freq,target_any_proximity,target_any_contact,target_celltype_proximity,target_celltype_contact)]