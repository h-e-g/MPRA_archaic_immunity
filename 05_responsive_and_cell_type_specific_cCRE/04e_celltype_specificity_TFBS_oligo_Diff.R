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
toc()

# source(sprintf("%s/scripts/%s/03_00_parameter_activity.R", MPRA_DIR, ANALYSIS_DIR))
source(sprintf("%s/scripts/%s/04_00_parameter_diff_activity.R", MPRA_DIR, ANALYSIS_DIR))

IN_DIR <- sprintf("%s/data/%s/03a_aggMPRA_analyse/%s/", MPRA_DIR, ANALYSIS_DIR, RUN_ID)
ACTIVE_DIR <- sprintf("%s/data/%s/03b_computeFDR/%s/Activity/%s/", MPRA_DIR, ANALYSIS_DIR, RUN_ID, CRITERIA_ACTIVE)
DIFF_DIR <- sprintf("%s/data/%s/03b_computeFDR/%s/Diff/%s/%s", MPRA_DIR, ANALYSIS_DIR, RUN_ID, CRITERIA_DIFF, CRITERIA_ACTIVE_OUT)

dir.create(DIFF_DIR, recursive = TRUE)
FIGURE_DIR <- sprintf("%s/figures/%s/%s/03_oligo_diff/%s/%s", MPRA_DIR, ANALYSIS_DIR, RUN_ID, CRITERIA_DIFF, CRITERIA_ACTIVE_OUT)
dir.create(FIGURE_DIR, recursive = TRUE)

# load oligo annotations
source(sprintf("%s/scripts/%s/02z__define_includedCRS.R", MPRA_DIR, ANALYSIS_DIR))


########################################################################################################
########################## Definition of cell type specific sites, same as Sup tables 2a-C #############
########################################################################################################


#### load differential activity between celltypes
oligo_activity_Diff_obs <- fread(sprintf("%s/all_oligos_diff_annotated__%s.tsv.gz", DIFF_DIR, CRITERIA_DIFF))
oligo_celltype_comp_NS_annot_obs <- oligo_activity_Diff_obs[grepl("celltype_comp", ANALYSIS_SUBTYPE) & grepl("_NS$", ANALYSIS_NAME), ]
oligo_celltype_comp_NS_annot_obs[, group1_labels := gsub("-ACE2", "", gsub("_NS_all", "", group1_labels))]
oligo_celltype_comp_NS_annot_obs[, group2_labels := gsub("-ACE2", "", gsub("_NS_all", "", group2_labels))]
oligo_celltype_comp_NS_annot_obs[, ANALYSIS_NAME := gsub("-ACE2", "", ANALYSIS_NAME)]
oligo_celltype_comp_NS_annot_obs_archaic <- oligo_celltype_comp_NS_annot_obs[oligo %chin% tested_and_ctrl_oligos_final_annot[type == "tested", oligo], ]

#### load activity stats
oligo_activity_obs <- fread(sprintf("%s/all_oligos_annotated__%s.tsv.gz", ACTIVE_DIR, CRITERIA_ACTIVE))
oligo_activity_obs_ctc <- oligo_activity_obs[power == 0 & boot == 0 & ANALYSIS_SUBTYPE == "celltype_cond", ]
oligo_activity_small <- oligo_activity_obs_ctc[condition == "NS", .(oligo, celline, oligo_class_CRITERIA, alpha_CRITERIA)]
oligo_activity_small <- dcast(oligo_activity_small, oligo ~ celline, value.var = c("oligo_class_CRITERIA", "alpha_CRITERIA"))

#### identify oligos that are cell type specific (upregulated)
oligo_K562_up <- intersect(
  oligo_celltype_comp_NS_annot_obs_archaic[ANALYSIS_NAME == "Celltype_HepG2_K562_NS" & oligo_diff_CRITERIA & log2FC_2vs1 > 0, oligo],
  oligo_celltype_comp_NS_annot_obs_archaic[ANALYSIS_NAME == "Celltype_K562_A549_NS" & oligo_diff_CRITERIA & log2FC_2vs1 < 0, oligo]
)
oligo_HepG2_up <- intersect(
  oligo_celltype_comp_NS_annot_obs_archaic[ANALYSIS_NAME == "Celltype_HepG2_K562_NS" & oligo_diff_CRITERIA & log2FC_2vs1 < 0, oligo],
  oligo_celltype_comp_NS_annot_obs_archaic[ANALYSIS_NAME == "Celltype_HepG2_A549_NS" & oligo_diff_CRITERIA & log2FC_2vs1 < 0, oligo]
)
oligo_A549_up <- intersect(
  oligo_celltype_comp_NS_annot_obs_archaic[ANALYSIS_NAME == "Celltype_K562_A549_NS" & oligo_diff_CRITERIA & log2FC_2vs1 > 0, oligo],
  oligo_celltype_comp_NS_annot_obs_archaic[ANALYSIS_NAME == "Celltype_HepG2_A549_NS" & oligo_diff_CRITERIA & log2FC_2vs1 > 0, oligo]
)

#### identify oligos that are cell type specific (DOWNregulated)
oligo_K562_lo <- intersect(
  oligo_celltype_comp_NS_annot_obs_archaic[ANALYSIS_NAME == "Celltype_HepG2_K562_NS" & oligo_diff_CRITERIA & log2FC_2vs1 < 0, oligo],
  oligo_celltype_comp_NS_annot_obs_archaic[ANALYSIS_NAME == "Celltype_K562_A549_NS" & oligo_diff_CRITERIA & log2FC_2vs1 > 0, oligo]
)
oligo_HepG2_lo <- intersect(
  oligo_celltype_comp_NS_annot_obs_archaic[ANALYSIS_NAME == "Celltype_HepG2_K562_NS" & oligo_diff_CRITERIA & log2FC_2vs1 > 0, oligo],
  oligo_celltype_comp_NS_annot_obs_archaic[ANALYSIS_NAME == "Celltype_HepG2_A549_NS" & oligo_diff_CRITERIA & log2FC_2vs1 > 0, oligo]
)
oligo_A549_lo <- intersect(
  oligo_celltype_comp_NS_annot_obs_archaic[ANALYSIS_NAME == "Celltype_K562_A549_NS" & oligo_diff_CRITERIA & log2FC_2vs1 < 0, oligo],
  oligo_celltype_comp_NS_annot_obs_archaic[ANALYSIS_NAME == "Celltype_HepG2_A549_NS" & oligo_diff_CRITERIA & log2FC_2vs1 < 0, oligo]
)


##############################################################################################################
###################### associate TFs with oligos #############################################################
##############################################################################################################

TF_score_annot <- fread(sprintf("%s/data/%s/00_TFscore/0_TFscore_all.tsv.gz", MPRA_DIR, ANALYSIS_DIR))

TF_score_annot <- TF_score_annot[oligo %chin% tested_and_ctrl_oligos_final_annot$oligo, ]
TF_score_annot <- merge(TF_score_annot, oligo_source[, .(oligo, GC)], by = "oligo", all.x = TRUE)
TF_score_annot[, score_GCadj := lm(score ~ GC)$residuals, by = .(matrix_id, names)]

TF_score_annot[, top10_GC := score_GCadj > quantile(score_GCadj, .9), by = .(matrix_id, names)]
TF_score_annot[, top10 := score > quantile(score, .9), by = .(matrix_id, names)]
TF_score_annot <- merge(TF_score_annot, tested_and_ctrl_oligos_final_annot[, .(oligo, posID, type, allele.label)], by = "oligo")

TF_score_annot[, K562_specific := oligo %chin% oligo_K562_up]
TF_score_annot[, HepG2_specific := oligo %chin% oligo_HepG2_up]
TF_score_annot[, A549_specific := oligo %chin% oligo_A549_up]


cre_TF_score <- TF_score_annot[, .(
  top10_GC = any(top10_GC),
  top10 = any(top10),
  K562_specific = any(K562_specific),
  HepG2_specific = any(HepG2_specific),
  A549_specific = any(A549_specific)
), by = .(posID, matrix_id, names)]

TF_enrich_CRS_K562 <- cre_TF_score[!is.na(matrix_id), as.list(unlist(fisher.test(table(K562_specific, top10_GC))[c("p.value", "estimate", "conf.int")])), by = .(matrix_id, names)]
TF_enrich_CRS_K562[, TF_ID := paste0(names, " (", matrix_id, ")")]
setnames(TF_enrich_CRS_K562, c("estimate.odds ratio", "conf.int1", "conf.int2"), c("OR", "OR_low", "OR_high"), skip_absent = TRUE)

TF_enrich_CRS_HepG2 <- cre_TF_score[!is.na(matrix_id), as.list(unlist(fisher.test(table(HepG2_specific, top10_GC))[c("p.value", "estimate", "conf.int")])), by = .(matrix_id, names)]
TF_enrich_CRS_HepG2[, TF_ID := paste0(names, " (", matrix_id, ")")]
setnames(TF_enrich_CRS_HepG2, c("estimate.odds ratio", "conf.int1", "conf.int2"), c("OR", "OR_low", "OR_high"), skip_absent = TRUE)

TF_enrich_CRS_A549 <- cre_TF_score[!is.na(matrix_id), as.list(unlist(fisher.test(table(A549_specific, top10_GC))[c("p.value", "estimate", "conf.int")])), by = .(matrix_id, names)]
TF_enrich_CRS_A549[, TF_ID := paste0(names, " (", matrix_id, ")")]
setnames(TF_enrich_CRS_A549, c("estimate.odds ratio", "conf.int1", "conf.int2"), c("OR", "OR_low", "OR_high"), skip_absent = TRUE)

TF_enrich_celltype_up <- list(K562 = TF_enrich_CRS_K562, HepG2 = TF_enrich_CRS_HepG2, A549 = TF_enrich_CRS_A549)
TF_enrich_celltype_up <- rbindlist(TF_enrich_celltype_up, idcol = "celltype")
TF_enrich_celltype_up[, FDR := p.adjust(p.value, method = "fdr"), by = .(celltype)]
TF_enrich_celltype_up[, type := "up"]



TF_score_annot[, K562_down := oligo %chin% oligo_K562_lo]
TF_score_annot[, HepG2_down := oligo %chin% oligo_HepG2_lo]
TF_score_annot[, A549_down := oligo %chin% oligo_A549_lo]


cre_TF_score <- TF_score_annot[, .(
  top10_GC = any(top10_GC),
  top10 = any(top10),
  K562_down = any(K562_down),
  HepG2_down = any(HepG2_down),
  A549_down = any(A549_down)
), by = .(posID, matrix_id, names)]


TF_enrich_CRS_K562_down <- cre_TF_score[!is.na(matrix_id), as.list(unlist(fisher.test(table(K562_down, top10_GC))[c("p.value", "estimate", "conf.int")])), by = .(matrix_id, names)]
TF_enrich_CRS_K562_down[, TF_ID := paste0(names, " (", matrix_id, ")")]
setnames(TF_enrich_CRS_K562_down, c("estimate.odds ratio", "conf.int1", "conf.int2"), c("OR", "OR_low", "OR_high"), skip_absent = TRUE)

TF_enrich_CRS_HepG2_down <- cre_TF_score[!is.na(matrix_id), as.list(unlist(fisher.test(table(HepG2_down, top10_GC))[c("p.value", "estimate", "conf.int")])), by = .(matrix_id, names)]
TF_enrich_CRS_HepG2_down[, TF_ID := paste0(names, " (", matrix_id, ")")]
setnames(TF_enrich_CRS_HepG2_down, c("estimate.odds ratio", "conf.int1", "conf.int2"), c("OR", "OR_low", "OR_high"), skip_absent = TRUE)

TF_enrich_CRS_A549_down <- cre_TF_score[!is.na(matrix_id), as.list(unlist(fisher.test(table(A549_down, top10_GC))[c("p.value", "estimate", "conf.int")])), by = .(matrix_id, names)]
TF_enrich_CRS_A549_down[, TF_ID := paste0(names, " (", matrix_id, ")")]
setnames(TF_enrich_CRS_A549_down, c("estimate.odds ratio", "conf.int1", "conf.int2"), c("OR", "OR_low", "OR_high"), skip_absent = TRUE)

TF_enrich_celltype_down <- list(K562 = TF_enrich_CRS_K562_down, HepG2 = TF_enrich_CRS_HepG2_down, A549 = TF_enrich_CRS_A549_down)
TF_enrich_celltype_down <- rbindlist(TF_enrich_celltype_down, idcol = "celltype")
TF_enrich_celltype_down[, FDR := p.adjust(p.value, method = "fdr"), by = .(celltype)]
TF_enrich_celltype_down[, type := "down"]

TF_enrich_celltype <- rbind(TF_enrich_celltype_up, TF_enrich_celltype_down)
TF_enrich_celltype <- TF_enrich_celltype[order(celltype, -type, p.value), ]
fwrite(TF_enrich_celltype, file = sprintf("%s/SupTable2d_celltype_specificity_TFBS_oligo_Diff__v%s.txt.gz", FIGURE_DIR, CRITERIA_DIFF), sep = "\t")
TF_enrich_celltype <- fread(sprintf("%s/SupTable2d_celltype_specificity_TFBS_oligo_Diff__v%s.txt.gz", FIGURE_DIR, CRITERIA_DIFF))
fwrite(TF_enrich_celltype[p.value<.05,],file=sprintf("%s/SupTable2d_celltype_specificity_TFBS_oligo_Diff__v%s_p05.txt", FIGURE_DIR, CRITERIA_DIFF),sep='\t')

plot_TFenrich_cellline <- function(TF_enrich_table, logged = TRUE, TF__order = NULL, CRS_group = "celltype-specific CRS", signif_Th = 0.05) {
  TF_enrich_table[, signif := FDR < signif_Th]
  if (!is.null(TF__order)) {
    TF_enrich_table[, TF_ID := factor(TF_ID, levels = TF__order, ordered = TRUE)]
  }

  if (logged == TRUE) {
    p <- ggplot(TF_enrich_table, aes(x = TF_ID, y = log2(OR), ymin = log2(OR_low), ymax = log2(OR_high), col = celltype, alpha = signif))
    OR_stat <- "log2(OR)"
  } else {
    p <- ggplot(TF_enrich_table, aes(x = TF_ID, y = OR, ymin = OR_low, ymax = OR_high, col = celltype, alpha = signif))
    OR_stat <- "OR"
  }
  p <- p + ylab(sprintf("Enrichment of TFBS in %s \n %s", CRS_group, OR_stat)) + xlab("TF")
  p <- p + geom_pointrange(size = .1) + geom_hline(yintercept = ifelse(logged, 0, 1), col = "grey")

  p <- p + theme_plot(lpos = "none", fontsize = 10) + coord_flip()
  p <- p + scale_color_manual(values = color_celline[c("A549", "HepG2", "K562")])
  p <- p + scale_alpha_manual(values = c("FALSE" = 0.5, "TRUE" = 1))
  p <- p + facet_grid(cols = vars(celltype), rows = vars(factor(type, c("up", "down"))), scales = "free_y", drop = TRUE)
  return(p)
}

pdf(sprintf("%s/Enrichment_TFBS_among_celltype_CRE_pointrange_%s.pdf", FIGURE_DIR, CRITERIA_DIFF), height = 3, width = 5)
TF__shown_up <- TF_enrich_celltype_up[order(-OR_low)][, .(TF_ID, names, specific_from = celltype)]
TF__shown_up <- TF__shown_up[!duplicated(names), head(.SD, 5), by = specific_from]
TF__shown_up <- TF__shown_up[order(factor(specific_from, levels = c("A549", "HepG2", "K562"))), ]
# TF__shown_lo= TF_enrich_celltype_down[order(OR_high)][FDR<.01,head(.SD,5),by=celltype][,unique(TF_ID)]
# TF__shown=c(TF__shown_up,TF__shown_lo)
TF__shown <- TF__shown_up$TF_ID
# p <- plot_TFenrich_cellline(TF_enrich_celltype[(TF_ID%in%TF__shown_up & type=="up") | (TF_ID%in%TF__shown_lo & type=="down"), ], logged=FALSE, TF__order=TF__shown,CRS_group="celltype specific CRS (FDR<1%)")
p <- plot_TFenrich_cellline(TF_enrich_celltype[(TF_ID %in% TF__shown_up$TF_ID & type == "up"), ], logged = TRUE, TF__order = rev(TF__shown), CRS_group = "celltype-specific CRE (up-regulated)", signif_Th = 0.05)
p <- p + theme(axis.text.y = element_text(colour = rep(color_celline[c("K562", "HepG2", "A549")], e = 5)))
print(p)
dev.off()

pdf(sprintf("%s/Enrichment_TFBS_among_celltype_CRE_downregulated_pointrange_%s.pdf", FIGURE_DIR, CRITERIA_DIFF), height = 3, width = 5)
# TF__shown_up <- TF_enrich_celltype_up[order(-OR_low)][, .(TF_ID, names, specific_from = celltype)]
# TF__shown_up <- TF__shown_up[!duplicated(name), head(.SD, 5), by = specific_from]
# TF__shown_up <- TF__shown_up[order(factor(specific_from, levels = c("A549", "HepG2", "K562"))), ]
TF__shown_lo <- TF_enrich_celltype_down[order(-OR_low)][, .(TF_ID, names, specific_from = celltype)]
TF__shown_lo <- TF__shown_lo[!duplicated(names), head(.SD, 5), by = specific_from]
TF__shown_lo <- TF__shown_lo[order(factor(specific_from, levels = c("A549", "HepG2", "K562"))), ]
# TF__shown=c(TF__shown_up,TF__shown_lo)
TF__shown <- TF__shown_lo$TF_ID
# p <- plot_TFenrich_cellline(TF_enrich_celltype[(TF_ID%in%TF__shown_lo & type=="up") | (TF_ID%in%TF__shown_lo & type=="down"), ], logged=FALSE, TF__order=TF__shown,CRS_group="celltype specific CRS (FDR<1%)")
p <- plot_TFenrich_cellline(TF_enrich_celltype[(TF_ID %in% TF__shown_lo$TF_ID & type == "up"), ], logged = TRUE, TF__order = rev(TF__shown), CRS_group = "celltype-specific CRE (down-regulated)", signif_Th = 0.05)
p <- p + theme(axis.text.y = element_text(colour = rep(color_celline[c("K562", "HepG2", "A549")], e = 5)))
print(p)
dev.off()

cat("All done!\n")
q("no")
