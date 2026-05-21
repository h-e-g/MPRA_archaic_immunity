MPRA_DIR <- "/pasteur/helix/projects/evo_immuno_pop/MPRA"
ANALYSIS_DIR <- "MPRA_count_exp6_analysis3"

source(sprintf("%s/scripts/%s/00_starter.R", MPRA_DIR, ANALYSIS_DIR))
source(sprintf("%s/scripts/%s/00_load_associations.R", MPRA_DIR, ANALYSIS_DIR))

RUN_ID <- "RUN3_Z2_nBC10"
# CRITERIA_ACTIVE <- "lfdr20_scrambled1pct_GCnorm"
#CRITERIA_ACTIVE <- "FDR5_scrambled5pct_FC0.2_GCnorm"
CRITERIA_ACTIVE <- "FDR5_FC1_FC0.2_GCnorm"

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
ACTIVITY_DIR <- sprintf("%s/02_oligo_activity/%s", FIGURE_DIR, CRITERIA_ACTIVE)
dir.create(ACTIVITY_DIR, showWarnings = FALSE, recursive = TRUE)

source(sprintf("%s/scripts/%s/02z__define_includedCRS.R", MPRA_DIR, ANALYSIS_DIR))


TFBM_annot <- fread(sprintf("%s/data/%s/00_TFBM_annot.tsv.gz", MPRA_DIR, ANALYSIS_DIR))
TF_score_annot <- fread(sprintf("%s/data/%s/00_TFscore/0_TFscore_all.tsv.gz", MPRA_DIR, ANALYSIS_DIR))
# TFBS_annot <- fread(sprintf("%s/data/%s/00_TFBM/0_TFBM_all.tsv.gz", MPRA_DIR, ANALYSIS_DIR))

oligo_activity_obs <- fread(sprintf("%s/all_oligos_annotated__%s.tsv.gz", ACTIVE_DIR, CRITERIA_ACTIVE))
oligo_activity_obs_ctc <- oligo_activity_obs[power == 0 & boot == 0 & ANALYSIS_SUBTYPE == "celltype_cond", ]
oligo_activity_obs_ctc_archaic <- oligo_activity_obs_ctc[oligo %in% tested_and_ctrl_oligos_final[type == "tested", oligo], ]

cat("ChromHMM, TFBM and distance , oligo activity ")

dir.create(sprintf("%s/ChromHMM_enrich", ACTIVITY_DIR), recursive = TRUE)
dir.create(sprintf("%s/TFBM_enrich", ACTIVITY_DIR), recursive = TRUE)
dir.create(sprintf("%s/Dist_enrich", ACTIVITY_DIR), recursive = TRUE)
############################################################################################
######## Enrichment of chromHMM states in enhancers and silencers (oligo level) : ##########
############################################################################################

# ################## all conditions ##################
# fig_data <- oligo_activity_obs_ctc_archaic[1:.N, ]
# fig_data[, chromHMM_full := ChromHMM_colors[match(ROADMAP_current_celline, paste(V1, V2, sep = "_")), V3]]
# fig_data[, chromHMM_full := factor(chromHMM_full, ChromHMM_colors$V3)]
# #  plot_data=fig_data[!is.na(chromHMM_full) & !is.na(CRS_class),.N,by=.(chromHMM_full,celline,condition,COND_ID,CRS_class,up=CRS_class=='enhancer',down=grepl('silencer',CRS_class),up_strong=CRS_class=='strong enhancer')]
# plot_data <- fig_data[!is.na(chromHMM_full) & !is.na(oligo_class_CRITERIA), .N, by = .(chromHMM_full, oligo_class_CRITERIA, up = oligo_class_CRITERIA == "enhancer", down = grepl("silencer", oligo_class_CRITERIA), up_strong = oligo_class_CRITERIA == "strong enhancer")]
# plot_data[, N_quies := sum(N[chromHMM_full == "Quiescent/Low"]), by = .(oligo_class_CRITERIA)]
# plot_data[, N_inact := sum(N[oligo_class_CRITERIA == "inactive"]), by = .(chromHMM_full)]
# plot_data[, N_quies_inact := sum(N[oligo_class_CRITERIA == "inactive" & chromHMM_full == "Quiescent/Low"])]
# plot_data[, c("OR", "P", "OR_lo", "OR_hi") := try(as.list(unlist(fisher.test(matrix(c(sum(N), sum(N_inact), sum(N_quies), sum(N_quies_inact)), 2))[c("estimate", "p.value", "conf.int")]))), by = .(chromHMM_full, oligo_class_CRITERIA)]
# plot_data[, FDR := p.adjust(P, "fdr")]


# ylim=range(plot_data[oligo_class_CRITERIA != "inactive" & oligo_class_CRITERIA != "strong silencer", log2(OR)])
# ylim[1]=min(floor(ylim[1]),-2)
# ylim[2]=max(ceiling(ylim[2]),4)
# pdf(sprintf("%s/ChromHMM_enrich/6a_Enrichment_chromHMM_per_oligo_class_pointrange__%s.pdf", ACTIVITY_DIR, CRITERIA_ACTIVE), height = 4, width = 2.5)
# p <- ggplot(plot_data[oligo_class_CRITERIA != "inactive" & oligo_class_CRITERIA != "strong silencer", ], aes(x = chromHMM_full, y = log2(OR), ymin = log2(OR_lo), ymax = log2(OR_hi), col = chromHMM_full, alpha = ifelse(FDR < 0.05, "signif", "ns")))
# p <- p + geom_pointrange(size = .6, fatten = 1) + geom_hline(yintercept = 0, col = "grey")
# p <- p + ylab("Enrichment in active CRS \n log2(OR)") + xlab("chromHmm class")
# p <- p + theme_plot(rotate.x = 90, lpos = "none") + scale_color_manual(values = setNames(ChromHMM_colors[, HEX], ChromHMM_colors[, V3]))
# p <- p + scale_alpha_manual(values = c(signif = 0.8, "ns" = 0.2))
# p <- p + facet_grid(rows = vars(factor(oligo_class_CRITERIA, c("strong enhancer", "enhancer", "silencer")))) + coord_cartesian(ylim = ylim)
# print(p)
# dev.off()

# fwrite(plot_data, file = sprintf("%s/6a_Enrichment_chromHMM_per_oligo_class_pointrange__%s__sourcedata.tsv", ACTIVITY_DIR, CRITERIA_ACTIVE), sep = "\t")

# ################## NS only (to avoid overconfidence) ##################
# fig_data <- oligo_activity_obs_ctc_archaic[condition == "NS", ]
# fig_data[, chromHMM_full := ChromHMM_colors[match(ROADMAP_current_celline, paste(V1, V2, sep = "_")), V3]]
# fig_data[, chromHMM_full := factor(chromHMM_full, ChromHMM_colors$V3)]
# #  plot_data=fig_data[!is.na(chromHMM_full) & !is.na(CRS_class),.N,by=.(chromHMM_full,celline,condition,COND_ID,CRS_class,up=CRS_class=='enhancer',down=grepl('silencer',CRS_class),up_strong=CRS_class=='strong enhancer')]
# plot_data <- fig_data[!is.na(chromHMM_full) & !is.na(oligo_class_CRITERIA), .N, by = .(chromHMM_full, oligo_class_CRITERIA, up = oligo_class_CRITERIA == "enhancer", down = grepl("silencer", oligo_class_CRITERIA), up_strong = oligo_class_CRITERIA == "strong enhancer")]
# plot_data[, N_quies := sum(N[chromHMM_full == "Quiescent/Low"]), by = .(oligo_class_CRITERIA)]
# plot_data[, N_inact := sum(N[oligo_class_CRITERIA == "inactive"]), by = .(chromHMM_full)]
# plot_data[, N_quies_inact := sum(N[oligo_class_CRITERIA == "inactive" & chromHMM_full == "Quiescent/Low"])]
# plot_data[, c("OR", "P", "OR_lo", "OR_hi") := try(as.list(unlist(fisher.test(matrix(c(sum(N), sum(N_inact), sum(N_quies), sum(N_quies_inact)), 2))[c("estimate", "p.value", "conf.int")]))), by = .(chromHMM_full, oligo_class_CRITERIA)]
# plot_data[, FDR := p.adjust(P, "fdr")]

# pdf(sprintf("%s/ChromHMM_enrich/6b_Enrichment_chromHMM_per_oligo_class_pointrange_NSonly__%s.pdf", ACTIVITY_DIR, CRITERIA_ACTIVE), height = 4, width = 2.5)
# p <- ggplot(plot_data[oligo_class_CRITERIA != "inactive" & oligo_class_CRITERIA != "strong silencer", ], aes(x = chromHMM_full, y = log2(OR), ymin = log2(OR_lo), ymax = log2(OR_hi), col = chromHMM_full, alpha = ifelse(FDR < 0.05, "signif", "ns")))
# p <- p + geom_pointrange(size = .6, fatten = 1) + geom_hline(yintercept = 0, col = "grey")
# p <- p + ylab("Enrichment in active CRS\n log2(OR)") + xlab("chromHmm class")
# p <- p + theme_plot(rotate.x = 90, lpos = "none") + scale_color_manual(values = setNames(ChromHMM_colors[, HEX], ChromHMM_colors[, V3]))
# p <- p + scale_alpha_manual(values = c(signif = 0.8, "ns" = 0.2))
# p <- p + facet_grid(rows = vars(factor(oligo_class_CRITERIA, c("strong enhancer", "enhancer", "silencer")))) + coord_cartesian(ylim = c(-2, 4))
# print(p)
# dev.off()
# fwrite(plot_data, file = sprintf("%s/ChromHMM_enrich/6b_Enrichment_chromHMM_per_oligo_class_pointrange_NSonly__%s__sourcedata.tsv", ACTIVITY_DIR, CRITERIA), sep = "\t")

# ################## by_celltype_condition (detail) ##################
# fig_data <- oligo_activity_obs_ctc_archaic
# fig_data[, chromHMM_full := ChromHMM_colors[match(ROADMAP_current_celline, paste(V1, V2, sep = "_")), V3]]
# fig_data[, chromHMM_full := factor(chromHMM_full, ChromHMM_colors$V3)]
# plot_data <- fig_data[!is.na(chromHMM_full) & !is.na(oligo_class_CRITERIA), .N, by = .(celline, condition, chromHMM_full, oligo_class_CRITERIA, up = oligo_class_CRITERIA == "enhancer", down = grepl("silencer", oligo_class_CRITERIA), up_strong = oligo_class_CRITERIA == "strong enhancer")]
# plot_data[, N_quies := sum(N[chromHMM_full == "Quiescent/Low"]), by = .(oligo_class_CRITERIA, celline, condition)]
# plot_data[, N_inact := sum(N[oligo_class_CRITERIA == "inactive"]), by = .(chromHMM_full, celline, condition)]
# plot_data[, N_quies_inact := sum(N[oligo_class_CRITERIA == "inactive" & chromHMM_full == "Quiescent/Low"]), by = .(celline, condition)]
# plot_data[, c("OR", "P", "OR_lo", "OR_hi") := try(as.list(unlist(fisher.test(matrix(c(sum(N), sum(N_inact), sum(N_quies), sum(N_quies_inact)), 2))[c("estimate", "p.value", "conf.int")]))), by = .(chromHMM_full, oligo_class_CRITERIA, celline, condition)]
# plot_data[, FDR := p.adjust(P, "fdr")]
# fwrite(plot_data, file = sprintf("%s/ChromHMM_enrich/6c_Enrichment_chromHMM_per_oligo_class_pointrange_detail__%s__sourcedata.tsv", ACTIVITY_DIR, CRITERIA), sep = "\t")

############################################################################################
######## Enrichment of chromHMM states in enhancers and silencers (cre level) : ############
############################################################################################

################## NS only (to avoid overconfidence) ##################
# preparte data
fig_data <- oligo_activity_obs_ctc_archaic[condition == "NS", ]
fig_data[, chromHMM_full := ChromHMM_colors[match(ROADMAP_current_celline, paste(V1, V2, sep = "_")), V3]]
fig_data[, chromHMM_full := factor(chromHMM_full, ChromHMM_colors$V3)]
# count crs, while keeping 0
unique_chromHMM_full <- unique(na.omit(fig_data$chromHMM_full))
unique_cre_class_CRITERIA <- unique(na.omit(fig_data$cre_class_CRITERIA))
all_combinations <- CJ(chromHMM_full = unique_chromHMM_full, cre_class_CRITERIA = unique_cre_class_CRITERIA)
agg_data <- fig_data[!is.na(chromHMM_full) & !is.na(cre_class_CRITERIA), .(N = length(unique(crsID))), by = .(chromHMM_full, cre_class_CRITERIA)]
plot_data <- merge(all_combinations, agg_data, by = c("chromHMM_full", "cre_class_CRITERIA"), all.x = TRUE)
plot_data[is.na(N), N := 0]
# compute odds ratios
plot_data[, N_quies := sum(N[chromHMM_full == "Quiescent/Low"]), by = .(cre_class_CRITERIA)]
plot_data[, N_inact := sum(N[cre_class_CRITERIA == "inactive"]), by = .(chromHMM_full)]
plot_data[, N_quies_inact := sum(N[cre_class_CRITERIA == "inactive" & chromHMM_full == "Quiescent/Low"])]
plot_data[, c("OR", "P", "OR_lo", "OR_hi") := try(as.list(unlist(fisher.test(matrix(c(sum(N), sum(N_inact), sum(N_quies), sum(N_quies_inact)), 2))[c("estimate", "p.value", "conf.int")]))), by = .(chromHMM_full, cre_class_CRITERIA)]
# add FDR and clean table
plot_data[, FDR := p.adjust(P, "fdr")]
plot_data <- plot_data[order(cre_class_CRITERIA, chromHMM_full)]
plot_data[FDR < .05, ]
pdf(sprintf("%s/ChromHMM_enrich/6d_Enrichment_chromHMM_per_cre_4class_pointrange_NSonly___%s.pdf", ACTIVITY_DIR, CRITERIA_ACTIVE), height = 5, width = 2.8)
p <- ggplot(plot_data[cre_class_CRITERIA != "inactive", ], aes(x = chromHMM_full, y = log2(OR), ymin = log2(OR_lo), ymax = log2(OR_hi), col = chromHMM_full, alpha = ifelse(FDR < 0.05, "signif", "ns")))
p <- p + geom_pointrange(size = .6, fatten = 1) + geom_hline(yintercept = 0, col = "grey")
p <- p + ylab("Enrichment in active CRE\n log2(OR)") + xlab("chromHmm class")
p <- p + theme_plot(rotate.x = 90, lpos = "none", fontsize = 11) + scale_color_manual(values = setNames(ChromHMM_colors[, HEX], ChromHMM_colors[, V3]))
p <- p + scale_alpha_manual(values = c(signif = 0.8, "ns" = 0.2))
p <- p + facet_grid(rows = vars(factor(cre_class_CRITERIA, c("strong enhancer", "enhancer", "silencer", "strong silencer"), labels = c("strong\nenhancer", "weak\nenhancer", "weak\nsilencer", "strong\nsilencer")))) + coord_cartesian(ylim = c(-2, 4))
print(p)
dev.off()
fwrite(plot_data, file = sprintf("%s/ChromHMM_enrich/6d_Enrichment_chromHMM_per_cre_4class_pointrange_NSonly__%s__sourcedata.tsv", ACTIVITY_DIR, CRITERIA_ACTIVE), sep = "\t")


### grouping weak and strong silencers
# prepare data
fig_data <- oligo_activity_obs_ctc_archaic[condition == "NS", ]
fig_data[, chromHMM_full := ChromHMM_colors[match(ROADMAP_current_celline, paste(V1, V2, sep = "_")), V3]]
fig_data[, chromHMM_full := factor(chromHMM_full, ChromHMM_colors$V3)]
fig_data[, cre_class3_CRITERIA := ifelse(cre_class_CRITERIA == "strong silencer", "silencer", cre_class_CRITERIA)]
fig_data[cre_class3_CRITERIA == "enhancer", cre_class3_CRITERIA := "weak enhancer"]
# count crs, while keeping 0
unique_chromHMM_full <- unique(na.omit(fig_data$chromHMM_full))
unique_cre_class3_CRITERIA <- unique(na.omit(fig_data$cre_class3_CRITERIA))
all_combinations <- CJ(chromHMM_full = unique_chromHMM_full, cre_class3_CRITERIA = unique_cre_class3_CRITERIA)
agg_data <- fig_data[!is.na(chromHMM_full) & !is.na(cre_class3_CRITERIA), .(N = length(unique(crsID))), by = .(chromHMM_full, cre_class3_CRITERIA)]
plot_data <- merge(all_combinations, agg_data, by = c("chromHMM_full", "cre_class3_CRITERIA"), all.x = TRUE)
plot_data[is.na(N), N := 0]
# compute odds ratios
plot_data[, N_quies := sum(N[chromHMM_full == "Quiescent/Low"]), by = .(cre_class3_CRITERIA)]
plot_data[, N_inact := sum(N[cre_class3_CRITERIA == "inactive"]), by = .(chromHMM_full)]
plot_data[, N_quies_inact := sum(N[cre_class3_CRITERIA == "inactive" & chromHMM_full == "Quiescent/Low"])]
plot_data[, c("OR", "P", "OR_lo", "OR_hi") := try(as.list(unlist(fisher.test(matrix(c(sum(N), sum(N_inact), sum(N_quies), sum(N_quies_inact)), 2))[c("estimate", "p.value", "conf.int")]))), by = .(chromHMM_full, cre_class3_CRITERIA)]
# add FDR and clean table
plot_data[, FDR := p.adjust(P, "fdr")]
plot_data <- plot_data[order(cre_class3_CRITERIA, chromHMM_full)]
plot_data[FDR < .05, ]
pdf(sprintf("%s/ChromHMM_enrich/6d_Enrichment_chromHMM_per_cre_3class_class_pointrange_NSonly__%s.pdf", ACTIVITY_DIR, CRITERIA_ACTIVE), height = 4, width = 2.8)
p <- ggplot(plot_data[cre_class3_CRITERIA != "inactive", ], aes(x = chromHMM_full, y = log2(OR), ymin = log2(OR_lo), ymax = log2(OR_hi), col = chromHMM_full, alpha = ifelse(FDR < 0.05, "signif", "ns")))
p <- p + geom_pointrange(size = .6, fatten = 1) + geom_hline(yintercept = 0, col = "grey")
p <- p + ylab("Enrichment in active CRE\n log2(OR)") + xlab("chromHmm class")
p <- p + theme_plot(rotate.x = 90, lpos = "none", fontsize = 10) + scale_color_manual(values = setNames(ChromHMM_colors[, HEX], ChromHMM_colors[, V3]))
p <- p + scale_alpha_manual(values = c(signif = 0.8, "ns" = 0.2))
p <- p + facet_grid(rows = vars(factor(cre_class3_CRITERIA, c("strong enhancer", "weak enhancer", "silencer"), labels = c("strong\nenhancer", "weak\nenhancer", "silencer")))) + coord_cartesian(ylim = c(-2, 4))
print(p)
dev.off()

fwrite(plot_data, file = sprintf("%s/ChromHMM_enrich/6d_Enrichment_chromHMM_per_cre_3class_pointrange_NSonly__%s__sourcedata.tsv", ACTIVITY_DIR, CRITERIA_ACTIVE), sep = "\t")
fwrite(plot_data, file = sprintf("%s/ChromHMM_enrich/SupTable1f_Enrichment_chromHMM_per_cre_3class_pointrange_NSonly__%s__sourcedata.tsv", ACTIVITY_DIR, CRITERIA_ACTIVE), sep = "\t")


### grouping weak and strong silencers & enhancers (2 class)
# prepare data
fig_data <- oligo_activity_obs_ctc_archaic[condition == "NS", ]
fig_data[, chromHMM_full := ChromHMM_colors[match(ROADMAP_current_celline, paste(V1, V2, sep = "_")), V3]]
fig_data[, chromHMM_full := factor(chromHMM_full, ChromHMM_colors$V3)]
fig_data[, cre_class3_CRITERIA := gsub('strong ','',cre_class_CRITERIA)]
# count crs, while keeping 0
unique_chromHMM_full <- unique(na.omit(fig_data$chromHMM_full))
unique_cre_class3_CRITERIA <- unique(na.omit(fig_data$cre_class3_CRITERIA))
all_combinations <- CJ(chromHMM_full = unique_chromHMM_full, cre_class3_CRITERIA = unique_cre_class3_CRITERIA)
agg_data <- fig_data[!is.na(chromHMM_full) & !is.na(cre_class3_CRITERIA), .(N = length(unique(crsID))), by = .(chromHMM_full, cre_class3_CRITERIA)]
plot_data <- merge(all_combinations, agg_data, by = c("chromHMM_full", "cre_class3_CRITERIA"), all.x = TRUE)
plot_data[is.na(N), N := 0]
# compute odds ratios
plot_data[, N_quies := sum(N[chromHMM_full == "Quiescent/Low"]), by = .(cre_class3_CRITERIA)]
plot_data[, N_inact := sum(N[cre_class3_CRITERIA == "inactive"]), by = .(chromHMM_full)]
plot_data[, N_quies_inact := sum(N[cre_class3_CRITERIA == "inactive" & chromHMM_full == "Quiescent/Low"])]
plot_data[, c("OR", "P", "OR_lo", "OR_hi") := try(as.list(unlist(fisher.test(matrix(c(sum(N), sum(N_inact), sum(N_quies), sum(N_quies_inact)), 2))[c("estimate", "p.value", "conf.int")]))), by = .(chromHMM_full, cre_class3_CRITERIA)]
# add FDR and clean table
plot_data[, FDR := p.adjust(P, "fdr")]
plot_data <- plot_data[order(cre_class3_CRITERIA, chromHMM_full)]
plot_data[FDR < .05, ]

pdf(sprintf("%s/ChromHMM_enrich/6d_Enrichment_chromHMM_per_cre_2class_class_pointrange_NSonly__%s.pdf", ACTIVITY_DIR, CRITERIA_ACTIVE), height = 4, width = 2.8)
p <- ggplot(plot_data[cre_class3_CRITERIA != "inactive", ], aes(x = chromHMM_full, y = log2(OR), ymin = log2(OR_lo), ymax = log2(OR_hi), col = chromHMM_full, alpha = ifelse(FDR < 0.05, "signif", "ns")))
p <- p + geom_pointrange(size = .6, fatten = 1) + geom_hline(yintercept = 0, col = "grey")
p <- p + ylab("Enrichment in active CRE\n log2(OR)") + xlab("chromHmm class")
p <- p + theme_plot(rotate.x = 90, lpos = "none", fontsize = 10) + scale_color_manual(values = setNames(ChromHMM_colors[, HEX], ChromHMM_colors[, V3]))
p <- p + scale_alpha_manual(values = c(signif = 0.8, "ns" = 0.2))
p <- p + facet_grid(rows = vars(factor(cre_class3_CRITERIA, c("enhancer", "silencer"), labels = c("positive CREs", "negative CREs")))) + coord_cartesian(ylim = c(-2, 4))
print(p)
dev.off()

fwrite(plot_data, file = sprintf("%s/ChromHMM_enrich/6d_Enrichment_chromHMM_per_cre_2class_pointrange_NSonly__%s__sourcedata.tsv", ACTIVITY_DIR, CRITERIA_ACTIVE), sep = "\t")
fwrite(plot_data, file = sprintf("%s/ChromHMM_enrich/SupTable1f_Enrichment_chromHMM_per_cre_2class_pointrange_NSonly__%s__sourcedata.tsv", ACTIVITY_DIR, CRITERIA_ACTIVE), sep = "\t")


################## by_celltype_condition (detail) ##################
fig_data <- oligo_activity_obs_ctc_archaic
fig_data[, chromHMM_full := ChromHMM_colors[match(ROADMAP_current_celline, paste(V1, V2, sep = "_")), V3]]
fig_data[, chromHMM_full := factor(chromHMM_full, ChromHMM_colors$V3)]
fig_data[, cre_class3_CRITERIA := ifelse(cre_class_CRITERIA == "strong silencer", "silencer", cre_class_CRITERIA)]
fig_data[cre_class3_CRITERIA == "enhancer", cre_class3_CRITERIA := "weak enhancer"]
# count crs, while keeping 0
unique_chromHMM_full <- unique(na.omit(fig_data$chromHMM_full))
unique_cre_class_CRITERIA <- unique(na.omit(fig_data$cre_class_CRITERIA))
unique_celline <- unique(na.omit(fig_data$celline))
unique_condition <- unique(na.omit(fig_data$condition))
all_combinations <- CJ(chromHMM_full = unique_chromHMM_full, cre_class_CRITERIA = unique_cre_class_CRITERIA, celline = unique_celline, condition = unique_condition)
agg_data <- fig_data[!is.na(chromHMM_full) & !is.na(cre_class_CRITERIA), .(N = length(unique(crsID))), by = .(chromHMM_full, cre_class_CRITERIA, celline, condition)]
plot_data <- merge(all_combinations, agg_data, by = c("chromHMM_full", "cre_class_CRITERIA", "celline", "condition"), all.x = TRUE)
plot_data[is.na(N), N := 0]
# compute odds ratios
plot_data <- fig_data[!is.na(chromHMM_full) & !is.na(cre_class_CRITERIA), .(N = length(unique(crsID))), by = .(celline, condition, chromHMM_full, cre_class_CRITERIA)]
plot_data[, N_quies := sum(N[chromHMM_full == "Quiescent/Low"]), by = .(cre_class_CRITERIA, celline, condition)]
plot_data[, N_inact := sum(N[cre_class_CRITERIA == "inactive"]), by = .(chromHMM_full, celline, condition)]
plot_data[, N_quies_inact := sum(N[cre_class_CRITERIA == "inactive" & chromHMM_full == "Quiescent/Low"]), by = .(celline, condition)]
plot_data[, c("OR", "P", "OR_lo", "OR_hi") := try(as.list(unlist(fisher.test(matrix(c(sum(N), sum(N_inact), sum(N_quies), sum(N_quies_inact)), 2))[c("estimate", "p.value", "conf.int")]))), by = .(chromHMM_full, cre_class_CRITERIA, celline, condition)]
plot_data <- plot_data[order(celline, relevel(factor(condition), "NS"), cre_class_CRITERIA, chromHMM_full)]
plot_data[, FDR := p.adjust(P, "fdr")]
fwrite(plot_data, file = sprintf("%s/ChromHMM_enrich/6e_Enrichment_chromHMM_per_cre_class4_pointrange_detail__%s__sourcedata.tsv", ACTIVITY_DIR, CRITERIA_ACTIVE), sep = "\t")

print(plot_data[FDR < .05, ])

###################################################################################
######## Enrichment of enhancers and silencers states in TFBM : ###################
###################################################################################

# # TFBS_count (not used)
# oligo_TFBS_count <- TFBS_annot[, .(n_TFBS = .N), by = .(oligo, matrix_id, names)]
# oligo_TFBS_count <- dcast(oligo_TFBS_count, oligo ~ matrix_id + names, value.var = "n_TFBS", fill = 0)
# oligo_TFBS_count <- melt(oligo_TFBS_count, id.vars = c("oligo"), value.name = "n_TFBS", variable.name = c("matrix_id.names"))
# oligo_TFBS_count <- oligo_TFBS_count[oligo %chin% tested_and_ctrl_oligos_final_annot$oligo, ]
# oligo_TFBS_count <- merge(oligo_TFBS_count, tested_and_ctrl_oligos_final_annot[, .(oligo, posID, type, allele.label)], by = "oligo")
# cre_TFBS_count <- oligo_TFBS_count[, .(n_TFBS = max(n_TFBS)), by = .(posID, matrix_id.names)]

# cre_class <- oligo_activity_obs_ctc_archaic[, .(
#   is_enhancer = any(grepl("enhancer", oligo_class_CRITERIA)),
#   is_silencer = any(grepl("silencer", oligo_class_CRITERIA)),
#   is_strong_enhancer = any(oligo_class_CRITERIA == "strong enhancer"),
#   is_strong_silencer = any(oligo_class_CRITERIA == "strong silencer"),
#   is_inactive = all(oligo_class_CRITERIA == "inactive")
# ), by = crsID]

# cre_class_long <- melt(cre_class, id.vars = c("crsID"), value.name = "is_active", variable.name = "class")
# fig_data <- merge(cre_class_long, cre_TFBS_count, by.x = "crsID", by.y = "posID", allow.cartesian = TRUE)
# fig_data <- fig_data[, .N, by = .(class, is_active, matrix_id.names, has_TFBM = n_TFBS > 0)]
# fig_data[, N_inact_noTFBM := N[is_active == TRUE & class == "is_inactive" & has_TFBM == FALSE], by = .(matrix_id.names)]
# fig_data[, N_inact_hasTFBM := N[is_active == TRUE & class == "is_inactive" & has_TFBM == TRUE], by = .(matrix_id.names)]
# fig_data[, N_active_noTFBM := N[is_active == TRUE & has_TFBM == FALSE], by = .(class, matrix_id.names)]
# fig_data[, N_active_hasTFBM := N[is_active == TRUE & has_TFBM == TRUE], by = .(class, matrix_id.names)]
# plot_data <- fig_data[is_active == TRUE & has_TFBM == TRUE & !is.na(N_inact_noTFBM) & !is.na(N_inact_hasTFBM) & !is.na(N_active_noTFBM) & !is.na(N_active_hasTFBM), ]
# plot_data[, c("OR", "P", "OR_lo", "OR_hi") := try(as.list(unlist(fisher.test(matrix(c(N_inact_noTFBM, N_inact_hasTFBM, N_active_noTFBM, N_active_hasTFBM), 2))[c("estimate", "p.value", "conf.int")]))), by = .(class, matrix_id.names)]
# plot_data[, FDR := p.adjust(P, "fdr")]
# fwrite(plot_data, file = sprintf("%s/6b_Enrichment_TFBM_per_class_pointrange_detail__%s__sourcedata.tsv", ACTIVITY_DIR, CRITERIA), sep = "\t")
# ##### results seem driven by GC content

# TF_score_annot=merge(TF_score_annot,CRS_TFBS_count,by=c('matrix_id','names'),all.x=TRUE)

######################################################################################################
######## Enrichment of enhancers and silencers states in high-TF affinity sequences ##################
######################################################################################################

TF_score_annot <- TF_score_annot[oligo %chin% tested_and_ctrl_oligos_final_annot$oligo, ]
TF_score_annot <- merge(TF_score_annot, oligo_source[, .(oligo, GC)], by = "oligo", all.x = TRUE)
TF_score_annot[, score_GCadj := lm(score ~ GC)$residuals, by = .(matrix_id, names)]

TF_score_annot[, top10_GC := score_GCadj > quantile(score_GCadj, .9), by = .(matrix_id, names)]
TF_score_annot[, top10 := score > quantile(score, .9), by = .(matrix_id, names)]
TF_score_annot <- merge(TF_score_annot, tested_and_ctrl_oligos_final_annot[, .(oligo, posID, type, allele.label)], by = "oligo")

# # assess enrichment, counting each oligo once
oligo_class <- oligo_activity_obs_ctc_archaic[, .(
  is_enhancer = any(grepl("enhancer", oligo_class_CRITERIA)),
  is_silencer = any(grepl("silencer", oligo_class_CRITERIA)),
  is_strong_enhancer = any(oligo_class_CRITERIA == "strong enhancer"),
  is_strong_silencer = any(oligo_class_CRITERIA == "strong silencer"),
  is_inactive = all(oligo_class_CRITERIA == "inactive")
), by = oligo]

TF_score_annot[, strong_enhancer := oligo %chin% oligo_class[is_strong_enhancer == TRUE, oligo]]
TF_score_annot[, weak_enhancer := oligo %chin% oligo_class[is_enhancer == TRUE & is_strong_enhancer == FALSE, oligo]]
TF_score_annot[, enhancer := oligo %chin% oligo_class[is_enhancer == TRUE, oligo]]
TF_score_annot[, silencer := oligo %chin% oligo_class[is_silencer == TRUE & is_strong_enhancer == FALSE & is_enhancer == FALSE, oligo]]
TF_score_annot[, weak_silencer := oligo %chin% oligo_class[is_silencer == TRUE & is_strong_enhancer == FALSE & is_enhancer == FALSE, oligo]]
TF_score_annot[, strong_silencer := oligo %chin% oligo_class[is_strong_silencer == TRUE, oligo]]

# ####### oligo level enrichment
# TF_score_annot[, TF_ID := paste0(names, " (", matrix_id, ")")]
# TF_enrich_CRS_enhancer <- TF_score_annot[!is.na(matrix_id), as.list(unlist(fisher.test(table(enhancer, top10_GC))[c("p.value", "estimate", "conf.int")])), by = .(matrix_id, names, TF_ID)]
# setnames(TF_enrich_CRS_enhancer, c("estimate.odds ratio", "conf.int1", "conf.int2"), c("OR", "OR_low", "OR_high"), skip_absent = TRUE)
# TF_enrich_CRS_enhancer[, type := "weak enhancer"]

# TF_enrich_CRS_enhancer_strong <- TF_score_annot[!is.na(matrix_id), as.list(unlist(fisher.test(table(strong_enhancer, top10_GC))[c("p.value", "estimate", "conf.int")])), by = .(matrix_id, names, TF_ID)]
# setnames(TF_enrich_CRS_enhancer_strong, c("estimate.odds ratio", "conf.int1", "conf.int2"), c("OR", "OR_low", "OR_high"), skip_absent = TRUE)
# TF_enrich_CRS_enhancer_strong[, type := "strong enhancer"]

# TF_enrich_CRS_silencer <- TF_score_annot[!is.na(matrix_id), as.list(unlist(fisher.test(table(silencer, top10_GC))[c("p.value", "estimate", "conf.int")])), by = .(matrix_id, names, TF_ID)]
# setnames(TF_enrich_CRS_silencer, c("estimate.odds ratio", "conf.int1", "conf.int2"), c("OR", "OR_low", "OR_high"), skip_absent = TRUE)
# TF_enrich_CRS_silencer[, type := "silencer"]

# TF_enrich_CRS_strong_silencer <- TF_score_annot[!is.na(matrix_id), as.list(unlist(fisher.test(table(strong_silencer, top10_GC))[c("p.value", "estimate", "conf.int")])), by = .(matrix_id, names, TF_ID)]
# setnames(TF_enrich_CRS_strong_silencer, c("estimate.odds ratio", "conf.int1", "conf.int2"), c("OR", "OR_low", "OR_high"), skip_absent = TRUE)
# TF_enrich_CRS_strong_silencer[, type := "strong silencer"]

# TF_enrich_CRS_weak_silencer <- TF_score_annot[!is.na(matrix_id), as.list(unlist(fisher.test(table(weak_silencer, top10_GC))[c("p.value", "estimate", "conf.int")])), by = .(matrix_id, names, TF_ID)]
# setnames(TF_enrich_CRS_weak_silencer, c("estimate.odds ratio", "conf.int1", "conf.int2"), c("OR", "OR_low", "OR_high"), skip_absent = TRUE)
# TF_enrich_CRS_weak_silencer[, type := "weak silencer"]

# TF_enrich_CRS_class <- rbind(TF_enrich_CRS_enhancer, TF_enrich_CRS_enhancer_strong, TF_enrich_CRS_silencer,TF_enrich_CRS_strong_silencer,TF_enrich_CRS_weak_silencer)
# TF_enrich_CRS_class[, FDR := p.adjust(p.value, method = "fdr"), by = type]

plot_TFenrich_class <- function(TF_enrich_table, logged = TRUE, TF__order = NULL, CRS_group = "", split_by_celltype = FALSE, include_strong_silencer = FALSE, include_strong = FALSE) {
  if (!is.null(TF__order)) {
    TF_enrich_table[, TF_ID := factor(TF_ID, levels = TF__order, ordered = TRUE)]
  }
  if (logged == TRUE) {
    p <- ggplot(TF_enrich_table, aes(x = TF_ID, y = log2(OR), ymin = log2(OR_low), ymax = log2(OR_high), col = type, alpha = FDR < 0.05))
    OR_stat <- "log2(OR)"
  } else {
    p <- ggplot(TF_enrich_table, aes(x = TF_ID, y = OR, ymin = OR_low, ymax = OR_high, col = type))
    OR_stat <- "OR"
  }
  p <- p + ylab(sprintf("Enrichment of TFBS in %s \n %s", CRS_group, OR_stat)) + xlab("TF")
  p <- p + geom_pointrange(size = .1) + geom_hline(yintercept = ifelse(logged, 0, 1), col = "grey")
  p <- p + theme_plot(rotate.x = 90, lpos = "none", fontsize = 10)
  p <- p + scale_color_manual(values = CRS_colors)
  p <- p + scale_alpha_manual(values = c("FALSE" = 0.5, "TRUE" = 1))

  if (include_strong_silencer) {
    crs_type <- c("strong enhancer", "weak enhancer", "weak silencer", "strong silencer")
  } else {
  	if (include_strong == TRUE) {
			crs_type <- c("strong enhancer", "weak enhancer", "silencer")
		}else{
			crs_type <- c("enhancer","silencer")
			}
  	}
  if (split_by_celltype == TRUE) {
    p <- p + facet_grid(rows = vars(factor(type, crs_type)), cols = vars(factor(celline, names(color_celline))), scales = "free_x", space = "free_x", drop = TRUE)
  } else {
    p <- p + facet_grid(rows = vars(factor(type, crs_type)), scales = "free_y", drop = TRUE)
  }
  return(p)
}


# pdf(sprintf("%s/TFBM_enrich/07a_Enrichment_top10TFscore_GCadj_3class_pointrange__oligo_level__%s.pdf", ACTIVITY_DIR, CRITERIA_ACTIVE), height = 4, width = 3)
# CRS_colors <- c("strong enhancer" = "#0f50b8", "weak enhancer" = "#4287f5", "silencer" = "#0fb88e")
# TF__shown <- TF_enrich_CRS_class[order(-OR_low)]
# TF__shown <- TF__shown[type%in%names(CRS_colors),][!duplicated(names)][FDR < .05, head(.SD, 5), by = type]
# TF__shown <- TF__shown[order(factor(type, levels = names(CRS_colors)), -OR_low), ]
# p <- plot_TFenrich_class(TF_enrich_CRS_class[(TF_ID %in% TF__shown$TF_ID & type%chin%names(CRS_colors)), ], logged = TRUE, TF__order = TF__shown$TF_ID, CRS_group = "active CRS")
# p <- p + theme(axis.text.x = element_text(colour = rep(CRS_colors, e = 5)))
# print(p)
# dev.off()


# pdf(sprintf("%s/TFBM_enrich/07b_Enrichment_top10TFscore_GCadj_4class_pointrange__oligo_level__%s.pdf", ACTIVITY_DIR, CRITERIA_ACTIVE), height = 4, width = 3)
# CRS_colors <- c("strong enhancer" = "#0f50b8", "weak enhancer" = "#4287f5", "weak silencer" = "#0fb88e",  "strong silencer" = "#047055")
# TF__shown <- TF_enrich_CRS_class[order(-OR_low)]
# TF__shown <- TF__shown[type%in%names(CRS_colors),][!duplicated(names)][, head(.SD, 5), by = type]
# TF__shown <- TF__shown[order(factor(type, levels = names(CRS_colors)), -OR_low), ]
# p <- plot_TFenrich_class(TF_enrich_CRS_class[(TF_ID %in% TF__shown$TF_ID & type%chin%names(CRS_colors)), ], logged = TRUE, TF__order = TF__shown$TF_ID, CRS_group = "active CRS", include_strong_silencer = TRUE)
# p <- p + theme(axis.text.x = element_text(colour = rep(CRS_colors, e = 5)))
# print(p)
# dev.off()

# fwrite( TF_enrich_CRS_class, file=sprintf("%s/TFBM_enrich/07ab_Enrichment_top10TFscore_GCadj_4class_pointrange__oligo_level__%s_sourcedata.tsv", ACTIVITY_DIR, CRITERIA), sep='\t')

############  assess enrichment, counting each cre once

# cre_class <- oligo_activity_obs_ctc_archaic[, .(
#   is_enhancer = any(grepl("enhancer", oligo_class_CRITERIA)),
#   is_silencer = any(grepl("silencer", oligo_class_CRITERIA)),
#   is_strong_enhancer = any(oligo_class_CRITERIA == "strong enhancer"),
#   is_strong_silencer = any(oligo_class_CRITERIA == "strong silencer"),
#   is_inactive = all(oligo_class_CRITERIA == "inactive")
# ),by = crsID]
# cre_class_long <- melt(cre_class, id.vars = c("crsID"), value.name = "is_active", variable.name = "class")
# fig_data <- merge(cre_class_long, cre_TF_score, by.x = "crsID", by.y = "posID", allow.cartesian = TRUE)
# fig_data <- fig_data[, .N, by = .(class, is_active, matrix_id, names, has_TFBM = top10_GC)]
# fig_data[, N_inact_noTFBM := N[is_active == TRUE & class == "is_inactive" & has_TFBM == FALSE], by = .(matrix_id, names)]
# fig_data[, N_inact_hasTFBM := N[is_active == TRUE & class == "is_inactive" & has_TFBM == TRUE], by = .(matrix_id, names)]
# fig_data[, N_active_noTFBM := N[is_active == TRUE & has_TFBM == FALSE], by = .(class, matrix_id, names)]
# fig_data[, N_active_hasTFBM := N[is_active == TRUE & has_TFBM == TRUE], by = .(class, matrix_id, names)]
# plot_data <- fig_data[is_active == TRUE & has_TFBM == TRUE & !is.na(N_inact_noTFBM) & !is.na(N_inact_hasTFBM) & !is.na(N_active_noTFBM) & !is.na(N_active_hasTFBM), ]
# plot_data[, c("OR", "P", "OR_lo", "OR_hi") := try(as.list(unlist(fisher.test(matrix(c(N_inact_noTFBM, N_inact_hasTFBM, N_active_noTFBM, N_active_hasTFBM), 2))[c("estimate", "p.value", "conf.int")]))), by = .(class, matrix_id, names)]
# plot_data[, FDR := p.adjust(P, "fdr")]
# fwrite(plot_data, file = sprintf("%s/6b_Enrichment_top10TFscore_GCadj_per_class_pointrange_detail__%s__sourcedata.tsv", ACTIVITY_DIR, CRITERIA_ACTIVE), sep = "\t")

cre_TF_score <- TF_score_annot[, .(
  top10_GC = any(top10_GC),
  top10 = any(top10),
  strong_enhancer = any(strong_enhancer),
	weak_enhancer = any(enhancer & !strong_enhancer),
  enhancer = any(enhancer),
  silencer = any(silencer & !enhancer & !strong_enhancer),
  weak_silencer = any(silencer & !enhancer & !strong_enhancer & !strong_silencer),
  strong_silencer = any(strong_silencer)
), by = .(posID, matrix_id, names)]


####### cre level enrichment
cre_TF_score[, TF_ID := paste0(names, " (", matrix_id, ")")]

TF_enrich_CRS_enhancer <- cre_TF_score[!is.na(matrix_id), as.list(unlist(fisher.test(table(enhancer, top10_GC))[c("p.value", "estimate", "conf.int")])), by = .(matrix_id, names, TF_ID)]
setnames(TF_enrich_CRS_enhancer, c("estimate.odds ratio", "conf.int1", "conf.int2"), c("OR", "OR_low", "OR_high"), skip_absent = TRUE)
TF_enrich_CRS_enhancer[, type := "enhancer"]

TF_enrich_CRS_enhancer_strong <- cre_TF_score[!is.na(matrix_id), as.list(unlist(fisher.test(table(strong_enhancer, top10_GC))[c("p.value", "estimate", "conf.int")])), by = .(matrix_id, names, TF_ID)]
setnames(TF_enrich_CRS_enhancer_strong, c("estimate.odds ratio", "conf.int1", "conf.int2"), c("OR", "OR_low", "OR_high"), skip_absent = TRUE)
TF_enrich_CRS_enhancer_strong[, type := "strong enhancer"]

TF_enrich_CRS_weak_enhancer <- cre_TF_score[!is.na(matrix_id), as.list(unlist(fisher.test(table(weak_enhancer, top10_GC))[c("p.value", "estimate", "conf.int")])), by = .(matrix_id, names, TF_ID)]
setnames(TF_enrich_CRS_weak_enhancer, c("estimate.odds ratio", "conf.int1", "conf.int2"), c("OR", "OR_low", "OR_high"), skip_absent = TRUE)
TF_enrich_CRS_weak_enhancer[, type := "weak enhancer"]

TF_enrich_CRS_silencer <- cre_TF_score[!is.na(matrix_id), as.list(unlist(fisher.test(table(silencer, top10_GC))[c("p.value", "estimate", "conf.int")])), by = .(matrix_id, names, TF_ID)]
setnames(TF_enrich_CRS_silencer, c("estimate.odds ratio", "conf.int1", "conf.int2"), c("OR", "OR_low", "OR_high"), skip_absent = TRUE)
TF_enrich_CRS_silencer[, type := "silencer"]

TF_enrich_CRS_strong_silencer <- cre_TF_score[!is.na(matrix_id), as.list(unlist(fisher.test(table(strong_silencer, top10_GC))[c("p.value", "estimate", "conf.int")])), by = .(matrix_id, names, TF_ID)]
setnames(TF_enrich_CRS_strong_silencer, c("estimate.odds ratio", "conf.int1", "conf.int2"), c("OR", "OR_low", "OR_high"), skip_absent = TRUE)
TF_enrich_CRS_strong_silencer[, type := "strong silencer"]

TF_enrich_CRS_weak_silencer <- cre_TF_score[!is.na(matrix_id), as.list(unlist(fisher.test(table(weak_silencer, top10_GC))[c("p.value", "estimate", "conf.int")])), by = .(matrix_id, names, TF_ID)]
setnames(TF_enrich_CRS_weak_silencer, c("estimate.odds ratio", "conf.int1", "conf.int2"), c("OR", "OR_low", "OR_high"), skip_absent = TRUE)
TF_enrich_CRS_weak_silencer[, type := "weak silencer"]

TF_enrich_CRS_class <- rbind(TF_enrich_CRS_enhancer, TF_enrich_CRS_weak_enhancer, TF_enrich_CRS_enhancer_strong, TF_enrich_CRS_silencer, TF_enrich_CRS_weak_silencer, TF_enrich_CRS_strong_silencer)
TF_enrich_CRS_class[, FDR := p.adjust(p.value, method = "fdr"), by = type]
TF_enrich_CRS_class <- merge(TF_enrich_CRS_class, unique(TFBM_annot[family != "", .(matrix_id, family)]), by = "matrix_id")
TF_enrich_CRS_class <- TF_enrich_CRS_class[order(type, -OR_low)]
fwrite(TF_enrich_CRS_class, file = sprintf("%s/TFBM_enrich/07cd_Enrichment_top10TFscore_GCadj_per_class_pointrange__cre_level__%s__sourcedata.tsv", ACTIVITY_DIR, CRITERIA_ACTIVE), sep = "\t")
fwrite(TF_enrich_CRS_class, file = sprintf("%s/TFBM_enrich/SupTable1g_Enrichment_top10TFscore_GCadj_per_class_pointrange__cre_level__%s__sourcedata.tsv", ACTIVITY_DIR, CRITERIA_ACTIVE), sep = "\t")
TF_enrich_CRS_class[FDR < .05, ][, .N, by = type]
TF_enrich_CRS_class <-fread(sprintf("%s/TFBM_enrich/07cd_Enrichment_top10TFscore_GCadj_per_class_pointrange__cre_level__%s__sourcedata.tsv", ACTIVITY_DIR, CRITERIA_ACTIVE), sep = "\t")

pdf(sprintf("%s/TFBM_enrich/07c_Enrichment_top10TFscore_GCadj_4class_pointrange__cre_level__%s.pdf", ACTIVITY_DIR, CRITERIA_ACTIVE), height = 5, width = 3)
CRS_colors <- c("strong enhancer" = "#0f50b8", "weak enhancer" = "#4287f5", "weak silencer" = "#0fb88e", "strong silencer" = "#047055")
TF__shown <- TF_enrich_CRS_class[order(-OR_low)]
TF__shown <- TF__shown[type %in% names(CRS_colors), ][!duplicated(names)][, head(.SD, 5), by = type]
TF__shown <- TF__shown[order(factor(type, levels = names(CRS_colors)), -OR_low), ]
p <- plot_TFenrich_class(TF_enrich_CRS_class[(TF_ID %in% TF__shown$TF_ID) & type %in% names(CRS_colors), ], logged = TRUE, TF__order = TF__shown$TF_ID, CRS_group = "active CRS", include_strong_silencer = TRUE)
p <- p + theme(axis.text.x = element_text(colour = rep(CRS_colors, e = 5)))
print(p)
dev.off()


pdf(sprintf("%s/TFBM_enrich/07d_Enrichment_top10TFscore_GCadj_3class_pointrange__cre_level__%s.pdf", ACTIVITY_DIR, CRITERIA_ACTIVE), height = 4, width = 3)
CRS_colors <- c("strong enhancer" = "#0f50b8", "weak enhancer" = "#4287f5", "silencer" = "#0fb88e")
TF__shown <- TF_enrich_CRS_class[order(-OR_low)]
TF__shown <- TF__shown[type %in% names(CRS_colors), ][!duplicated(names)]
TF__shown[, count_family := 0]
# TF__shown[, count_family := seq_len(.N),by=.(type,family)]
TF__shown <- TF__shown[count_family <= 2, head(.SD, 5), by = type]
TF__shown <- TF__shown[order(factor(type, levels = names(CRS_colors)), -OR_low), ]
p <- plot_TFenrich_class(TF_enrich_CRS_class[(TF_ID %in% TF__shown$TF_ID) & type %in% names(CRS_colors), ], logged = TRUE, TF__order = TF__shown$TF_ID, CRS_group = "active CRS", include_strong = TRUE)
p <- p + theme(axis.text.x = element_text(colour = rep(CRS_colors, e = 5)))
print(p)
dev.off()

pdf(sprintf("%s/TFBM_enrich/07d_Enrichment_top10TFscore_GCadj_3class_pointrange__cre_level__%s__top2_per_family.pdf", ACTIVITY_DIR, CRITERIA_ACTIVE), height = 4, width = 3)
CRS_colors <- c("strong enhancer" = "#0f50b8", "weak enhancer" = "#4287f5", "silencer" = "#0fb88e")
TF__shown <- TF_enrich_CRS_class[order(-OR_low)]
TF__shown <- TF__shown[type %in% names(CRS_colors), ][!duplicated(names)]
TF__shown[, count_family := seq_len(.N), by = .(type, family)]
TF__shown <- TF__shown[count_family <= 2, head(.SD, 5), by = type]
TF__shown <- TF__shown[order(factor(type, levels = names(CRS_colors)), -OR_low), ]
p <- plot_TFenrich_class(TF_enrich_CRS_class[(TF_ID %in% TF__shown$TF_ID) & type %in% names(CRS_colors), ], logged = TRUE, TF__order = TF__shown$TF_ID, CRS_group = "active CRS", include_strong = TRUE)
p <- p + theme(axis.text.x = element_text(colour = rep(CRS_colors, e = 5)))
print(p)
dev.off()


pdf(sprintf("%s/TFBM_enrich/07d_Enrichment_top10TFscore_GCadj_2class_pointrange__cre_level__%s__top2_per_family.pdf", ACTIVITY_DIR, CRITERIA_ACTIVE), height = 4, width = 3)
CRS_colors <- c("enhancer" = "#4287f5", "silencer" = "#0fb88e")
TF__shown <- TF_enrich_CRS_class[order(-OR_low)]
TF__shown <- TF__shown[type %in% names(CRS_colors), ][!duplicated(names)]
TF__shown[, count_family := seq_len(.N), by = .(type, family)]
TF__shown <- TF__shown[count_family <= 2, head(.SD, 5), by = type]
TF__shown <- TF__shown[order(factor(type, levels = names(CRS_colors)), -OR_low), ]
p <- plot_TFenrich_class(TF_enrich_CRS_class[(TF_ID %in% TF__shown$TF_ID) & type %in% names(CRS_colors), ], logged = TRUE, TF__order = TF__shown$TF_ID, CRS_group = "active CRS")
p <- p + theme(axis.text.x = element_text(colour = rep(CRS_colors, e = 5)))
print(p)
dev.off()


selected_TFs <- NULL
# Initialize an empty list to store the results for each type
current_data <- TF_enrich_CRS_class[order(-OR_low)]
# Loop over each type
for (myType in names(CRS_colors)) {
  # Select the top 5 TFs by score
  top_TFs <- head(current_data[type == myType, .(TF_ID, names)], 5)
  # Update the list of selected TFs
  selected_TFs <- rbind(selected_TFs, top_TFs)
  # Subset data for the current type, excluding already selected TFs
  current_data <- TF_enrich_CRS_class[!names %in% selected_TFs$names, ][order(-OR_low)]
}

pdf(sprintf("%s/TFBM_enrich/07d_Enrichment_top10TFscore_GCadj_3class_pointrange__cre_level__%s__with_dups.pdf", ACTIVITY_DIR, CRITERIA_ACTIVE), height = 4, width = 3)
p <- plot_TFenrich_class(TF_enrich_CRS_class[(TF_ID %in% selected_TFs$TF_ID) & type %in% names(CRS_colors), ], logged = TRUE, TF__order = selected_TFs$TF_ID, CRS_group = "active CRS")
p <- p + theme(axis.text.x = element_text(colour = rep(CRS_colors, e = 5)))
print(p)
dev.off()

########################################################
########## enrichment by class and cell line ###########
########################################################

################ assess enrichment separately per cell type ####################

cre_class_ct <- oligo_activity_obs_ctc_archaic[condition == "NS", .(
  is_enhancer = any(grepl("enhancer", oligo_class_CRITERIA)),
  is_strong_enhancer = any(oligo_class_CRITERIA == "strong enhancer"),
  is_weak_enhancer = any(grepl("enhancer", oligo_class_CRITERIA)) & !any(oligo_class_CRITERIA == "strong enhancer"),
  is_silencer = any(grepl("silencer", oligo_class_CRITERIA)) & !any(grepl("enhancer", oligo_class_CRITERIA)),
  is_strong_silencer = any(oligo_class_CRITERIA == "strong silencer"),
  is_weak_silencer = any(grepl("silencer", oligo_class_CRITERIA)) & !any(grepl("enhancer", oligo_class_CRITERIA)) & !any(oligo_class_CRITERIA == "strong silencer"),
  is_inactive = all(oligo_class_CRITERIA == "inactive")
), by = .(posID, celline)]

cre_TF_score <- TF_score_annot[, .(
  top10_GC = any(top10_GC),
  top10 = any(top10)
), by = .(posID, matrix_id, names)]

cre_TF_score_ct <- merge(cre_TF_score, cre_class_ct, by = c("posID"), allow.cartesian = TRUE)

cre_TF_score_ct[, TF_ID := paste0(names, " (", matrix_id, ")")]
TF_enrich_CRS_enhancer <- cre_TF_score_ct[!is.na(matrix_id), as.list(unlist(fisher.test(table(is_enhancer, top10_GC))[c("p.value", "estimate", "conf.int")])), by = .(matrix_id, names, TF_ID, celline)]
setnames(TF_enrich_CRS_enhancer, c("estimate.odds ratio", "conf.int1", "conf.int2"), c("OR", "OR_low", "OR_high"), skip_absent = TRUE)
TF_enrich_CRS_enhancer[, type := "enhancer"]

TF_enrich_CRS_weak_enhancer <- cre_TF_score_ct[!is.na(matrix_id), as.list(unlist(fisher.test(table(is_weak_enhancer, top10_GC))[c("p.value", "estimate", "conf.int")])), by = .(matrix_id, names, TF_ID, celline)]
setnames(TF_enrich_CRS_weak_enhancer, c("estimate.odds ratio", "conf.int1", "conf.int2"), c("OR", "OR_low", "OR_high"), skip_absent = TRUE)
TF_enrich_CRS_weak_enhancer[, type := "weak enhancer"]

TF_enrich_CRS_enhancer_strong <- cre_TF_score_ct[!is.na(matrix_id), as.list(unlist(fisher.test(table(is_strong_enhancer, top10_GC))[c("p.value", "estimate", "conf.int")])), by = .(matrix_id, names, TF_ID, celline)]
setnames(TF_enrich_CRS_enhancer_strong, c("estimate.odds ratio", "conf.int1", "conf.int2"), c("OR", "OR_low", "OR_high"), skip_absent = TRUE)
TF_enrich_CRS_enhancer_strong[, type := "strong enhancer"]

TF_enrich_CRS_silencer <- cre_TF_score_ct[!is.na(matrix_id), as.list(unlist(fisher.test(table(is_silencer, top10_GC))[c("p.value", "estimate", "conf.int")])), by = .(matrix_id, names, TF_ID, celline)]
setnames(TF_enrich_CRS_silencer, c("estimate.odds ratio", "conf.int1", "conf.int2"), c("OR", "OR_low", "OR_high"), skip_absent = TRUE)
TF_enrich_CRS_silencer[, type := "silencer"]

TF_enrich_CRS_weak_silencer <- cre_TF_score_ct[!is.na(matrix_id), as.list(unlist(fisher.test(table(is_weak_silencer, top10_GC))[c("p.value", "estimate", "conf.int")])), by = .(matrix_id, names, TF_ID, celline)]
setnames(TF_enrich_CRS_weak_silencer, c("estimate.odds ratio", "conf.int1", "conf.int2"), c("OR", "OR_low", "OR_high"), skip_absent = TRUE)
TF_enrich_CRS_weak_silencer[, type := "weak silencer"]

TF_enrich_CRS_strong_silencer <- cre_TF_score_ct[!is.na(matrix_id), as.list(unlist(fisher.test(table(factor(is_strong_silencer, levels = c(FALSE, TRUE)), top10_GC))[c("p.value", "estimate", "conf.int")])), by = .(matrix_id, names, TF_ID, celline)]
setnames(TF_enrich_CRS_strong_silencer, c("estimate.odds ratio", "conf.int1", "conf.int2"), c("OR", "OR_low", "OR_high"), skip_absent = TRUE)
TF_enrich_CRS_strong_silencer[, type := "strong silencer"]

TF_enrich_CRS_class <- rbind(TF_enrich_CRS_enhancer, TF_enrich_CRS_weak_enhancer, TF_enrich_CRS_enhancer_strong, TF_enrich_CRS_silencer, TF_enrich_CRS_weak_silencer, TF_enrich_CRS_strong_silencer)
TF_enrich_CRS_class[, FDR := p.adjust(p.value, method = "fdr"), by = type]

TF_enrich_CRS_class <- merge(TF_enrich_CRS_class, unique(TFBM_annot[family != "", .(matrix_id, family)]), by = "matrix_id")
TF_enrich_CRS_class <- TF_enrich_CRS_class[order(type, celline, -OR_low)]
fwrite(TF_enrich_CRS_class, file = sprintf("%s/TFBM_enrich/07ef_Enrichment_top10TFscore_GCadj_per_class_pointrange_by_ct__cre_level__%s__sourcedata.tsv", ACTIVITY_DIR, CRITERIA_ACTIVE), sep = "\t")
TF_enrich_CRS_class <- fread(sprintf("%s/TFBM_enrich/07ef_Enrichment_top10TFscore_GCadj_per_class_pointrange_by_ct__cre_level__%s__sourcedata.tsv", ACTIVITY_DIR, CRITERIA_ACTIVE), sep = "\t")

pdf(sprintf("%s/TFBM_enrich/07e_Enrichment_top10TFscore_GCadj_4class_pointrange_by_ct__cre_level__%s.pdf", ACTIVITY_DIR, CRITERIA_ACTIVE), height = 4, width = 6)
CRS_colors <- c("strong enhancer" = "#0f50b8", "weak enhancer" = "#4287f5", "weak silencer" = "#0fb88e", "strong silencer" = "#047055")
TF__shown <- TF_enrich_CRS_class[order(-OR_low)]
TF__shown <- TF__shown[type %in% names(CRS_colors), ][!duplicated(names)][, head(.SD, 5), by = .(celline, type)]
TF__shown <- TF__shown[order(factor(type, levels = names(CRS_colors)), celline, -OR_low), ]
p <- plot_TFenrich_class(TF_enrich_CRS_class[(paste(TF_ID, celline) %in% TF__shown[, paste(TF_ID, celline)] & type %in% names(CRS_colors)), ], logged = TRUE, TF__order = TF__shown$TF_ID, CRS_group = "active CRS", include_strong_silencer = TRUE, split_by_celltype = TRUE)
p <- p + theme(axis.text.x = element_text(colour = rep(CRS_colors, e = 5)))
print(p)
dev.off()

pdf(sprintf("%s/TFBM_enrich/07f_Enrichment_top10TFscore_GCadj_3class_pointrange_by_ct__cre_level__%s.pdf", ACTIVITY_DIR, CRITERIA_ACTIVE), height = 4, width = 6)
CRS_colors <- c("strong enhancer" = "#0f50b8", "weak enhancer" = "#4287f5", "silencer" = "#0fb88e")
TF__shown <- TF_enrich_CRS_class[order(-OR_low)]
TF__shown <- TF__shown[type %in% names(CRS_colors), ][!duplicated(names)][, head(.SD, 5), by = .(celline, type)]
TF__shown <- TF__shown[order(factor(type, levels = names(CRS_colors)), celline, -OR_low), ]
p <- plot_TFenrich_class(TF_enrich_CRS_class[(paste(TF_ID, celline) %in% TF__shown[, paste(TF_ID, celline)] & type %in% names(CRS_colors)), ], logged = TRUE, TF__order = TF__shown$TF_ID, CRS_group = "active CRS", split_by_celltype = TRUE, include_strong = TRUE)
p <- p + theme(axis.text.x = element_text(colour = rep(CRS_colors, e = 5)))
print(p)
dev.off()

pdf(sprintf("%s/TFBM_enrich/07f_Enrichment_top10TFscore_GCadj_3class_pointrange_by_ct__cre_level__%s__top2_per_family.pdf", ACTIVITY_DIR, CRITERIA_ACTIVE), height = 4, width = 6)
CRS_colors <- c("strong enhancer" = "#0f50b8", "weak enhancer" = "#4287f5", "silencer" = "#0fb88e")
TF__shown <- TF_enrich_CRS_class[order(-OR_low)]
TF__shown <- TF__shown[type %in% names(CRS_colors), ][!duplicated(names)]
TF__shown[, count_family := seq_len(.N), by = .(celline, type, family)]
TF__shown <- TF__shown[count_family <= 2, ][, head(.SD, 5), by = .(celline, type)]
TF__shown <- TF__shown[order(factor(type, levels = names(CRS_colors)), celline, -OR_low), ]
p <- plot_TFenrich_class(TF_enrich_CRS_class[(paste(TF_ID, celline) %in% TF__shown[, paste(TF_ID, celline)] & type %in% names(CRS_colors)), ], logged = TRUE, TF__order = TF__shown$TF_ID, CRS_group = "active CRS", split_by_celltype = TRUE, include_strong = TRUE)
p <- p + theme(axis.text.x = element_text(colour = rep(CRS_colors, e = 5)))
print(p)
dev.off()

pdf(sprintf("%s/TFBM_enrich/07f_Enrichment_top10TFscore_GCadj_2class_pointrange_by_ct__cre_level__%s__top2_per_family.pdf", ACTIVITY_DIR, CRITERIA_ACTIVE), height = 4, width = 6)
CRS_colors <- c("enhancer" = "#4287f5", "silencer" = "#0fb88e")
TF__shown <- TF_enrich_CRS_class[order(-OR_low)]
TF__shown <- TF__shown[type %in% names(CRS_colors), ][!duplicated(names)]
TF__shown[, count_family := seq_len(.N), by = .(celline, type, family)]
TF__shown <- TF__shown[count_family <= 2, ][, head(.SD, 5), by = .(celline, type)]
TF__shown <- TF__shown[order(factor(type, levels = names(CRS_colors)), celline, -OR_low), ]
p <- plot_TFenrich_class(TF_enrich_CRS_class[(paste(TF_ID, celline) %in% TF__shown[, paste(TF_ID, celline)] & type %in% names(CRS_colors)), ], logged = TRUE, TF__order = TF__shown$TF_ID, CRS_group = "active CRS", split_by_celltype = TRUE)
p <- p + theme(axis.text.x = element_text(colour = rep(CRS_colors, e = 5)))
print(p)
dev.off()

###########################################################################################
######## Enrichment of enhancers and silencers based oin distance to TSS ##################
###########################################################################################

oligo_source_v3 <- fread(sprintf("%s/data/%s/00_oligo_annot_v3_withTSS.txt", MPRA_DIR, ANALYSIS_DIR))
TSS_dist <- melt(oligo_source_v3[, .(oligo, crsID, Dist_NearestTSS_K562, Dist_NearestTSS_A549, Dist_NearestTSS_HepG2)], id.vars = c("oligo", "crsID"), value.name = "Dist_NearestTSS")
TSS_dist[, celline := gsub("Dist_NearestTSS_", "", variable)]
TSS_dist[, variable := NULL]
dist_th_kb_10bin <- c(0, 1, 2, 5, 10, 20, 50, 100, 200, 500, Inf)
dist_th_kb_5bin <- c(0, 2, 10, 50, 100, Inf)
TSS_dist[oligo %chin% tested_and_ctrl_oligos_final_annot$oligo, .(nOligo = .N, A549 = sum(celline == "A549"), HepG2 = sum(celline == "HepG2"), K562 = sum(celline == "K562")), keyby = cut(Dist_NearestTSS / 1000, dist_th_kb_10bin, include.lowest = TRUE)]
#           cut nOligo  A549 HepG2  K562
#        <fctr>  <int> <int> <int> <int>
#  1:     [0,1]    667   198   261   208
#  2:     (1,2]    503   165   174   164
#  3:     (2,5]   1036   311   419   306
#  4:    (5,10]   1398   444   540   414
#  5:   (10,20]   2071   669   841   561
#  6:   (20,50]   4224  1434  1475  1315
#  7:  (50,100]   3735  1259  1360  1116
#  8: (100,200]   4073  1336  1060  1677
#  9: (200,500]   4647  1698  1439  1510
# 10: (500,Inf]   1277   363   308   606


# oligo_Dist_ct <- merge(oligo_activity_obs_ctc_archaic, TSS_dist, by = c("oligo", "celline", "crsID"), allow.cartesian = TRUE)
# oligo_Dist_ct_count <- oligo_Dist_ct[, .N, by = .(oligo_class_CRITERIA = gsub("strong silencer", "silencer", oligo_class_CRITERIA), dist_bin = cut(Dist_NearestTSS / 1000, dist_th_kb_10bin, include.lowest = TRUE))]
# oligo_Dist_ct_count[, Ntot_bin := sum(N), by = .(dist_bin)]
# oligo_Dist_ct_count[, as.list(c(Ntot_bin, unlist(binom.test(N, Ntot_bin)[c("estimate", "conf.int")]))), by = .(oligo_class_CRITERIA, dist_bin)]
# oligo_Dist_ct_Pct <- oligo_Dist_ct_count[, as.list(unlist(binom.test(N, Ntot_bin)[c("estimate", "conf.int")])), by = .(oligo_class_CRITERIA, dist_bin)]
# setnames(oligo_Dist_ct_Pct, c("estimate.probability of success", "conf.int1", "conf.int2"), c("Pct", "Pct_lo", "Pct_hi"), skip_absent = TRUE)


# logged <- TRUE
# CRS_colors <- c("strong enhancer" = "#0f50b8", "enhancer" = "#4287f5", "silencer" = "#0fb88e")
# pdf(sprintf("%s/Dist_enrich/08a_Enrichment_Dist10Bin_pointrange__all_conditions__oligo_level__%s.pdf", ACTIVITY_DIR, CRITERIA_ACTIVE), height = 4, width = 3)
# p <- ggplot(oligo_Dist_ct_Pct[oligo_class_CRITERIA != "inactive"], aes(x = dist_bin, y = Pct, ymin = Pct_lo, ymax = Pct_hi, col = oligo_class_CRITERIA)) +
#   geom_pointrange()
# p <- p + ylab("Pct of oligos") + xlab("distance to nearest TSS")
# p <- p + geom_pointrange(size = .1, position = position_dodge(0.2)) + geom_hline(yintercept = ifelse(logged, 0, 1), col = "grey")
# p <- p + theme_plot(rotate.x = 90, lpos = "none")
# p <- p + scale_color_manual(values = CRS_colors)
# # p <- p + facet_grid(rows = vars(factor(celline,names(color_celline))), scales = "free_y", drop = TRUE)
# print(p)
# dev.off()


oligo_Dist_ct <- merge(oligo_activity_obs_ctc_archaic, TSS_dist, by = c("oligo", "celline", "crsID"), allow.cartesian = TRUE)
cre_Dist_ct_count <- oligo_Dist_ct[, .(N = length(unique(paste(celline, crsID)))), by = .(oligo_class_CRITERIA = gsub("strong silencer", "silencer", oligo_class_CRITERIA), dist_bin = cut(Dist_NearestTSS / 1000, dist_th_kb_10bin, include.lowest = TRUE))]
cre_Dist_ct_count[, Ntot_bin := sum(N), by = .(dist_bin)]
cre_Dist_ct_count[, as.list(c(Ntot_bin, unlist(binom.test(N, Ntot_bin)[c("estimate", "conf.int")]))), by = .(oligo_class_CRITERIA, dist_bin)]
cre_Dist_ct_Pct <- cre_Dist_ct_count[, as.list(c(N, Ntot_bin, unlist(binom.test(N, Ntot_bin)[c("estimate", "conf.int")]))), by = .(oligo_class_CRITERIA, dist_bin)]
setnames(cre_Dist_ct_Pct, c("estimate.probability of success", "conf.int1", "conf.int2"), c("Pct", "Pct_lo", "Pct_hi"), skip_absent = TRUE)

logged <- TRUE
CRS_colors <- c("strong enhancer" = "#0f50b8", "enhancer" = "#4287f5", "silencer" = "#0fb88e")
pdf(sprintf("%s/Dist_enrich/08b_Enrichment_Dist10Bin_pointrange__all_conditions__cre_level__%s.pdf", ACTIVITY_DIR, CRITERIA_ACTIVE), height = 4, width = 4)
p <- ggplot(cre_Dist_ct_Pct[oligo_class_CRITERIA != "inactive"], aes(x = dist_bin, y = Pct, ymin = Pct_lo, ymax = Pct_hi, col = oligo_class_CRITERIA))
p <- p + ylab("Percentage of active CRE (across celltypes)") + xlab("distance to nearest TSS  (kb)")
p <- p + geom_pointrange(size = .5, position = position_dodge(0.6)) + geom_hline(yintercept = ifelse(logged, 0, 1), col = "grey")
p <- p + theme_plot(rotate.x = 90, lpos = "right", fontsize = 12)
p <- p + scale_color_manual(values = CRS_colors)
# p <- p + facet_grid(rows = vars(factor(celline,names(color_celline))), scales = "free_y", drop = TRUE)
print(p)
dev.off()

logged <- TRUE
CRS_colors <- c("strong enhancer" = "#0f50b8", "enhancer" = "#4287f5", "silencer" = "#0fb88e")
pdf(sprintf("%s/Dist_enrich/08b_Enrichment_Dist10Bin_pointrange__all_conditions__cre_level_split__%s.pdf", ACTIVITY_DIR, CRITERIA_ACTIVE), height = 4, width = 3)
p <- ggplot(cre_Dist_ct_Pct[oligo_class_CRITERIA != "inactive"], aes(x = dist_bin, y = Pct, ymin = Pct_lo, ymax = Pct_hi, col = oligo_class_CRITERIA))
p <- p + ylab("Percentage of active CRE (across celltypes)") + xlab("distance to nearest TSS (kb)")
p <- p + geom_pointrange(size = .3, position = position_dodge(0.2)) + geom_hline(yintercept = ifelse(logged, 0, 1), col = "grey")
p <- p + theme_plot(rotate.x = 90, lpos = "none", fontsize = 12)
p <- p + scale_color_manual(values = CRS_colors)
p <- p + facet_grid(rows = vars(oligo_class_CRITERIA), drop = TRUE) # , scales = "free_y"
print(p)
dev.off()


oligo_Dist_ct <- merge(oligo_activity_obs_ctc_archaic, TSS_dist, by = c("oligo", "celline", "crsID"), allow.cartesian = TRUE)
cre_Dist_ct_count <- oligo_Dist_ct[, .(N = length(unique(paste(celline, crsID)))), by = .(oligo_class_CRITERIA = gsub("strong silencer", "silencer", oligo_class_CRITERIA), dist_bin = cut(Dist_NearestTSS / 1000, dist_th_kb_5bin, include.lowest = TRUE))]
cre_Dist_ct_count[, Ntot_bin := sum(N), by = .(dist_bin)]
cre_Dist_ct_count[, as.list(c(Ntot_bin, unlist(binom.test(N, Ntot_bin)[c("estimate", "conf.int")]))), by = .(oligo_class_CRITERIA, dist_bin)]
cre_Dist_ct_Pct <- cre_Dist_ct_count[, as.list(unlist(binom.test(N, Ntot_bin)[c("estimate", "conf.int")])), by = .(oligo_class_CRITERIA, dist_bin)]
setnames(cre_Dist_ct_Pct, c("estimate.probability of success", "conf.int1", "conf.int2"), c("Pct", "Pct_lo", "Pct_hi"), skip_absent = TRUE)

logged <- TRUE
CRS_colors <- c("strong enhancer" = "#0f50b8", "enhancer" = "#4287f5", "silencer" = "#0fb88e")
pdf(sprintf("%s/Dist_enrich/08c_Enrichment_Dist5Bin_pointrange__all_conditions__cre_level__%s.pdf", ACTIVITY_DIR, CRITERIA_ACTIVE), height = 4, width = 2.5)
p <- ggplot(cre_Dist_ct_Pct[oligo_class_CRITERIA != "inactive"], aes(x = dist_bin, y = Pct, ymin = Pct_lo, ymax = Pct_hi, col = oligo_class_CRITERIA))
p <- p + ylab("Percentage of active CRE (across celltypes)") + xlab("distance to nearest TSS (kb)")
p <- p + geom_pointrange(size = .5, position = position_dodge(0.7)) + geom_hline(yintercept = ifelse(logged, 0, 1), col = "grey")
p <- p + theme_plot(rotate.x = 90, lpos = "none", fontsize = 12)
p <- p + scale_color_manual(values = CRS_colors)
# p <- p + facet_grid(rows = vars(factor(celline,names(color_celline))), scales = "free_y", drop = TRUE)
print(p)
dev.off()

logged <- TRUE
CRS_colors <- c("strong enhancer" = "#0f50b8", "enhancer" = "#4287f5", "silencer" = "#0fb88e")
pdf(sprintf("%s/Dist_enrich/08c_Enrichment_Dist5Bin_pointrange__all_conditions__cre_level__split__%s.pdf", ACTIVITY_DIR, CRITERIA_ACTIVE), height = 5, width = 2)
p <- ggplot(cre_Dist_ct_Pct[oligo_class_CRITERIA != "inactive"], aes(x = dist_bin, y = Pct, ymin = Pct_lo, ymax = Pct_hi, col = oligo_class_CRITERIA))
p <- p + ylab("Percentage of active CRE (across celltypes)") + xlab("distance to\nnearest TSS (kb)")
p <- p + geom_pointrange(size = .3, position = position_dodge(0.7)) + geom_hline(yintercept = ifelse(logged, 0, 1), col = "grey")
p <- p + theme_plot(rotate.x = 90, lpos = "none", fontsize = 12)
p <- p + scale_color_manual(values = CRS_colors)
p <- p + facet_grid(rows = vars(factor(oligo_class_CRITERIA, c("strong enhancer", "enhancer", "silencer"), labels = c("strong enhancer", "weak enhancer", "silencer"))), drop = TRUE) # scales = "free_y"
print(p)
dev.off()


logged <- TRUE
CRS_colors <- c("enhancer" = "#4287f5", "silencer" = "#0fb88e")
cre_Dist_ct_count <- oligo_Dist_ct[, .(N = length(unique(paste(celline, crsID)))), by = .(oligo_class_CRITERIA = gsub("strong ", "", oligo_class_CRITERIA), dist_bin = cut(Dist_NearestTSS / 1000, dist_th_kb_5bin, include.lowest = TRUE))]
cre_Dist_ct_count[, Ntot_bin := sum(N), by = .(dist_bin)]
cre_Dist_ct_count[, as.list(c(Ntot_bin, unlist(binom.test(N, Ntot_bin)[c("estimate", "conf.int")]))), by = .(oligo_class_CRITERIA, dist_bin)]
cre_Dist_ct_Pct <- cre_Dist_ct_count[, as.list(unlist(binom.test(N, Ntot_bin)[c("estimate", "conf.int")])), by = .(oligo_class_CRITERIA, dist_bin)]
setnames(cre_Dist_ct_Pct, c("estimate.probability of success", "conf.int1", "conf.int2"), c("Pct", "Pct_lo", "Pct_hi"), skip_absent = TRUE)
pdf(sprintf("%s/Dist_enrich/08d_Enrichment_Dist5Bin_pointrange__all_conditions__2cre_level__split__%s.pdf", ACTIVITY_DIR, CRITERIA_ACTIVE), height = 5, width = 2)
p <- ggplot(cre_Dist_ct_Pct[oligo_class_CRITERIA != "inactive"], aes(x = dist_bin, y = Pct, ymin = Pct_lo, ymax = Pct_hi, col = oligo_class_CRITERIA))
p <- p + ylab("Percentage of active CRE (across celltypes)") + xlab("distance to\nnearest TSS (kb)")
p <- p + geom_pointrange(size = .3, position = position_dodge(0.7)) + geom_hline(yintercept = ifelse(logged, 0, 1), col = "grey")
p <- p + theme_plot(rotate.x = 90, lpos = "none", fontsize = 12)
p <- p + scale_color_manual(values = CRS_colors)
p <- p + facet_grid(rows = vars(factor(oligo_class_CRITERIA, c("enhancer", "silencer"), labels = c("positive CREs", "negative CREs"))), drop = TRUE) # scales = "free_y"
print(p)
dev.off()
# cre_Dist_ct_count <- oligo_Dist_ct[, .(N = length(unique(paste(crsID)))), by = .(oligo_class_CRITERIA = gsub("strong silencer", "silencer", oligo_class_CRITERIA), dist_bin = cut(Dist_NearestTSS / 1000, dist_th_kb_5bin, include.lowest = TRUE), celline)]
# cre_Dist_ct_count[, Ntot_bin := sum(N), by = .(dist_bin, celline)]
# cre_Dist_ct_count[, as.list(c(Ntot_bin, unlist(binom.test(N, Ntot_bin)[c("estimate", "conf.int")]))), by = .(oligo_class_CRITERIA, dist_bin, celline)]
# cre_Dist_ct_Pct <- cre_Dist_ct_count[, as.list(unlist(binom.test(N, Ntot_bin)[c("estimate", "conf.int")])), by = .(oligo_class_CRITERIA, dist_bin, celline)]
# setnames(cre_Dist_ct_Pct, c("estimate.probability of success", "conf.int1", "conf.int2"), c("Pct", "Pct_lo", "Pct_hi"), skip_absent = TRUE)

# logged <- TRUE
# CRS_colors <- c("strong enhancer" = "#0f50b8", "enhancer" = "#4287f5", "silencer" = "#0fb88e")
# pdf(sprintf("%s/08d_Enrichment_Dist5Bin_pointrange_by_ct__all_conditions__cre_level__%s.pdf", ACTIVITY_DIR, CRITERIA_ACTIVE), height = 4, width = 3)
# p <- ggplot(cre_Dist_ct_Pct[oligo_class_CRITERIA != "inactive"], aes(x = dist_bin, y = Pct, ymin = Pct_lo, ymax = Pct_hi, col = oligo_class_CRITERIA)) +
#   geom_pointrange()
# p <- p + ylab("Pct of cre ") + xlab("distance to nearest TSS")
# p <- p + geom_pointrange(size = .1, position = position_dodge(0.2)) + geom_hline(yintercept = ifelse(logged, 0, 1), col = "grey")
# p <- p + theme_plot(rotate.x = 90, lpos = "none")
# p <- p + scale_color_manual(values = CRS_colors)
# p <- p + facet_grid(rows = vars(factor(celline, names(color_celline))), scales = "free_y", drop = TRUE)
# print(p)
# dev.off()

# oligo_Dist_any <- merge(oligo_activity_obs_ctc_archaic, oligo_source_v3[, .(oligo, Dist_NearestGeneStart)], by = c("oligo"), allow.cartesian = TRUE)
# cre_Dist_any_count <- oligo_Dist_any[, .(N = length(unique(paste(crsID)))), by = .(oligo_class_CRITERIA = gsub("strong silencer", "silencer", oligo_class_CRITERIA), dist_bin = cut(Dist_NearestGeneStart / 1000, dist_th_kb_10bin, include.lowest = TRUE))]
# cre_Dist_any_count[, Ntot_bin := sum(N), by = .(dist_bin)]
# cre_Dist_any_count[, as.list(c(Ntot_bin, unlist(binom.test(N, Ntot_bin)[c("estimate", "conf.int")]))), by = .(oligo_class_CRITERIA, dist_bin)]
# cre_Dist_any_Pct <- cre_Dist_any_count[, as.list(unlist(binom.test(N, Ntot_bin)[c("estimate", "conf.int")])), by = .(oligo_class_CRITERIA, dist_bin)]
# setnames(cre_Dist_any_Pct, c("estimate.probability of success", "conf.int1", "conf.int2"), c("Pct", "Pct_lo", "Pct_hi"), skip_absent = TRUE)

# logged <- TRUE
# CRS_colors <- c("strong enhancer" = "#0f50b8", "enhancer" = "#4287f5", "silencer" = "#0fb88e")
# pdf(sprintf("%s/08e_Enrichment_Dist10Bin_pointrange__cre_level__%s.pdf", ACTIVITY_DIR, CRITERIA_ACTIVE), height = 4, width = 3)
# p <- ggplot(cre_Dist_any_Pct[oligo_class_CRITERIA != "inactive"], aes(x = dist_bin, y = Pct, ymin = Pct_lo, ymax = Pct_hi, col = oligo_class_CRITERIA)) +
#   geom_pointrange()
# p <- p + ylab("Pct of cres") + xlab("distance to nearest TSS")
# p <- p + geom_pointrange(size = .1, position = position_dodge(0.2)) + geom_hline(yintercept = ifelse(logged, 0, 1), col = "grey")
# p <- p + theme_plot(rotate.x = 90, lpos = "none")
# p <- p + scale_color_manual(values = CRS_colors)
# # p <- p + facet_grid(rows = vars(factor(celline,names(color_celline))), scales = "free_y", drop = TRUE)
# print(p)
# dev.off()



# oligo_Dist_ct <- merge(oligo_activity_obs_ctc_archaic, TSS_dist, by = c("oligo", "celline", "crsID"), allow.cartesian = TRUE)
# oligo_Dist_ct_count <- oligo_Dist_ct[condition == "NS", .N, by = .(oligo_class_CRITERIA = gsub("strong silencer", "silencer", oligo_class_CRITERIA), dist_bin = cut(Dist_NearestTSS / 1000, dist_th_kb_10bin, include.lowest = TRUE))]
# oligo_Dist_ct_count[, Ntot_bin := sum(N), by = .(dist_bin)]
# oligo_Dist_ct_count[, as.list(c(Ntot_bin, unlist(binom.test(N, Ntot_bin)[c("estimate", "conf.int")]))), by = .(oligo_class_CRITERIA, dist_bin)]
# oligo_Dist_ct_Pct <- oligo_Dist_ct_count[, as.list(unlist(binom.test(N, Ntot_bin)[c("estimate", "conf.int")])), by = .(oligo_class_CRITERIA, dist_bin)]
# setnames(oligo_Dist_ct_Pct, c("estimate.probability of success", "conf.int1", "conf.int2"), c("Pct", "Pct_lo", "Pct_hi"), skip_absent = TRUE)


# logged <- TRUE
# CRS_colors <- c("strong enhancer" = "#0f50b8", "enhancer" = "#4287f5", "silencer" = "#0fb88e")
# pdf(sprintf("%s/08f_Enrichment_Dist10Bin_pointrange__NSonly__oligo_level__%s.pdf", ACTIVITY_DIR, CRITERIA_ACTIVE), height = 4, width = 3)
# p <- ggplot(oligo_Dist_ct_Pct[oligo_class_CRITERIA != "inactive"], aes(x = dist_bin, y = Pct, ymin = Pct_lo, ymax = Pct_hi, col = oligo_class_CRITERIA)) +
#   geom_pointrange()
# p <- p + ylab("Pct of oligos") + xlab("distance to nearest TSS")
# p <- p + geom_pointrange(size = .1, position = position_dodge(0.2)) + geom_hline(yintercept = ifelse(logged, 0, 1), col = "grey")
# p <- p + theme_plot(rotate.x = 90, lpos = "none")
# p <- p + scale_color_manual(values = CRS_colors)
# # p <- p + facet_grid(rows = vars(factor(celline,names(color_celline))), scales = "free_y", drop = TRUE)
# print(p)
# dev.off()


# oligo_Dist_ct <- merge(oligo_activity_obs_ctc_archaic, TSS_dist, by = c("oligo", "celline", "crsID"), allow.cartesian = TRUE)
# oligo_Dist_ct_count <- oligo_Dist_ct[, .N, by = .(oligo_class_CRITERIA = gsub("strong silencer", "silencer", oligo_class_CRITERIA), celline, condition, dist_bin = cut(Dist_NearestTSS / 1000, dist_th_kb_10bin, include.lowest = TRUE))]
# oligo_Dist_ct_count[, Ntot_bin := sum(N), by = .(celline, condition, dist_bin)]
# oligo_Dist_ct_count[, as.list(c(Ntot_bin, unlist(binom.test(N, Ntot_bin)[c("estimate", "conf.int")]))), by = .(oligo_class_CRITERIA, celline, condition, dist_bin)]
# oligo_Dist_ct_Pct <- oligo_Dist_ct_count[, as.list(unlist(binom.test(N, Ntot_bin)[c("estimate", "conf.int")])), by = .(oligo_class_CRITERIA, celline, condition, dist_bin)]
# setnames(oligo_Dist_ct_Pct, c("estimate.probability of success", "conf.int1", "conf.int2"), c("Pct", "Pct_lo", "Pct_hi"), skip_absent = TRUE)


# logged <- TRUE
# CRS_colors <- c("strong enhancer" = "#0f50b8", "enhancer" = "#4287f5", "silencer" = "#0fb88e")
# pdf(sprintf("%s/08g_Enrichment_Dist10Bin_pointrange_by_ct__oligo_level__%s.pdf", ACTIVITY_DIR, CRITERIA_ACTIVE), height = 4, width = 3)
# p <- ggplot(oligo_Dist_ct_Pct[oligo_class_CRITERIA != "inactive" & condition == "NS"], aes(x = dist_bin, y = Pct, ymin = Pct_lo, ymax = Pct_hi, col = oligo_class_CRITERIA)) +
#   geom_pointrange()
# p <- p + ylab("Pct of oligos") + xlab("distance to nearest TSS")
# p <- p + geom_pointrange(size = .1, position = position_dodge(0.2)) + geom_hline(yintercept = ifelse(logged, 0, 1), col = "grey")
# p <- p + theme_plot(rotate.x = 90, lpos = "none")
# p <- p + scale_color_manual(values = CRS_colors)
# p <- p + facet_grid(rows = vars(factor(celline, names(color_celline))), scales = "free_y", drop = TRUE)
# p <- p + theme(axis.text.x = element_text(colour = rep(CRS_colors, e = 5)))
# print(p)
# dev.off()

# oligo_Dist_ct <- merge(oligo_activity_obs_ctc_archaic, TSS_dist, by = c("oligo", "celline", "crsID"), allow.cartesian = TRUE)
# oligo_Dist_ct_count <- oligo_Dist_ct[, .N, by = .(oligo_class_CRITERIA = gsub("strong silencer", "silencer", oligo_class_CRITERIA), celline, condition, dist_bin = cut(Dist_NearestTSS / 1000, dist_th_kb_5bin, include.lowest = TRUE))]
# oligo_Dist_ct_count[, Ntot_bin := sum(N), by = .(celline, condition, dist_bin)]
# oligo_Dist_ct_count[, as.list(c(Ntot_bin, unlist(binom.test(N, Ntot_bin)[c("estimate", "conf.int")]))), by = .(oligo_class_CRITERIA, celline, condition, dist_bin)]
# oligo_Dist_ct_Pct <- oligo_Dist_ct_count[, as.list(unlist(binom.test(N, Ntot_bin)[c("estimate", "conf.int")])), by = .(oligo_class_CRITERIA, celline, condition, dist_bin)]
# setnames(oligo_Dist_ct_Pct, c("estimate.probability of success", "conf.int1", "conf.int2"), c("Pct", "Pct_lo", "Pct_hi"), skip_absent = TRUE)


# logged <- TRUE
# CRS_colors <- c("strong enhancer" = "#0f50b8", "enhancer" = "#4287f5", "silencer" = "#0fb88e")
# pdf(sprintf("%s/08h_Enrichment_Dist5Bin_pointrange_by_ct__oligo_level__%s.pdf", ACTIVITY_DIR, CRITERIA_ACTIVE), height = 4, width = 3)
# p <- ggplot(oligo_Dist_ct_Pct[oligo_class_CRITERIA != "inactive" & condition == "NS"], aes(x = dist_bin, y = Pct, ymin = Pct_lo, ymax = Pct_hi, col = oligo_class_CRITERIA)) +
#   geom_pointrange()
# p <- p + ylab("Pct of oligos") + xlab("distance to nearest TSS")
# p <- p + geom_pointrange(size = .1, position = position_dodge(0.2)) + geom_hline(yintercept = ifelse(logged, 0, 1), col = "grey")
# p <- p + theme_plot(rotate.x = 90, lpos = "none")
# p <- p + scale_color_manual(values = CRS_colors)
# p <- p + facet_grid(rows = vars(factor(celline, names(color_celline))), scales = "free_y", drop = TRUE)
# print(p)
# dev.off()

cat("All done !")
q("no")

# pdf(sprintf("%s/10h_Enrichment_TFBS_among_CRS_class_by_celline_pointrange_%s.pdf", ACTIVITY_DIR, CRITERIA_ACTIVE), height = 4, width = 3)
# TF__shown <- TF_enrich_CRS_class[order(-OR_low)]
# TF__shown <- TF__shown[!duplicated(paste(names,celline))][FDR < .05, head(.SD, 5), by = .(type, celline)]
# TF__shown <- TF__shown[order(factor(type, levels = c("strong enhancer", "weak enhancer", "silencer")), -OR_low), ]
# CRS_colors <- c("strong enhancer" = "#0f50b8", "weak enhancer" = "#4287f5", "silencer" = "#0fb88e")
# p <- plot_TFenrich_class(TF_enrich_CRS_class[(TF_ID %in% TF__shown$TF_ID), ], logged = TRUE, TF__order = TF__shown$TF_ID, CRS_group = "active CRS")
# p <- p + theme(axis.text.x = element_text(colour = rep(CRS_colors, e = 5)))
# print(p)
# dev.off()


# pdf(sprintf("%s/10h_Enrichment_TFBS_among_CRS_class_by_celline_pointrange_%s.pdf", ACTIVITY_DIR, CRITERIA_ACTIVE), height = 4, width = 3)
# TF__shown <- TF_enrich_CRS_class[order(-OR_low)]
# TF__shown <- TF__shown[!duplicated(paste(names,celline))][FDR < .05, head(.SD, 5), by = .(type, celline)]
# TF__shown <- TF__shown[order(factor(type, levels = c("strong enhancer", "weak enhancer", "silencer")), -OR_low), ]
# CRS_colors <- c("strong enhancer" = "#0f50b8", "weak enhancer" = "#4287f5", "silencer" = "#0fb88e")
# p <- plot_TFenrich_class(TF_enrich_CRS_class[(TF_ID %in% TF__shown$TF_ID), ], logged = TRUE, TF__order = TF__shown$TF_ID, CRS_group = "active CRS")
# p <- p + theme(axis.text.x = element_text(colour = rep(CRS_colors, e = 5)))
# print(p)
# dev.off()


# oligo_activity_DT_TFBS=merge(oligo_activity_obs_ctc,TF_score_annot,by='oligo',all.x=TRUE,allow.cartesian=TRUE)
# oligo_activity_DT_TFBS[,sig_up:=grepl('enhancer',oligo_class_CRITERIA),by=.(matrix_id,names)]
# oligo_activity_DT_TFBS[,sig_down:=grepl('silencer',oligo_class_CRITERIA),by=.(matrix_id,names)]


# # check most active TFsn& plot enrichment
# TF_enrich_active_top10=oligo_activity_DT_TFBS[!is.na(matrix_id),as.list(unlist(fisher.test(table(sig_up,top10))[c('p.value','estimate','conf.int')])),by=.(analysis_name,matrix_id,names)]
# TF_enrich_active_top10[,COND_ID:=gsub('-ACE2','',gsub('_all','',analysis_name))]
# TF_enrich_active_top10[,TF_ID:=paste0(names," (",matrix_id,")")]
# setnames(TF_enrich_active_top10,c('estimate.odds ratio','conf.int1','conf.int2'),c('OR','OR_low','OR_high'),skip_absent = TRUE)


### removed (error management too complex)

# TF_enrich_active_top10=CRS_activity_DT_TFBS[!is.na(matrix_id),{
#     tab = table(sig_up,top10)
#     if(sum(top10,na.rm=T)>0 & sum(sig_up,na.rm=T)>0){
#       Fisher_test=unlist(fisher.test(tab)[c('p.value','estimate','conf.int')])
#       res =as.list(Fisher_test)
#     }else{
#       if(sum(top10,na.rm=T)==0){cat('no CRS with TF binding in condition\n')}
#       if(sum(sig_up,na.rm=T)==0){cat('no active CRS with in condition \n')}
#       res = as.list(c(-999,-999,-999,-999))
#     }
#     },by=.(analysis_name, matrix_id, names)]

# TF_enrich_repressor_top10=CRS_activity_DT_TFBS[!is.na(matrix_id),{
#     tab = table(down,top10)
#     if(sum(top10,na.rm=T)>0 & sum(down,na.rm=T)>0){
#       Fisher_test=unlist(fisher.test(tab)[c('p.value','estimate','conf.int')])
#       res =as.list(Fisher_test)
#     }else{
#       if(sum(top10,na.rm=T)==0){cat('no CRS with TF binding in condition\n')}
#       if(sum(sig_down,na.rm=T)==0){cat('no repressor CRS with in condition \n')}
#       res = as.list(c(-999,-999,-999,-999))
#     }
#     },by=.(analysis_name, matrix_id, names)]


# plot_TFenrich <- function(TF_enrich_table, logged = TRUE, TF__order = NULL, CRS_group = "active CRS (alpha>1, FDR<1%)") {
#   if (!is.null(TF__order)) {
#     TF_enrich_table[, TF_ID := factor(TF_ID, levels = TF__order, ordered = TRUE)]
#   }

#   if (logged == TRUE) {
#     p <- ggplot(TF_enrich_table, aes(x = TF_ID, y = log2(OR), ymin = log2(OR_low), ymax = log2(OR_high), col = COND_ID))
#     OR_stat <- "log2(OR)"
#   } else {
#     p <- ggplot(TF_enrich_table, aes(x = TF_ID, y = OR, ymin = OR_low, ymax = OR_high, col = COND_ID))
#     OR_stat <- "OR"
#   }
#   p <- p + ylab(sprintf("Enrichment of TFBS in %s \n %s", CRS_group, OR_stat)) + xlab("TF")
#   p <- p + geom_pointrange(size = .1) + geom_hline(yintercept = ifelse(logged, 0, 1), col = "grey")

#   p <- p + theme_plot(rotate.x = 90, lpos = "none") + coord_flip()
#   p <- p + scale_color_manual(values = color_setup_simplified_norep)
#   p <- p + facet_grid(cols = vars(COND_ID), scales = "free_y", drop = TRUE)
#   return(p)
# }

# plot_TFenrich_byTF <- function(TF_enrich_table, logged = TRUE, TF__order = NULL, CRS_group = "active CRS (alpha>1, FDR<1%)") {
#   if (logged == TRUE) {
#     p <- ggplot(TF_enrich_table, aes(x = COND_ID, y = log2(OR), ymin = log2(OR_low), ymax = log2(OR_high), col = COND_ID))
#     OR_stat <- "log2(OR)"
#   } else {
#     p <- ggplot(TF_enrich_table, aes(x = COND_ID, y = OR, ymin = OR_low, ymax = OR_high, col = COND_ID))
#     OR_stat <- "OR"
#   }
#   p <- p + ylab(sprintf("Enrichment of TFBS in %s \n %s", CRS_group, OR_stat)) + xlab("TF")
#   p <- p + geom_pointrange(size = .1) + geom_hline(yintercept = ifelse(logged, 0, 1), col = "grey")
#   p <- p + theme_plot(rotate.x = 90, lpos = "none")
#   p <- p + scale_color_manual(values = color_setup_simplified_norep)
#   if (!is.null(TF__order)) {
#     p <- p + facet_wrap(~ factor(TF_ID, TF__order), scales = "free_y", drop = TRUE)
#   } else {
#     p <- p + facet_wrap(~TF_ID, scales = "free_y", drop = TRUE)
#   }
#   return(p)
# }


# plot_TFenrich_cellline <- function(TF_enrich_table, logged = TRUE, TF__order = NULL, CRS_group = "celltype-specific CRS (FDR<1%)") {
#   if (!is.null(TF__order)) {
#     TF_enrich_table[, TF_ID := factor(TF_ID, levels = TF__order, ordered = TRUE)]
#   }

#   if (logged == TRUE) {
#     p <- ggplot(TF_enrich_table, aes(x = TF_ID, y = log2(OR), ymin = log2(OR_low), ymax = log2(OR_high), col = celltype))
#     OR_stat <- "log2(OR)"
#   } else {
#     p <- ggplot(TF_enrich_table, aes(x = TF_ID, y = OR, ymin = OR_low, ymax = OR_high, col = celltype))
#     OR_stat <- "OR"
#   }
#   p <- p + ylab(sprintf("Enrichment of TFBS in %s \n %s", CRS_group, OR_stat)) + xlab("TF")
#   p <- p + geom_pointrange(size = .1) + geom_hline(yintercept = ifelse(logged, 0, 1), col = "grey")

#   p <- p + theme_plot(rotate.x = 90, lpos = "none") + coord_flip()
#   p <- p + scale_color_manual(values = color_celline[c("A549", "HepG2", "K562")])
#   p <- p + facet_grid(cols = vars(celltype), rows = vars(factor(type, c("up", "down"))), scales = "free_y", drop = TRUE)
#   return(p)
# }



# plot_TFenrich_class <- function(TF_enrich_table, logged = TRUE, TF__order = NULL, CRS_group = "") {
#   if (!is.null(TF__order)) {
#     TF_enrich_table[, TF_ID := factor(TF_ID, levels = TF__order, ordered = TRUE)]
#   }
#   if (logged == TRUE) {
#     p <- ggplot(TF_enrich_table, aes(x = TF_ID, y = log2(OR), ymin = log2(OR_low), ymax = log2(OR_high), col = type))
#     OR_stat <- "log2(OR)"
#   } else {
#     p <- ggplot(TF_enrich_table, aes(x = TF_ID, y = OR, ymin = OR_low, ymax = OR_high, col = type))
#     OR_stat <- "OR"
#   }
#   p <- p + ylab(sprintf("Enrichment of TFBS in %s \n %s", CRS_group, OR_stat)) + xlab("TF")
#   p <- p + geom_pointrange(size = .1) + geom_hline(yintercept = ifelse(logged, 0, 1), col = "grey")
#   p <- p + theme_plot(rotate.x = 90, lpos = "none")
#   p <- p + scale_color_manual(values = CRS_colors)
#   p <- p + facet_grid(rows = vars(factor(type, c("strong enhancer", "weak enhancer", "silencer"))), scales = "free_y", drop = TRUE)
#   return(p)
# }

# pdf(sprintf("%s/10e_Enrichment_TFBS_among_active_CRS_per_condition_pointrange_%s_V1_noGCadj.pdf", FIGURE_DIR, "FDR1_norm"), height = 3, width = 7)
# TF__shown <- TF_enrich_active_top10[order(-OR_low), ][, head(.SD, 10), by = COND_ID][, unique(TF_ID)]
# p <- plot_TFenrich(TF_enrich_active_top10[TF_ID %in% TF__shown])
# print(p)
# dev.off()

# results are suspiciously driven by GC content, so we adjust scores on GC content of CRS and repeat enrichment

# TF_enrich_active_top10GC=CRS_activity_DT_TFBS[!is.na(matrix_id),as.list(unlist(fisher.test(table(sig_up,top10_GC))[c('p.value','estimate','conf.int')])),by=.(analysis_name,matrix_id,names)]
# TF_enrich_active_top10GC=merge(TF_enrich_active_top10GC,condition_summary,by='analysis_name')
# TF_enrich_active_top10GC[,TF_ID:=paste0(names," (",matrix_id,")")]
# setnames(TF_enrich_active_top10GC,c('estimate.odds ratio','conf.int1','conf.int2'),c('OR','OR_low','OR_high'),skip_absent = TRUE)
# TF_enrich_active_top10GC[,Cell_specificity:=summary(lm((log2(OR)/(log2(OR_high)-log2(OR_low))*3.92)~celline,data=.SD))$r.squared,by=.(TF_ID)]
# TF_enrich_active_top10GC[order(-Cell_specificity),]

# TF_enrich_active_top10GC[,Cell_specificity:=summary(lm((OR/(OR_high-OR_low)*3.92)~celline,data=.SD))$r.squared,by=.(TF_ID)]
# TF_enrich_active_top10GC[,mean_HepG2:=mean((log2(OR)/(log2(OR_high)-log2(OR_low))*3.92)[celline=='HepG2']),by=.(TF_ID)]
# TF_enrich_active_top10GC[,mean_A549:=mean((log2(OR)/(log2(OR_high)-log2(OR_low))*3.92)[celline=='A549']),by=.(TF_ID)]
# TF_enrich_active_top10GC[,mean_K562:=mean((log2(OR)/(log2(OR_high)-log2(OR_low))*3.92)[celline=='K562']),by=.(TF_ID)]

# pdf(sprintf('%s/10e_Enrichment_TFBS_among_active_CRS_per_condition_pointrange_%s.pdf',FIGURE_DIR,'FDR1_norm'),height=3,width=7)
# TF__shown=TF_enrich_active_top10GC[order(-OR_low),][,head(.SD,10),by=COND_ID][,unique(TF_ID)]
# p <- plot_TFenrich(TF_enrich_active_top10GC[TF_ID%in%TF__shown],logged=FALSE,TF__shown)
# print(p)
# dev.off()

# pdf(sprintf('%s/10e_Enrichment_TFBS_among_active_CRS_per_condition_pointrange_%s_cellSpecific.pdf',FIGURE_DIR,'FDR1_norm'),height=5,width=7)
# # TF__shown_K562=TF_enrich_active_top10GC[mean_K562>0,][order(-Cell_specificity),unique(TF_ID)][1:4]
# # TF__shown_HepG2=TF_enrich_active_top10GC[mean_HepG2>0,][order(-Cell_specificity),unique(TF_ID)][1:4]
# # TF__shown_A549=TF_enrich_active_top10GC[mean_A549>0,][order(-Cell_specificity),unique(TF_ID)][1:4]
# TF__shown_K562=TF_enrich_active_top10GC[mean_K562>0,][order(-mean_K562+mean_HepG2+mean_A549),unique(TF_ID)][1:4]
# TF__shown_HepG2=TF_enrich_active_top10GC[mean_HepG2>0,][order(-mean_HepG2+mean_K562+mean_A549),unique(TF_ID)][1:4]
# TF__shown_A549=TF_enrich_active_top10GC[mean_A549>0,][order(-mean_A549+mean_HepG2+mean_K562),unique(TF_ID)][1:4]

# TF__shown=c(TF__shown_K562,TF__shown_HepG2,TF__shown_A549)
# p <- plot_TFenrich_byTF(TF_enrich_active_top10GC[TF_ID%in%TF__shown],logged=FALSE,TF__order=TF__shown)
# print(p)
# dev.off()

# ##### TF associated with silencer CRS #####
# TF_enrich_repressor_top10GC=CRS_activity_DT_TFBS[!is.na(matrix_id),as.list(unlist(fisher.test(table(sig_down,top10_GC))[c('p.value','estimate','conf.int')])),by=.(analysis_name,matrix_id,names)]
# TF_enrich_repressor_top10GC=merge(TF_enrich_repressor_top10GC,condition_summary,by='analysis_name')
# TF_enrich_repressor_top10GC[,TF_ID:=paste0(names," (",matrix_id,")")]
# setnames(TF_enrich_repressor_top10GC,c('estimate.odds ratio','conf.int1','conf.int2'),c('OR','OR_low','OR_high'),skip_absent = TRUE)
# # TF_enrich_repressor_top10GC[,Cell_specificity:=summary(lm((log2(OR)/(log2(OR_high)-log2(OR_low))*3.92)~celline,data=.SD))$r.squared,by=.(TF_ID)]
# # TF_enrich_repressor_top10GC[order(-Cell_specificity),]

# # TF_enrich_repressor_top10GC[,Cell_specificity:=summary(lm((OR/(OR_high-OR_low)*3.92)~celline,data=.SD))$r.squared,by=.(TF_ID)]
# TF_enrich_repressor_top10GC[,mean_HepG2:=mean((log2(OR)/(log2(OR_high)-log2(OR_low))*3.92)[celline=='HepG2']),by=.(TF_ID)]
# TF_enrich_repressor_top10GC[,mean_A549:=mean((log2(OR)/(log2(OR_high)-log2(OR_low))*3.92)[celline=='A549']),by=.(TF_ID)]
# TF_enrich_repressor_top10GC[,mean_K562:=mean((log2(OR)/(log2(OR_high)-log2(OR_low))*3.92)[celline=='K562']),by=.(TF_ID)]

# pdf(sprintf('%s/10f_Enrichment_TFBS_among_repressor_CRS_per_condition_pointrange_%s.pdf',FIGURE_DIR,'FDR1_norm'),height=3,width=7)
# TF__shown=TF_enrich_repressor_top10GC[order(-OR_low),][,head(.SD,10),by=COND_ID][,unique(TF_ID)]
# p <- plot_TFenrich(TF_enrich_repressor_top10GC[TF_ID%in%TF__shown],logged=FALSE,TF__shown,CRS_group="silencer CRS (alpha<1, FDR<1%)")
# print(p)
# dev.off()

# pdf(sprintf('%s/10f_Enrichment_TFBS_among_repressor_CRS_per_condition_pointrange_%s_cellSpecific.pdf',FIGURE_DIR,'FDR1_norm'),height=5,width=7)
# # TF__shown_K562=TF_enrich_active_top10GC[mean_K562>0,][order(-Cell_specificity),unique(TF_ID)][1:4]
# # TF__shown_HepG2=TF_enrich_active_top10GC[mean_HepG2>0,][order(-Cell_specificity),unique(TF_ID)][1:4]
# # TF__shown_A549=TF_enrich_active_top10GC[mean_A549>0,][order(-Cell_specificity),unique(TF_ID)][1:4]
# TF__shown_K562=TF_enrich_repressor_top10GC[mean_K562>0,][order(-mean_K562+mean_HepG2+mean_A549),unique(TF_ID)][1:4]
# TF__shown_HepG2=TF_enrich_repressor_top10GC[mean_HepG2>0,][order(-mean_HepG2+mean_K562+mean_A549),unique(TF_ID)][1:4]
# TF__shown_A549=TF_enrich_repressor_top10GC[mean_A549>0,][order(-mean_A549+mean_HepG2+mean_K562),unique(TF_ID)][1:4]

# TF__shown=c(TF__shown_K562,TF__shown_HepG2,TF__shown_A549)
# p <- plot_TFenrich_byTF(TF_enrich_repressor_top10GC[TF_ID%in%TF__shown], logged=FALSE, TF__order=TF__shown,CRS_group="silencer CRS (alpha<1, FDR<1%)")
# print(p)
# dev.off()


# TF_enrich_enhancer_any=CRS_activity_DT_TFBS[!is.na(matrix_id),as.list(unlist(fisher.test(table(unique(CRS)%in%CRS[sig_up],unique(CRS)%in%CRS[top10_GC]))[c('p.value','estimate','conf.int')])),by=.(matrix_id,names)]
# TF_enrich_silencer_any=CRS_activity_DT_TFBS[!is.na(matrix_id),as.list(unlist(fisher.test(table(unique(CRS)%in%CRS[sig_down],unique(CRS)%in%CRS[top10_GC]))[c('p.value','estimate','conf.int')])),by=.(matrix_id,names)]
# TF_enrich_any=list(enhancer=TF_enrich_enhancer_any,silencer=TF_enrich_silencer_any)
# TF_enrich_any=rbindlist(TF_enrich_any,idcol='type')
# setnames(TF_enrich_any,c('estimate.odds ratio','conf.int1','conf.int2'),c('OR','OR_low','OR_high'),skip_absent = TRUE)
# TF_enrich_any[,TF_ID:=paste0(names," (",matrix_id,")")]

# mean_activity_CRS=CRS_activity_DT_TFBS[,.(TF_affinity=mean(score_GCadj),
#                                           mean_CRS_activity_allconditions=mean(alpha.GC_norm),
#                                           mean_CRS_activity_HepG2=mean(alpha.GC_norm[celline=='HepG2']),
#                                           mean_CRS_activity_K562=mean(alpha.GC_norm[celline=='K562']),
#                                           mean_CRS_activity_A549=mean(alpha.GC_norm[celline=='A549'])
#                                           ),by=.(CRS,matrix_id,names)]

# #corTF_activity=CRS_activity_DT_TFBS[,.(cor_TF_activity=cor(TF_affinity,mean_CRS_activity_allconditions,use='p')),by=.(matrix_id,names)]
# lasso_data=dcast(mean_activity_CRS,CRS+mean_CRS_activity_allconditions+mean_CRS_activity_K562+mean_CRS_activity_A549+mean_CRS_activity_HepG2~matrix_id+names,value.var='TF_affinity')

# lasso_data_y=lasso_data[,mean_CRS_activity_allconditions]
# lasso_data_x=lasso_data[,mget(colnames(lasso_data)[grepl('^MA',colnames(lasso_data))])]
# run_lasso=function(lasso_data_y,lasso_data_x,suffix='_all'){
#   w.na=which(apply(is.na(lasso_data_x),1,any))
#   cv_lasso=cv.glmnet(as.matrix(lasso_data_x[-w.na,]), lasso_data_y[-w.na], alpha=1)
#   lasso_fit=glmnet(as.matrix(lasso_data_x[-w.na,]), lasso_data_y[-w.na], alpha=1)
#   coef.min=coef(lasso_fit, s = cv_lasso$lambda.min)
#   coef.1se=coef(lasso_fit, s = cv_lasso$lambda.1se)
#   DT_lasso=data.table(rownames(coef.1se),coeff.1se=as.matrix(coef.1se),coeff.min=as.matrix(coef.min))
#   setnames(DT_lasso,c('V1','coeff.1se.s1','coeff.min.s1'),c('TF.matrix',paste0(c('coeff.1se','coeff.min'),suffix)))
#   DT_lasso
# }
# DT_lasso_all=run_lasso(lasso_data[,mean_CRS_activity_allconditions],lasso_data_x,'.all')
# DT_lasso_A549=run_lasso(lasso_data[,mean_CRS_activity_A549],lasso_data_x,'.A549')
# DT_lasso_HepG2=run_lasso(lasso_data[,mean_CRS_activity_HepG2],lasso_data_x,'.HepG2')
# DT_lasso_K562=run_lasso(lasso_data[,mean_CRS_activity_K562],lasso_data_x,'.K562')

# DT_lasso_all=merge(DT_lasso_all,DT_lasso_A549,by='TF.matrix')
# DT_lasso_all=merge(DT_lasso_all,DT_lasso_HepG2,by='TF.matrix')
# DT_lasso_all=merge(DT_lasso_all,DT_lasso_K562,by='TF.matrix')

# fwrite(DT_lasso_all,file=sprintf('%s/lasso_coefficients_CRS_activity_vs_TF_bycelltype.txt',FIGURE_DIR),sep='\t')
