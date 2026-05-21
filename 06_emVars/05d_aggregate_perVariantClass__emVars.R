MPRA_DIR <- "/pasteur/helix/projects/evo_immuno_pop/MPRA"
ANALYSIS_DIR <- "MPRA_count_exp6_analysis3"

source(sprintf("%s/scripts/%s/00_starter.R", MPRA_DIR, ANALYSIS_DIR))
source(sprintf("%s/scripts/%s/00_load_associations.R", MPRA_DIR, ANALYSIS_DIR))
source(sprintf("%s/scripts/misc_plots.R", MPRA_DIR))

RUN_ID <- "RUN2_Z2_nBC10"
CRITERIA_ACTIVE <- "lfdr20_scrambled1pct_GCnorm"
FILTER_ACTIVE <- TRUE
CRITERIA_EMVARS <- "EmVar_lfdr20_FC.2"
EMVAR_LIST <- 'significant'
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
  if (cmd[i] == "--emVar_list" || cmd[i] == "-l") {
    EMVAR_LIST <- cmd[i + 1] # either 'suggestive' or 'significant'
  }
}

if(!EMVAR_LIST %in% c('suggestive','significant')) {
  cat("WARNING: --emVar_list must be either 'suggestive' or 'significant', setting to 'significant' and continuing\n")
  EMVAR_LIST <- 'significant'
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
#FIGURE_DIR <- sprintf("%s/figures/%s/04_emVars/%s/%s", MPRA_DIR, ANALYSIS_DIR, CRITERIA_EMVARS, CRITERIA_ACTIVE_OUT)

dir.create(FIGURE_DIR, recursive = TRUE)
dir.create(sprintf("%s/SupTables/", FIGURE_DIR), recursive = TRUE)


# load tested oligos
source(sprintf("%s/scripts/%s/02z__define_includedCRS.R", MPRA_DIR, ANALYSIS_DIR))
dir.create(sprintf("%s/Pct_emVar_by_pop/", FIGURE_DIR))
dir.create(sprintf("%s/Pct_emVar_by_pop_x_adaptive_status/", FIGURE_DIR))
dir.create(sprintf("%s/Pct_emVar_by_source/", FIGURE_DIR))
dir.create(sprintf("%s/Pct_emVar_by_match/", FIGURE_DIR))
dir.create(sprintf("%s/Pct_emVar_by_scenario/", FIGURE_DIR))

# emVARs_all_results <- fread(sprintf("%s/all_emVars_results.tsv.gz", IN_DIR), sep = "\t")

# oligo_activity_file <- sprintf("%s/oligo_activity__all.tsv.gz", OUT_DIR)
# oligo_activity <- fread(oligo_activity_file)
# oligo_activity_obs <- fread(sprintf("%s/all_oligos_annotated.tsv.gz", OUT_DIR)) # where does this come from	?


########## load emVar data  ##########
emVARs_obs_ctc <- fread(file = sprintf("%s/all_emVars_annotated_celltype__%s.tsv", EMVAR_DIR, CRITERIA_ACTIVE_OUT), sep = "\t")

emVARs_obs_ctc_full <- fread(sprintf("%s/all_emVars_annotated_full_celltype__%s.tsv", EMVAR_DIR, CRITERIA_ACTIVE_OUT))

## ARCHAIC ALLELES ARE ASSOCIATED WITH LOWER ENHANCERS ACTIVITY AND STRONGER REPRESSOR ACTIVITY ?

############################################################################################################
################################## Effect of archaic alleles on expression: POPULATION #####################
############################################################################################################

minSNP <- 5
tic("effects of archaic alleles on expression: by Population\n")
cat("\nby Population  (all conditions merged)\n")

Pct_emVar_byPop <- get_Pct_emVars(emVARs_obs_ctc_full, split_by = c("POP_adaptive"), total_by = NULL, emVarlist = EMVAR_LIST)
Pct_emVar_byPop <- make_it_long(Pct_emVar_byPop, split_by = "POP_adaptive")
Pct_emVar_byPop[N_test > minSNP, FDR := p.adjust(Pvalue, "fdr"), by = measure]
Pct_emVar_byPop[N_test > minSNP][order(Pvalue), head(.SD, 5), by = measure]
Pct_emVar_byPop[, POP_adaptive := toupper(gsub("SGDP_", "", POP_adaptive))]
fwrite(Pct_emVar_byPop, sprintf("%s/Pct_emVar_by_pop/01a_Pct_emVars_%s_by_population.tsv", FIGURE_DIR, EMVAR_LIST), sep = "\t")

for (i_MEASURE in seq_len(Measure_table[, .N])) {
  MEASURE <- Measure_table[i_MEASURE, measure]
  MEASURE_LABEL <- Measure_table[i_MEASURE, measure_label]
  MEASURE_NULL <- Measure_table[i_MEASURE, null]

  FigData <- Pct_emVar_byPop[measure == MEASURE & N_test > minSNP, ]
  p <- ggplot(FigData)
  p <- p + geom_hline(yintercept = MEASURE_NULL, col = "lightgrey", linetype = "dashed")
  p <- p + geom_pointrange(aes(x = factor(POP_adaptive, names(color_populations_MPRA)), y = Pct * 100, ymin = Pct_lo * 100, ymax = Pct_hi * 100, col = POP_adaptive), size = .3)
  # p <- p + facet_grid(rows=vars(measure)) + ylab("Percentage of alleles") + xlab("Population")
  if (is.na(MEASURE_NULL)) {
    MEASURE_NULL <- as.numeric(Pct_emVar_byPop[measure == MEASURE, .(Pct = median(Pct))])
  }
  p <- p + geom_hline(yintercept = MEASURE_NULL * 100, linetype = 2)
  p <- p + theme_plot(rotate.x = 90) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + guides(col = "none") + xlab("")
  p <- p + scale_color_manual(values = color_populations_MPRA)
  p <- p + ylab(MEASURE_LABEL) + xlab("Population")
  pdf(sprintf("%s/Pct_emVar_by_pop/01a_Pct_%s_%s_by_population.pdf", FIGURE_DIR, MEASURE, EMVAR_LIST), height = 2.3, width = 1.5)
  print(p)
  dev.off()
}

cat("\nby Population  (each condition separately\n")
Pct_emVar_byPop <- get_Pct_emVars(emVARs_obs_ctc_full, split_by = c("POP_adaptive", "COND_ID"), total_by = c("COND_ID"), emVarlist = EMVAR_LIST)
Pct_emVar_byPop <- make_it_long(Pct_emVar_byPop, split_by = c("POP_adaptive", "COND_ID"))
Pct_emVar_byPop[N_test > minSNP, FDR := p.adjust(Pvalue, "fdr"), by = measure]
Pct_emVar_byPop[N_test > minSNP][order(Pvalue), head(.SD, 5), by = measure]
Pct_emVar_byPop[, POP_adaptive := toupper(gsub("SGDP_", "", POP_adaptive))]
Pct_emVar_byPop[, COND_ID := factor(COND_ID, levels = names(color_setup_simplified_norep))]

fwrite(Pct_emVar_byPop, sprintf("%s/Pct_emVar_by_pop/01b_Pct_emVars_%s_by_population_byCondition.tsv", FIGURE_DIR, EMVAR_LIST), sep = "\t")

for (i_MEASURE in seq_len(Measure_table[, .N])) {
  MEASURE <- Measure_table[i_MEASURE, measure]
  MEASURE_LABEL <- Measure_table[i_MEASURE, measure_label]
  MEASURE_NULL <- Measure_table[i_MEASURE, null]

  FigData <- Pct_emVar_byPop[measure == MEASURE & N_test > minSNP, ]
  p <- ggplot(FigData)
  p <- p + geom_hline(yintercept = MEASURE_NULL, col = "lightgrey", linetype = "dashed")
  p <- p + geom_pointrange(aes(x = factor(POP_adaptive, names(color_populations_MPRA)), y = Pct * 100, ymin = Pct_lo * 100, ymax = Pct_hi * 100, col = COND_ID, alpha = Pvalue < 0.01), size = .3)
	p <- p + geom_pointrange(data = FigData[Pvalue < 0.01, ], aes(x = factor(POP_adaptive, names(color_populations_MPRA)), y = Pct * 100, ymin = Pct_lo * 100, ymax = Pct_hi * 100), col = "red", size = .3, alpha = .5)
  # p <- p + facet_grid(rows=vars(measure)) + ylab("Percentage of alleles") + xlab("Population")
  if (is.na(MEASURE_NULL)) {
    MEASURE_NULL <- Pct_emVar_byPop[measure == MEASURE, .(Pct = median(Pct)), by = COND_ID]
    p <- p + geom_hline(data = MEASURE_NULL, aes(yintercept = Pct * 100), linetype = 2)
  } else {
    p <- p + geom_hline(yintercept = MEASURE_NULL * 100, linetype = 2)
  }
  p <- p + theme_plot(rotate.x = 90) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + guides(col = "none") + xlab("")
  p <- p + scale_color_manual(values = color_setup_simplified_norep)
  p <- p + scale_alpha_manual(values = alpha_TRUEFALSE)
  p <- p + ylab(MEASURE_LABEL) + xlab("Population")
	p <- p + facet_wrap(~factor(COND_ID,names(color_setup_simplified_norep)),ncol=4)
  p <- p + theme(legend.box = "vertical")
  pdf(sprintf("%s/Pct_emVar_by_pop/01b_Pct_%s_%s_by_population_byCondition.pdf", FIGURE_DIR, MEASURE, EMVAR_LIST), height = 4, width = 4)
  print(p)
  dev.off()
}

cat("\nby Population  (each celltype separately\n")
Pct_emVar_byPop <- get_Pct_emVars(emVARs_obs_ctc_full[condition == "NS", ], split_by = c("POP_adaptive", "celltype"), total_by = c("celltype"), emVarlist = EMVAR_LIST)
Pct_emVar_byPop <- make_it_long(Pct_emVar_byPop, split_by = c("POP_adaptive", "celltype"))
Pct_emVar_byPop[N_test > minSNP, FDR := p.adjust(Pvalue, "fdr"), by = measure]
Pct_emVar_byPop[N_test > minSNP][order(Pvalue), head(.SD, 5), by = measure]
Pct_emVar_byPop[, POP_adaptive := toupper(gsub("SGDP_", "", POP_adaptive))]
fwrite(Pct_emVar_byPop, sprintf("%s/Pct_emVar_by_pop/01c_Pct_emVars_%s_by_population_byCelltypeNS.tsv", FIGURE_DIR, EMVAR_LIST), sep = "\t")

for (i_MEASURE in seq_len(Measure_table[, .N])) {
  MEASURE <- Measure_table[i_MEASURE, measure]
  MEASURE_LABEL <- Measure_table[i_MEASURE, measure_label]
  MEASURE_NULL <- Measure_table[i_MEASURE, null]

  FigData <- Pct_emVar_byPop[measure == MEASURE & N_test > minSNP, ]
  p <- ggplot(FigData)
  p <- p + geom_hline(yintercept = MEASURE_NULL, col = "lightgrey", linetype = "dashed")
  p <- p + geom_pointrange(aes(x = factor(POP_adaptive, names(color_populations_MPRA)), y = Pct * 100, ymin = Pct_lo * 100, ymax = Pct_hi * 100, col = celltype), size = .3, alpha = .5)
  p <- p + geom_pointrange(data = FigData[Pvalue < 0.01, ], aes(x = factor(POP_adaptive, names(color_populations_MPRA)), y = Pct * 100, ymin = Pct_lo * 100, ymax = Pct_hi * 100), col = "red", size = .3, alpha = .5)
  # p <- p + facet_grid(rows=vars(measure)) + ylab("Percentage of alleles") + xlab("Population")
  if (is.na(MEASURE_NULL)) {
    MEASURE_NULL <- Pct_emVar_byPop[measure == MEASURE, .(Pct = median(Pct)), by = celltype]
    p <- p + geom_hline(data = MEASURE_NULL, aes(yintercept = Pct * 100), linetype = 2)
  } else {
    p <- p + geom_hline(yintercept = MEASURE_NULL * 100, linetype = 2)
  }
  p <- p + theme_plot(rotate.x = 90) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + guides(col = "none") + xlab("")
  p <- p + scale_color_manual(values = color_celline)
  # p <- p + scale_alpha_manual(values=alpha_TRUEFALSE)
  p <- p + ylab(MEASURE_LABEL) + xlab("Population")
  p <- p + facet_grid(rows = vars(celltype))
  p <- p + theme(legend.box = "vertical")
  pdf(sprintf("%s/Pct_emVar_by_pop/01c_Pct_%s_%s_by_population_byCelltypeNS.pdf", FIGURE_DIR, MEASURE, EMVAR_LIST), height = 2.3, width = 1.7)
  print(p)
  dev.off()
}
toc()


############################################################################################################
################################## Effect of archaic alleles on expression: ARCHAIC SOURCE  ################
############################################################################################################

color_archaic_MPRA["undetermined"] <- "lightgrey"

tic("effects of archaic alleles on expression: by Archaic source\n")
cat("\nby Archaic source  (all conditions merged)\n")

Pct_emVar_bySource <- get_Pct_emVars(emVARs_obs_ctc[Introgression_source_top != "", ], split_by = c("Introgression_source_top"), total_by = NULL, emVarlist = EMVAR_LIST)
Pct_emVar_bySource <- make_it_long(Pct_emVar_bySource, split_by = "Introgression_source_top")
Pct_emVar_bySource[N_test > minSNP, FDR := p.adjust(Pvalue, "fdr"), by = measure]
Pct_emVar_bySource[N_test > minSNP][order(Pvalue), head(.SD, 5), by = measure]
fwrite(Pct_emVar_bySource, sprintf("%s/Pct_emVar_by_source/02a_Pct_emVars_%s_by_source.tsv", FIGURE_DIR, EMVAR_LIST), sep = "\t")

for (i_MEASURE in seq_len(Measure_table[, .N])) {
  MEASURE <- Measure_table[i_MEASURE, measure]
  MEASURE_LABEL <- Measure_table[i_MEASURE, measure_label]
  MEASURE_NULL <- Measure_table[i_MEASURE, null]

  FigData <- Pct_emVar_bySource[measure == MEASURE & N_test > minSNP, ]
  p <- ggplot(FigData)
  p <- p + geom_hline(yintercept = MEASURE_NULL, col = "lightgrey", linetype = "dashed")
  p <- p + geom_pointrange(aes(x = factor(Introgression_source_top, names(color_archaic_MPRA)), y = Pct * 100, ymin = Pct_lo * 100, ymax = Pct_hi * 100, col = Introgression_source_top), size = .5)
  # p <- p + facet_grid(rows=vars(measure)) + ylab("Percentage of alleles") + xlab("Archaic source")
  if (is.na(MEASURE_NULL)) {
    MEASURE_NULL <- as.numeric(Pct_emVar_bySource[measure == MEASURE, .(Pct = median(Pct))])
  }
  p <- p + geom_hline(yintercept = MEASURE_NULL * 100, linetype = 2)
  p <- p + theme_plot(rotate.x = 90) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + guides(col = "none") + xlab("")
  p <- p + scale_color_manual(values = color_archaic_MPRA)
  p <- p + ylab(MEASURE_LABEL) + xlab("Introgression source")
  pdf(sprintf("%s/Pct_emVar_by_source/02a_Pct_%s_%s_by_source.pdf", FIGURE_DIR, MEASURE, EMVAR_LIST), height = 2.6, width = 1.3)
  print(p)
  dev.off()
}

cat("\nby Archaic source  (each condition separately\n")
Pct_emVar_bySource <- get_Pct_emVars(emVARs_obs_ctc[Introgression_source_top != "", ], split_by = c("Introgression_source_top", "COND_ID"), total_by = c("COND_ID"), emVarlist = EMVAR_LIST)
Pct_emVar_bySource <- make_it_long(Pct_emVar_bySource, split_by = c("Introgression_source_top", "COND_ID"))
Pct_emVar_bySource[N_test > minSNP, FDR := p.adjust(Pvalue, "fdr"), by = measure]
Pct_emVar_bySource[N_test > minSNP][order(Pvalue), head(.SD, 5), by = measure]
Pct_emVar_bySource[, COND_ID := factor(COND_ID, levels = names(color_setup_simplified_norep))]

fwrite(Pct_emVar_bySource, sprintf("%s/Pct_emVar_by_source/02b_Pct_emVars_%s_by_source_byCondition.tsv", FIGURE_DIR, EMVAR_LIST), sep = "\t")

for (i_MEASURE in seq_len(Measure_table[, .N])) {
  MEASURE <- Measure_table[i_MEASURE, measure]
  MEASURE_LABEL <- Measure_table[i_MEASURE, measure_label]
  MEASURE_NULL <- Measure_table[i_MEASURE, null]

  FigData <- Pct_emVar_bySource[measure == MEASURE & N_test > minSNP, ]
  p <- ggplot(FigData)
  p <- p + geom_hline(yintercept = MEASURE_NULL, col = "lightgrey", linetype = "dashed")
  p <- p + geom_pointrange(aes(x = factor(Introgression_source_top, names(color_archaic_MPRA)), y = Pct * 100, ymin = Pct_lo * 100, ymax = Pct_hi * 100, col = COND_ID, alpha = Pvalue < 0.01), size = .3)
	p <- p + geom_pointrange(data = FigData[Pvalue < 0.01, ], aes(x = factor(Introgression_source_top, names(color_archaic_MPRA)), y = Pct * 100, ymin = Pct_lo * 100, ymax = Pct_hi * 100), col = "red", size = .3, alpha = .5)
  
  # p <- p + facet_grid(rows=vars(measure)) + ylab("Percentage of alleles") + xlab("Archaic source")
  if (is.na(MEASURE_NULL)) {
    MEASURE_NULL <- Pct_emVar_bySource[measure == MEASURE, .(Pct = median(Pct)), by = COND_ID]
    p <- p + geom_hline(data = MEASURE_NULL, aes(yintercept = Pct * 100), linetype = 2)
  } else {
    p <- p + geom_hline(yintercept = MEASURE_NULL * 100, linetype = 2)
  }
  p <- p + theme_plot(rotate.x = 90) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + guides(col = "none") + xlab("")
  p <- p + scale_color_manual(values = color_setup_simplified_norep)
  p <- p + scale_alpha_manual(values = alpha_TRUEFALSE)
  p <- p + ylab(MEASURE_LABEL) + xlab("Introgression source")
  p <- p + facet_wrap(~factor(COND_ID,names(color_setup_simplified_norep)),ncol=4)
  p <- p + theme(legend.box = "vertical")
  pdf(sprintf("%s/Pct_emVar_by_source/02b_Pct_%s_%s_by_source_byCondition.pdf", FIGURE_DIR, MEASURE, EMVAR_LIST), height = 4, width = 4)
  print(p)
  dev.off()
}

cat("\nby Archaic source  (each celltype separately\n")
Pct_emVar_bySource <- get_Pct_emVars(emVARs_obs_ctc[condition == "NS" & Introgression_source_top != "", ], split_by = c("Introgression_source_top", "celltype"), total_by = c("celltype"), emVarlist = EMVAR_LIST)
Pct_emVar_bySource <- make_it_long(Pct_emVar_bySource, split_by = c("Introgression_source_top", "celltype"))
Pct_emVar_bySource[N_test > minSNP, FDR := p.adjust(Pvalue, "fdr"), by = measure]
Pct_emVar_bySource[N_test > minSNP][order(Pvalue), head(.SD, 5), by = measure]
fwrite(Pct_emVar_bySource, sprintf("%s/Pct_emVar_by_source/02c_Pct_emVars_%s_by_source_byCelltypeNS.tsv", FIGURE_DIR, EMVAR_LIST), sep = "\t")

for (i_MEASURE in seq_len(Measure_table[, .N])) {
  MEASURE <- Measure_table[i_MEASURE, measure]
  MEASURE_LABEL <- Measure_table[i_MEASURE, measure_label]
  MEASURE_NULL <- Measure_table[i_MEASURE, null]

  FigData <- Pct_emVar_bySource[measure == MEASURE & N_test > minSNP, ]
  p <- ggplot(FigData)
  p <- p + geom_hline(yintercept = MEASURE_NULL, col = "lightgrey", linetype = "dashed")
  p <- p + geom_pointrange(aes(x = factor(Introgression_source_top, names(color_archaic_MPRA)), y = Pct * 100, ymin = Pct_lo * 100, ymax = Pct_hi * 100, col = celltype), size = .3, alpha = .5)
  p <- p + geom_pointrange(data = FigData[Pvalue < 0.01, ], aes(x = factor(Introgression_source_top, names(color_archaic_MPRA)), y = Pct * 100, ymin = Pct_lo * 100, ymax = Pct_hi * 100), col = "red", size = .3, alpha = .5)
  # p <- p + facet_grid(rows=vars(measure)) + ylab("Percentage of alleles") + xlab("Archaic source")
  if (is.na(MEASURE_NULL)) {
    MEASURE_NULL <- Pct_emVar_bySource[measure == MEASURE, .(Pct = median(Pct)), by = celltype]
    p <- p + geom_hline(data = MEASURE_NULL, aes(yintercept = Pct * 100), linetype = 2)
  } else {
    p <- p + geom_hline(yintercept = MEASURE_NULL * 100, linetype = 2)
  }
  p <- p + theme_plot(rotate.x = 90) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + guides(col = "none") + xlab("")
  p <- p + scale_color_manual(values = color_celline)
  # p <- p + scale_alpha_manual(values=alpha_TRUEFALSE)
  p <- p + ylab(MEASURE_LABEL) + xlab("Introgression source")
  p <- p + facet_grid(rows = vars(celltype))
  p <- p + theme(legend.box = "vertical")
  pdf(sprintf("%s/Pct_emVar_by_source/02c_Pct_%s_%s_by_source_byCelltypeNS.pdf", FIGURE_DIR, MEASURE, EMVAR_LIST), height = 2.3, width = 1.7)
  print(p)
  dev.off()
}
toc()

############################################################################################################
################################## Effect of archaic alleles on expression: ALLELE MATCH  ##################
############################################################################################################


tic("effects of archaic alleles on expression: by Allele match\n")
cat("\nby Allele match (all conditions merged)\n")

color_archaic_MPRA["no match"] <- "lightgrey"

Pct_emVar_byMatch <- get_Pct_emVars(emVARs_obs_ctc[allele_match != "", ], split_by = c("allele_match"), total_by = NULL, emVarlist = EMVAR_LIST)
Pct_emVar_byMatch <- make_it_long(Pct_emVar_byMatch, split_by = "allele_match")
Pct_emVar_byMatch[N_test > minSNP, FDR := p.adjust(Pvalue, "fdr"), by = measure]
Pct_emVar_byMatch[N_test > minSNP][order(Pvalue), head(.SD, 5), by = measure]
fwrite(Pct_emVar_byMatch, sprintf("%s/Pct_emVar_by_match/03a_Pct_emVars_%s_by_match.tsv", FIGURE_DIR, EMVAR_LIST), sep = "\t")

for (i_MEASURE in seq_len(Measure_table[, .N])) {
  MEASURE <- Measure_table[i_MEASURE, measure]
  MEASURE_LABEL <- Measure_table[i_MEASURE, measure_label]
  MEASURE_NULL <- Measure_table[i_MEASURE, null]

  FigData <- Pct_emVar_byMatch[measure == MEASURE & N_test > minSNP, ]
  p <- ggplot(FigData)
  p <- p + geom_hline(yintercept = MEASURE_NULL, col = "lightgrey", linetype = "dashed")
  p <- p + geom_pointrange(aes(x = factor(allele_match, names(color_archaic_MPRA)), y = Pct * 100, ymin = Pct_lo * 100, ymax = Pct_hi * 100, col = allele_match), size = .5)
  # p <- p + facet_grid(rows=vars(measure)) + ylab("Percentage of alleles") + xlab("Archaic source")
  if (is.na(MEASURE_NULL)) {
    MEASURE_NULL <- as.numeric(Pct_emVar_byMatch[measure == MEASURE, .(Pct = median(Pct))])
  }
  p <- p + geom_hline(yintercept = MEASURE_NULL * 100, linetype = 2)
  p <- p + theme_plot(rotate.x = 90) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + guides(col = "none") + xlab("")
  p <- p + scale_color_manual(values = color_archaic_MPRA)
  p <- p + ylab(MEASURE_LABEL) + xlab("Allele match")
  pdf(sprintf("%s/Pct_emVar_by_match/03a_Pct_%s_%s_by_match.pdf", FIGURE_DIR, MEASURE, EMVAR_LIST), height = 2.6, width = 1.3)
  print(p)
  dev.off()
}

cat("\nby Allele match  (each condition separately\n")
Pct_emVar_byMatch <- get_Pct_emVars(emVARs_obs_ctc[allele_match != "", ], split_by = c("allele_match", "COND_ID"), total_by = c("COND_ID"), emVarlist = EMVAR_LIST)
Pct_emVar_byMatch <- make_it_long(Pct_emVar_byMatch, split_by = c("allele_match", "COND_ID"))
Pct_emVar_byMatch[N_test > minSNP, FDR := p.adjust(Pvalue, "fdr"), by = measure]
Pct_emVar_byMatch[N_test > minSNP][order(Pvalue), head(.SD, 5), by = measure]
Pct_emVar_byMatch[, COND_ID := factor(COND_ID, levels = names(color_setup_simplified_norep))]

fwrite(Pct_emVar_byMatch, sprintf("%s/Pct_emVar_by_match/03b_Pct_emVars_%s_by_match_byCondition.tsv", FIGURE_DIR, EMVAR_LIST), sep = "\t")

for (i_MEASURE in seq_len(Measure_table[, .N])) {
  MEASURE <- Measure_table[i_MEASURE, measure]
  MEASURE_LABEL <- Measure_table[i_MEASURE, measure_label]
  MEASURE_NULL <- Measure_table[i_MEASURE, null]

  FigData <- Pct_emVar_byMatch[measure == MEASURE & N_test > minSNP, ]
  p <- ggplot(FigData)
  p <- p + geom_hline(yintercept = MEASURE_NULL, col = "lightgrey", linetype = "dashed")
  p <- p + geom_pointrange(aes(x = factor(allele_match, names(color_archaic_MPRA)), y = Pct * 100, ymin = Pct_lo * 100, ymax = Pct_hi * 100, col = COND_ID, alpha = Pvalue < 0.01), size = .3)
  p <- p + geom_pointrange(data = FigData[Pvalue < 0.01, ], aes(x = factor(allele_match, names(color_archaic_MPRA)), y = Pct * 100, ymin = Pct_lo * 100, ymax = Pct_hi * 100), col = "red", size = .3, alpha = .5)

	# p <- p + facet_grid(rows=vars(measure)) + ylab("Percentage of alleles") + xlab("Archaic source")
  if (is.na(MEASURE_NULL)) {
    MEASURE_NULL <- Pct_emVar_byMatch[measure == MEASURE, .(Pct = median(Pct)), by = COND_ID]
    p <- p + geom_hline(data = MEASURE_NULL, aes(yintercept = Pct * 100), linetype = 2)
  } else {
    p <- p + geom_hline(yintercept = MEASURE_NULL * 100, linetype = 2)
  }
  p <- p + theme_plot(rotate.x = 90) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + guides(col = "none") + xlab("")
  p <- p + scale_color_manual(values = color_setup_simplified_norep)
  p <- p + scale_alpha_manual(values = alpha_TRUEFALSE)
  p <- p + ylab(MEASURE_LABEL) + xlab("Allele match")
  p <- p + facet_wrap(~factor(COND_ID,names(color_setup_simplified_norep)),ncol=4)
  p <- p + theme(legend.box = "vertical")
  pdf(sprintf("%s/Pct_emVar_by_match/03b_Pct_%s_%s_by_match_byCondition.pdf", FIGURE_DIR, MEASURE, EMVAR_LIST), height = 4, width = 2.1)
  print(p)
  dev.off()
}

cat("\nby Allele match  (each celltype separately\n")
Pct_emVar_byMatch <- get_Pct_emVars(emVARs_obs_ctc[condition == "NS" & allele_match != "", ], split_by = c("allele_match", "celltype"), total_by = c("celltype"), emVarlist = EMVAR_LIST)
Pct_emVar_byMatch <- make_it_long(Pct_emVar_byMatch, split_by = c("allele_match", "celltype"))
Pct_emVar_byMatch[N_test > minSNP, FDR := p.adjust(Pvalue, "fdr"), by = measure]
Pct_emVar_byMatch[N_test > minSNP][order(Pvalue), head(.SD, 5), by = measure]
fwrite(Pct_emVar_byMatch, sprintf("%s/Pct_emVar_by_match/03c_Pct_emVars_%s_by_match_byCelltypeNS.tsv", FIGURE_DIR, EMVAR_LIST), sep = "\t")

for (i_MEASURE in seq_len(Measure_table[, .N])) {
  MEASURE <- Measure_table[i_MEASURE, measure]
  MEASURE_LABEL <- Measure_table[i_MEASURE, measure_label]
  MEASURE_NULL <- Measure_table[i_MEASURE, null]

  FigData <- Pct_emVar_byMatch[measure == MEASURE & N_test > minSNP, ]
  p <- ggplot(FigData)
  p <- p + geom_hline(yintercept = MEASURE_NULL, col = "lightgrey", linetype = "dashed")
  p <- p + geom_pointrange(aes(x = factor(allele_match, names(color_archaic_MPRA)), y = Pct * 100, ymin = Pct_lo * 100, ymax = Pct_hi * 100, col = celltype), size = .3, alpha = .5)
  p <- p + geom_pointrange(data = FigData[Pvalue < 0.01, ], aes(x = factor(allele_match, names(color_archaic_MPRA)), y = Pct * 100, ymin = Pct_lo * 100, ymax = Pct_hi * 100), col = "red", size = .3, alpha = .5)
  # p <- p + facet_grid(rows=vars(measure)) + ylab("Percentage of alleles") + xlab("Archaic source")
  if (is.na(MEASURE_NULL)) {
    MEASURE_NULL <- Pct_emVar_byMatch[measure == MEASURE, .(Pct = median(Pct)), by = celltype]
    p <- p + geom_hline(data = MEASURE_NULL, aes(yintercept = Pct * 100), linetype = 2)
  } else {
    p <- p + geom_hline(yintercept = MEASURE_NULL * 100, linetype = 2)
  }
  p <- p + theme_plot(rotate.x = 90) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + guides(col = "none") + xlab("")
  p <- p + scale_color_manual(values = color_celline)
  # p <- p + scale_alpha_manual(values=alpha_TRUEFALSE)
  p <- p + ylab(MEASURE_LABEL) + xlab("Allele match")
  p <- p + facet_grid(rows = vars(celltype))
  p <- p + theme(legend.box = "vertical")
  pdf(sprintf("%s/Pct_emVar_by_match/03c_Pct_%s_%s_by_match_byCelltypeNS.pdf", FIGURE_DIR, MEASURE, EMVAR_LIST), height = 2.3, width = 1.7)
  print(p)
  dev.off()
}
toc()


############################################################################################################
################################## Effect of archaic alleles on expression: INTROGRESSION SCENARIO  ########
############################################################################################################


tic("effects of archaic alleles on expression: by Introgression scenario\n")
cat("\nby Introgression scenario (all conditions merged)\n")

# color_archaic_MPRA["no match"] <- "lightgrey"

Pct_emVar_byScenario <- get_Pct_emVars(emVARs_obs_ctc[Introgression_scenario_top != "", ], split_by = c("Introgression_scenario_top"), total_by = NULL, emVarlist = EMVAR_LIST)
Pct_emVar_byScenario <- make_it_long(Pct_emVar_byScenario, split_by = "Introgression_scenario_top")
Pct_emVar_byScenario[N_test > minSNP, FDR := p.adjust(Pvalue, "fdr"), by = measure]
Pct_emVar_byScenario[N_test > minSNP][order(Pvalue), head(.SD, 5), by = measure]
fwrite(Pct_emVar_byScenario, sprintf("%s/Pct_emVar_by_scenario/04a_Pct_emVars_%s_by_scenario.tsv", FIGURE_DIR, EMVAR_LIST), sep = "\t")

for (i_MEASURE in seq_len(Measure_table[, .N])) {
  MEASURE <- Measure_table[i_MEASURE, measure]
  MEASURE_LABEL <- Measure_table[i_MEASURE, measure_label]
  MEASURE_NULL <- Measure_table[i_MEASURE, null]

  FigData <- Pct_emVar_byScenario[measure == MEASURE & N_test > minSNP, ]
  p <- ggplot(FigData)
  p <- p + geom_hline(yintercept = MEASURE_NULL, col = "lightgrey", linetype = "dashed")
  p <- p + geom_pointrange(aes(x = Introgression_scenario_top, y = Pct * 100, ymin = Pct_lo * 100, ymax = Pct_hi * 100, col = Introgression_scenario_top), size = .5)
  # p <- p + facet_grid(rows=vars(measure)) + ylab("Percentage of alleles") + xlab("Archaic source")
  if (is.na(MEASURE_NULL)) {
    MEASURE_NULL <- as.numeric(Pct_emVar_byScenario[measure == MEASURE, .(Pct = median(Pct))])
  }
  p <- p + geom_hline(yintercept = MEASURE_NULL * 100, linetype = 2)
  p <- p + theme_plot(rotate.x = 90) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + guides(col = "none") + xlab("")
  # p <- p + scale_color_manual(values=color_archaic_MPRA)
  p <- p + ylab(MEASURE_LABEL) + xlab("Introgression scenario")
  pdf(sprintf("%s/Pct_emVar_by_scenario/04a_Pct_%s_%s_by_scenario.pdf", FIGURE_DIR, MEASURE, EMVAR_LIST), height = 2.6, width = 1.3)
  print(p)
  dev.off()
}

cat("\nby Introgression scenario  (each condition separately\n")
Pct_emVar_byScenario <- get_Pct_emVars(emVARs_obs_ctc[Introgression_scenario_top != "", ], split_by = c("Introgression_scenario_top", "COND_ID"), total_by = c("COND_ID"), emVarlist = EMVAR_LIST)
Pct_emVar_byScenario <- make_it_long(Pct_emVar_byScenario, split_by = c("Introgression_scenario_top", "COND_ID"))
Pct_emVar_byScenario[N_test > minSNP, FDR := p.adjust(Pvalue, "fdr"), by = measure]
Pct_emVar_byScenario[N_test > minSNP][order(Pvalue), head(.SD, 5), by = measure]
Pct_emVar_byScenario[, COND_ID := factor(COND_ID, levels = names(color_setup_simplified_norep))]

fwrite(Pct_emVar_byScenario, sprintf("%s/Pct_emVar_by_scenario/04b_Pct_emVars_%s_by_scenario_byCondition.tsv", FIGURE_DIR, EMVAR_LIST), sep = "\t")

for (i_MEASURE in seq_len(Measure_table[, .N])) {
  MEASURE <- Measure_table[i_MEASURE, measure]
  MEASURE_LABEL <- Measure_table[i_MEASURE, measure_label]
  MEASURE_NULL <- Measure_table[i_MEASURE, null]

  FigData <- Pct_emVar_byScenario[measure == MEASURE & N_test > minSNP, ]
  p <- ggplot(FigData)
  p <- p + geom_hline(yintercept = MEASURE_NULL, col = "lightgrey", linetype = "dashed")
  p <- p + geom_pointrange(aes(x = Introgression_scenario_top, y = Pct * 100, ymin = Pct_lo * 100, ymax = Pct_hi * 100, col = COND_ID, alpha = Pvalue < 0.01), size = .3)
  p <- p + geom_pointrange(data = FigData[Pvalue < 0.01, ], aes(x = Introgression_scenario_top, y = Pct * 100, ymin = Pct_lo * 100, ymax = Pct_hi * 100), col = "red", size = .3, alpha = .5)

	# p <- p + facet_grid(rows=vars(measure)) + ylab("Percentage of alleles") + xlab("Archaic source")
  if (is.na(MEASURE_NULL)) {
    MEASURE_NULL <- Pct_emVar_byScenario[measure == MEASURE, .(Pct = median(Pct)), by = COND_ID]
    p <- p + geom_hline(data = MEASURE_NULL, aes(yintercept = Pct * 100), linetype = 2)
  } else {
    p <- p + geom_hline(yintercept = MEASURE_NULL * 100, linetype = 2)
  }
  p <- p + theme_plot(rotate.x = 90) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + guides(col = "none") + xlab("")
  p <- p + scale_color_manual(values = color_setup_simplified_norep)
  p <- p + scale_alpha_manual(values = alpha_TRUEFALSE)
  p <- p + ylab(MEASURE_LABEL) + xlab("Introgression scenario")
  p <- p + facet_wrap(~factor(COND_ID,names(color_setup_simplified_norep)),ncol=4)
  p <- p + theme(legend.box = "vertical")
  pdf(sprintf("%s/Pct_emVar_by_scenario/04b_Pct_%s_%s_by_scenario_byCondition.pdf", FIGURE_DIR, MEASURE, EMVAR_LIST), height = 4, width = 4)
  print(p)
  dev.off()
}

cat("\nby Introgression scenario  (each celltype separately\n")
Pct_emVar_byScenario <- get_Pct_emVars(emVARs_obs_ctc[condition == "NS" & Introgression_scenario_top != "", ], split_by = c("Introgression_scenario_top", "celltype"), total_by = c("celltype"), emVarlist = EMVAR_LIST)
Pct_emVar_byScenario <- make_it_long(Pct_emVar_byScenario, split_by = c("Introgression_scenario_top", "celltype"))
Pct_emVar_byScenario[N_test > minSNP, FDR := p.adjust(Pvalue, "fdr"), by = measure]
Pct_emVar_byScenario[N_test > minSNP][order(Pvalue), head(.SD, 5), by = measure]
fwrite(Pct_emVar_byScenario, sprintf("%s/Pct_emVar_by_scenario/04c_Pct_emVars_%s_by_scenario_byCelltypeNS.tsv", FIGURE_DIR, EMVAR_LIST), sep = "\t")

for (i_MEASURE in seq_len(Measure_table[, .N])) {
  MEASURE <- Measure_table[i_MEASURE, measure]
  MEASURE_LABEL <- Measure_table[i_MEASURE, measure_label]
  MEASURE_NULL <- Measure_table[i_MEASURE, null]

  FigData <- Pct_emVar_byScenario[measure == MEASURE & N_test > minSNP, ]
  p <- ggplot(FigData)
  p <- p + geom_hline(yintercept = MEASURE_NULL, col = "lightgrey", linetype = "dashed")
  p <- p + geom_pointrange(aes(x = Introgression_scenario_top, y = Pct * 100, ymin = Pct_lo * 100, ymax = Pct_hi * 100, col = celltype), size = .3, alpha = .5)
  p <- p + geom_pointrange(data = FigData[Pvalue < 0.01, ], aes(x = Introgression_scenario_top, y = Pct * 100, ymin = Pct_lo * 100, ymax = Pct_hi * 100), col = "red", size = .3, alpha = .5)
  # p <- p + facet_grid(rows=vars(measure)) + ylab("Percentage of alleles") + xlab("Archaic source")
  if (is.na(MEASURE_NULL)) {
    MEASURE_NULL <- Pct_emVar_byScenario[measure == MEASURE, .(Pct = median(Pct)), by = celltype]
    p <- p + geom_hline(data = MEASURE_NULL, aes(yintercept = Pct * 100), linetype = 2)
  } else {
    p <- p + geom_hline(yintercept = MEASURE_NULL * 100, linetype = 2)
  }
  p <- p + theme_plot(rotate.x = 90) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + guides(col = "none") + xlab("")
  p <- p + scale_color_manual(values = color_celline)
  # p <- p + scale_alpha_manual(values=alpha_TRUEFALSE)
  p <- p + ylab(MEASURE_LABEL) + xlab("Introgression scenario")
  p <- p + facet_grid(rows = vars(celltype))
  p <- p + theme(legend.box = "vertical")
  pdf(sprintf("%s/Pct_emVar_by_scenario/04c_Pct_%s_%s_by_scenario_byCelltypeNS.pdf", FIGURE_DIR, MEASURE, EMVAR_LIST), height = 2.3, width = 1.7)
  print(p)
  dev.off()
}
toc()


############################################################################################################
################################## Effect of archaic alleles on expression: POP x ADAPTIVE  ################
############################################################################################################


minSNP <- 5
tic("effects of archaic alleles on expression: by Population x is adaptive\n")
cat("\nby Population  (all conditions merged)\n")

Pct_emVar_byPop <- get_Pct_emVars(emVARs_obs_ctc_full, split_by = c("POP_adaptive", "is_adaptive"), total_by = NULL, emVarlist = EMVAR_LIST)
Pct_emVar_byPop <- make_it_long(Pct_emVar_byPop, split_by = c("POP_adaptive", "is_adaptive"))
Pct_emVar_byPop[N_test > minSNP, FDR := p.adjust(Pvalue, "fdr"), by = measure]
Pct_emVar_byPop[N_test > minSNP][order(Pvalue), head(.SD, 5), by = measure]
Pct_emVar_byPop[, POP_adaptive := toupper(gsub("SGDP_", "", POP_adaptive))]
fwrite(Pct_emVar_byPop, sprintf("%s/Pct_emVar_by_pop_x_adaptive_status/05a_Pct_emVars_%s_by_population_x_isAdaptive.tsv", FIGURE_DIR, EMVAR_LIST), sep = "\t")

for (i_MEASURE in seq_len(Measure_table[, .N])) {
  MEASURE <- Measure_table[i_MEASURE, measure]
  MEASURE_LABEL <- Measure_table[i_MEASURE, measure_label]
  MEASURE_NULL <- Measure_table[i_MEASURE, null]

  FigData <- Pct_emVar_byPop[measure == MEASURE & N_test > minSNP, ]
  p <- ggplot(FigData)
  p <- p + geom_hline(yintercept = MEASURE_NULL, col = "lightgrey", linetype = "dashed")
  p <- p + geom_pointrange(aes(x = is_adaptive, y = Pct * 100, ymin = Pct_lo * 100, ymax = Pct_hi * 100, col = POP_adaptive), size = .3)
  # p <- p + facet_grid(rows=vars(measure)) + ylab("Percentage of alleles") + xlab("Population")
  if (is.na(MEASURE_NULL)) {
    MEASURE_NULL <- as.numeric(Pct_emVar_byPop[measure == MEASURE, .(Pct = median(Pct))])
  }
  p <- p + geom_hline(yintercept = MEASURE_NULL * 100, linetype = 2)
  p <- p + theme_plot(rotate.x = 90) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + guides(col = "none") + xlab("")
  p <- p + scale_color_manual(values = color_populations_MPRA)
  p <- p + ylab(MEASURE_LABEL) + xlab("is adaptive?") + facet_grid(cols = vars(factor(POP_adaptive, names(color_populations_MPRA))))
  pdf(sprintf("%s/Pct_emVar_by_pop_x_adaptive_status/05a_Pct_%s_%s_by_population_x_isAdaptive.pdf", FIGURE_DIR, MEASURE, EMVAR_LIST), height = 2.3, width = 4)
  print(p)
  dev.off()
}

# cat("\nby Population  (each condition separately\n")
# Pct_emVar_byPop <- get_Pct_emVars(emVARs_obs_ctc_full, split_by = c( "POP_adaptive", "COND_ID"), total_by = c( "COND_ID"))
# Pct_emVar_byPop <- make_it_long(Pct_emVar_byPop, split_by = c("POP_adaptive","COND_ID"))
# Pct_emVar_byPop[N_test > minSNP, FDR := p.adjust(Pvalue, "fdr"), by = measure]
# Pct_emVar_byPop[N_test > minSNP][order(Pvalue), head(.SD, 5), by = measure]
# Pct_emVar_byPop[,POP_adaptive:=toupper(gsub('SGDP_','',POP_adaptive))]
# Pct_emVar_byPop[, COND_ID := factor(COND_ID, levels = names(color_setup_simplified_norep))]

# fwrite(Pct_emVar_byPop, sprintf("%s/Pct_emVar_by_pop/01b_Pct_emVars_by_population_byCondition.tsv", FIGURE_DIR), sep = "\t")

# for (i_MEASURE in seq_len(Measure_table[,.N])) {
#   MEASURE <- Measure_table[i_MEASURE, measure]
#   MEASURE_LABEL <- Measure_table[i_MEASURE, measure_label]
#   MEASURE_NULL <- Measure_table[i_MEASURE, null]

#   FigData <- Pct_emVar_byPop[measure==MEASURE & N_test > minSNP,]
#   p <- ggplot(FigData)
#   p <- p + geom_hline(yintercept = MEASURE_NULL, col = "lightgrey", linetype = "dashed")
#   p <- p + geom_pointrange(aes(x = COND_ID, y = Pct * 100, ymin = Pct_lo * 100, ymax = Pct_hi * 100, col = COND_ID, alpha=Pvalue<0.01), size = .3)
#   #p <- p + facet_grid(rows=vars(measure)) + ylab("Percentage of alleles") + xlab("Population")
#   if (is.na(MEASURE_NULL)) {
#     MEASURE_NULL <- Pct_emVar_byPop[measure == MEASURE, .(Pct=median(Pct)), by = POP_adaptive]
#     p <- p + geom_hline(data=MEASURE_NULL, aes(yintercept = Pct*100), linetype=2)
#   }else{
#     p <- p + geom_hline(yintercept=MEASURE_NULL*100, linetype=2)
#   }
#   p <- p + theme_plot(rotate.x = 90) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + guides(col = "none") + xlab("")
#   p <- p + scale_color_manual(values=color_setup_simplified_norep)
#   p <- p + scale_alpha_manual(values=alpha_TRUEFALSE)
#   p <- p + ylab(MEASURE_LABEL) + xlab('Condition')
#   p <- p + facet_grid(rows = vars(factor(POP_adaptive,names(color_populations_MPRA))))
#  p <- p + theme(legend.box = "vertical")
#   pdf(sprintf("%s/Pct_emVar_by_pop/01b_Pct_%s_by_population_byCondition.pdf", FIGURE_DIR, MEASURE), height = 4, width = 2.1)
#   print(p)
#   dev.off()
# }

cat("\nby Population x is adaptive (each celltype separately\n")
Pct_emVar_byPop <- get_Pct_emVars(emVARs_obs_ctc_full[condition == "NS", ], split_by = c("POP_adaptive", "is_adaptive", "celltype"), total_by = c("celltype"), emVarlist = EMVAR_LIST)
Pct_emVar_byPop <- make_it_long(Pct_emVar_byPop, split_by = c("POP_adaptive", "is_adaptive", "celltype"))
Pct_emVar_byPop[N_test > minSNP, FDR := p.adjust(Pvalue, "fdr"), by = measure]
Pct_emVar_byPop[N_test > minSNP][order(Pvalue), head(.SD, 5), by = measure]
Pct_emVar_byPop[, POP_adaptive := toupper(gsub("SGDP_", "", POP_adaptive))]
fwrite(Pct_emVar_byPop, sprintf("%s/Pct_emVar_by_pop_x_adaptive_status/05c_Pct_emVars_%s_by_population_x_isAdaptive_byCelltypeNS.tsv", FIGURE_DIR, EMVAR_LIST), sep = "\t")

for (i_MEASURE in seq_len(Measure_table[, .N])) {
  MEASURE <- Measure_table[i_MEASURE, measure]
  MEASURE_LABEL <- Measure_table[i_MEASURE, measure_label]
  MEASURE_NULL <- Measure_table[i_MEASURE, null]

  FigData <- Pct_emVar_byPop[measure == MEASURE & N_test > minSNP, ]
  p <- ggplot(FigData)
  p <- p + geom_hline(yintercept = MEASURE_NULL, col = "lightgrey", linetype = "dashed")
  p <- p + geom_pointrange(aes(x = is_adaptive, y = Pct * 100, ymin = Pct_lo * 100, ymax = Pct_hi * 100, col = celltype), size = .3, alpha = .5)
  p <- p + geom_pointrange(data = FigData[Pvalue < 0.01, ], aes(x = is_adaptive, y = Pct * 100, ymin = Pct_lo * 100, ymax = Pct_hi * 100), col = "red", size = .3, alpha = .5)
  # p <- p + facet_grid(rows=vars(measure)) + ylab("Percentage of alleles") + xlab("Population")
  if (is.na(MEASURE_NULL)) {
    MEASURE_NULL <- Pct_emVar_byPop[measure == MEASURE, .(Pct = median(Pct)), by = celltype]
    p <- p + geom_hline(data = MEASURE_NULL, aes(yintercept = Pct * 100), linetype = 2)
  } else {
    p <- p + geom_hline(yintercept = MEASURE_NULL * 100, linetype = 2)
  }
  p <- p + theme_plot(rotate.x = 90) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + guides(col = "none") + xlab("")
  p <- p + scale_color_manual(values = color_celline)
  # p <- p + scale_alpha_manual(values=alpha_TRUEFALSE)
  p <- p + ylab(MEASURE_LABEL) + xlab("is adaptive?")
  p <- p + facet_grid(rows = vars(celltype), cols = vars(factor(POP_adaptive, names(color_populations_MPRA))))
  p <- p + theme(legend.box = "vertical")
  pdf(sprintf("%s/Pct_emVar_by_pop_x_adaptive_status/05c_Pct_%s_%s_by_population_x_isAdaptive_byCelltypeNS.pdf", FIGURE_DIR, MEASURE, EMVAR_LIST), height = 2.3, width = 4)
  print(p)
  dev.off()
}
toc()




############################################################################################################
############################################# END  #########################################################
############################################################################################################


# emVARs_obs_ctc[is_emVar_CRITERIA == TRUE, .(.N, mean(ifelse(grepl("enhancer", cre_class_CRITERIA), log2FC_archaic_vs_modern > 0, log2FC_archaic_vs_modern < 0))), keyby = .(grepl("enhancer", cre_class_CRITERIA), ANALYSIS_NAME)][grepl == TRUE, ]





# emVar_byPop <- emVARs_obs_ctc_full[, .(
#   N_expression_increasing_allele = length(unique(ID[is_emVar_CRITERIA == TRUE & log2FC_archaic_vs_modern > 0])),
#   N_activity_increasing_allele = length(unique(ID[is_emVar_CRITERIA == TRUE & ifelse(grepl("enhancer", cre_class_CRITERIA), log2FC_archaic_vs_modern > 0, log2FC_archaic_vs_modern < 0)])),
#   N_activity_altering_allele = length(unique(ID[is_emVar_CRITERIA == TRUE])),
#   N_CRE_overlaping_allele = length(unique(ID[is_tested_cre==TRUE])),
#   N_introgressed_allele = length(unique(ID))
# ), keyby = .(POP_adaptive, is_enhancer = grepl("enhancer", cre_class_CRITERIA))]
# emVar_byPop[, Pct_CREoverlap := N_CRE_overlaping_allele/N_introgressed_allele, by = .(is_enhancer, POP_adaptive)]
# emVar_byPop[, Pct_emVar := N_activity_altering_allele/N_introgressed_allele, by = .(is_enhancer, POP_adaptive)]
# emVar_byPop[, Pct_IncreasingActivity := N_activity_increasing_allele/N_introgressed_allele, by = .(is_enhancer, POP_adaptive)]
# emVar_byPop[, Pct_IncreasingExpression := N_expression_increasing_allele/N_introgressed_allele, by = .(is_enhancer, POP_adaptive)]
# emVar_byPop[, P_overlap := binom.test(N_CRE_overlaping_allele, N_introgressed_allele, p = emVar_byPop[, sum(N_CRE_overlaping_allele) / sum(N_introgressed_allele)])$p.value, by = .(is_enhancer, POP_adaptive)]
# emVar_byPop[, P_emVar := binom.test(N_activity_altering_allele, N_introgressed_allele, p = emVar_byPop[, sum(N_activity_altering_allele) / sum(N_introgressed_allele)])$p.value, by = .(is_enhancer, POP_adaptive)]
# emVar_byPop[, P_sign := binom.test(N_activity_increasing_allele, N_activity_altering_allele)$p.value, by = .(POP_adaptive, is_enhancer)]
# emVar_byPop[, P_sign_expression := binom.test(N_expression_increasing_allele, N_activity_altering_allele)$p.value, by = .(POP_adaptive, is_enhancer)]



# emVar_byPop <- emVARs_obs_ctc_full[, .(
#   N_activity_increasing_allele = length(unique(ID[is_emVar_CRITERIA == TRUE & ifelse(grepl("enhancer", cre_class_CRITERIA), log2FC_archaic_vs_modern > 0, log2FC_archaic_vs_modern < 0)])),
#   N_activity_altering_allele = length(unique(ID[is_emVar_CRITERIA == TRUE])),
#   N_introgressed_allele = length(unique(ID))
# ), keyby = .(POP_adaptive, is_enhancer = grepl("enhancer", cre_class_CRITERIA, COND_ID))]
# emVar_byPop[, P_emVar := binom.test(N_activity_altering_allele, N_introgressed_allele, p = emVar_byPop[, sum(N_activity_altering_allele) / sum(N_introgressed_allele)])$p.value, by = .(is_enhancer, POP_adaptive)]
# emVar_byPop[, P_sign := binom.test(N_activity_increasing_allele, N_activity_altering_allele)$p.value, by = .(POP_adaptive, is_enhancer)]

# emVARs_obs_ctc[is_emVar_CRITERIA == TRUE, .(.N, mean(ifelse(grepl("enhancer", cre_class_CRITERIA), log2FC_archaic_vs_modern > 0, log2FC_archaic_vs_modern < 0))), keyby = .(Introgression_scenario_top)][grepl == TRUE, ]


# # ARCHAIC DERIVED ALLELES ARE ASSOCIATED WITH LOWER ENHANCERS ACTIVITY
# FigData <- emVARs_obs_ctc[is_emVar_CRITERIA == TRUE, ]
# FigData <- FigData[Introgression_scenario_top != "", .(
#   N_emVars = .N,
#   N_higher_activy_archaic = sum(ifelse(grepl("enhancer", cre_class_CRITERIA), log2FC_archaic_vs_modern > 0, log2FC_archaic_vs_modern < 0)),
#   N_higher_activy_modern = sum(ifelse(grepl("enhancer", cre_class_CRITERIA), log2FC_archaic_vs_modern < 0, log2FC_archaic_vs_modern > 0))
# ), keyby = .(crs_type = ifelse(grepl("enhancer", cre_class_CRITERIA), "enhancer", "silencer"), Introgression_scenario_top)]
# FigData[, pct_higher_activy_archaic_inf := binom.test(N_higher_activy_archaic, N_emVars)$conf.int[1], keyby = .(crs_type, Introgression_scenario_top)]
# FigData[, pct_higher_activy_archaic_sup := binom.test(N_higher_activy_archaic, N_emVars)$conf.int[2], keyby = .(crs_type, Introgression_scenario_top)]
# FigData[, binom_P := binom.test(N_higher_activy_archaic, N_emVars)$p.value, keyby = .(crs_type, Introgression_scenario_top)]

# p <- ggplot(FigData) +
#   geom_hline(yintercept = 50, col = "lightgrey", linetype = "dashed")
# p <- p + geom_pointrange(aes(x = Introgression_scenario_top, y = N_higher_activy_archaic / N_emVars * 100, ymin = pct_higher_activy_archaic_inf * 100, ymax = pct_higher_activy_archaic_sup * 100, col = Introgression_scenario_top), size = .5)
# p <- p + facet_grid(cols = vars(crs_type)) + ylab("Percentage of emVars that have\nhigher activity in archaics\n compared to modern humans")
# p <- p + theme_plot(rotate.x = 90) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + guides(col = "none") + xlab("")

# pdf(sprintf("%s/Pct_emVar_by_scenario/Pct_emVars_by_introgression_scenario.pdf", FIGURE_DIR), height = 3.5, width = 3)
# print(p)
# dev.off()


# # emVARs_annot_obs[is_introgressed & posID %chin% posID[FDR<.05 & abs(log2FC_archaic_vs_modern)>.2] & ANALYSIS_SUBTYPE=='celltype_cond',]
# FigData <- emVARs_obs_ctc[is_emVar_CRITERIA == TRUE, ]
# FigData <- merge(FigData, SNP_annot_emVars, by = "posID")
# FigData <- FigData[, .(
#   N_emVars = .N,
#   N_higher_activy_archaic = sum(ifelse(grepl("enhancer", oligo_class_loose), log2FC_archaic_vs_modern > 0, log2FC_archaic_vs_modern < 0)),
#   N_higher_activy_modern = sum(ifelse(grepl("enhancer", oligo_class_loose), log2FC_archaic_vs_modern < 0, log2FC_archaic_vs_modern > 0))
# ), keyby = .(crs_type = ifelse(grepl("enhancer", oligo_class_loose), "enhancer", "silencer"), cut(best_expFC, quantile(best_expFC, seq(0, 1, .2))))]
# FigData[, pct_higher_activy_archaic_inf := binom.test(N_higher_activy_archaic, N_emVars)$conf.int[1], keyby = .(crs_type, cut)]
# FigData[, pct_higher_activy_archaic_sup := binom.test(N_higher_activy_archaic, N_emVars)$conf.int[2], keyby = .(crs_type, cut)]
# FigData[, binom_P := binom.test(N_higher_activy_archaic, N_emVars)$p.value, keyby = .(crs_type, cut)]
# # there is an issue with that: observations are not independant

# FigData <- emVARs_obs_ctc[is_emVar_CRITERIA == TRUE, ]
# FigData <- merge(FigData, SNP_annot_emVars, by = "posID")
# FigData <- FigData[, .(
#   N_emVars = .N,
#   N_higher_activy_archaic = sum(ifelse(grepl("enhancer", oligo_class_loose), log2FC_archaic_vs_modern > 0, log2FC_archaic_vs_modern < 0)),
#   N_higher_activy_modern = sum(ifelse(grepl("enhancer", oligo_class_loose), log2FC_archaic_vs_modern < 0, log2FC_archaic_vs_modern > 0))
# ), keyby = .(crs_type = ifelse(grepl("enhancer", oligo_class_loose), "enhancer", "silencer"), cut(best_context_specificty_BoW / best_expFC, quantile(best_context_specificty_BoW / best_expFC, seq(0, 1, .2))))]
# FigData[, pct_higher_activy_archaic_inf := binom.test(N_higher_activy_archaic, N_emVars)$conf.int[1], keyby = .(crs_type, cut)]
# FigData[, pct_higher_activy_archaic_sup := binom.test(N_higher_activy_archaic, N_emVars)$conf.int[2], keyby = .(crs_type, cut)]
# FigData[, binom_P := binom.test(N_higher_activy_archaic, N_emVars)$p.value, keyby = .(crs_type, cut)]



# # ARCHAIC DERIVED ALLELES ARE ASSOCIATED WITH LOWER ENHANCERS ACTIVITY
# FigData <- emVARs_annot_active_ok[grepl("all", ANALYSIS_NAME) & FDR < .05 & abs(log2FC_archaic_vs_modern) > .2 & !ctrl & is_introgressed, ]
# FigData <- FigData[, .(
#   N_emVars = .N,
#   N_higher_activy_archaic = sum(ifelse(grepl("enhancer", oligo_class_loose), log2FC_archaic_vs_modern > 0, log2FC_archaic_vs_modern < 0)),
#   N_higher_activy_modern = sum(ifelse(grepl("enhancer", oligo_class_loose), log2FC_archaic_vs_modern < 0, log2FC_archaic_vs_modern > 0))
# ), keyby = .(crs_type = ifelse(grepl("enhancer", oligo_class_loose), "enhancer", "silencer"), Adaptive_from)]
# FigData[, pct_higher_activy_archaic_inf := binom.test(N_higher_activy_archaic, N_emVars)$conf.int[1], keyby = .(crs_type, Adaptive_from)]
# FigData[, pct_higher_activy_archaic_sup := binom.test(N_higher_activy_archaic, N_emVars)$conf.int[2], keyby = .(crs_type, Adaptive_from)]
# FigData[, binom_P := binom.test(N_higher_activy_archaic, N_emVars)$p.value, keyby = .(crs_type, Adaptive_from)]

# p <- ggplot(FigData) +
#   geom_hline(yintercept = 50, col = "lightgrey", linetype = "dashed")
# p <- p + geom_pointrange(aes(x = Adaptive_from, y = N_higher_activy_archaic / N_emVars * 100, ymin = pct_higher_activy_archaic_inf * 100, ymax = pct_higher_activy_archaic_sup * 100, col = Adaptive_from), size = .5)
# p <- p + facet_grid(cols = vars(crs_type)) + ylab("Percentage of emVars that have\nhigher activity in archaics\n compared to modern humans")
# p <- p + theme_plot(rotate.x = 90) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + guides(col = "none") + xlab("")

# pdf(sprintf("%s/Agg_emVars_by_archaic_origin_v2.pdf", FIGURE_DIR), height = 3.5, width = 3)
# print(p)
# dev.off()



# # ARCHAIC DERIVED ALLELES ARE ASSOCIATED WITH LOWER ENHANCERS ACTIVITY
# FigData <- emVARs_annot_active_ok[grepl("all", ANALYSIS_NAME) & FDR < .05 & abs(log2FC_archaic_vs_modern) > .2 & !ctrl & is_introgressed, ]
# FigData <- FigData[, .(
#   N_emVars = .N,
#   N_higher_activy_archaic = sum(ifelse(grepl("enhancer", oligo_class_loose), log2FC_archaic_vs_modern > 0, log2FC_archaic_vs_modern < 0)),
#   N_higher_activy_modern = sum(ifelse(grepl("enhancer", oligo_class_loose), log2FC_archaic_vs_modern < 0, log2FC_archaic_vs_modern > 0))
# ), keyby = .(crs_type = ifelse(grepl("enhancer", oligo_class_loose), "enhancer", "silencer"), COND_ID = gsub("(.*_.*)_all", "\\1", ANALYSIS_NAME))]

# FigData[, pct_higher_activy_archaic_inf := binom.test(N_higher_activy_archaic, N_emVars)$conf.int[1], keyby = .(crs_type, COND_ID)]
# FigData[, pct_higher_activy_archaic_sup := binom.test(N_higher_activy_archaic, N_emVars)$conf.int[2], keyby = .(crs_type, COND_ID)]
# FigData[, binom_P := binom.test(N_higher_activy_archaic, N_emVars)$p.value, keyby = .(crs_type, COND_ID)]
# FigData[, COND_ID := gsub("-ACE2", "", COND_ID)]

# p <- ggplot(FigData) +
#   geom_hline(yintercept = 50, col = "lightgrey", linetype = "dashed")
# p <- p + geom_pointrange(aes(x = COND_ID, y = N_higher_activy_archaic / N_emVars * 100, ymin = pct_higher_activy_archaic_inf * 100, ymax = pct_higher_activy_archaic_sup * 100, col = COND_ID), size = .5)
# p <- p + facet_grid(cols = vars(crs_type)) + ylab("Percentage of emVars that have\nhigher activity in archaics\n compared to modern humans")
# p <- p + theme_plot(rotate.x = 90) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + guides(col = "none") + xlab("")
# p <- p + scale_color_manual(values = c(color_setup_simplified_norep, HepG2_all = color_celline["HepG2"], K562_all = color_celline["K562"], A459_all = color_celline["A459"]))
# pdf(sprintf("%s/Agg_emVars_by_COND_ID_v2.pdf", FIGURE_DIR), height = 3.5, width = 3)
# print(p)
# dev.off()

# # Denisovan alleles are associated with stronger activity of immune-related CRS

# FigData <- emVARs_annot_active_ok[grepl("all", ANALYSIS_NAME) & FDR < .05 & abs(log2FC_archaic_vs_modern) > .2 & !ctrl & is_introgressed, ]
# FigData <- FigData[, .(
#   N_emVars = .N,
#   N_higher_activy_archaic = sum(ifelse(grepl("enhancer", oligo_class_loose), log2FC_archaic_vs_modern > 0, log2FC_archaic_vs_modern < 0)),
#   N_higher_activy_modern = sum(ifelse(grepl("enhancer", oligo_class_loose), log2FC_archaic_vs_modern < 0, log2FC_archaic_vs_modern > 0))
# ),
# keyby = .(Adaptive_from, ANALYSIS_NAME, crs_type = ifelse(grepl("enhancer", oligo_class_loose), "enhancer", "silencer"))
# ]
# FigData[, pct_higher_activy_archaic_inf := binom.test(N_higher_activy_archaic, N_emVars)$conf.int[1], keyby = .(crs_type, Adaptive_from, ANALYSIS_NAME)]
# FigData[, pct_higher_activy_archaic_sup := binom.test(N_higher_activy_archaic, N_emVars)$conf.int[2], keyby = .(crs_type, Adaptive_from, ANALYSIS_NAME)]
# FigData[, binom_P := binom.test(N_higher_activy_archaic, N_emVars)$p.value, keyby = .(crs_type, Adaptive_from, ANALYSIS_NAME)]
# FigData <- merge(FigData, condition_summary, by.x = "ANALYSIS_NAME", by.y = "analysis_name")
# FigData[, COND_ID := factor(COND_ID, levels = names(color_setup_simplified_norep))]
# p <- ggplot(FigData) +
#   geom_hline(yintercept = 50, col = "lightgrey")
# p <- p + geom_pointrange(aes(x = COND_ID, y = N_higher_activy_archaic / N_emVars * 100, ymin = pct_higher_activy_archaic_inf * 100, ymax = pct_higher_activy_archaic_sup * 100, col = COND_ID))
# p <- p + facet_grid(rows = vars(crs_type), cols = vars(Adaptive_from))
# p <- p + ylab("Percentage of emVars that have\nhigher activity in archaics compared to AMH")
# p <- p + theme_plot(rotate.x = 90) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + scale_color_manual(values = color_setup_simplified_norep) + guides(col = "none")
# p <- p + scale_color_manual(values = color_setup_simplified_norep)

# pdf(sprintf("%s/Agg_emVars_by_COND_origin.pdf", FIGURE_DIR), height = 7, width = 7)
# print(p)
# dev.off()


# FigData <- emVARs_annot_active_ok[grepl("all", ANALYSIS_NAME) & FDR < .05 & abs(log2FC_archaic_vs_modern) > .2 & !ctrl & is_introgressed, ]
# FigData <- FigData[, .(
#   N_emVars = .N,
#   N_higher_activy_archaic = sum(ifelse(grepl("enhancer", oligo_class_loose), log2FC_archaic_vs_modern > 0, log2FC_archaic_vs_modern < 0)),
#   N_higher_activy_modern = sum(ifelse(grepl("enhancer", oligo_class_loose), log2FC_archaic_vs_modern < 0, log2FC_archaic_vs_modern > 0))
# ),
# keyby = .(Adaptive_from, ANALYSIS_NAME)
# ]
# FigData[, pct_higher_activy_archaic_inf := binom.test(N_higher_activy_archaic, N_emVars)$conf.int[1], keyby = .(Adaptive_from, ANALYSIS_NAME)]
# FigData[, pct_higher_activy_archaic_sup := binom.test(N_higher_activy_archaic, N_emVars)$conf.int[2], keyby = .(Adaptive_from, ANALYSIS_NAME)]
# FigData[, binom_P := binom.test(N_higher_activy_archaic, N_emVars)$p.value, keyby = .(Adaptive_from, ANALYSIS_NAME)]
# FigData <- merge(FigData, condition_summary, by.x = "ANALYSIS_NAME", by.y = "analysis_name")
# FigData[, COND_ID := factor(COND_ID, levels = names(color_setup_simplified_norep))]
# p <- ggplot(FigData) +
#   geom_hline(yintercept = 50, col = "lightgrey")
# p <- p + geom_pointrange(aes(x = COND_ID, y = N_higher_activy_archaic / N_emVars * 100, ymin = pct_higher_activy_archaic_inf * 100, ymax = pct_higher_activy_archaic_sup * 100, col = COND_ID))
# p <- p + facet_grid(cols = vars(Adaptive_from))
# p <- p + theme_plot(rotate.x = 90) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + scale_color_manual(values = color_setup_simplified_norep) + guides(col = "none")
# p <- p + scale_color_manual(values = color_setup_simplified_norep)
# p <- p + ylab("Percentage of emVars that have\nhigher activity in archaics compared to AMH")

# pdf(sprintf("%s/Agg_emVars_by_COND_origin_v2.pdf", FIGURE_DIR), height = 7, width = 7)
# print(p)
# dev.off()

# # ARCHAIC DERIVED ALLELES ARE ASSOCIATED WITH LOWER ENHANCERS ACTIVITY

# emVARs_annot_active_ok[grepl("all", ANALYSIS_NAME) & FDR < .05 & abs(log2FC_archaic_vs_modern) > .2 & !ctrl & is_introgressed, .(.N, mean(ifelse(grepl("enhancer", oligo_class_loose), log2FC_archaic_vs_modern > 0, log2FC_archaic_vs_modern < 0))), keyby = .(Adaptive_from, ANALYSIS_NAME, grepl("enhancer", oligo_class_loose))]

# OAS_SNP1 <- "12:113364366"
# OAS_SNP2 <- "12:113364375"
# #############################################################################
# #############################################################################
# #############################################################################

# # A note on the analysis startegy for power related analyses :
# # For each CRS, we are going to :
# # 1. report the number of barcode in the library and
# # 2. estimate the power to detect activity at the reported level of activity
# #    2a. this can be computed for the Winner's curse adjusted beta values
# #        (use FIQT: F DR I nverse Q uantile T ransformation).
# #         doi.org/10.1093/bioinformatics/btw303
# #    2b. perhaps use local FDR from fdrtool rather than FDR ?

# # 3. test whether probability of having an emVar varies across tissue/conditionb, adjiusting for difference in power.
