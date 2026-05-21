# 06c_upset_plot_emVars_Diff.R

# running: sbatch -p geh,common --mem=30G  00_Rscript.sh MPRA_count_exp6_analysisZ/06c_upset_plot_emVars_Diff.R

MPRA_DIR <- "/pasteur/helix/projects/evo_immuno_pop/MPRA"
ANALYSIS_DIR <- "MPRA_count_exp6_analysis3"

source(sprintf("%s/scripts/%s/00_starter.R", MPRA_DIR, ANALYSIS_DIR))
source(sprintf("%s/scripts/%s/00_load_associations.R", MPRA_DIR, ANALYSIS_DIR))


RUN_ID <- "RUN3_Z2_nBC10"
CRITERIA_ACTIVE <- "FDR5_scrambled5pct_FC0.2_GCnorm"
CRITERIA_EMVARS <- "EmVar_FDR5_FC.2"
CRITERIA_DIFF <- "Diff_FDR5_FC.2"
CRITERIA_EMVARS_DIFF <- "EmVarDiff_FDR5_FC.2"
FILTER_ACTIVE <- TRUE
FILTER_DIFF <- FALSE
FILTER_EMVAR <- TRUE
STIM_ONLY <- FALSE

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
  if (cmd[i] == "--criteria_emvar" || cmd[i] == "-e") {
    CRITERIA_EMVARS <- cmd[i + 1]
  }
  if (cmd[i] == "--criteria_emvar_diff" || cmd[i] == "-ed") {
    CRITERIA_EMVARS_DIFF <- cmd[i + 1]
  }
  if (cmd[i] == "--filter_active" || cmd[i] == "-fa") {
    FILTER_ACTIVE <- as.logical(cmd[i + 1])
  }
  if (cmd[i] == "--filter_diff" || cmd[i] == "-fd") {
    FILTER_DIFF <- as.logical(cmd[i + 1])
  }
	if(cmd[i] == "--filter_emvar" || cmd[i] == "-fe") {
    FILTER_EMVAR <- as.logical(cmd[i + 1])
	}
  if(cmd[i] == "--stim_only" || cmd[i] == "-s") {
    STIM_ONLY <- as.logical(cmd[i + 1])
  }
}

if (!FILTER_ACTIVE) {
  CRITERIA_ACTIVE_OUT <- paste0("noActivityFilter_", CRITERIA_ACTIVE)
} else {
  CRITERIA_ACTIVE_OUT <- CRITERIA_ACTIVE
}

if (!FILTER_DIFF) {
  CRITERIA_DIFF_OUT <- paste0("noDiffFilter_", CRITERIA_DIFF)
} else {
  CRITERIA_DIFF_OUT <- CRITERIA_DIFF
}

if (!FILTER_EMVAR) {
  CRITERIA_EMVARS_OUT <- paste0("noEmVarFilter_", CRITERIA_EMVARS)
} else {
  CRITERIA_EMVARS_OUT <- CRITERIA_EMVARS
}
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
SNP_annot_v5 <- merge(SNP_annot_v4[,-"Introgression_scenario"],selected_annot_wide,by=c('ID','posID'))
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

FIGURE_DIR <- sprintf("%s/figures/%s/%s/05_emVarDiff/%s/%s/%s/%s/%s", MPRA_DIR, ANALYSIS_DIR, RUN_ID, CRITERIA_EMVARS_DIFF, CRITERIA_EMVARS_OUT, CRITERIA_DIFF_OUT, CRITERIA_ACTIVE_OUT, STIM_ONLY_OUT)

dir.create(FIGURE_DIR, recursive = TRUE)
dir.create(sprintf("%s/SupTables/", FIGURE_DIR), recursive = TRUE)

# load tested oligos
source(sprintf("%s/scripts/%s/02z__define_includedCRS.R", MPRA_DIR, ANALYSIS_DIR))
dir.create(sprintf("%s/upset_plots/", FIGURE_DIR))

####################################################################
################# compute FDR, emVars diff ######################
####################################################################

# emVARs_Diff_obs_response = fread(file = sprintf("%s/all_emVars_diff_obs_response.tsv", EMVAR_DIFF_DIR), sep = "\t")


# emVARs_Diff_any <- emVARs_Diff_obs_response[is_emVar_Diff_CRITERIA == TRUE, .(posID, celline, condition, Introgression_source = Introgression_source_top, Introgression_scenario = Introgression_scenario_top, allele_match)]
# emVars_Diff_any_increased <- emVARs_Diff_obs_response[is_emVar_Diff_CRITERIA == TRUE & log2FC_archaic_vs_modern > 0, .(posID, celline, condition, Introgression_source = Introgression_source_top, Introgression_scenario = Introgression_scenario_top, allele_match)]
# emVars_any_decreased <- emVARs_annot_ctc[is_emVar_CRITERIA == TRUE & log2FC_archaic_vs_modern < 0, .(posID, celline, condition, Introgression_source = Introgression_source_top, Introgression_scenario = Introgression_scenario_top, allele_match)]

# emVars_NS <- emVARs_annot_ctc[is_emVar_CRITERIA == TRUE & condition == "NS", .(posID, celline, condition, Introgression_source = Introgression_source_top, Introgression_scenario = Introgression_scenario_top, allele_match)]
# emVars_NS_increased <- emVARs_annot_ctc[is_emVar_CRITERIA == TRUE & condition == "NS" & log2FC_archaic_vs_modern > 0, .(posID, celline, condition, Introgression_source = Introgression_source_top, Introgression_scenario = Introgression_scenario_top, allele_match)]
# emVars_NS_decreased <- emVARs_annot_ctc[is_emVar_CRITERIA == TRUE & condition == "NS" & log2FC_archaic_vs_modern < 0, .(posID, celline, condition, Introgression_source = Introgression_source_top, Introgression_scenario = Introgression_scenario_top, allele_match)]


# oligo_response <- merge(oligo_response, SNP_annot_v5[, .(crsID, POP_adaptive, Introgressed_from = Adaptive_from, Introgression_scenario)], by = "crsID")
# oligo_response_any <- oligo_response[FDR < .05 & abs(log2FC_2vs1) > logFC_TH & excluded == FALSE & ctrl == FALSE, .(crsID, celline, condition, type, POP_adaptive, Introgressed_from)]
# oligo_response_increased <- oligo_response[FDR < .05 & log2FC_2vs1 > logFC_TH & excluded == FALSE & ctrl == FALSE, .(crsID, celline, condition, type, POP_adaptive, Introgressed_from)]
# oligo_response_decreased <- oligo_response[FDR < .05 & log2FC_2vs1 < -logFC_TH & excluded == FALSE & ctrl == FALSE, .(crsID, celline, condition, type, POP_adaptive, Introgressed_from)]

# N_Context <- oligo_response_any[, .(nContext = .N), by = crsID][, .N, keyby = nContext][, Pct := N / sum(N)][1:.N]

response_emVars_obs <- fread(sprintf("%s/all_emVars_diff_obs_response.tsv", EMVARDIFF_DIR))
response_emVars_any <- response_emVars_obs[is_emVar_Diff_CRITERIA==TRUE, .(crsID, celline, condition, Introgression_source = Introgression_source_top, Introgression_scenario = Introgression_scenario_top, allele_match)]
response_emVars_stim <- response_emVars_obs[is_emvar_cre_group2 & !is_emvar_cre_group1 & is_emVar_Diff_CRITERIA, .(crsID, celline, condition, Introgression_source = Introgression_source_top, Introgression_scenario = Introgression_scenario_top, allele_match)]
response_emVars_ns <- response_emVars_obs[is_emvar_cre_group1 & !is_emvar_cre_group2 & is_emVar_Diff_CRITERIA, .(crsID, celline, condition, Introgression_source = Introgression_source_top, Introgression_scenario = Introgression_scenario_top, allele_match)]
#response_emVars_any <- merge(response_emVars_any, condition_summary, by = "COND_ID")

response_emVars_any_suggestive <- response_emVars_obs[crsID %in% response_emVars_any[,crsID] & pval_emp<0.05, .(crsID, celline, condition, Introgression_source = Introgression_source_top, Introgression_scenario = Introgression_scenario_top, allele_match)]
# emVars_any <- emVARs_annot_ctc[is_emVar_CRITERIA == TRUE, .(posID, celline, condition, Introgression_source = Introgression_source_top, Introgression_scenario = Introgression_scenario_top, allele_match)]
# emVars_any_increased <- emVARs_annot_ctc[is_emVar_CRITERIA == TRUE & log2FC_archaic_vs_modern > 0, .(posID, celline, condition, Introgression_source = Introgression_source_top, Introgression_scenario = Introgression_scenario_top, allele_match)]
# emVars_any_decreased <- emVARs_annot_ctc[is_emVar_CRITERIA == TRUE & log2FC_archaic_vs_modern < 0, .(posID, celline, condition, Introgression_source = Introgression_source_top, Introgression_scenario = Introgression_scenario_top, allele_match)]

# emVars_NS <- emVARs_annot_ctc[is_emVar_CRITERIA == TRUE & condition == "NS", .(posID, celline, condition, Introgression_source = Introgression_source_top, Introgression_scenario = Introgression_scenario_top, allele_match)]
# emVars_NS_increased <- emVARs_annot_ctc[is_emVar_CRITERIA == TRUE & condition == "NS" & log2FC_archaic_vs_modern > 0, .(posID, celline, condition, Introgression_source = Introgression_source_top, Introgression_scenario = Introgression_scenario_top, allele_match)]
# emVars_NS_decreased <- emVARs_annot_ctc[is_emVar_CRITERIA == TRUE & condition == "NS" & log2FC_archaic_vs_modern < 0, .(posID, celline, condition, Introgression_source = Introgression_source_top, Introgression_scenario = Introgression_scenario_top, allele_match)]


####################################################################
######################## uspet plots, oligo Diff ###################
####################################################################

# library(MPRAnalyze)
library(ComplexUpset)
library(ggplot2)
library(dplyr)
library(data.table)

names_stimulation <- names(color_setup_simplified_norep)[!grepl("_NS$", names(color_setup_simplified_norep))]
#names_stimulation <- setdiff(names_stimulation, "HepG2_DEX")

create_upset <- function(fig_data, ...) {
  Archaic_origin <- case_when(
    fig_data$Introgression_source == "Vindija" ~ "Vindija",
    fig_data$Introgression_source == "Denisova" ~ "Denisova",
    fig_data$Introgression_source == "Vindija/Denisova" ~ "Vindija/Denisova",
    fig_data$Introgression_source == "undetermined" ~ "Undetermined",
		fig_data$Introgression_source == "Undetermined" ~ "Undetermined",
    fig_data$Introgression_source == "" ~ "Undetermined",
    TRUE ~ "WTF"
  )
  archaic_color <- c("AMH" = grey(0.8), "Denisova" = "#80CDC1", "Vindija" = "#FDAE61", "Vindija/Denisova" = "#FEE08B", "Undetermined" = grey(0.5))

  Archaic_origin <- factor(Archaic_origin, c("Denisova", "Vindija/Denisova", "Vindija", "Undetermined"))
  # my_base_annotations <- list(
  #   "Intersection size" = intersection_size(
  #     mapping = aes(fill = Archaic_origin),
  #     text = list(
  #       vjust = .5,
  #       hjust = -0.1,
  #       angle = 90,
  #       size = 2.5
  #     ),
  #     bar_number_threshold = 1, # show all numbers on top of bars
  #     width = 0.8
  #   ) # reduce width of the bars
  #   + scale_y_continuous(expand = expansion(mult = c(0, 0.1)))
  #     + scale_fill_manual(values = archaic_color)
  #     + guides(fill = guide_legend(title = "Archaic allele observed in:"))
  #     + theme(
  #       # hide grid lines
  #       panel.grid.major = element_blank(),
  #       panel.grid.minor = element_blank(),
  #     )
  # )

  my_base_annotations <- list(
    "Intersection size" = intersection_size(
      mapping = aes(fill = Archaic_origin),
      text = list(
        vjust = .5,
        hjust = -0.1,
        angle = 90,
        size = 2.5
      ),
      bar_number_threshold = 1, # show all numbers on top of bars
      width = 0.8
    ) # reduce width of the bars
    + scale_y_continuous(expand = expansion(mult = c(0, 0.2)))
      + scale_fill_manual(values = archaic_color)
      + guides(fill = guide_legend(title = "Introgressed from:"))
      + theme(
        # hide grid lines
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
      )
  )

  my_queries <- lapply(seq_along(color_setup_simplified_norep[names_stimulation]), function(i) {
    upset_query(set = names_stimulation[i], fill = color_setup_simplified_norep[names_stimulation[i]])
  })

  valid_sets <- names_stimulation[colSums(fig_data[, mget(names_stimulation)]) > 0]
  valid_colors <- color_setup_simplified_norep[valid_sets]
  my_queries <- lapply(seq_along(valid_sets), function(i) {
    upset_query(set = valid_sets[i], fill = valid_colors[[valid_sets[i]]])
  })
  #
  p <- upset(as.data.frame(fig_data),
    colnames(fig_data),
    name = "celltype & condition",
    height_ratio = 1.25,
    width_ratio = 0.25,
    n_intersections = 50,
    queries = my_queries,
    base_annotations = my_base_annotations,
    sort_sets = FALSE,
    set_sizes = (upset_set_size()
    + geom_text(aes(label = ..count..), hjust = 1.1, stat = "count", size = 2.5)
      + scale_y_continuous(expand = expansion(mult = c(0.2, 0)), trans = "reverse")
      + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())),
    matrix = intersection_matrix(geom = geom_point(
      shape = "circle filled",
      size = 3.5,
      stroke = 0.45
    )) + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()),
    ...
  )
  # p <- p + theme_plot(rotate.x=90)
  return(p)
}


data_names <- c("response_emVars_any", "response_emVars_stim", "response_emVars_ns", "response_emVars_any_suggestive")

data_to_plot <- lapply(seq_along(data_names), function(i) {
  data_name <- data_names[i]
  data <- get(data_name)
  #  data <- data[!duplicated(crsID),]
  data <- dcast(data, crsID + Introgression_scenario + Introgression_source ~ factor(paste(celline, condition, sep = "_"), names_stimulation), fill = 0, fun.aggregate = length, drop = c(TRUE, FALSE))
  data <- cbind(data[, c("crsID", "Introgression_scenario", "Introgression_source")], data[, mget(rev(names_stimulation))] > 0)
  rownames(data) <- data$crsID
  data
})

for (i in seq_along(data_to_plot)) {
  fig_data <- data_to_plot[[i]]
  cat(data_names[i], "\n")
  p <- create_upset(unique(fig_data))
  
  pdf(sprintf("%s/upset_plots/01a_upsetPlot_responseEmVars_%s.pdf", FIGURE_DIR, data_names[i]), height = 3, width = ifelse(i==4,8,5))
  print(p)
  dev.off()
  fwrite(unique(fig_data), sprintf("%s/SupTables/SupTable_01a_upsetPlot_responseEmVars_%s_sourcedata.tsv.gz", FIGURE_DIR, data_names[i]),sep='\t')
  #fig_data <- fread(sprintf("%s/SupTables/SupTable_01a_upsetPlot_responseEmVars_%s_sourcedata.tsv.gz", FIGURE_DIR, "response_emVars_any_suggestive"),sep='\t')
}

# data_to_plot <- lapply(seq_along(data_names), function(i) {
#   data_name <- data_names[i]
#   data <- get(data_name)
#   data <- dcast(data, crsID + POP_adaptive + Adaptive_from ~ COND_ID, fill = 0, fun.aggregate = length)
#   data <- cbind(data[, c("crsID", "POP_adaptive", "Adaptive_from")], data[, mget(rev(names_stimulation))] > 0)
#   rownames(data) <- data$crsID
#   data
# })


# for (i in seq_along(data_to_plot)) {
#   fig_data <- data_to_plot[[i]]
#   cat(data_names[i], "\n")
#   p <- create_upset(fig_data)

#   pdf(sprintf("%s/13a_upsetPlot_cellines_%s.pdf", FIGURE_DIR, data_names[i]), height = 4, width = 5)
#   print(p)
#   dev.off()

#   p <- create_upset(fig_data, sort_intersections_by = "degree")

#   pdf(sprintf("%s/13b_upsetPlot_cellines_%s_notsorted.pdf", FIGURE_DIR, data_names[i]), height = 4, width = 5.5)
#   print(p)
#   dev.off()
# }

#################################################################
####### piechart of context-specific responses ##################
#################################################################

N_Context <- response_emVars_any[, .(nContext = .N), by = crsID][, .N, keyby = nContext][, Pct := N / sum(N)][seq_len(.N)]
N_Context

# make a pie chart
N_Context[, label := paste0("  ", round(Pct * 100, ifelse(Pct > 0.02, 0, 1)), "%")]
#N_Context[, label := c(paste0("  ", round(Pct * 100), "%")[Pct > .03], "", "2.4%", "")]

N_Context[, x_pos := 1.65]
p <- ggplot(N_Context, aes(x = "", y = Pct, fill = factor(nContext))) +
  geom_bar(stat = "identity", width = 1, col = "black")
p <- p + coord_polar("y", start = 0) + scale_fill_brewer(palette = "Blues") + theme_void() + labs(fill = "Number of contexts\nwhere CRS is responsive")
p <- p + geom_text(aes(x = x_pos, label = label), position = position_stack(vjust = 0.5)) + theme(legend.position = "bottom")
pdf(sprintf("%s/upset_plots/01b_piechart_responseEmVars_context_number.pdf", FIGURE_DIR), height = 5, width = 4)
print(p)
dev.off()

N_Context_sugg <- response_emVars_any_suggestive[, .(nContext = .N), by = crsID][, .N, keyby = nContext][, Pct := N / sum(N)][seq_len(.N)]
N_Context_sugg

# make a pie chart
N_Context_sugg[, label := paste0("  ", round(Pct * 100, ifelse(Pct > 0.02, 0, 1)), "%")]
#N_Context[, label := c(paste0("  ", round(Pct * 100), "%")[Pct > .03], "", "2.4%", "")]

N_Context_sugg[, x_pos := 1.65]
p <- ggplot(N_Context_sugg, aes(x = "", y = Pct, fill = factor(nContext))) +
  geom_bar(stat = "identity", width = 1, col = "black")
p <- p + coord_polar("y", start = 0) + scale_fill_brewer(palette = "Blues") + theme_void() + labs(fill = "Number of contexts\nwhere CRS is responsive")
p <- p + geom_text(aes(x = x_pos, label = label), position = position_stack(vjust = 0.5)) + theme(legend.position = "bottom")
pdf(sprintf("%s/upset_plots/01c_piechart_responseEmVars_context_number_Pemp05.pdf", FIGURE_DIR), height = 5, width = 4)
print(p)
dev.off()


cat('\nAll done!\n')
q("no")