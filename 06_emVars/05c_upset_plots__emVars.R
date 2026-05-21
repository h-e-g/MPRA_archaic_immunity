# running: sbatch -p geh,common --mem=30G  00_Rscript.sh MPRA_count_exp6_analysisZ/03b_aggMPRAnalyse_oligo_activity.R

MPRA_DIR <- "/pasteur/helix/projects/evo_immuno_pop/MPRA"
ANALYSIS_DIR <- "MPRA_count_exp6_analysis3"

source(sprintf("%s/scripts/%s/00_starter.R", MPRA_DIR, ANALYSIS_DIR))
source(sprintf("%s/scripts/%s/00_load_associations.R", MPRA_DIR, ANALYSIS_DIR))
source(sprintf("%s/scripts/misc_plots.R", MPRA_DIR))

RUN_ID <- "RUN2_Z2_nBC10"
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
# load annotation of oligos (before: CRS)
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
dir.create(sprintf("%s/upset_plots/", FIGURE_DIR))

# selected_annot

# IN_DIR <- sprintf("%s/data/%s/03a_aggMPRA_analyse/%s/", MPRA_DIR, ANALYSIS_DIR, RUN_ID)
# OUT_DIR <- sprintf("%s/data/%s/03b_computeFDR/%s/", MPRA_DIR, ANALYSIS_DIR, RUN_ID)
# dir.create(OUT_DIR, recursive = TRUE)

# emVARs_all_results <- fread(sprintf("%s/all_emVars_results.tsv.gz", IN_DIR), sep = "\t")
# emVARs_all_results[, crsID := posID]

# CRS_activity_file <- sprintf("%s/CRS_activity__all.tsv.gz", OUT_DIR)
# CRS_activity <- fread(CRS_activity_file)

# emVARs_all_active <- fread(sprintf("%s/all_emVars_activeCRS_annotated.tsv.gz", OUT_DIR), sep = "\t")

# emVARs_celltype_active <- fread(sprintf("%s/all_emVars_activeCRS_annotated_celltype.tsv", OUT_DIR))
# emVARs_celltype_active <- merge(emVARs_celltype_active, condition_summary, by.x = "ANALYSIS_NAME", by.y = "analysis_name")

# emVARs_reps_active <- fread(sprintf("%s/all_emVars_by_sample_activeCRS_annotated.tsv", OUT_DIR))
# emVARs_reps_active <- merge(emVARs_reps_active, condition_summary_reps, by.x = "ANALYSIS_NAME", by.y = "analysis_name")

emVARs_annot_ctc <- fread(sprintf("%s/all_emVars_annotated_celltype__%s.tsv", EMVAR_DIR, CRITERIA_ACTIVE_OUT))

emVars_any <- emVARs_annot_ctc[is_emVar_CRITERIA == TRUE, .(posID, celline, condition, Introgression_source = Introgression_source_top, Introgression_scenario = Introgression_scenario_top, allele_match)]
emVars_any_increased <- emVARs_annot_ctc[is_emVar_CRITERIA == TRUE & log2FC_archaic_vs_modern > 0, .(posID, celline, condition, Introgression_source = Introgression_source_top, Introgression_scenario = Introgression_scenario_top, allele_match)]
emVars_any_decreased <- emVARs_annot_ctc[is_emVar_CRITERIA == TRUE & log2FC_archaic_vs_modern < 0, .(posID, celline, condition, Introgression_source = Introgression_source_top, Introgression_scenario = Introgression_scenario_top, allele_match)]

emVars_NS <- emVARs_annot_ctc[is_emVar_CRITERIA == TRUE & condition == "NS", .(posID, celline, condition, Introgression_source = Introgression_source_top, Introgression_scenario = Introgression_scenario_top, allele_match)]
emVars_NS_increased <- emVARs_annot_ctc[is_emVar_CRITERIA == TRUE & condition == "NS" & log2FC_archaic_vs_modern > 0, .(posID, celline, condition, Introgression_source = Introgression_source_top, Introgression_scenario = Introgression_scenario_top, allele_match)]
emVars_NS_decreased <- emVARs_annot_ctc[is_emVar_CRITERIA == TRUE & condition == "NS" & log2FC_archaic_vs_modern < 0, .(posID, celline, condition, Introgression_source = Introgression_source_top, Introgression_scenario = Introgression_scenario_top, allele_match)]

emVars_any_suggestive <- emVARs_annot_ctc[posID %in% emVars_any[,posID]  & pval_emp<0.05 , .(posID, celline, condition, Introgression_source = Introgression_source_top, Introgression_scenario = Introgression_scenario_top, allele_match)]

####################################################################
######################## uspet plots, emVars #######################
####################################################################

# library(MPRAnalyze)
library(ComplexUpset)
library(ggplot2)
library(dplyr)
library(data.table)

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
  # Allele_match <- case_when(
  #   fig_data$allele_match == "Vindija" ~ "Vindija",
  #   fig_data$allele_match == "Denisova" ~ "Denisova",
  #   fig_data$allele_match == "Vindija/Denisova" ~ "Vindija/Denisova",
  #   fig_data$allele_match == "no match" ~ "no match",
  #   TRUE ~ "WTF"
  # )
  archaic_color <- c("AMH" = grey(0.8), "Denisova" = "#80CDC1", "Vindija" = "#FDAE61", "Vindija/Denisova" = "#FEE08B", "Undetermined" = grey(0.5))

  Archaic_origin <- factor(Archaic_origin, c("Denisova", "Vindija/Denisova", "Vindija", "Undetermined"))
  # Allele_match <- factor(Allele_match, c("Vindija", "Denisova", "Vindija/Denisova", "no match"))

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
  my_queries <- lapply(seq_along(color_celline), function(i) {
    upset_query(set = names(color_celline)[i], fill = color_celline[i])
  })

  p <- upset(as.data.frame(fig_data),
    colnames(fig_data),
    name = "celltype",
    height_ratio = 0.25,
    width_ratio = 0.5,
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


data_names <- c(
  "emVars_any", "emVars_any_increased", "emVars_any_decreased", "emVars_NS", "emVars_NS_increased", "emVars_NS_decreased", "emVars_any_suggestive"
)

data_to_plot <- lapply(seq_along(data_names), function(i) {
  data_name <- data_names[i]
  data <- get(data_name)[celline != "Calu3", ]
  data <- dcast(data, posID + Introgression_scenario + Introgression_source ~ celline, fill = 0, fun.aggregate = length)
  data <- cbind(data[, c("posID", "Introgression_scenario", "Introgression_source")], data[, mget(rev(names(color_celline[names(color_celline) != "Calu3"])))] > 0)
  rownames(data) <- data$posID
  data
})


for (i in seq_along(data_to_plot)) {
  fig_data <- data_to_plot[[i]]
  cat(data_names[i], "\n")
  p <- create_upset(fig_data)

  pdf(sprintf("%s/upset_plots/upsetPlot_%s.pdf", FIGURE_DIR, data_names[i]), height = 4, width = 5)
  print(p)
  dev.off()
  fwrite(unique(fig_data), sprintf("%s/SupTables/SupTable_01a_upsetPlot_responseEmVars_%s_sourcedata.tsv.gz", FIGURE_DIR, data_names[i]),sep='\t')
  
  p <- create_upset(fig_data, sort_intersections_by = "degree")

  pdf(sprintf("%s/upset_plots/upsetPlot_%s_notsorted.pdf", FIGURE_DIR, data_names[i]), height = 4, width = 5.5)
  print(p)
  dev.off()
}


#################################################################
####### piechart of context-specific responses ##################
#################################################################

N_Context <- emVars_any[, .(nContext = .N), by = posID][, .N, keyby = nContext][, Pct := N / sum(N)][seq_len(.N)]
N_Context

# make a pie chart
N_Context[, label := paste0("  ", round(Pct * 100, ifelse(Pct > 0.02, 0, 1)), "%")]
#N_Context[, label := c(paste0("  ", round(Pct * 100), "%")[Pct > .03], "", "2.4%", "")]

N_Context[, x_pos := 1.65]
p <- ggplot(N_Context, aes(x = "", y = Pct, fill = factor(nContext))) +
  geom_bar(stat = "identity", width = 1, col = "black")
p <- p + coord_polar("y", start = 0) + scale_fill_brewer(palette = "Blues") + theme_void() + labs(fill = "Number of contexts\nwhere CRS is responsive")
p <- p + geom_text(aes(x = x_pos, label = label), position = position_stack(vjust = 0.5)) + theme(legend.position = "bottom")
pdf(sprintf("%s/upset_plots/01b_piechart_EmVars_context_number.pdf", FIGURE_DIR), height = 5, width = 4)
print(p)
dev.off()

N_Context_sugg <- emVars_any_suggestive[, .(nContext = .N), by = posID][, .N, keyby = nContext][, Pct := N / sum(N)][seq_len(.N)]
N_Context_sugg

# make a pie chart
N_Context_sugg[, label := paste0("  ", round(Pct * 100, ifelse(Pct > 0.02, 0, 1)), "%")]
#N_Context[, label := c(paste0("  ", round(Pct * 100), "%")[Pct > .03], "", "2.4%", "")]

N_Context_sugg[, x_pos := 1.65]
p <- ggplot(N_Context_sugg, aes(x = "", y = Pct, fill = factor(nContext))) +
  geom_bar(stat = "identity", width = 1, col = "black")
p <- p + coord_polar("y", start = 0) + scale_fill_brewer(palette = "Blues") + theme_void() + labs(fill = "Number of contexts\nwhere CRS is responsive")
p <- p + geom_text(aes(x = x_pos, label = label), position = position_stack(vjust = 0.5)) + theme(legend.position = "bottom")
pdf(sprintf("%s/upset_plots/01c_piechart_EmVars_context_number_Pemp05.pdf", FIGURE_DIR), height = 5, width = 4)
print(p)
dev.off()



N_celline <- emVars_any[, .(nCelltype = length(unique(celline))), by = posID][, .N, keyby = nCelltype][, Pct := N / sum(N)][seq_len(.N)]
N_celline

# make a pie chart
N_celline[, label := paste0("  ", round(Pct * 100, ifelse(Pct > 0.02, 0, 1)), "%")]
#N_Context[, label := c(paste0("  ", round(Pct * 100), "%")[Pct > .03], "", "2.4%", "")]

N_celline[, x_pos := 1.65]
p <- ggplot(N_celline, aes(x = "", y = Pct, fill = factor(nCelltype))) +
  geom_bar(stat = "identity", width = 1, col = "black")
p <- p + coord_polar("y", start = 0) + scale_fill_brewer(palette = "Blues") + theme_void() + labs(fill = "Number of contexts\nwhere CRS is responsive")
p <- p + geom_text(aes(x = x_pos, label = label), position = position_stack(vjust = 0.5)) + theme(legend.position = "bottom")
pdf(sprintf("%s/upset_plots/01d_piechart_EmVars_cellline_number.pdf", FIGURE_DIR), height = 5, width = 4)
print(p)
dev.off()

N_celline_sugg <- emVars_any_suggestive[, .(nCelltype = length(unique(celline))), by = posID][, .N, keyby = nCelltype][, Pct := N / sum(N)][seq_len(.N)]
N_celline_sugg

# make a pie chart
N_celline_sugg[, label := paste0("  ", round(Pct * 100, ifelse(Pct > 0.02, 0, 1)), "%")]
#N_Context[, label := c(paste0("  ", round(Pct * 100), "%")[Pct > .03], "", "2.4%", "")]

N_celline_sugg[, x_pos := 1.65]
p <- ggplot(N_celline_sugg, aes(x = "", y = Pct, fill = factor(nCelltype))) +
  geom_bar(stat = "identity", width = 1, col = "black")
p <- p + coord_polar("y", start = 0) + scale_fill_brewer(palette = "Blues") + theme_void() + labs(fill = "Number of contexts\nwhere CRS is responsive")
p <- p + geom_text(aes(x = x_pos, label = label), position = position_stack(vjust = 0.5)) + theme(legend.position = "bottom")

pdf(sprintf("%s/upset_plots/01e_piechart_EmVars_cellline_number_Pemp05.pdf", FIGURE_DIR), height = 5, width = 4)
print(p)
dev.off()

cat('\nAll done !\n')
q('no')
#############################################################################
#############################################################################
#############################################################################
