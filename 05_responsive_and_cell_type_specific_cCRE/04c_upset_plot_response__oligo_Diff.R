# 04b_upset_plot_oligo_Diff.R

# running: sbatch -p geh,common --mem=30G  00_Rscript.sh MPRA_count_exp6_analysisZ/04b_upset_plot_oligo_Diff.R

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


####################################################################
################# compute FDR, CRS difference ######################
####################################################################

oligo_activity_Diff_obs <- fread(sprintf("%s/all_oligos_diff_annotated__%s.tsv.gz", DIFF_DIR, CRITERIA_DIFF))

oligo_response <- oligo_activity_Diff_obs[ANALYSIS_SUBTYPE == "response", ]
oligo_response[, COND_ID := gsub("-ACE2", "", gsub("response_(.*)", "\\1", ANALYSIS_NAME))]

oligo_response <- merge(oligo_response, condition_summary, by = "COND_ID")
oligo_response <- merge(oligo_response, SNP_annot[, .(posID, POP_adaptive, Introgressed_from = Adaptive_from, Introgression_scenario_v2)], by = "posID")

oligo_response_any <- unique(oligo_response[oligo_diff_CRITERIA == TRUE, .(posID, celline, condition, type, POP_adaptive, Introgressed_from)])
oligo_response_increased <- unique(oligo_response[oligo_diff_CRITERIA == TRUE & log2FC_2vs1 > 0, .(posID, celline, condition, type, POP_adaptive, Introgressed_from)])
oligo_response_decreased <- unique(oligo_response[oligo_diff_CRITERIA == TRUE & log2FC_2vs1 < 0, .(posID, celline, condition, type, POP_adaptive, Introgressed_from)])

N_Context <- oligo_response_any[, .(nContext = .N), by = posID][, .N, keyby = nContext][, Pct := N / sum(N)][seq_len(.N)]
N_Context

oligo_response_suggestive <- unique(oligo_response[posID %in% oligo_response_any$posID & pval_emp < 0.05, .(posID, celline, condition, type, POP_adaptive, Introgressed_from)])
N_Context_Pemp <- oligo_response_suggestive[, .(nContext = .N), by = posID][, .N, keyby = nContext][, Pct := N / sum(N)][seq_len(.N)]
N_Context_Pemp

oligo_response_noViral <- unique(oligo_response[oligo_diff_CRITERIA == TRUE & !condition %in% c("SARS", "IAV"), .(posID, celline, condition, type, POP_adaptive, Introgressed_from)])
N_Context_noViral <- oligo_response_noViral[, .(nContext = .N), by = posID][, .N, keyby = nContext][, Pct := N / sum(N)][seq_len(.N)]
N_Context_noViral

oligo_response_suggestive_noViral <- unique(oligo_response[posID %in% oligo_response_noViral[, posID] & pval_emp < 0.05 & !condition %in% c("SARS", "IAV"), .(posID, celline, condition, type, POP_adaptive, Introgressed_from)])
N_Context_Pemp_noViral <- oligo_response_suggestive_noViral[, .(nContext = .N), by = posID][, .N, keyby = nContext][, Pct := N / sum(N)][seq_len(.N)]
N_Context_Pemp_noViral

# oligo_response_strong <- unique(oligo_response[oligo_diff_CRITERIA == TRUE & lfdr < 0.01, .(posID, celline, condition, type, POP_adaptive, Introgressed_from)])
# oligo_response_suggestive_strong <- unique(oligo_response[posID %in% oligo_response_strong$posID & pval.LRT < 0.05, .(posID, celline, condition, type, POP_adaptive, Introgressed_from)])
# N_Context_lfdr1 <- oligo_response_suggestive_strong[, .(nContext = .N), by = posID][, .N, keyby = nContext][, Pct := N / sum(N)][seq_len(.N)]
# N_Context_lfdr1


# oligo_response_strong1 <- unique(oligo_response[oligo_diff_CRITERIA == TRUE & FDR < 0.01, .(posID, celline, condition, type, POP_adaptive, Introgressed_from)])
# N_Context_FDR1 <- oligo_response_strong1[, .(nContext = .N), by = posID][, .N, keyby = nContext][, Pct := N / sum(N)][seq_len(.N)]
# N_Context_FDR1

# oligo_response_strong01 <- unique(oligo_response[oligo_diff_CRITERIA == TRUE & FDR < 0.001, .(posID, celline, condition, type, POP_adaptive, Introgressed_from)])
# N_Context_FDR01 <- oligo_response_strong01[, .(nContext = .N), by = posID][, .N, keyby = nContext][, Pct := N / sum(N)][seq_len(.N)]
# N_Context_FDR01

# oligo_response_suggestive_strong1 <- unique(oligo_response[posID %in% oligo_response_strong1$posID & pval_emp < 0.05, .(posID, celline, condition, type, POP_adaptive, Introgressed_from)])
# N_Context_FDR1_Pemp <- oligo_response_suggestive_strong1[, .(nContext = .N), by = posID][, .N, keyby = nContext][, Pct := N / sum(N)][seq_len(.N)]
# N_Context_FDR1_Pemp

# oligo_response_suggestive_strong1 <- unique(oligo_response[posID %in% sample(posID, length(oligo_response_strong1$posID)) & pval_emp < 0.05, .(posID, celline, condition, type, POP_adaptive, Introgressed_from)])
# N_Context_FDR1_Pemp <- oligo_response_suggestive_strong1[, .(nContext = .N), by = posID][, .N, keyby = nContext][, Pct := N / sum(N)][seq_len(.N)]
# N_Context_FDR1_Pemp


# oligo_response_strong <- unique(oligo_response[oligo_diff_CRITERIA == TRUE & FDR < 0.001, .(posID, celline, condition, type, POP_adaptive, Introgressed_from)])
# oligo_response_suggestive_strong <- unique(oligo_response[posID %in% oligo_response_strong$posID & pval_emp < 0.05, .(posID, celline, condition, type, POP_adaptive, Introgressed_from)])


# oligo_response_suggestive <- unique(oligo_response[posID %in% oligo_response_any$posID & pval.LRT < 0.1, .(posID, celline, condition, type, POP_adaptive, Introgressed_from)])
# N_Context_suggestive <- oligo_response_suggestive[, .(nContext = .N), by = posID][, .N, keyby = nContext][, Pct := N / sum(N)][seq_len(.N)]
# N_Context_suggestive

####################################################################
######################## uspet plots, oligo Diff ###################
####################################################################

# library(MPRAnalyze)
library(ComplexUpset)
library(ggplot2)
library(dplyr)
library(data.table)

names_stimulation <- names(color_setup_simplified_norep)[!grepl("_NS$", names(color_setup_simplified_norep))]

create_upset <- function(fig_data, ...) {
  fig_data$Archaic_origin <- case_when(
    fig_data$Introgressed_from == "Vindija" ~ "Vindija",
    fig_data$Introgressed_from == "Denisova" ~ "Denisova",
    fig_data$Introgressed_from == "Archaic (shared)" ~ "Both genomes",
    fig_data$Introgressed_from == "Archaic (unknown)" ~ "Not found in either genome\n (indels or unknown archaic)",
    TRUE ~ "WTF"
  )
  archaic_color <- c("AMH" = grey(0.8), "Denisova" = "#80CDC1", "Vindija" = "#FDAE61", "Both genomes" = "#FEE08B", "Not found in either genome\n (indels or unknown archaic)" = grey(0.5))

  fig_data$Archaic_origin <- factor(fig_data$Archaic_origin, c("Vindija", "Denisova", "Both genomes", "Not found in either genome\n (indels or unknown archaic)"))
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
    "Intersection size" = intersection_size() # REMOVE mapping
  )
  my_queries <- lapply(seq_along(color_setup_simplified_norep[names_stimulation]), function(i) {
    upset_query(set = names_stimulation[i], fill = color_setup_simplified_norep[names_stimulation[i]])
  })

  valid_sets <- names_stimulation[colSums(fig_data[, mget(names_stimulation)]) > 0]
  valid_colors <- color_setup_simplified_norep[valid_sets]
  my_queries <- lapply(seq_along(valid_sets), function(i) {
    upset_query(set = valid_sets[i], fill = valid_colors[[valid_sets[i]]])
  })

  p <- upset(as.data.frame(fig_data),
    colnames(fig_data),
    name = "celltype & condition",
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


data_names <- c("oligo_response_any", "oligo_response_increased", "oligo_response_decreased", "oligo_response_suggestive", "oligo_response_noViral", "oligo_response_suggestive_noViral")

data_to_plot <- lapply(seq_along(data_names), function(i) {
  data_name <- data_names[i]
  data <- get(data_name)
  #  data <- data[!duplicated(posID),]
  data <- dcast(data, posID + Introgressed_from ~ factor(paste(celline, condition, sep = "_"), names_stimulation), fill = 0, fun.aggregate = length, drop = c(TRUE, FALSE))
  data <- cbind(data[, c("posID", "Introgressed_from")], data[, mget(rev(names_stimulation))] > 0)
  rownames(data) <- data$posID
  data
})

for (i in seq_along(data_to_plot)) {
  fig_data <- data_to_plot[[i]]
  cat(data_names[i], "\n")
  p <- create_upset(unique(fig_data))

  pdf(sprintf("%s/02a_upsetPlot_responseCRS_%s.pdf", FIGURE_DIR, data_names[i]), height = 5, width = 8)
  print(p)
  dev.off()
}

data_to_plot <- lapply(seq_along(data_names), function(i) {
  data_name <- data_names[i]
  data <- get(data_name)[celline != "Calu3", ]
  data <- dcast(data, posID + POP_adaptive + Introgressed_from ~ celline, fill = 0, fun.aggregate = length)
  data <- cbind(data[, c("posID", "POP_adaptive", "Introgressed_from")], data[, mget(rev(names(color_celline[names(color_celline) != "Calu3"])))] > 0)
  rownames(data) <- data$posID
  data
})


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


# make a pie chart
N_Context[, label := paste0("  ", round(Pct * 100, ifelse(Pct > 0.02, 0, 1)), "%")]
# N_Context[, label := c(paste0("  ", round(Pct * 100), "%")[Pct > .03], "", "2.4%", "")]

N_Context[, x_pos := 1.65 + rep(0, .N)]
p <- ggplot(N_Context, aes(x = "", y = Pct, fill = factor(nContext))) +
  geom_bar(stat = "identity", width = 1, col = "black")
p <- p + coord_polar("y", start = 0) + scale_fill_brewer(palette = "Blues") + theme_void() + labs(fill = "Number of contexts\nwhere CRS is responsive")
p <- p + geom_text(aes(x = x_pos, label = label), position = position_stack(vjust = 0.5)) + theme(legend.position = "bottom")
pdf(sprintf("%s/02b_piechart_responseCRS_context_number.pdf", FIGURE_DIR), height = 5, width = 4)
print(p)
  dev.off()

  N_Context[, label := paste0("  ", round(Pct * 100, ifelse(Pct > 0.02, 0, 1)), "%")]
  # N_Context[, label := c(paste0("  ", round(Pct * 100), "%")[Pct > .03], "", "2.4%", "")]
  N_Context[, x_pos := 1.65 + rep(0, .N)]
  p <- ggplot(N_Context, aes(x = "", y = Pct, fill = factor(nContext))) +
    geom_bar(stat = "identity", width = 1, col = "black")
  p <- p + coord_polar("y", start = 0) + scale_fill_brewer(palette = "Blues") + theme_void() + labs(fill = "Number of\ncontexts\nwhere CRE\nis responsive")
  p <- p + geom_text(aes(x = x_pos, label = label), position = position_stack(vjust = 0.8), size = 5) + theme(legend.position = "right")
  p <- p + theme( legend.text = element_text(size = 12),  legend.title = element_text(size = 12, hjust = 0), plot.margin = margin(1, 1, 1, 1, "cm"))
  pdf(sprintf("%s/02b_piechart_responseCRS_context_number_larger.pdf", FIGURE_DIR), height = 4, width = 6)
  print(p)
  dev.off()


# make a pie chart
N_Context_Pemp[, label := paste0("  ", round(Pct * 100, ifelse(Pct > 0.02, 0, 1)), "%")]
# N_Context[, label := c(paste0("  ", round(Pct * 100), "%")[Pct > .03], "", "2.4%", "")]
N_Context_Pemp[, x_pos := 1.65 + rep(0, .N)]
p <- ggplot(N_Context_Pemp, aes(x = "", y = Pct, fill = factor(nContext))) +
  geom_bar(stat = "identity", width = 1, col = "black")
p <- p + coord_polar("y", start = 0) + scale_fill_brewer(palette = "Blues") + theme_void() + labs(fill = "Number of contexts\nwhere CRS is responsive\n (Pemp < 0.05)")
p <- p + geom_text(aes(x = x_pos, label = label), position = position_stack(vjust = 0.5)) + theme(legend.position = "bottom")
pdf(sprintf("%s/02b_piechart_responseCRS_context_number_suggestive.pdf", FIGURE_DIR), height = 5, width = 4)
print(p)
dev.off()



# make a pie chart
N_Context_Pemp[, label := paste0("  ", round(Pct * 100, ifelse(Pct > 0.02, 0, 1)), "%")]
# N_Context[, label := c(paste0("  ", round(Pct * 100), "%")[Pct > .03], "", "2.4%", "")]
N_Context_Pemp[, x_pos := 1.65 + rep(0, .N)]
p <- ggplot(N_Context_Pemp, aes(x = "", y = Pct, fill = factor(nContext))) +
  geom_bar(stat = "identity", width = 1, col = "black")
p <- p + coord_polar("y", start = 0) + scale_fill_brewer(palette = "Blues") + theme_void() + labs(fill = "Number of\ncontexts\nwhere CRE\nis responsive\n(Pemp < 0.05)")
p <- p + geom_text(aes(x = x_pos, label = label), position = position_stack(vjust = 0.5), size = 5) + theme(legend.position = "right")
p <- p + theme( legend.text = element_text(size = 12),  legend.title = element_text(size = 12, hjust = 0), plot.margin = margin(1, 1, 1, 1, "cm"))
pdf(sprintf("%s/02b_piechart_responseCRS_context_number_suggestivelarger.pdf", FIGURE_DIR), height = 4, width = 6)
print(p)
dev.off()


# make a pie chart
N_Context_noViral[, label := paste0("  ", round(Pct * 100, ifelse(Pct > 0.02, 0, 1)), "%")]
# N_Context[, label := c(paste0("  ", round(Pct * 100), "%")[Pct > .03], "", "2.4%", "")]

N_Context_noViral[, x_pos := 1.65 + rep(0, .N)]
p <- ggplot(N_Context_noViral, aes(x = "", y = Pct, fill = factor(nContext))) +
  geom_bar(stat = "identity", width = 1, col = "black")
p <- p + coord_polar("y", start = 0) + scale_fill_brewer(palette = "Blues") + theme_void() + labs(fill = "Number of contexts\nwhere CRS is responsive")
p <- p + geom_text(aes(x = x_pos, label = label), position = position_stack(vjust = 0.5)) + theme(legend.position = "bottom")
pdf(sprintf("%s/02b_piechart_responseCRS_context_number_noViral.pdf", FIGURE_DIR), height = 5, width = 4)
print(p)
dev.off()



# make a pie chart
N_Context_Pemp_noViral[, label := paste0("  ", round(Pct * 100, ifelse(Pct > 0.02, 0, 1)), "%")]
# N_Context[, label := c(paste0("  ", round(Pct * 100), "%")[Pct > .03], "", "2.4%", "")]
N_Context_Pemp_noViral[, x_pos := 1.65 + rep(0, .N)]
p <- ggplot(N_Context_Pemp_noViral, aes(x = "", y = Pct, fill = factor(nContext))) +
  geom_bar(stat = "identity", width = 1, col = "black")
p <- p + coord_polar("y", start = 0) + scale_fill_brewer(palette = "Blues") + theme_void() + labs(fill = "Number of contexts\nwhere CRS is responsive\n (Pemp < 0.05)")
p <- p + geom_text(aes(x = x_pos, label = label), position = position_stack(vjust = 0.5)) + theme(legend.position = "bottom")
pdf(sprintf("%s/02b_piechart_responseCRS_context_number_suggestive_noViral.pdf", FIGURE_DIR), height = 5, width = 4)
print(p)
dev.off()



# # make a pie chart
# N_Context_FDR1[, label := paste0("  ", round(Pct * 100, ifelse(Pct > 0.02, 0, 1)), "%")]
# # N_Context[, label := c(paste0("  ", round(Pct * 100), "%")[Pct > .03], "", "2.4%", "")]

# N_Context_FDR1[, x_pos := 1.65 + rep(0, .N)]
# p <- ggplot(N_Context_FDR1, aes(x = "", y = Pct, fill = factor(nContext))) +
#   geom_bar(stat = "identity", width = 1, col = "black")
# p <- p + coord_polar("y", start = 0) + scale_fill_brewer(palette = "Blues") + theme_void() + labs(fill = "Number of contexts\nwhere CRS is responsive")
# p <- p + geom_text(aes(x = x_pos, label = label), position = position_stack(vjust = 0.5)) + theme(legend.position = "bottom")
# pdf(sprintf("%s/02b_piechart_responseCRS_context_number__FDR1.pdf", FIGURE_DIR), height = 5, width = 4)
# print(p)
# dev.off()


# # make a pie chart
# N_Context_lfdr1[, label := paste0("  ", round(Pct * 100, ifelse(Pct > 0.02, 0, 1)), "%")]
# # N_Context[, label := c(paste0("  ", round(Pct * 100), "%")[Pct > .03], "", "2.4%", "")]

# N_Context_lfdr1[, x_pos := 1.65 + rep(0, .N)]
# p <- ggplot(N_Context_lfdr1, aes(x = "", y = Pct, fill = factor(nContext))) +
#   geom_bar(stat = "identity", width = 1, col = "black")
# p <- p + coord_polar("y", start = 0) + scale_fill_brewer(palette = "Blues") + theme_void() + labs(fill = "Number of contexts\nwhere CRS is responsive")
# p <- p + geom_text(aes(x = x_pos, label = label), position = position_stack(vjust = 0.5)) + theme(legend.position = "bottom")
# pdf(sprintf("%s/02b_piechart_responseCRS_context_number__lfdr1.pdf", FIGURE_DIR), height = 5, width = 4)
# print(p)
# dev.off()


# library(ggplot2)
# library(dplyr)
# library(data.table)

# make a pie chart
# # N_Context[,label:=paste0("  ",round(Pct*100,ifelse(Pct>0.02,0,1)),"%")]
# N_Context[, label := c(paste0("  ", round(Pct * 100), "%")[Pct > .03], "", "2.4%", "")]

# N_Context[, x_pos := 1.65 + c(0, 0, 0, 0, 0, 0, 0)]
# p <- ggplot(N_Context, aes(x = "", y = Pct, fill = factor(nContext))) +
#   geom_bar(stat = "identity", width = 1, col = "black")
# p <- p + coord_polar("y", start = 0) + scale_fill_brewer(palette = "Blues") + theme_void() + labs(fill = "Number of contexts\nwhere CRS is responsive")
# p <- p + geom_text(aes(x = x_pos, label = label), position = position_stack(vjust = 0.5)) + theme(legend.position = "bottom")
# pdf(sprintf("%s/02b_piechart_responseCRS_context_number.pdf", FIGURE_DIR), height = 5, width = 4)
# print(p)
# dev.off()
