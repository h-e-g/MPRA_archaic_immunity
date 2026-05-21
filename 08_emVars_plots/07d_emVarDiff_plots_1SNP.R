MPRA_DIR <- "/pasteur/helix/projects/evo_immuno_pop/MPRA"
ANALYSIS_DIR <- "MPRA_count_exp6_analysis3"
SCRIPT_DIR <- sprintf("%s/scripts/%s/", MPRA_DIR, ANALYSIS_DIR)

source(sprintf("%s/scripts/%s/00_starter.R", MPRA_DIR, ANALYSIS_DIR))

RUN_ID <- "RUN3_Z2_nBC10"
CRITERIA_ACTIVE_OUT <- "FDR5_scrambled5pct_FC0.2_GCnorm"
CRITERIA_DIFF_OUT <- "noDiffFilter_Diff_FDR5_FC.2"
CRITERIA_EMVARS_OUT <- "EmVar_FDR5_FC.2"
CRITERIA_EMVARS_DIFF <- "EmVarDiff_FDR5_FC.2"
STIM_ONLY <- FALSE

cmd <- commandArgs(trailingOnly = TRUE)
print(cmd)

for (i in seq_along(cmd)) {
  if (cmd[i] == "--run_dir" || cmd[i] == "-r") {
    RUN_ID <- cmd[i + 1]
  }
  if (cmd[i] == "--criteria_active" || cmd[i] == "-a") {
    CRITERIA_ACTIVE_OUT <- cmd[i + 1]
  }
  if (cmd[i] == "--criteria_diff" || cmd[i] == "-d") {
    CRITERIA_DIFF_OUT <- cmd[i + 1]
  }
  if (cmd[i] == "--criteria_emvar" || cmd[i] == "-e") {
    CRITERIA_EMVARS_OUT <- cmd[i + 1]
  }
  if (cmd[i] == "--criteria_emvar_diff" || cmd[i] == "-ed") {
    CRITERIA_EMVARS_DIFF <- cmd[i + 1]
  }
  if (cmd[i] == "--stim_only" || cmd[i] == "-s") {
    STIM_ONLY <- as.logical(cmd[i + 1])
  }
}

CRITERIA_ACTIVE <- gsub("noActivityFilter_", "", CRITERIA_ACTIVE_OUT)
CRITERIA_DIFF <- gsub("noDiffFilter_", "", CRITERIA_DIFF_OUT)
CRITERIA_EMVARS <- gsub("noEmVarFilter_", "", CRITERIA_EMVARS_OUT)
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
SNP_annot_v5 <- merge(SNP_annot_v4[, -"Introgression_scenario"], selected_annot_wide, by = c("ID", "posID"))
toc()
# load annotations of SNPs
SNP_annot <- SNP_annot_v5
SNP_annot[, NON_INTROGRESSED := ifelse(INTROGRESSED.allele == ANCESTRAL, DERIVED, ANCESTRAL)]


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

library(ggplot2)

FIGURE_DIR <- sprintf("%s/figures/%s/%s", MPRA_DIR, ANALYSIS_DIR, RUN_ID)




oligo_activity_file <- sprintf("%s/oligo_activity__all__%s.tsv.gz", ACTIVE_DIR, CRITERIA_ACTIVE)
oligo_activity_class <- fread(oligo_activity_file)
oligo_activity_file <- sprintf("%s/all_oligos_annotated__%s.tsv.gz", ACTIVE_DIR, CRITERIA_ACTIVE)
oligo_activity_obs <- fread(oligo_activity_file)

emVARs_obs <- fread(file = sprintf("%s/all_emVars_annotated_celltype__%s.tsv", EMVAR_DIR, CRITERIA_ACTIVE_OUT))
emVARs_obs[, posID := crsID]
emVARs_obs <- emVARs_obs[is_emVar_CRITERIA == TRUE, ]

SNP_annot_emVars <- fread(file = sprintf("%s/SNP_annot_emVars__%s.tsv", EMVAR_DIR, CRITERIA_ACTIVE_OUT))

response_emVars <- fread(sprintf("%s/all_emVars_diff_obs_response.tsv", EMVARDIFF_DIR))
response_emVars[, posID := crsID]
response_emVars <- response_emVars[is_emVar_Diff_CRITERIA == TRUE]
response_emVars[,length(unique(posID))] # 94
# myPOSID <- response_emVars[1, posID]

celltype_emVars <- fread(sprintf("%s/all_emVars_diff_obs_celltype.tsv", EMVARDIFF_DIR))
celltype_emVars[, posID := crsID]
celltype_emVars <- celltype_emVars[is_emVar_Diff_CRITERIA == TRUE]


################################################################################################################################
################################ make emVar plot across stimuli/celltypes ######################################################
################################################################################################################################

plot_SNP <- function(POSID, select = c("all", "replicates"), GC_correct = TRUE, normalize = TRUE, target_cellines = NULL, target_conditions = NULL, ALLELES = NULL, use_introgressed_allele = TRUE, highlight_conditions = NULL, highlight_conditions_red = NULL) {
  # when specificied alleles should be fo the form 'ANC/DER'
  select <- match.arg(select)
  SNP_DT <- oligo_activity_obs[posID == POSID, ]
  SNP_DT_core <- SNP_DT[, .(oligo, ANALYSIS_NAME, allele, alpha, alpha.se, alpha.se_norm, alpha.GC, alpha.GC_norm, pval.SE = pval.SE_norm)]
  cat(".")
  # SNP_DT_core <- melt(SNP_DT_core, measure.vars = list(activity = c("activity_allele.1", "activity_allele.2"), allele = c("allele.1", "allele.2")), id.vars = c("analysis_name", "logFC.se"), variable.name = "allele.num")
  if (!GC_correct) {
    SNP_DT_core[, logactivity := log2(alpha)]
    SNP_DT_core[, logactivity.se := alpha.se]
  } else {
    if (normalize) {
      SNP_DT_core[, logactivity := log2(alpha.GC_norm)]
      SNP_DT_core[, logactivity.se := alpha.se_norm / log(2)]
    } else {
      SNP_DT_core[, logactivity := log2(alpha.GC)]
      SNP_DT_core[, logactivity.se := alpha.se / log(2)]
    }
  }
  SNP_DT_core[, logAct_hi := logactivity + 1.96 * logactivity.se]
  SNP_DT_core[, logAct_low := logactivity - 1.96 * logactivity.se]

  if (is.null(ALLELES)) {
    if (is.na(SNP_annot[posID == POSID, .(INTROGRESSED.allele)]) || SNP_annot[posID == POSID, .(INTROGRESSED.allele)] == "" || !use_introgressed_allele) {
      if (is.na(SNP_annot[posID == POSID, .(ANCESTRAL)]) || SNP_annot[posID == POSID, .(ANCESTRAL)] == "") {
        ALLELES <- unlist(unique(SNP_annot[posID == POSID, .(allele.1, allele.2)]))
        Allele_labels <- "Allele"
      } else {
        ALLELES <- unlist(unique(SNP_annot[posID == POSID & !is.na(ANCESTRAL), .(ANCESTRAL, DERIVED)]))
        Allele_labels <- "Allele (Ancestral/Derived)"
      }
    } else {
      ALLELES <- unlist(unique(SNP_annot[posID == POSID & !is.na(INTROGRESSED.allele), .(NON_INTROGRESSED, INTROGRESSED.allele)]))
      Allele_labels <- "Allele (Modern/Introgressed)"
    }
  } else {
    ALLELES <- unlist(str_split(ALLELES, "/"))
    Allele_labels <- "Allele"
  }

  if (select == "all") {
    SNP_DT_core <- SNP_DT_core[grepl("_all", ANALYSIS_NAME), ]
    SNP_DT_core <- merge(SNP_DT_core, condition_summary, by.x = "ANALYSIS_NAME", by.y = "analysis_name")
    SNP_DT_core[, ANALYSIS_NAME := gsub("_all", "", ANALYSIS_NAME)]
  } else {
    SNP_DT_core <- SNP_DT_core[grepl("_R", ANALYSIS_NAME), ]
    SNP_DT_core <- merge(SNP_DT_core, condition_summary_reps, by.x = "ANALYSIS_NAME", by.y = "analysis_name")
    SNP_DT_core[, ANALYSIS_NAME := gsub("-ACE2", "", ANALYSIS_NAME)]
  }
  if (!is.null(target_cellines)) {
    SNP_DT_core <- SNP_DT_core[celline %in% target_cellines, ]
  }
  if (!is.null(target_conditions)) {
    SNP_DT_core <- SNP_DT_core[condition %in% target_conditions, ]
  }

  if (any(nchar(ALLELES) > 3)) {
    max_nchar <- max(nchar(ALLELES))
    ALLELES[nchar(ALLELES) > 3] <- paste0(substr(ALLELES[nchar(ALLELES) > 3], 1, 1), "..", substr(ALLELES[nchar(ALLELES) > 3], max_nchar, max_nchar))
    SNP_DT_core[allele == ALLELES[nchar(ALLELES) > 3], allele := paste0(substr(allele, 1, 1), "..", substr(allele, nchar(allele), nchar(allele)))]
  }
  # print(SNP_DT_core)
  # print(ALLELES)

  if (select == "all") {
    p <- ggplot(SNP_DT_core, aes(
      x = factor(allele, ALLELES), y = logactivity,
      ymin = logAct_low,
      ymax = logAct_hi,
      col = COND_ID
    ))
    p <- p + geom_pointrange(alpha = .5)
    p <- p + scale_color_manual(values = color_setup_simplified_norep)
  } else {
    p <- ggplot(SNP_DT_core, aes(
      x = factor(allele, ALLELES), y = logactivity,
      ymin = logAct_low,
      ymax = logAct_hi,
      col = ANALYSIS_NAME
    ))
    p <- p + geom_pointrange(alpha = .5, position = position_jitter(height = 0, width = .2), size = .3)
    p <- p + scale_color_manual(values = color_setup_simplified)
  }
  p <- p + facet_grid(cols = vars(celline, factor(condition, c("NS", "IAV", "SARS", "IFNA2b", "DEX", "TNFa"))))
  p <- p + theme_plot(rotate.x = ifelse(max(nchar(ALLELES)) > 1, 90, 0)) + xlab(Allele_labels)
  p <- p + geom_hline(yintercept = 0, linetype = 2, col = "lightgrey") + guides(color = "none")
  p <- p + ylab("log2( CRS activity )")

  if (!is.null(highlight_conditions)) {
    stars_data <- SNP_DT_core[COND_ID %in% highlight_conditions, # Filter for specific facets
      .(x = 1.5, y = max(logactivity) + 0.1, label = "*", is_red = COND_ID %in% highlight_conditions_red), # Summarize with max(logactivity)
      by = .(COND_ID, condition, celline)
    ] # Group by condition and celline

    # add stars to significant facets
    p <- p + geom_text(data = stars_data[is_red == FALSE, ], aes(x = x, y = y, label = label), inherit.aes = FALSE, size = 6, col = "black")
    p <- p + geom_text(data = stars_data[is_red == TRUE, ], aes(x = x, y = y, label = label), inherit.aes = FALSE, size = 6, col = "red")
  }
  p
}

FIGURE_EMVARDIFF_DIR <- sprintf("%s/figures/%s/%s/05_emVarDiff/%s/%s/%s/%s/%s", MPRA_DIR, ANALYSIS_DIR, RUN_ID, CRITERIA_EMVARS_DIFF, CRITERIA_EMVARS_OUT, CRITERIA_DIFF_OUT, CRITERIA_ACTIVE_OUT, STIM_ONLY_OUT)

dir.create(sprintf("%s/response_emVar_plots/all_conditions", FIGURE_EMVARDIFF_DIR), showWarnings = FALSE, recursive = TRUE)
dir.create(sprintf("%s/response_emVar_plots/HepG2_only", FIGURE_EMVARDIFF_DIR), showWarnings = FALSE, recursive = TRUE)
dir.create(sprintf("%s/response_emVar_plots/A549_only", FIGURE_EMVARDIFF_DIR), showWarnings = FALSE, recursive = TRUE)
dir.create(sprintf("%s/response_emVar_plots/K562_only", FIGURE_EMVARDIFF_DIR), showWarnings = FALSE, recursive = TRUE)

tic("plotting response emVars")
for (myPOSID in response_emVars[, unique(posID)]) {
  print(myPOSID)
  cat("\n", myPOSID, "\n")
  myRSID <- try(SNP_annot[posID == myPOSID, get("rsID")])
  highlight_conds <- emVARs_obs[posID == myPOSID, unique(COND_ID)]
  highlight_conds_red <- response_emVars[posID == myPOSID, unique(COND_ID)]
  if (class(myRSID) == "try-error" || is.na(myRSID) || myRSID == "") {
    myVariantID <- myPOSID
    rsID_print <- gsub(":", "-", myPOSID)
  } else {
    myVariantID <- paste0(unique(myRSID), " (", myPOSID, ")")
    rsID_print <- paste(gsub(":", "-", myPOSID), unique(myRSID), sep = "_")
  }
  # criteria_list <- unlist(str_split(SNP_annot[posID == myPOSID, selection_criteria], ","))
  # for (type.CRS in criteria_list) {
  # dir.create(sprintf("%s/2_emVar_plots/01_response_emVars_active_signif/%s", FIGURE_DIR, type.CRS), showWarnings = FALSE, recursive = TRUE)
  # all conditions
  pdf(sprintf("%s/response_emVar_plots/all_conditions/variant_%s_01a_SNPeffect.pdf", FIGURE_EMVARDIFF_DIR, rsID_print), height = 2, width = 4)
  p <- try(plot_SNP(POSID = myPOSID, GC_correct = TRUE, select = "all", highlight_conditions = highlight_conds, highlight_conditions_red = highlight_conds_red))
  try(print(p))
  dev.off()

  pdf(sprintf("%s/response_emVar_plots/all_conditions/variant_%s_01b_SNPeffect_reps.pdf", FIGURE_EMVARDIFF_DIR, rsID_print), height = 2, width = 4)
  p <- try(plot_SNP(POSID = myPOSID, GC_correct = TRUE, select = "replicates", highlight_conditions = highlight_conds, highlight_conditions_red = highlight_conds_red))
  try(print(p))
  dev.off()

  if ("HepG2" %in% response_emVars[posID == myPOSID, celline]) {
    pdf(sprintf("%s/response_emVar_plots/HepG2_only/variant_%s_02a_SNPeffect.pdf", FIGURE_EMVARDIFF_DIR, rsID_print), height = 2, width = 1.8)
    p <- try(plot_SNP(POSID = myPOSID, GC_correct = TRUE, select = "all", highlight_conditions = highlight_conds, highlight_conditions_red = highlight_conds_red, target_cellines = "HepG2"))
    try(print(p))
    dev.off()

    pdf(sprintf("%s/response_emVar_plots/HepG2_only/variant_%s_02b_SNPeffect_reps.pdf", FIGURE_EMVARDIFF_DIR, rsID_print), height = 2, width = 1.8)
    p <- try(plot_SNP(POSID = myPOSID, GC_correct = TRUE, select = "replicates", highlight_conditions = highlight_conds, highlight_conditions_red = highlight_conds_red, target_cellines = "HepG2"))
    try(print(p))
    dev.off()
  }

  # A549
  if ("A549" %in% response_emVars[posID == myPOSID, celline]) {
    pdf(sprintf("%s/response_emVar_plots/A549_only/variant_%s_02a_SNPeffect.pdf", FIGURE_EMVARDIFF_DIR, rsID_print), height = 2, width = 1.8)
    p <- try(plot_SNP(POSID = myPOSID, GC_correct = TRUE, select = "all", highlight_conditions = highlight_conds, highlight_conditions_red = highlight_conds_red, target_cellines = "A549"))
    try(print(p))
    dev.off()


    pdf(sprintf("%s/response_emVar_plots/A549_only/variant_%s_02b_SNPeffect_reps.pdf", FIGURE_EMVARDIFF_DIR, rsID_print), height = 2, width = 1.8)
    p <- try(plot_SNP(POSID = myPOSID, GC_correct = TRUE, select = "replicates", highlight_conditions = highlight_conds, highlight_conditions_red = highlight_conds_red, target_cellines = "A549"))
    try(print(p))
    dev.off()
  }
  # K562
  if ("K562" %in% response_emVars[posID == myPOSID, celline]) {
    pdf(sprintf("%s/response_emVar_plots/K562_only/variant_%s_02a_SNPeffect.pdf", FIGURE_EMVARDIFF_DIR, rsID_print), height = 2, width = 1.8)
    p <- try(plot_SNP(POSID = myPOSID, GC_correct = TRUE, select = "all", highlight_conditions = highlight_conds, highlight_conditions_red = highlight_conds_red, target_cellines = "K562"))
    try(print(p))
    dev.off()

    pdf(sprintf("%s/response_emVar_plots/K562_only/variant_%s_02b_SNPeffect_reps.pdf", FIGURE_EMVARDIFF_DIR, rsID_print), height = 2, width = 1.8)
    p <- try(plot_SNP(POSID = myPOSID, GC_correct = TRUE, select = "replicates", highlight_conditions = highlight_conds, highlight_conditions_red = highlight_conds_red, target_cellines = "K562"))
    try(print(p))
    dev.off()
  }
}
toc()

dir.create(sprintf("%s/celltype_emVar_plots/all_conditions", FIGURE_EMVARDIFF_DIR), showWarnings = FALSE, recursive = TRUE)
dir.create(sprintf("%s/celltype_emVar_plots/NS_only", FIGURE_EMVARDIFF_DIR), showWarnings = FALSE, recursive = TRUE)

tic("plotting celltype emVars")
for (myPOSID in celltype_emVars[, unique(posID)]) {
  print(myPOSID)
  cat("\n", myPOSID, "\n")
  myRSID <- try(SNP_annot[posID == myPOSID, get("rsID")])
  highlight_conds <- emVARs_obs[posID == myPOSID, unique(COND_ID)]
  if (class(myRSID) == "try-error" || is.na(myRSID) || myRSID == "") {
    myVariantID <- myPOSID
    rsID_print <- gsub(":", "-", myPOSID)
  } else {
    myVariantID <- paste0(unique(myRSID), " (", myPOSID, ")")
    rsID_print <- paste(gsub(":", "-", myPOSID), unique(myRSID), sep = "_")
  }
  # all conditions
  pdf(sprintf("%s/celltype_emVar_plots/all_conditions/variant_%s_01a_SNPeffect.pdf", FIGURE_EMVARDIFF_DIR, rsID_print), height = 2, width = 4)
  p <- try(plot_SNP(POSID = myPOSID, GC_correct = TRUE, select = "all", highlight_conditions = highlight_conds))
  try(print(p))
  dev.off()

  pdf(sprintf("%s/celltype_emVar_plots/all_conditions/variant_%s_01b_SNPeffect_reps.pdf", FIGURE_EMVARDIFF_DIR, rsID_print), height = 2, width = 4)
  p <- try(plot_SNP(POSID = myPOSID, GC_correct = TRUE, select = "replicates", highlight_conditions = highlight_conds))
  try(print(p))
  dev.off()

  # NS only
  if (any(celltype_emVars[posID == myPOSID, grepl("NS", group1_labels) & grepl("NS", group2_labels)])) {
    pdf(sprintf("%s/celltype_emVar_plots/NS_only/variant_%s_02a_SNPeffect.pdf", FIGURE_EMVARDIFF_DIR, rsID_print), height = 2, width = 4)
    p <- try(plot_SNP(POSID = myPOSID, GC_correct = TRUE, select = "all", highlight_conditions = highlight_conds, target_conditions = "NS"))
    try(print(p))
    dev.off()

    pdf(sprintf("%s/celltype_emVar_plots/NS_only/variant_%s_02b_SNPeffect_reps.pdf", FIGURE_EMVARDIFF_DIR, rsID_print), height = 2, width = 4)
    p <- try(plot_SNP(POSID = myPOSID, GC_correct = TRUE, select = "all", highlight_conditions = highlight_conds, target_conditions = "NS"))
    try(print(p))
    dev.off()
  }
}
toc()

cat("\nAll done!\n")
q("no")
