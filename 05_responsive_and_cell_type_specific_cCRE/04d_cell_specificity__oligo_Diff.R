# 04b_upset_plot_oligo_Diff.R

# running: sbatch -p geh,common --mem=30G  00_Rscript.sh MPRA_count_exp6_analysisZ/04b_upset_plot_oligo_Diff.R

MPRA_DIR <- "/pasteur/helix/projects/evo_immuno_pop/MPRA"
ANALYSIS_DIR <- "MPRA_count_exp6_analysis3"

source(sprintf("%s/scripts/%s/00_starter.R", MPRA_DIR, ANALYSIS_DIR))
source(sprintf("%s/scripts/%s/00_load_associations.R", MPRA_DIR, ANALYSIS_DIR))
source(sprintf("%s/scripts/misc_plots.R", MPRA_DIR))

RUN_ID <- "RUN3_Z2_nBC10"
#CRITERIA_ACTIVE <- "FDR5_scrambled5pct_FC0.2_GCnorm"
CRITERIA_ACTIVE <- "FDR5_FC1_FC0.2_GCnorm"
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
   CRITERIA_ACTIVE_OUT <- paste0('noActivityFilter_',CRITERIA_ACTIVE)
 }else{
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
FIGURE_DIR <- sprintf("%s/figures/%s/%s/03_oligo_diff/%s/%s", MPRA_DIR, ANALYSIS_DIR, RUN_ID, CRITERIA_DIFF, CRITERIA_ACTIVE_OUT )
dir.create(FIGURE_DIR, recursive = TRUE)

# load oligo annotations
source(sprintf("%s/scripts/%s/02z__define_includedCRS.R", MPRA_DIR, ANALYSIS_DIR))


########################################################################################################
########################## Definition of cell type specific sites, Sup tables 2a-C #####################
########################################################################################################
# create clean tables of cell type specific oligos

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
oligo_activity_small <- oligo_activity_obs_ctc[condition == "NS", .(oligo, crsID, celline, oligo_class_CRITERIA, alpha_CRITERIA)]
oligo_activity_small <- dcast(oligo_activity_small, crsID + oligo ~ celline, value.var = c("oligo_class_CRITERIA", "alpha_CRITERIA"))

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

cre_A549_up= oligo_activity_obs_ctc[oligo%in%oligo_A549_up,unique(crsID)]
cre_HepG2_up= oligo_activity_obs_ctc[oligo%in%oligo_HepG2_up,unique(crsID)]
cre_K562_up= oligo_activity_obs_ctc[oligo%in%oligo_K562_up,unique(crsID)]

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

cre_A549_lo= oligo_activity_obs_ctc[oligo%in%oligo_A549_lo,unique(crsID)]
cre_HepG2_lo= oligo_activity_obs_ctc[oligo%in%oligo_HepG2_lo,unique(crsID)]
cre_K562_lo= oligo_activity_obs_ctc[oligo%in%oligo_K562_lo,unique(crsID)]

#### identify oligos that are cell type specific (upregulated)
oligo_K562_up_sugg <- intersect(
  oligo_celltype_comp_NS_annot_obs_archaic[ANALYSIS_NAME == "Celltype_HepG2_K562_NS" & pval_emp < 0.05 & log2FC_2vs1 > 0, oligo],
  oligo_celltype_comp_NS_annot_obs_archaic[ANALYSIS_NAME == "Celltype_K562_A549_NS" & pval_emp < 0.05 & log2FC_2vs1 < 0, oligo]
)
oligo_HepG2_up_sugg <- intersect(
  oligo_celltype_comp_NS_annot_obs_archaic[ANALYSIS_NAME == "Celltype_HepG2_K562_NS" & pval_emp < 0.05 & log2FC_2vs1 < 0, oligo],
  oligo_celltype_comp_NS_annot_obs_archaic[ANALYSIS_NAME == "Celltype_HepG2_A549_NS" & pval_emp < 0.05 & log2FC_2vs1 < 0, oligo]
)
oligo_A549_up_sugg <- intersect(
  oligo_celltype_comp_NS_annot_obs_archaic[ANALYSIS_NAME == "Celltype_K562_A549_NS" & pval_emp < 0.05 & log2FC_2vs1 > 0, oligo],
  oligo_celltype_comp_NS_annot_obs_archaic[ANALYSIS_NAME == "Celltype_HepG2_A549_NS" & pval_emp < 0.05 & log2FC_2vs1 > 0, oligo]
)

cre_A549_up_sugg= oligo_activity_obs_ctc[oligo%in%oligo_A549_up_sugg,unique(crsID)]
cre_HepG2_up_sugg= oligo_activity_obs_ctc[oligo%in%oligo_HepG2_up_sugg,unique(crsID)]
cre_K562_up_sugg= oligo_activity_obs_ctc[oligo%in%oligo_K562_up_sugg,unique(crsID)]

### EXTRACT STATS FOR CELL TYPE SPECIFIC OLIGOS (UPREGULATED)
oligo_K562_specific_up <- oligo_celltype_comp_NS_annot_obs_archaic[crsID %in% cre_K562_up, .(oligo, crsID, allele.label, group1_labels, group2_labels, FDR, pval_emp, log2FC_2vs1, log2FC.se, pval.LRT, oligo_diff_CRITERIA)]
oligo_HepG2_specific_up <- oligo_celltype_comp_NS_annot_obs_archaic[crsID %in% cre_HepG2_up, .(oligo, crsID, allele.label, group1_labels, group2_labels, FDR, pval_emp, log2FC_2vs1, log2FC.se, pval.LRT, oligo_diff_CRITERIA)]
oligo_A549_specific_up <- oligo_celltype_comp_NS_annot_obs_archaic[crsID %in% cre_A549_up, .(oligo, crsID, allele.label, group1_labels, group2_labels, FDR, pval_emp, log2FC_2vs1, log2FC.se, pval.LRT, oligo_diff_CRITERIA)]

oligo_K562_specific_up <- merge(oligo_K562_specific_up, oligo_activity_small[crsID %in% cre_K562_up, ], by =  c("oligo",'crsID'), allow.cartesian = TRUE)
oligo_HepG2_specific_up <- merge(oligo_HepG2_specific_up, oligo_activity_small[crsID %in% cre_HepG2_up, ], by =  c("oligo",'crsID'), allow.cartesian = TRUE)
oligo_A549_specific_up <- merge(oligo_A549_specific_up, oligo_activity_small[crsID %in% cre_A549_up, ], by =  c("oligo",'crsID'), allow.cartesian = TRUE)

### EXTRACT STATS FOR CELL TYPE SPECIFIC OLIGOS (DOWNREGULATED)

oligo_K562_specific_down <- oligo_celltype_comp_NS_annot_obs_archaic[crsID %in% cre_K562_lo, .(oligo, crsID, allele.label, group1_labels, group2_labels, FDR, pval_emp, log2FC_2vs1, log2FC.se, pval.LRT, oligo_diff_CRITERIA)]
oligo_HepG2_specific_down <- oligo_celltype_comp_NS_annot_obs_archaic[crsID %in% cre_HepG2_lo, .(oligo, crsID, allele.label, group1_labels, group2_labels, FDR, pval_emp, log2FC_2vs1, log2FC.se, pval.LRT, oligo_diff_CRITERIA)]
oligo_A549_specific_down <- oligo_celltype_comp_NS_annot_obs_archaic[crsID %in% cre_A549_lo, .(oligo, crsID, allele.label, group1_labels, group2_labels, FDR, pval_emp, log2FC_2vs1, log2FC.se, pval.LRT, oligo_diff_CRITERIA)]

oligo_K562_specific_down <- merge(oligo_K562_specific_down, oligo_activity_small[crsID %in% cre_K562_lo, ], by = c("oligo",'crsID'), allow.cartesian = TRUE)
oligo_HepG2_specific_down <- merge(oligo_HepG2_specific_down, oligo_activity_small[crsID %in% cre_HepG2_lo, ], by =  c("oligo",'crsID'), allow.cartesian = TRUE)
oligo_A549_specific_down <- merge(oligo_A549_specific_down, oligo_activity_small[crsID %in% cre_A549_lo, ], by =  c("oligo",'crsID'), allow.cartesian = TRUE)

##########################################  A549-specific ##########################################
###### reshape data to wide format (A549-specific upregulated)
oligo_A549_specific_up[, log2FC := ifelse(group2_labels == "A549", log2FC_2vs1, -log2FC_2vs1)]
oligo_A549_specific_up[, comparison_label := ifelse(group2_labels == "A549", group1_labels, group2_labels)]
oligo_A549_specific_up <- dcast(oligo_A549_specific_up[group2_labels == "A549" | group1_labels == "A549", ], oligo + crsID + allele.label + oligo_class_CRITERIA_A549 + alpha_CRITERIA_A549 + alpha_CRITERIA_HepG2 + oligo_class_CRITERIA_HepG2 + alpha_CRITERIA_K562 + oligo_class_CRITERIA_K562~ comparison_label, value.var = c("log2FC", "log2FC.se", "FDR", "pval_emp"))
oligo_A549_specific_up[, target_cellline := "A549"]
oligo_A549_specific_up <- oligo_A549_specific_up[, .(`CRE ID`=crsID, allele.label, `unique CRS ID`=oligo, 
  target_cellline,
  alpha_A549 = alpha_CRITERIA_A549,
  crs_activity_class_target = oligo_class_CRITERIA_A549,
  specificity = "higher transcription in A549",
  # alpha_K562 = alpha_CRITERIA_K562,
  #  crs_activity_class_K562 = oligo_class_CRITERIA_K562,
  log2FC_K562, log2FC.se_K562, pval_emp_K562, FDR_K562,
  # alpha_HepG2 = alpha_CRITERIA_HepG2,
	# crs_activity_class_HepG2 = oligo_class_CRITERIA_HepG2,
  log2FC_HepG2, log2FC.se_HepG2, pval_emp_HepG2, FDR_HepG2,
	is_A549_specific=ifelse(oligo%in%c(oligo_A549_up),'x','')
)]
fwrite(oligo_A549_specific_up, file = sprintf("%s/celltype_specificity_oligo_A549_specific_up.txt", FIGURE_DIR), sep = "\t")

###### reshape data to wide format (A549-specific downregulated)
oligo_A549_specific_down[, log2FC := ifelse(group2_labels == "A549", log2FC_2vs1, -log2FC_2vs1)]
oligo_A549_specific_down[, comparison_label := ifelse(group2_labels == "A549", group1_labels, group2_labels)]
oligo_A549_specific_down <- dcast(oligo_A549_specific_down[group2_labels == "A549" | group1_labels == "A549", ], oligo + crsID + allele.label +  oligo_class_CRITERIA_A549 + alpha_CRITERIA_A549 + alpha_CRITERIA_HepG2 + oligo_class_CRITERIA_HepG2 + alpha_CRITERIA_K562 + oligo_class_CRITERIA_K562 ~ comparison_label, value.var = c("log2FC", "log2FC.se", "FDR", "pval_emp"))
oligo_A549_specific_down[, target_cellline := "A549"]
oligo_A549_specific_down <- oligo_A549_specific_down[, .(`CRE ID`=crsID, allele.label, `unique CRS ID`=oligo, 
  target_cellline,
  alpha_A549 = alpha_CRITERIA_A549,
  crs_activity_class_target = oligo_class_CRITERIA_A549,
  specificity = "lower transcription in A549",
  # alpha_K562 = alpha_CRITERIA_K562,
	# crs_activity_class_K562 = oligo_class_CRITERIA_K562,
  log2FC_K562, log2FC.se_K562, pval_emp_K562, FDR_K562,
  # alpha_HepG2 = alpha_CRITERIA_HepG2,
	# crs_activity_class_HepG2 = oligo_class_CRITERIA_HepG2,
  log2FC_HepG2, log2FC.se_HepG2, pval_emp_HepG2, FDR_HepG2,
	is_A549_specific=ifelse(oligo%in%c(oligo_A549_lo),'x','')
)]
fwrite(oligo_A549_specific_down, file = sprintf("%s/celltype_specificity_oligo_A549_specific_down.txt", FIGURE_DIR), sep = "\t")

###### generate SUpTable 2A ######
SupTable2a_oligo_A549_specific_all <- rbind(oligo_A549_specific_up, oligo_A549_specific_down)
SupTable2a_oligo_A549_specific_all[,score:=max(sign(log2FC_HepG2+log2FC_K562)*abs(log2FC_HepG2+log2FC_K562)),by=`CRE ID`]
SupTable2a_oligo_A549_specific_all <- SupTable2a_oligo_A549_specific_all[order(-score,allele.label),-'score']
fwrite(SupTable2a_oligo_A549_specific_all, file = sprintf("%s/SupTables/SupTable2a_oligo_A549_specific_all.txt", FIGURE_DIR), sep = "\t") #  fwrite SupTable S2a-c


##########################################  HEPG2-specific ##########################################
###### reshape data to wide format (HepG2-specific upregulated)
oligo_HepG2_specific_up[, log2FC := ifelse(group2_labels == "HepG2", log2FC_2vs1, -log2FC_2vs1)]
oligo_HepG2_specific_up[, comparison_label := ifelse(group2_labels == "HepG2", group1_labels, group2_labels)]
oligo_HepG2_specific_up <- dcast(oligo_HepG2_specific_up[group2_labels == "HepG2" | group1_labels == "HepG2", ], oligo + crsID + allele.label +  oligo_class_CRITERIA_A549 + alpha_CRITERIA_A549 + alpha_CRITERIA_HepG2 + oligo_class_CRITERIA_HepG2 + alpha_CRITERIA_K562 + oligo_class_CRITERIA_K562 ~ comparison_label, value.var = c("log2FC", "log2FC.se", "FDR", "pval_emp"))
oligo_HepG2_specific_up[, target_cellline := "HepG2"]
oligo_HepG2_specific_up <- oligo_HepG2_specific_up[, .(`CRE ID`=crsID, allele.label, `unique CRS ID`=oligo, 
  target_cellline,
  alpha_HepG2 = alpha_CRITERIA_HepG2,
  crs_activity_class_target = oligo_class_CRITERIA_HepG2,
  specificity = "higher transcription in HepG2",
  # alpha_K562 = alpha_CRITERIA_K562,
	# crs_activity_class_K562 = oligo_class_CRITERIA_K562,
  log2FC_K562, log2FC.se_K562, pval_emp_K562, FDR_K562,
  # alpha_A549 = alpha_CRITERIA_A549,
	# crs_activity_class_A549 = oligo_class_CRITERIA_A549,
  log2FC_A549, log2FC.se_A549,  pval_emp_A549, FDR_A549,
	is_HepG2_specific=ifelse(oligo%in%c(oligo_HepG2_up),'x','')
)]
fwrite(oligo_HepG2_specific_up, file = sprintf("%s/celltype_specificity_oligo_HepG2_specific_up.txt", FIGURE_DIR), sep = "\t")

###### reshape data to wide format (HepG2-specific downregulated)
oligo_HepG2_specific_down[, log2FC := ifelse(group2_labels == "HepG2", log2FC_2vs1, -log2FC_2vs1)]
oligo_HepG2_specific_down[, comparison_label := ifelse(group2_labels == "HepG2", group1_labels, group2_labels)]
oligo_HepG2_specific_down <- dcast(oligo_HepG2_specific_down[group2_labels == "HepG2" | group1_labels == "HepG2", ], oligo + crsID + allele.label +  oligo_class_CRITERIA_A549 + alpha_CRITERIA_A549 + alpha_CRITERIA_HepG2 + oligo_class_CRITERIA_HepG2 + alpha_CRITERIA_K562 + oligo_class_CRITERIA_K562 ~ comparison_label, value.var = c("log2FC", "log2FC.se", "FDR", "pval_emp"))
oligo_HepG2_specific_down[, target_cellline := "HepG2"]
oligo_HepG2_specific_down <- oligo_HepG2_specific_down[, .(`CRE ID`=crsID, allele.label, `unique CRS ID`=oligo, 
  target_cellline,
  alpha_HepG2 = alpha_CRITERIA_HepG2,
  crs_activity_class_target = oligo_class_CRITERIA_HepG2,
  specificity = "lower transcription in HepG2",
  # alpha_K562 = alpha_CRITERIA_K562,
	# crs_activity_class_K562 = oligo_class_CRITERIA_K562,
  log2FC_K562, log2FC.se_K562, pval_emp_K562,FDR_K562, 
  # alpha_A549 = alpha_CRITERIA_A549,
	# crs_activity_class_A549 = oligo_class_CRITERIA_A549,
  log2FC_A549, log2FC.se_A549, pval_emp_A549, FDR_A549,
	is_HepG2_specific=ifelse(oligo%in%c(oligo_HepG2_lo),'x','')
)]
fwrite(oligo_HepG2_specific_down, file = sprintf("%s/celltype_specificity_oligo_HepG2_specific_down.txt", FIGURE_DIR), sep = "\t")

###### generate SUpTable 2B ######
SupTable2b_oligo_HepG2_specific_all <- rbind(oligo_HepG2_specific_up, oligo_HepG2_specific_down)

SupTable2b_oligo_HepG2_specific_all[,score:=max(sign(log2FC_A549+log2FC_K562)*abs(log2FC_A549+log2FC_K562)),by=`CRE ID`]
SupTable2b_oligo_HepG2_specific_all <- SupTable2b_oligo_HepG2_specific_all[order(-score,allele.label),-'score']
fwrite(SupTable2b_oligo_HepG2_specific_all, file = sprintf("%s/SupTables/SupTable2b_oligo_HepG2_specific_all.txt", FIGURE_DIR), sep = "\t") #  fwrite SupTable S2a-c

##########################################  K562-specific ##########################################
###### reshape data to wide format (K562-specific upregulated)
oligo_K562_specific_up[, log2FC := ifelse(group2_labels == "K562", log2FC_2vs1, -log2FC_2vs1)]
oligo_K562_specific_up[, comparison_label := ifelse(group2_labels == "K562", group1_labels, group2_labels)]
oligo_K562_specific_up <- dcast(oligo_K562_specific_up[group2_labels == "K562" | group1_labels == "K562", ], oligo + crsID + allele.label + oligo_class_CRITERIA_A549 + alpha_CRITERIA_A549 + alpha_CRITERIA_HepG2 + oligo_class_CRITERIA_HepG2 + alpha_CRITERIA_K562 + oligo_class_CRITERIA_K562 ~ comparison_label, value.var = c("log2FC", "log2FC.se", "FDR", "pval_emp"))
oligo_K562_specific_up[, target_cellline := "K562"]
oligo_K562_specific_up <- oligo_K562_specific_up[, .(`CRE ID`=crsID, allele.label, `unique CRS ID`=oligo, 
  target_cellline,
  alpha_K562 = alpha_CRITERIA_K562,
  crs_activity_class_target = oligo_class_CRITERIA_K562,
  specificity = "higher transcription in K562",
  # alpha_A549 = alpha_CRITERIA_A549,
	# crs_activity_class_A549 = oligo_class_CRITERIA_A549,
  log2FC_A549, log2FC.se_A549, pval_emp_A549, FDR_A549, 
  # alpha_HepG2 = alpha_CRITERIA_HepG2,
	# crs_activity_class_HepG2 = oligo_class_CRITERIA_HepG2,
  log2FC_HepG2, log2FC.se_HepG2,  pval_emp_HepG2, FDR_HepG2,
	is_K562_specific=ifelse(oligo%in%c(oligo_K562_up),'x','')
)]
fwrite(oligo_K562_specific_up, file = sprintf("%s/celltype_specificity_oligo_K562_specific_up.txt", FIGURE_DIR), sep = "\t")

###### reshape data to wide format (K562-specific downregulated)
oligo_K562_specific_down[, log2FC := ifelse(group2_labels == "K562", log2FC_2vs1, -log2FC_2vs1)]
oligo_K562_specific_down[, comparison_label := ifelse(group2_labels == "K562", group1_labels, group2_labels)]
oligo_K562_specific_down <- dcast(oligo_K562_specific_down[group2_labels == "K562" | group1_labels == "K562", ], oligo + crsID + allele.label +  oligo_class_CRITERIA_A549 + alpha_CRITERIA_A549 + alpha_CRITERIA_HepG2 + oligo_class_CRITERIA_HepG2 + alpha_CRITERIA_K562 + oligo_class_CRITERIA_K562 ~ comparison_label, value.var = c("log2FC", "log2FC.se", "FDR", "pval_emp"))
oligo_K562_specific_down[, target_cellline := "K562"]
oligo_K562_specific_down <- oligo_K562_specific_down[, .(`CRE ID`=crsID, allele.label, `unique CRS ID`=oligo, 
  target_cellline,
  alpha_K562 = alpha_CRITERIA_K562,
  crs_activity_class_target = oligo_class_CRITERIA_K562,
  specificity = "lower transcription in K562",
  # alpha_A549 = alpha_CRITERIA_A549,
	# crs_activity_class_A549 = oligo_class_CRITERIA_A549,
  log2FC_A549, log2FC.se_A549, pval_emp_A549, FDR_A549,
  # alpha_HepG2 = alpha_CRITERIA_HepG2,
	# crs_activity_class_HepG2 = oligo_class_CRITERIA_HepG2,
  log2FC_HepG2, log2FC.se_HepG2, pval_emp_HepG2, FDR_HepG2,
	is_K562_specific=ifelse(oligo%in%c(oligo_K562_lo),'x','')
)]
fwrite(oligo_K562_specific_down, file = sprintf("%s/celltype_specificity_oligo_K562_specific_down.txt", FIGURE_DIR), sep = "\t")

###### generate SUpTable 2C ######
SupTable2c_oligo_K562_specific_all <- rbind(oligo_K562_specific_up, oligo_K562_specific_down)
SupTable2c_oligo_K562_specific_all[,score:=max(sign(log2FC_A549+log2FC_HepG2)*abs(log2FC_A549+log2FC_HepG2)),by=`CRE ID`]
SupTable2c_oligo_K562_specific_all <- SupTable2c_oligo_K562_specific_all[order(-score,allele.label),-'score']
fwrite(SupTable2c_oligo_K562_specific_all, file = sprintf("%s/SupTables/SupTable2c_oligo_K562_specific_all.txt", FIGURE_DIR), sep = "\t") #  fwrite SupTable S2a-c

####################################################################################################################
####################################################################################################################
####################################################################################################################

celltype_specificity_count <- rbind(
  data.table(type = "oligo up", K562 = length(oligo_K562_up), HepG2 = length(oligo_HepG2_up), A549 = length(oligo_A549_up)),
  data.table(type = "oligo down", K562 = length(oligo_K562_lo), HepG2 = length(oligo_HepG2_lo), A549 = length(oligo_A549_lo)),
  data.table(
    type = "cre up", 
		K562 = length(SupTable2c_oligo_K562_specific_all[grep("higher", specificity), unique(`CRE ID`)]),
    HepG2 = length(SupTable2b_oligo_HepG2_specific_all[grep("higher", specificity), unique(`CRE ID`)]),
    A549 = length(SupTable2a_oligo_A549_specific_all[grep("higher", specificity), unique(`CRE ID`)])
  ),
  data.table(
    type = "cre down", 
		K562 = length(SupTable2c_oligo_K562_specific_all[grep("low", specificity), unique(`CRE ID`)]),
    HepG2 = length(SupTable2b_oligo_HepG2_specific_all[grep("low", specificity), unique(`CRE ID`)]),
    A549 = length(SupTable2a_oligo_A549_specific_all[grep("low", specificity), unique(`CRE ID`)])
  )
)

celltype_specificity_count
#          type  K562 HepG2  A549
#        <char> <int> <int> <int>
# 1:   oligo up   816   343   299
# 2: oligo down   739   140   160
# 3:     cre up   505   220   236
# 4:   cre down   447   114   109
fwrite(celltype_specificity_count, file = sprintf("%s/celltype_specificity_oligo_counts.txt", FIGURE_DIR))


#################################################################
####### pairwise comparison between cell types  #################
#################################################################

# TODO: add THESE plots


cat('All done!')
q('no')