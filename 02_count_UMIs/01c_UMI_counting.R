MPRA_DIR <- "/pasteur/helix/projects/evo_immuno_pop/MPRA"
ANALYSIS_DIR <- "MPRA_count_exp6_analysis3"

source(sprintf("%s/scripts/%s/00_starter.R", MPRA_DIR, ANALYSIS_DIR))
source(sprintf("%s/scripts/%s/00_load_associations.R", MPRA_DIR, ANALYSIS_DIR))
source(sprintf("%s/scripts/misc_plots.R", MPRA_DIR))

cmd <- commandArgs()
print(cmd)

VERSION_IN <- "1"
VERSION_OUT <- "1"
SUB_SAMPLE <- -1
SUB_SAMPLE_Pct <- FALSE
for (i in seq_along(cmd)) {
  if (cmd[i] == "--sample" | cmd[i] == "-s") {
    j <- as.numeric(cmd[i + 1])
  } # ID of the set to test
  if (cmd[i] == "--version_in" | cmd[i] == "-vi") {
    VERSION_IN <- cmd[i + 1]
  } # version of input files
  if (cmd[i] == "--version_out" | cmd[i] == "-vo") {
    VERSION_OUT <- cmd[i + 1]
  } # version of output files
  if (cmd[i] == "--sub_sample" | cmd[i] == "-u") {
    if (grepl("%", cmd[i + 1])) {
      SUB_SAMPLE_Pct <- TRUE
      SUB_SAMPLE <- as.numeric(gsub("%", "", cmd[i + 1])) / 100
    } else {
      SUB_SAMPLE <- as.numeric(cmd[i + 1])
    }
  } # subsampling parameters
}

mySETUP_ID <- Indexes[, unique(SETUP_ID)[j]]
Indexes[SETUP_ID == mySETUP_ID, ]

# load oligo annotations
oligo_source <- fread(sprintf("%s/data/%s/00_oligo_annot_v1.txt", MPRA_DIR, ANALYSIS_DIR))
oligo_type <- oligo_source[, .(type = paste(unique(type), collapse = "\\")), by = oligo]
# list oligos associated with >1 source
dup_oligo <- oligo_source[duplicated(oligo), unique(oligo)]


if (SUB_SAMPLE > 0) {
  # define subsampling directory
  if (SUB_SAMPLE_Pct == TRUE) {
    ANALYSIS_DIR <- sprintf("%s/01d_subsampling/SUB_%sPct", ANALYSIS_DIR, 100 * SUB_SAMPLE)
  } else {
    ANALYSIS_DIR <- sprintf("%s/01d_subsampling/SUB_%s", ANALYSIS_DIR, SUB_SAMPLE)
  }
}
# Load UMI counts per barcode
#  - for each sample (ID_rep) with material (RNA or DNA)
UMI_perBC_annot <- list()
for (MATERIAL in c("RNA", "DNA")) {
  # get sample ID
  SAMPLE_ID <- Indexes[SETUP_ID == mySETUP_ID & material == MATERIAL, unique(Sample_ID)]
  SID <- Indexes[SETUP_ID == mySETUP_ID & material == MATERIAL, unique(SID)]
  # read UMI counts
  UMI_perBC_annot[[MATERIAL]] <- fread(sprintf("%s/data/%s/01b_MPRA_results_counts_BC/MPRA_results_counts_BC__%s__v%s.txt.gz", MPRA_DIR, ANALYSIS_DIR, SAMPLE_ID, VERSION_IN))
}
# merge data from RNA and DNA samples & add library info
UMI_perBC_annot <- rbindlist(UMI_perBC_annot, idcol = "material")
UMI_perBC_annot <- merge(UMI_perBC_annot[BC_type == "1a_assoc_unique", ], LIB_summary, by = c("library", "material"))

# add BC_class (low complexity, known, unknown)
UMI_perBC_annot[, BC_class := ifelse(shannon < 10, "low complexity", ifelse(known_BC1, "known", "unknown"))]

# cast to get nUMI_perBC for each oligo, barcode1, known_BC1, shannon, celltype, condition, replicate, material (RNA or DNA)
UMI_perBC_annot_DNA_vs_RNA <- dcast(UMI_perBC_annot, barcode1 + known_BC1 + shannon + BC_type + BC_class + celltype + condition + replicate + SETUP_ID + barcode_lib + celline + COND_ID + experiment + oligo ~ material, value.var = list("nUMI_perBC"), fill = 0, fun.agg = sum)

#  assign outlier Z scores & to filter out barcodes with extreme activity.
UMI_perBC_annot_DNA_vs_RNA[RNA > 0 & DNA > 0, Ratio := RNA / DNA]
UMI_perBC_annot_DNA_vs_RNA[RNA > 0 & DNA > 0, outlier_Z := (Ratio - mean(Winsorize(Ratio))) / sd(Winsorize(Ratio)), by = .(oligo, celltype, condition, replicate)]
# UMI_perBC_annot_DNA_vs_RNA[,mean(abs(outlier_Z)>3,na.rm=T),keyby=.(cut(DNA,c(0:10,20,30,50,Inf),include.lowest=T))]


# save to file
dir.create(sprintf("%s/data/%s/01c_MPRA_results_BC_activity", MPRA_DIR, ANALYSIS_DIR))

fwrite(UMI_perBC_annot_DNA_vs_RNA, file = sprintf("%s/data/%s/01c_MPRA_results_BC_activity/01c_MPRA_results_BC_activity__%s__v%s.txt.gz", MPRA_DIR, ANALYSIS_DIR, mySETUP_ID, VERSION_OUT), sep = "\t")
# UMI_perBC_annot_DNA_vs_RNA <- fread(file = sprintf("%s/data/%s/01c_MPRA_results_BC_activity/01c_MPRA_results_BC_activity__%s__v%s.txt.gz", MPRA_DIR, ANALYSIS_DIR, mySETUP_ID, VERSION_OUT))

##### TODO: reassess if the following part is useful in the end.
# TODO: remove if not needed

### calculate Enrichment of RNA vs DNA, for each oligo
# -------------------------------------------------------------------------------
# UMI counts per oligo, for each library (ID_rep), material (RNA or DNA)
UMI_per_oligo <- UMI_perBC_annot[BC_type == "1a_assoc_unique" & BC_class == "known",
  .(nUMI = sum(nUMI_perBC), nBC = length(unique(barcode1))),
  by = .(oligo, library, celltype, material, condition, replicate)
]
# merge counts per oligo and material (RNA or DNA)
UMI_per_oligo_DNA_vs_RNA <- dcast(UMI_per_oligo, oligo + celltype + condition + replicate ~ material, value.var = list("nUMI", "nBC"))
# calculate total number of UMIs per library per material
UMI_per_oligo_DNA_vs_RNA[!is.na(oligo), nUMI_DNA_tot := sum(nUMI_DNA, na.rm = TRUE), by = .(celltype, condition, replicate)]
UMI_per_oligo_DNA_vs_RNA[!is.na(oligo), nUMI_RNA_tot := sum(nUMI_RNA, na.rm = TRUE), by = .(celltype, condition, replicate)]
# filter out oligos with missing data
UMI_per_oligo_DNA_vs_RNA <- UMI_per_oligo_DNA_vs_RNA[!is.na(oligo) & !is.na(nUMI_RNA) & !is.na(nUMI_DNA)]
# perform Fisher's exact test: RNA vs DNA expression
UMI_per_oligo_DNA_vs_RNA <- UMI_per_oligo_DNA_vs_RNA[,
  fisher.test(matrix(c(nUMI_DNA_tot, nUMI_RNA_tot, nUMI_DNA, nUMI_RNA), 2))[c("estimate", "p.value")],
  by = .(
    oligo, celltype, condition, replicate,
    nUMI_DNA_tot, nUMI_RNA_tot, nUMI_DNA, nUMI_RNA,
    nBC_DNA, nBC_RNA
  )
]
# p-values correction: false discovery rate (FDR)
UMI_per_oligo_DNA_vs_RNA[, FDR := p.adjust(p.value, "fdr")]
# compute fold-enrichment (FE manually)
UMI_per_oligo_DNA_vs_RNA[, FE := (1 + nUMI_RNA) / (1 + nUMI_DNA) / (median((1 + nUMI_RNA) / (1 + nUMI_DNA)))]

# TODO: make sure there is no duplicates in the oligo_type (one interpretation per oligo)
# annotate oligo with type and GC content
UMI_per_oligo_DNA_vs_RNA <- merge(UMI_per_oligo_DNA_vs_RNA, oligo_type, by = "oligo", all.x = TRUE)
UMI_per_oligo_DNA_vs_RNA[, type_simple := case_when(
  grepl("Promoter", type) ~ "Promoter",
  grepl("scrambled", type) ~ "scrambled",
  TRUE ~ "other"
)]
UMI_per_oligo_DNA_vs_RNA <- merge(UMI_per_oligo_DNA_vs_RNA, oligo_source[, .(oligo, GC, crsID)], all.x = TRUE)

# compute adjusted fold-enrichment
setnames(UMI_per_oligo_DNA_vs_RNA, "estimate", "OR")
UMI_per_oligo_DNA_vs_RNA[, OR_GCadj := 2^(log2(OR) - lm(log2(OR) ~ GC)$coef[2] * GC - lm(log2(OR) ~ GC)$coef[1])]

# compute GC-content corrected activity P-values
UMI_per_oligo_DNA_vs_RNA[, P_scrambled := 2 * pnorm(abs(log2(OR_GCadj)),
                        median(log2(OR_GCadj[type == "scrambled"]), na.rm = TRUE),
                        IQR(log2(OR_GCadj[type == "scrambled"]), na.rm = TRUE) / 1.35,
                        low = FALSE
                      ), by = celltype]

# compute FDR on GC-corrected activity P-values
UMI_per_oligo_DNA_vs_RNA[, FDR_scrambled := p.adjust(P_scrambled, "fdr"), by = celltype]

# write output
dir.create(sprintf("%s/data/%s/01c_MPRA_results_oligo_activity", MPRA_DIR, ANALYSIS_DIR), recursive = TRUE)
fwrite(UMI_per_oligo_DNA_vs_RNA, file = sprintf("%s/data/%s/01c_MPRA_results_oligo_activity/01c_MPRA_results_oligo_activity__%s__v%s.txt.gz", MPRA_DIR, ANALYSIS_DIR, mySETUP_ID, VERSION_OUT), sep = "\t")


cat("\nAll done\n")
