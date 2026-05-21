####### DONE : implement subsmapling plots


MPRA_DIR <- "/pasteur/helix/projects/evo_immuno_pop/MPRA"
ANALYSIS_DIR <- "MPRA_count_exp6_analysis3"
source(sprintf("%s/scripts/misc_plots.R", MPRA_DIR))
source(sprintf("%s/scripts/%s/00_starter.R", MPRA_DIR, ANALYSIS_DIR))
source(sprintf("%s/scripts/%s/00_load_associations.R", MPRA_DIR, ANALYSIS_DIR))



FIGURE_DIR <- sprintf("%s/figures/%s", MPRA_DIR, ANALYSIS_DIR)
dir.create(FIGURE_DIR, showWarnings = FALSE, recursive = TRUE)

################################################################################
########    QC plots                                                     #######
################################################################################

VERSION_IN <- "1"
VERSION_OUT <- "1"

# load oligo annotations
oligo_source <- fread(sprintf("%s/data/%s/00_oligo_annot_v1.txt", MPRA_DIR, ANALYSIS_DIR))
#oligo_source[, allele.num2 := paste0("A", 1:.N, sep = ""), by = .(posID, shift, strand)]
oligo_type <- oligo_source[, .(type = paste(unique(type), collapse = "\\")), by = oligo]
# list oligos associated with >1 source
dup_oligo <- oligo_source[duplicated(oligo), unique(oligo)]

# load available subsamplings
SUB_DIRs <- list.dirs(sprintf("%s/data/%s/01d_subsampling", MPRA_DIR, ANALYSIS_DIR), full.names = FALSE, recursive = FALSE)
SUB_DIRs <- c("SUB_100Pct", SUB_DIRs)
###############################################################################
###########  reading/formatting read/UMI counts             ###################
###############################################################################

FAILED_SAMPLES <- c()

SUB_Counts_by_lib <- list()
SUB_UMI_per_BC <- list()
SUB_BC_recovery <- list()
for (i_SUB in seq_along(SUB_DIRs)) {
  SUB_DIR <- SUB_DIRs[i_SUB]
  cat("\n", SUB_DIR, ":", LIB_summary[, .N], "samples.\n")
  Counts_by_lib <- list()
  UMI_per_BC <- list()
  BC_recovery <- list()
  for (i in seq_len(LIB_summary[, .N])) {
    SID <- LIB_summary[i, SID]
    SAMPLE_ID <- LIB_summary[i, library]
    cat(i, "")
    # load UMI per BC
    if (SUB_DIR == "SUB_100Pct") {
      Counts_by_lib[[SAMPLE_ID]] <- fread(sprintf("%s/data/%s/01b_Count_reads_UMI_type_frequency/Count_reads_UMI_type_frequency__%s__v%s.tsv", MPRA_DIR, ANALYSIS_DIR, SAMPLE_ID, VERSION_IN))


      UMI_per_BC[[SAMPLE_ID]] <- fread(sprintf("%s/data/%s/01b_MPRA_results_counts_BC/MPRA_results_counts_BC__%s__v%s.txt.gz", MPRA_DIR, ANALYSIS_DIR, SAMPLE_ID, VERSION_IN))
      BC_recovery[[SAMPLE_ID]] <- fread(sprintf("%s/data/%s/01b_BC_recovery/BC_recovery__%s__v%s.txt", MPRA_DIR, ANALYSIS_DIR, SAMPLE_ID, VERSION_IN))
    } else {
      Counts_by_lib[[SAMPLE_ID]] <- fread(sprintf("%s/data/%s/01d_subsampling/%s/01b_Count_reads_UMI_type_frequency/Count_reads_UMI_type_frequency__%s__v%s.tsv", MPRA_DIR, ANALYSIS_DIR, SUB_DIR, SAMPLE_ID, VERSION_IN))
      UMI_per_BC[[SAMPLE_ID]] <- fread(sprintf("%s/data/%s/01d_subsampling/%s/01b_MPRA_results_counts_BC/MPRA_results_counts_BC__%s__v%s.txt.gz", MPRA_DIR, ANALYSIS_DIR, SUB_DIR, SAMPLE_ID, VERSION_IN))
      BC_recovery[[SAMPLE_ID]] <- try(fread(sprintf("%s/data/%s/01d_subsampling/%s/01b_BC_recovery/BC_recovery__%s__v%s.txt", MPRA_DIR, ANALYSIS_DIR, SUB_DIR, SAMPLE_ID, VERSION_IN)))
      if (any(class(BC_recovery[[SAMPLE_ID]]) == "try-error")) {
        FAILED_SAMPLES <- c(FAILED_SAMPLES, paste(gsub("SUB_(.*)Pct", "\\1%", SUB_DIR), Indexes[, which(unique(Sample_ID) == SAMPLE_ID)], sep = "__"))
        BC_recovery[[SAMPLE_ID]] <- NULL
      }
    }
  }
  Counts_by_lib <- rbindlist(Counts_by_lib, idcol = "library")
  Counts_by_lib <- merge(Counts_by_lib, LIB_summary, by = "library")
  UMI_per_BC <- rbindlist(UMI_per_BC, idcol = "library")
  BC_recovery <- rbindlist(BC_recovery, idcol = "library")
  SUB_Counts_by_lib[[SUB_DIR]] <- Counts_by_lib
  SUB_UMI_per_BC[[SUB_DIR]] <- UMI_per_BC
  SUB_BC_recovery[[SUB_DIR]] <- BC_recovery
}
cat("\n'")
cat(FAILED_SAMPLES, sep = "' '")
cat("'\n\n")


SUB_Counts_by_lib <- rbindlist(SUB_Counts_by_lib, idcol = "SUB")
SUB_Counts_by_lib[, sub_pct := as.numeric(gsub("SUB_(.*)Pct", "\\1", SUB))]
SUB_UMI_per_BC <- rbindlist(SUB_UMI_per_BC, idcol = "SUB")
SUB_UMI_per_BC[, sub_pct := as.numeric(gsub("SUB_(.*)Pct", "\\1", SUB))]
SUB_BC_recovery <- rbindlist(SUB_BC_recovery, idcol = "SUB")
SUB_BC_recovery[, sub_pct := as.numeric(gsub("SUB_(.*)Pct", "\\1", SUB))]

stats_sub <- SUB_UMI_per_BC[BC_type == "1a_assoc_unique", .(
  nBC = .N,
  nUMI = sum(nUMI_perBC),
  med_UMI_perBC = as.numeric(median(nUMI_perBC)),
  rank_mean = mean(quantile(nUMI_perBC, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))),
  winsor_mean = mean(Winsorize(nUMI_perBC))
), keyby = .(sub_pct, library)]

read_count <- SUB_Counts_by_lib[, .(million_read = sum(N) / 1e6, pass_reads = sum(N[known_BC1 & sameBC & known_BC2])), keyby = .(library, sub_pct)]
read_count[, sequencing_efficiency := 100 * (pass_reads / (million_read * 1e6))]

stats_sub <- merge(stats_sub, read_count, by = c("library", "sub_pct"))

stats_sub <- merge(stats_sub, SUB_BC_recovery, by = c("library", "sub_pct"))
stats_sub <- merge(stats_sub, LIB_summary, by = c("library"))
stats_sub <- merge(stats_sub, Associations_Filtered[, .(expected_barcodes = .N), by = .(barcode_lib = barcode_library)], by = c("barcode_lib"))
#stats_sub <- merge(stats_sub, unique(Indexes[, .(experiment, SETUP_ID)]), by = "SETUP_ID")
stats_sub[, sequencing_saturation := 100 * (1 - nUMI / (pass_reads))]
stats_sub <- stats_sub[order(library, sub_pct)]
stats_sub[, sequencing_saturation_new := 100 * (1 - diff(c(0, nUMI)) / diff(c(0, pass_reads))), by = .(library)]

stats_sub[, AdditionalUMI_for20pctMoreReads := (100 - sequencing_saturation_new) / 100 * sequencing_efficiency / 100 * 0.2 * million_read * 1e6]
stats_sub[, PctAdditionalUMI_for20pctMoreReads := AdditionalUMI_for20pctMoreReads / nUMI]
stats_sub[, NewUMIper100read := (100 - sequencing_saturation_new) / 100 * sequencing_efficiency / 100 * 100]

stats_sub[, group := paste(experiment, barcode_lib, sep = "-")]
stats_sub_bckgd <- stats_sub[, mget(setdiff(colnames(stats_sub), c("group")))]



dir.create(sprintf("%s/01d_subsampling/01_per_million_read", FIGURE_DIR), showWarnings = FALSE, recursive = TRUE)
dir.create(sprintf("%s/01d_subsampling/02_percentage_total_read", FIGURE_DIR), showWarnings = FALSE, recursive = TRUE)

decision_making <- stats_sub[sub_pct == 100, ][order(material, experiment, barcode_lib, million_read), .(material, experiment, barcode_lib, library, million_read, sequencing_saturation_new, NewUMIper100read, PctAdditionalUMI_for20pctMoreReads, Pct_BC_over3UMI, nUMI_perBC_winsorizedmean = winsor_mean)]
fwrite(decision_making, sprintf("%s/01d_subsampling/stats_subsampling_for_decision_making.txt", FIGURE_DIR), sep = "\t")

#############################################
######         per million read        ######
#############################################

addScales <- function(p, y_pct = FALSE, x_pct = FALSE) {
  p <- p + scale_color_manual(values = color_setup, name = "sample")
  if (y_pct) {
    p <- p + scale_y_continuous(labels = scales::percent, limit = c(0, 1))
  } else {
    p <- p + scale_y_continuous(trans = "log10")
  }
  if (x_pct) {
    p <- p + scale_x_continuous(labels = scales::percent, trans = "log10")
  } else {
    p <- p + scale_x_continuous(trans = "log10")
  }
  p
}

### total number of UMI
# raw
p <- ggplot(stats_sub, aes(x = million_read, y = nUMI, col = SETUP_ID))
p <- p + geom_line(data = stats_sub_bckgd, aes(group = paste(experiment, barcode_lib, SETUP_ID)), col = "lightgrey", alpha = .2) + geom_line(data = stats_sub) + geom_point()
p <- p + facet_grid(cols = vars(group), rows = vars(material))
p <- p + theme_plot() + ylab("Total Number of UMIs")
p <- addScales(p)
pdf(sprintf("%s/01d_subsampling/01_per_million_read/01a_saturation_nUMI_vs_million_read.pdf", FIGURE_DIR), height = 6, width = 7.5)
print(p)
dev.off()

# per 1M expected barcode
p <- ggplot(stats_sub, aes(x = million_read, y = nUMI / expected_barcodes * 1e6, col = SETUP_ID))
p <- p + geom_line(data = stats_sub_bckgd, aes(group = paste(experiment, barcode_lib, SETUP_ID)), col = "lightgrey", alpha = .2) + geom_line(data = stats_sub) + geom_point()
p <- p + facet_grid(cols = vars(group), rows = vars(material))
p <- p + theme_plot() + ylab("Total Number of UMIs per 1M expected barcodes")
p <- addScales(p)
pdf(sprintf("%s/01d_subsampling/01_per_million_read/01b_saturation_nUMI_per1MBC_vs_million_read.pdf", FIGURE_DIR), height = 6, width = 7.5)
print(p)
dev.off()

# raw - linear
p <- ggplot(stats_sub, aes(x = million_read, y = nUMI, col = SETUP_ID))
p <- p + geom_line(data = stats_sub_bckgd, aes(group = paste(experiment, barcode_lib, SETUP_ID)), col = "lightgrey", alpha = .2) + geom_line(data = stats_sub) + geom_point()
p <- p + facet_grid(cols = vars(group), rows = vars(material))
p <- p + theme_plot() + ylab("Total Number of UMIs")
p <- p + scale_x_continuous(trans = "log10") + scale_y_continuous() + scale_color_manual(values = color_setup, name = "sample")
pdf(sprintf("%s/01d_subsampling/01_per_million_read/01c_saturation_nUMI_linear_vs_million_read.pdf", FIGURE_DIR), height = 6, width = 7.5)
print(p)
dev.off()

# per 1M expected barcode
p <- ggplot(stats_sub, aes(x = million_read, y = nUMI / expected_barcodes * 1e6, col = SETUP_ID))
p <- p + geom_line(data = stats_sub_bckgd, aes(group = paste(experiment, barcode_lib, SETUP_ID)), col = "lightgrey", alpha = .2) + geom_line(data = stats_sub) + geom_point()
p <- p + facet_grid(cols = vars(group), rows = vars(material))
p <- p + theme_plot() + ylab("Total Number of UMIs per 1M expected barcodes")
p <- p + scale_x_continuous(trans = "log10") + scale_y_continuous() + scale_color_manual(values = color_setup, name = "sample")
pdf(sprintf("%s/01d_subsampling/01_per_million_read/01d_saturation_nUMI_per1MBC_linear_vs_million_read.pdf", FIGURE_DIR), height = 6, width = 7.5)
print(p)
dev.off()
### Mean UMI per BC
# raw
p <- ggplot(, aes(x = million_read, y = nUMI / nBC, col = SETUP_ID))
p <- p + geom_line(data = stats_sub_bckgd, aes(group = paste(experiment, barcode_lib, SETUP_ID)), col = "lightgrey", alpha = .2) + geom_line(data = stats_sub) + geom_point()
p <- p + facet_grid(cols = vars(group), rows = vars(material))
p <- p + theme_plot() + ylab("Mean number of UMIs per barcodes")
p <- addScales(p) + scale_y_continuous(breaks = c(1, 3, 5, 10, 15, 20), limit = c(0.1, 20))
pdf(sprintf("%s/01d_subsampling/01_per_million_read/02a_saturation_meanUMI_perBC_vs_million_read.pdf", FIGURE_DIR), height = 6, width = 7.5)
print(p)
dev.off()

# winsorized
p <- ggplot(stats_sub, aes(x = million_read, y = winsor_mean, col = SETUP_ID))
p <- p + geom_line(data = stats_sub_bckgd, aes(group = paste(experiment, barcode_lib, SETUP_ID)), col = "lightgrey", alpha = .2) + geom_line(data = stats_sub) + geom_point()
p <- p + facet_grid(cols = vars(group), rows = vars(material))
p <- p + theme_plot() + ylab("Mean number of UMIs per barcodes (winsorized)")
p <- addScales(p) + scale_y_continuous(breaks = c(1, 3, 5, 10, 15, 20), limit = c(0.1, 20))
pdf(sprintf("%s/01d_subsampling/01_per_million_read/02b_saturation_winsorizedUMI_perBC_vs_million_read.pdf", FIGURE_DIR), height = 6, width = 7.5)
print(p)
dev.off()

### barcode recovery
# any
p <- ggplot(stats_sub, aes(x = million_read, y = Pct_BC_recovered, col = SETUP_ID))
p <- p + geom_line(data = stats_sub_bckgd, aes(group = paste(experiment, barcode_lib, SETUP_ID)), col = "lightgrey", alpha = .2) + geom_line(data = stats_sub) + geom_point()
p <- p + facet_grid(cols = vars(group), rows = vars(material))
p <- p + theme_plot() + ylab("Percentage of BC recovered")
p <- addScales(p, y_pct = TRUE)
pdf(sprintf("%s/01d_subsampling/01_per_million_read/03a_BC_recovery_vs_million_read.pdf", FIGURE_DIR), height = 6, width = 7.5)
print(p)
dev.off()

# >=3 UMI
p <- ggplot(stats_sub, aes(x = million_read, y = Pct_BC_over3UMI, col = SETUP_ID))
p <- p + geom_line(data = stats_sub_bckgd, aes(group = paste(experiment, barcode_lib, SETUP_ID)), col = "lightgrey", alpha = .2) + geom_line(data = stats_sub) + geom_point()
p <- p + facet_grid(cols = vars(group), rows = vars(material))
p <- p + theme_plot() + ylab("Percentage of BC recovered (3+ UMIs)")
p <- addScales(p, y_pct = TRUE)

pdf(sprintf("%s/01d_subsampling/01_per_million_read/03b_BC_recovery_3UMI_vs_million_read.pdf", FIGURE_DIR), height = 6, width = 7.5)
print(p)
dev.off()

### oligo recovery
# any
p <- ggplot(stats_sub, aes(x = million_read, y = Pct_oligo_recovered, col = SETUP_ID))
p <- p + geom_line(data = stats_sub_bckgd, aes(group = paste(experiment, barcode_lib, SETUP_ID)), col = "lightgrey", alpha = .2) + geom_line(data = stats_sub) + geom_point()
p <- p + facet_grid(cols = vars(group), rows = vars(material))
p <- p + theme_plot() + ylab("Percentage of oligo recovered")
p <- addScales(p, y_pct = TRUE)

pdf(sprintf("%s/01d_subsampling/01_per_million_read/04a_oligo_recovery_vs_million_read.pdf", FIGURE_DIR), height = 6, width = 7.5)
print(p)
dev.off()

# >10 BC
p <- ggplot(stats_sub, aes(x = million_read, y = Pct_oligo_over10BC, col = SETUP_ID))
p <- p + geom_line(data = stats_sub_bckgd, aes(group = paste(experiment, barcode_lib, SETUP_ID)), col = "lightgrey", alpha = .2) + geom_line(data = stats_sub) + geom_point()
p <- p + facet_grid(cols = vars(group), rows = vars(material))
p <- p + theme_plot() + ylab("Percentage of oligo recovered (>10BC)")
p <- addScales(p)
pdf(sprintf("%s/01d_subsampling/01_per_million_read/04b_oligo_recovery_10BC_vs_million_read.pdf", FIGURE_DIR), height = 6, width = 7.5)
print(p)
dev.off()

# >30 BC
p <- ggplot(stats_sub, aes(x = million_read, y = Pct_oligo_over30BC, col = SETUP_ID))
p <- p + geom_line(data = stats_sub_bckgd, aes(group = paste(experiment, barcode_lib, SETUP_ID)), col = "lightgrey", alpha = .2) + geom_line(data = stats_sub) + geom_point()
p <- p + facet_grid(cols = vars(group), rows = vars(material))
p <- p + theme_plot() + ylab("Percentage of oligo recovered (>30BC)")
p <- addScales(p, y_pct = TRUE)

pdf(sprintf("%s/01d_subsampling/01_per_million_read/04c_oligo_recovery_30BC_vs_million_read.pdf", FIGURE_DIR), height = 6, width = 7.5)
print(p)
dev.off()

# >100 BC
p <- ggplot(stats_sub, aes(x = million_read, y = Pct_oligo_over100BC, col = SETUP_ID))
p <- p + geom_line(data = stats_sub_bckgd, aes(group = paste(experiment, barcode_lib, SETUP_ID)), col = "lightgrey", alpha = .2) + geom_line(data = stats_sub) + geom_point()
p <- p + facet_grid(cols = vars(group), rows = vars(material))
p <- p + theme_plot() + ylab("Percentage of oligo recovered (>100BC)")
p <- addScales(p, y_pct = TRUE)

pdf(sprintf("%s/01d_subsampling/01_per_million_read/04d_oligo_recovery_100BC_vs_million_read.pdf", FIGURE_DIR), height = 6, width = 7.5)
print(p)
dev.off()

# >30 BC w 3UMI
p <- ggplot(stats_sub, aes(x = million_read, y = Pct_oligo_over30BCw3UMI, col = SETUP_ID))
p <- p + geom_line(data = stats_sub_bckgd, aes(group = paste(experiment, barcode_lib, SETUP_ID)), col = "lightgrey", alpha = .2) + geom_line(data = stats_sub) + geom_point()
p <- p + facet_grid(cols = vars(group), rows = vars(material))
p <- p + theme_plot() + ylab("Percentage of oligo recovered (>30BC with 3+ UMI)")
p <- addScales(p, y_pct = TRUE)

pdf(sprintf("%s/01d_subsampling/01_per_million_read/04e_oligo_recovery_30BCw3UMI_vs_million_read.pdf", FIGURE_DIR), height = 6, width = 7.5)
print(p)
dev.off()

### sequencing saturation
p <- ggplot(stats_sub, aes(x = million_read, y = sequencing_saturation / 100, col = SETUP_ID))
p <- p + geom_line(data = stats_sub_bckgd, aes(group = paste(experiment, barcode_lib, SETUP_ID)), col = "lightgrey", alpha = .2) + geom_line(data = stats_sub) + geom_point()
p <- p + facet_grid(cols = vars(group), rows = vars(material))
p <- p + theme_plot() + ylab("sequencing_saturation (% of duplicate reads)")
p <- addScales(p, y_pct = TRUE)

pdf(sprintf("%s/01d_subsampling/01_per_million_read/05_sequencing_saturation_vs_million_read.pdf", FIGURE_DIR), height = 6, width = 7.5)
print(p)
dev.off()

### sequencing saturation NEW READS
p <- ggplot(stats_sub, aes(x = million_read, y = sequencing_saturation_new / 100, col = SETUP_ID))
p <- p + geom_line(data = stats_sub_bckgd, aes(group = paste(experiment, barcode_lib, SETUP_ID)), col = "lightgrey", alpha = .2) + geom_line(data = stats_sub) + geom_point()
p <- p + facet_grid(cols = vars(group), rows = vars(material))
p <- p + theme_plot() + ylab("sequencing_saturation (% of duplicate reads)")
p <- addScales(p, y_pct = TRUE)

pdf(sprintf("%s/01d_subsampling/01_per_million_read/05b_sequencing_saturation_lastReads_vs_million_read.pdf", FIGURE_DIR), height = 6, width = 7.5)
print(p)
dev.off()



### newUMIs per 100 reads
p <- ggplot(stats_sub, aes(x = million_read, y = NewUMIper100read, col = SETUP_ID))
p <- p + geom_line(data = stats_sub_bckgd, aes(group = paste(experiment, barcode_lib, SETUP_ID)), col = "lightgrey", alpha = .2) + geom_line(data = stats_sub) + geom_point()
p <- p + facet_grid(cols = vars(group), rows = vars(material))
p <- p + theme_plot() + ylab("(Number of new UMI expected per 100 new reads")
p <- p + scale_x_continuous(trans = "log10") + scale_y_continuous() + scale_color_manual(values = color_setup, name = "sample")
pdf(sprintf("%s/01d_subsampling/01_per_million_read/06a_newUMI_per100reads_vs_million_read.pdf", FIGURE_DIR), height = 6, width = 7.5)
print(p)
dev.off()

# PctAdditionalUMI_for20pctMoreReads
p <- ggplot(stats_sub, aes(x = million_read, y = PctAdditionalUMI_for20pctMoreReads, col = SETUP_ID))
p <- p + geom_line(data = stats_sub_bckgd, aes(group = paste(experiment, barcode_lib, SETUP_ID)), col = "lightgrey", alpha = .2) + geom_line(data = stats_sub) + geom_point()
p <- p + facet_grid(cols = vars(group), rows = vars(material))
p <- p + theme_plot() + ylab("(Number of new UMI expected per 100 new reads")
p <- addScales(p, y_pct = TRUE)
pdf(sprintf("%s/01d_subsampling/01_per_million_read/06b_PctAdditionalUMI_for20pctMoreReads_vs_million_read.pdf", FIGURE_DIR), height = 6, width = 7.5)
print(p)
dev.off()

#############################################
######     percentage of total reads   ######
#############################################

### total number of UMI
# raw
p <- ggplot(stats_sub, aes(x = sub_pct / 100, y = nUMI, col = SETUP_ID))
p <- p + geom_line(data = stats_sub_bckgd, aes(group = paste(experiment, barcode_lib, SETUP_ID)), col = "lightgrey", alpha = .2) + geom_line(data = stats_sub) + geom_point()
p <- p + facet_grid(cols = vars(group), rows = vars(material))
p <- p + theme_plot() + ylab("Total Number of UMIs") + xlab("Percentage of reads")
p <- addScales(p, x_pct = TRUE)
pdf(sprintf("%s/01d_subsampling/02_percentage_total_read/01a_saturation_nUMI_vs_percentage_read.pdf", FIGURE_DIR), height = 6, width = 7.5)
print(p)
dev.off()

# per 1M expected barcode
p <- ggplot(stats_sub, aes(x = sub_pct / 100, y = nUMI / expected_barcodes * 1e6, col = SETUP_ID))
p <- p + geom_line(data = stats_sub_bckgd, aes(group = paste(experiment, barcode_lib, SETUP_ID)), col = "lightgrey", alpha = .2) + geom_line(data = stats_sub) + geom_point()
p <- p + facet_grid(cols = vars(group), rows = vars(material))
p <- p + theme_plot() + ylab("Total Number of UMIs per 1M expected barcodes") + xlab("Percentage of reads")
p <- addScales(p, x_pct = TRUE)
pdf(sprintf("%s/01d_subsampling/02_percentage_total_read/01b_saturation_nUMI_per1MBC_vs_percentage_read.pdf", FIGURE_DIR), height = 6, width = 7.5)
print(p)
dev.off()

### Mean UMI per BC
# raw
p <- ggplot(stats_sub, aes(x = sub_pct / 100, y = nUMI / nBC, col = SETUP_ID))
p <- p + geom_line(data = stats_sub_bckgd, aes(group = paste(experiment, barcode_lib, SETUP_ID)), col = "lightgrey", alpha = .2) + geom_line(data = stats_sub) + geom_point()
p <- p + facet_grid(cols = vars(group), rows = vars(material))
p <- p + theme_plot() + ylab("Mean number of UMIs per barcodes") + xlab("Percentage of reads")
p <- addScales(p, x_pct = TRUE) + scale_y_continuous(breaks = c(1, 3, 5, 10, 15, 20), limit = c(0.1, 20))
pdf(sprintf("%s/01d_subsampling/02_percentage_total_read/02a_saturation_meanUMI_perBC_vs_percentage_read.pdf", FIGURE_DIR), height = 6, width = 7.5)
print(p)
dev.off()

# winsorized
p <- ggplot(stats_sub, aes(x = sub_pct / 100, y = winsor_mean, col = SETUP_ID))
p <- p + geom_line(data = stats_sub_bckgd, aes(group = paste(experiment, barcode_lib, SETUP_ID)), col = "lightgrey", alpha = .2) + geom_line(data = stats_sub) + geom_point()
p <- p + facet_grid(cols = vars(group), rows = vars(material))
p <- p + theme_plot() + ylab("Mean number of UMIs per barcodes (winsorized)") + xlab("Percentage of reads")
p <- addScales(p, x_pct = TRUE) + scale_y_continuous(breaks = c(1, 3, 5, 10, 15, 20), limit = c(0.1, 20))
pdf(sprintf("%s/01d_subsampling/02_percentage_total_read/02b_saturation_winsorizedUMI_perBC_vs_percentage_read.pdf", FIGURE_DIR), height = 6, width = 7.5)
print(p)
dev.off()

### barcode recovery
# any
p <- ggplot(stats_sub, aes(x = sub_pct / 100, y = Pct_BC_recovered, col = SETUP_ID))
p <- p + geom_line(data = stats_sub_bckgd, aes(group = paste(experiment, barcode_lib, SETUP_ID)), col = "lightgrey", alpha = .2) + geom_line(data = stats_sub) + geom_point()
p <- p + facet_grid(cols = vars(group), rows = vars(material))
p <- p + theme_plot() + ylab("Percentage of BC recovered") + xlab("Percentage of reads")
p <- addScales(p, x_pct = TRUE, y_pct = TRUE)
pdf(sprintf("%s/01d_subsampling/02_percentage_total_read/03a_BC_recovery_vs_percentage_read.pdf", FIGURE_DIR), height = 6, width = 7.5)
print(p)
dev.off()

# >=3 UMI
p <- ggplot(stats_sub, aes(x = sub_pct / 100, y = Pct_BC_over3UMI, col = SETUP_ID))
p <- p + geom_line(data = stats_sub_bckgd, aes(group = paste(experiment, barcode_lib, SETUP_ID)), col = "lightgrey", alpha = .2) + geom_line(data = stats_sub) + geom_point()
p <- p + facet_grid(cols = vars(group), rows = vars(material))
p <- p + theme_plot() + ylab("Percentage of BC recovered (3+ UMIs)") + xlab("Percentage of reads")
p <- addScales(p, x_pct = TRUE, y_pct = TRUE)
pdf(sprintf("%s/01d_subsampling/02_percentage_total_read/03b_BC_recovery_3UMI_vs_percentage_read.pdf", FIGURE_DIR), height = 6, width = 7.5)
print(p)
dev.off()

### oligo recovery
# any
p <- ggplot(stats_sub, aes(x = sub_pct / 100, y = Pct_oligo_recovered, col = SETUP_ID))
p <- p + geom_line(data = stats_sub_bckgd, aes(group = paste(experiment, barcode_lib, SETUP_ID)), col = "lightgrey", alpha = .2) + geom_line(data = stats_sub) + geom_point()
p <- p + facet_grid(cols = vars(group), rows = vars(material))
p <- p + theme_plot() + ylab("Percentage of oligo recovered") + xlab("Percentage of reads")
p <- addScales(p, x_pct = TRUE, y_pct = TRUE)
pdf(sprintf("%s/01d_subsampling/02_percentage_total_read/04a_oligo_recovery_vs_percentage_read.pdf", FIGURE_DIR), height = 6, width = 7.5)
print(p)
dev.off()

# >10 BC
p <- ggplot(stats_sub, aes(x = sub_pct / 100, y = Pct_oligo_over10BC, col = SETUP_ID))
p <- p + geom_line(data = stats_sub_bckgd, aes(group = paste(experiment, barcode_lib, SETUP_ID)), col = "lightgrey", alpha = .2) + geom_line(data = stats_sub) + geom_point()
p <- p + facet_grid(cols = vars(group), rows = vars(material))
p <- p + theme_plot() + ylab("Percentage of oligo recovered (>10BC)") + xlab("Percentage of reads")
p <- addScales(p, x_pct = TRUE, y_pct = TRUE)
pdf(sprintf("%s/01d_subsampling/02_percentage_total_read/04b_oligo_recovery_10BC_vs_percentage_read.pdf", FIGURE_DIR), height = 6, width = 7.5)
print(p)
dev.off()

# >30 BC
p <- ggplot(stats_sub, aes(x = sub_pct / 100, y = Pct_oligo_over30BC, col = SETUP_ID))
p <- p + geom_line(data = stats_sub_bckgd, aes(group = paste(experiment, barcode_lib, SETUP_ID)), col = "lightgrey", alpha = .2) + geom_line(data = stats_sub) + geom_point()
p <- p + facet_grid(cols = vars(group), rows = vars(material))
p <- p + theme_plot() + ylab("Percentage of oligo recovered (>30BC)") + xlab("Percentage of reads")
p <- addScales(p, x_pct = TRUE, y_pct = TRUE)
pdf(sprintf("%s/01d_subsampling/02_percentage_total_read/04c_oligo_recovery_30BC_vs_percentage_read.pdf", FIGURE_DIR), height = 6, width = 7.5)
print(p)
dev.off()

# >100 BC
p <- ggplot(stats_sub, aes(x = sub_pct / 100, y = Pct_oligo_over100BC, col = SETUP_ID))
p <- p + geom_line(data = stats_sub_bckgd, aes(group = paste(experiment, barcode_lib, SETUP_ID)), col = "lightgrey", alpha = .2) + geom_line(data = stats_sub) + geom_point()
p <- p + facet_grid(cols = vars(group), rows = vars(material))
p <- p + theme_plot() + ylab("Percentage of oligo recovered (>100BC)") + xlab("Percentage of reads")
p <- addScales(p, x_pct = TRUE, y_pct = TRUE)
pdf(sprintf("%s/01d_subsampling/02_percentage_total_read/04d_oligo_recovery_100BC_vs_percentage_read.pdf", FIGURE_DIR), height = 6, width = 7.5)
print(p)
dev.off()

# >30 BC w 3UMI
p <- ggplot(stats_sub, aes(x = sub_pct / 100, y = Pct_oligo_over30BCw3UMI, col = SETUP_ID))
p <- p + geom_line(data = stats_sub_bckgd, aes(group = paste(experiment, barcode_lib, SETUP_ID)), col = "lightgrey", alpha = .2) + geom_line(data = stats_sub) + geom_point()
p <- p + facet_grid(cols = vars(group), rows = vars(material))
p <- p + theme_plot() + ylab("Percentage of oligo recovered (>30BC with 3 UMI or more)") + xlab("Percentage of reads")
p <- addScales(p, x_pct = TRUE, y_pct = TRUE)
pdf(sprintf("%s/01d_subsampling/02_percentage_total_read/04e_oligo_recovery_30BCw3UMI_vs_percentage_read.pdf", FIGURE_DIR), height = 6, width = 7.5)
print(p)
dev.off()

### sequencing saturation
p <- ggplot(stats_sub, aes(x = sub_pct / 100, y = sequencing_saturation / 100, col = SETUP_ID))
p <- p + geom_line(data = stats_sub_bckgd, aes(group = paste(experiment, barcode_lib, SETUP_ID)), col = "lightgrey", alpha = .2) + geom_line(data = stats_sub) + geom_point()
p <- p + facet_grid(cols = vars(group), rows = vars(material))
p <- p + theme_plot() + ylab("sequencing_saturation (% of duplicate reads)") + xlab("Percentage of reads")
p <- addScales(p, y_pct = TRUE, x_pct = TRUE)
pdf(sprintf("%s/01d_subsampling/02_percentage_total_read/05_sequencing_saturation_vs_percentage_read.pdf", FIGURE_DIR), height = 6, width = 7.5)
print(p)
dev.off()

### sequencing saturation
p <- ggplot(stats_sub, aes(x = sub_pct / 100, y = sequencing_saturation / 100, col = SETUP_ID))
p <- p + geom_line(data = stats_sub_bckgd, aes(group = paste(experiment, barcode_lib, SETUP_ID)), col = "lightgrey", alpha = .2) + geom_line(data = stats_sub) + geom_point()
p <- p + facet_grid(cols = vars(group), rows = vars(material))
p <- p + theme_plot() + ylab("sequencing_saturation (% of duplicate reads)") + xlab("Percentage of reads")
p <- addScales(p, y_pct = TRUE, x_pct = TRUE)
pdf(sprintf("%s/01d_subsampling/02_percentage_total_read/05_sequencing_saturation_vs_percentage_read.pdf", FIGURE_DIR), height = 6, width = 7.5)
print(p)
dev.off()


### sequencing saturation
p <- ggplot(stats_sub, aes(x = sub_pct / 100, y = sequencing_saturation_new / 100, col = SETUP_ID))
p <- p + geom_line(data = stats_sub_bckgd, aes(group = paste(experiment, barcode_lib, SETUP_ID)), col = "lightgrey", alpha = .2) + geom_line(data = stats_sub) + geom_point()
p <- p + facet_grid(cols = vars(group), rows = vars(material))
p <- p + theme_plot() + ylab("sequencing_saturation (% of duplicate reads)") + xlab("Percentage of reads")
p <- addScales(p, y_pct = TRUE, x_pct = TRUE)
pdf(sprintf("%s/01d_subsampling/02_percentage_total_read/05b_sequencing_saturation_lastReads_vs_percentage_read.pdf", FIGURE_DIR), height = 6, width = 7.5)
print(p)
dev.off()

### newUMIs per 100 reads
p <- ggplot(stats_sub, aes(x = sub_pct / 100, y = NewUMIper100read, col = SETUP_ID))
p <- p + geom_line(data = stats_sub_bckgd, aes(group = paste(experiment, barcode_lib, SETUP_ID)), col = "lightgrey", alpha = .2) + geom_line(data = stats_sub) + geom_point()
p <- p + facet_grid(cols = vars(group), rows = vars(material))
p <- p + theme_plot() + ylab("(Number of new UMI expected per 100 new reads") + xlab("Percentage of reads")
p <- p + scale_x_continuous(labels = scales::percent, trans = "log10") + scale_y_continuous() + scale_color_manual(values = color_setup, name = "sample")
pdf(sprintf("%s/01d_subsampling/02_percentage_total_read/06a_newUMI_per100reads_vs_percentage_read.pdf", FIGURE_DIR), height = 6, width = 7.5)
print(p)
dev.off()

# PctAdditionalUMI_for20pctMoreReads
p <- ggplot(stats_sub, aes(x = sub_pct / 100, y = PctAdditionalUMI_for20pctMoreReads, col = SETUP_ID))
p <- p + geom_line(data = stats_sub_bckgd, aes(group = paste(experiment, barcode_lib, SETUP_ID)), col = "lightgrey", alpha = .2) + geom_line(data = stats_sub) + geom_point()
p <- p + facet_grid(cols = vars(group), rows = vars(material))
p <- p + theme_plot() + ylab("(Number of new UMI expected per 100 new reads") + xlab("Percentage of reads")
p <- addScales(p, y_pct = TRUE, x_pct = TRUE)
pdf(sprintf("%s/01d_subsampling/02_percentage_total_read/06b_PctAdditionalUMI_for20pctMoreReads_vs_percentage_read.pdf", FIGURE_DIR), height = 6, width = 7.5)
print(p)
dev.off()


########################################
######   AGGREGATED RESULTS       ######
########################################

### newUMIs per 100 reads
p <- ggplot(stats_sub, aes(x = million_read, y = NewUMIper100read, col = SETUP_ID))
p <- p + geom_point() + geom_line(alpha = .2)
p <- p + facet_grid(rows = vars(material))
p <- p + theme_plot() + ylab("(Number of new UMI expected per 100 new reads") + xlab("Number of reads (millions)")
p <- p + scale_x_continuous() + scale_y_continuous() + scale_color_manual(values = color_setup, name = "sample")
p <- p + geom_hline(yintercept = 10, col = "darkred", linetype = 2, alpha = .5) + geom_hline(yintercept = 5, col = "red", linetype = 1)
pdf(sprintf("%s/01d_subsampling/06a_newUMI_per100reads_vs_millions_of_reads_aggregated.pdf", FIGURE_DIR), height = 8, width = 6)
print(p)
dev.off()

### sequencing saturation NEW READS
p <- ggplot(stats_sub, aes(x = million_read, y = sequencing_saturation_new / 100, col = SETUP_ID))
p <- p + geom_point() + geom_line(alpha = .2)
p <- p + facet_grid(rows = vars(material))
p <- p + theme_plot() + ylab("sequencing_saturation (% of duplicate reads)")
p <- p + scale_x_continuous() + scale_y_continuous() + scale_color_manual(values = color_setup, name = "sample")
p <- p + geom_hline(yintercept = .8, col = "darkred", linetype = 2, alpha = .5) + geom_hline(yintercept = .9, col = "red", linetype = 1)

pdf(sprintf("%s/01d_subsampling/05b_sequencing_saturation_lastReads_vs_million_read_aggregated.pdf", FIGURE_DIR), height = 8, width = 6)
print(p)
dev.off()

cat("All Done!\n")
q('no')
