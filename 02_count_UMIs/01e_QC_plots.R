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

LIB_summary[, SETUP_ID := factor(SETUP_ID, names(color_setup))]
LIB_summary <- LIB_summary[order(SETUP_ID, material)]
LIB_summary <- LIB_summary[!celltype%chin%c('Calu3','A549'),]
color_setup <- color_setup[!names(color_setup)%chin%c("Calu3_NS_R1","A549_NS_R1")]

# load oligo annotations
oligo_source <- fread(sprintf("%s/data/%s/00_oligo_annot_v1.txt", MPRA_DIR, ANALYSIS_DIR))
oligo_source[, allele.num := paste0("A", seq_len(.N), sep = ""), by = .(crsID, shift, strand)]
oligo_type <- oligo_source[, .(type = paste(unique(type), collapse = "\\")), by = oligo]
# list oligos associated with >1 source
dup_oligo <- oligo_source[duplicated(oligo), unique(oligo)]

###############################################################################
###########  reading/formatting read/UMI counts             ###################
###############################################################################

Counts_by_lib <- list()
for (i in seq_len(LIB_summary[, .N])) {
  SID <- LIB_summary[i, SID]
  SAMPLE_ID <- LIB_summary[i, library]
  # load UMI per BC
  Counts_by_lib[[SAMPLE_ID]] <- fread(sprintf("%s/data/%s/01b_Count_reads_UMI_type_frequency/Count_reads_UMI_type_frequency__%s__v%s.tsv", MPRA_DIR, ANALYSIS_DIR, SAMPLE_ID, VERSION_IN))
}
Counts_by_lib <- rbindlist(Counts_by_lib, idcol = "library")
Counts_by_lib <- merge(Counts_by_lib, LIB_summary, by = "library")

Counts_by_lib[, read_class := case_when(
  !sameBC & !known_BC1 & !known_BC2 ~ "discordant barcodes, unknown",
  !sameBC & (known_BC1 | known_BC2) ~ "discordant barcodes, partially known",
  sameBC & known_BC1 & shannon & shannon_UMI ~ "known BC",
  sameBC & known_BC1 & shannon & !shannon_UMI ~ "known BC, low complexity UMI",
  sameBC & !known_BC1 & shannon ~ "unknown BC",
  sameBC & !known_BC1 & !shannon ~ "unknown BC, low complexity"
)]
Counts_by_lib <- Counts_by_lib[, .(N = sum(N) / 1e6, nUMI = sum(nUMI) / 1e6), by = .(library, read_class, experiment, barcode_lib)]
Counts_by_lib[, library := factor(library, LIB_summary$library)]
Counts_by_lib <- Counts_by_lib[order(experiment, barcode_lib, library, read_class), ]
###############################################################################
###########            plotting read/UMI counts             ###################
###############################################################################
dir.create(sprintf("%s/01b_UMI_filtering_and_recovery/",FIGURE_DIR))
pdf(sprintf("%s/01b_UMI_filtering_and_recovery/00_01_detail_reads_perlib.pdf", FIGURE_DIR), height = 5, width = 7)
p <- ggplot(Counts_by_lib, aes(x = library, y = N, fill = factor(read_class, names(color_readclass))))
p <- p + geom_bar(stat = "Identity", position = "stack") + ylab("read counts (million)")
p <- p + theme_plot(rotate.x = 90) + scale_fill_manual(values = color_readclass) + facet_grid(col = vars(paste(experiment, "-", barcode_lib)), space = "free", scales = "free")
print(p)
dev.off()

pdf(sprintf("%s/01b_UMI_filtering_and_recovery/00_02_detail_UMI_perlib.pdf", FIGURE_DIR), height = 5, width = 7)
p <- ggplot(Counts_by_lib, aes(x = library, y = nUMI, fill = factor(read_class, names(color_readclass))))
p <- p + geom_bar(stat = "Identity", position = "stack") + ylab("UMI counts (million)")
p <- p + theme_plot(rotate.x = 90) + scale_fill_manual(values = color_readclass) + facet_grid(col = vars(paste(experiment, "-", barcode_lib)), space = "free", scales = "free")
print(p)
dev.off()
###############################################################################
###########     reading & formatting PctDup UMIs       ########################
###############################################################################

UMI_perBC <- list()
Pct_dupUMI <- list()
for (i in seq_len(LIB_summary[, .N])) {
  SID <- LIB_summary[i, SID]
  SAMPLE_ID <- LIB_summary[i, library]
  # load UMI per BC
  UMI_perBC[[i]] <- fread(sprintf("%s/data/%s/01b_MPRA_results_counts_BC/MPRA_results_counts_BC__%s__v%s.txt.gz", MPRA_DIR, ANALYSIS_DIR, SAMPLE_ID, VERSION_IN))
  # load stats on pct duplicated UMIs
  Pct_dupUMI[[i]] <- fread(sprintf("%s/data/%s/01b_Pct_dupUMI/Pct_dupUMI__%s__v%s.tsv", MPRA_DIR, ANALYSIS_DIR, SAMPLE_ID, VERSION_IN))
}
UMI_perBC <- rbindlist(UMI_perBC)
UMI_perBC <- merge(UMI_perBC, LIB_summary, by = "library")
UMI_perBC[, BC_class := ifelse(shannon < 10, "low complexity", ifelse(known_BC1, "known", "unknown"))]
UMI_perBC[, SETUP_ID := factor(SETUP_ID, names(color_setup))]

Pct_dupUMI <- rbindlist(Pct_dupUMI)
Pct_dupUMI <- merge(Pct_dupUMI, LIB_summary, by = "library")
Pct_dupUMI[, complexity := factor(complexity, unique(complexity))]
Pct_dupUMI[, SETUP_ID := factor(SETUP_ID, names(color_setup))]

###############################################################################
###########            plotting PctDup UMIs            ########################
###############################################################################

### duplicated index
p <- ggplot(Pct_dupUMI, aes(
  x = complexity,
  y = Pct_dup_index,
  ymin = Pct_dup_index_low,
  ymax = Pct_dup_index_high,
  col = SETUP_ID
))
p <- p + geom_pointrange() + theme_plot(rotate.x = 90, lpos = "right") + facet_grid(cols = vars(SETUP_ID), rows = vars(material))
p <- p + scale_color_manual(values = color_setup) + ylab("Pct of UMIs associated with >1 index")

pdf(sprintf("%s/01b_UMI_filtering_and_recovery/01_01_Pct_dupUMI.pdf", FIGURE_DIR), width = 2 + 2 * N_sample, height = 4)
print(p)
dev.off()

### duplicated BC
p <- ggplot(Pct_dupUMI, aes(
  x = complexity,
  y = Pct_dup_BC,
  ymin = Pct_dup_BC_low,
  ymax = Pct_dup_BC_high,
  col = SETUP_ID
))
p <- p + geom_pointrange() + theme_plot(rotate.x = 90, lpos = "right") + facet_grid(cols = vars(SETUP_ID), rows = vars(material))
p <- p + scale_color_manual(values = color_setup) + ylab("Pct of UMIs associated with >1 BC")
pdf(sprintf("%s/01b_UMI_filtering_and_recovery/01_02_Pct_dupBC_perUMI.pdf", FIGURE_DIR), width = 2 + 2 * N_sample, height = 4)
print(p)
dev.off()

### number unique/duplicated BC per complexity BIN
FigData <- melt(Pct_dupUMI[, .(library, complexity, Nb_dup_BC, Nb_uniq_BC, material, SETUP_ID)], id.vars = c("library", "complexity", "SETUP_ID", "material"))

p <- ggplot(FigData, aes(
  x = complexity, y = value,
  fill = ifelse(grepl("uniq", variable), "unique", "duplicated")
))
p <- p + geom_bar(stat = "Identity", position = "stack") + theme_plot(rotate.x = 90, lpos = "right")
p <- p + facet_grid(cols = vars(SETUP_ID), rows = vars(material))
p <- p + ylab("Nb of UMIs")
pdf(sprintf("%s/01b_UMI_filtering_and_recovery/01_03_Nb_dupBC_perUMI.pdf", FIGURE_DIR), width = 2 + 1.5 * N_sample, height = 4)
print(p)
dev.off()

###############################################################################
###########            plotting UMI counts             ########################
###############################################################################
UMI_perBC[, library := factor(library, LIB_summary$library)]

UMI_perBC_mean <- UMI_perBC[, .(nUMI_perBC = mean(nUMI_perBC)), keyby = library]
fwrite(UMI_perBC_mean, file = sprintf("%s/01b_UMI_filtering_and_recovery/02_01_UMI_perBC_mean_by_lib.tsv", FIGURE_DIR), sep = "\t")

p <- ggplot(UMI_perBC, aes(x = library, y = nUMI_perBC, fill = SETUP_ID, alpha = material)) +
  scale_y_continuous(trans = "log10")
p <- p + rasterize(geom_violin(scale = "width"), dpi = 200) + rasterize(geom_boxplot(notch = TRUE, fill = "#FFFFFF", alpha = .5), dpi = 200)
p <- p + ylab("Number of UMIs per barcode") + xlab("library")
p <- p + scale_fill_manual(values = color_setup) + scale_alpha_manual(values = c(DNA = 1, RNA = 0.5), guide = "none")
p <- p + theme_plot(rotate.x = 90) + guides(fill = 'none')
pdf(sprintf("%s/01b_UMI_filtering_and_recovery/02_02_UMI_perBC_by_lib.pdf", FIGURE_DIR), height = 5, width = 2 + N_sample / 3)
print(p)
dev.off()

UMI_per_oligo <- UMI_perBC[BC_type == "1a_assoc_unique" & BC_class == "known",
  .(nUMI = sum(nUMI_perBC), nBC = length(unique(barcode1))),
  keyby = .(oligo, library, celltype, material, condition, replicate, SETUP_ID)
]
LIB_summary <- LIB_summary[order(factor(celline, c("Calu3", "A549", "HepG2", "K562")),celltype, factor(condition, c("NS", "IFNA2b", "DEX", "SARS", "IAV",'TNFa')), replicate, material)]
UMI_per_oligo[, library := factor(library, LIB_summary$library)]

p <- ggplot(UMI_per_oligo, aes(x = library, y = nUMI, fill = SETUP_ID, alpha = material)) +
  scale_y_continuous(trans = "log10")
p <- p + rasterize(geom_violin(scale = "width"), dpi = 200) + rasterize(geom_boxplot(notch = TRUE, fill = "#FFFFFF", alpha = .5), dpi = 200)
p <- p + ylab("Number of UMIs per oligo") + xlab("library")
p <- p + scale_fill_manual(values = color_setup) + scale_alpha_manual(values = c(DNA = 1, RNA = 0.5), guide = "none")
p <- p + theme_plot(rotate.x = 90) + guides(fill = 'none')

pdf(sprintf("%s/01b_UMI_filtering_and_recovery/02_03_UMI_per_oligo_by_lib.pdf", FIGURE_DIR), height = 5, width = 2 + N_sample / 3)
print(p)
dev.off()

library(DescTools)
UMI_per_oligo_mean <- UMI_per_oligo[, .(nUMI_per_oligo = mean(Winsorize(nUMI)), nBC_per_oligo = mean(Winsorize(nBC))), keyby = library]
UMI_per_oligo_mean <- merge(UMI_per_oligo_mean, LIB_summary, by = "library")

# numbers for Labmeeting 07/10/2024
UMI_per_oligo_mean[, .(min(nUMI_per_oligo), mean(nUMI_per_oligo), max(nUMI_per_oligo)), keyby = .(barcode_lib, material)]
UMI_per_oligo_mean[, .(min(nBC_per_oligo), mean(nBC_per_oligo), max(nBC_per_oligo)), keyby = .(barcode_lib)]
UMI_per_oligo_mean[, .(min(nUMI_per_oligo / nBC_per_oligo), mean(nUMI_per_oligo / nBC_per_oligo), max(nUMI_per_oligo / nBC_per_oligo)), keyby = .(barcode_lib, material)]


UMI_per_oligo_mean <- UMI_per_oligo[, .(nUMI_per_oligo = mean(nUMI), nBC_per_oligo = mean(nBC)), keyby = library]
UMI_per_oligo_mean <- merge(UMI_per_oligo_mean, LIB_summary, by = "library")
# for future refrence
UMI_per_oligo_mean[, .(min(nUMI_per_oligo), mean(nUMI_per_oligo), max(nUMI_per_oligo)), keyby = .(barcode_lib, material)]
UMI_per_oligo_mean[, .(min(nBC_per_oligo), mean(nBC_per_oligo), max(nBC_per_oligo)), keyby = .(barcode_lib)]
UMI_per_oligo_mean[, .(min(nUMI_per_oligo / nBC_per_oligo), mean(nUMI_per_oligo / nBC_per_oligo), max(nUMI_per_oligo / nBC_per_oligo)), keyby = .(barcode_lib, material)]

fwrite(UMI_per_oligo_mean, file = sprintf("%s/01b_UMI_filtering_and_recovery/02_04_UMI_BC_per_oligo_mean_by_lib.tsv", FIGURE_DIR), sep = "\t")



p <- ggplot(UMI_per_oligo, aes(x = library, y = nBC, fill = SETUP_ID, alpha = material)) +
  scale_y_continuous(trans = "log10")
p <- p + rasterize(geom_violin(scale = "width"), dpi = 200) + rasterize(geom_boxplot(notch = TRUE, fill = "#FFFFFF", alpha = .5), dpi = 200)
p <- p + ylab("Number of barcodes per oligo") + xlab("library")
p <- p + scale_fill_manual(values = color_setup) + scale_alpha_manual(values = c(DNA = 1, RNA = 0.5), guide = "none")
p <- p + theme_plot(rotate.x = 90) +guides(fill="none")

pdf(sprintf("%s/01b_UMI_filtering_and_recovery/02_05_BC_per_oligo_by_lib.pdf", FIGURE_DIR), height = 5, width = 2 + N_sample / 3)
print(p)
dev.off()


###############################################################################
###########            reading BC recovery                  ###################
###############################################################################

BC_recovery <- list()
for (i in seq_len(LIB_summary[, .N])) {
  SID <- LIB_summary[i, as.character(SID)]
  SAMPLE_ID <- LIB_summary[i, as.character(library)]
  # load oligo activity
  BC_recovery[[SAMPLE_ID]] <- fread(sprintf("%s/data/%s/01b_BC_recovery/BC_recovery__%s__v%s.txt", MPRA_DIR, ANALYSIS_DIR, SAMPLE_ID, VERSION_IN))
}
BC_recovery <- rbindlist(BC_recovery, idcol = "library")
BC_recovery <- merge(BC_recovery, unique(LIB_summary[, .(library, celltype, condition, replicate, SETUP_ID)]), by = "library")
BC_recovery_long <- melt(BC_recovery, id.vars = c("library", "celltype", "condition", "replicate", "SETUP_ID"))
BC_recovery_long[, type := case_when(
  grepl("Pct_oligo", variable) ~ "Pct_oligo",
  grepl("Pct_BC", variable) ~ "Pct_BC", TRUE ~ NA
)]
BC_recovery_long[, group := gsub("Pct_(oligo|BC)_", "", variable)]

p <- ggplot(BC_recovery_long, aes(x = library, fill = group, y = value)) +
  geom_bar(stat = "Identity", position = "dodge")
p <- p + theme_plot(rotate.x = 90) + facet_grid(cols = vars(type))
pdf(sprintf("%s/01b_UMI_filtering_and_recovery/03_01_recovery_BC_and_oligo.pdf", FIGURE_DIR), height = 3, width = 1 + N_sample / 3)
print(p)
dev.off()

cat("All Done!\n")
q('no')



BC_perCRS <- UMI_perBC[BC_type == "1a_assoc_unique", .(nBC = length(unique(barcode1))), keyby = .(CRS, SETUP_ID, COND_ID, experiment, celline, condition, material, replicate, BC_class, barcode_library = barcode_lib)]
BC_perCRS <- BC_perCRS[CRS != "", ]
BC_perCRS_SETUP <- UMI_perBC[BC_type == "1a_assoc_unique", .(nIntegratedBC = length(unique(barcode1))), keyby = .(CRS, celline, condition, replicate, BC_class, barcode_library = barcode_lib)]
BC_perCRS_SETUP <- BC_perCRS_SETUP[CRS != "", ]
BC_perCRS_total <- Associations_Filtered[, .(nKnownBC = length(unique(BC))), by = .(CRS, barcode_library)]

BC_perCRS <- merge(BC_perCRS, BC_perCRS_total, by = c("CRS", "barcode_library"))
BC_perCRS <- merge(BC_perCRS, BC_perCRS_SETUP, by = c("CRS", "barcode_library","celline","condition","replicate","BC_class"))
BC_perCRS[,PCT_Integrated := nIntegratedBC/nKnownBC]
BC_perCRS[,PCT_recovered := nBC/nIntegratedBC]
