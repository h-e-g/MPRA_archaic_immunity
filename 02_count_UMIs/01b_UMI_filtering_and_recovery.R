MPRA_DIR <- "/pasteur/helix/projects/evo_immuno_pop/MPRA"
ANALYSIS_DIR <- "MPRA_count_exp6_analysis3"

source(sprintf("%s/scripts/%s/00_starter.R", MPRA_DIR, ANALYSIS_DIR))
source(sprintf("%s/scripts/%s/00_load_associations.R", MPRA_DIR, ANALYSIS_DIR))

VERSION_IN <- "1"
VERSION_OUT <- "1"
SUB_SAMPLE <- -1
SUB_SAMPLE_Pct <- FALSE

cmd <- commandArgs()
print(cmd)
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

SAMPLE_ID <- Indexes[, unique(Sample_ID)[j]]
SID <- Indexes[Sample_ID == SAMPLE_ID, head(SID, 1)]
BC_LIB <- Indexes[Sample_ID == SAMPLE_ID, head(barcode_lib, 1)]

################################################################################
##############                  BC per UMI  plots                    ###########
################################################################################

# load MPRA_reads (recommended: >400 G)
MPRA_reads <- sprintf("%s/data/%s/01a_MPRA_results/MPRA_results__%s__v%s.txt.gz", MPRA_DIR, ANALYSIS_DIR, SAMPLE_ID, VERSION_IN)
DT_MPRA <- fread(MPRA_reads)
DT_MPRA[, library := SAMPLE_ID]
DT_MPRA[, .N, by = cut(shannon_UMI, 10)]

if (SUB_SAMPLE > 0) {
  # sub sample the reads
  if (SUB_SAMPLE_Pct == TRUE) {
    DT_MPRA <- DT_MPRA[sample(seq_len(.N), SUB_SAMPLE * .N, replace = FALSE), ]
  } else {
    DT_MPRA <- DT_MPRA[sample(seq_len(.N), min(SUB_SAMPLE, .N), replace = FALSE), ]
  }
  # define subsampling directory & update analysis directory
  if (SUB_SAMPLE_Pct == TRUE) {
    ANALYSIS_DIR <- sprintf("%s/01d_subsampling/SUB_%sPct", ANALYSIS_DIR, 100 * SUB_SAMPLE)
  } else {
    ANALYSIS_DIR <- sprintf("%s/01d_subsampling/SUB_%s", ANALYSIS_DIR, SUB_SAMPLE)
  }
  dir.create(ANALYSIS_DIR, showWarnings = FALSE, recursive = TRUE)
}
# find duplicated UMIs (per library)
BC_perUMI_strict <- DT_MPRA[sameBC == TRUE, .(nBC = length(unique(barcode1))), by = .(library, UMI)]
dup_UMIs <- BC_perUMI_strict[nBC > 1, unique(UMI)]

# count barcodes per UMIs (per library)
BC_perUMI <- DT_MPRA[sameBC == TRUE, .(.N,
  nBC = length(unique(barcode1)),
  nIndex = length(unique(index)),
  nBC_index = length(unique(paste(index, barcode1)))
), by = .(library, UMI, cut(shannon_UMI, 10))]

dir.create(sprintf("%s/data/%s/01b_BC_perUMI", MPRA_DIR, ANALYSIS_DIR), showWarnings = FALSE, recursive = TRUE)
fwrite(BC_perUMI, file = sprintf("%s/data/%s/01b_BC_perUMI/BC_perUMI__%s__v%s.tsv.gz", MPRA_DIR, ANALYSIS_DIR, SAMPLE_ID, VERSION_OUT), sep = "\t")


Pct_dupUMI <- BC_perUMI[, .(.N,
  Nb_dup_BC = sum(nBC > 1),
  Nb_dup_index = sum(nIndex > 1)
), keyby = .(library, complexity = cut)]

Pct_dupUMI[, Pct_dup_index_low := binom.test(Nb_dup_index, N)$conf.int[1], by = .(library, complexity)]
Pct_dupUMI[, Pct_dup_index_high := binom.test(Nb_dup_index, N)$conf.int[2], by = .(library, complexity)]
Pct_dupUMI[, Pct_dup_index := Nb_dup_index / N]

Pct_dupUMI[, Pct_dup_BC_low := binom.test(Nb_dup_BC, N)$conf.int[1], by = .(library, complexity)]
Pct_dupUMI[, Pct_dup_BC_high := binom.test(Nb_dup_BC, N)$conf.int[2], by = .(library, complexity)]
Pct_dupUMI[, Pct_dup_BC := Nb_dup_BC / N]
Pct_dupUMI[, Nb_uniq_BC := N - Nb_dup_BC]

dir.create(sprintf("%s/data/%s/01b_Pct_dupUMI", MPRA_DIR, ANALYSIS_DIR), showWarnings = FALSE, recursive = TRUE)
fwrite(Pct_dupUMI, file = sprintf("%s/data/%s/01b_Pct_dupUMI/Pct_dupUMI__%s__v%s.tsv", MPRA_DIR, ANALYSIS_DIR, SAMPLE_ID, VERSION_OUT), sep = "\t")

###### alignment of unmatched barcodes : in most cases there is low quality alignment between the two barcodes.
# BC1=DT_MPRA[sameBC==FALSE,barcode1]
# BC2=DT_MPRA[sameBC==FALSE,barcode2]
# aln_scores=c()
# samp=sample(1:length(BC1),1000)
# tic()
# for (i in 1: 1000){
#   cat(i,'')
#   aln_scores[i]=pairwiseAlignment(BC1[samp[i]],reverseComplement(DNAString(BC2[samp[i]])))@score
# }
# toc()
##### DONE: find the cause of duplicated UMIs across indexes, libraries and BC, and remove them early on
##### duplicated UMIs are mainly due to low complexity. However filtering these low complexity UMIs is unnecessarily stringent and removes a high number of UMIs
# we have thus decided to keep low complexity UMIs as long as they are not in duplicate for the library under consideration

Count_reads_UMI_type_frequency <- DT_MPRA[, .(.N,
  nUMI = length(unique(UMI)),
  nUMI_BC = length(unique(paste(barcode1, UMI))),
  nUMI_BC_index = length(unique(paste(barcode1, UMI, index)))
), by = .(known_BC1, sameBC, known_BC2, shannon >= 11, shannon_UMI >= 11)][order(-N)]
dir.create(sprintf("%s/data/%s/01b_Count_reads_UMI_type_frequency", MPRA_DIR, ANALYSIS_DIR), showWarnings = FALSE, recursive = TRUE)
fwrite(Count_reads_UMI_type_frequency, file = sprintf("%s/data/%s/01b_Count_reads_UMI_type_frequency/Count_reads_UMI_type_frequency__%s__v%s.tsv", MPRA_DIR, ANALYSIS_DIR, SAMPLE_ID, VERSION_OUT), sep = "\t")

# Count_reads_type_frequency[,Pct_total:=round(N/sum(N)*100,1)]
# Count_reads_type_frequency[,sum(Pct_total),by=.(knownIndex)]
# Count_reads_type_frequency[,sum(Pct_total),by=.(knownIndex_proxy)]
# Count_reads_type_frequency[,sum(Pct_total),by=.(known_BC1)]
# Count_reads_type_frequency[,sum(Pct_total),by=.(sameBC)]
# Count_reads_type_frequency[,sum(Pct_total),by=.(known_BC2)]
# Count_reads_type_frequency[,sum(Pct_total),by=.(shannon)]
# Count_reads_type_frequency[,sum(Pct_total),by=.(shannon_UMI)]


##### NOT DONE: filter out low complexity UMIs and low complexity barcodes (too much filtering)
# filter out duplicated UMIs (if the same UMI is found for multiples barcodes, we filter the UMI)

DT_MPRA_clean <- DT_MPRA[sameBC & known_BC1 & known_BC2 & !UMI %chin% dup_UMIs, ]
dir.create(sprintf("%s/data/%s/01b_MPRA_results_clean", MPRA_DIR, ANALYSIS_DIR), showWarnings = FALSE, recursive = TRUE)
fwrite(DT_MPRA_clean, sprintf("%s/data/%s/01b_MPRA_results_clean/MPRA_results_clean__%s__v%s.txt.gz", MPRA_DIR, ANALYSIS_DIR, SAMPLE_ID, VERSION_OUT), sep = "\t")

UMI_perBC <- DT_MPRA_clean[, .(nUMI_perBC = length(unique(UMI))), by = .(barcode1, known_BC1, shannon, library)]

UMI_perBC_annot <- merge(UMI_perBC, Associations_raw[barcode_library == BC_LIB, .(barcode1 = BC, oligo, N_assoc = N, BC_type = type)], by = c("barcode1"), all.x = TRUE)
UMI_perBC_annot[is.na(N_assoc), N_assoc := 0]
UMI_perBC_annot[is.na(N_assoc), BC_type := "0_noCRS"]
dir.create(sprintf("%s/data/%s/01b_MPRA_results_counts_BC", MPRA_DIR, ANALYSIS_DIR), showWarnings = FALSE, recursive = TRUE)
fwrite(UMI_perBC_annot, file = sprintf("%s/data/%s/01b_MPRA_results_counts_BC/MPRA_results_counts_BC__%s__v%s.txt.gz", MPRA_DIR, ANALYSIS_DIR, SAMPLE_ID, VERSION_OUT), sep = "\t")
# UMI_perBC_annot=fread(sprintf("%s/data/%s/01b_MPRA_results_counts_BC/MPRA_results_counts_BC__%s__v%s.txt.gz", MPRA_DIR, ANALYSIS_DIR, SAMPLE_ID, VERSION_OUT))

UMI_perBC_stats_BCtype <- UMI_perBC_annot[, .(.N, nUMI = sum(nUMI_perBC)), by = BC_type][order(BC_type), .(BC_type, N, Pct_BC = N / sum(N), nUMI, Pct_UMI = nUMI / sum(nUMI))]
dir.create(sprintf("%s/data/%s/01b_UMI_perBC_stats_BCtype", MPRA_DIR, ANALYSIS_DIR), showWarnings = FALSE, recursive = TRUE)
fwrite(UMI_perBC_stats_BCtype, file = sprintf("%s/data/%s/01b_UMI_perBC_stats_BCtype/UMI_perBC_stats_BCtype__%s__v%s.txt.gz", MPRA_DIR, ANALYSIS_DIR, SAMPLE_ID, VERSION_OUT), sep = "\t")

BC_per_oligo <- UMI_perBC_annot[, .(
  nBC = length(unique(barcode1)),
  nBCover3 = length(unique(barcode1[nUMI_perBC >= 3]))
), by = .(oligo, library)]

## compute percentage of barcodes from the library that are recovered in the sample
BC_recovery <- Associations_Filtered[barcode_library == BC_LIB, .(
  Pct_oligo_recovered = mean(unique(oligo) %chin% UMI_perBC_annot[, oligo]),
  Pct_BC_recovered = mean(unique(BC) %chin% UMI_perBC_annot[, barcode1]),
  Pct_BC_over3UMI = mean(unique(BC) %chin% UMI_perBC_annot[nUMI_perBC >= 3, barcode1]),
  Pct_oligo_over10BC = mean(unique(oligo) %chin% BC_per_oligo[nBC >= 10, oligo]),
  Pct_oligo_over30BC = mean(unique(oligo) %chin% BC_per_oligo[nBC >= 30, oligo]),
  Pct_oligo_over100BC = mean(unique(oligo) %chin% BC_per_oligo[nBC >= 100, oligo]),
  Pct_oligo_over30BCw3UMI = mean(unique(oligo) %chin% BC_per_oligo[nBCover3 >= 30, oligo])
)]
dir.create(sprintf("%s/data/%s/01b_BC_recovery", MPRA_DIR, ANALYSIS_DIR), showWarnings = FALSE, recursive = TRUE)
fwrite(BC_recovery, file = sprintf("%s/data/%s/01b_BC_recovery/BC_recovery__%s__v%s.txt", MPRA_DIR, ANALYSIS_DIR, SAMPLE_ID, VERSION_OUT), sep = "\t")

UMI_perBC_full <- DT_MPRA[sameBC & !UMI %chin% dup_UMIs, .(nUMI_perBC = length(unique(UMI))), by = .(barcode1, known_BC1, shannon, library)]
UMI_perBC_full <- merge(UMI_perBC_full, Associations_raw[barcode_library == BC_LIB, .(barcode1 = BC, oligo, N_assoc = N, BC_type = type)], by = c("barcode1"), all.x = TRUE)
UMI_perBC_full[is.na(N_assoc), N_assoc := 0]
UMI_perBC_full[is.na(N_assoc), BC_type := "0_noCRS"]
dir.create(sprintf("%s/data/%s/01b_MPRA_results_counts_BC_FULL", MPRA_DIR, ANALYSIS_DIR), showWarnings = FALSE, recursive = TRUE)
fwrite(UMI_perBC_full, file = sprintf("%s/data/%s/01b_MPRA_results_counts_BC_FULL/MPRA_results_counts_BC_FULL__%s__v%s.txt.gz", MPRA_DIR, ANALYSIS_DIR, SAMPLE_ID, VERSION_OUT), sep = "\t")

# cat("\nAll done\n")
