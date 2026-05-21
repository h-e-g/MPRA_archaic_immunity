#!/usr/bin/env Rscript

# DONE: implement the possibility to set boot parameter >0 to enable bootsraping of barcodes
# DONE implement the possibility to set perm_oligo parameter >0 to enable the permutation of crsID between oligo (for power calculation)
# DONE: implement annotation of expected difference in activity between TESTED oligo (WHEN boot>0)
# DONE: run emvars with perm 0, boot1 & perm_oligo 1 to assess power as a function of effect size
# DONE: run emvars with perm 1, boot1 & perm_oligo 1 to assess false positives in our power caluclations.
# DONE :find a way to normalize by oligo activity between differnt experiments (so that median and IQR of oligo activity are the same across experiments/conditions)
# idea: compute library size factos in a more robust manner than total depth.
# TODO: simplify.  try to merge different branches (eg. oligo and emVars, or replace oligo with new approach)

MPRA_DIR <- "/pasteur/helix/projects/evo_immuno_pop/MPRA"
ANALYSIS_DIR <- "MPRA_count_exp6_analysis3"

source(sprintf("%s/scripts/misc_plots.R", MPRA_DIR))
source(sprintf("%s/scripts/%s/00_starter.R", MPRA_DIR, ANALYSIS_DIR))
source(sprintf("%s/scripts/%s/00_load_associations.R", MPRA_DIR, ANALYSIS_DIR))


library(MPRAnalyze)

cmd <- commandArgs(trailingOnly = TRUE)
print(cmd)
# parameters defaults
VERSION_IN <- 1
VERSION_OUT <- 1
instruction_file <- sprintf("%s/scripts/%s/02a_MPRA_analyse_instructions_full_sorted.txt", MPRA_DIR, ANALYSIS_DIR)
OUT_DIR <- sprintf("%s/data/%s/02b_runMPRAnalyze/", MPRA_DIR, ANALYSIS_DIR)
sub_sample <- 2000
min_nBC <- 10
Outlier_Z_threshold <- 3
perm_barcode <- 0
perm_allele <- 0
perm_group <- 0
perm_oligo <- 0
boot <- 0
ADJ_LIB <- FALSE
ADJ_REP <- FALSE
# instruction_line <- 1
# chunk_nb <- 1
my_oligo_num='0_0'

for (i in seq_along(cmd)) {
  if (cmd[i] == "--version_in" || cmd[i] == "-vi") {
    VERSION_IN <- cmd[i + 1]
  } # version of input files
  if (cmd[i] == "--version_out" | cmd[i] == "-vo") {
    VERSION_OUT <- cmd[i + 1]
  } # version of output files
  if (cmd[i] == "--sub_sample" || cmd[i] == "-n") {
    sub_sample <- as.numeric(cmd[i + 1])
  } # max number of barcodes to use per oligo
  if (cmd[i] == "--min_nBC" || cmd[i] == "-n") {
    min_nBC <- as.numeric(cmd[i + 1])
  } # min number of barcodes to use per oligo
  if (cmd[i] == "--outlier_threshold" || cmd[i] == "-z") {
    Outlier_Z_threshold <- as.numeric(cmd[i + 1])
  } # max Outlier Z value to include Barcode
  if (cmd[i] == "--instruction_file" || cmd[i] == "-i") {
    instruction_file <- cmd[i + 1]
  } # instruction_file containing instructions for MPRAnalyze
  if (cmd[i] == "--instruction_line" || cmd[i] == "-il") {
    instruction_line <- as.numeric(cmd[i + 1])
  }
  if (cmd[i] == "--chunk_file" || cmd[i] == "-c") {
    chunk_file <- cmd[i + 1]
  }
  if (cmd[i] == "--perm_allele") {
    perm_allele <- as.numeric(cmd[i + 1])
  }
  if (cmd[i] == "--perm_group") {
    perm_group <- as.numeric(cmd[i + 1])
  }
  if (cmd[i] == "--perm_barcode") {
    perm_barcode <- as.numeric(cmd[i + 1])
  }
  if (cmd[i] == "--boot" || cmd[i] == "-b") {
    boot <- as.numeric(cmd[i + 1])
  }
  if (cmd[i] == "--power" || cmd[i] == "--perm_oligo") {
    perm_oligo <- as.numeric(cmd[i + 1])
  }
  if (cmd[i] == "--adj_lib" || cmd[i] == "-l") {
    ADJ_LIB <- as.logical(cmd[i + 1])
  }
  if (cmd[i] == "--adj_rep" || cmd[i] == "-r") {
    ADJ_REP <- as.logical(cmd[i + 1])
  }
  if (cmd[i] == "--output_dir" || cmd[i] == "-o") {
    OUT_DIR <- cmd[i + 1]
  }
}

tic("reading instructions")
# loading instruction file
instruction_DT <- fread(instruction_file, header = TRUE, fill = TRUE)
if (all(is.na(instruction_DT[, group2_samples]))) {
  instruction_DT[, group2_samples := NULL]
  instruction_DT[, group2_samples := ""]
}
ANALYSIS_NAME <- instruction_DT[instruction_line, analysis_name]
ANALYSIS_TYPE <- instruction_DT[instruction_line, analysis_type]
ANALYSIS_SUBTYPE <- instruction_DT[instruction_line, analysis_subtype]
# loading chunk file
test_list <- fread(chunk_file)
chunk_file <- gsub("/+", "/", chunk_file)
chunk_ANALYSIS_TYPE <- gsub(".*/(.*)/(.*)/(.*)/chunk_([0-9]+)_([0-9]+).txt", "\\1", chunk_file)
chunk_ANALYSIS_SUBTYPE <- gsub(".*/(.*)/(.*)/(.*)/chunk_([0-9]+)_([0-9]+).txt", "\\2", chunk_file)
chunk_ANALYSIS_NAME <- gsub(".*/(.*)/(.*)/(.*)/chunk_([0-9]+)_([0-9]+).txt", "\\3", chunk_file)
chunk_instruction_line <- gsub(".*/(.*)/(.*)/(.*)/chunk_([0-9]+)_([0-9]+).txt", "\\4", chunk_file)
chunk_nb <- gsub(".*/(.*)/(.*)/(.*)/chunk_([0-9]+)_([0-9]+).txt", "\\5", chunk_file)

STOP <- FALSE
if (ANALYSIS_NAME != chunk_ANALYSIS_NAME) {
  STOP <- TRUE
  cat(sprintf("\n%s does not match line %s of %s: wrong analysis_name\n", chunk_file, instruction_line, instruction_file))
}
if (ANALYSIS_TYPE != chunk_ANALYSIS_TYPE) {
  STOP <- TRUE
  cat(sprintf("\n%s does not match, line %s of %s: wrong analysis_type\n", chunk_file, instruction_line, instruction_file))
}
if (ANALYSIS_SUBTYPE != chunk_ANALYSIS_SUBTYPE) {
  STOP <- TRUE
  cat(sprintf("\n%s does not match, line %s of %s: wrong analysis_subtype\n", chunk_file, instruction_line, instruction_file))
}
# if (instruction_line != chunk_instruction_line) {
#   STOP <- TRUE
#   cat(sprintf("\n%s does not match, line %s of %s: wrong instruction_line\n", chunk_file, instruction_line, instruction_file))
# }
if (STOP) {
  stop("Chunk file error. job cancelled\n")
}

# extracting analysis parameters
cat("\n", instruction_line, ":", ANALYSIS_NAME, "-", ANALYSIS_TYPE, ":")

# extracting group 1 samples
group_1_samples <- instruction_DT[instruction_line, group1_samples]
if (group_1_samples != "") {
  group_1_samples <- as.character(str_split(group_1_samples, ",", simplify = TRUE))
  if (any(!group_1_samples %in% Indexes[, unique(SETUP_ID)])) {
    stop(sprintf("samples %s not found in Indexes", paste(group_1_samples[!group_1_samples %in% Indexes[, unique(SETUP_ID)]], collapse = ", ")))
  }
} else {
  stop("no samples specified")
}

# extracting group 2 samples
group_2_samples <- instruction_DT[instruction_line, group2_samples]
if (group_2_samples != "") {
  group_2_samples <- str_split(group_2_samples, ",", simplify = TRUE)
  if (any(!group_2_samples %in% Indexes[, unique(SETUP_ID)])) {
    stop(sprintf("samples %s not found in Indexes", paste(group_2_samples[!group_2_samples %in% Indexes[, unique(SETUP_ID)]], collapse = ", ")))
  }
} else {
  group_2_samples <- c()
  cat("no samples specified for group 2\n")
}
toc()

tic("loading oligo & SNP annotations")
# load annotation of oligos
oligo_source <- fread(sprintf("%s/data/%s/00_oligo_annot_v2.txt", MPRA_DIR, ANALYSIS_DIR))
# define order of priority and assign each oligo to a single annotation
oligo_source <- oligo_source[order(factor(type, c(
  "archaic", "COVID", "Expecto", "MPRA", "eQTLs", "eQTL_repeat", "Purged", "randomSNP",
  "Promoter_Lung", "Promoter_Tcell", "Promoter_Mono_NS", "Promoter_Mono_STIM", "scrambled"
))), ]
oligo_source <- oligo_source[!duplicated(oligo), ]
# load annotations of SNPs
SNP_annot <- fread(sprintf("%s/data/%s/00_SNP_annot_v1.txt", MPRA_DIR, ANALYSIS_DIR))
toc()

tic("loading data")
# load UMI per oligo
#  - for each sample (ID_rep) with material (RNA or DNA)
load_time <- c()
UMI_perBC_annot <- list()
for (mySETUP_ID in unique(c(group_1_samples, group_2_samples))) {
  cat("reading UMI counts for", mySETUP_ID, "\n")
  # read UMI counts
  UMI_perBC_annot[[mySETUP_ID]] <- fread(sprintf("%s/data/%s/01c_MPRA_results_BC_activity/01c_MPRA_results_BC_activity__%s__v%s.txt.gz", MPRA_DIR, ANALYSIS_DIR, mySETUP_ID, VERSION_IN))
}

# merge data from RNA and DNA samples & add library info
UMI_perBC_annot <- rbindlist(UMI_perBC_annot, use.names = TRUE)
# annotate BCs
UMI_per_BC_annot_oligo <- merge(UMI_perBC_annot, unique(oligo_source)[oligo != "", .(oligo, source, type, shift, SNP, strand, crsID, allele, GC)], by = "oligo", allow.cartesian = TRUE)
# filter out outlier BCs
UMI_per_BC_annot_oligo <- UMI_per_BC_annot_oligo[abs(outlier_Z) < Outlier_Z_threshold, ]
# subsampling Barcodes per oligo (ie. CRS + allele)
set.seed(sub_sample)
UMI_per_BC_annot_oligo <- UMI_per_BC_annot_oligo[, .SD[sample(seq_len(.N), min(.N, sub_sample)), ], by = .(oligo, SETUP_ID, allele)]
# UMI_per_BC_annot_oligo <- UMI_per_BC_annot_oligo[,.(oligo,barcode_lib,barcode1,celltype,replicate,SETUP_ID,condition,material,nUMI_perBC,source,type,shift,SNP,strand,crsID,allele,GC)]
# convert to long format
UMI_perBC_annot <- melt(UMI_per_BC_annot_oligo, measure.vars = c("RNA", "DNA"), value.name = c("nUMI_perBC"), variable.name = "material")
toc()

if (boot > 0) {
  tic("boostraping")
  # bootstrap barcodes for each oligo & library
  set.seed(boot)
  UMI_perBC_annot[, nUMI_perBC := 1 + rbinom(.N,
    size = (sum(nUMI_perBC) - .N) * boot,
    prob = nUMI_perBC / sum(nUMI_perBC)
  ), by = .(SETUP_ID, material)]
  toc()
}

# create random combinations of BC to define NULL oligo activity
if (perm_barcode > 0) {
  tic("permuting barcode between oligos (for NULL oligo activity)")
  set.seed(perm_barcode)
  UMI_per_BC_annot_oligo_perm <- UMI_per_BC_annot_oligo[seq_len(.N), ]
  setnames(UMI_per_BC_annot_oligo_perm, "oligo", "oligo_true")
  UMI_per_BC_annot_oligo_perm[, c("source", "type", "shift", "SNP", "strand", "crsID", "allele", "GC") := NULL]
  Associations_Filtered_perm <- Associations_Filtered[seq_len(.N), ]
  Associations_Filtered_perm[, oligo := sample(oligo), by = barcode_library]
  UMI_per_BC_annot_oligo_perm <- merge(UMI_per_BC_annot_oligo_perm, Associations_Filtered_perm[, .(barcode1 = BC, barcode_lib = barcode_library, oligo)], by = c("barcode_lib", "barcode1"))
  UMI_per_BC_annot_oligo_perm <- merge(UMI_per_BC_annot_oligo_perm, unique(oligo_source)[oligo != "", .(oligo, source, type, shift, SNP, strand, crsID, allele, GC)], by = "oligo")
  UMI_per_BC_annot_oligo_perm <- UMI_per_BC_annot_oligo_perm[, mget(colnames(UMI_per_BC_annot_oligo))]
  UMI_per_BC_annot_oligo_perm <- merge(UMI_per_BC_annot_oligo_perm, technical_covariates[, .(SETUP_ID, plasmid_library, lentivirus_preparation)], by = "SETUP_ID")
  toc()
}

if (perm_oligo > 0) {
  tic("assigning pairs of oligos with differential activity to each crsID for power analysis")
  # assign random pair of oligos to each SNP (to allow for power analysis)
  UMI_per_BC_annot_oligo[, crsID := NULL]
  crsID_perm <- unique(oligo_source)[oligo != "", .(oligo, crsID)]
  nBC_mean <- Associations_Filtered[, .N, by = .(oligo, barcode_library)][, .(nBC = mean(N * ifelse(barcode_library %in% c("LIB1", "LIB2"), 2, 1))), by = oligo]
  crsID_perm <- merge(crsID_perm, nBC_mean, by = "oligo")
  set.seed(perm_oligo)
  crsID_perm[, crsID := sample(crsID), by = cut(nBC, c(0, 30, 100, 200, 500, 1000, 2000))]
  # keep only crsID with 2 oligo exactly
  nOligo_by_crsID <- crsID_perm[, .N, by = crsID]
  crsID_over2oligo <- nOligo_by_crsID[N > 2, crsID]
  UMI_per_BC_annot_oligo <- merge(UMI_per_BC_annot_oligo, crsID_perm, by = "oligo")
  UMI_per_BC_annot_oligo[, nBC := NULL]
  toc()
}

UMI_per_BC_annot_oligo <- merge(UMI_per_BC_annot_oligo, technical_covariates[, .(SETUP_ID, plasmid_library, lentivirus_preparation)], by = "SETUP_ID")

tic("computing depth per library")
####### compute total UMI per oligo
DEPTH <- UMI_perBC_annot[, .(total_UMI = sum(nUMI_perBC)), by = .(celltype, material, replicate, SETUP_ID, condition)]
BC_meanDEPTH <- merge(UMI_perBC_annot, DEPTH, by = c("celltype", "material", "replicate", "SETUP_ID", "condition"))
BC_meanDEPTH[, mean_logUPM := mean(log(1 + (nUMI_perBC / total_UMI * 1e6))), by = .(barcode1, barcode_lib, oligo, material)]

## adjust depth per library to ensure similar RNA/DNA ratio across libraries (similar to DESeq)
## DE seq finds the ratio of each read count to the geometric mean of all read counts for that gene across all samples (the denominator serving as a pseudo-reference sample [24]). The median of these ratios for a sample, called the 'size factor', is used to scale that sample.
## here it would be : size factor = median( across all BCs of a library of UMI_per_million / exp(mean(log(UMI_per_million) across all libraries))

Size_factors <- BC_meanDEPTH[, .(size_factor = median(nUMI_perBC / total_UMI * 1e6 / (exp(mean_logUPM) - 1), na.rm = TRUE)), by = .(celltype, material, replicate, SETUP_ID, condition)]

DEPTH <- merge(DEPTH, Size_factors, by = c("celltype", "material", "replicate", "SETUP_ID", "condition"))
DEPTH[, total_UMI_sizeCorrected := total_UMI * size_factor]
#
DEPTH <- dcast(DEPTH, celltype + condition + replicate + SETUP_ID ~ material, value.var = "total_UMI_sizeCorrected")
toc()


tic("initializing")
OUTS_DIR <- sprintf("%s/%s/%s/%s/perm_%s_%s_%s", OUT_DIR, ANALYSIS_TYPE, ANALYSIS_SUBTYPE, ANALYSIS_NAME, perm_barcode, perm_group, perm_allele)
if (boot > 0 || perm_oligo > 0) {
  OUTS_DIR <- sprintf("%s/power%s_boot%s", OUTS_DIR, perm_oligo, boot)
}
dir.create(OUTS_DIR, recursive = TRUE, showWarnings = FALSE)

if (chunk_nb == 1 & perm_oligo > 0) {
  # for first chunk save the crsID-oligo correspondance
  fwrite(crsID_perm, file = sprintf("%s/crsID_perm.txt", OUTS_DIR), sep = "\t")
}
toc()

############ choose oligo to test
counter <- 0
alphas_and_pvalues <- list()
for (test in test_list$test) {
  counter <- counter + 1
  my_oligo_num <- paste(chunk_nb, counter, sep = "_")
  # for (instruction_line in seq_len(instruction_DT[, .N])) {
  # does oligo number exist for this comparison

  if (grepl("EMVAR", toupper(ANALYSIS_TYPE))) {
    analysis_unit <- "SNPs"
    # get position of target SNP
    tic(paste("shaping data for crsID number", counter, ":", test))
    if (perm_oligo > 0 && test %chin% crsID_over2oligo) {
      cat("over 2 oligo for this position, skipping\n")
      next
    }
    SNP_data <- UMI_per_BC_annot_oligo[crsID == test, ]
    # SNP_data <- SNP_data[SETUP_ID %chin% c(group_1_samples, group_2_samples), .(nUMI_perBC = sum(nUMI_perBC)), by = .(oligo, barcode1, barcode_lib, celltype, condition, replicate, SETUP_ID, material, type, crsID, shift, strand, allele, GC)]
    SNP_data <- SNP_data[SETUP_ID %chin% c(group_1_samples, group_2_samples), .(DNA = sum(DNA), RNA = sum(RNA)), by = .(oligo, barcode1, barcode_lib, celltype, condition, replicate, SETUP_ID, type, crsID, shift, strand, allele, GC, plasmid_library, lentivirus_preparation)]
    if (SNP_data[, length(unique(barcode1))] == 0) {
      cat("no barcode in the current celltype for this position, skipping\n")
      next
    }
    # SNP_data <- dcast(SNP_data, oligo + barcode1 + barcode_lib + celltype + condition + replicate + SETUP_ID + type + strand + crsID + shift + allele + GC ~ factor(material,c('DNA','RNA')), value.var = "nUMI_perBC", fill = 0, drop= c(TRUE, FALSE))
    SNP_data <- merge(SNP_data, DEPTH, by = c("celltype", "condition", "replicate", "SETUP_ID"), suffix = c("", "_depth"))

    # add info on library and remove non-unique/low confidence oligo-barcodes associations
    # SNP_data <- merge(SNP_data, Associations_Filtered[type == "1a_assoc_unique", .(barcode1 = BC, oligo, barcode_lib = barcode_library, BC_type = type)], by = c("barcode1", "oligo", "barcode_lib"))

    if (perm_oligo == 0 & test %in% SNP_annot[,crsID]) {
      # get numeric allele to harmonize easily across SNPs
      allele_annot <- melt(SNP_annot[crsID == test, .(crsID, allele.1, allele.2)], id.vars = "crsID")
      setnames(allele_annot, c("variable", "value"), c("allele", "allele.char"))
      allele_annot <- unique(allele_annot) # collapse duplicated annotation (eg. random SNP & archaic SNP)
      allele_annot[, allele.num := as.numeric(gsub("allele.", "", allele))]
      allele_annot[, allele := gsub("allele.", "a", allele)]
      SNP_data <- merge(allele_annot[, .(allele.char, allele)], SNP_data, by.x = "allele.char", by.y = "allele")
      oligo_list <- SNP_data[order(allele), unique(oligo)]
    } else {
      oligo_list <- SNP_data[, unique(oligo)]
      SNP_data[, allele := setNames(c("a1", "a2"), oligo_list)[oligo]]
    }
    # too_many_BC <- SNP_data[, .(N_BC = length(unique(paste(barcode1, barcode_lib)))), by = .(oligo, SETUP_ID, allele)][, N_BC] > sub_sample
    # if (sub_sample > 0 && any(too_many_BC)) {
    #   SNP_data <- SNP_data[, .SD[sample(seq_len(.N), min(.N, sub_sample)), ], by = .(oligo, SETUP_ID, allele)]
    #   # SNP_data <- SNP_data[order(SETUP_ID,oligo, -DNA),]
    #   # SNP_data <- SNP_data[,head(.SD, min(.N, sub_sample/2)), by = .(oligo,SETUP_ID)]
    # }
    SNP_data[, group := ifelse(SETUP_ID %chin% group_2_samples, "2", "1")]
    nBCs <- SNP_data[, .N, by = .(group, allele)]
    SNP_stats <- SNP_data[, .(crsID = test, oligo_num = my_oligo_num, ANALYSIS_NAME, ANALYSIS_SUBTYPE, ANALYSIS_TYPE, nBCs = .N, nUMI_DNA = sum(DNA), nUMI_RNA = sum(RNA), nUMI_DNA_perBC = mean(DNA), nUMI_RNA_perBC = mean(RNA)), keyby = .(group, allele)]
    SNP_stats <- dcast(SNP_stats, crsID + oligo_num + ANALYSIS_NAME + ANALYSIS_SUBTYPE + ANALYSIS_TYPE ~ factor(paste0("g", group, allele), levels = c("g1a1", "g1a2", "g2a1", "g2a2")), fill = 0, value.var = c("nBCs", "nUMI_DNA", "nUMI_RNA", "nUMI_DNA_perBC", "nUMI_RNA_perBC"), drop = FALSE)

    nLIB <- SNP_data[, .(nLIB = length(unique(paste(plasmid_library, lentivirus_preparation)))), by = .(group, allele)][, max(nLIB)]
    nRep <- SNP_data[, .(nRep = length(unique(SETUP_ID))), by = .(group, allele)][, max(nRep)]
    # SNP_print <- SNP_data[, .(.N,
    #   nBC_RNA = sum(RNA > 0), nUMI_RNA = sum(RNA), RNA_depth = unique(RNA_depth),
    #   nBC_DNA = sum(DNA > 0), nUMI_DNA = sum(DNA), DNA_depth = unique(DNA_depth)
    # ), keyby = .(celltype, condition, barcode_lib, allele, group, replicate, SETUP_ID, oligo, crsID)]
    # SNP_print[, OR := (nUMI_RNA / nBC_RNA) / (nUMI_DNA / nBC_DNA) / (RNA_depth / DNA_depth)]
    # print(SNP_print)
    # fwrite(SNP_print, file = sprintf("%s/stats_v%s_SNP%s_max%sBC.txt", OUTS_DIR, VERSION_OUT, my_oligo_num, sub_sample), sep = "\t")

    if (perm_allele > 0) {
      set.seed(perm_allele)
      SNP_data[, allele := sample(allele)]
    }
    if (SNP_data[allele == "a1", length(unique(barcode1)) == 0]) {
      cat("no barcode for allele 1, skipping\n", sep = "")
      next
    }
    if (SNP_data[allele == "a2", length(unique(barcode1)) == 0]) {
      cat("no barcode for allele 2, skipping\n", sep = "")
      next
    }
    if (length(group_2_samples) > 0) {
      if (perm_group > 0) {
        set.seed(perm_group)
        SNP_data[, group := sample(group)]
      }
      rownames(SNP_data) <- SNP_data[, paste(barcode1, barcode_lib, celltype, condition, replicate, group, crsID, shift, allele, plasmid_library, lentivirus_preparation, sep = "_")]
      # prepare barcode annotation for MPRAAnalyze
      my.colAnnot <- SNP_data[, .(barcode1, barcode_lib, celltype, condition, replicate, strand, group, crsID, shift, allele, oligo, plasmid_library, lentivirus_preparation)]
      rownames(my.colAnnot) <- rownames(SNP_data)
    } else {
      rownames(SNP_data) <- SNP_data[, paste(barcode1, barcode_lib, celltype, condition, replicate, crsID, shift, allele, plasmid_library, lentivirus_preparation, sep = "_")]
      # prepare barcode annotation for MPRAAnalyze
      my.colAnnot <- SNP_data[, .(barcode1, barcode_lib, celltype, condition, replicate, strand, crsID, shift, allele, oligo, plasmid_library, lentivirus_preparation)]
      rownames(my.colAnnot) <- rownames(SNP_data)
    }
    # prepare RNA counts for MPRAAnalyze
    rna.counts <- matrix(SNP_data[, RNA], 1, SNP_data[, .N])
    dimnames(rna.counts) <- list(test, rownames(SNP_data))
    # prepare DNA counts for MPRAAnalyze
    dna.counts <- matrix(SNP_data[, DNA], 1, SNP_data[, .N])
    dimnames(dna.counts) <- list(test, rownames(SNP_data))

    obj <- MpraObject(
      dnaCounts = dna.counts, rnaCounts = rna.counts,
      dnaAnnot = my.colAnnot, rnaAnnot = my.colAnnot
    )
    obj <- setDepthFactors(obj, dnaDepth = SNP_data$DNA_depth, rnaDepth = SNP_data$RNA_depth)
    toc()

    if (length(group_2_samples) == 0) {
      tic("assessing SNP effect only")
      if (nLIB > 1 && ADJ_LIB) {
        obj_quantif <- analyzeQuantification(obj = obj, dnaDesign = ~ 0 + allele + paste(plasmid_library, lentivirus_preparation), rnaDesign = ~ 0 + allele + paste(plasmid_library, lentivirus_preparation))
        if (nRep > 1 && ADJ_REP) {
          obj_test <- analyzeComparative(obj = obj, dnaDesign = ~ allele + paste(replicate, plasmid_library, lentivirus_preparation), rnaDesign = ~ 1 + allele + paste(replicate, plasmid_library, lentivirus_preparation), reducedDesign = ~ 1 + paste(replicate, plasmid_library, lentivirus_preparation), fit.se = TRUE)
        } else {
          obj_test <- analyzeComparative(obj = obj, dnaDesign = ~ allele + paste(plasmid_library, lentivirus_preparation), rnaDesign = ~ 1 + allele + paste(plasmid_library, lentivirus_preparation), reducedDesign = ~ 1 + paste(plasmid_library, lentivirus_preparation), fit.se = TRUE)
        }
      } else {
        obj_quantif <- analyzeQuantification(obj = obj, dnaDesign = ~ 0 + allele, rnaDesign = ~ 0 + allele)
        if (nRep > 1 && ADJ_REP) {
          obj_test <- analyzeComparative(obj = obj, dnaDesign = ~ replicate + allele, rnaDesign = ~ 1 + replicate + allele, reducedDesign = ~ 1 + replicate, fit.se = TRUE) # with replicate
          # obj_test <- analyzeComparative(obj = obj, dnaDesign = ~ replicate + barcode1 + allele, rnaDesign = ~ 1 + replicate + allele, reducedDesign = ~1 +replicate, fit.se = TRUE) # with barcode in DNA, no change in pvalue
        } else {
          obj_test <- analyzeComparative(obj = obj, dnaDesign = ~allele, rnaDesign = ~ 1 + allele, reducedDesign = ~1, fit.se = TRUE) # simple, no replicate
        }
      }
      toc()
      alpha_quantif <- getAlpha(obj_quantif, by.factor = "all")[, c("allelea1", "allelea2")]
      colnames(alpha_quantif) <- c("a1", "a2")
      LRT <- as.data.table(testLrt(obj_test))[, .(logFC, statistic, df.test, pval.LRT = pval)]
      logFC.se <- obj_test@modelFits$r.se[1 + sum(apply(obj_test@designs@rnaFull, 2, sum) != 0)]
      if (is.null(logFC.se)) {
        logFC.se <- NA
      }
      # get # of UMIs and barcodes
      alphas_and_pvalues[[test]] <- cbind(data.table(crsID = test, oligo1 = oligo_list[1], oligo2 = oligo_list[2], GC = SNP_data[, mean(GC)], alpha_quantif, LRT, logFC.se, SNP_stats[, -"crsID"]))
      # SNP_data[,.(sum(DNA>0),sum(RNA>0)),by=.(oligo,celltype,condition)]
      # fwrite(alphas_and_pvalues, file = sprintf("%s/alphas_and_pvalues_v%s_SNP%s_max%sBC_Z%s.txt", OUTS_DIR, VERSION_OUT, my_oligo_num, sub_sample, Outlier_Z_threshold), sep = "\t")
      cat(sprintf("\ncrsID number %s processed (%s)\n", counter, test))
    }
    if (length(group_2_samples) > 0) {
      if (SNP_data[allele == "a1", length(unique(group)) < 2]) {
        cat("no barcode for allele 1 in group ", setdiff(1:2, SNP_data[allele == "a1", unique(group)]), ", skipping\n", sep = "")
        next
      }
      if (SNP_data[allele == "a2", length(unique(group)) < 2]) {
        cat("no barcode for allele 2 in group ", setdiff(1:2, SNP_data[allele == "a2", unique(group)]), ", skipping\n", sep = "")
        next
      }
      if (nLIB > 1 && ADJ_LIB) {
        tic("assessing difference in SNP effect between group")
        obj_quantif <- analyzeQuantification(obj = obj, dnaDesign = ~ 0 + paste(allele, group, sep = "_") + paste(plasmid_library, lentivirus_preparation), rnaDesign = ~ 0 + paste(allele, group, sep = "_") + paste(plasmid_library, lentivirus_preparation))
        obj_diffSNP <- analyzeComparative(obj = obj, dnaDesign = ~ allele * group + paste(plasmid_library, lentivirus_preparation), rnaDesign = ~ 1 + allele * group + paste(plasmid_library, lentivirus_preparation), reducedDesign = ~ 1 + allele + group + paste(plasmid_library, lentivirus_preparation), fit.se = TRUE)
        toc()
        tic("assessing SNP effect jointly across groups")
        obj_joinSNP <- analyzeComparative(obj = obj, dnaDesign = ~ allele * group + paste(plasmid_library, lentivirus_preparation), rnaDesign = ~ 1 + allele * group + paste(plasmid_library, lentivirus_preparation), reducedDesign = ~ 1 + group + paste(plasmid_library, lentivirus_preparation), fit.se = TRUE)
        toc()
      } else {
        tic("assessing difference in SNP effect between group")
        obj_quantif <- analyzeQuantification(obj = obj, dnaDesign = ~ 0 + paste(allele, group, sep = "_"), rnaDesign = ~ 0 + paste(allele, group, sep = "_"))
        obj_diffSNP <- analyzeComparative(obj = obj, dnaDesign = ~ allele * group, rnaDesign = ~ 1 + allele * group, reducedDesign = ~ 1 + allele + group, fit.se = TRUE)
        toc()
        tic("assessing SNP effect jointly across groups")
        obj_joinSNP <- analyzeComparative(obj = obj, dnaDesign = ~ allele * group, rnaDesign = ~ 1 + allele * group, reducedDesign = ~ 1 + group, fit.se = TRUE)
        toc()
      }
      alpha_quantif <- getAlpha(obj_quantif, by.factor = "all")[, c(1, 3, 2, 4)] # c('group1_a1','group1_a2','group2_a1','group2_a2')
      colnames(alpha_quantif) <- c("group1_a1", "group1_a2", "group2_a1", "group2_a2")
      delta_logFC.se <- obj_diffSNP@modelFits$r.se[1 + sum(apply(obj_diffSNP@designs@rnaFull, 2, sum) != 0)]
      if (is.null(delta_logFC.se)) {
        delta_logFC.se <- NA
      }
      interaction_results <- cbind(as.data.table(testLrt(obj_diffSNP))[, .(delta_logFC = logFC, statistic_int = statistic, df.test_int = df.test, pval_int = pval)], delta_logFC.se = delta_logFC.se, Pval_joint = testLrt(obj_joinSNP)$pval)

      alphas_and_pvalues[[test]] <- data.table(crsID = test, oligo1 = oligo_list[1], oligo2 = oligo_list[2], GC = SNP_data[, mean(GC)], alpha_quantif, interaction_results, SNP_stats[, -"crsID"])
      # fwrite(alphas_and_pvalues, file = sprintf("%s/alphas_and_pvalues_v%s_SNP%s_max%sBC_Z%s.txt", OUTS_DIR, VERSION_OUT, my_oligo_num, sub_sample, Outlier_Z_threshold), sep = "\t")
      cat(sprintf("\ncrsID number %s processed (%s)\n", counter, test))
    }
  }



  if (grepl("OLIGO", toupper(ANALYSIS_TYPE))) {
    analysis_unit <- "oligos"
    tic("shaping data")
    if (perm_barcode == 0 || length(group_2_samples) > 0) {
      # use actual BCs for each oligo
      oligo_data <- UMI_per_BC_annot_oligo[oligo == test, ]
    } else {
      # to estimate NULL oligo activity, use random combination of BCs
      oligo_data <- UMI_per_BC_annot_oligo_perm[oligo == test, ]
    }
    # oligo_data <- oligo_data[SETUP_ID %chin% c(group_1_samples, group_2_samples), .(nUMI_perBC = sum(nUMI_perBC)), by = .(oligo, barcode1, barcode_lib, celltype, condition, replicate, SETUP_ID, type, crsID, shift, strand, allele, GC)]
    # oligo_data <- dcast(oligo_data, oligo + barcode1 + barcode_lib + celltype + condition + replicate + SETUP_ID + type + strand + crsID + shift + allele + GC ~ factor(material,c('DNA','RNA')), value.var = "nUMI_perBC", fill = 0, drop=c(TRUE,FALSE))
    oligo_data <- merge(oligo_data, DEPTH, by = c("celltype", "condition", "SETUP_ID"), suffix = c("", "_depth"))
    #    oligo_data <- merge(oligo_data, DEPTH_0, by = c("celltype", "condition", "SETUP_ID"), suffix = c("", "_depth"))

    oligo_data[, group := ifelse(SETUP_ID %chin% group_2_samples, "2", "1")]

    if (oligo_data[, .N] < min_nBC) {
      print(sprintf("skipping %s for oligo %s, only %s observation(s)", ANALYSIS_NAME, test, oligo_data[, .N]))
      next
    }
    # if(oligo_data[DNA>0,.N]==0){
    #   print(sprintf("skipping %s for oligo %s, no observation(s) with DNA",ANALYSIS_NAME, test))
    #     next
    # }
    # if(oligo_data[RNA>0,.N]==0){
    #   print(sprintf("skipping %s for oligo %s, no observation(s) with RNA",ANALYSIS_NAME, test))
    #     next
    # }
    # too_many_BC <- oligo_data[, .(N_BC = length(paste(unique(barcode1, barcode_lib)))), by = SETUP_ID][, N_BC] > sub_sample
    # if (sub_sample > 0 && any(too_many_BC)) {
    #   if (oligo_data[, .N] > sub_sample && sub_sample > 0) {
    #     oligo_data <- oligo_data[, .SD[sample(seq_len(.N), min(.N, sub_sample))], by = SETUP_ID]
    #   }
    # }
    if (length(group_2_samples) > 0) {
      if (perm_group > 0) {
        set.seed(perm_group)
        oligo_data[, group := sample(group)]
      }
      if (oligo_data[, length(unique(group))] < 2) {
        print(sprintf("skipping %s for oligo %s, one group has no barcodes", ANALYSIS_NAME, test, oligo_data[, .N]))
        next
      } else {
        if (oligo_data[, .N, by = group][, any(N < min_nBC)]) {
          print(sprintf("skipping %s for oligo %s, at least one group has only %s observation(s)", ANALYSIS_NAME, test, oligo_data[, .N, by = group][, min(N)]))
          next
        }
        # if(oligo_data[,.(N=sum(DNA>0)),by=group][,any(N==0)]){
        #    print(sprintf("skipping %s for oligo %s,  at least one group has no observation(s) with DNA",ANALYSIS_NAME, test))
        #   next
        # }
        # if(oligo_data[,.(N=sum(RNA>0)),by=group][,any(N==0)]){
        #   print(sprintf("skipping %s for oligo %s, at least one group has no observation(s) with RNA",ANALYSIS_NAME, test))
        #   next
        #    }
      }
      rownames(oligo_data) <- oligo_data[, paste(barcode1, barcode_lib, celltype, condition, replicate, group, crsID, shift, allele, plasmid_library, lentivirus_preparation, sep = "_")]
      # prepare barcode annotation for MPRAAnalyze
      my.colAnnot <- oligo_data[, .(barcode1, barcode_lib, celltype, strand, condition, replicate, group, crsID, shift, allele, oligo, SETUP_ID, plasmid_library, lentivirus_preparation)]
    } else {
      rownames(oligo_data) <- oligo_data[, paste(barcode1, barcode_lib, celltype, condition, replicate, crsID, shift, allele, plasmid_library, lentivirus_preparation, sep = "_")]
      # prepare barcode annotation for MPRAAnalyze
      my.colAnnot <- oligo_data[, .(barcode1, barcode_lib, celltype, strand, condition, replicate, crsID, shift, allele, oligo, SETUP_ID, plasmid_library, lentivirus_preparation)]
    }
    rownames(my.colAnnot) <- rownames(oligo_data)

    # prepare RNA counts for MPRAAnalyze
    rna.counts <- matrix(oligo_data[, RNA], 1, oligo_data[, .N])
    dimnames(rna.counts) <- list(test, rownames(oligo_data))
    # prepare DNA counts for MPRAAnalyze
    dna.counts <- matrix(oligo_data[, DNA], 1, oligo_data[, .N])
    dimnames(dna.counts) <- list(test, rownames(oligo_data))

    obj <- MpraObject(
      dnaCounts = dna.counts, rnaCounts = rna.counts,
      dnaAnnot = my.colAnnot, rnaAnnot = my.colAnnot
    )
    obj <- setDepthFactors(obj, dnaDepth = oligo_data$DNA_depth, rnaDepth = oligo_data$RNA_depth)
    toc()

    oligo_stats <- oligo_data[, .(oligo = test, oligo_num = my_oligo_num, ANALYSIS_NAME, ANALYSIS_SUBTYPE, ANALYSIS_TYPE, nBC = .N, nUMI_DNA = sum(DNA), nUMI_RNA = sum(RNA), nUMI_DNA_perBC = mean(DNA), nUMI_RNA_perBC = mean(RNA)), by = group]
    oligo_stats <- dcast(oligo_stats, oligo + oligo_num + ANALYSIS_NAME + ANALYSIS_SUBTYPE + ANALYSIS_TYPE ~ factor(paste0("g", group), levels = c("g1", "g2")), fill = 0, value.var = c("nBC", "nUMI_DNA", "nUMI_RNA", "nUMI_DNA_perBC", "nUMI_RNA_perBC"), drop = FALSE)

    nLIB <- oligo_data[, .(nLIB = length(unique(plasmid_library, lentivirus_preparation))), by = .(group)][, max(nLIB)]
    nRep <- oligo_data[, .(nRep = length(unique(replicate))), by = .(group)][, max(nRep)]
    # oligo_print <- oligo_data[, .(.N,
    #   nBC_RNA = sum(RNA > 0), nUMI_RNA = sum(RNA), RNA_depth = unique(RNA_depth),
    #   nBC_DNA = sum(DNA > 0), nUMI_DNA = sum(DNA), DNA_depth = unique(DNA_depth)
    # ), by = .(celltype, condition, barcode_lib, replicate, group, SETUP_ID)]
    # oligo_print[, OR := (nUMI_RNA / nBC_RNA) / (nUMI_DNA / nBC_DNA) / (RNA_depth / DNA_depth)]
    # oligo_print[, OR_adj := (nUMI_RNA /nUMI_DNA) / (RNA_depth / DNA_depth)]
    # oligo_print[, oligo := test]
    # print(oligo_print)
    # fwrite(oligo_print, file = sprintf("%s/stats_v%s_oligo%s_max%sBC.txt", OUTS_DIR, VERSION_OUT, my_oligo_num, sub_sample), sep = "\t")
    if (length(group_2_samples) == 0) {
      tic("assessing oligo activity only")
      if (nLIB > 1 && ADJ_LIB) {
        obj_quantif <- analyzeQuantification(obj = obj, dnaDesign = ~ 1 + paste(plasmid_library, lentivirus_preparation), rnaDesign = ~ 1 + paste(plasmid_library, lentivirus_preparation))
        obj_test <- analyzeComparative(obj = obj, dnaDesign = ~ 1 + paste(plasmid_library, lentivirus_preparation), rnaDesign = ~ 1 + paste(plasmid_library, lentivirus_preparation), reducedDesign = ~1 + paste(plasmid_library, lentivirus_preparation), fit.se = TRUE)
      } else {
        obj_quantif <- analyzeQuantification(obj = obj, dnaDesign = ~1, rnaDesign = ~1)
        obj_test <- analyzeComparative(obj = obj, dnaDesign = ~1, rnaDesign = ~1, reducedDesign = ~1, fit.se = TRUE)
      }
      toc()
      alpha_quantif <- getAlpha(obj_quantif)
      alpha.se <- obj_test@modelFits$r.se[1 + sum(apply(obj_test@designs@rnaFull, 2, sum) != 0)]
      if (is.null(alpha.se)) {
        alpha.se <- NA
      }

      alphas_and_pvalues[[test]] <- data.table(oligo = test, GC = oligo_data[, mean(GC)], alpha_quantif, alpha.se, oligo_stats[, -"oligo"])
      # fwrite(alphas_and_pvalues, file = sprintf("%s/alphas_and_pvalues_v%s_oligo%s_max%sBC_Z%s.txt", OUTS_DIR, VERSION_OUT, my_oligo_num, sub_sample, Outlier_Z_threshold), sep = "\t")
      cat(sprintf("\n oligo number %s processed (%s)\n", counter, test))
    }
    if (length(group_2_samples) > 0) {
      tic("assessing oligo activity difference between group")
      if (nLIB > 1 && ADJ_LIB) {
        obj_quantif <- analyzeQuantification(obj = obj, dnaDesign = ~ 0 + group + paste(plasmid_library, lentivirus_preparation), rnaDesign = ~ 0 + group + paste(plasmid_library, lentivirus_preparation))
      } else {
        obj_quantif <- analyzeQuantification(obj = obj, dnaDesign = ~ 0 + group, rnaDesign = ~ 0 + group)
        # obj_quantif <- analyzeQuantification(obj = obj, dnaDesign = ~ 0 + replicate + group, rnaDesign = ~ 0 + replicate+ group)
      }

      alpha_quantif <- getAlpha(obj_quantif, by.factor = "all")[, c("group1", "group2")]
      if (nLIB > 1 && ADJ_LIB) {
        obj_diffgroup <- analyzeComparative(obj = obj, dnaDesign = ~ group + paste(plasmid_library, lentivirus_preparation), rnaDesign = ~ 1 + group + paste(plasmid_library, lentivirus_preparation), reducedDesign = ~ 1 + paste(plasmid_library, lentivirus_preparation), fit.se = TRUE)
      } else {
        obj_diffgroup <- analyzeComparative(obj = obj, dnaDesign = ~group, rnaDesign = ~ 1 + group, reducedDesign = ~1, fit.se = TRUE)
        # obj_diffgroup <- analyzeComparative(obj = obj, dnaDesign = ~ 1 + group + replicate, rnaDesign = ~ 1 + replicate + group, reducedDesign = ~1 + replicate, fit.se = TRUE)
      }
      toc()
      LRT <- as.data.table(testLrt(obj_diffgroup))[, .(logFC, statistic, df.test, pval.LRT = pval)]
      logFC.se <- obj_diffgroup@modelFits$r.se[1 + sum(apply(obj_diffgroup@designs@rnaFull, 2, sum) != 0)]
      if (is.null(logFC.se)) {
        logFC.se <- NA
      }

      toc()
      alphas_and_pvalues[[test]] <- data.table(oligo = test, GC = oligo_data[, mean(GC)], alpha_quantif, LRT, logFC.se, oligo_stats[, -"oligo"])
      # fwrite(alphas_and_pvalues, file = sprintf("%s/alphas_and_pvalues_v%s_oligo%s_max%sBC_Z%s.txt", OUTS_DIR, VERSION_OUT, my_oligo_num, sub_sample, Outlier_Z_threshold), sep = "\t")
      cat(sprintf("\n oligo number %s processed (%s)\n", counter, test))
    }
  }
}
alphas_and_pvalues <- rbindlist(alphas_and_pvalues)
output_name <- sprintf("%s/alphas_and_pvalues_v%s_%s_chunk%s_max%sBC_Z%s.txt", OUTS_DIR, VERSION_OUT, analysis_unit, chunk_nb, sub_sample, Outlier_Z_threshold)

cat(nrow(alphas_and_pvalues), "test performed")

cat("\nprinting to", output_name)
fwrite(alphas_and_pvalues, file = output_name, sep = "\t")

cat("\nAll done")
q("no")
