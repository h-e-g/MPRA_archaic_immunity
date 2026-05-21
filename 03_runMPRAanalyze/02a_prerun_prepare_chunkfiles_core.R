MPRA_DIR <- "/pasteur/helix/projects/evo_immuno_pop/MPRA"
ANALYSIS_DIR <- "MPRA_count_exp6_analysis3"

source(sprintf("%s/scripts/misc_plots.R", MPRA_DIR))
source(sprintf("%s/scripts/%s/00_starter.R", MPRA_DIR, ANALYSIS_DIR))
source(sprintf("%s/scripts/%s/00_load_associations.R", MPRA_DIR, ANALYSIS_DIR))


library(MPRAnalyze)

FIGURE_DIR <- sprintf("%s/figures/%s", MPRA_DIR, ANALYSIS_DIR)

# parameters defaults
VERSION_IN <- 1
VERSION_OUT <- 1
instruction_file <- sprintf("%s/scripts/%s/02a_MPRA_analyse_instructions_full_sorted.txt", MPRA_DIR, ANALYSIS_DIR)
sub_sample <- 2000
min_nBC <- 10 # not used so far
outlier_threshold <- 3
target_time_s=3600 # target of 1 hour per job (set limit between 2 and 4 h)
cst_time=6 # values based on empirical test + some margin
linear_complexity=0.004 # values based on empirical test + some margin
squared_complexity=0 # values based on empirical test + some margin

cmd <- commandArgs()
print(cmd)

for (i in seq_along(cmd)) {
  if (cmd[i] == "--version_in" | cmd[i] == "-vi") {
    VERSION_IN <- as.numeric(cmd[i + 1])
  } # version of input files
  if (cmd[i] == "--version_out" | cmd[i] == "-vo") {
    VERSION_OUT <- as.numeric(cmd[i + 1])
  } # version of output files
  if (cmd[i] == "--sub_sample" | cmd[i] == "-u") {
    if (grepl("%", cmd[i + 1])) {
      SUB_SAMPLE_Pct <- TRUE
      SUB_SAMPLE <- as.numeric(gsub("%", "", cmd[i + 1])) / 100
    } else {
      SUB_SAMPLE <- as.numeric(cmd[i + 1])
    }
  } # subsampling parameters
  if (cmd[i] == "--instruction_file" | cmd[i] == "-i") {
    instruction_file <- cmd[i + 1]
  } # instruction file
  if (cmd[i] == "--instruction_line" | cmd[i] == "-l") {
    instruction_line <- as.numeric(cmd[i + 1])
  } # instruction file
  if (cmd[i] == "--min_nBC" | cmd[i] == "-n") {
    min_nBC <- as.numeric(cmd[i + 1])
  } # minimum number of BCs in a oligo
  if (cmd[i] == "--outlier_threshold" | cmd[i] == "-z") {
    outlier_threshold <- as.numeric(cmd[i + 1])
  } # threshold for outlier detection
  if (cmd[i] == "--target_time_s" | cmd[i] == "-t") {
    target_time_s <- as.numeric(cmd[i + 1])
  } # time target for each analysis
  if (cmd[i] == "--cst_time" | cmd[i] == "-c") {
    cst_time <- as.numeric(cmd[i + 1])
  } # time constant for time estimation
  if (cmd[i] == "--linear_complexity" | cmd[i] == "-lc") {
    linear_complexity <- as.numeric(cmd[i + 1])
  } # linear complexity of the analysis
  if (cmd[i] == "--squared_complexity" | cmd[i] == "-sc") {
    squared_complexity <- as.numeric(cmd[i + 1])
  } # squared complexity of the analysis
  if (cmd[i] == "--output_dir" || cmd[i] == "-o") {
     OUT_DIR <- cmd[i + 1]
  }
}


time_estimate <- function(nBC){ cst_time + nBC * linear_complexity}

instruction_DT <- fread(instruction_file, header = TRUE, fill = TRUE)
if(all(is.na(instruction_DT[, group2_samples]))){
  instruction_DT[, group2_samples:=NULL]
  instruction_DT[, group2_samples:='']
}

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
load_time=c()
UMI_perBC_annot <- list()
for (mySETUP_ID in unique(LIB_summary[, SETUP_ID])) {
#  for (MATERIAL in c("RNA", "DNA")) {
    # get sample ID
 #   SAMPLE_ID <- Indexes[SETUP_ID == mySETUP_ID & material == MATERIAL, unique(Sample_ID)]
    # SID <- Indexes[SETUP_ID == mySETUP_ID & material == MATERIAL, unique(SID)]
    cat("reading UMI counts for", mySETUP_ID, "\n")
    # read UMI counts
    tim=Sys.time()
    UMI_perBC_annot[[mySETUP_ID]] <- fread(sprintf("%s/data/%s/01c_MPRA_results_BC_activity/01c_MPRA_results_BC_activity__%s__v%s.txt.gz", MPRA_DIR, ANALYSIS_DIR, mySETUP_ID, VERSION_IN))
    load_time[mySETUP_ID] <- Sys.time() - tim
  }

# merge data from RNA and DNA samples & add library info
UMI_perBC_annot <- rbindlist(UMI_perBC_annot, use.names=TRUE)
toc()
# annotate BCs
UMI_per_BC_annot_oligo <- merge(UMI_perBC_annot, unique(oligo_source)[oligo != "", .(oligo, source, type, shift, SNP, strand, crsID, allele, GC)], by = "oligo", allow.cartesian = TRUE)
# filter out outlier BCs 
UMI_per_BC_annot_oligo <- UMI_per_BC_annot_oligo[abs(outlier_Z)<outlier_threshold,]
# subsampling Barcodes per oligo (ie. CRS + allele)
set.seed(sub_sample)
UMI_per_BC_annot_oligo <- UMI_per_BC_annot_oligo[, .SD[sample(seq_len(.N), min(.N, sub_sample)), ], by = .(oligo, SETUP_ID, allele)]
UMI_per_BC_long <- melt(UMI_per_BC_annot_oligo, measure.vars=c('RNA','DNA'), value.name=c('nUMI_perBC'), variable.name='material',id.vars=c('barcode1','barcode_lib','oligo','celltype', 'replicate', 'SETUP_ID', 'condition'))

#for (instruction_line in seq_len(instruction_DT[, .N])) {
  cat("\n", instruction_line, ":")
  ANALYSIS_NAME <- instruction_DT[instruction_line, analysis_name]
  ANALYSIS_TYPE <- instruction_DT[instruction_line, analysis_type]
  ANALYSIS_SUBTYPE <- instruction_DT[instruction_line, analysis_subtype]
  group_1_samples <- instruction_DT[instruction_line, group1_samples]
  if (group_1_samples != "") {
    group_1_samples <- as.character(str_split(group_1_samples, ",", simplify = TRUE))
    if (any(!group_1_samples %in% Indexes[, unique(SETUP_ID)])) {
      stop(sprintf("samples %s not found in Indexes", paste(group_1_samples[!group_1_samples %in% Indexes[, unique(SETUP_ID)]], collapse = ", ")))
    }
  } else {
    stop("no samples specified")
  }

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

  if (grepl("EMVAR", toupper(ANALYSIS_TYPE))) {
    Ntest <- UMI_per_BC_annot_oligo[SETUP_ID %chin% c(group_1_samples, group_2_samples), ][crsID %chin% oligo_source[,unique(crsID)], ][!duplicated(crsID), .N]
  } else {
    Ntest <- UMI_per_BC_annot_oligo[SETUP_ID %chin% c(group_1_samples, group_2_samples), ][!duplicated(oligo), .N]
  }
  instruction_DT[instruction_line, N_test := Ntest]  
  # compute data loading time
  testdata_load_time <- sum(load_time[ c(group_1_samples, group_2_samples)])
  # count BC per test (SNP or oligo)
  if (grepl("EMVAR", toupper(ANALYSIS_TYPE))) {
    nBC_per_test <- UMI_per_BC_annot_oligo[ SETUP_ID %chin% c(group_1_samples, group_2_samples), ][crsID %chin% oligo_source[,unique(crsID)], ][,.N,by=crsID]
  } else {
    nBC_per_test <- UMI_per_BC_annot_oligo[ SETUP_ID %chin% c(group_1_samples, group_2_samples), ][,.N,by=oligo]
  }
  # estimate running time per test & order by decreasing time
  nBC_per_test[,time_per_test_s:=time_estimate(N)]
  nBC_per_test=nBC_per_test[order(-time_per_test_s)]
  
  # estimate total running time & deduce number of chunks
  cat(sprintf("Estimated time to run the full analysis: %s min\n", round(nBC_per_test[,sum(time_per_test_s)]/60, 2)))
  N_CHUNK=ceiling(nBC_per_test[,sum(time_per_test_s)/(target_time_s-testdata_load_time)])
  instruction_DT[instruction_line, N_chunk := N_CHUNK]
  # define optimal chunks
  cat(sprintf("splitting in %s chunks\n", N_CHUNK))
  CHUNK_DT=data.table(chunk_nb=seq_len(N_CHUNK), chunk_time_s=testdata_load_time)
  # initialize chunks
  CHUNK_tests=list()
  for (chunk_i in seq_len(N_CHUNK)){
    if(grepl("EMVAR", toupper(ANALYSIS_TYPE))){
      test_ID=nBC_per_test[chunk_i,crsID]
    }else{
      test_ID=nBC_per_test[chunk_i,oligo]
    }
    CHUNK_tests[[chunk_i]]=test_ID
  }
  # progressively add each SNP/oligo to the chunk with lowest running time
  for (test_i in (N_CHUNK+1):nrow(nBC_per_test)) {
    if(test_i%%100==0){cat(sprintf("%s/%s\n",test_i,nrow(nBC_per_test)))}
    # find chunk with lowest running time
    chunk_number <- CHUNK_DT[,which.min(chunk_time_s)]
    if(grepl("EMVAR", toupper(ANALYSIS_TYPE))){
      test_ID=nBC_per_test[test_i,crsID]
    }else{
      test_ID=nBC_per_test[test_i,oligo]
    }
    # update chunk & associated running time
    CHUNK_tests[[chunk_number]]=c(CHUNK_tests[[chunk_number]],test_ID)
    CHUNK_DT[chunk_number,chunk_time_s:=chunk_time_s+nBC_per_test[test_i,time_per_test_s]]
  }
  # print chunks 
  dir.create(sprintf("%s/chunk_files/%s/%s/%s/%s/%s",MPRA_DIR,ANALYSIS_DIR,OUT_DIR,ANALYSIS_TYPE,ANALYSIS_SUBTYPE,ANALYSIS_NAME), recursive = TRUE, showWarnings = FALSE)
  for (chunk_number in seq_len(N_CHUNK)){
    chunk_file=sprintf("%s/chunk_files/%s/%s/%s/%s/%s/chunk_%s_%s.txt",MPRA_DIR, ANALYSIS_DIR,OUT_DIR,ANALYSIS_TYPE,ANALYSIS_SUBTYPE,ANALYSIS_NAME,instruction_line, chunk_number)
    fwrite(data.table(instruction_line,ANALYSIS_TYPE,ANALYSIS_SUBTYPE,ANALYSIS_NAME,test=CHUNK_tests[[chunk_number]]),file=chunk_file,sep='\t')
  }
#}

# fwrite(instruction_DT, file = instruction_file, sep = "\t")
# fwrite(instruction_DT[grepl('TNFa|NS_R3|NS_all',analysis_name) & !grepl('A549_NS_R1',group1_samples)], file = sprintf("%s/scripts/%s/02a_MPRA_analyse_instructions_short.txt", MPRA_DIR, ANALYSIS_DIR), sep = "\t")
cat('All done!')
q('no')


##TODO: genrate all chunk files
# instruction line and chunk number are argumens of core script
# only the BC/UMI data relevant to the instruction line are read
# add reading the chunk file, & looping through chunk file into core script
# generate laucher script that runs through instruction lines and runs all chunks
# permute etc are arguments of launcher script