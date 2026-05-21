MPRA_DIR <- "/pasteur/helix/projects/evo_immuno_pop/MPRA"
ANALYSIS_DIR <- "MPRA_count_exp6_analysis3"


source(sprintf("%s/scripts/%s/00_starter.R", MPRA_DIR, ANALYSIS_DIR))
source(sprintf("%s/scripts/%s/00_load_associations.R", MPRA_DIR, ANALYSIS_DIR))

VERSION_OUT <- 1

cmd <- commandArgs()
print(cmd)
for (i in seq_along(cmd)) {
  if (cmd[i] == "--sample" || cmd[i] == "-s") {
    j <- as.numeric(cmd[i + 1])
  } # ID of the set to test
  if (cmd[i] == "--version" || cmd[i] == "-v") {
    VERSION_OUT <- as.numeric(cmd[i + 1])
  } # ID of the set to test
}


# deactivate next line after HN00196048 data has been reformatted
# Indexes <- Indexes[SeqRUN != "HN00196048", ]
SAMPLE_ID <- Indexes[, unique(Sample_ID)[j]]
###### random subset of reads (for testeing purposes)
# FASTQ_1_file=sprintf('%s/outs/test_UMItools/sub1.fastq.gz',MPRA_DIR)
# FASTQ_2_file=sprintf('%s/outs/test_UMItools/sub2.fastq.gz',MPRA_DIR)

# tim=Sys.time()
# FASTQ_1=readDNAStringSet(FASTQ_1_file, format="fastq")
# print(Sys.time()-tim)
# FASTQ_2=readDNAStringSet(FASTQ_2_file, format="fastq")
# print(Sys.time()-tim)


        #  barcode1   barcode2              UMI      index     order known_BC1 known_BC2 sameBC nT_BC nA_BC nG_BC nC_BC   shannon nT_UMI nA_UMI nG_UMI nC_UMI shannon_UMI               library
        #             <char>     <char>           <char>     <char>     <int>    <lgcl>    <lgcl> <lgcl> <int> <int> <int> <int>     <num>  <int>  <int>  <int>  <int>       <num>                <char>
        # 1: NGGGTTAGGGGGCCG ATGCCTCAAA AAATCATGCCTCAAAA TTTGAGGCAT         1     FALSE     FALSE  FALSE     2     1     9     2 13.005380      3      8      1      4   15.493841 A549-ACE2_TNFa_R1_RNA
        # 2: GGTCCTCTTGGTAAA TTATCGATTA ACTCTTTATCGATTAA TAATCGATAA         2     FALSE     FALSE  FALSE     5     3     4     3 17.874516      7      5      1      3   15.994501 A549-ACE2_TNFa_R1_RNA
        # 3: GTCTTTAGTGTGATT ACTCTAACCC ACCTCACTCTAACCCA GGGTTAGAGT         3     FALSE     FALSE  FALSE     8     2     4     1 14.682911      3      5      0      8   13.665943 A549-ACE2_TNFa_R1_RNA
        # 4: GTTGTTTTTGTTTTT CACAAAAAAA AACCCCACAAAAAAAA TTTTTTTGTG         4     FALSE     FALSE  FALSE    12     0     3     0  6.650555      0     11      0      5    8.068649 A549-ACE2_TNFa_R1_RNA
        # 5: GAACCACCGTTCGCG CTTAAGAACA AAAACCTTAAGAACAA TGTTCTTAAG


myIndexes <- Indexes[Sample_ID == SAMPLE_ID, ]
DT_MPRA <- list()
for (i in myIndexes[, seq_len(.N)]) {
  SAMPLE_ID <- myIndexes[i, Sample_ID]
  SEQRUN <- myIndexes[i, SeqRUN]
  EXP <- myIndexes[i, experiment]
  BC_LIB <- myIndexes[i, barcode_lib]
  LANE <- myIndexes[i, Lane]
  SID <- myIndexes[i, SID]
  REP <- myIndexes[i, replicate]
  if(EXP == 'exp5' & SEQRUN=='B7075-9'){
    FASTQ_R1_file <- sprintf("%s/data/fastq/%s/%s_%s_L%s_R2_001.fastq.gz", MPRA_DIR, SEQRUN, SAMPLE_ID, SID, str_pad(LANE, 3, "left", "0"))
    FASTQ_UMI_file <- sprintf("%s/data/fastq/%s/%s_%s_L%s_R1_001.fastq.gz", MPRA_DIR, SEQRUN, SAMPLE_ID, SID, str_pad(LANE, 3, "left", "0"))
    FASTQ_R2_file <- sprintf("%s/data/fastq/%s/%s_%s_L%s_R3_001.fastq.gz", MPRA_DIR, SEQRUN, SAMPLE_ID, SID, str_pad(LANE, 3, "left", "0"))
    FASTQ_INDEX_file <- sprintf("%s/data/fastq/%s/%s_%s_L%s_I1_001.fastq.gz", MPRA_DIR, SEQRUN, SAMPLE_ID, SID, str_pad(LANE, 3, "left", "0"))
  }else{
    FASTQ_R1_file <- sprintf("%s/data/fastq/%s/%s_%s_L%s_R1_001.fastq.gz", MPRA_DIR, SEQRUN, SAMPLE_ID, SID, str_pad(LANE, 3, "left", "0"))
    FASTQ_UMI_file <- sprintf("%s/data/fastq/%s/%s_%s_L%s_R2_001.fastq.gz", MPRA_DIR, SEQRUN, SAMPLE_ID, SID, str_pad(LANE, 3, "left", "0"))
    FASTQ_R2_file <- sprintf("%s/data/fastq/%s/%s_%s_L%s_R3_001.fastq.gz", MPRA_DIR, SEQRUN, SAMPLE_ID, SID, str_pad(LANE, 3, "left", "0"))
    FASTQ_INDEX_file <- sprintf("%s/data/fastq/%s/%s_%s_L%s_I1_001.fastq.gz", MPRA_DIR, SEQRUN, SAMPLE_ID, SID, str_pad(LANE, 3, "left", "0"))
  }
  FASTQ_R1 <- readDNAStringSet(FASTQ_R1_file, format = "fastq")
  FASTQ_UMI <- readDNAStringSet(FASTQ_UMI_file, format = "fastq")
  FASTQ_R2 <- readDNAStringSet(FASTQ_R2_file, format = "fastq")
  FASTQ_INDEX <- readDNAStringSet(FASTQ_R2_file, format = "fastq")

  INDEX_length <- 10
  UMI_length <- 16
  BC1_length <- 15
  BC2_length <- ifelse(EXP=='exp5' & SEQRUN=='B7075-9',10,15)

  DT_MPRA[[i]] <- data.table(
    barcode1 = as.character(subseq(FASTQ_R1, 1, BC1_length)),
    barcode2 = as.character(reverseComplement(subseq(FASTQ_R2, 1, BC2_length))),
    UMI = as.character(subseq(FASTQ_UMI,1,UMI_length)),
    index = as.character(subseq(FASTQ_INDEX,1,INDEX_length))
  )
  # this is step isdone for each index to account for the possibility that BC2_length may vary for some indexes
  DT_MPRA[[i]][, known_BC1 := barcode1 %chin% Associations_Filtered$BC]
  DT_MPRA[[i]][, known_BC2 := substr(barcode2, 1, BC2_length) %chin% substr(Associations_Filtered$BC,1,BC2_length)]
  # if barcode2 is shorter than barcode 1, we compare it to the last bases of barcode 1
  DT_MPRA[[i]][, sameBC := substr(barcode2, 1, BC2_length) == substr(barcode1, 1+BC1_length-BC2_length, BC1_length)]

}
DT_MPRA <- rbindlist(DT_MPRA)
DT_MPRA[, order := seq_len(.N)]

# R1primer_seq='GCAAAGTGAACACATCGCTAAGCGAAAGCTAAG'
# R1primer_revC='CTTAGCTTTCGCTTAGCGATGTGTTCACTTTGC'
#
# R2primer_seq='GCTCCTCGCCCTTGCTCACCATGGTGGCGACCGGT'
# R2primer_revC='ACCGGTCGCCACCATGGTGAGCAAGGGCGAGGAGC'
#
# P7_seq='CAAGCAGAAGACGGCATACGAGAT'
# P7_revC='ATCTCGTATGCCGTCTTCTGCTTG'
#
# P5_seq='AATGATACGGCGACCACCGAGATCTACAC'
# P5_revC='GTGTAGATCTCGGTGGTCGCCGTATCATT'


tim <- Sys.time()
DT_MPRA[, nT_BC := str_count(barcode1, "T")]
DT_MPRA[, nA_BC := str_count(barcode1, "A")]
DT_MPRA[, nG_BC := str_count(barcode1, "G")]
DT_MPRA[, nC_BC := str_count(barcode1, "C")]
DT_MPRA[, shannon := -(nT_BC * log((1 + nT_BC) / 16) + nC_BC * log((1 + nC_BC) / 16) + nG_BC * log((1 + nG_BC) / 16) + nA_BC * log((1 + nA_BC) / 16))]
print(Sys.time() - tim)

tim <- Sys.time()
DT_MPRA[, nT_UMI := str_count(UMI, "T")]
DT_MPRA[, nA_UMI := str_count(UMI, "A")]
DT_MPRA[, nG_UMI := str_count(UMI, "G")]
DT_MPRA[, nC_UMI := str_count(UMI, "C")]
DT_MPRA[, shannon_UMI := -(nT_UMI * log((1 + nT_UMI) / 16) + nC_UMI * log((1 + nC_UMI) / 16) + nG_UMI * log((1 + nG_UMI) / 16) + nA_UMI * log((1 + nA_UMI) / 16))]
print(Sys.time() - tim)

dir.create(sprintf("%s/data/%s/01a_MPRA_results/", MPRA_DIR, ANALYSIS_DIR))
fwrite(DT_MPRA, file = sprintf("%s/data/%s/01a_MPRA_results/MPRA_results__%s__v%s.txt.gz", MPRA_DIR, ANALYSIS_DIR, SAMPLE_ID, VERSION_OUT), sep = "\t")
# DT_MPRA=fread(sprintf('%s/data/%s/01a_MPRA_results__%s__%s__v%s.txt.gz',MPRA_DIR,ANALYSIS_DIR,SID,SAMPLE_ID,VERSION_OUT))
cat("\nAll done\n")