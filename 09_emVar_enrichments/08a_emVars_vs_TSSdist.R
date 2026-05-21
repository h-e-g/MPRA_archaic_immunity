MPRA_DIR <- "/pasteur/helix/projects/evo_immuno_pop/MPRA"
ANALYSIS_DIR <- "MPRA_count_exp6_analysis3"
EIP = "/pasteur/helix/projects/evo_immuno_pop"
RESOURCE_DIR = sprintf("%s/single_cell/resources", EIP)


source(sprintf("%s/scripts/%s/00_starter.R", MPRA_DIR, ANALYSIS_DIR))
source(sprintf("%s/scripts/misc_plots.R", MPRA_DIR))
#source(sprintf("%s/scripts/%s/00_load_associations.R", MPRA_DIR, ANALYSIS_DIR))

RUN_ID <- "RUN3_Z2_nBC10"
CRITERIA_ACTIVE_OUT <- "FDR5_scrambled5pct_FC0.2_GCnorm"
CRITERIA_EMVARS_OUT <- "EmVar_FDR5_FC.2"
EMVAR_LIST <- "significant"
cmd <- commandArgs(trailingOnly = TRUE)
print(cmd)

for (i in seq_along(cmd)) {
  if (cmd[i] == "--run_dir" || cmd[i] == "-r") {
    RUN_ID <- cmd[i + 1]
  }
  if (cmd[i] == "--criteria_active" || cmd[i] == "-a") {
    CRITERIA_ACTIVE_OUT <- cmd[i + 1]
  }
  if (cmd[i] == "--criteria_emvar" || cmd[i] == "-e") {
    CRITERIA_EMVARS_OUT <- cmd[i + 1]
  }
  if (cmd[i] == "--emVar_list" || cmd[i] == "-l") {
    EMVAR_LIST <- cmd[i + 1] # either 'suggestive' or 'significant'
  }
}

if(!EMVAR_LIST %in% c('suggestive','significant')) {
  cat("WARNING: --emVar_list must be either 'suggestive' or 'significant', setting to 'significant' and continuing\n")
  EMVAR_LIST <- 'significant'
}

CRITERIA_ACTIVE <- gsub("noActivityFilter_", "", CRITERIA_ACTIVE_OUT)
CRITERIA_EMVARS <- gsub("noEmVarFilter_", "", CRITERIA_EMVARS_OUT)


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
##### emVars parameters
source(sprintf("%s/scripts/%s/05_00_parameter_emVars.R", MPRA_DIR, ANALYSIS_DIR))

IN_DIR <- sprintf("%s/data/%s/03a_aggMPRA_analyse/%s/", MPRA_DIR, ANALYSIS_DIR, RUN_ID)
ACTIVE_DIR <- sprintf("%s/data/%s/03b_computeFDR/%s/Activity/%s/", MPRA_DIR, ANALYSIS_DIR, RUN_ID, CRITERIA_ACTIVE)
EMVAR_DIR <- sprintf("%s/data/%s/03b_computeFDR/%s/emVars/%s/%s", MPRA_DIR, ANALYSIS_DIR, RUN_ID, CRITERIA_EMVARS, CRITERIA_ACTIVE_OUT)

library(ggplot2)

FIGURE_DIR <- sprintf("%s/figures/%s/%s/04_emVars/%s/%s/", MPRA_DIR, ANALYSIS_DIR, RUN_ID, CRITERIA_EMVARS, CRITERIA_ACTIVE_OUT)
dir.create(sprintf("%s/Pct_emVar_by_Dist", FIGURE_DIR), recursive = TRUE)
dir.create(sprintf("%s/Pct_emVar_by_ChromHMM", FIGURE_DIR), recursive = TRUE)

##### load oligo activity
oligo_activity_file <- sprintf("%s/oligo_activity__all__%s.tsv.gz", ACTIVE_DIR, CRITERIA_ACTIVE)
oligo_activity_class <- fread(oligo_activity_file)

oligo_activity_file <- sprintf("%s/all_oligos_annotated__%s.tsv.gz", ACTIVE_DIR, CRITERIA_ACTIVE)
oligo_activity_obs <- fread(oligo_activity_file)

emVARs_obs <- fread(file = sprintf("%s/all_emVars_annotated_celltype__%s.tsv", EMVAR_DIR, CRITERIA_ACTIVE_OUT))
emVARs_obs[, posID := crsID]
# emVARs_obs <- emVARs_obs[is_emVar_CRITERIA==TRUE,]


crs_targets <- fread(sprintf("%s/data/%s/00_oligo_targets_v4.txt", MPRA_DIR, ANALYSIS_DIR))

SNP_annot_emVars <- fread(file = sprintf("%s/SNP_annot_emVars__%s.tsv", EMVAR_DIR, CRITERIA_ACTIVE_OUT))
SNP_annot_emVars <- merge(SNP_annot_emVars, unique(crs_targets[method == "NearestTSS" & celltype == "any", .(posID = crsID, DistTSS_kb = -Score1, NearestTSS = TargetGene)]), by = c("posID"), allow.cartesian = TRUE, all.x = TRUE)

SNP_annot_emVars[nCelltype_NS > 0, mean(nCelltype_NS == 1)]
SNP_annot_emVars[nCelltype_NS == 1, mean(nCelltype_P05_NS > 1)]
SNP_annot_emVars[nCelltype>0,wilcox.test(DistTSS_kb[nCelltype_P05==3],DistTSS_kb[nCelltype_P05==1])]
# W = 23062, p-value = 0.9123
# alternative hypothesis: true location shift is not equal to 0

SNP_annot_emVars <- merge(SNP_annot_emVars, SNP_annot[, .(ID, crsID, gwas, LD_block_AFR, introgression_locus)], by = "ID")
SNP_annot_gwas <- SNP_annot_emVars[, .(nSNP = .N, 
                                      has_gwas = any(gwas != ""), 
                                      nTraits = length(unique(na.omit(gsub("(.*) - ", "", unlist(strsplit(gwas, "/")))))),
                                      traits = paste(unlist(unique(na.omit(gsub("(.*) - ", "", unlist(strsplit(gwas, "/")))))), collapse = "/"),
                                      nEmVar = sum(nCelltype_NS != 0), 
                                      nShared_EmVars = sum((nCelltype_NS > 0) & (nCelltype_P05_NS > 1))
                                      ), by = introgression_locus]
SNP_annot_gwas[, .(introgression_locus, trait = unique(gsub("(.*) - ", "", unlist(strsplit(traits, "/")))))]

# 
set.seed(42)
mean(replicate(1000, {
   SNP_annot_emVars[sample(seq_len(.N), replace = F), ][!duplicated(introgression_locus)][nCelltype_NS > 0 & nCelltype_P05_NS == 1, mean(gwas != "")]
 }), na.rm = T)
 # [1] 0.333908

set.seed(42)
 mean(replicate(1000, {
   SNP_annot_emVars[sample(seq_len(.N), replace = F), ][!duplicated(introgression_locus)][nCelltype_NS > 0 & nCelltype_P05_NS > 1, mean(gwas != "")]
 }), na.rm = T)
 # [1] 0.1988697

set.seed(42)
mean(replicate(1000, {
   SNP_annot_emVars[sample(seq_len(.N), replace = F), ][!duplicated(introgression_locus)][nCelltype > 0 & nCelltype_P05== 1, mean(gwas != "")]
 }), na.rm = T)
# 0.3657831
set.seed(42)
mean(replicate(1000, {
   SNP_annot_emVars[sample(seq_len(.N), replace = F), ][!duplicated(introgression_locus)][nCelltype > 0 & nCelltype_P05== 3, mean(gwas != "")]
 }), na.rm = T)
# 0.3297125


annot_crs <- unique(oligo_source[, .(crsID, E114, E118, E123)])
annot_crs[, E114_short := case_when(grepl("1_|2_", E114) ~ 3, grepl("6_|7_", E114) ~ 2, TRUE ~ 0)]
annot_crs[, E118_short := case_when(grepl("1_|2_", E118) ~ 3, grepl("6_|7_", E118) ~ 2, TRUE ~ 0)]
annot_crs[, E123_short := case_when(grepl("1_|2_", E123) ~ 3, grepl("6_|7_", E123) ~ 2, TRUE ~ 0)]
SNP_annot_emVars <- merge(SNP_annot_emVars, annot_crs, by = "crsID")

SNP_annot_emVars[,chromHmm_annot:=pmax(E114_short,E118_short,E123_short)]
SNP_annot[,POS_B38:=as.numeric(gsub('.*_([0-9]+)_.*_.*','\\1',variantId_hg38))]


CHROMHMM_FULLSTACK = fread(sprintf("%s/references/ATAC/ChromHMM/hg38_genome_100_segments.bed.gz", RESOURCE_DIR), col.names = c("CHR", "START", "END", "ANNOTATION"))
CHROMHMM_FULLSTACK = makeGRangesFromDataFrame(CHROMHMM_FULLSTACK, keep.extra.columns = TRUE)
seqlevelsStyle(CHROMHMM_FULLSTACK)='NCBI'

SNP_annot_GR=makeGRangesFromDataFrame(SNP_annot[,.(posID,POS_B38,CHROM)], keep.extra.columns = TRUE,start.field='POS_B38', end.field='POS_B38')
oo <- findOverlaps(SNP_annot_GR,CHROMHMM_FULLSTACK)
SNP_annot_GR$ANNOTATION=NA
SNP_annot_GR$ANNOTATION[queryHits(oo)]=CHROMHMM_FULLSTACK[subjectHits(oo),]$ANNOTATION
CHROMHMM_FULLSTACK_METADATA = fread(sprintf("%s/references/ATAC/ChromHMM/state_annotations_processed.csv", RESOURCE_DIR))
CHROMHMM_FULLSTACK_METADATA[,ANNOTATION:=paste(state_order_by_group ,mneumonics,sep='_')]
SNP_annot_GR <- merge(as.data.table(SNP_annot_GR),CHROMHMM_FULLSTACK_METADATA[,.(ANNOTATION, Group)],by='ANNOTATION')
SNP_annot_emVars <- merge(SNP_annot_emVars, SNP_annot_GR, by = "posID")

#SNP_annot_emVars[nCelltype>0 & nCelltype_P05!=2,chisq.test(table(chromHmm_annot,nCelltype_P05))]
SNP_annot_emVars[,chromHmm_diff:=pmax(E114_short,E118_short,E123_short)-pmin(E114_short,E118_short,E123_short)]


# emVARs_annot <- fread(file = sprintf("%s/all_emVars_annotated.tsv.gz", IN_DIR), sep = "\t")
# emVARs_annot <- emVARs_annot[power == 0 & perm == "perm_0_0_0" & grepl("celltype", ANALYSIS_SUBTYPE), ]


# # from 00f
# introgressed_SNPs_annot <- list()
# for (CHR in 1:22) {
#   cat("CHR", CHR, "\n")
#   introgressed_SNPs_annot[[CHR]] <- try(fread(sprintf("%s/data/vcfs/annotate_introgressed/full_annotation/haplotypes_aSNP_with_introgressed_and_frequency_%s.tsv.gz", MPRA_DIR, CHR)))
# }
# # introgressed_SNPs_annot[[6]]=NULL
# introgressed_SNPs_annot <- rbindlist(introgressed_SNPs_annot)

# INTRO.freqs <- introgressed_SNPs_annot[, .(ID, POP, freq_introgressed_allele=target_SNP_freq.fixed, freq_introgressed_haplotype=target_SNP_on_archaic_hap_freq.fixed,  GUESSED_INTROGRESSED.fixed, GUESSED_INTROGRESSED.nonfixed)]

emVARs_annot <- merge(emVARs_obs, unique(crs_targets[method == "NearestTSS" & celltype == "any", .(crsID, DistTSS_kb = -Score1, NearestTSS = TargetGene)]), by = c("crsID"), allow.cartesian = TRUE, all.x = TRUE)
emVARs_annot <- merge(emVARs_annot, unique(crs_targets[method == "NearestTSS" & celltype != "any", .(crsID, celltype, DistActiveTSS_kb = -Score1, NearestActiveTSS = TargetGene)]), by = c("crsID", "celltype"), allow.cartesian = TRUE, all.x = TRUE)
emVARs_annot <- merge(emVARs_annot, unique(oligo_activity_obs[, .(crsID, celline, ROADMAP_current_celline)]), by = c("crsID", "celline"), allow.cartesian = TRUE, all.x = TRUE)



Dist_breaks <- c(-1, .5, 1, 2, 5, 10, 20, 30, 40, 50, 100, 200, 500, Inf)
Dist_breaks <- c(-1, 1, 2, 5, 10, 20, 30, 40, 50, 100, 200, 500, Inf)
# Dist_breaks= c(-1,2,5,10,20,50,100,200,500,Inf)
# Dist_breaks= c(-1,10,50,100,500,Inf)
emVARs_annot[, distTSS_bin := cut(DistTSS_kb, Dist_breaks)]
emVARs_annot[, distActiveTSS_bin := cut(DistActiveTSS_kb, Dist_breaks)]

#####################################################################################################################
####################################### Pct emVar by Distance to TSS ################################################
#####################################################################################################################

minSNP <- 5

Pct_emVar_byDistTSS <- get_Pct_emVars(emVARs_annot, split_by = c("distTSS_bin"), total_by = NULL, emVarlist = EMVAR_LIST)
Pct_emVar_byDistTSS <- make_it_long(Pct_emVar_byDistTSS, split_by = c("distTSS_bin"))
Pct_emVar_byDistTSS[N_test > minSNP, FDR := p.adjust(Pvalue, "fdr"), by = measure]
Pct_emVar_byDistTSS[N_test > minSNP][order(Pvalue), head(.SD, 5), by = measure]
fwrite(Pct_emVar_byDistTSS, sprintf("%s/Pct_emVar_by_Dist/01a_Pct_emVar_%s_byDistTSS.tsv", FIGURE_DIR, EMVAR_LIST), sep = "\t", quote = FALSE)

for (i_MEASURE in seq_len(Measure_table[, .N])) {
  MEASURE <- Measure_table[i_MEASURE, measure]
  MEASURE_LABEL <- Measure_table[i_MEASURE, measure_label]
  MEASURE_NULL <- Measure_table[i_MEASURE, null]

  FigData <- Pct_emVar_byDistTSS[measure == MEASURE & N_test > minSNP, ]
  p <- ggplot(FigData)
  p <- p + geom_hline(yintercept = MEASURE_NULL, col = "lightgrey", linetype = "dashed")
  p <- p + geom_pointrange(aes(x = distTSS_bin, y = Pct * 100, ymin = Pct_lo * 100, ymax = Pct_hi * 100, alpha = ifelse(FDR < 0.05,'FDR<0.05','ns')), size = .3)
  # p <- p + facet_grid(rows=vars(measure)) + ylab("Percentage of alleles") + xlab("Population")
  if (is.na(MEASURE_NULL)) {
    MEASURE_NULL <- as.numeric(Pct_emVar_byDistTSS[measure == MEASURE, .(Pct = median(Pct))])
  }
  p <- p + geom_hline(yintercept = MEASURE_NULL * 100, linetype = 2)
  p <- p + theme_plot(rotate.x = 90, fontsize = 9) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + guides(col = "none") + xlab("")
  p <- p + scale_alpha_manual(values = setNames(alpha_TRUEFALSE,c('ns','FDR<0.05')))
  # p <- p + scale_color_manual(values=color_populations_MPRA)
  p <- p + ylab(MEASURE_LABEL) + xlab("Distance to Nearest TSS (kb)")
  pdf(sprintf("%s/Pct_emVar_by_Dist/01a_Pct_%s_%s_by_DistTSS.pdf", FIGURE_DIR, MEASURE, EMVAR_LIST), height = 2.3, width = 1.5)
  print(p)
  dev.off()
}


Pct_emVar_byDistActiveTSS <- get_Pct_emVars(emVARs_annot, split_by = c("distActiveTSS_bin", "celline"), total_by = "celline", emVarlist = EMVAR_LIST)
Pct_emVar_byDistActiveTSS <- make_it_long(Pct_emVar_byDistActiveTSS, split_by = c("distActiveTSS_bin", "celline"))
Pct_emVar_byDistActiveTSS[N_test > minSNP, FDR := p.adjust(Pvalue, "fdr"), by = measure]
Pct_emVar_byDistActiveTSS[N_test > minSNP][order(Pvalue), head(.SD, 5), by = measure]
fwrite(Pct_emVar_byDistActiveTSS, sprintf("%s/Pct_emVar_by_Dist/01b_Pct_emVar_%s_byDistActiveTSS_byCelltype.tsv", FIGURE_DIR, EMVAR_LIST), sep = "\t", quote = FALSE)

for (i_MEASURE in seq_len(Measure_table[, .N])) {
  MEASURE <- Measure_table[i_MEASURE, measure]
  MEASURE_LABEL <- Measure_table[i_MEASURE, measure_label]
  MEASURE_NULL <- Measure_table[i_MEASURE, null]

  FigData <- Pct_emVar_byDistActiveTSS[measure == MEASURE & N_test > minSNP, ]
  p <- ggplot(FigData)
  p <- p + geom_hline(yintercept = MEASURE_NULL, col = "lightgrey", linetype = "dashed")
  p <- p + geom_pointrange(aes(x = distActiveTSS_bin, y = Pct * 100, ymin = Pct_lo * 100, ymax = Pct_hi * 100, col = celline), size = .3, alpha = .5)
  p <- p + geom_pointrange(data = FigData[Pvalue < 0.01, ], aes(x = distActiveTSS_bin, y = Pct * 100, ymin = Pct_lo * 100, ymax = Pct_hi * 100), col = "red", size = .3, alpha = .5)
  # p <- p + facet_grid(rows=vars(measure)) + ylab("Percentage of alleles") + xlab("Population")
  if (is.na(MEASURE_NULL)) {
    MEASURE_NULL <- Pct_emVar_byDistActiveTSS[measure == MEASURE, .(Pct = median(Pct)), by = celline]
    p <- p + geom_hline(data = MEASURE_NULL, aes(yintercept = Pct * 100), linetype = 2)
  } else {
    p <- p + geom_hline(yintercept = MEASURE_NULL * 100, linetype = 2)
  }
  p <- p + theme_plot(rotate.x = 90) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + guides(col = "none") + xlab("")
  p <- p + scale_color_manual(values = color_celline)
  # p <- p + scale_alpha_manual(values=alpha_TRUEFALSE)
  p <- p + ylab(MEASURE_LABEL) + xlab("Distance to Nearest Active TSS (kb)")
  p <- p + facet_grid(rows = vars(celline))
  p <- p + theme(legend.box = "vertical")
  pdf(sprintf("%s/Pct_emVar_by_Dist/01b_Pct_%s_%s_by_Dist_byCelltype.pdf", FIGURE_DIR, MEASURE, EMVAR_LIST), height = 2.3, width = 1.7)
  print(p)
  dev.off()
}
write('In Fig1a, solid dots correspond to FDR<5%\nIn Fig1b, dots in red correspond to P<0.01', file = sprintf("%s/Pct_emVar_by_Dist/README.md", FIGURE_DIR), append = TRUE)


Pct_emVar_byChromHMM <- get_Pct_emVars(emVARs_annot, split_by = c("ROADMAP_current_celline"), total_by = NULL, emVarlist = EMVAR_LIST)
Pct_emVar_byChromHMM <- make_it_long(Pct_emVar_byChromHMM, split_by = "ROADMAP_current_celline")
Pct_emVar_byChromHMM[N_test > minSNP, FDR := p.adjust(Pvalue, "fdr"), by = measure]
Pct_emVar_byChromHMM[N_test > minSNP][order(Pvalue), head(.SD, 5), by = measure]
Pct_emVar_byChromHMM <- Pct_emVar_byChromHMM[order(measure, ROADMAP_current_celline), ]
Pct_emVar_byChromHMM[, ROADMAP_current_celline := factor(ROADMAP_current_celline, levels = ChromHMM_colors[, paste(V1, V2, sep = "_")])]
fwrite(Pct_emVar_byChromHMM, sprintf("%s/Pct_emVar_by_ChromHMM/01c_Pct_emVar_%s_by_ChromHMM.tsv", FIGURE_DIR, EMVAR_LIST), sep = "\t", quote = FALSE)

for (i_MEASURE in seq_len(Measure_table[, .N])) {
  MEASURE <- Measure_table[i_MEASURE, measure]
  MEASURE_LABEL <- Measure_table[i_MEASURE, measure_label]
  MEASURE_NULL <- Measure_table[i_MEASURE, null]

  FigData <- Pct_emVar_byChromHMM[measure == MEASURE & N_test > minSNP, ]
  p <- ggplot(FigData)
  p <- p + geom_hline(yintercept = MEASURE_NULL, col = "lightgrey", linetype = "dashed")
  p <- p + geom_pointrange(aes(x = ROADMAP_current_celline, y = Pct * 100, ymin = Pct_lo * 100, ymax = Pct_hi * 100, col = ROADMAP_current_celline, alpha = ifelse(FDR < 0.05,'FDR<0.05','ns')), size = .3)
  # p <- p + geom_pointrange(data= FigData[Pvalue<0.01,],aes(x = ROADMAP_current_celline, y = Pct * 100, ymin = Pct_lo * 100, ymax = Pct_hi * 100),col = 'red', size = .3, alpha=.5)
  # p <- p + facet_grid(rows=vars(measure)) + ylab("Percentage of alleles") + xlab("Population")
  if (is.na(MEASURE_NULL)) {
    MEASURE_NULL <- Pct_emVar_byChromHMM[measure == MEASURE, .(Pct = median(Pct))]
    p <- p + geom_hline(data = MEASURE_NULL, aes(yintercept = Pct * 100), linetype = 2)
  } else {
    p <- p + geom_hline(yintercept = MEASURE_NULL * 100, linetype = 2)
  }
  p <- p + theme_plot(rotate.x = 90) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + guides(col = "none") + xlab("")
  p <- p + scale_color_manual(values = setNames(ChromHMM_colors[, HEX], ChromHMM_colors[, paste(V1, V2, sep = "_")]))
  p <- p + scale_alpha_manual(values = setNames(alpha_TRUEFALSE,c('ns','FDR<0.05')))
  p <- p + ylab(MEASURE_LABEL) + xlab("ChromHMM state")
  # p <- p + facet_grid(rows = vars(celline))
  p <- p + theme(legend.box = "vertical")
  pdf(sprintf("%s/Pct_emVar_by_ChromHMM/01c_Pct_%s_%s_by_ChromHMM.pdf", FIGURE_DIR, MEASURE, EMVAR_LIST), height = 2.3, width = 1.7)
  print(p)
  dev.off()
}




minSNP <- 5
Pct_emVar_byChromHMM <- get_Pct_emVars(emVARs_annot, split_by = c("ROADMAP_current_celline", "celline"), total_by = "celline", emVarlist = EMVAR_LIST)
Pct_emVar_byChromHMM <- make_it_long(Pct_emVar_byChromHMM, split_by = c("ROADMAP_current_celline", "celline"))
Pct_emVar_byChromHMM[N_test > minSNP, FDR := p.adjust(Pvalue, "fdr"), by = measure]
Pct_emVar_byChromHMM[N_test > minSNP][order(Pvalue), head(.SD, 5), by = measure]
Pct_emVar_byChromHMM <- Pct_emVar_byChromHMM[order(measure, celline, ROADMAP_current_celline), ]
Pct_emVar_byChromHMM[, ROADMAP_current_celline := factor(ROADMAP_current_celline, levels = ChromHMM_colors[, paste(V1, V2, sep = "_")])]
fwrite(Pct_emVar_byChromHMM, sprintf("%s/Pct_emVar_by_ChromHMM/01d_Pct_emVar_%s_by_ChromHMM_byCelltype.tsv", FIGURE_DIR, EMVAR_LIST), sep = "\t", quote = FALSE)

for (i_MEASURE in seq_len(Measure_table[, .N])) {
  MEASURE <- Measure_table[i_MEASURE, measure]
  MEASURE_LABEL <- Measure_table[i_MEASURE, measure_label]
  MEASURE_NULL <- Measure_table[i_MEASURE, null]

  FigData <- Pct_emVar_byChromHMM[measure == MEASURE & N_test > minSNP, ]
  p <- ggplot(FigData)
  p <- p + geom_hline(yintercept = MEASURE_NULL, col = "lightgrey", linetype = "dashed")
  p <- p + geom_pointrange(aes(x = ROADMAP_current_celline, y = Pct * 100, ymin = Pct_lo * 100, ymax = Pct_hi * 100, col = celline), size = .3, alpha = .5)
  p <- p + geom_pointrange(data = FigData[Pvalue < 0.01, ], aes(x = ROADMAP_current_celline, y = Pct * 100, ymin = Pct_lo * 100, ymax = Pct_hi * 100), col = "red", size = .3, alpha = .5)
  # p <- p + facet_grid(rows=vars(measure)) + ylab("Percentage of alleles") + xlab("Population")
  if (is.na(MEASURE_NULL)) {
    MEASURE_NULL <- Pct_emVar_byChromHMM[measure == MEASURE, .(Pct = median(Pct)), by = celline]
    p <- p + geom_hline(data = MEASURE_NULL, aes(yintercept = Pct * 100), linetype = 2)
  } else {
    p <- p + geom_hline(yintercept = MEASURE_NULL * 100, linetype = 2)
  }
  p <- p + theme_plot(rotate.x = 90) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + guides(col = "none") + xlab("")
  p <- p + scale_color_manual(values = color_celline)
  # p <- p + scale_alpha_manual(values=alpha_TRUEFALSE)
  p <- p + ylab(MEASURE_LABEL) + xlab("ChromHMM state")
  p <- p + facet_grid(rows = vars(celline))
  p <- p + theme(legend.box = "vertical")
  pdf(sprintf("%s/Pct_emVar_by_ChromHMM/01d_Pct_%s_%s_by_ChromHMM_byCelltype.pdf", FIGURE_DIR, MEASURE, EMVAR_LIST), height = 2.3, width = 1.7)
  print(p)
  dev.off()
}
write('In Fig1c, solid dots correspond to FDR<5%\nIn Fig1d, dots in red correspond to P<0.01', file = sprintf("%s/Pct_emVar_by_ChromHMM/README.md", FIGURE_DIR), append = TRUE)


cat("All done\n")
q("no")
  
# emVARs_annot_v2 <- emVARs_annot_v2[method=='NearestTSS',]

# emVARs_annot_v2 <- emVARs_annot_v2[method=='NearestGene',]

# N_crs_diff=emVARs_annot_v2[,.(N_crs=length(unique(crsID)),N_crs_sig=length(unique(crsID[FDR<.05])),N_crs_up=length(unique(crsID[FDR<.05 & log2FC_archaic_vs_modern>.2])),N_crs_down=length(unique(crsID[FDR<.05 & log2FC_archaic_vs_modern< -.2]))),,keyby=cut(-Score1,c(-1,.5,1,2,5,10,50,100,200,500,Inf))]

# N_crs_diff=emVARs_annot_v2[,.(N_crs=length(unique(crsID)),
# 															N_crs_sig=length(unique(crsID[FDR<.05])),
# 															N_crs_up=length(unique(crsID[FDR<.05 & log2FC_archaic_vs_modern>.2])),
# 															N_crs_down=length(unique(crsID[FDR<.05 & log2FC_archaic_vs_modern< -.2])),
# 															N_crs_up_stg=length(unique(crsID[FDR<.05 & log2FC_archaic_vs_modern>.5])),
# 															N_crs_down_stg=length(unique(crsID[FDR<.05 & log2FC_archaic_vs_modern< -.5])),
# 															N_crs_up_wk=length(unique(crsID[FDR<.05 & log2FC_archaic_vs_modern>0])),
# 															N_crs_down_wk=length(unique(crsID[FDR<.05 & log2FC_archaic_vs_modern< 0]))),
# 															keyby=cut(-Score1,c(-1,.5,1,2,5,10,50,100,200,500,Inf))]
# N_crs_diff[,N_crs_perkb:=N_crs/c(0.5,0.5,1,3,5,40,50,100,300,500)]
# N_crs_diff[,Pct_sig:=N_crs_sig/N_crs]
# N_crs_diff[,Pct_down:=N_crs_down/N_crs]
# N_crs_diff[,Pct_up:=N_crs_up/N_crs]
# N_crs_diff[,Pct_down_wk:=N_crs_down_wk/N_crs]
# N_crs_diff[,Pct_up_wk:=N_crs_up_wk/N_crs]
# N_crs_diff[,Pct_down_stg:=N_crs_down_stg/N_crs]
# N_crs_diff[,Pct_up_stg:=N_crs_up_stg/N_crs]


# N_crs_diff_ct=emVARs_annot_v2[,.(N_crs=length(unique(crsID)),
# 															N_crs_sig=length(unique(crsID[FDR<.05])),
# 															N_crs_up=length(unique(crsID[FDR<.05 & log2FC_archaic_vs_modern>.2])),
# 															N_crs_down=length(unique(crsID[FDR<.05 & log2FC_archaic_vs_modern< -.2])),
# 															N_crs_up_stg=length(unique(crsID[FDR<.05 & log2FC_archaic_vs_modern>.5])),
# 															N_crs_down_stg=length(unique(crsID[FDR<.05 & log2FC_archaic_vs_modern< -.5])),
# 															N_crs_up_wk=length(unique(crsID[FDR<.05 & log2FC_archaic_vs_modern>0])),
# 															N_crs_down_wk=length(unique(crsID[FDR<.05 & log2FC_archaic_vs_modern< 0]))),
# 															keyby=.(celltype,cut(-Score1,c(-1,.5,1,2,5,10,50,100,200,500,Inf)))]


# N_crs_diff_ct[,N_crs_perkb:=N_crs/c(0.5,0.5,1,3,5,40,50,100,300,500)]
# N_crs_diff_ct[,Pct_sig:=N_crs_sig/N_crs]
# N_crs_diff_ct[,Pct_down:=N_crs_down/N_crs]
# N_crs_diff_ct[,Pct_up:=N_crs_up/N_crs]
# N_crs_diff_ct[,Pct_down_wk:=N_crs_down_wk/N_crs]
# N_crs_diff_ct[,Pct_up_wk:=N_crs_up_wk/N_crs]
# N_crs_diff_ct[,Pct_down_stg:=N_crs_down_stg/N_crs]
# N_crs_diff_ct[,Pct_up_stg:=N_crs_up_stg/N_crs]

# N_crs_diff_ct=emVARs_annot_v2[,.(N_crs=length(unique(crsID)),
# 															N_crs_sig=length(unique(crsID[FDR<.05])),
# 															N_crs_up=length(unique(crsID[FDR<.05 & log2FC_archaic_vs_modern>.2])),
# 															N_crs_down=length(unique(crsID[FDR<.05 & log2FC_archaic_vs_modern< -.2])),
# 															N_crs_up_stg=length(unique(crsID[FDR<.05 & log2FC_archaic_vs_modern>.5])),
# 															N_crs_down_stg=length(unique(crsID[FDR<.05 & log2FC_archaic_vs_modern< -.5])),
# 															N_crs_up_wk=length(unique(crsID[FDR<.05 & log2FC_archaic_vs_modern>0])),
# 															N_crs_down_wk=length(unique(crsID[FDR<.05 & log2FC_archaic_vs_modern< 0]))),
# 															keyby=.(celltype,cut(-Score1,c(-1,2,10,100,Inf)))]


# N_crs_diff[,N_crs_perkb:=N_crs/c(0.5,0.5,1,3,5,40,50,100,300,500)]
# N_crs_diff_ct[,Pct_sig:=N_crs_sig/N_crs]
# N_crs_diff_ct[,Pct_down:=N_crs_down/N_crs]
# N_crs_diff_ct[,Pct_up:=N_crs_up/N_crs]
# N_crs_diff_ct[,Pct_down_wk:=N_crs_down_wk/N_crs]
# N_crs_diff_ct[,Pct_up_wk:=N_crs_up_wk/N_crs]
# N_crs_diff_ct[,Pct_down_stg:=N_crs_down_stg/N_crs]
# N_crs_diff_ct[,Pct_up_stg:=N_crs_up_stg/N_crs]
