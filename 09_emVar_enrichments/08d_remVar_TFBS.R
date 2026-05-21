

MPRA_DIR <- "/pasteur/helix/projects/evo_immuno_pop/MPRA"
ANALYSIS_DIR <- "MPRA_count_exp6_analysis3"

source(sprintf("%s/scripts/%s/00_starter.R", MPRA_DIR, ANALYSIS_DIR))
source(sprintf("%s/scripts/misc_plots.R", MPRA_DIR))



RUN_ID <- "RUN3_Z2_nBC10"
CRITERIA_ACTIVE <- "FDR5_scrambled5pct_FC0.2_GCnorm"
CRITERIA_EMVARS <- "EmVar_FDR5_FC.2"
CRITERIA_DIFF <- "Diff_FDR5_FC.2"
CRITERIA_EMVARS_DIFF <- "EmVarDiff_FDR5_FC.2"
FILTER_ACTIVE <- TRUE
FILTER_DIFF <- FALSE
FILTER_EMVAR <- TRUE
STIM_ONLY <- FALSE

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
  if (cmd[i] == "--criteria_emvar" || cmd[i] == "-e") {
    CRITERIA_EMVARS <- cmd[i + 1]
  }
  if (cmd[i] == "--criteria_emvar_diff" || cmd[i] == "-ed") {
    CRITERIA_EMVARS_DIFF <- cmd[i + 1]
  }
  if (cmd[i] == "--filter_active" || cmd[i] == "-fa") {
    FILTER_ACTIVE <- as.logical(cmd[i + 1])
  }
  if (cmd[i] == "--filter_diff" || cmd[i] == "-fd") {
    FILTER_DIFF <- as.logical(cmd[i + 1])
  }
	if(cmd[i] == "--filter_emvar" || cmd[i] == "-fe") {
    FILTER_EMVAR <- as.logical(cmd[i + 1])
	}
  if(cmd[i] == "--stim_only" || cmd[i] == "-s") {
    STIM_ONLY <- as.logical(cmd[i + 1])
  }
}

if (!FILTER_ACTIVE) {
  CRITERIA_ACTIVE_OUT <- paste0("noActivityFilter_", CRITERIA_ACTIVE)
} else {
  CRITERIA_ACTIVE_OUT <- CRITERIA_ACTIVE
}

if (!FILTER_DIFF) {
  CRITERIA_DIFF_OUT <- paste0("noDiffFilter_", CRITERIA_DIFF)
} else {
  CRITERIA_DIFF_OUT <- CRITERIA_DIFF
}

if (!FILTER_EMVAR) {
  CRITERIA_EMVARS_OUT <- paste0("noEmVarFilter_", CRITERIA_EMVARS)
} else {
  CRITERIA_EMVARS_OUT <- CRITERIA_EMVARS
}
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
SNP_annot_v5 <- merge(SNP_annot_v4[,-"Introgression_scenario"],selected_annot_wide,by=c('ID','posID'))
toc()


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

FIGURE_DIR <- sprintf("%s/figures/%s/%s/05_emVarDiff/%s/%s/%s/%s/%s", MPRA_DIR, ANALYSIS_DIR, RUN_ID, CRITERIA_EMVARS_DIFF, CRITERIA_EMVARS_OUT, CRITERIA_DIFF_OUT, CRITERIA_ACTIVE_OUT, STIM_ONLY_OUT)

dir.create(FIGURE_DIR, recursive = TRUE)
dir.create(sprintf("%s/SupTables/", FIGURE_DIR), recursive = TRUE)

# load tested oligos
source(sprintf("%s/scripts/%s/02z__define_includedCRS.R", MPRA_DIR, ANALYSIS_DIR))
dir.create(sprintf("%s/remVar_TFBS/", FIGURE_DIR))

####################################################################
################# compute FDR, emVars diff ######################
####################################################################

# emVARs_Diff_obs_response = fread(file = sprintf("%s/all_emVars_diff_obs_response.tsv", EMVAR_DIFF_DIR), sep = "\t")


# emVARs_Diff_any <- emVARs_Diff_obs_response[is_emVar_Diff_CRITERIA == TRUE, .(posID, celline, condition, Introgression_source = Introgression_source_top, Introgression_scenario = Introgression_scenario_top, allele_match)]
# emVars_Diff_any_increased <- emVARs_Diff_obs_response[is_emVar_Diff_CRITERIA == TRUE & log2FC_archaic_vs_modern > 0, .(posID, celline, condition, Introgression_source = Introgression_source_top, Introgression_scenario = Introgression_scenario_top, allele_match)]
# emVars_any_decreased <- emVARs_annot_ctc[is_emVar_CRITERIA == TRUE & log2FC_archaic_vs_modern < 0, .(posID, celline, condition, Introgression_source = Introgression_source_top, Introgression_scenario = Introgression_scenario_top, allele_match)]

# emVars_NS <- emVARs_annot_ctc[is_emVar_CRITERIA == TRUE & condition == "NS", .(posID, celline, condition, Introgression_source = Introgression_source_top, Introgression_scenario = Introgression_scenario_top, allele_match)]
# emVars_NS_increased <- emVARs_annot_ctc[is_emVar_CRITERIA == TRUE & condition == "NS" & log2FC_archaic_vs_modern > 0, .(posID, celline, condition, Introgression_source = Introgression_source_top, Introgression_scenario = Introgression_scenario_top, allele_match)]
# emVars_NS_decreased <- emVARs_annot_ctc[is_emVar_CRITERIA == TRUE & condition == "NS" & log2FC_archaic_vs_modern < 0, .(posID, celline, condition, Introgression_source = Introgression_source_top, Introgression_scenario = Introgression_scenario_top, allele_match)]


# oligo_response <- merge(oligo_response, SNP_annot_v5[, .(crsID, POP_adaptive, Introgressed_from = Adaptive_from, Introgression_scenario)], by = "crsID")
# oligo_response_any <- oligo_response[FDR < .05 & abs(log2FC_2vs1) > logFC_TH & excluded == FALSE & ctrl == FALSE, .(crsID, celline, condition, type, POP_adaptive, Introgressed_from)]
# oligo_response_increased <- oligo_response[FDR < .05 & log2FC_2vs1 > logFC_TH & excluded == FALSE & ctrl == FALSE, .(crsID, celline, condition, type, POP_adaptive, Introgressed_from)]
# oligo_response_decreased <- oligo_response[FDR < .05 & log2FC_2vs1 < -logFC_TH & excluded == FALSE & ctrl == FALSE, .(crsID, celline, condition, type, POP_adaptive, Introgressed_from)]

# N_Context <- oligo_response_any[, .(nContext = .N), by = crsID][, .N, keyby = nContext][, Pct := N / sum(N)][1:.N]
emVARs_annot_ctc <- fread(sprintf("%s/all_emVars_annotated_celltype__%s.tsv", EMVAR_DIR, CRITERIA_ACTIVE_OUT))

emVars_any <- emVARs_annot_ctc[is_emVar_CRITERIA == TRUE, .(crsID=posID, ID, celline, condition, Introgression_source = Introgression_source_top, Introgression_scenario = Introgression_scenario_top, allele_match)]
emVars_any_suggestive <- emVARs_annot_ctc[crsID %in% emVars_any[,crsID]  & pval_emp<0.05 , .(crsID, ID, celline, condition, Introgression_source = Introgression_source_top, Introgression_scenario = Introgression_scenario_top, allele_match)]

response_emVars_obs <- fread(sprintf("%s/all_emVars_diff_obs_response.tsv", EMVARDIFF_DIR))
response_emVars_any <- response_emVars_obs[is_emVar_Diff_CRITERIA==TRUE, .(crsID, ID,celline, condition, Introgression_source = Introgression_source_top, Introgression_scenario = Introgression_scenario_top, allele_match)]
response_emVars_stim <- response_emVars_obs[is_emvar_cre_group2 & !is_emvar_cre_group1 & is_emVar_Diff_CRITERIA, .(crsID, ID, celline, condition, Introgression_source = Introgression_source_top, Introgression_scenario = Introgression_scenario_top, allele_match)]
response_emVars_ns <- response_emVars_obs[is_emvar_cre_group1 & !is_emvar_cre_group2 & is_emVar_Diff_CRITERIA, .(crsID, ID, celline, condition, Introgression_source = Introgression_source_top, Introgression_scenario = Introgression_scenario_top, allele_match)]
#response_emVars_any <- merge(response_emVars_any, condition_summary, by = "COND_ID")

response_emVars_any_suggestive <- response_emVars_obs[crsID %in% response_emVars_any[,crsID] & pval_emp<0.05, .(crsID, ID, celline, condition, Introgression_source = Introgression_source_top, Introgression_scenario = Introgression_scenario_top, allele_match)]

TFBM_SNP_overlap=fread(sprintf("%s/data/%s/00_TFBM/1a_TFBM_diff_all.tsv.gz", MPRA_DIR, ANALYSIS_DIR))
TFBM_SNP_overlap = TFBM_SNP_overlap[crsID %chin% SNPs_with_more_than_one_introgressed_allele,]

TFBM_SNP_overlap[,is_remVar:=crsID %chin% response_emVars_any_suggestive[,crsID]]
TFBM_SNP_overlap[,is_emVar:=crsID %chin% emVars_any[,crsID]]

TFBM_SNP_overlap[,DIFF_AvsM_pct:=case_when(!is_introgressed~NA,
                                    allele.2==INTROGRESSED.allele~(exp(A2)-exp(A1))/exp(A1),
                                    allele.1==INTROGRESSED.allele~(exp(A1)-exp(A2))/exp(A2),
                                    TRUE~Inf)]

best_TFBM_diff= TFBM_SNP_overlap[pmax(A1,A2)>5 & !is.infinite(DIFF_AvsM),][order(-abs(DIFF_AvsM)),head(.SD,1),by=crsID]

getScoresTFBM <- function(snp_list_1, snp_list_2, TF_scores=TFBM_SNP_overlap[,.(snps=crsID,TF_score=DIFF_AvsM)]){
  DT1=data.table(snps=snp_list_1)
  DT1= merge(DT1,TF_scores,by='snps',all.x=TRUE)
  DT1[ is.na(TF_score),TF_score:=0]

  DT2=data.table(snps=snp_list_2)
  DT2= merge(DT2,TF_scores,by='snps',all.x=TRUE)
  DT2[ is.na(TF_score),TF_score:=0]
  
  results <- list(N_1=DT1[,.N],N_tfbs_1=DT1[,sum(TF_score!=0)],median_1=DT1[TF_score!=0,median(TF_score)],  N_2=DT2[,.N],N_tfbs_2=DT2[,sum(TF_score!=0)], median_2=DT2[TF_score!=0,median(TF_score)], P_wilcox=wilcox.test(DT1[,TF_score],DT2[,TF_score])$p.value)
  return(results)
}

TFBM_SNP_overlap_AvsMpct_comparison_remVar_vs_emVar <- TFBM_SNP_overlap[,getScoresTFBM(unique(response_emVars_any_suggestive$crsID),
                                                                               setdiff(emVars_any_suggestive$crsID,unique(response_emVars_any_suggestive$crsID)),
                                                                               .SD[,.(snps=crsID,TF_score=DIFF_AvsM_pct)]),
                    by=.(matrix_id,names)]

TFBM_SNP_overlap_AvsMpct_comparison_remVar_vs_emVar[,FDR:=p.adjust(P_wilcox,method='fdr')]

TFBM_SNP_overlap_AvsM_comparison_remVar_vs_emVar <- TFBM_SNP_overlap[,getScoresTFBM(unique(response_emVars_any_suggestive$crsID),
                                                                               setdiff(emVars_any_suggestive$crsID,unique(response_emVars_any_suggestive$crsID)),
                                                                               .SD[,.(snps=crsID,TF_score=DIFF_AvsM)]),
                    by=.(matrix_id,names)]
TFBM_SNP_overlap_AvsM_comparison_remVar_vs_emVar[,FDR:=p.adjust(P_wilcox,method='fdr')]


TFBM_SNP_overlap_absAvsM_comparison_remVar_vs_emVar <- TFBM_SNP_overlap[,getScoresTFBM(unique(response_emVars_any_suggestive$crsID),
                                                                               setdiff(emVars_any_suggestive$crsID,unique(response_emVars_any_suggestive$crsID)),
                                                                               .SD[,.(snps=crsID,TF_score=abs(DIFF_AvsM))]),
                    by=.(matrix_id,names)]
TFBM_SNP_overlap_absAvsM_comparison_remVar_vs_emVar[,FDR:=p.adjust(P_wilcox,method='fdr')]
fwrite(TFBM_SNP_overlap_absAvsM_comparison_emVar_vs_all, sprintf("%s/remVar_TFBS/TFBM_SNP_overlap_absAvsM_comparison_emVar_vs_all.tsv", FIGURE_DIR), sep = "\t", quote = FALSE)



######
TFBM_SNP_overlap_absAvsM_comparison_emVar_vs_all=list()

for(i in seq_len(condition_summary[,.N])){
  myCL=condition_summary[i,celline]
  myCOND=condition_summary[i,condition]
  myCOND_ID=condition_summary[i,COND_ID]
  cat(sprintf("%s %s\n",myCL,myCOND))
  TFBM_SNP_overlap_absAvsM_comparison_emVar_vs_all[[myCOND_ID]] <- TFBM_SNP_overlap[,getScoresTFBM(emVars_any_suggestive[celline==myCL & condition==myCOND,crsID],
                                                                               setdiff(unique(emVARs_annot_ctc$crsID),emVars_any_suggestive[celline==myCL & condition==myCOND,crsID]),
                                                                               .SD[,.(snps=crsID,TF_score=abs(DIFF_AvsM))]),
                    by=.(matrix_id,names)]
}
TFBM_SNP_overlap_absAvsM_comparison_emVar_vs_all=rbindlist(TFBM_SNP_overlap_absAvsM_comparison_emVar_vs_all,idcol='COND_ID')
TFBM_SNP_overlap_absAvsM_comparison_emVar_vs_all[,FDR:=p.adjust(P_wilcox,method='fdr'),by=COND_ID]

dir.create(sprintf("%s/remVar_TFBS", FIGURE_DIR), recursive = TRUE)
fwrite(TFBM_SNP_overlap_absAvsM_comparison_emVar_vs_all, sprintf("%s/remVar_TFBS/TFBM_SNP_overlap_absAvsM_comparison_emVar_vs_all.tsv", FIGURE_DIR), sep = "\t", quote = FALSE)


TFBM_SNP_overlap_absAvsM_comparison_remVar_vs_all=list()

for(i in seq_len(condition_summary[condition!='NS',.N])){
  myCL=condition_summary[condition!='NS',][i,celline]
  myCOND=condition_summary[condition!='NS',][i,condition]
  myCOND_ID=condition_summary[condition!='NS',][i,COND_ID]
  cat(sprintf("%s %s\n",myCL,myCOND))
  my_remVars=intersect(emVars_any_suggestive[celline==myCL & condition==myCOND,crsID],response_emVars_any_suggestive[celline==myCL & condition==myCOND,crsID])
  TFBM_SNP_overlap_absAvsM_comparison_remVar_vs_all[[myCOND_ID]] <- TFBM_SNP_overlap[,getScoresTFBM(my_remVars,
                                                                               setdiff(unique(emVARs_annot_ctc$crsID),my_remVars),
                                                                               .SD[,.(snps=crsID,TF_score=abs(DIFF_AvsM))]),
                    by=.(matrix_id,names)]
}
TFBM_SNP_overlap_absAvsM_comparison_remVar_vs_all=rbindlist(TFBM_SNP_overlap_absAvsM_comparison_remVar_vs_all,idcol='COND_ID')
TFBM_SNP_overlap_absAvsM_comparison_remVar_vs_all[,FDR:=p.adjust(P_wilcox,method='fdr'),by=COND_ID]

dir.create(sprintf("%s/remVar_TFBS", FIGURE_DIR), recursive = TRUE)
fwrite(TFBM_SNP_overlap_absAvsM_comparison_remVar_vs_all, sprintf("%s/remVar_TFBS/TFBM_SNP_overlap_absAvsM_comparison_remVar_vs_all.tsv", FIGURE_DIR), sep = "\t", quote = FALSE)



feat_Zdiff=fread(sprintf('%s/data/%s/00_DeepSEA__effect_5682_emVar_testedSNP__Beluga__dfcc57f5-f2c1-4a35-812a-56477b236645/DeepSEA__5682snps__Beluga__effects.tsv.gz',MPRA_DIR, ANALYSIS_DIR))
# feat_Zdiff[grepl('H3K27ac',variable) & grepl('K562|HepG2|A549',variable),] 
# TODO : compare this with emVars systematically 

my_emVars=emVARs_annot_ctc[is_emVar_CRITERIA==TRUE,unique(ID)]

emVar_beluga=merge(emVARs_annot_ctc,feat_Zdiff[assay_simple%in%c('H3K27ac','Accessibility') & grepl('K562|HepG2|A549',variable),.(ID,assay_simple,variable,Z_introgressed,diff_introgressed)],by='ID',allow.cartesian=TRUE)
emVar_beluga[,same_celltype:=grepl(unique(celline),variable),by=celline]
emVar_beluga[!ID %chin% my_emVars & condition=='NS',any(abs(Z_introgressed)>2),by=.(ID,assay_simple)][,mean(V1),keyby=.(assay_simple)]
# [1] 0.1149571
emVar_beluga[ID %chin% my_emVars & condition=='NS',any(abs(Z_introgressed)>2),by=.(ID,assay_simple)][,mean(V1),keyby=.(assay_simple)]
# [1] 0.2278665

emVar_beluga[!ID %chin% my_emVars & condition=='NS' & same_celltype,any(abs(Z_introgressed)>3),by=.(ID,assay_simple)][,mean(V1),keyby=.(assay_simple)]
# [1] 0.06668784
emVar_beluga[ID %chin% my_emVars & condition=='NS' & same_celltype,any(abs(Z_introgressed)>3),by=.(ID,assay_simple)][,mean(V1),keyby=.(assay_simple)]
# [1] 0.1465893
emVar_beluga[condition=='NS' & same_celltype, .(Beluga=any(abs(Z_introgressed)>3)),by=.(ID,assay_simple,celline)][,unlist(fisher.test(ID %chin% my_emVars,Beluga)[c('estimate','p.value')]),keyby=.(assay_simple,celline)]

emVar_beluga[ID %chin% my_emVars & condition=='NS' & same_celltype & is_emVar_CRITERIA,any(abs(Z_introgressed)>3),by=.(ID, celline)][,mean(V1),by= celline]
#    celline        V1
# 1:   HepG2 0.1280000
# 2:    A549 0.2258065
# 3:    K562 0.2692308
#### 13-27% of detectable emVars (FDR<5%) are predicted by beluga to impact chromatin accessibility or H3K27Ac in the same cell line (|Z| > 3)
emVar_beluga[!ID %chin% my_emVars & condition=='NS' & same_celltype,any(abs(Z_introgressed)>3),by=.(ID, celline)][,mean(V1),by= celline]
#   celline         V1
# 1:    A549 0.03222719
# 2:   HepG2 0.04668149
# 3:    K562 0.04636393
#### This percentage falls to <4.6% for variants that are not emVars


emVar_beluga[same_celltype & condition=='NS' & abs(Z_introgressed)>3,mean(is_emVar_CRITERIA),by=celline]
#    celline        V1
# 1:    K562 0.1773256
# 2:    A549 0.1340996
# 3:   HepG2 0.1771772
#### 13-18% of variants predicted by beluga to impact chromatin accessibility of H3K27Ac (|Z| > 3) have detectable emVars in the same cell line

emVar_beluga[same_celltype & condition=='NS' & abs(Z_introgressed)<3,mean(is_emVar_CRITERIA),by=celline]
#    celline         V1
# 1:    A549 0.01348335
# 2:   HepG2 0.06180127
# 3:    K562 0.02247090
#### This percentage falls to 1-6% for variants that are not predicted to as functional by beluga (|Z| < 3)



emVar_beluga[same_celltype & condition=='NS' & abs(Z_introgressed)<2,mean(is_emVar_CRITERIA),by=celline]


emVar_beluga[same_celltype & condition=='NS' & abs(Z_introgressed)>3,mean(pval_emp<0.05),by=celline]
#    celline        V1
# 1:    K562 0.1616679	
# 2:    A549 0.1331418
# 3:   HepG2 0.1565087
#### 13-16% of variants predicted by beluga to impact chromatin accessibility of H3K27Ac have detectable emVars in the same cell line
emVar_beluga[same_celltype & abs(Z_introgressed)<3,mean(pval_emp<0.05),by=celline]

p <- ggplot(emVar_beluga[(abs(Z_introgressed)>2 & (pval_emp<0.05 & ID %chin% my_emVars)) & condition=='NS',], aes(x=log2FC_archaic_vs_modern,y=diff_introgressed,col=COND_ID)) + geom_point() + facet_grid(variable~celline)
p <- p + scale_color_manual(values=color_setup_simplified_norep) + theme_plot()

pdf(sprintf('%s/remVar_TFBS/compared_effect_size.pdf',FIGURE_DIR))
print(p)
dev.off()

# p <- ggplot(emVar_beluga[(abs(Z_introgressed)>2 &(pval_emp<0.05 & ID %chin% my_emVars)) & condition=='NS',], aes(x=log2FC_archaic_vs_modern,y=Z_introgressed,col=COND_ID)) + geom_point() + facet_grid(variable~celline)
# p <- p + scale_color_manual(values=color_setup_simplified_norep) + theme_plot()

# pdf(sprintf('%s/remVar_TFBS/compared_effect_size_Z.pdf',FIGURE_DIR))
# print(p)
# dev.off()

Beluga_diff_comparison_emVar_vs_all=list()
Beluga_diff_vs_emVar_correlation=list()

Test_cor_emVar <- function(snplist, emVar_table, TF_scores){
		emVar_TF=merge(emVar_table[ID%in%snplist,.(ID,beta)],TF_scores[,.(ID,diff_score)],by='ID')
		emVar_TF[,.(r=cor(beta,diff_score),rho=cor(beta,diff_score,method='s'),p.cor=cor.test(beta,diff_score)$p.value,p.rho=cor.test(beta,diff_score,method='s')$p.value)]
}

for(i in seq_len(condition_summary[,.N])){
  myCL=condition_summary[i,celline]
  myCOND=condition_summary[i,condition]
  myCOND_ID=condition_summary[i,COND_ID]
  cat(sprintf("%s %s\n",myCL,myCOND))
  my_emVars=emVars_any_suggestive[celline==myCL & condition==myCOND,ID]
  Beluga_diff_comparison_emVar_vs_all[[myCOND_ID]] <- feat_Zdiff[,getScoresTFBM(my_emVars,
                                                                               setdiff(unique(emVARs_annot_ctc$ID),my_emVars),
                                                                               .SD[,.(snps=ID,TF_score=abs(diff_introgressed))]),
                    by=.(variable)]
	Beluga_diff_vs_emVar_correlation[[myCOND_ID]] <- feat_Zdiff[ID %chin% my_emVars,Test_cor_emVar(my_emVars, 
																							emVARs_annot_ctc[celline==myCL & condition==myCOND , .(ID, beta=log2FC_archaic_vs_modern)],
																								 .SD[,.(ID, diff_score = diff_introgressed)] )
																								 ,by=variable]
}
Beluga_diff_comparison_emVar_vs_all=rbindlist(Beluga_diff_comparison_emVar_vs_all,idcol='COND_ID')
Beluga_diff_comparison_emVar_vs_all[,FDR:=p.adjust(P_wilcox,method='fdr'),by=COND_ID]


Beluga_diff_vs_emVar_correlation=rbindlist(Beluga_diff_vs_emVar_correlation,idcol='COND_ID')
Beluga_diff_vs_emVar_correlation[,FDR.rho:=p.adjust(p.rho,method='fdr'),by=COND_ID]
Beluga_diff_vs_emVar_correlation[,FDR.pearson:=p.adjust(p.r,method='fdr'),by=COND_ID]
fwrite(Beluga_diff_vs_emVar_correlation, sprintf("%s/remVar_TFBS/Beluga_diff_vs_emVar_correlation.tsv", FIGURE_DIR), sep = "\t", quote = FALSE)

Beluga_diff_comparison_emVar_vs_all=merge(Beluga_diff_comparison_emVar_vs_all,Beluga_diff_vs_emVar_correlation,by=c('COND_ID','variable'))
colnames(Beluga_diff_comparison_emVar_vs_all)=c('COND_ID','predicted_epigenetic_mark','N_emVars','N_emVars_withmark','median_emVars','N_non_emVars','N_non_emVars_withmark','median_non_emVars','P_wilcox','FDR','Pearson_cor_effect_size', 'Spearman_cor_effect_size', 'Pearson_pval', 'Spearman_pval', 'Spearman_FDR')
fwrite(Beluga_diff_comparison_emVar_vs_all, sprintf("%s/remVar_TFBS/Beluga_diff_comparison_emVar_vs_all.tsv", FIGURE_DIR), sep = "\t", quote = FALSE)


for(i in seq_len(condition_summary[condition!='NS',.N])){
  myCL=condition_summary[condition!='NS',][i,celline]
  myCOND=condition_summary[condition!='NS',][i,condition]
  myCOND_ID=condition_summary[condition!='NS',][i,COND_ID]
  cat(sprintf("%s %s\n",myCL,myCOND))
  my_remVars=intersect(emVars_any_suggestive[celline==myCL & condition==myCOND,ID],response_emVars_any_suggestive[celline==myCL & condition==myCOND,ID])
  Beluga_diff_comparison_remVar_vs_all[[myCOND_ID]] <- feat_Zdiff[,getScoresTFBM(my_remVars,
                                                                               setdiff(unique(emVARs_annot_ctc$ID),my_remVars),
                                                                               .SD[,.(snps=ID,TF_score=abs(diff_introgressed))]),
                    by=.(variable)]
}
Beluga_diff_comparison_remVar_vs_all=rbindlist(Beluga_diff_comparison_remVar_vs_all,idcol='COND_ID')
Beluga_diff_comparison_remVar_vs_all[,FDR:=p.adjust(P_wilcox,method='fdr'),by=COND_ID]
fwrite(Beluga_diff_comparison_remVar_vs_all, sprintf("%s/remVar_TFBS/Beluga_diff_comparison_remVar_vs_all.tsv", FIGURE_DIR), sep = "\t", quote = FALSE)




Beluga_Z_comparison_remVar_vs_all=list()

for(i in seq_len(condition_summary[condition!='NS',.N])){
  myCL=condition_summary[condition!='NS',][i,celline]
  myCOND=condition_summary[condition!='NS',][i,condition]
  myCOND_ID=condition_summary[condition!='NS',][i,COND_ID]
  cat(sprintf("%s %s\n",myCL,myCOND))
  my_remVars=intersect(emVars_any_suggestive[celline==myCL & condition==myCOND,ID],response_emVars_any_suggestive[celline==myCL & condition==myCOND,ID])
  Beluga_Z_comparison_remVar_vs_all[[myCOND_ID]] <- feat_Zdiff[,getScoresTFBM(my_remVars,
                                                                               setdiff(unique(emVARs_annot_ctc$ID),my_remVars),
                                                                               .SD[,.(snps=ID,TF_score=abs(Z_introgressed))]),
                    by=.(variable)]
}
Beluga_Z_comparison_remVar_vs_all=rbindlist(Beluga_Z_comparison_remVar_vs_all,idcol='COND_ID')
Beluga_Z_comparison_remVar_vs_all[,FDR:=p.adjust(P_wilcox,method='fdr'),by=COND_ID]
fwrite(Beluga_Z_comparison_remVar_vs_all, sprintf("%s/remVar_TFBS/Beluga_Z_comparison_remVar_vs_all.tsv", FIGURE_DIR), sep = "\t", quote = FALSE)
