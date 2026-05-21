# note : requires both filtered and non filtered emVars have been run.
MPRA_DIR <- "/pasteur/helix/projects/evo_immuno_pop/MPRA"
ANALYSIS_DIR <- "MPRA_count_exp6_analysis3"

source(sprintf("%s/scripts/%s/00_starter.R", MPRA_DIR, ANALYSIS_DIR))
source(sprintf("%s/scripts/%s/00_load_associations.R", MPRA_DIR, ANALYSIS_DIR))
source(sprintf("%s/scripts/misc_plots.R", MPRA_DIR))

RUN_ID <- "RUN3_Z2_nBC10"
#CRITERIA_ACTIVE <- "FDR5_scrambled5pct_FC0.2_GCnorm"
CRITERIA_ACTIVE <- "FDR5_FC1_FC0.2_GCnorm"
CRITERIA_EMVARS <- "EmVar_FDR5_FC.2"

cmd <- commandArgs(trailingOnly = TRUE)
print(cmd)

for (i in seq_along(cmd)) {
  if (cmd[i] == "--run_dir" || cmd[i] == "-r") {
    RUN_ID <- cmd[i + 1]
  }
  if (cmd[i] == "--criteria_active" || cmd[i] == "-a") {
    CRITERIA_ACTIVE <- cmd[i + 1]
  }
  if (cmd[i] == "--filter_active" || cmd[i] == "-f") {
    FILTER_ACTIVE <- as.logical(cmd[i + 1])
  }
  if (cmd[i] == "--criteria_emvar" || cmd[i] == "-e") {
    CRITERIA_EMVARS <- cmd[i + 1]
  }
}


CRITERIA_ACTIVE_NOFILTER <- paste0("noActivityFilter_", CRITERIA_ACTIVE)
#CRITERIA_ACTIVE_OUT <- paste0("merged_", CRITERIA_ACTIVE)
CRITERIA_ACTIVE_OUT <- CRITERIA_ACTIVE

tic("loading oligo & SNP annotations")
# load annotation of oligos (beforz: CRS)
oligo_source <- fread(sprintf("%s/data/%s/00_oligo_annot_v2.txt", MPRA_DIR, ANALYSIS_DIR))

# load annotations of SNPs
SNP_annot_v4 <- fread(sprintf("%s/data/%s/00_SNP_annot_v4.txt", MPRA_DIR, ANALYSIS_DIR))
SNP_annot <- SNP_annot_v4
# load annotation per POP
selected_annot <- fread(sprintf("%s/data/%s/selected_snp_annotation.tsv.gz", MPRA_DIR, ANALYSIS_DIR))

# load clean annotation of introgression
selected_annot_wide <- fread(sprintf("%s/data/%s/selected_snp_annotation_wide.tsv.gz", MPRA_DIR, ANALYSIS_DIR))
toc()

##### active parameters
source(sprintf("%s/scripts/%s/03_00_parameter_activity.R", MPRA_DIR, ANALYSIS_DIR))
##### emVars parameters
source(sprintf("%s/scripts/%s/05_00_parameter_emVars.R", MPRA_DIR, ANALYSIS_DIR))

IN_DIR <- sprintf("%s/data/%s/03a_aggMPRA_analyse/%s/", MPRA_DIR, ANALYSIS_DIR, RUN_ID)
ACTIVE_DIR <- sprintf("%s/data/%s/03b_computeFDR/%s/Activity/%s/", MPRA_DIR, ANALYSIS_DIR, RUN_ID, CRITERIA_ACTIVE)
EMVAR_DIR <- sprintf("%s/data/%s/03b_computeFDR/%s/emVars/%s/%s", MPRA_DIR, ANALYSIS_DIR, RUN_ID, CRITERIA_EMVARS, CRITERIA_ACTIVE)
EMVAR_NOFILTER_DIR <- sprintf("%s/data/%s/03b_computeFDR/%s/emVars/%s/%s", MPRA_DIR, ANALYSIS_DIR, RUN_ID, CRITERIA_EMVARS, CRITERIA_ACTIVE_NOFILTER)


dir.create(EMVAR_DIR, recursive = TRUE)
FIGURE_DIR <- sprintf("%s/figures/%s/%s/04_emVars/%s/%s", MPRA_DIR, ANALYSIS_DIR, RUN_ID, CRITERIA_EMVARS, CRITERIA_ACTIVE_OUT)
#FIGURE_DIR <- sprintf("%s/figures/%s/04_emVars/%s/%s", MPRA_DIR, ANALYSIS_DIR, CRITERIA_EMVARS, CRITERIA_ACTIVE_OUT)

dir.create(FIGURE_DIR, recursive = TRUE)
dir.create(sprintf("%s/SupTables/", FIGURE_DIR), recursive = TRUE)
cat("\nprinting output in:", FIGURE_DIR)


# load tested oligos
source(sprintf("%s/scripts/%s/02z__define_includedCRS.R", MPRA_DIR, ANALYSIS_DIR))

library(ggrepel)

########## load emVar data  ##########

# unfiltered
emVARs_annot_noFilter <- fread(file = sprintf("%s/all_emVars_annotated__%s.tsv.gz", EMVAR_NOFILTER_DIR, CRITERIA_ACTIVE_NOFILTER), sep = "\t")
emVARs_annot_noFilter <- emVARs_annot_noFilter[ANALYSIS_SUBTYPE == "celltype_cond"]
emVARs_annot_noFilter <- merge(emVARs_annot_noFilter, condition_summary[, -"celltype"], by.x = "ANALYSIS_NAME", by.y = "analysis_name")

emVARs_annot_noFilter[is_emVar_CRITERIA==TRUE,length(unique(crsID))]
CRE_emVars_unfiltered <-emVARs_annot_noFilter[is_emVar_CRITERIA==TRUE,unique(crsID)]
length(CRE_emVars_unfiltered)
# 628
active_CRE_emVars_unfiltered <-emVARs_annot_noFilter[is_emVar_CRITERIA==TRUE & cre_class_CRITERIA!='inactive',unique(crsID)]
length(active_CRE_emVars_unfiltered)
# 421
inactive_CRE_emVars_unfiltered <- setdiff(CRE_emVars_unfiltered, active_CRE_emVars_unfiltered)
length(inactive_CRE_emVars_unfiltered)
# 207

emVARs_annot_noFilter[is_emVar_CRITERIA==TRUE,length(unique(crsID)),by=cre_class_CRITERIA]
# WARNING: one emVars may be counted multiple times here
#    cre_class_CRITERIA    V1
#                <char> <int>
# 1:           silencer   164
# 2:    strong enhancer   102
# 3:           enhancer   192
# 4:           inactive   284
# 5:    strong silencer    18

####### double check activity of inactive CRE that contain an emVars

oligo_activity_obs <- fread(sprintf("%s/all_oligos_annotated__%s.tsv.gz", ACTIVE_DIR, CRITERIA_ACTIVE))

emVARs_inactive_CRE <- emVARs_annot_noFilter[ is_emVar_CRITERIA==TRUE & crsID%in%inactive_CRE_emVars_unfiltered,.(crsID,celline,condition)]
oligo_activity_emVars_inactiveCRE <- 	 paste(crsID, celline, condition)%in%emVARs_inactive_CRE[,paste(crsID, celline, condition)],]

oligo_activity_emVars_inactiveCRE[,.(suggestive_activity=any(pval_emp<.05)),by=crsID][,mean(suggestive_activity)]
# [1] 0.52657
oligo_activity_emVars_inactiveCRE[,.(suggestive_activity=any(pval_emp<.05)),by=crsID][,sum(suggestive_activity)]
# [1] 109

# activity filtered
emVARs_annot_obs <- fread(file = sprintf("%s/all_emVars_annotated__%s.tsv.gz", EMVAR_DIR, CRITERIA_ACTIVE), sep = "\t")
emVARs_annot_obs_ctc <- emVARs_annot_obs[ANALYSIS_SUBTYPE == "celltype_cond"]
emVARs_annot_obs_ctc <- merge(emVARs_annot_obs_ctc, condition_summary[, -"celltype"], by.x = "ANALYSIS_NAME", by.y = "analysis_name")

CRE_emVars_filtered <-emVARs_annot_obs_ctc[is_emVar_CRITERIA==TRUE,unique(crsID)]
length(CRE_emVars_filtered)
# 689
length(CRE_emVars_filtered)/length(active_CRE_emVars_unfiltered)
# 1.63658

emVARs_annot_obs_ctc[is_emVar_CRITERIA==TRUE,mean(log2FC_archaic_vs_modern>0)]
# [1] 0.4970258

emVARs_annot_obs_ctc[is_emVar_CRITERIA==TRUE,length(unique(crsID)),by=cre_class_CRITERIA]
#    cre_class_CRITERIA    V1
#                <char> <int>
# 1:           silencer   291
# 2:    strong enhancer   146
# 3:           enhancer   301
# 4:    strong silencer    34

emVARs_annot_obs_ctc[is_emVar_CRITERIA==TRUE & abs(log2FC_archaic_vs_modern)>0.5 ,length(unique(crsID))]
# 90
emVARs_annot_obs_ctc[is_emVar_CRITERIA==TRUE & abs(log2FC_archaic_vs_modern)>0.5 ,][order(CHROM,POS_b37,-abs(log2FC_archaic_vs_modern))][!duplicated(posID),][,sum(c(0,diff(POS_b37))>0 & c(0,diff(POS_b37))<1e5)]
# 21
# 21 snps are <100kb from their predecessor, so 69 independent signals

replicate(100,{emVARs_annot_obs_ctc[!duplicated(posID),][sample(1:.N,90) ,][order(CHROM,POS_b37)][,sum(c(0,diff(POS_b37))>0 & c(0,diff(POS_b37))<1e5)]})
# we are well within the range of what is expected by chance, so we cannot conclude that emVARs with high log2FC are more clustered than expected by chance

emVARs_annot_obs_ctc[is_emVar_CRITERIA==TRUE & abs(log2FC_archaic_vs_modern)>1 ,length(unique(crsID))]
#9
emVARs_annot_obs_ctc[is_emVar_CRITERIA==TRUE & log2FC_archaic_vs_modern>0.5 ,length(unique(crsID))]
#50
emVARs_annot_obs_ctc[is_emVar_CRITERIA==TRUE & log2FC_archaic_vs_modern< -0.5 ,length(unique(crsID))]
#41 
emVARs_annot_obs_ctc[is_emVar_CRITERIA==TRUE & abs(log2FC_archaic_vs_modern)>1 ,length(unique(crsID))]

count_strong_emVars_by_POP <- emVARs_annot_obs_ctc[,.(N=length(unique(ID))),keyby=.(is_emVar_CRITERIA == TRUE & abs(log2FC_archaic_vs_modern) > 0.5, POP_adaptive_top)]
# is_emVar_CRITERIA POP_adaptive_top    N
#                <lgcl>           <char> <int>
#  1:             FALSE                    172
#  2:             FALSE             Agta   871
#  3:             FALSE              CHB   397
#  4:             FALSE              GBR   228
#  5:             FALSE              IBS   199
#  6:             FALSE              JPT   178
#  7:             FALSE              PJL   544
#  8:             FALSE           Papuan   840
#  9:             FALSE              STU   408
# 10:              TRUE                      1
# 11:              TRUE             Agta    22
# 12:              TRUE              CHB     8
# 13:              TRUE              GBR     3
# 14:              TRUE              IBS     5
# 15:              TRUE              JPT     1
# 16:              TRUE              PJL    15
# 17:              TRUE           Papuan    23
# 18:              TRUE              STU    12
count_strong_emVars_in_underRepresented_POP <- count_strong_emVars_by_POP[,.(N=sum(N)),keyby=.(is_emVar_CRITERIA == TRUE   , POP_adaptive_top%in% c("Papuan", "Agta", "PJL", "STU"))]

count_strong_emVars_in_underRepresented_POP[,fisher.test(matrix(N,2,2))]


# one crsID is strong in both directions
intersect(emVARs_annot_obs_ctc[is_emVar_CRITERIA==TRUE & log2FC_archaic_vs_modern< -0.5,crsID], emVARs_annot_obs_ctc[is_emVar_CRITERIA==TRUE & log2FC_archaic_vs_modern> 0.5,crsID])
# "1:57365276" (C8A)

CS_results_finemap_eQTL_tested <- fread(sprintf("%s/data/%s/CredibleSets/CredibleSets_eQTLs.tsv.gz", MPRA_DIR, ANALYSIS_DIR))
CS_results_finemap_eQTL_tested <- merge(SNP_annot_v4[, .(ID, variantId_hg38, crsID, rsID, REF, ALT, INTROGRESSED.allele, allele2_is_REF)], CS_results_finemap_eQTL_tested, all.x = TRUE, by.x = "variantId_hg38", by.y = "variantId")



########################################################################################################
######################################## eQTL catalog  only ############################################
########################################################################################################

eQTL_list <- CS_results_finemap_eQTL_tested[(studyType == "eqtl" | studyType == "sceqtl") & eQTL_type == 'ge', ][order(variantId_hg38, eQTL_gene, eQTL_tissue, pValueExponent)]
eQTL_list[, INTROGRESSED_BETA := ifelse(INTROGRESSED.allele == ALT, beta, ifelse(INTROGRESSED.allele == REF, -beta, NA))]
eQTL_list <- eQTL_list[!is.na(INTROGRESSED_BETA), ]

# doesn't change the pattern
#eQTL_list <- eQTL_list[abs(INTROGRESSED_BETA)>quantile(abs(INTROGRESSED_BETA),.25),]


# emVARs_obs_ctc[,.(.N,Pct_blood=mean(!is.na(blood_eQTL) & blood_eQTL),Pct_lung=mean(!is.na(lung_eQTL) & lung_eQTL),Pct_liver=mean(!is.na(liver_eQTL) & liver_eQTL)),keyby=.(celline,is_emVar_CRITERIA)]
# emVARs_obs_ctc[!is.na(lung_eQTL),.(.N,Pct_K562=mean(is_emVar_CRITERIA[celline=='K562']),Pct_HEPG2=mean(is_emVar_CRITERIA[celline=='HepG2']),Pct_A549=mean(is_emVar_CRITERIA[celline=='A549'])),keyby=.(blood_eQTL,liver_eQTL,lung_eQTL)]

eQTL_vs_emVar <- merge(eQTL_list, emVARs_annot_obs_ctc[, .(ID, is_emVar_CRITERIA, log2FC_archaic_vs_modern, celline, condition, COND_ID, Nearest)], by = "ID", allow.cartesian = TRUE, all.y = TRUE)
# eQTL_vs_emVar <- merge(eQTL_list, emVARs_annot_obs_ctc[, .(ID, is_emVar_CRITERIA, log2FC_archaic_vs_modern, celline, condition, COND_ID, Nearest)], by.x = c("ID",'eQTL_gene'), by.y = c("ID","Nearest"), allow.cartesian = TRUE, all.y = TRUE)


eQTL_vs_emVar[is_emVar_CRITERIA==TRUE, mean(sign(INTROGRESSED_BETA) == sign(log2FC_archaic_vs_modern),na.rm=T), by = COND_ID]
eQTL_vs_emVar[is_emVar_CRITERIA==TRUE,][order(Nearest, -abs(log2FC_archaic_vs_modern))][!duplicated(paste(Nearest, COND_ID))][, cor(INTROGRESSED_BETA, log2FC_archaic_vs_modern,use='p'), by = COND_ID]

has_eQTL <- eQTL_vs_emVar[is_emVar_CRITERIA == TRUE, .(
  has_eQTL = any(!is.na(eQTL_study)),
  has_same_direction_eQTL = any(!is.na(eQTL_study) & sign(INTROGRESSED_BETA) == sign(log2FC_archaic_vs_modern)),
  has_opposite_direction_eQTL = any(!is.na(eQTL_study) & sign(INTROGRESSED_BETA) != sign(log2FC_archaic_vs_modern))
), by = .(ID, COND_ID)]
has_eQTL_all <- eQTL_vs_emVar[is_emVar_CRITERIA == TRUE, .(
  has_eQTL = any(!is.na(eQTL_study)),
  has_same_direction_eQTL = any(!is.na(eQTL_study) & sign(INTROGRESSED_BETA) == sign(log2FC_archaic_vs_modern)),
  has_opposite_direction_eQTL = any(!is.na(eQTL_study) & sign(INTROGRESSED_BETA) != sign(log2FC_archaic_vs_modern))
), by = .(ID)]

# has_eQTL <- eQTL_vs_emVar[is_emVar_CRITERIA == TRUE & abs(log2FC_archaic_vs_modern)>.5, .(
#   has_eQTL = any(!is.na(eQTL_study)),
#   has_same_direction_eQTL = any(!is.na(eQTL_study) & sign(INTROGRESSED_BETA) == sign(log2FC_archaic_vs_modern)),
#   has_opposite_direction_eQTL = any(!is.na(eQTL_study) & sign(INTROGRESSED_BETA) != sign(log2FC_archaic_vs_modern))
# ), by = .(ID, COND_ID)]
# has_eQTL_all <- eQTL_vs_emVar[is_emVar_CRITERIA == TRUE & abs(log2FC_archaic_vs_modern)>.5, .(
#   has_eQTL = any(!is.na(eQTL_study)),
#   has_same_direction_eQTL = any(!is.na(eQTL_study) & sign(INTROGRESSED_BETA) == sign(log2FC_archaic_vs_modern)),
#   has_opposite_direction_eQTL = any(!is.na(eQTL_study) & sign(INTROGRESSED_BETA) != sign(log2FC_archaic_vs_modern))
# ), by = .(ID)]
cat("Percentage emVars with an eQTL (any eQTL catalog):\n")

count_eQTL_emVar <- has_eQTL_all[, .(
  has_same_direction_eQTL = sum(has_same_direction_eQTL),
  has_opposite_direction_eQTL = sum(has_opposite_direction_eQTL),
  has_same_direction_eQTL_only = sum(has_same_direction_eQTL & !has_opposite_direction_eQTL),
  has_opposite_direction_eQTL_only = sum(has_opposite_direction_eQTL & !has_same_direction_eQTL),
  has_same_and_opposite_direction_eQTL = sum(has_opposite_direction_eQTL & has_same_direction_eQTL),
  has_any_eQTL = sum(has_eQTL)
)]


cat("alltogether\n")
print(has_eQTL_all[, mean(has_eQTL)])
# 0.3613933

cat("detail by condition\n")
print(has_eQTL[, mean(has_eQTL), by = COND_ID])

#          COND_ID        V1
#  1:     HepG2_NS 0.3360000
#  2:     A549_IAV 0.3894737
#  3:      A549_NS 0.3387097
#  4:    A549_SARS 0.3281250
#  5:    A549_TNFa 0.3818182
#  6:    HepG2_DEX 0.3673469
#  7: HepG2_IFNA2b 0.3926941
#  8:   HepG2_TNFa 0.3589744
#  9:     K562_DEX 0.3980583
# 10:  K562_IFNA2b 0.4571429
# 11:      K562_NS 0.3750000
# 12:    K562_TNFa 0.4318182


cat("count by direction\n")
count_eQTL_emVar <- has_eQTL_all[, .(
  has_same_direction_eQTL = sum(has_same_direction_eQTL),
  has_opposite_direction_eQTL = sum(has_opposite_direction_eQTL),
  has_same_direction_eQTL_only = sum(has_same_direction_eQTL & !has_opposite_direction_eQTL),
  has_opposite_direction_eQTL_only = sum(has_opposite_direction_eQTL & !has_same_direction_eQTL),
  has_same_and_opposite_direction_eQTL = sum(has_opposite_direction_eQTL & has_same_direction_eQTL),
  has_any_eQTL = sum(has_eQTL)
)]
print(count_eQTL_emVar)

cat("Percentages\n")

print(count_eQTL_emVar[,.(has_same_direction_eQTL=has_same_direction_eQTL/has_any_eQTL,has_opposite_direction_eQTL_only=has_opposite_direction_eQTL_only/has_any_eQTL,has_same_and_opposite_direction_eQTL=has_same_and_opposite_direction_eQTL/has_any_eQTL)])




emVARs_annot_obs_ctc[is_emVar_CRITERIA==TRUE & condition=='NS',length(unique(crsID))]
# [1] 342


eQTL_vs_emVar <- merge(eQTL_list, emVARs_annot_obs_ctc[, .(ID, is_emVar_CRITERIA, log2FC_archaic_vs_modern, celline, condition, COND_ID, Nearest)], by = "ID", allow.cartesian = TRUE, all.y = TRUE)
# eQTL_vs_emVar[is_emVar_CRITERIA==TRUE, mean(sign(INTROGRESSED_BETA) == sign(log2FC_archaic_vs_modern),na.rm=T), by = COND_ID]
# eQTL_vs_emVar[is_emVar_CRITERIA==TRUE,][order(Nearest, -abs(log2FC_archaic_vs_modern))][!duplicated(paste(Nearest, COND_ID))][, cor(INTROGRESSED_BETA, log2FC_archaic_vs_modern,use='p'), by = COND_ID]

has_eQTL <- eQTL_vs_emVar[is_emVar_CRITERIA == TRUE, .(
  has_eQTL = any(!is.na(eQTL_study)),
  has_same_direction_eQTL = any(!is.na(eQTL_study) & sign(INTROGRESSED_BETA) == sign(log2FC_archaic_vs_modern)),
  has_opposite_direction_eQTL = any(!is.na(eQTL_study) & sign(INTROGRESSED_BETA) != sign(log2FC_archaic_vs_modern))
), by = .(ID, COND_ID)]

has_eQTL_all <- eQTL_vs_emVar[is_emVar_CRITERIA == TRUE, .(
  has_eQTL = any(!is.na(eQTL_study)),
  has_same_direction_eQTL = any(!is.na(eQTL_study) & sign(INTROGRESSED_BETA) == sign(log2FC_archaic_vs_modern)),
  has_opposite_direction_eQTL = any(!is.na(eQTL_study) & sign(INTROGRESSED_BETA) != sign(log2FC_archaic_vs_modern))
), by = .(ID)]

has_eQTL_NS <- eQTL_vs_emVar[is_emVar_CRITERIA == TRUE & condition=='NS', .(
  has_eQTL = any(!is.na(eQTL_study)),
  has_same_direction_eQTL = any(!is.na(eQTL_study) & sign(INTROGRESSED_BETA) == sign(log2FC_archaic_vs_modern)),
  has_opposite_direction_eQTL = any(!is.na(eQTL_study) & sign(INTROGRESSED_BETA) != sign(log2FC_archaic_vs_modern))
), by = .(ID)]

cat("Percentage emVars with an eQTL (any GTEx only):\n")
cat("alltogether\n")
has_eQTL_all[, mean(has_eQTL)] #

cat("by condition\n")
has_eQTL[, mean(has_eQTL), by = COND_ID]

cat("NS only\n")
has_eQTL_NS[, mean(has_eQTL)]
has_eQTL_NS[, mean(has_eQTL),by=COND_ID]

cat("detail by direction\n")
count_eQTL_emVar <- has_eQTL_all[, .(
  has_same_direction_eQTL = sum(has_same_direction_eQTL),
  has_opposite_direction_eQTL = sum(has_opposite_direction_eQTL),
  has_same_direction_eQTL_only = sum(has_same_direction_eQTL & !has_opposite_direction_eQTL),
  has_opposite_direction_eQTL_only = sum(has_opposite_direction_eQTL & !has_same_direction_eQTL),
  has_same_and_opposite_direction_eQTL = sum(has_opposite_direction_eQTL & has_same_direction_eQTL),
  has_any_eQTL = sum(has_eQTL)
)]
count_eQTL_emVar

cat("Percentages\n")
count_eQTL_emVar[,.(has_same_direction_eQTL=has_same_direction_eQTL/has_any_eQTL,has_opposite_direction_eQTL_only=has_opposite_direction_eQTL_only/has_any_eQTL,has_same_and_opposite_direction_eQTL=has_same_and_opposite_direction_eQTL/has_any_eQTL)]

cat("detail by direction (NS only)\n")
count_eQTL_emVar_NS <- has_eQTL_NS[, .(
  has_same_direction_eQTL = sum(has_same_direction_eQTL),
  has_opposite_direction_eQTL = sum(has_opposite_direction_eQTL),
  has_same_direction_eQTL_only = sum(has_same_direction_eQTL & !has_opposite_direction_eQTL),
  has_opposite_direction_eQTL_only = sum(has_opposite_direction_eQTL & !has_same_direction_eQTL),
  has_same_and_opposite_direction_eQTL = sum(has_opposite_direction_eQTL & has_same_direction_eQTL),
  has_any_eQTL = sum(has_eQTL)
)]
count_eQTL_emVar_NS

cat("Percentages\n")
count_eQTL_emVar_NS[,.(has_same_direction_eQTL=has_same_direction_eQTL/has_any_eQTL,has_opposite_direction_eQTL_only=has_opposite_direction_eQTL_only/has_any_eQTL,has_same_and_opposite_direction_eQTL=has_same_and_opposite_direction_eQTL/has_any_eQTL)]

# best_eQTL_list <- eQTL_list[,head(.SD,1),by=.(variantId_hg38)]
# ctr_emVARs_eQTL <- emVARs_all_results[grepl("eQTL", oligo1) & power == 0 & pval.LRT < .01 & ANALYSIS_SUBTYPE == "celltype_cond" & perm == "perm_0_0_0", ]
# # ctr_emVARs_eQTL=merge(ctr_emVARs_eQTL,SNP_annot_v4[,.(ID,variantId_hg38,crsID, rsID, REF, ALT,INTROGRESSED.allele)],by='crsID',allow.cartesian=TRUE)
# eQTL_vs_emVar <- merge(eQTL_list, ctr_emVARs_eQTL, by = "crsID", allow.cartesian = TRUE)

# ctr_emVARs_eQTL_2 <- ctr_emVARs_eQTL[, .(Pos_emVar = sum(logFC_2vs1 > 0), Neg_emVar = sum(logFC_2vs1 < 0)), keyby = .(crsID, oligo1, oligo2)]
# eQTL_list_2 <- eQTL_list[crsID %in% ctr_emVARs_eQTL_2$crsID, .(Pos_eQTL = sum(beta > 0), Neg_eQTL = sum(beta < 0)), keyby = .(crsID, REF, ALT)]

# eQTL_vs_emVar[, .(Pos_emVar = sum((logFC_2vs1 * ifelse(allele2_is_REF, -1, 1)) > 0), Neg_emVar = sum((logFC_2vs1 * ifelse(allele2_is_REF, -1, 1)) < 0), Pos_eQTL = sum(beta > 0), Neg_eQTL = sum(beta < 0)), keyby = .(crsID, REF, ALT)] <- merge(eQTL_list, ctr_emVARs_eQTL, by = "crsID", allow.cartesian = TRUE)

# Pct_concordant_GTex <- eQTL_vs_emVar[is_emVar_CRITERIA == TRUE, .(nSNP = .N, nCOncordant_effect = sum(sign(INTROGRESSED_BETA) == sign(log2FC_archaic_vs_modern)), Pct_concondant = mean(sign(INTROGRESSED_BETA) == sign(log2FC_archaic_vs_modern))), by = COND_ID]
# Pct_concordant_GTex[, P := binom.test(nCOncordant_effect, nSNP)$p.value, by = COND_ID]

# Pct_concordant_GTex

###########################
SNP_annot_emVars <- fread(sprintf("%s/SNP_annot_emVars__%s.tsv", EMVAR_DIR, CRITERIA_ACTIVE))

SNP_annot_emVars[,emVar_A549:=posID %in% emVARs_annot_obs_ctc[is_emVar_CRITERIA==TRUE & celline=='A549',unique(crsID)]]
SNP_annot_emVars[,emVar_K562:=posID %in% emVARs_annot_obs_ctc[is_emVar_CRITERIA==TRUE & celline=='K562',unique(crsID)]]
SNP_annot_emVars[,emVar_HepG2:=posID %in% emVARs_annot_obs_ctc[is_emVar_CRITERIA==TRUE & celline=='HepG2',unique(crsID)]]

SNP_annot_emVars[,fisher.test(table(emVar_A549,emVar_HepG2))[c('estimate','p.value')]]
#    estimate      p.value
# 1: 7.544687 3.265507e-35
SNP_annot_emVars[,fisher.test(table(emVar_K562,emVar_HepG2))[c('estimate','p.value')]]
#    estimate      p.value
# 1:  4.89408 1.518064e-23
SNP_annot_emVars[,fisher.test(table(emVar_A549,emVar_K562))[c('estimate','p.value')]]
#    estimate      p.value
# 1: 10.10471 3.412669e-34

SNP_annot_emVars[nCelltype>=1,.N,keyby=nCelltype][,Pct:=N/sum(N)][seq_len(.N),]
#    nCelltype     N        Pct
# 1:         1   537 0.77939042
# 2:         2   107 0.15529753
# 3:         3    45 0.06531205

SNP_annot_emVars[nCelltype==1,.N,keyby=nCelltype_P05>1][,Pct:=N/sum(N)][seq_len(.N),]
#    nCelltype_P05     N       Pct
#           <lgcl> <int>     <num>
# 1:         FALSE   191 0.3556797
# 2:          TRUE   346 0.6443203

any_emVars <- SNP_annot_emVars[nCelltype>0,posID]
celline_specific_emVars <- SNP_annot_emVars[nCelltype==1,posID]

SNP_annot_emVars <- merge(SNP_annot_emVars, SNP_annot[, .(ID, crsID, gwas, LD_block_AFR, introgression_locus)], by = "ID")
# SNP_annot_gwas <- SNP_annot_emVars[, .(nSNP = .N, 
#                                       has_gwas = any(gwas != ""), 
#                                       nTraits = length(unique(na.omit(gsub("(.*) - ", "", unlist(strsplit(gwas, "/")))))),
#                                       traits = paste(unlist(unique(na.omit(gsub("(.*) - ", "", unlist(strsplit(gwas, "/")))))), collapse = "/"),
#                                       nEmVar = sum(nCelltype_NS != 0), 
#                                       nShared_EmVars = sum((nCelltype_NS > 0) & (nCelltype_P05_NS > 1))
#                                       ), by = introgression_locus]
# SNP_annot_gwas[, .(introgression_locus, trait = unique(gsub("(.*) - ", "", unlist(strsplit(traits, "/")))))]


SNP_annot_emVars[nCelltype > 0,mean(gwas != "")] #0.3323657
SNP_annot_emVars[nCelltype > 0 & nCelltype_P05 == 1, mean(gwas != "")]
#0.3769634
SNP_annot_emVars[nCelltype > 0 & nCelltype_P05 > 1, mean(gwas != "")]
#0.315261
# SNP_annot_emVars[nCelltype > 0,.(.N, mean(gwas != "")),keyby=nContext>1]


SNP_annot_emVars[nCelltype > 0,mean(gwas != "")] #0.3323657
SNP_annot_emVars[nCelltype > 0 & nCelltype_P05_NS == 0,mean(gwas != "")] #0.3055556
SNP_annot_emVars[nCelltype > 0 & nCelltype_P05_NS > 0,mean(gwas != "")] #0.3418468



eQTL_catalog_Full <- CS_results_finemap_eQTL_tested[!is.na(studyType),]
SNP_annot_emVars[, has_eQTL := crsID %in% eQTL_catalog_Full$crsID]

SNP_annot_emVars[gwas != "" & nCelltype > 0,mean(!has_eQTL)]
#[1] 0.1222707
SNP_annot_emVars[gwas != "" & nCelltype > 0,sum(!has_eQTL)]
# 28
SNP_annot_emVars[gwas != "" & nCelltype > 0,sum(!has_eQTL & nCelltype_NS > 0)] # 13
SNP_annot_emVars[gwas != "" & nCelltype > 0,sum(!has_eQTL & nCelltype_NS==0 & nCelltype_P05_NS > 0)] # 10
SNP_annot_emVars[gwas != "" & nCelltype > 0,sum(!has_eQTL & nCelltype_NS==0 & nCelltype_P05_NS == 0)] # 5


SNP_annot_emVars[gwas != "" & nCelltype_NS > 0,mean(!has_eQTL)]
SNP_annot_emVars[gwas != "" & nCelltype_NS == 0 & nCelltype> 0,mean(!has_eQTL)]
SNP_annot_emVars[gwas != "" & nCelltype_P05_NS == 0 & nCelltype> 0,mean(!has_eQTL)]
SNP_annot_emVars[gwas != "" & nCelltype_P05_NS > 0 & nCelltype> 0,mean(!has_eQTL)]

cat("All done !\n")
q("no")