# running: sbatch -p geh,common --mem=30G  00_Rscript.sh MPRA_count_exp6_analysisZ/03b_aggMPRAnalyse_oligo_activity.R

MPRA_DIR <- "/pasteur/helix/projects/evo_immuno_pop/MPRA"
ANALYSIS_DIR <- "MPRA_count_exp6_analysis3"
SCRIPT_DIR <- sprintf("%s/scripts/%s/", MPRA_DIR, ANALYSIS_DIR)


source(sprintf("%s/scripts/%s/00_starter.R", MPRA_DIR, ANALYSIS_DIR))
source(sprintf("%s/scripts/%s/00_load_associations.R", MPRA_DIR, ANALYSIS_DIR))
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
selected_annot_wide[,Introgression_source_top:=Introgression_source_top_initial]
selected_annot_wide[Introgression_source_top=='both',Introgression_source_top:='Vindija/Denisova']
selected_annot_wide[Introgression_source_top=='',Introgression_source_top:='Undetermined']

SNP_annot_v5 <- merge(SNP_annot_v4[,-"Introgression_scenario"],selected_annot_wide,by=c('ID','posID'))
toc()

# load oligo annotations
source(sprintf("%s/scripts/%s/02z__define_includedCRS.R", MPRA_DIR, ANALYSIS_DIR))


##### active parameters
source(sprintf("%s/scripts/%s/03_00_parameter_activity.R", MPRA_DIR, ANALYSIS_DIR))
##### diff parameters
source(sprintf("%s/scripts/%s/04_00_parameter_diff_activity.R", MPRA_DIR, ANALYSIS_DIR))
##### emVars parameters
source(sprintf("%s/scripts/%s/05_00_parameter_emVars.R", MPRA_DIR, ANALYSIS_DIR))
##### emVars diff parameters
source(sprintf("%s/scripts/%s/06_00_parameter_emVars_Diff.R", MPRA_DIR, ANALYSIS_DIR))

IN_DIR <- sprintf("%s/data/%s/03a_aggMPRA_analyse/%s/", MPRA_DIR, ANALYSIS_DIR, RUN_ID)
ACTIVE_DIR <- sprintf("%s/data/%s/03b_computeFDR/%s/Activity/%s/", MPRA_DIR, ANALYSIS_DIR, RUN_ID, CRITERIA_ACTIVE)
DIFF_DIR <- sprintf("%s/data/%s/03b_computeFDR/%s/Diff/%s/%s/", MPRA_DIR, ANALYSIS_DIR, RUN_ID, CRITERIA_DIFF, CRITERIA_ACTIVE_OUT)
EMVAR_DIR <- sprintf("%s/data/%s/03b_computeFDR/%s/emVars/%s/%s", MPRA_DIR, ANALYSIS_DIR, RUN_ID, CRITERIA_EMVARS, CRITERIA_ACTIVE_OUT)
EMVARDIFF_DIR <- sprintf("%s/data/%s/03b_computeFDR/%s/EmVar_Diff/%s/%s/%s/%s/%s", MPRA_DIR, ANALYSIS_DIR, RUN_ID, CRITERIA_EMVARS_DIFF, CRITERIA_EMVARS_OUT, CRITERIA_DIFF_OUT, CRITERIA_ACTIVE_OUT, STIM_ONLY_OUT)

dir.create(EMVARDIFF_DIR, recursive = TRUE)

FIGURE_DIR <- sprintf("%s/figures/%s/%s/05_emVarDiff/%s/%s/%s/%s/%s", MPRA_DIR, ANALYSIS_DIR, RUN_ID, CRITERIA_EMVARS_DIFF, CRITERIA_EMVARS_OUT, CRITERIA_DIFF_OUT, CRITERIA_ACTIVE_OUT, STIM_ONLY_OUT)
dir.create(FIGURE_DIR, recursive = TRUE)
dir.create(sprintf('%s/SupTables/', FIGURE_DIR, recursive = TRUE))

emVARs_all_results <- fread(sprintf("%s/all_emVars_results.tsv.gz", IN_DIR), sep = "\t")
emVARs_all_results <- emVARs_all_results[crsID %in% tested_and_ctrl_oligos_final[type=='tested',posID], ]

emVARs_Diff_all_results <- fread(sprintf("%s/all_emVars_Diff_results.tsv.gz", IN_DIR), sep = "\t")
emVARs_Diff_all_results <- emVARs_Diff_all_results[crsID %in% tested_and_ctrl_oligos_final[type=='tested',posID], ]
# oligo_activity_file <- sprintf("%s/oligo_activity__all.tsv.gz", OUT_DIR)
# oligo_activity <- fread(oligo_activity_file)

# emVARs_annot_active_ok <- fread(sprintf("%s/all_emVars_activeCRS_ok_annotated_celltype_v2.tsv", OUT_DIR))


# emVARs_Diff_annot <- merge(emVARs_Diff_all_results[df.int > 0, ], SNP_annot_v5, by = c("crsID"), allow.cartesian = TRUE)
# emVARs_Diff_annot[, celltype := str_split(ANALYSIS_NAME, "_|-", simplify = T)[, 1]]
# emVARs_Diff_annot <- merge(emVARs_Diff_annot, crs_targets, by = c("crsID", "celltype"), all.x = TRUE)

oligo_activity_obs <- fread(sprintf("%s/all_oligos_annotated__%s.tsv.gz", ACTIVE_DIR, CRITERIA_ACTIVE))

all_emVars_obs <- fread(sprintf("%s/all_emVars_annotated__%s.tsv.gz", EMVAR_DIR, CRITERIA_ACTIVE_OUT), sep = "\t")
all_Diff_obs <- fread(sprintf("%s/all_oligos_diff_annotated__%s.tsv.gz", DIFF_DIR, CRITERIA_DIFF), sep = "\t")

comparison_labels <- fread(sprintf("%s/04a_comparison_labels.txt", SCRIPT_DIR))

oligo_target <- fread(sprintf("%s/data/%s/00_oligo_targets_v4.txt", MPRA_DIR, ANALYSIS_DIR), sep = "\t")
oligo_target_collapsed <- fread(sprintf("%s/data/%s/00_oligo_targets_collapsed_v4.txt", MPRA_DIR, ANALYSIS_DIR), sep = "\t")
crs_targets <- unique(oligo_target_collapsed[, -"oligo"])
setnames(crs_targets, "celltype", "celline")

# emVARs_annot <- merge(emVARs_annot, unique(oligo_activity_obs[, .(cre_class_CRITERIA, posID = crsID, ANALYSIS_NAME)])
#   , by = c("posID", "ANALYSIS_NAME"), all.x = TRUE, suffix = c(".emVAR", ".oligo"))
# # emVARs_annot <- emVARs_annot[is_tested_cre == TRUE, ]

# if (any(grepl(".x$", colnames(emVARs_annot)))) {
#   stop("merge issue")
# }



####################################################################
################# compute FDR,  response emVars ####################
####################################################################
cat('\nchecking number of tested SNP per condition/permutation\n')
emVARs_Diff_all_results[, .(.N, sum(df.test_int > 0)), keyby = .(power, boot, ANALYSIS_NAME, Z_th, perm)]

cat('\ncomputing FDR\n')
emVARs_Diff_annot <- merge(emVARs_Diff_all_results[df.test_int > 0, ], SNP_annot_v5, by = "crsID", allow.cartesian = TRUE)
emVARs_Diff_annot[, Delta_log2FC_der_vs_anc := ifelse(allele.1 == ANCESTRAL, delta_logFC_2vs1 / log(2), -delta_logFC_2vs1 / log(2))]
emVARs_Diff_annot[, Delta_log2FC_archaic_vs_modern := ifelse(allele.2 == INTROGRESSED.allele, delta_logFC_2vs1 / log(2), -delta_logFC_2vs1 / log(2))]
emVARs_Diff_annot[, delta_log2FC.se := delta_logFC.se / log(2)]

emVARs_Diff_annot <- merge(emVARs_Diff_annot, comparison_labels, by.x = "ANALYSIS_NAME", by.y = "analysis_name")

########## ########## ########## ########## ########## ########## ########## ########## 
########## manage filtering criteria to handle filtering on emvars and/or diff ######## 
########## ########## ########## ########## ########## ########## ########## ########## 

# subset to cres that are active/emVars/Diff in the current condition/celltype
  all_active_cres <- oligo_activity_obs[oligo_active_CRITERIA == TRUE, .(oligo, crsID, ANALYSIS_NAME, ANALYSIS_SUBTYPE)]
	all_diff_cres <- all_Diff_obs[oligo_diff_CRITERIA == TRUE, .(crsID, ANALYSIS_NAME, ANALYSIS_SUBTYPE)]
	all_emvar_cres <- all_emVars_obs[is_emVar_CRITERIA == TRUE, .(crsID, ANALYSIS_NAME, ANALYSIS_SUBTYPE)]

# active filter
 if (FILTER_ACTIVE==TRUE) {
	 emVARs_Diff_annot[, is_active_cre_group1 := paste(crsID, group1_labels) %chin% all_active_cres[, paste(crsID, ANALYSIS_NAME)]]
   emVARs_Diff_annot[, is_active_cre_group2 := paste(crsID, group2_labels) %chin% all_active_cres[, paste(crsID, ANALYSIS_NAME)]]
   emVARs_Diff_annot[, is_active_cre := is_active_cre_group1 | is_active_cre_group2 ]
 } else {
   emVARs_Diff_annot[, is_active_cre := TRUE]
 }

# diff filter
 if (FILTER_DIFF==TRUE) {
	 emVARs_Diff_annot[, is_diff_cre := paste(crsID, ANALYSIS_NAME) %chin% all_diff_cres[, paste(crsID, ANALYSIS_NAME)]]
	} else {
   emVARs_Diff_annot[, is_diff_cre := TRUE]
 }

# emVar filter
  if (FILTER_EMVAR==TRUE) {
	 	emVARs_Diff_annot[, is_emvar_cre_group1 := paste(crsID, group1_labels) %chin% all_emvar_cres[, paste(crsID, ANALYSIS_NAME)]]
    emVARs_Diff_annot[, is_emvar_cre_group2 := paste(crsID, group2_labels) %chin% all_emvar_cres[, paste(crsID, ANALYSIS_NAME)]]
    if(STIM_ONLY==TRUE) {
      emVARs_Diff_annot[, is_emvar_cre := is_emvar_cre_group2 ]
    }else{
      emVARs_Diff_annot[, is_emvar_cre := is_emvar_cre_group1 | is_emvar_cre_group2 ]
    }
  } else {
   emVARs_Diff_annot[, is_emvar_cre := TRUE]
 }
# combine filters
emVARs_Diff_annot[, is_tested_cre := is_active_cre & is_emvar_cre & is_diff_cre ]

########## ########## ########## ########## ########## ########## ########## ########## 
########## compute local and global FDR on choson CREs/conditions/celltypes ########### 
########## ########## ########## ########## ########## ########## ########## ########## 
# emVARs_Diff_annot[FDR < .05 & power == 0 & perm == "perm_0_0_0" & grepl("celltype", ANALYSIS_SUBTYPE), .N, by = ANALYSIS_NAME]
emVARs_Diff_annot[, bin_logP := cut(-log10(pval_int), breaks = c(seq(0, 5, by = 0.25), Inf))]

lfdr_estim <- emVARs_Diff_annot[, .(n0 = sum(perm != "perm_0_0_0" & df.test_int > 0 & is_tested_cre), n1 = sum(perm == "perm_0_0_0" & df.test_int > 0  & is_tested_cre)), keyby = .(Z_th, power, boot, ANALYSIS_NAME, bin_logP)]
lfdr_estim[, lfdr := pmin(1, n0 / n1), keyby = .(Z_th, power, boot, ANALYSIS_NAME, bin_logP)]
lfdr_estim[, monotonic_bins := cumsum(!duplicated(rev(cummax(rev(lfdr))))), keyby = .(Z_th, power, boot, ANALYSIS_NAME)]
lfdr_estim[, lfdr := pmin(1, sum(n0) / sum(n1)), keyby = .(Z_th, power, boot, ANALYSIS_NAME, monotonic_bins)]
emVARs_Diff_annot <- merge(emVARs_Diff_annot, lfdr_estim[, .(Z_th, power, boot, ANALYSIS_NAME, bin_logP, lfdr)], by = c("Z_th", "power", "boot", "ANALYSIS_NAME", "bin_logP"), all.x = TRUE)
emVARs_Diff_annot[ is_tested_cre == FALSE, lfdr := 1]

emVARs_Diff_annot <- emVARs_Diff_annot[order(perm, pval_int), .SD, keyby = .(Z_th, power, boot, ANALYSIS_NAME)]

# emVARs_Diff_annot <- merge(all_emVars_annot[FDR <= 0.05, .(Z_th, power, boot, ANALYSIS_NAME, crsID)], emVARs_Diff_annot, by = c("Z_th", "power", "boot", "ANALYSIS_NAME", "crsID"), all = TRUE)
emVARs_Diff_annot <- emVARs_Diff_annot[order(pval_int), .SD, by = .(Z_th, power, boot, ANALYSIS_NAME)]
emVARs_Diff_annot[, N_FP := (.1 + cumsum(perm != "perm_0_0_0" & df.test_int > 0 & is_tested_cre)), by = .(Z_th, power, boot, ANALYSIS_NAME)]
emVARs_Diff_annot[, N_POS := (.1 + cumsum(perm == "perm_0_0_0" & df.test_int > 0 & is_tested_cre)), by = .(Z_th, power, boot, ANALYSIS_NAME)]
emVARs_Diff_annot[, pval_emp := N_FP/(.1 + sum(perm == "perm_0_0_0" & df.test_int > 0 & is_tested_cre)), by = .(Z_th, power, boot, ANALYSIS_NAME)]
emVARs_Diff_annot[ is_tested_cre == TRUE , FDR := pmin(1, rev(cummin(rev(N_FP / N_POS)))), by = .(Z_th, power, boot, ANALYSIS_NAME)]
emVARs_Diff_annot[ is_tested_cre == FALSE , FDR := 1]

emVARs_Diff_annot[, is_emVar_Diff_CRITERIA := FDR < FDR_TH_EMVARDIFF & lfdr < LFDR_TH_EMVARDIFF & abs(Delta_log2FC_archaic_vs_modern) > LOG2FC_TH_EMVARDIFF ]

emVARs_Diff_obs <- emVARs_Diff_annot[power == 0 & perm == "perm_0_0_0", ]
emVARs_Diff_perm <- emVARs_Diff_annot[power == 0 & perm == "perm_0_1_0", ]
emVARs_Diff_power <- emVARs_Diff_annot[power == 1, ]

fwrite(emVARs_Diff_obs, file = sprintf("%s/all_emVars_diff_annotated.tsv.gz", EMVARDIFF_DIR), sep = "\t")
fwrite(emVARs_Diff_perm, file = sprintf("%s/all_emVars_diff_annotated_perm.tsv.gz", EMVARDIFF_DIR), sep = "\t")
fwrite(emVARs_Diff_power, file = sprintf("%s/all_emVars_diff_power.tsv.gz", EMVARDIFF_DIR), sep = "\t")


####################################################################
################# extract final sets of emVars  ####################
####################################################################

emVARs_Diff_obs_response <- emVARs_Diff_obs[grepl("response", ANALYSIS_SUBTYPE),]
COND_SUMMARY <- condition_summary[, .(celline, condition, COND_ID, ANALYSIS_NAME=gsub('(.*)_all','response_\\1',analysis_name))]
emVARs_Diff_obs_response <- merge(emVARs_Diff_obs_response, COND_SUMMARY, by= "ANALYSIS_NAME")

fwrite(emVARs_Diff_obs_response, file = sprintf("%s/all_emVars_diff_obs_response.tsv", EMVARDIFF_DIR), sep = "\t")

N_ResponseEmvar <- emVARs_Diff_obs_response[is_emVar_Diff_CRITERIA==TRUE,length(unique(crsID))]
N_ResponseEmvar_x_cond <- emVARs_Diff_obs_response[is_emVar_Diff_CRITERIA==TRUE,length(unique(paste(crsID, celline, condition)))]

cat(sprintf('%s response emVars identified, corresponding to a total of %s significant interactions', N_ResponseEmvar, N_ResponseEmvar_x_cond))

emVARs_Diff_obs_celltype <- emVARs_Diff_obs[grepl("celltype_comp", ANALYSIS_SUBTYPE),]
fwrite(emVARs_Diff_obs_celltype, file = sprintf("%s/all_emVars_diff_obs_celltype.tsv", EMVARDIFF_DIR), sep = "\t")
emVARs_Diff_obs_celltype <- emVARs_Diff_obs_celltype[is_emVar_Diff_CRITERIA==TRUE & ANALYSIS_SUBTYPE=='celltype_comp',]
emVARs_Diff_obs_celltype <- emVARs_Diff_obs_celltype[posID %in% tested_and_ctrl_oligos_final_annot[type=='tested',posID],]
emVARs_Diff_obs_celltype[is_emVar_Diff_CRITERIA==TRUE,length(unique(posID))] # 282
emVARs_Diff_obs_celltype[, celline1_labels :=  gsub("-ACE2", "", gsub("Celltype_(.*)_(.*)_(.*)", "\\1", ANALYSIS_NAME))]
emVARs_Diff_obs_celltype[, celline2_labels := gsub("-ACE2", "", gsub("Celltype_(.*)_(.*)_(.*)", "\\2", ANALYSIS_NAME))]
emVARs_Diff_obs_celltype[, condition :=  gsub("-ACE2", "", gsub("Celltype_(.*)_(.*)_(.*)", "\\3", ANALYSIS_NAME))]

setnames(crs_targets,'celltype','celline',skip_absent=TRUE)
emVARs_Diff_obs_celltype <- merge(emVARs_Diff_obs_celltype, unique(crs_targets[,.(crsID,target_any_proximity,target_any_contact)]), by = c("crsID"), all.x = TRUE)

SupTable3d <- emVARs_Diff_obs_celltype[,.(celline_1=celline1_labels,celline_2=celline2_labels,condition,crsID, ID,variantId_hg38,rsID,CHROM,POS_b37,
INTROGRESSED=INTROGRESSED.allele, `NON-INTROGRESSED`=ifelse(INTROGRESSED.allele==REF,ALT,REF),
introgression_locus, 
nObs_allele1_Cell_line_1 = nBCs_g1a1, nObs_allele2_Cell_line_1 = nBCs_g1a2, nObs_allele1_Cell_line_2 = nBCs_g2a1, nObs_allele2_Cell_line_2 = nBCs_g2a2,
activity_allele1_Cell_line_1=group1_a1, activity_allele2Cell_line_1=group1_a2, activity_allele1_Cell_line_2=group2_a1, activity_allele2_Cell_line_2=group2_a2,
is_emvar_Cell_line_1=is_emvar_cre_group1,is_emvar_Cell_line_2=is_emvar_cre_group2,
Delta_log2FC_archaic_vs_modern, delta_log2FC.se, pval_int, pval_emp, lfdr, FDR, is_emVar_Diff_CRITERIA,
eQTL_in=effect_in,
gwas_trait=gwas,
POP_introgressed, Source_introgression, allele_match, POP_adaptive_inital_def, POP_adaptive_top, Introgression_scenario_top, Introgression_source_top, max_intro_allele_freq, max_intro_haplo_freq, 
target_any_proximity,target_any_contact, Nearest, Lung, Mono, Tcell)]

fwrite(SupTable3d , file=sprintf("%s/SupTables/SupTable3d_celllineDiff_emVars.tsv", FIGURE_DIR),sep='\t')


N_CelltypeDiff_Emvar <- emVARs_Diff_obs_celltype[ is_emVar_Diff_CRITERIA==TRUE ,length(unique(crsID))]
N_CelltypeDiff_Emvar_x_comp <- emVARs_Diff_obs_celltype[ is_emVar_Diff_CRITERIA==TRUE , length(unique(paste(crsID, ANALYSIS_NAME)))]

cat(sprintf('%s celltype-dependent emVars identified, corresponding to a total of %s significant interactions', N_CelltypeDiff_Emvar, N_CelltypeDiff_Emvar_x_comp))

response_emVars_obs <- fread(sprintf("%s/all_emVars_diff_obs_response.tsv", EMVARDIFF_DIR))
response_emVars_obs <- response_emVars_obs[posID %in% tested_and_ctrl_oligos_final_annot[type=='tested',posID],]
response_emVars_obs[is_emVar_Diff_CRITERIA==TRUE,length(unique(posID))] # 94

response_emVars_obs <- merge(response_emVars_obs, crs_targets, by = c("crsID", "celline"), all.x = TRUE)

CS_results_finemap_gwas_tested <- fread(sprintf("%s/data/%s/CredibleSets/CredibleSets_gwas.tsv.gz", MPRA_DIR, ANALYSIS_DIR))
trait_annot=fread(sprintf("%s/data/%s/CredibleSets/gwas_trait_annotation.txt", MPRA_DIR, ANALYSIS_DIR))
trait_annot[,studyId:=gsub('(.*) - (.*)','\\1',V1)]
trait_annot[studyId=="FINNGEN_R12_E4_HYTHY_AI_STRICT",is_immune:=0]
CS_results_finemap_gwas_tested <- merge(CS_results_finemap_gwas_tested,trait_annot,by='studyId',all.x=TRUE)
response_emVars_IDs <- response_emVars_obs[is_emVar_Diff_CRITERIA==TRUE,unique(variantId_hg38)]
CS_results_finemap_gwas_tested[variantId%in%response_emVars_IDs & is_immune==1,length(unique(variantId))] 
#17
CS_results_finemap_gwas_tested[,full_name:=gsub('(.*) - NA', '\\1',full_name)]

CS_results_finemap_gwas_tested[variantId%in%response_emVars_IDs & is_immune==1,.(variantId,trait=full_name,finemappingMethod)]

SupTable4f <- merge(CS_results_finemap_gwas_tested[variantId%in%response_emVars_IDs & is_immune==1,.(variantId_hg38=variantId,trait=full_name,finemappingMethod)], 
      response_emVars_obs[is_emVar_Diff_CRITERIA==TRUE, .(celline,condition,crsID, variantId_hg38,rsID)], by='variantId_hg38')

fwrite(SupTable4f, file=sprintf("%s/SupTables/SupTable4f_response_emVars_gwas_immune.tsv", FIGURE_DIR),sep='\t')

SupTable4b <- response_emVars_obs[,.(celline,condition,crsID, ID,variantId_hg38,rsID,CHROM,POS_b37,
INTROGRESSED=INTROGRESSED.allele, `NON-INTROGRESSED`=ifelse(INTROGRESSED.allele==REF,ALT,REF),
introgression_locus, 
nObs_allele1_ns = nBCs_g1a1, nObs_allele2_ns = nBCs_g1a2, nObs_allele1_stim = nBCs_g2a1, nObs_allele2_stim = nBCs_g2a2,
activity_allele1_NS=group1_a1, activity_allele2_NS=group1_a2, activity_allele1_STIM=group2_a1, activity_allele2_STIM=group2_a2,
is_emvar_NS=is_emvar_cre_group1,is_emvar_STIM=is_emvar_cre_group2,
Delta_log2FC_archaic_vs_modern, delta_log2FC.se, pval_int, pval_emp, lfdr, FDR, is_emVar_Diff_CRITERIA,
eQTL_in=effect_in,
gwas_trait=gwas,
POP_introgressed, Source_introgression, allele_match, POP_adaptive_inital_def, POP_adaptive_top, Introgression_scenario_top, Introgression_source_top, max_intro_allele_freq, max_intro_haplo_freq, 
target_any_proximity,target_any_contact,target_celltype_proximity,target_celltype_contact, Nearest, Lung, Mono, Tcell)]

fwrite(SupTable4b , file=sprintf("%s/SupTables/SupTable4b_response_emVars.tsv", FIGURE_DIR),sep='\t')

cat("\nAll done\n")
q("no")

# ####################################################################
# ################# compute FDR,  response emVars ####################
# #################### significant emVars only #######################
# ####################################################################

# # extract significant emVars at 5% FDR
# emVars_signif <- all_emVars_annot[power == 0 & perm == "perm_0_0_0" & grepl("celltype_cond", ANALYSIS_SUBTYPE) & FDR < 0.05, .(crsID, ANALYSIS_NAME, oligo1, oligo2, GC, pval.LRT, log2FC_archaic_vs_modern, log2FC.se,
#   FDR_emVAR = FDR, oligo_class = oligo_class_loose, selection_criteria, selection_criteria.simple, Lung, Mono, Nearest, Tcell, CHROM, POS_b37, allele.1, allele.2, allele2_is_REF, rsID,
#   Introgressed_in, Introgressed_from, Introgression_scenario, POP_adaptive, Adaptive_from,
#   ANCESTRAL, DERIVED, INTROGRESSED.allele, REF, ALT, Vindija.der, Chagyrskaya.der, Altai.der, Denisova.der, MaxPosterior_D_VY, MaxPosterior_V_DY,
#   YRI.der, SGDP_African.der, GBR.der, IBS.der, SGDP_WestEurasian.der, CHB.der, JPT.der, SGDP_Eastasian.der, PJL.der, STU.der, SGDP_Agta.der, SGDP_Papuan.der, ctrl, excluded
# )]
# emVars_signif <- merge(emVars_signif, condition_summary, by.x = "ANALYSIS_NAME", by.y = "analysis_name")

# # compute FDR for response emVars, focusing on significant emVars
# response_emVars <- emVARs_Diff_annot[power == 0 & grepl("response", ANALYSIS_SUBTYPE), ]
# response_emVars[, analysis_name := gsub("response_(.*)", "\\1_all", ANALYSIS_NAME)]
# response_emVars <- merge(response_emVars, condition_summary, by = "analysis_name")
# response_emVars <- merge(emVars_signif, response_emVars, by = c("crsID", "COND_ID"), suffix = c(".emVar", ".response"))
# response_emVars <- response_emVars[paste(posID, COND_ID) %chin% response_emVars[perm == "perm_0_0_0", paste(posID, COND_ID)], ]
# response_emVars <- response_emVars[order(pval_int), .SD, by = ANALYSIS_NAME.response]
# response_emVars[, N_FP := (.1 + cumsum(perm != "perm_0_0_0" & df.test_int > 0)), by = .(ANALYSIS_NAME.response)]
# response_emVars[, N_POS := (.1 + cumsum(perm == "perm_0_0_0" & df.test_int > 0)), by = .(ANALYSIS_NAME.response)]
# response_emVars[, FDR := pmin(1, rev(cummin(rev(N_FP / N_POS)))), by = .(ANALYSIS_NAME.response)]

# lfdr_estim <- response_emVars[, .(n0 = sum(perm != "perm_0_0_0" & df.test_int > 0), n1 = sum(perm == "perm_0_0_0" & df.test_int > 0)), keyby = .(Z_th, power, boot, ANALYSIS_NAME.response, bin_logP)]
# lfdr_estim[, lfdr := pmin(1, n0 / n1), keyby = .(Z_th, power, boot, ANALYSIS_NAME.response, bin_logP)]
# lfdr_estim[, monotonic_bins := cumsum(!duplicated(rev(cummax(rev(lfdr))))), keyby = .(Z_th, power, boot, ANALYSIS_NAME.response)]
# lfdr_estim[, lfdr := pmin(1, sum(n0) / sum(n1)), keyby = .(Z_th, power, boot, ANALYSIS_NAME.response, monotonic_bins)]
# response_emVars <- merge(response_emVars, lfdr_estim[, .(Z_th, power, boot, ANALYSIS_NAME.response, bin_logP, lfdr)], by = c("Z_th", "power", "boot", "ANALYSIS_NAME.response", "bin_logP"), all.x = TRUE)
# response_emVars <- response_emVars[order(perm, pval_int), .SD, keyby = .(Z_th, power, boot, ANALYSIS_NAME.response)]
# fwrite(response_emVars[power == 0 & perm == "perm_0_0_0", ], file = sprintf("%s/all_response_emVars.tsv", OUT_DIR), sep = "\t")

# # emVars_signif[!grepl("weak|inactive", oligo_class) & type.oligo %in% c("archaic", "COVID"), length(unique(posID))]
# ####################################################################
# ################# compute FDR,  response emVars ####################
# ############ active CRS with significant emVars only ###############
# ####################################################################

# # extract significant emVars at 5% FDR
# emVars_signif_active <- emVARs_annot_active_ok[power == 0 & perm == "perm_0_0_0" & grepl("celltype_cond", ANALYSIS_SUBTYPE) & FDR < 0.05, .(crsID, ANALYSIS_NAME, oligo1, oligo2, GC, pval.LRT, log2FC_archaic_vs_modern, log2FC.se,
#   FDR_emVAR = FDR, oligo_class_loose = oligo_class_loose
#   # ,selection_criteria, selection_criteria.simple, Lung, Mono, Nearest, Tcell, CHROM, POS_b37, allele.1, allele.2, allele2_is_REF, rsID,
#   # Introgressed_in, Introgressed_from, Introgression_scenario, POP_adaptive, Adaptive_from,
#   # ANCESTRAL, DERIVED, INTROGRESSED.allele, REF, ALT, Vindija.der, Chagyrskaya.der, Altai.der, Denisova.der, MaxPosterior_D_VY, MaxPosterior_V_DY,
#   # YRI.der, SGDP_African.der, GBR.der, IBS.der, SGDP_WestEurasian.der, CHB.der, JPT.der, SGDP_Eastasian.der, PJL.der, STU.der, SGDP_Agta.der, SGDP_Papuan.der, ctrl, excluded
# )]
# emVars_signif_active <- merge(emVars_signif_active, condition_summary, by.x = "ANALYSIS_NAME", by.y = "analysis_name")

# # compute FDR for response emVars, focusing on significant emVars
# response_emVars_active <- emVARs_Diff_annot[power == 0 & grepl("response", ANALYSIS_SUBTYPE), ]
# response_emVars_active[, analysis_name := gsub("response_(.*)", "\\1_all", ANALYSIS_NAME)]
# response_emVars_active <- merge(response_emVars_active, condition_summary, by = "analysis_name")

# response_emVars_active <- merge(emVars_signif_active, response_emVars_active, by = c("crsID", "COND_ID"), suffix = c(".emVar", ".response"))
# response_emVars_active <- response_emVars_active[paste(posID, COND_ID) %chin% response_emVars[perm == "perm_0_0_0", paste(posID, COND_ID)], ]
# response_emVars_active <- response_emVars_active[order(pval_int), .SD, by = ANALYSIS_NAME.response]
# response_emVars_active[, N_FP := (.1 + cumsum(perm != "perm_0_0_0" & df.test_int > 0)), by = .(ANALYSIS_NAME.response)]
# response_emVars_active[, N_POS := (.1 + cumsum(perm == "perm_0_0_0" & df.test_int > 0)), by = .(ANALYSIS_NAME.response)]
# response_emVars_active[, FDR := pmin(1, rev(cummin(rev(N_FP / N_POS)))), by = .(ANALYSIS_NAME.response)]

# lfdr_estim <- response_emVars_active[, .(n0 = sum(perm != "perm_0_0_0" & df.test_int > 0), n1 = sum(perm == "perm_0_0_0" & df.test_int > 0)), keyby = .(Z_th, power, boot, ANALYSIS_NAME.response, bin_logP)]
# lfdr_estim[, lfdr := pmin(1, n0 / n1), keyby = .(Z_th, power, boot, ANALYSIS_NAME.response, bin_logP)]
# lfdr_estim[, monotonic_bins := cumsum(!duplicated(rev(cummax(rev(lfdr))))), keyby = .(Z_th, power, boot, ANALYSIS_NAME.response)]
# lfdr_estim[, lfdr := pmin(1, sum(n0) / sum(n1)), keyby = .(Z_th, power, boot, ANALYSIS_NAME.response, monotonic_bins)]
# response_emVars_active <- merge(response_emVars_active, lfdr_estim[, .(Z_th, power, boot, ANALYSIS_NAME.response, bin_logP, lfdr)], by = c("Z_th", "power", "boot", "ANALYSIS_NAME.response", "bin_logP"), all.x = TRUE)
# response_emVars_active <- response_emVars_active[order(perm, pval_int), .SD, keyby = .(Z_th, power, boot, ANALYSIS_NAME.response)]
# fwrite(response_emVars_active[power == 0 & perm == "perm_0_0_0", ], file = sprintf("%s/all_response_emVars_active.tsv", OUT_DIR), sep = "\t")

# response_emVars_active_signif <- response_emVars_active[FDR < .05 & selection_criteria.simple %in% c("top5pct_archaic_ImmVIP", "COVID_LD") & power == 0 & perm == "perm_0_0_0", ][order(crsID), .(crsID, rsID, COND_ID, oligo_class_loose, selection_criteria, pval.LRT, log2FC_archaic_vs_modern, log2FC.se, FDR_emVAR, Nearest, Lung, Mono, Tcell, Introgressed_in, Introgressed_from, Introgression_scenario, POP_adaptive, Adaptive_from, ANCESTRAL, DERIVED, INTROGRESSED.allele, Delta_log2FC_archaic_vs_modern, delta_log2FC.se, pval_int, FDR)]
# fwrite(response_emVars_active_signif, file = sprintf("%s/all_response_emVars_active_signif.tsv", OUT_DIR), sep = "\t")

# cat("\nAll done\n")
# q("no")
