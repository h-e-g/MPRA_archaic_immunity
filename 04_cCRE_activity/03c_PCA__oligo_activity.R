
MPRA_DIR <- "/pasteur/helix/projects/evo_immuno_pop/MPRA"
ANALYSIS_DIR <- "MPRA_count_exp6_analysis3"

source(sprintf("%s/scripts/%s/00_starter.R", MPRA_DIR, ANALYSIS_DIR))
source(sprintf("%s/scripts/%s/00_load_associations.R", MPRA_DIR, ANALYSIS_DIR))




RUN_ID <- "RUN3_Z2_nBC10"
CRITERIA_ACTIVE <- "FDR5_scrambled5pct_FC0.2_GCnorm"

cmd <- commandArgs(trailingOnly = TRUE)
print(cmd)

for (i in seq_along(cmd)) {
  if (cmd[i] == "--run_dir" || cmd[i] == "-r") {
    RUN_ID <- cmd[i + 1]
  }
	if (cmd[i] == "--criteria" || cmd[i] == "-r") {
    CRITERIA_ACTIVE <- cmd[i + 1]
  }
}
source(sprintf("%s/scripts/%s/03_00_parameter_activity.R", MPRA_DIR, ANALYSIS_DIR))


IN_DIR <- sprintf("%s/data/%s/03a_aggMPRA_analyse/%s/", MPRA_DIR, ANALYSIS_DIR, RUN_ID)
ACTIVE_DIR <- sprintf("%s/data/%s/03b_computeFDR/%s/Activity/%s/", MPRA_DIR, ANALYSIS_DIR, RUN_ID, CRITERIA_ACTIVE)
dir.create(ACTIVE_DIR, recursive = TRUE)

FIGURE_DIR <- sprintf("%s/figures/%s/%s", MPRA_DIR, ANALYSIS_DIR, RUN_ID)

ACTIVITY_DIR <- sprintf("%s/02_oligo_activity/%s", FIGURE_DIR, CRITERIA_ACTIVE)
dir.create(sprintf("%s/PCs_activity/", ACTIVITY_DIR), showWarnings = FALSE, recursive = TRUE) 
dir.create(sprintf("%s/SupTables/", ACTIVITY_DIR), showWarnings = FALSE, recursive = TRUE) 

tic("loading oligo & SNP annotations")
# load annotation of oligos 
oligo_source <- fread(sprintf("%s/data/%s/00_oligo_annot_v2.txt", MPRA_DIR, ANALYSIS_DIR))

# load annotations of SNPs
SNP_annot <- fread(sprintf("%s/data/%s/00_SNP_annot_v3.txt", MPRA_DIR, ANALYSIS_DIR))
toc()
# required libraries
library(lme4)
library(umap)
library(pcaMethods)

source(sprintf("%s/scripts/%s/02z__define_includedCRS.R", MPRA_DIR, ANALYSIS_DIR))



cat('PCA , oligo activity\n')
################################################################################
########   PCA                                                     #######
################################################################################

oligo_activity_perm <- fread(sprintf("%s/all_oligos_annotated_perm__%s.tsv.gz", ACTIVE_DIR, CRITERIA_ACTIVE))
oligo_activity_perm_ctc <- oligo_activity_perm[power == 0 & boot == 0 & ANALYSIS_SUBTYPE == "celltype_cond", ]

oligo_activity_obs <- fread(sprintf("%s/all_oligos_annotated__%s.tsv.gz", ACTIVE_DIR, CRITERIA_ACTIVE))
# oligo_activity_obs <- oligo_activity_obs[oligo %chin% tested_and_ctrl_oligos_final$oligo, ]
oligo_activity_obs_ctc <- oligo_activity_obs[power == 0 & boot == 0 & ANALYSIS_SUBTYPE == "celltype_cond", ]
oligo_activity_obs[, logActivity := log2(alpha_CRITERIA)]

####  create per sample activity matrix
activity_matrix <- dcast(oligo_activity_obs[ANALYSIS_SUBTYPE == "sample" & excluded == FALSE, ], oligo ~ ANALYSIS_NAME, value.var = "logActivity")
oligo_names <- activity_matrix$oligo
activity_matrix <- as.matrix(activity_matrix[, mget(condition_summary_reps$analysis_name)])
rownames(activity_matrix) <- oligo_names

####  remove NAs
library(pcaMethods)
na.remove <- apply(activity_matrix, 1, function(x) mean(is.na(x)) > .9)
sum(na.remove) # 11 oligo with >90% NA were removed from data
PCs_activity <- pca(t(activity_matrix[!na.remove, ]), scale = "uv", center = TRUE, nPcs = ncol(activity_matrix))

# generate observed and permuted Data 
activity_matrix_noNA <- PCs_activity@completeObs@.Data

activity_matrix_noNA_null <- apply(activity_matrix_noNA, 2, sample)
PCs_activity_null <- pca(apply(activity_matrix_noNA_null, 2, sample), scale = "uv", center = TRUE, nPcs = ncol(activity_matrix))

# remove/condition cell type effect
activity_matrix_res <- apply(activity_matrix_noNA, 2, function(x) {
  lm(x ~ condition_summary_reps$COND_ID)$res
})
PCs_activity_residuals <- pca(activity_matrix_res, scale = "uv", center = TRUE, nPcs = ncol(activity_matrix))


activity_matrix_null_res <- apply(activity_matrix_noNA_null, 2, function(x) {
  lm(x ~ condition_summary_reps$COND_ID)$res
})
PCs_activity_residuals_null <- pca(activity_matrix_null_res, scale = "uv", center = TRUE, nPcs = ncol(activity_matrix))


activity_matrix_noNA_order <- activity_matrix_noNA[technical_covariates$SETUP_ID,]
technical_covariates[,COND_ID:=paste(gsub('-ACE2','',celltype),condition)]

########@ test Batch correction 
#x=activity_matrix_noNA_order[,1]
#summary(lmer(x ~ (1|celltype) + (1|COND_ID) + (1|plasmid_library) + (1|lentivirus_preparation),data=technical_covariates))
#as.data.table(summary(lmer(x ~ (1|celltype) + (1|COND_ID) + (1|plasmid_library) + (1|lentivirus_preparation),data=technical_covariates))$varcor)[,vcov/sum(vcov)]

# correct batch effect and repeat correction
correctBatch=function(x){
  mod=lmer(x ~ (1|celltype) + (1|COND_ID) + (1|plasmid_library) + (1|lentivirus_preparation),data=technical_covariates); x-ranef(mod)$plasmid_library
  lib_effects <- setNames(ranef(mod)$plasmid_library$'(Intercept)',rownames(ranef(mod)$plasmid_library))
  batch_effect <- setNames(ranef(mod)$lentivirus_preparation$'(Intercept)',rownames(ranef(mod)$lentivirus_preparation))
  x-lib_effects[technical_covariates$plasmid_library]-batch_effect[technical_covariates$lentivirus_preparation]
}
Batch_corrected <- apply(activity_matrix_noNA_order, 2, correctBatch)

PCs_activity_corrected <- pca(Batch_corrected, scale = "uv", center = TRUE, nPcs = ncol(activity_matrix))

## estimate varianbce explained by Batch, cell type and condition
Pct_variance <- apply(activity_matrix_noNA_order, 2, function(alpha){as.data.table(summary(lmer(alpha ~ (1|celltype) + (1|COND_ID) + (1|plasmid_library) + (1|lentivirus_preparation),data=technical_covariates))$varcor)[,setNames(vcov/sum(vcov),grp)]})
Pct_variance <- data.table(oligo=colnames(Pct_variance),t(Pct_variance))[,.(oligo,lentivirus_preparation,plasmid_library,celltype,condition=COND_ID,Residual)]
fwrite(Pct_variance,file=sprintf("%s/PCs_activity/SupTable1e_Variance_decomposition.tsv", ACTIVITY_DIR),sep='\t')

oligo_order=fread(sprintf("%s/figures/%s/Z_figures/SupTable_1/SupTable1c_tested_oligos.tsv", MPRA_DIR, ANALYSIS_DIR),sep='\t')
oligo_order[ , oligo_order:=1:.N]
Pct_variance=merge(oligo_order[,.(crsID,oligo_order)],Pct_variance,by.x="crsID",by.y='oligo')[order(oligo_order),-'oligo_order']
fwrite(Pct_variance,file=sprintf("%s/SupTables/SupTable1e_Variance_decomposition.tsv", ACTIVITY_DIR),sep='\t')



PC_data <- data.table(k = 1:42, observed = PCs_activity@R2, null = PCs_activity_null@R2, resid = PCs_activity_residuals@R2, null_resid = PCs_activity_residuals_null@R2)
PC_scores <- as.data.table(melt(PCs_activity@scores))
colnames(PC_scores) <- c("sample", "PC_k", "PCscore")
PC_scores <- merge(PC_scores, condition_summary_reps, by.x = "sample", by.y = "analysis_name")
PC_score_wide <- dcast(PC_scores[as.numeric(PC_k) < 10, ], sample + celltype + condition + COND_ID + experiment + barcode_lib ~ PC_k, value.var = "PCscore")
fwrite(PC_data,file=sprintf("%s/PCs_activity/PC_R2.tsv", ACTIVITY_DIR),sep='\t')
fwrite(PC_score_wide,file=sprintf("%s/PCs_activity/PC_score_wide.tsv", ACTIVITY_DIR),sep='\t')

##################### plot PCA results #####################
# plot R2
pdf(sprintf("%s/PCs_activity/PC_activity_R2.pdf", ACTIVITY_DIR), height = 5, width = 4)
p <- ggplot(melt(PC_data, id.vars = "k"), aes(x = k, y = value, group = variable, colour = variable)) +
  geom_line() +
  theme_plot()
print(p)
dev.off()

# plot PC score by condID
pdf(sprintf("%s/PCs_activity/PC_activity_scores_by_condID.pdf", ACTIVITY_DIR), height = 7, width = 7)
p <- ggplot(PC_scores[as.numeric(PC_k) < 10, ], aes(x = COND_ID, y = PCscore, fill = COND_ID)) +
  geom_boxplot() +
  facet_wrap(~PC_k) +
  theme_plot(rotate.x = 90)
p <- p + scale_fill_manual(values = color_setup_simplified_norep)
print(p)
dev.off()

# plot PC scores PCs 1-2
pdf(sprintf("%s/PCs_activity/PC_activity_all_PC12.pdf", ACTIVITY_DIR), height = 7, width = 7)
p <- ggplot(PC_score_wide, aes(x = PC1, y = PC2, col = COND_ID)) +
  geom_point() +
  theme_plot()
p <- p + scale_color_manual(values = color_setup_simplified_norep)
print(p)
dev.off()


##################### plot PCA results, batch corrected #####################
PC_scores_corrected <- as.data.table(melt(PCs_activity_corrected@scores))
colnames(PC_scores_corrected ) <- c("sample", "PC_k", "PCscore")
PC_scores_corrected  <- merge(PC_scores_corrected , condition_summary_reps, by.x = "sample", by.y = "analysis_name")
PC_score_corrected_wide <- dcast(PC_scores_corrected [as.numeric(PC_k) < 10, ], sample + celltype + condition + COND_ID + experiment + barcode_lib ~ PC_k, value.var = "PCscore")

fwrite(PC_score_corrected_wide,file=sprintf("%s/PCs_activity/PC_score_corrected_wide.tsv", ACTIVITY_DIR),sep='\t')

pdf(sprintf("%s/PCs_activity/PC_activity_scores_by_condID_batchCorrected.pdf", ACTIVITY_DIR), height = 7, width = 7)
p <- ggplot(PC_scores_corrected[as.numeric(PC_k) < 10, ], aes(x = COND_ID, y = PCscore, fill = COND_ID)) +
  geom_boxplot() +
  facet_wrap(~PC_k) +
  theme_plot(rotate.x = 90)
p <- p + scale_fill_manual(values = color_setup_simplified_norep)
print(p)
dev.off()

# plot PC scores corrected PCs 1-2
pdf(sprintf("%s/PCs_activity/PC_activity_all_PC12_batchCorrected.pdf", ACTIVITY_DIR), height = 7, width = 7)
p <- ggplot(PC_score_corrected_wide, aes(x = PC1, y = PC2, col = COND_ID)) +
  geom_point() +
  theme_plot()
p <- p + scale_color_manual(values = color_setup_simplified_norep)
print(p)
dev.off()
# plot PC scores PCs 3-4
pdf(sprintf("%s/PCs_activity/PC_activity_all_PC34_batchCorrected.pdf", ACTIVITY_DIR), height = 7, width = 7)
p <- ggplot(PC_score_corrected_wide, aes(x = PC3, y = PC4, col = COND_ID)) +
  geom_point() +
  theme_plot()
p <- p + scale_color_manual(values = color_setup_simplified_norep)
print(p)
dev.off()
# plot PC scores PCs 5-6
pdf(sprintf("%s/PCs_activity/PC_activity_all_PC56_batchCorrected.pdf", ACTIVITY_DIR), height = 7, width = 7)
p <- ggplot(PC_score_corrected_wide, aes(x = PC5, y = PC6, col = COND_ID)) +
  geom_point() +
  theme_plot()
p <- p + scale_color_manual(values = color_setup_simplified_norep)
print(p)
dev.off()
# plot PC scores PCs 7-8
pdf(sprintf("%s/PCs_activity/PC_activity_all_PC78_batchCorrected.pdf", ACTIVITY_DIR), height = 7, width = 7)
p <- ggplot(PC_score_corrected_wide, aes(x = PC7, y = PC8, col = COND_ID)) +
  geom_point() +
  theme_plot()
p <- p + scale_color_manual(values = color_setup_simplified_norep)
print(p)
dev.off()

##################### plot UMPA, no batch correction #####################
set.seed(12345)

activity.umap <- umap(activity_matrix_noNA, n_neighbors = 15, min_dist = .2, n_components = 2)
UMAP <- as.data.table(activity.umap$layout)
colnames(UMAP) <- c("UMAP1", "UMAP2")
UMAP[,sample :=rownames(activity_matrix_noNA)]
UMAP <- merge(UMAP, condition_summary_reps, by.x = "sample", by.y = "analysis_name")


# plot umap scores 
pdf(sprintf("%s/PCs_activity/umap_all_PC12_raw.pdf", ACTIVITY_DIR), height = 4, width = 3)
p <- ggplot(UMAP, aes(x = UMAP1, y = UMAP2, col = COND_ID)) +
  geom_point(alpha=.8) +
  theme_plot(fontsize=12)
p <- p + scale_color_manual(values = color_setup_simplified_norep) + guides(colour = guide_legend(override.aes = list(size = 2), ncol=3))
print(p)
dev.off()

##################### plot UMAP, batch corrected #####################
set.seed(123)
activity.umap <- umap(Batch_corrected, n_neighbors = 15, min_dist = .2, n_components = 2)
UMAP <- as.data.table(activity.umap$layout)
colnames(UMAP) <- c("UMAP1", "UMAP2")
UMAP[,sample :=rownames(Batch_corrected)]
UMAP <- merge(UMAP, condition_summary_reps, by.x = "sample", by.y = "analysis_name")

# plot umap scores 
pdf(sprintf("%s/PCs_activity/umap_all_PC12_batchcorrected.pdf", ACTIVITY_DIR), height = 4, width = 3)
p <- ggplot(UMAP, aes(x = UMAP1, y = UMAP2, col = COND_ID)) +
  geom_point(alpha=.8) +
  theme_plot(fontsize=12)
p <- p + scale_color_manual(values = color_setup_simplified_norep) + guides(colour = guide_legend(override.aes = list(size = 2), ncol=3))
print(p)
dev.off()

##################### plot UMAP, batch corrected, higher min dist #####################
set.seed(123)
activity.umap <- umap(Batch_corrected, n_neighbors = 15, min_dist = .5, n_components = 2)
UMAP <- as.data.table(activity.umap$layout)
colnames(UMAP) <- c("UMAP1", "UMAP2")
UMAP[,sample :=rownames(Batch_corrected)]
UMAP <- merge(UMAP, condition_summary_reps, by.x = "sample", by.y = "analysis_name")

# plot umap scores 
pdf(sprintf("%s/PCs_activity/umap_all_PC12_batchcorrected_minDist.5.pdf", ACTIVITY_DIR), height = 4, width = 3)
p <- ggplot(UMAP, aes(x = UMAP1, y = UMAP2, col = COND_ID)) +
  geom_point(alpha=.8) +
  theme_plot(fontsize=12)
p <- p + scale_color_manual(values = color_setup_simplified_norep) + guides(colour = guide_legend(override.aes = list(size = 2), ncol=3))
print(p)
dev.off()

###################### plot UMPA per cellline, batch corrected, higher min dist #####################

UMAP=list()
for(celline in condition_summary_reps[,unique(celline)]){
	set.seed(123)
	activity.umap <- umap(Batch_corrected[grepl(celline,rownames(Batch_corrected)),], n_neighbors = 5,  n_components = 2)
	UMAP[[celline]] <- as.data.table(activity.umap$layout)
colnames(UMAP[[celline]]) <- c("UMAP1", "UMAP2")
UMAP[[celline]][,sample :=rownames(Batch_corrected)[grepl(celline,rownames(Batch_corrected))]]
UMAP[[celline]] <- merge(UMAP[[celline]], condition_summary_reps, by.x = "sample", by.y = "analysis_name")
}
UMAP=rbindlist(UMAP)

# plot umap scores 
pdf(sprintf("%s/PCs_activity/umap_all_PC12_batchcorrected_minDist.5_percelline.pdf", ACTIVITY_DIR), height = 4, width = 7)
p <- ggplot(UMAP, aes(x = UMAP1, y = UMAP2, col = COND_ID))
p <- p + geom_point(alpha=.8) #+ stat_ellipse(geom = "polygon", aes(fill = COND_ID), alpha = 0.25)
p <- p + theme_plot(fontsize=12) + facet_grid(cols=vars(celline))
p <- p + scale_color_manual(values = color_setup_simplified_norep) + guides(colour = guide_legend(override.aes = list(size = 2), ncol=3))
print(p)
dev.off()

#############################################################################
#############################################################################
#############################################################################

#### possible interpretations of variability between residuals

PC_scores_resid <- as.data.table(melt(PCs_activity_residuals@scores))
colnames(PC_scores_resid) <- c("sample", "PC_k", "PCscore")
PC_scores_resid <- merge(PC_scores_resid, condition_summary_reps, by.x = "sample", by.y = "analysis_name")

PC_score_resid_wide <- dcast(PC_scores_resid[as.numeric(PC_k) < 10, ], sample + celltype + condition + COND_ID + experiment + barcode_lib ~ PC_k, value.var = "PCscore")


pdf(sprintf("%s/PCs_activity/PC_activity_residual_PC12.pdf", ACTIVITY_DIR), height = 7, width = 7)
p <- ggplot(PC_score_resid_wide, aes(x = PC1, y = PC2, col = COND_ID)) +
  geom_point() +
  theme_plot()
p <- p + scale_color_manual(values = color_setup_simplified_norep)
print(p)
dev.off()
# plot PC scores PCs 3-4
pdf(sprintf("%s/PCs_activity/PC_activity_residual_PC34.pdf", ACTIVITY_DIR), height = 7, width = 7)
p <- ggplot(PC_score_resid_wide, aes(x = PC3, y = PC4, col = COND_ID)) +
  geom_point() +
  theme_plot()
p <- p + scale_color_manual(values = color_setup_simplified_norep)
print(p)
dev.off()


loading_residuals_PCs <- data.table(oligo = rownames(PCs_activity_residuals@loadings), PCs_activity_residuals@loadings[, 1:10])
TFscores <- fread(sprintf("%s/data/%s/00_TFscore/0_TFscore_all.tsv.gz", MPRA_DIR, ANALYSIS_DIR))
loading_residuals_PCs_withTF <- merge(loading_residuals_PCs, TFscores, by = "oligo")

loading_residuals_PCs_withTF[, .(COR = cor(PC1, score, method = "s"), P = cor.test(PC1, score, method = "s")$p.value), by = .(matrix_id, names)][order(P)][1:20, ]

loading_residuals_PCs_withTF[, .(COR = cor(PC1, score, method = "s"), P = cor.test(PC1, score, method = "s")$p.value), by = .(matrix_id, names)][order(P)][1:20, ]
#     matrix_id       names         COR            P
#        <char>      <char>       <num>        <num>
#  1:  MA0033.1       FOXL1 -0.07317933 5.636319e-13
#  2:  MA0670.1        NFIA  0.07282960 7.259405e-13
#  3:  MA1540.1       NR5A1  0.07088947 2.892581e-12
#  4:  MA0025.1       NFIL3 -0.07071121 3.278354e-12
#  5:  MA0161.2        NFIC  0.06991791 5.701769e-12
#  6:  MA0843.1         TEF -0.06780790 2.412596e-11
#  7:  MA0748.2         YY2  0.06578919 9.213762e-11
#  8:  MA0043.2         HLF -0.06504262 1.497445e-10
#  9:  MA0639.1         DBP -0.06490552 1.636182e-10
# 10:  MA0592.1       ESRRA  0.06389120 3.133929e-10
# 11:  MA0713.1      PHOX2A -0.06380860 3.302809e-10
# 12:  MA0130.1     ZNF354C  0.06208660 9.720151e-10
# 13:  MA0106.1        TP53  0.06165438 1.268827e-09
# 14:  MA1579.1      ZBTB26  0.06104672 1.839910e-09
# 15:  MA0025.2       NFIL3 -0.06096804 1.930107e-09
# 16:  MA0703.1       LMX1B -0.06084321 2.082120e-09
# 17:  MA0681.2      PHOX2B -0.06067247 2.309058e-09
# 18:  MA0892.1        GSX1 -0.05940098 4.945490e-09
# 19:  MA0259.1 ARNT::HIF1A  0.05933039 5.156778e-09
# 20:  MA0592.3       ESRRA  0.05932214 5.182034e-09

##### chat GPT interpretations of PC1 of residuals
# * Many of these transcription factors (e.g. TEF, HLF, DBP, NFIL3) are closely linked to circadian rhythm regulation — this suggests a potential temporal component to the variability you're seeing in your MPRA replicates.
# * Several factors (e.g. ARNT::HIF1A, ESRRA) regulate metabolic pathways and cellular stress responses, which can fluctuate depending on cellular state or environmental conditions.
# * Some, like TP53, are key stress response regulators, which might be activated differently across replicates due to variable experimental stress or subtle differences in conditions.

loading_residuals_PCs_withTF[, .(COR = cor(PC2, score, method = "s"), P = cor.test(PC2, score, method = "s")$p.value), by = .(matrix_id, names)][order(P)][1:20, ]
#     matrix_id  names         COR            P
#        <char> <char>       <num>        <num>
#  1:  MA0080.5   SPI1 -0.06802565 2.083030e-11
#  2:  MA0762.1   ETV2 -0.06129203 1.584263e-09
#  3:  MA0081.2   SPIB -0.06092325 1.983364e-09
#  4:  MA0156.2    FEV -0.05867358 7.593284e-09
#  5:  MA0098.3   ETS1 -0.05758808 1.426435e-08
#  6:  MA0475.2   FLI1 -0.05708746 1.900616e-08
#  7:  MA0764.1   ETV4 -0.05649365 2.663158e-08
#  8:  MA0080.2   SPI1 -0.05633336 2.915359e-08
#  9:  MA0474.2    ERG -0.05517782 5.557081e-08
# 10:  MA0473.1   ELF1 -0.05488708 6.523279e-08
# 11:  MA0761.1   ETV1 -0.05477726 6.929052e-08
# 12:  MA0765.1   ETV5 -0.05223314 2.715779e-07
# 13:  MA0641.1   ELF4 -0.05151983 3.939534e-07
# 14:  MA0156.1    FEV -0.05142818 4.130958e-07
# 15:  MA0759.1   ELK3 -0.05094353 5.301994e-07
# 16:  MA0598.3    EHF -0.05058985 6.352305e-07
# 17:  MA0745.2  SNAI2  0.04964105 1.025570e-06
# 18:  MA0080.4   SPI1 -0.04871399 1.624329e-06
# 19:  MA0028.1   ELK1 -0.04865262 1.674061e-06
# 20:  MA0136.2   ELF5 -0.04838490 1.908629e-06
#     matrix_id  names         COR            P

cat(loading_residuals_PCs_withTF[, .(COR = cor(PC2, score, method = "s"), P = cor.test(PC2, score, method = "s")$p.value), by = .(matrix_id, names)][order(P)][1:20, names])
# SPI1 ETV2 SPIB FEV ETS1 FLI1 ETV4 SPI1 ERG ELF1 ETV1 ETV5 ELF4 FEV ELK3 EHF SNAI2 SPI1 ELK1 ELF5

##### chat GPT interpretations of PC2 of residuals
# * ETS Family Dominance: Most of these TFs (SPI1, ETVs, ERG, ELF1, etc.) belong to the ETS family, which are major regulators of hematopoiesis, immune cell function, and vascular development.
# * Redundancy and Variability: The repeated appearance of SPI1 (PU.1) suggests it might be a particularly strong driver of variability between replicates, perhaps because it can respond dynamically to subtle differences in cell state (like differentiation status or immune activation).
# * Development and Plasticity: Factors like SNAI2, ETV4/5, and ELK3 are involved in cell migration, EMT, and plasticity, which can contribute to variable regulatory outcomes across replicates in cell-based assays.
# * Stress and Signaling: ELK1 (and related ELK factors) are immediate early response genes activated by MAPK pathways—making them sensitive to subtle differences in culture conditions, stress, or serum factors.

loading_residuals_PCs_withTF[, .(COR = cor(PC3, score, method = "s"), P = cor.test(PC3, score, method = "s")$p.value), by = .(matrix_id, names)][order(P)][1:20, ]
#     matrix_id  names        COR            P
#        <char> <char>      <num>        <num>
#  1:  MA0599.1   KLF5 0.08251681 4.228850e-16
#  2:  MA0516.2    SP2 0.07800287 1.523001e-14
#  3:  MA0162.2   EGR1 0.07528965 1.193718e-13
#  4:  MA1653.1 ZNF148 0.07413677 2.802104e-13
#  5:  MA0079.3    SP1 0.07352878 4.371834e-13
#  6:  MA0516.1    SP2 0.07345027 4.629095e-13
#  7:  MA1522.1    MAZ 0.07328351 5.225839e-13
#  8:  MA0753.1 ZNF740 0.07227510 1.081696e-12
#  9:  MA0753.2 ZNF740 0.07151906 1.854324e-12
# 10:  MA0528.1 ZNF263 0.06876872 1.257532e-11
# 11:  MA0039.4   KLF4 0.06856073 1.449107e-11
# 12:  MA0471.2   E2F6 0.06819982 1.851515e-11
# 13:  MA0746.1    SP3 0.06803525 2.069550e-11
# 14:  MA0079.2    SP1 0.06787800 2.301287e-11
# 15:  MA1578.1  VEZF1 0.06787508 2.305823e-11
# 16:  MA0746.2    SP3 0.06741886 3.132991e-11
# 17:  MA0079.4    SP1 0.06621888 6.950099e-11
# 18:  MA0471.1   E2F6 0.06550100 1.112064e-10
# 19:  MA0741.1  KLF16 0.06541758 1.174116e-10
# 20:  MA0685.1    SP4 0.06309056 5.198275e-10
##### chat GPT interpretations of PC3 of residuals

# * GC-Rich Promoter Regulation
#  A significant portion of your list (SP1, SP2, SP3, SP4, MAZ, ZNF263, KLFs) preferentially bind GC-rich regions, which are common in CpG islands.
#  If your MPRA library is enriched for GC-rich sequences, variability could be driven by differences in chromatin accessibility, DNA methylation, or competition between these factors.

# * Immediate Early Genes and Stress Response: EGR1 and some KLFs respond rapidly to environmental changes (serum, stress), so even small differences in experimental conditions can affect their activity.

# * Cell Cycle and Differentiation Control :E2F6, KLF4, and KLF5 regulate cell cycle genes and can affect proliferation, potentially introducing variability between replicates depending on subtle differences in cell confluence or passage number.



### reading outputs:

# MPRA_DIR <- "/pasteur/helix/projects/evo_immuno_pop/MPRA"
# ANALYSIS_DIR <- "MPRA_count_exp6_analysis2"
# RUN_ID <- "RUN2_Z2_nBC10"
# OUT_DIR <- sprintf("%s/data/%s/03b_computeFDR/%s/", MPRA_DIR, ANALYSIS_DIR, RUN_ID)
# oligo_activity_obs <- fread(sprintf("%s/all_oligos_annotated.tsv.gz", OUT_DIR))
# CRS_activity <- fread(sprintf("%s/CRS_activity__celltype_cond.tsv.gz", OUT_DIR))

# A note on the analysis startegy for power related analyses :
# For each CRS, we are going to :
# 1. report the number of barcode in the library and
# 2. estimate the power to detect activity at the reported level of activity
#    2a. this can be computed for the Winner's curse adjusted beta values
#        (use FIQT: F DR I nverse Q uantile T ransformation).
#         doi.org/10.1093/bioinformatics/btw303
#    2b. perhaps use local FDR from fdrtool rather than FDR ?

# 3. test whether probability of having an emVar varies across tissue/conditionb, adjiusting for difference in power.
cat('\nAll done!\n')
q('no')