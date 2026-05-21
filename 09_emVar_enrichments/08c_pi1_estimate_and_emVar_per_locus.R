
MPRA_DIR <- "/pasteur/helix/projects/evo_immuno_pop/MPRA"
ANALYSIS_DIR <- "MPRA_count_exp6_analysis3"

source(sprintf("%s/scripts/%s/00_starter.R", MPRA_DIR, ANALYSIS_DIR))
source(sprintf("%s/scripts/%s/00_load_associations.R", MPRA_DIR, ANALYSIS_DIR))
source(sprintf("%s/scripts/misc_plots.R", MPRA_DIR))

RUN_ID <- "RUN3_Z2_nBC10"
CRITERIA_ACTIVE <- "FDR5_scrambled5pct_FC0.2_GCnorm"
CRITERIA_EMVARS <- "EmVar_FDR5_FC.2"
FILTER_ACTIVE <- TRUE

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

if (!FILTER_ACTIVE) {
  CRITERIA_ACTIVE_OUT <- paste0("noActivityFilter_", CRITERIA_ACTIVE)
} else {
  CRITERIA_ACTIVE_OUT <- CRITERIA_ACTIVE
}

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

# load gene targets
oligo_target_collapsed <- fread(sprintf("%s/data/%s/00_oligo_targets_collapsed_v4.txt", MPRA_DIR, ANALYSIS_DIR), sep = "\t")
oligo_target <- fread(sprintf("%s/data/%s/00_oligo_targets_v4.txt", MPRA_DIR, ANALYSIS_DIR))
crs_targets <- unique(oligo_target_collapsed[, -"oligo"])
NearestGene <- unique(oligo_target[method == "NearestGene" & gene_type != "", .(posID = crsID, TargetGene)])
locus_genes <- merge(SNP_annot[,.(posID, introgression_locus )],NearestGene, by="posID")
locus_genes <- locus_genes[introgression_locus!="",.(locus_genes=paste(unique(TargetGene),collapse='/')),by=introgression_locus]
# load genes <100kb of each crs
NearbyGenes <- fread(sprintf("%s/data/%s/00_oligo_nearbyGenes_100kb.txt.gz", MPRA_DIR, ANALYSIS_DIR))

toc()

library(ggrepel)

##### active parameters
source(sprintf("%s/scripts/%s/03_00_parameter_activity.R", MPRA_DIR, ANALYSIS_DIR))
##### emVars parameters
source(sprintf("%s/scripts/%s/05_00_parameter_emVars.R", MPRA_DIR, ANALYSIS_DIR))

IN_DIR <- sprintf("%s/data/%s/03a_aggMPRA_analyse/%s/", MPRA_DIR, ANALYSIS_DIR, RUN_ID)
ACTIVE_DIR <- sprintf("%s/data/%s/03b_computeFDR/%s/Activity/%s/", MPRA_DIR, ANALYSIS_DIR, RUN_ID, CRITERIA_ACTIVE)
EMVAR_DIR <- sprintf("%s/data/%s/03b_computeFDR/%s/emVars/%s/%s", MPRA_DIR, ANALYSIS_DIR, RUN_ID, CRITERIA_EMVARS, CRITERIA_ACTIVE_OUT)

dir.create(EMVAR_DIR, recursive = TRUE)
FIGURE_DIR <- sprintf("%s/figures/%s/%s/04_emVars/%s/%s", MPRA_DIR, ANALYSIS_DIR, RUN_ID, CRITERIA_EMVARS, CRITERIA_ACTIVE_OUT)
#FIGURE_DIR <- sprintf("%s/figures/%s/04_emVars/%s/%s", MPRA_DIR, ANALYSIS_DIR, CRITERIA_EMVARS, CRITERIA_ACTIVE_OUT)

dir.create(FIGURE_DIR, recursive = TRUE)
dir.create(sprintf("%s/SupTables/", FIGURE_DIR), recursive = TRUE)
cat("\nprinting output in:", FIGURE_DIR)


emVARs_obs_ctc <- fread(sprintf("%s/all_emVars_annotated_celltype__%s.tsv", EMVAR_DIR, CRITERIA_ACTIVE_OUT))

#######################################################################################################
########################### pi1_estim by COND ID & Source  ############################################
#######################################################################################################

pi1_estim <- emVARs_obs_ctc[Introgression_source_top%in%c("Vindija",'Denisova'),
      .(nMarg=sum(pval_emp<.05),
        nNS=sum(pval_emp>.05), 
        pi1=1-mean(pval_emp>.05)/0.95, 
        pi1.se=sd(replicate(1000,1-mean(sample(pval_emp,.N,replace=TRUE)>.05)/.95)))
        , keyby=.(COND_ID=factor(COND_ID,names(color_setup_simplified_norep)),Introgression_source_top)]

fwrite(pi1_estim, file = sprintf("%s/SupTables/pi1_estim_emVar_per_condition_Vindija_vs_Denisova.tsv", FIGURE_DIR), sep = "\t")

Fisher_results=pi1_estim[,fisher.test(.SD[,.(nNS,nMarg)]),keyby=.(COND_ID)]
fwrite(Fisher_results, file = sprintf("%s/SupTables/Fisher_results_percent_emVars_Vindija_vs_Denisova.tsv", FIGURE_DIR), sep = "\t")

pi1_estim <- fread(sprintf("%s/SupTables/pi1_estim_emVar_per_condition_Vindija_vs_Denisova.tsv", FIGURE_DIR))
Fisher_results <- fread(sprintf("%s/SupTables/Fisher_results_percent_emVars_Vindija_vs_Denisova.tsv", FIGURE_DIR))

Fisher_results <- Fisher_results[ ,.(conf.int=paste0('[',paste(round(conf.int,2),collapse=','),']')),keyby=.(COND_ID,OR=estimate,  pval=p.value)]
Fisher_results[,FDR:=p.adjust(pval,"fdr")]
pi1_estim_WithTest = merge(pi1_estim, Fisher_results, by='COND_ID')
pi1_estim_WithTest[!duplicated(COND_ID),FDR:=NA]
pi1_estim_WithTest[!duplicated(COND_ID),OR:=NA]
pi1_estim_WithTest[!duplicated(COND_ID),pval:=NA]
pi1_estim_WithTest[!duplicated(COND_ID),conf.int:=NA]

########################### Table 5a
fwrite(pi1_estim_WithTest, file = sprintf("%s/SupTables/SupTable5a_pi1_estim_emVar_per_condition_POP_x_source_signed.tsv", FIGURE_DIR), sep = "\t")


p <- ggplot(pi1_estim,aes(x=Introgression_source_top,y=pi1, fill=COND_ID, alpha=Introgression_source_top))
p <- p + geom_bar(stat='Identity',position='dodge') + ylab('Percentage of regulatory alleles (pi1 estimate)')
p <- p + geom_errorbar(aes(x=Introgression_source_top, ymin=pi1-1.96*pi1.se,ymax=pi1+1.96*pi1.se), width=.5) +facet_grid(cols=vars(COND_ID))
p <- p + scale_fill_manual(values=color_setup_simplified_norep)
p <- p + scale_alpha_manual(values=c('Vindija'=1,'Denisova'=0.7))+theme_plot(rotate.x=90,fontsize=11)
pdf(sprintf("%s/08a_pi1_estim_emVar_per_condition_Vindija_vs_Denisova.pdf", FIGURE_DIR), height = 4, width = 4)
print(p)
dev.off()

########################### Figure 5a
pi1_estim <-merge(pi1_estim,condition_summary,by='COND_ID')
  p <- ggplot(pi1_estim,aes(x=factor(ifelse(Introgression_source_top=='Vindija','N','D'),levels=c('N','D')),y=pi1, fill=COND_ID, alpha=ifelse(Introgression_source_top=='Vindija','N','D')))
  p <- p + geom_bar(stat='Identity',position='dodge') + ylab('Percentage of regulatory \nalleles (pi1 estimate)')
  p <- p + geom_errorbar(aes(x=factor(ifelse(Introgression_source_top=='Vindija','N','D'),levels=c('N','D')), ymin=pi1-1.96*pi1.se,ymax=pi1+1.96*pi1.se), width=.5)
  p <- p + facet_grid(cols=vars(celline, factor(condition, levels=condition_order)))
  p <- p + scale_fill_manual(values=color_setup_simplified_norep) + xlab('Archaic source')
  p <- p + scale_alpha_manual(values=c('N'=1,'D'=0.7)) + theme_plot(rotate.x=0,fontsize=9) + guides(fill=FALSE,alpha=FALSE)
  pdf(sprintf("%s/08a_pi1_estim_emVar_per_condition_Vindija_vs_Denisova_v2.pdf", FIGURE_DIR), height = 2.5, width = 5)
  print(p)
  dev.off()


pi1_estim_archaic_origin <- emVARs_obs_ctc[,
      .(nMarg=sum(pval_emp<.05),
        nNS=sum(pval_emp>.05), 
        pi1=1-mean(pval_emp>.05)/0.95, 
        pi1.se=sd(replicate(1000,1-mean(sample(pval_emp,.N,replace=TRUE)>.05)/.95)))
        , keyby=.(COND_ID=factor(COND_ID,names(color_setup_simplified_norep)),
								SOURCE=ifelse(grepl('archaic-derived',Introgression_scenario),
																	ifelse(Introgression_source_top%in%c('Vindija','Denisova'),
																						Introgression_source_top,'undetermined/Shared'),'reintrogressed'))]


########################### Figure 5a
pi1_estim_archaic_origin <-merge(pi1_estim_archaic_origin,condition_summary,by='COND_ID')
pi1_estim_archaic_origin[,SOURCE:=factor(ifelse(toupper(substr(SOURCE,1,1))=='V','N',toupper(substr(SOURCE,1,1))),levels=c('N','D','R','U'))]
  p <- ggplot(pi1_estim_archaic_origin,aes(x=SOURCE,y=pi1, fill=COND_ID, alpha=SOURCE))
  p <- p + geom_bar(stat='Identity',position='dodge') + ylab('Percentage of regulatory \nalleles (pi1 estimate)')
  p <- p + geom_errorbar(aes(x=SOURCE, ymin=pi1-1.96*pi1.se,ymax=pi1+1.96*pi1.se), width=.5)
  p <- p + facet_grid(cols=vars(celline, factor(condition, levels=condition_order)))
  p <- p + scale_fill_manual(values=color_setup_simplified_norep) + xlab('Archaic source')
  p <- p + scale_alpha_manual(values=c('N'=1,'D'=0.7,'R'=0.4)) + theme_plot(rotate.x=0,fontsize=9) + guides(fill=FALSE,alpha=FALSE)
  pdf(sprintf("%s/08a_pi1_estim_emVar_per_condition_Vindija_vs_Denisova_V2.1_origin.pdf", FIGURE_DIR), height = 2.5, width = 5)
  print(p)
  dev.off()


pi1_estim_archaic_origin <- emVARs_obs_ctc[,
      .(nMarg=sum(pval_emp<.05),
        nNS=sum(pval_emp>.05), 
        pi1=1-mean(pval_emp>.05)/0.95, 
        pi1.se=sd(replicate(1000,1-mean(sample(pval_emp,.N,replace=TRUE)>.05)/.95)))
        , keyby=.(COND_ID=factor(COND_ID,names(color_setup_simplified_norep)),
								SOURCE=ifelse(Introgression_source_top%in%c('Vindija','Denisova'),
																						Introgression_source_top,'undetermined/Shared'))]


pi1_estim_archaic_origin <- emVARs_obs_ctc[,
      .(nMarg=sum(pval_emp<.05),
        nNS=sum(pval_emp>.05), 
        pi1=1-mean(pval_emp>.05)/0.95, 
        pi1.se=sd(replicate(1000,1-mean(sample(pval_emp,.N,replace=TRUE)>.05)/.95)))
        , keyby=.(COND_ID=factor(COND_ID,names(color_setup_simplified_norep)),
								SOURCE=ifelse(Introgression_source_top%in%c('Vindija','Denisova'),
																						Introgression_source_top,'undetermined/Shared'),
								DERIVED=ifelse(grepl('archaic-derived',Introgression_scenario),'','_R'))]

########################### Figure 5a
pi1_estim_archaic_origin <-merge(pi1_estim_archaic_origin,condition_summary,by='COND_ID')
pi1_estim_archaic_origin[,SOURCE:=factor(paste0(ifelse(toupper(substr(SOURCE,1,1))=='V','N',toupper(substr(SOURCE,1,1))),DERIVED),levels=c('N','D','U','N_R','D_R','U_R'))]
  p <- ggplot(pi1_estim_archaic_origin,aes(x=SOURCE,y=pi1, fill=COND_ID, alpha=SOURCE))
  p <- p + geom_bar(stat='Identity',position='dodge') + ylab('Percentage of regulatory \nalleles (pi1 estimate)')
  p <- p + geom_errorbar(aes(x=SOURCE, ymin=pi1-1.96*pi1.se,ymax=pi1+1.96*pi1.se), width=.5)
  p <- p + facet_grid(cols=vars(celline, factor(condition, levels=condition_order)))
  p <- p + scale_fill_manual(values=color_setup_simplified_norep) + xlab('Archaic source')
  p <- p + scale_alpha_manual(values=c('N'=1,'D'=0.7,'R'=0.4)) + theme_plot(rotate.x=90,fontsize=9) + guides(fill=FALSE,alpha=FALSE)
  pdf(sprintf("%s/08a_pi1_estim_emVar_per_condition_Vindija_vs_Denisova_V2.2_origin.pdf", FIGURE_DIR), height = 2.5, width = 8)
  print(p)
  dev.off()

fwrite(pi1_estim, file = sprintf("%s/SupTables/pi1_estim_emVar_per_condition_Vindija_vs_Denisova.tsv", FIGURE_DIR), sep = "\t")

Fisher_results=pi1_estim[,fisher.test(.SD[,.(nNS,nMarg)]),keyby=.(COND_ID)]
fwrite(Fisher_results, file = sprintf("%s/SupTables/Fisher_results_percent_emVars_Vindija_vs_Denisova.tsv", FIGURE_DIR), sep = "\t")

pi1_estim <- fread(sprintf("%s/SupTables/pi1_estim_emVar_per_condition_Vindija_vs_Denisova.tsv", FIGURE_DIR))
Fisher_results <- fread(sprintf("%s/SupTables/Fisher_results_percent_emVars_Vindija_vs_Denisova.tsv", FIGURE_DIR))

Fisher_results <- Fisher_results[ ,.(conf.int=paste0('[',paste(round(conf.int,2),collapse=','),']')),keyby=.(COND_ID,OR=estimate,  pval=p.value)]
Fisher_results[,FDR:=p.adjust(pval,"fdr")]
pi1_estim_WithTest = merge(pi1_estim, Fisher_results, by='COND_ID')
pi1_estim_WithTest[!duplicated(COND_ID),FDR:=NA]
pi1_estim_WithTest[!duplicated(COND_ID),OR:=NA]
pi1_estim_WithTest[!duplicated(COND_ID),pval:=NA]
pi1_estim_WithTest[!duplicated(COND_ID),conf.int:=NA]

########################### Table 5a
fwrite(pi1_estim_WithTest, file = sprintf("%s/SupTables/SupTable5a_pi1_estim_emVar_per_condition_POP_x_source_signed.tsv", FIGURE_DIR), sep = "\t")







#######################################################################################################
########################### pi1_estim_signed by COND ID & Source  #####################################
#######################################################################################################
TH=.05
pi1_estim_signed <- emVARs_obs_ctc[Introgression_source_top%in%c("Vindija",'Denisova'),
      .(.N, 
        nMarg_pos=sum(pval_emp<TH & log2FC_archaic_vs_modern>0),
        nMarg_neg=sum(pval_emp<TH & log2FC_archaic_vs_modern<0),
        nNS=sum(pval_emp>TH), 
        pi1_pos=1-mean(pval_emp>TH | log2FC_archaic_vs_modern<0)/(1-TH/2), 
        pi1_pos.se=sd(replicate(1000,1-mean(sample(pval_emp>TH | log2FC_archaic_vs_modern<0,.N,replace=TRUE))/(1-TH/2))),
        pi1_neg=1-mean(pval_emp>TH | log2FC_archaic_vs_modern>0)/(1-TH/2), 
        pi1_neg.se=sd(replicate(1000,1-mean(sample(pval_emp>TH | log2FC_archaic_vs_modern>0,.N,replace=TRUE))/(1-TH/2)))),
        keyby=.(COND_ID=factor(COND_ID,names(color_setup_simplified_norep)), Introgression_source_top)]

total=emVARs_obs_ctc[,.(.N, 
        nMarg_pos=sum(pval_emp<TH & log2FC_archaic_vs_modern>0),
        nMarg_neg=sum(pval_emp<TH & log2FC_archaic_vs_modern<0),
        nNS=sum(pval_emp>TH)),by=.(COND_ID=factor(COND_ID,names(color_setup_simplified_norep)))]

pi1_estim_signed_v2=merge(pi1_estim_signed,total,by=c('COND_ID'),suffix=c('','_total'))
pi1_estim_signed_v2[, pval_pos:=fisher.test(matrix(c(nMarg_neg_total+nNS_total-nMarg_neg-nNS,nMarg_pos_total-nMarg_pos, nMarg_neg+nNS,nMarg_pos),2,2))$p.value, by=.(COND_ID, Introgression_source_top)]
pi1_estim_signed_v2[, pval_neg:=fisher.test(matrix(c(nMarg_pos_total+nNS_total-nMarg_pos-nNS,nMarg_neg_total-nMarg_neg, nMarg_pos+nNS,nMarg_neg),2,2))$p.value, by=.(COND_ID, Introgression_source_top)]

pi1_estim_signed_v3 <- pi1_estim_signed_v2[order(pmin(pval_pos,pval_neg))]
pi1_estim_signed_v3[,FDR_pos:=p.adjust(pval_pos,'fdr')]
pi1_estim_signed_v3[,FDR_neg:=p.adjust(pval_neg,'fdr')]
pi1_estim_signed_v3

Figdata <- pi1_estim_signed_v3
#Figdata[,class:=factor(paste(region, Source_introgression.byPOP), levels=c('EUR Vindija','EAS Vindija','SAS Vindija','Papuan Vindija','Papuan Denisova','Agta Vindija','Agta Denisova'))]
Figdata <- melt(Figdata,measure.vars=list(c('pi1_pos','pi1_neg'),c('pi1_pos.se','pi1_neg.se')),value.name=c('pi1','pi1.se'))
levels(Figdata$variable)=c('pos','neg')
Figdata[,class:=factor(paste0(ifelse(Introgression_source_top=='Vindija','N','D'),ifelse(variable=='pos','+','-')),levels=c('N+','N-','D+', 'D-'),)]

p <- ggplot(Figdata,aes(x=class,y=pi1, fill=COND_ID, alpha=Introgression_source_top))
p <- p + geom_bar(stat='Identity',position='dodge') + ylab('Percentage of regulatory alleles (pi1 estimate)')
p <- p + geom_errorbar(aes(x=class, ymin=pi1-1.96*pi1.se,ymax=pi1+1.96*pi1.se), width=.2) +facet_wrap(~COND_ID)
p <- p + scale_fill_manual(values=color_setup_simplified_norep)
p <- p + scale_alpha_manual(values=c('Vindija'=1,'Denisova'=0.7))+theme_plot(rotate.x=90,fontsize=11)
p <- p + guides(fill=FALSE) + xlab('')
pdf(sprintf("%s/08b_pi1_estim_emVar_per_condition_Vindija_vs_Denisova_signed.pdf", FIGURE_DIR), height = 8, width = 4)
print(p)
dev.off()

############################################################################################################
########################### pi1_estim by COND ID, POP & Source  ############################################
############################################################################################################

emVARs_obs_ctc_full <- merge(emVARs_obs_ctc, selected_annot,by=c('ID','posID'),allow.cartesian=TRUE,all.x=TRUE, suffix=c("",".byPOP"))
emVARs_obs_ctc_full[,region:=case_when(POP%in%c('GBR','IBS')~'EUR',
                                       POP%in%c('CHB','JPT')~'EAS',
                                       POP%in%c('STU','PJL')~'SAS',
                                       POP%in%c('Agta')~'AGT',
																			 POP%in%c('Papuan')~'PAP',
																			 TRUE~POP)]
emVARs_obs_ctc_full <- emVARs_obs_ctc_full[!duplicated(paste(posID,COND_ID,region)) & POP!='',]
pi1_estim <- emVARs_obs_ctc_full[,
      .(.N, 
        nMarg=sum(pval_emp<.05),
        nNS=sum(pval_emp>.05), 
        pi1=1-mean(pval_emp>.05)/0.95, 
        pi1.se=sd(replicate(1000,1-mean(sample(pval_emp,.N,replace=TRUE)>.05)/.95))),
        keyby=.(COND_ID=factor(COND_ID,names(color_setup_simplified_norep)), region, Source_introgression.byPOP)]

fwrite(pi1_estim, file = sprintf("%s/SupTables/pi1_estim_emVar_per_condition_POP_x_source.tsv", FIGURE_DIR), sep = "\t")

Figdata <- pi1_estim[N>100 & Source_introgression.byPOP%in%c('Vindija','Denisova'),]
Figdata[,class:=factor(paste(region, Source_introgression.byPOP), levels=c('EUR Vindija','EAS Vindija','SAS Vindija','PAP Vindija','PAP Denisova','AGT Vindija','AGT Denisova'))]
p <- ggplot(Figdata,aes(x=class,y=pi1, fill=COND_ID, alpha=Source_introgression.byPOP))
p <- p + geom_bar(stat='Identity',position='dodge') + ylab('Percentage of regulatory alleles (pi1 estimate)')
p <- p + geom_errorbar(aes(x=class, ymin=pi1-1.96*pi1.se,ymax=pi1+1.96*pi1.se), width=.5) +facet_wrap(~factor(COND_ID,levels=names(color_setup_simplified_norep)))
p <- p + scale_fill_manual(values=color_setup_simplified_norep)
p <- p + scale_alpha_manual(values=c('Vindija'=1,'Denisova'=0.7))+theme_plot(rotate.x=90,fontsize=11)
p <- p + guides(fill=FALSE) + xlab('')
pdf(sprintf("%s/08c_pi1_estim_emVar_per_condition_POP_x_source.pdf", FIGURE_DIR), height = 8, width = 4)
print(p)
dev.off()

################################################################################################################
########################### pi1_estim_signed by COND ID, POP & Source  ############################################
########################### Fig 5b and S16 ############################################
################################################################################################################

TH=.05
pi1_estim_signed <- emVARs_obs_ctc_full[,
      .(.N, 
        nMarg_pos=sum(pval_emp<TH & log2FC_archaic_vs_modern>0),
        nMarg_neg=sum(pval_emp<TH & log2FC_archaic_vs_modern<0),
        nNS=sum(pval_emp>TH), 
        pi1_pos=1-mean(pval_emp>TH | log2FC_archaic_vs_modern<0)/(1-TH/2), 
        pi1_pos.se=sd(replicate(1000,1-mean(sample(pval_emp>TH | log2FC_archaic_vs_modern<0,.N,replace=TRUE))/(1-TH/2))),
        pi1_neg=1-mean(pval_emp>TH | log2FC_archaic_vs_modern>0)/(1-TH/2), 
        pi1_neg.se=sd(replicate(1000,1-mean(sample(pval_emp>TH | log2FC_archaic_vs_modern>0,.N,replace=TRUE))/(1-TH/2)))),
        keyby=.(COND_ID=factor(COND_ID,names(color_setup_simplified_norep)), region, Source_introgression.byPOP)]

total=emVARs_obs_ctc[,.(.N, 
        nMarg_pos=sum(pval_emp<TH & log2FC_archaic_vs_modern>0),
        nMarg_neg=sum(pval_emp<TH & log2FC_archaic_vs_modern<0),
        nNS=sum(pval_emp>TH)),by=.(COND_ID=factor(COND_ID,names(color_setup_simplified_norep)))]

Figdata <- pi1_estim_signed[N>100 & Source_introgression.byPOP%in%c('Vindija','Denisova'),]
Figdata[,class:=factor(paste(region, ifelse(Source_introgression.byPOP=='Vindija','N','D'),sep=' - '), levels=c('EUR - N','EAS - N','SAS - N','PAP - N','PAP - D','AGT - N','AGT - D'))]
Figdata <- melt(Figdata,measure.vars=list(c('pi1_pos','pi1_neg'),c('pi1_pos.se','pi1_neg.se')),value.name=c('pi1','pi1.se'))
levels(Figdata$variable)=c('pos','neg')

########################### Fig S16
p <- ggplot(Figdata,aes(x=class,y=pi1, fill=COND_ID, alpha=Source_introgression.byPOP))
p <- p + geom_bar(stat='Identity',position='dodge') + ylab('Percentage of regulatory alleles (pi1 estimate)')
p <- p + geom_errorbar(aes(x=class, ymin=pi1-1.96*pi1.se,ymax=pi1+1.96*pi1.se), width=.5) 
p <- p + facet_grid(cols=vars(factor(COND_ID,names(color_setup_simplified_norep))),rows=vars(variable))
p <- p + scale_fill_manual(values=color_setup_simplified_norep)
p <- p + scale_alpha_manual(values=c('Vindija'=1,'Denisova'=0.7))+theme_plot(rotate.x=90,fontsize=10)
p <- p + guides(fill=FALSE) + xlab('')
pdf(sprintf("%s/08d_pi1_estim_emVar_per_condition_POP_x_source_signed.pdf", FIGURE_DIR), height = 5, width = 8)
print(p)
dev.off()

Figdata <- merge(Figdata, condition_summary[,.(COND_ID,condition)], by='COND_ID')
p <- ggplot(Figdata[COND_ID%in%c('A549_IAV','A549_TNFa')],aes(x=class,y=pi1, fill=COND_ID, alpha=Source_introgression.byPOP))
p <- p + geom_bar(stat='Identity',position='dodge') + ylab('Percentage of regulatory\n alleles (pi1 estimate)')
p <- p + geom_errorbar(aes(x=class, ymin=pi1-1.96*pi1.se,ymax=pi1+1.96*pi1.se), width=.5) + facet_grid(cols=vars(factor(ifelse(variable=='pos','increasing','decreasing'),c('increasing','decreasing'))),rows=vars(condition))
p <- p + scale_fill_manual(values=color_setup_simplified_norep)
p <- p + scale_alpha_manual(values=c('Vindija'=1,'Denisova'=0.7))+theme_plot(rotate.x=90,fontsize=10)
p <- p + guides(fill=FALSE, alpha=FALSE) + xlab('')
pdf(sprintf("%s/08d_pi1_estim_emVar_per_condition_POP_x_source_signed_subset.pdf", FIGURE_DIR), height = 2.5, width = 2.5)
print(p)
dev.off()

########################### Fig 5b
p <- p + facet_grid(cols=vars(condition,factor(ifelse(variable=='pos','increasing','decreasing'),c('increasing','decreasing'))))
pdf(sprintf("%s/08d_pi1_estim_emVar_per_condition_POP_x_source_signed_subset_v2.pdf", FIGURE_DIR), height = 2.5, width = 3.3)
print(p)
dev.off()

pi1_estim_signed_v2=merge(pi1_estim_signed,total,by=c('COND_ID'),suffix=c('','_total'))
pi1_estim_signed_v2[, pval_pos:=fisher.test(matrix(c(nMarg_neg_total+nNS_total-nMarg_neg-nNS,nMarg_pos_total-nMarg_pos, nMarg_neg+nNS,nMarg_pos),2,2))$p.value, by=.(COND_ID, region, Source_introgression.byPOP)]
pi1_estim_signed_v2[, pval_neg:=fisher.test(matrix(c(nMarg_pos_total+nNS_total-nMarg_pos-nNS,nMarg_neg_total-nMarg_neg, nMarg_pos+nNS,nMarg_neg),2,2))$p.value, by=.(COND_ID, region, Source_introgression.byPOP)]

pi1_estim_signed_v3 <- pi1_estim_signed_v2[order(pmin(pval_pos,pval_neg))][Source_introgression.byPOP%in%c('Denisova','Vindija') & N>100]
pi1_estim_signed_v3[,FDR_pos:=p.adjust(pval_pos,'fdr')]
pi1_estim_signed_v3[,FDR_neg:=p.adjust(pval_neg,'fdr')]
pi1_estim_signed_v3

########################### Table S5b
fwrite(pi1_estim_signed_v3, file = sprintf("%s/SupTables/SupTable5b_pi1_estim_emVar_per_condition_POP_x_source_signed.tsv", FIGURE_DIR), sep = "\t")

pi1_estim_signed_v3 <- pi1_estim_signed_v3[COND_ID%in%c('A549_IAV','A549_TNFa'),][order(pmin(pval_pos,pval_neg))][Source_introgression.byPOP%in%c('Denisova','Vindija') & N>100]
pi1_estim_signed_v3[,FDR_pos:=p.adjust(pval_pos,'fdr')]
pi1_estim_signed_v3[,FDR_neg:=p.adjust(pval_neg,'fdr')]
pi1_estim_signed_v3

################################################################################################################
########################### pi1_estim_signed by COND ID, POP & Source (adaptive only)  ############################################
################################################################################################################

TH=.05
pi1_estim_signed <- emVARs_obs_ctc_full[is_adaptive==TRUE,
      .(.N, 
        nMarg_pos=sum(pval_emp<TH & log2FC_archaic_vs_modern>0),
        nMarg_neg=sum(pval_emp<TH & log2FC_archaic_vs_modern<0),
        nNS=sum(pval_emp>TH), 
        pi1_pos=1-mean(pval_emp>TH | log2FC_archaic_vs_modern<0)/(1-TH/2), 
        pi1_pos.se=sd(replicate(1000,1-mean(sample(pval_emp>TH | log2FC_archaic_vs_modern<0,.N,replace=TRUE))/(1-TH/2))),
        pi1_neg=1-mean(pval_emp>TH | log2FC_archaic_vs_modern>0)/(1-TH/2), 
        pi1_neg.se=sd(replicate(1000,1-mean(sample(pval_emp>TH | log2FC_archaic_vs_modern>0,.N,replace=TRUE))/(1-TH/2)))),
        keyby=.(COND_ID=factor(COND_ID,names(color_setup_simplified_norep)), region, Source_introgression.byPOP)]

total=emVARs_obs_ctc[,.(.N, 
        nMarg_pos=sum(pval_emp<TH & log2FC_archaic_vs_modern>0),
        nMarg_neg=sum(pval_emp<TH & log2FC_archaic_vs_modern<0),
        nNS=sum(pval_emp>TH)),by=.(COND_ID=factor(COND_ID,names(color_setup_simplified_norep)))]

Figdata <- pi1_estim_signed[N>100 & Source_introgression.byPOP%in%c('Vindija','Denisova'),]
Figdata[,class:=factor(paste(region, Source_introgression.byPOP), levels=c('EUR Vindija','EAS Vindija','SAS Vindija','Papuan Vindija','Papuan Denisova','Agta Vindija','Agta Denisova'))]
Figdata <- melt(Figdata,measure.vars=list(c('pi1_pos','pi1_neg'),c('pi1_pos.se','pi1_neg.se')),value.name=c('pi1','pi1.se'))
levels(Figdata$variable)=c('pos','neg')

p <- ggplot(Figdata,aes(x=class,y=pi1, fill=COND_ID, alpha=Source_introgression.byPOP))
p <- p + geom_bar(stat='Identity',position='dodge') + ylab('Percentage of regulatory alleles (pi1 estimate)')
p <- p + geom_errorbar(aes(x=class, ymin=pi1-1.96*pi1.se,ymax=pi1+1.96*pi1.se), width=.5) +facet_grid(cols=vars(COND_ID),rows=vars(variable))
p <- p + scale_fill_manual(values=color_setup_simplified_norep)
p <- p + scale_alpha_manual(values=c('Vindija'=1,'Denisova'=0.7))+theme_plot(rotate.x=90,fontsize=10)
p <- p + guides(fill=FALSE) + xlab('')


pdf(sprintf("%s/08e_pi1_estim_emVar_per_condition_POP_x_source_signed_adapative.pdf", FIGURE_DIR), height = 5, width = 8)
print(p)
dev.off()

pi1_estim_signed_v2=merge(pi1_estim_signed,total,by=c('COND_ID'),suffix=c('','_total'))
pi1_estim_signed_v2[, pval_pos:=fisher.test(matrix(c(nMarg_neg_total+nNS_total-nMarg_neg-nNS,nMarg_pos_total-nMarg_pos, nMarg_neg+nNS,nMarg_pos),2,2))$p.value, by=.(COND_ID, region, Source_introgression.byPOP)]
pi1_estim_signed_v2[, pval_neg:=fisher.test(matrix(c(nMarg_pos_total+nNS_total-nMarg_pos-nNS,nMarg_neg_total-nMarg_neg, nMarg_pos+nNS,nMarg_neg),2,2))$p.value, by=.(COND_ID, region, Source_introgression.byPOP)]

pi1_estim_signed_v3 <- pi1_estim_signed_v2[order(pmin(pval_pos,pval_neg))][Source_introgression.byPOP%in%c('Denisova','Vindija') & N>100]
pi1_estim_signed_v3[,FDR_pos:=p.adjust(pval_pos,'fdr')]
pi1_estim_signed_v3[,FDR_neg:=p.adjust(pval_neg,'fdr')]
pi1_estim_signed_v3

################################################################################################################
########################### comparison of Pct of emVars per locus.  ############################################
################################################################################################################

dir.create(sprintf("%s/Pct_emVar_per_locus/", FIGURE_DIR), recursive = TRUE)
cat("\ntop 5 introgressed loci with highest Pct of emVars (any condition)\n")

minSNP=5
Pct_emVars_per_locus <- get_Pct_emVars(emVARs_obs_ctc, split_by = "introgression_locus")
Pct_emVars_per_locus <- make_it_long(Pct_emVars_per_locus, split_by = "introgression_locus")
Pct_emVars_per_locus <- merge(Pct_emVars_per_locus, locus_genes,by='introgression_locus')
Pct_emVars_per_locus[N_test > minSNP, FDR := p.adjust(Pvalue, "fdr"), by = measure]
Pct_emVars_per_locus[N_test > minSNP][order(Pvalue), head(.SD, 5), by = measure]
Pct_emVars_per_locus <- add_manhattan_coord(Pct_emVars_per_locus)
fwrite(Pct_emVars_per_locus[order(measure, Pvalue),],file=sprintf("%s/Pct_emVar_per_locus/01a_Pct_emVar_per_introgressed_locus.tsv", FIGURE_DIR),sep='\t')
print(Pct_emVars_per_locus[order(Pvalue),][Pvalue<0.01,])

for (i_MEASURE in seq_len(Measure_table[,.N])) {
  MEASURE <- Measure_table[i_MEASURE, measure]
  MEASURE_LABEL <- Measure_table[i_MEASURE, measure_label]
  MEASURE_NULL <- Measure_table[i_MEASURE, null]
  if (is.na(MEASURE_NULL)) {
    MEASURE_NULL <- Pct_emVars_per_locus[measure == MEASURE, median(Pct)]
  }
  pdf(sprintf("%s/Pct_emVar_per_locus/01b_Pct_%s_per_introgressed_locus.pdf", FIGURE_DIR, MEASURE), height = 3.5, width = 4.5)
  FigData <- Pct_emVars_per_locus[N_test > minSNP & measure == MEASURE,]
  p <- ggplot(FigData) + theme_plot(fontsize = 12) + theme(rect=element_blank(),line=element_line())
  #p <- p + geom_rect(aes(xmin = start_pos, xmax = end_pos, ymin = Pct_lo, ymax = Pct_hi,fill = FDR<0.05))
  #p <- p + geom_segment(aes(x = start_pos, xend = end_pos, y = Pct, yend = Pct), linewidth = 2, col='black')
  p <- p + geom_pointrange(aes(x = mid_pos, ymin = Pct_lo, ymax = Pct_hi, y=Pct, col=Pvalue<0.01, alpha=N_test))
  p <- p + geom_text_repel(data = FigData[Pvalue<0.01,], aes(x = mid_pos, y = Pct, label = locus_genes), size = 1, segment.size = 0.1)

  #p <- p + scale_x_continuous(label = '', breaks = sort(c(chrom_length$chrom_start,chrom_length$chrom_end)))
  p <- p + geom_vline(col='lightgrey', xintercept = chrom_length$chrom_end+5e5, linetype=3)
  p <- p + geom_hline(col='red', yintercept = MEASURE_NULL, linetype=2)
  p <- p + scale_x_continuous(label = chrom_length$chromosome , breaks = chrom_length$chrom_mid)
  p <- p + scale_color_manual(values=color_TRUEFALSE)
  p <- p + ylab(paste(MEASURE_LABEL, "per locus")) + xlab('Chromosome')
  print(p)
  dev.off()
}

cat("\ntop 5 introgressed loci with highest Pct of emVars, per cell line (any stimulus)\n")
Pct_emVars_per_locus_celltype <- get_Pct_emVars(emVARs_obs_ctc, split_by = c("introgression_locus", "celltype"), total_by = "celltype")
Pct_emVars_per_locus_celltype <- make_it_long(Pct_emVars_per_locus_celltype, split_by = c("introgression_locus", "celltype"))
Pct_emVars_per_locus_celltype <- merge(Pct_emVars_per_locus_celltype, locus_genes,by='introgression_locus')
Pct_emVars_per_locus_celltype[N_test > minSNP, FDR := p.adjust(Pvalue, "fdr"), by = measure]
Pct_emVars_per_locus_celltype[N_test > minSNP][order(Pvalue), head(.SD, 5), by = measure]
Pct_emVars_per_locus_celltype <- add_manhattan_coord(Pct_emVars_per_locus_celltype)
fwrite(Pct_emVars_per_locus_celltype[order(measure, celltype, Pvalue),],file=sprintf("%s/Pct_emVar_per_locus/01c_Pct_emVar_per_introgressed_locus__byCelltype.tsv", FIGURE_DIR),sep='\t')
print(Pct_emVars_per_locus_celltype[order(Pvalue),][Pvalue<0.01,])


for (i_MEASURE in seq_len(Measure_table[,.N])) {
  MEASURE <- Measure_table[i_MEASURE, measure]
  MEASURE_LABEL <- Measure_table[i_MEASURE, measure_label]
  MEASURE_NULL <- Measure_table[i_MEASURE, null]
  
  pdf(sprintf("%s/Pct_emVar_per_locus/01d_Pct_%s_per_introgressed_locus__byCelltype.pdf", FIGURE_DIR, MEASURE), height = 4, width = 6)
  FigData <- Pct_emVars_per_locus_celltype[measure == MEASURE & N_test > minSNP, ]
  p <- ggplot(FigData) 
  p <- p + theme_plot(fontsize = 12) + theme(rect=element_blank(),line=element_line())
  # p <- p + geom_rect(aes(xmin = start_pos, xmax = end_pos, ymin = Pct_lo, ymax = Pct_hi, fill = celltype))
  # p <- p + geom_segment(aes(x = start_pos, xend = end_pos, y = Pct, yend = Pct), linewidth = 2)
  p <- p + geom_pointrange(aes(x = mid_pos, ymin = Pct_lo, ymax = Pct_hi, y=Pct, col=celltype, alpha=N_test))
  p <- p + geom_pointrange(data=FigData[Pvalue<0.01,],aes(x = mid_pos, ymin = Pct_lo, ymax = Pct_hi, y=Pct, alpha=N_test),col='red')
  p <- p + geom_text_repel(data = FigData[Pvalue<0.01,], aes(x = mid_pos, y = Pct, label = locus_genes), size = 1, segment.size = 0.1)
  p <- p + geom_vline(col='lightgrey', xintercept = chrom_length$chrom_end+5e5, linetype=3)
  if (is.na(MEASURE_NULL)) {
    MEASURE_NULL <- Pct_emVars_per_locus_celltype[measure == MEASURE, .(Pct=median(Pct)), by = celltype]
    p <- p + geom_hline(data=MEASURE_NULL, aes(yintercept = Pct), linetype=2)
   }else{
    p <- p + geom_hline(yintercept=MEASURE_NULL, linetype=2)
   }
  p <- p + scale_x_continuous(label = chrom_length$chromosome , breaks = chrom_length$chrom_mid)
  p <- p + scale_color_manual(values = color_celline)
  p <- p + ylab(paste(MEASURE_LABEL, "per locus")) + xlab('Chromosome')
  p <- p + facet_grid(rows = vars(celltype))
  print(p)
  dev.off()
}

cat("\ntop 5 introgressed loci with highest Pct of emVars, per cell type (NS only)\n")
Pct_emVars_per_locus_celltypeNS <- get_Pct_emVars(emVARs_obs_ctc[condition == "NS", ], split_by = c("introgression_locus", "celltype"), total_by = "celltype")
Pct_emVars_per_locus_celltypeNS <- make_it_long(Pct_emVars_per_locus_celltypeNS, split_by = c("introgression_locus", "celltype"))
Pct_emVars_per_locus_celltypeNS <- merge(Pct_emVars_per_locus_celltypeNS, locus_genes,by='introgression_locus')
Pct_emVars_per_locus_celltypeNS[N_test > minSNP, FDR := p.adjust(Pvalue, "fdr"), by = measure]
Pct_emVars_per_locus_celltypeNS[N_test > minSNP][order(Pvalue), head(.SD, 5), by = measure]
Pct_emVars_per_locus_celltypeNS <- add_manhattan_coord(Pct_emVars_per_locus_celltypeNS)
fwrite(Pct_emVars_per_locus_celltypeNS[order(measure, celltype, Pvalue),],file=sprintf("%s/Pct_emVar_per_locus/01e_Pct_emVar_per_introgressed_locus__byCelltypeNS.tsv", FIGURE_DIR),sep='\t')
print(Pct_emVars_per_locus_celltypeNS[order(Pvalue),][Pvalue<0.01,])


for (i_MEASURE in seq_len(Measure_table[,.N])) {

  MEASURE <- Measure_table[i_MEASURE, measure]
  MEASURE_LABEL <- Measure_table[i_MEASURE, measure_label]
  MEASURE_NULL <- Measure_table[i_MEASURE, null]
  pdf(sprintf("%s/Pct_emVar_per_locus/01f_Pct_%s_per_introgressed_locus__byCelltypeNS.pdf", FIGURE_DIR, MEASURE), height = 4, width = 6)
  FigData <- Pct_emVars_per_locus_celltypeNS[measure == MEASURE & N_test > minSNP, ]
  p <- ggplot(FigData) 
  p <- p + theme_plot(fontsize = 12) + theme(rect=element_blank(),line=element_line())
  p <- p + geom_pointrange(data=FigData,aes(x = mid_pos, ymin = Pct_lo, ymax = Pct_hi, y=Pct, alpha=N_test, col=celltype))
  p <- p + geom_pointrange(data=FigData[Pvalue<0.01,],aes(x = mid_pos, ymin = Pct_lo, ymax = Pct_hi, y=Pct, alpha=N_test),col='red')
  p <- p + geom_text_repel(data = FigData[Pvalue<0.01,], aes(x = mid_pos, y = Pct, label = locus_genes), size = 1, segment.size = 0.1)
  # p <- p + geom_rect(aes(xmin = start_pos, xmax = end_pos, ymin = Pct_lo, ymax = Pct_hi, fill = celltype))
  # p <- p + geom_segment(aes(x = start_pos, xend = end_pos, y = Pct, yend = Pct), linewidth = 2)
  p <- p + geom_vline(col='lightgrey', xintercept = chrom_length$chrom_end+5e5, linetype=3)
  if (is.na(MEASURE_NULL)) {
    MEASURE_NULL <- Pct_emVars_per_locus_celltypeNS[measure == MEASURE, .(Pct=median(Pct)), by = celltype]
    p <- p + geom_hline(data=MEASURE_NULL, aes(yintercept = Pct), linetype=2)
   }else{
    p <- p + geom_hline(yintercept=MEASURE_NULL, linetype=2)
   }
  p <- p + scale_x_continuous(label = chrom_length$chromosome , breaks = chrom_length$chrom_mid)
  p <- p + scale_color_manual(values = color_celline)
  p <- p + ylab(paste(MEASURE_LABEL, "per locus")) + facet_grid(rows = vars(celltype))
  
  print(p)
  dev.off()
}

cat("\ntop 5 introgressed loci with highest Pct of emVars, per condition\n")
Pct_emVars_per_locus_condition <- get_Pct_emVars(emVARs_obs_ctc, split_by = c("introgression_locus", "COND_ID"), total_by = "COND_ID")
Pct_emVars_per_locus_condition <- make_it_long(Pct_emVars_per_locus_condition, split_by = c("introgression_locus", "COND_ID"))
Pct_emVars_per_locus_condition <- merge(Pct_emVars_per_locus_condition, locus_genes,by='introgression_locus')
Pct_emVars_per_locus_condition[N_test > 10, FDR := p.adjust(Pvalue, "fdr"), by = measure]
Pct_emVars_per_locus_condition[N_test > 10][order(Pvalue), head(.SD, 5), by = measure]
Pct_emVars_per_locus_condition <- add_manhattan_coord(Pct_emVars_per_locus_condition)
fwrite(Pct_emVars_per_locus_condition[order(measure,COND_ID,Pvalue),],file=sprintf("%s/Pct_emVar_per_locus/01g_Pct_emVar_per_introgressed_locus__byCondition.tsv", FIGURE_DIR),sep='\t')

for (i_MEASURE in seq_len(Measure_table[,.N])) {
  MEASURE <- Measure_table[i_MEASURE, measure]
  MEASURE_LABEL <- Measure_table[i_MEASURE, measure_label]
  MEASURE_NULL <- Measure_table[i_MEASURE, null]

  pdf(sprintf("%s/Pct_emVar_per_locus/01h_Pct_%s_per_introgressed_locus__byCondition.pdf", FIGURE_DIR, MEASURE), height = 9, width = 6)
  FigData <- Pct_emVars_per_locus_condition[measure == MEASURE & N_test > minSNP, ]
  p <- ggplot(FigData) 
  p <- p + theme_plot(fontsize = 12) + theme(rect=element_blank(),line=element_line())
  p <- p + geom_pointrange(data=FigData,aes(x = mid_pos, ymin = Pct_lo, ymax = Pct_hi, y=Pct, alpha=N_test, col=COND_ID))
  p <- p + geom_pointrange(data=FigData[Pvalue<0.01,],aes(x = mid_pos, ymin = Pct_lo, ymax = Pct_hi, y=Pct, alpha=N_test),col='red')
  p <- p + geom_text_repel(data = FigData[Pvalue<0.01,], aes(x = mid_pos, y = Pct, label = locus_genes), size = 1, segment.size = 0.1)
  # p <- p + geom_rect(aes(xmin = start_pos, xmax = end_pos, ymin = Pct_lo, ymax = Pct_hi, fill = celltype))
  # p <- p + geom_segment(aes(x = start_pos, xend = end_pos, y = Pct, yend = Pct), linewidth = 2)
  p <- p + geom_vline(col='lightgrey', xintercept = chrom_length$chrom_end+5e5, linetype=3)
  if (is.na(MEASURE_NULL)) {
    MEASURE_NULL <- Pct_emVars_per_locus_condition[measure == MEASURE, .(Pct=median(Pct)), by = COND_ID]
    p <- p + geom_hline(data=MEASURE_NULL, aes(yintercept = Pct), linetype=2)
   }else{
    p <- p + geom_hline(yintercept=MEASURE_NULL, linetype=2)
  }
  p <- p + scale_x_continuous(label = chrom_length$chromosome , breaks = chrom_length$chrom_mid)
  p <- p + scale_color_manual(values = color_setup_simplified_norep)
  p <- p + ylab(paste(MEASURE_LABEL, "per locus")) + xlab('Chromosome')
  p <- p + facet_grid(rows = vars(COND_ID))
  p <- p + theme(legend.box = "vertical")
  print(p)
  dev.off()
}

################################################################################################################
########################### comparison of Pct of emVars per locus (suggestive)  ################################
########################### Figure 5c & Table S5c ##############################################################
################################################################################################################

dir.create(sprintf("%s/Pct_emVar_per_locus/", FIGURE_DIR), recursive = TRUE)
cat("\ntop 5 introgressed loci with highest Pct of emVars (any condition)\n")

minSNP=5
Pct_emVars_per_locus <- get_Pct_emVars(emVARs_obs_ctc, split_by = "introgression_locus", emVarlist='suggestive')
Pct_emVars_per_locus <- make_it_long(Pct_emVars_per_locus, split_by = "introgression_locus")
Pct_emVars_per_locus <- merge(Pct_emVars_per_locus, locus_genes,by='introgression_locus')
Pct_emVars_per_locus[N_test > minSNP, FDR := p.adjust(Pvalue, "fdr"), by = measure]
Pct_emVars_per_locus[N_test > minSNP][order(Pvalue), head(.SD, 5), by = measure]
Pct_emVars_per_locus <- add_manhattan_coord(Pct_emVars_per_locus)
fwrite(Pct_emVars_per_locus[order(measure, Pvalue),],file=sprintf("%s/Pct_emVar_per_locus/02a_Pct_emVar_per_introgressed_locus_suggestive.tsv", FIGURE_DIR),sep='\t')
print(Pct_emVars_per_locus[order(Pvalue),][Pvalue<0.01,])

for (i_MEASURE in seq_len(Measure_table[,.N])) {
  MEASURE <- Measure_table[i_MEASURE, measure]
  MEASURE_LABEL <- Measure_table[i_MEASURE, measure_label]
  MEASURE_NULL <- Measure_table[i_MEASURE, null]
  if (is.na(MEASURE_NULL)) {
    MEASURE_NULL <- Pct_emVars_per_locus[measure == MEASURE, median(Pct)]
  }
  pdf(sprintf("%s/Pct_emVar_per_locus/02b_Pct_%s_per_introgressed_locus_suggestive.pdf", FIGURE_DIR, MEASURE), height = 3.5, width = 4.5)
  FigData <- Pct_emVars_per_locus[N_test > minSNP & measure == MEASURE,]
  p <- ggplot(FigData) + theme_plot(fontsize = 12) + theme(rect=element_blank(),line=element_line())
  #p <- p + geom_rect(aes(xmin = start_pos, xmax = end_pos, ymin = Pct_lo, ymax = Pct_hi,fill = FDR<0.05))
  #p <- p + geom_segment(aes(x = start_pos, xend = end_pos, y = Pct, yend = Pct), linewidth = 2, col='black')
  p <- p + geom_pointrange(aes(x = mid_pos, ymin = Pct_lo, ymax = Pct_hi, y=Pct, col=Pvalue<0.01, alpha=N_test))
  p <- p + geom_text_repel(data = FigData[Pvalue<0.01,], aes(x = mid_pos, y = Pct, label = locus_genes), size = 1, segment.size = 0.1)

  #p <- p + scale_x_continuous(label = '', breaks = sort(c(chrom_length$chrom_start,chrom_length$chrom_end)))
  p <- p + geom_vline(col='lightgrey', xintercept = chrom_length$chrom_end+5e5, linetype=3)
  p <- p + geom_hline(col='red', yintercept = MEASURE_NULL, linetype=2)
  p <- p + scale_x_continuous(label = chrom_length$chromosome , breaks = chrom_length$chrom_mid)
  p <- p + scale_color_manual(values=color_TRUEFALSE)
  p <- p + ylab(paste(MEASURE_LABEL, "per locus")) + xlab('Chromosome')
  print(p)
  dev.off()
}

cat("\ntop 5 introgressed loci with highest Pct of emVars, per cell line (any stimulus)\n")
Pct_emVars_per_locus_celltype <- get_Pct_emVars(emVARs_obs_ctc, split_by = c("introgression_locus", "celltype"), total_by = "celltype", emVarlist='suggestive')
Pct_emVars_per_locus_celltype <- make_it_long(Pct_emVars_per_locus_celltype, split_by = c("introgression_locus", "celltype"))
Pct_emVars_per_locus_celltype <- merge(Pct_emVars_per_locus_celltype, locus_genes,by='introgression_locus')
Pct_emVars_per_locus_celltype[N_test > minSNP, FDR := p.adjust(Pvalue, "fdr"), by = measure]
Pct_emVars_per_locus_celltype[N_test > minSNP][order(Pvalue), head(.SD, 5), by = measure]
Pct_emVars_per_locus_celltype <- add_manhattan_coord(Pct_emVars_per_locus_celltype)
fwrite(Pct_emVars_per_locus_celltype[order(measure, celltype, Pvalue),],file=sprintf("%s/Pct_emVar_per_locus/02c_Pct_emVar_per_introgressed_locus__byCelltype_suggestive.tsv", FIGURE_DIR),sep='\t')
print(Pct_emVars_per_locus_celltype[order(Pvalue),][Pvalue<0.01,])


for (i_MEASURE in seq_len(Measure_table[,.N])) {
  MEASURE <- Measure_table[i_MEASURE, measure]
  MEASURE_LABEL <- Measure_table[i_MEASURE, measure_label]
  MEASURE_NULL <- Measure_table[i_MEASURE, null]
  
  pdf(sprintf("%s/Pct_emVar_per_locus/02d_Pct_%s_per_introgressed_locus__byCelltype_suggestive.pdf", FIGURE_DIR, MEASURE), height = 4, width = 6)
  FigData <- Pct_emVars_per_locus_celltype[measure == MEASURE & N_test > minSNP, ]
  p <- ggplot(FigData) 
  p <- p + theme_plot(fontsize = 12) + theme(rect=element_blank(),line=element_line())
  # p <- p + geom_rect(aes(xmin = start_pos, xmax = end_pos, ymin = Pct_lo, ymax = Pct_hi, fill = celltype))
  # p <- p + geom_segment(aes(x = start_pos, xend = end_pos, y = Pct, yend = Pct), linewidth = 2)
  p <- p + geom_pointrange(aes(x = mid_pos, ymin = Pct_lo, ymax = Pct_hi, y=Pct, col=celltype, alpha=N_test))
  p <- p + geom_pointrange(data=FigData[Pvalue<0.01,],aes(x = mid_pos, ymin = Pct_lo, ymax = Pct_hi, y=Pct, alpha=N_test),col='red')
  p <- p + geom_text_repel(data = FigData[Pvalue<0.01,], aes(x = mid_pos, y = Pct, label = locus_genes), size = 1, segment.size = 0.1)
  p <- p + geom_vline(col='lightgrey', xintercept = chrom_length$chrom_end+5e5, linetype=3)
  if (is.na(MEASURE_NULL)) {
    MEASURE_NULL <- Pct_emVars_per_locus_celltype[measure == MEASURE, .(Pct=median(Pct)), by = celltype]
    p <- p + geom_hline(data=MEASURE_NULL, aes(yintercept = Pct), linetype=2)
   }else{
    p <- p + geom_hline(yintercept=MEASURE_NULL, linetype=2)
   }
  p <- p + scale_x_continuous(label = chrom_length$chromosome , breaks = chrom_length$chrom_mid)
  p <- p + scale_color_manual(values = color_celline)
  p <- p + ylab(paste(MEASURE_LABEL, "per locus")) + xlab('Chromosome')
  p <- p + facet_grid(rows = vars(celltype))
  print(p)
  dev.off()
}


cat("\ntop 5 introgressed loci with highest Pct of emVars, per cell type (NS only)\n")
Pct_emVars_per_locus_celltypeNS <- get_Pct_emVars(emVARs_obs_ctc[condition == "NS", ], split_by = c("introgression_locus", "celltype"), total_by = "celltype", emVarlist='suggestive')
Pct_emVars_per_locus_celltypeNS <- make_it_long(Pct_emVars_per_locus_celltypeNS, split_by = c("introgression_locus", "celltype"))
Pct_emVars_per_locus_celltypeNS <- merge(Pct_emVars_per_locus_celltypeNS, locus_genes,by='introgression_locus')
Pct_emVars_per_locus_celltypeNS[N_test > minSNP, FDR := p.adjust(Pvalue, "fdr"), by = measure]
Pct_emVars_per_locus_celltypeNS[N_test > minSNP][order(Pvalue), head(.SD, 5), by = measure]
Pct_emVars_per_locus_celltypeNS <- add_manhattan_coord(Pct_emVars_per_locus_celltypeNS)
fwrite(Pct_emVars_per_locus_celltypeNS[order(measure, celltype, Pvalue),],file=sprintf("%s/Pct_emVar_per_locus/02e_Pct_emVar_per_introgressed_locus__byCelltypeNS_suggestive.tsv", FIGURE_DIR),sep='\t')
print(Pct_emVars_per_locus_celltypeNS[order(Pvalue),][Pvalue<0.01,])


for (i_MEASURE in seq_len(Measure_table[,.N])) {

  MEASURE <- Measure_table[i_MEASURE, measure]
  MEASURE_LABEL <- Measure_table[i_MEASURE, measure_label]
  MEASURE_NULL <- Measure_table[i_MEASURE, null]
  pdf(sprintf("%s/Pct_emVar_per_locus/02f_Pct_%s_per_introgressed_locus__byCelltypeNS_suggestive.pdf", FIGURE_DIR, MEASURE), height = 4, width = 6)
  FigData <- Pct_emVars_per_locus_celltypeNS[measure == MEASURE & N_test > minSNP, ]
  p <- ggplot(FigData) 
  p <- p + theme_plot(fontsize = 12) + theme(rect=element_blank(),line=element_line())
  p <- p + geom_pointrange(data=FigData,aes(x = mid_pos, ymin = Pct_lo, ymax = Pct_hi, y=Pct, alpha=N_test, col=celltype))
  p <- p + geom_pointrange(data=FigData[Pvalue<0.01,],aes(x = mid_pos, ymin = Pct_lo, ymax = Pct_hi, y=Pct, alpha=N_test),col='red')
  p <- p + geom_text_repel(data = FigData[Pvalue<0.01,], aes(x = mid_pos, y = Pct, label = locus_genes), size = 1, segment.size = 0.1)
  # p <- p + geom_rect(aes(xmin = start_pos, xmax = end_pos, ymin = Pct_lo, ymax = Pct_hi, fill = celltype))
  # p <- p + geom_segment(aes(x = start_pos, xend = end_pos, y = Pct, yend = Pct), linewidth = 2)
  p <- p + geom_vline(col='lightgrey', xintercept = chrom_length$chrom_end+5e5, linetype=3)
  if (is.na(MEASURE_NULL)) {
    MEASURE_NULL <- Pct_emVars_per_locus_celltypeNS[measure == MEASURE, .(Pct=median(Pct)), by = celltype]
    p <- p + geom_hline(data=MEASURE_NULL, aes(yintercept = Pct), linetype=2)
   }else{
    p <- p + geom_hline(yintercept=MEASURE_NULL, linetype=2)
   }
  p <- p + scale_x_continuous(label = chrom_length$chromosome , breaks = chrom_length$chrom_mid)
  p <- p + scale_color_manual(values = color_celline)
  p <- p + ylab(paste(MEASURE_LABEL, "per locus")) + facet_grid(rows = vars(celltype))
  
  print(p)
  dev.off()
}

########################### Figure 5c 
cat("\ntop 5 introgressed loci with highest Pct of emVars, per condition\n")
# Pct_emVars_per_locus_condition <- emVARs_obs_ctc[, .(N_introgressed_allele = length(unique(ID)), N_activity_altering_allele = length(unique(ID[is_emVar_CRITERIA]))), by = .(introgression_locus, COND_ID)]
# Pct_emVars_per_locus_condition[, N_activity_altering_allele_tot := sum(N_activity_altering_allele), by = COND_ID]
# Pct_emVars_per_locus_condition[, N_tested_tot := sum(N_introgressed_allele), by = COND_ID]
# Pct_emVars_per_locus_condition[, P := binom.test(N_activity_altering_allele, N_introgressed_allele, p = N_activity_altering_allele_tot / N_tested_tot, alternative = "greater")$p.value, by = .(introgression_locus, COND_ID)]
# Pct_emVars_per_locus_condition[, Pct := binom.test(N_activity_altering_allele, N_introgressed_allele, p = N_activity_altering_allele_tot / N_tested_tot, alternative = "greater")$estimate, by = .(introgression_locus, COND_ID)]
# Pct_emVars_per_locus_condition[, Pct_lo := binom.test(N_activity_altering_allele, N_introgressed_allele, p = N_activity_altering_allele_tot / N_tested_tot, alternative = "greater")$conf.int[1], by = .(introgression_locus, COND_ID)]
# Pct_emVars_per_locus_condition[, Pct_hi := binom.test(N_activity_altering_allele, N_introgressed_allele, p = N_activity_altering_allele_tot / N_tested_tot, alternative = "greater")$conf.int[2], by = .(introgression_locus, COND_ID)]
# Pct_emVars_per_locus_condition[N_introgressed_allele > 10, FDR := p.adjust(P, "fdr")]
# Pct_emVars_per_locus_condition[N_introgressed_allele > 10][order(P), head(.SD, 3), by = COND_ID]

Pct_emVars_per_locus_condition <- get_Pct_emVars(emVARs_obs_ctc, split_by = c("introgression_locus", "COND_ID"), total_by = "COND_ID", emVarlist='suggestive')
Pct_emVars_per_locus_condition <- make_it_long(Pct_emVars_per_locus_condition, split_by = c("introgression_locus", "COND_ID"))
Pct_emVars_per_locus_condition <- merge(Pct_emVars_per_locus_condition, locus_genes,by='introgression_locus')
Pct_emVars_per_locus_condition[N_test > 10, FDR := p.adjust(Pvalue, "fdr"), by = measure]
Pct_emVars_per_locus_condition[N_test > 10][order(Pvalue), head(.SD, 5), by = measure]
Pct_emVars_per_locus_condition <- add_manhattan_coord(Pct_emVars_per_locus_condition)
fwrite(Pct_emVars_per_locus_condition[order(measure,COND_ID,Pvalue),],file=sprintf("%s/Pct_emVar_per_locus/02g_Pct_emVar_per_introgressed_locus__byCondition_suggestive.tsv", FIGURE_DIR),sep='\t')

for (i_MEASURE in seq_len(Measure_table[,.N])) {
  MEASURE <- Measure_table[i_MEASURE, measure]
  MEASURE_LABEL <- Measure_table[i_MEASURE, measure_label]
  MEASURE_NULL <- Measure_table[i_MEASURE, null]

  pdf(sprintf("%s/Pct_emVar_per_locus/02h_Pct_%s_per_introgressed_locus__byCondition_suggestive.pdf", FIGURE_DIR, MEASURE), height = 9, width = 6)
  FigData <- Pct_emVars_per_locus_condition[measure == MEASURE & N_test > minSNP, ]
  p <- ggplot(FigData) 
  p <- p + theme_plot(fontsize = 12) + theme(rect=element_blank(),line=element_line())
  p <- p + geom_pointrange(data=FigData,aes(x = mid_pos, ymin = Pct_lo, ymax = Pct_hi, y=Pct, alpha=N_test, col=COND_ID))
  p <- p + geom_pointrange(data=FigData[Pvalue<0.01,],aes(x = mid_pos, ymin = Pct_lo, ymax = Pct_hi, y=Pct, alpha=N_test),col='red')
  p <- p + geom_text_repel(data = FigData[Pvalue<0.01,], aes(x = mid_pos, y = Pct, label = locus_genes), size = 1, segment.size = 0.1)
  # p <- p + geom_rect(aes(xmin = start_pos, xmax = end_pos, ymin = Pct_lo, ymax = Pct_hi, fill = celltype))
  # p <- p + geom_segment(aes(x = start_pos, xend = end_pos, y = Pct, yend = Pct), linewidth = 2)
  p <- p + geom_vline(col='lightgrey', xintercept = chrom_length$chrom_end+5e5, linetype=3)
  if (is.na(MEASURE_NULL)) {
    MEASURE_NULL <- Pct_emVars_per_locus_condition[measure == MEASURE, .(Pct=median(Pct)), by = COND_ID]
    p <- p + geom_hline(data=MEASURE_NULL, aes(yintercept = Pct), linetype=2)
   }else{
    p <- p + geom_hline(yintercept=MEASURE_NULL, linetype=2)
  }
  p <- p + scale_x_continuous(label = chrom_length$chromosome , breaks = chrom_length$chrom_mid)
  p <- p + scale_color_manual(values = color_setup_simplified_norep)
  p <- p + ylab(paste(MEASURE_LABEL, "per locus")) + xlab('Chromosome')
  p <- p + facet_grid(rows = vars(COND_ID))
  p <- p + theme(legend.box = "vertical")
  print(p)
  dev.off()
}

# FigData <- Pct_emVars_per_locus_condition[grepl('Expression',measure)]
#  p <- ggplot(FigData[FDR>0.01],aes(x= mid_pos,y=-log10(FDR)))+geom_point(col='lightgrey',size=.2)
#  p <- p + geom_point(data=FigData[FDR<=0.01,], aes(x= mid_pos,y=-log10(FDR),fill=COND_ID,shape=measure))
#  p <- p + scale_shape_manual(values=c('ExpressionIncreasing'=24,'ExpressionDecreasing'=25))
#  p <- p + scale_fill_manual(values=color_setup_simplified_norep)
#  p <- p + scale_color_manual(values=color_setup_simplified_norep)+theme_plot(fontsize=11)
#  fig_annotations <- FigData[FDR<=0.01,][order(Pvalue)][!duplicated(locus_genes)]
#  p <- p + geom_text_repel(data = fig_annotations, aes(x = mid_pos, y = -log10(FDR), label = locus_genes, col=COND_ID), size = 2, segment.size = 0.1)
#    # p <- p + geom_rect(aes(xmin = start_pos, xmax = end_pos, ymin = Pct_lo, ymax = Pct_hi, fill = celltype))
#   # p <- p + geom_segment(aes(x = start_pos, xend = end_pos, y = Pct, yend = Pct), linewidth = 2)
#   p <- p + geom_vline(col='lightgrey', xintercept = chrom_length$chrom_end+5e5, linetype=3)
#   p <- p + scale_x_continuous(label = chrom_length$chromosome , breaks = chrom_length$chrom_mid)
#   p <- p + ylab('emVar enrichment at the locus [log10 (FDR)]') + xlab('Chromosome')
#  # p <- p + facet_grid(rows = vars(COND_ID))
#   p <- p + theme(legend.box = "vertical")

#   pdf(sprintf("%s/Pct_emVar_per_locus/03a_EmvarEnrichent_of_introgressed_loci_FDR1__byCondition_suggestive.pdf", FIGURE_DIR), height = 4, width = 9)
#     print(p)
#   dev.off()




# FigData <- Pct_emVars_per_locus_condition[grepl('Expression',measure)]
# FigData <-merge(FigData, condition_summary[,.(COND_ID,celline)], by='COND_ID')
#  p <- ggplot(FigData[FDR>0.05],aes(x= mid_pos,y=-log10(FDR)))+geom_point(col='lightgrey',size=.2)
#  p <- p + geom_point(data=FigData[FDR<=0.05,], aes(x= mid_pos,y=-log10(FDR),fill=COND_ID,shape=measure))
#  p <- p + scale_shape_manual(values=c('ExpressionIncreasing'=24,'ExpressionDecreasing'=25))
#  p <- p + scale_fill_manual(values=color_setup_simplified_norep)
#  p <- p + scale_color_manual(values=color_setup_simplified_norep)+theme_plot(fontsize=11)
#  fig_annotations <- FigData[FDR<=0.05,][order(Pvalue)][!duplicated(locus_genes)]
#  p <- p + geom_text_repel(data = fig_annotations, aes(x = mid_pos, y = -log10(FDR), label = locus_genes, col=COND_ID), size = 2, segment.size = 0.1)
#    # p <- p + geom_rect(aes(xmin = start_pos, xmax = end_pos, ymin = Pct_lo, ymax = Pct_hi, fill = celltype))
#   # p <- p + geom_segment(aes(x = start_pos, xend = end_pos, y = Pct, yend = Pct), linewidth = 2)
#   p <- p + geom_vline(col='lightgrey', xintercept = chrom_length$chrom_end+5e5, linetype=3)
#   p <- p + scale_x_continuous(label = chrom_length$chromosome , breaks = chrom_length$chrom_mid)
#   p <- p + ylab('emVar enrichment at the locus [log10 (FDR)]') + xlab('Chromosome')
#   p <- p + facet_grid(rows = vars(celline))
#   p <- p + theme(legend.box = "vertical")

#   pdf(sprintf("%s/Pct_emVar_per_locus/03a_EmvarEnrichent_of_introgressed_loci_FDR5__byCondition_suggestive.pdf", FIGURE_DIR), height = 6, width = 9)
#     print(p)
#   dev.off()

minSNP=5
TableS5c <- Pct_emVars_per_locus_condition[grepl('Expression',measure),]
total_condition <- get_Pct_emVars(emVARs_obs_ctc, split_by = c("COND_ID"), total_by = "COND_ID", emVarlist='suggestive')
total_condition <- make_it_long(total_condition, split_by = c("COND_ID"))
total_condition <- total_condition[grepl('Expression',measure),]
TableS5c <- merge(TableS5c, total_condition[,.(N_obs_tot=N_obs,N_test_tot=N_test, COND_ID,measure)], by=c('measure','COND_ID'), allow.cartesian=TRUE)
#TableS5c[,pct_global:=N_obs_tot/N_test_tot]
#TableS5c[,pi1_global:=pmax(0,1-1-(1-pct_global)/0.975)]
#TableS5c[,pi1_local:=pmax(0,1-(1-N_obs/N_test)/0.975)]
#TableS5c[,FE:=Pct/pct_global]
TableS5c[,c('OR','Fisher_P'):=fisher.test(matrix(c(N_obs,N_obs_tot-N_obs,N_test-N_obs,N_test_tot-N_obs_tot-N_test+N_obs),2,2), alternative = "greater")[c('estimate','p.value')],by=.(measure,COND_ID,introgression_locus)]
#TableS5c[,Binom_P:=binom.test(N_obs,N_test,p=pct_global, alternative = "greater")$p.value,by=.(measure,COND_ID,introgression_locus)]
TableS5c[N_test >= minSNP,FDR:=p.adjust(Fisher_P,"fdr"),by=.(measure,COND_ID)]
# CCL17 & CCL22:10.1111/j.1398-9995.2009.02095.x

TableS5c_print <- TableS5c[N_test >= minSNP,.(introgression_locus,COND_ID,
      Effect_on_expression=gsub('Expression','',measure),
      N_regtulatory_allele=N_obs,
      N_tested_allele=N_test,
      Pct_regulatory=Pct*100,
      Pct_regulatory_95ci=sprintf('[%s-%s]',round(Pct_lo*100,1),round(Pct_hi*100,1)),
      Pct_regulatory_genome_wide=N_obs_tot/N_test_tot,
      OR, Fisher_P, FDR, locus_genes)]
fwrite(TableS5c_print,file=sprintf("%s/SupTables/SupTable5c_Pct_emVar_per_introgressed_locus__byCondition_suggestive.tsv", FIGURE_DIR),sep='\t')
# emVARs_obs_ctc_gene <- merge(emVARs_obs_ctc, unique(NearbyGenes[distance_to_gene<1e4,.(crsID,gene_name)]), by='crsID', allow.cartesian=TRUE)

# emVARs_obs_ctc_gene <- merge(emVARs_obs_ctc, unique(NearbyGenes[distance_to_gene<1e4,.(crsID,gene_name,gene_type)]), by='crsID', allow.cartesian=TRUE)

# minSNP=5
# Pct_emVars_per_gene_condition <- get_Pct_emVars(emVARs_obs_ctc_gene, split_by = c("gene_name",'gene_type',"introgression_locus","COND_ID"), total_by = "COND_ID", emVarlist='suggestive')
# Pct_emVars_per_gene_condition <- make_it_long(Pct_emVars_per_gene_condition, split_by = c("gene_name",'gene_type', "introgression_locus","COND_ID"))
# Pct_emVars_per_gene_condition[N_test >= minSNP, FDR := p.adjust(Pvalue, "fdr"), by = measure]
# Pct_emVars_per_gene_condition[N_test >= minSNP][order(Pvalue), head(.SD, 5), by = measure]
# Pct_emVars_per_gene_condition <- add_manhattan_coord(Pct_emVars_per_gene_condition)
# fwrite(Pct_emVars_per_gene_condition[order(measure,COND_ID,Pvalue),],file=sprintf("%s/Pct_emVar_per_locus/04a_Pct_emVar_per_introgressed_gene100kb__byCondition_suggestive.tsv", FIGURE_DIR),sep='\t')


# FigData <- Pct_emVars_per_locus_condition[grepl('Expression',measure)]
# FigData <-merge(FigData, condition_summary[,.(COND_ID,celline)], by='COND_ID')
#  p <- ggplot(FigData[FDR>0.05],aes(x= mid_pos,y=-log10(FDR)))+geom_point(col='lightgrey',size=.2)
#  p <- p + geom_point(data=FigData[FDR<=0.05,], aes(x= mid_pos,y=-log10(FDR),fill=COND_ID,shape=measure))
#  p <- p + scale_shape_manual(values=c('ExpressionIncreasing'=24,'ExpressionDecreasing'=25))
#  p <- p + scale_fill_manual(values=color_setup_simplified_norep)
#  p <- p + scale_color_manual(values=color_setup_simplified_norep)+theme_plot(fontsize=11)
#  fig_annotations <- Pct_emVars_per_gene_condition[Pvalue<.01  & paste(COND_ID,introgression_locus, measure) %in% FigData[FDR<.05,paste(COND_ID,introgression_locus, measure)]  ,.(introgression_locus,  COND_ID, measure, gene_name, Pvalue)]
#  fig_annotations <- fig_annotations[order(introgression_locus,COND_ID,Pvalue), .(locus_genes=paste(unique(gene_name),collapse='/')), by=.(introgression_locus,COND_ID,measure)]
#  fig_annotations=merge(fig_annotations,FigData[FDR<=0.05,.(mid_pos,FDR,introgression_locus,COND_ID,celline,measure)],by=c('introgression_locus','COND_ID','measure'),allow.cartesian=TRUE)
#   p <- p + geom_text_repel(data = fig_annotations, aes(x = mid_pos, y = -log10(FDR), label = locus_genes, col=COND_ID), size = 2, segment.size = 0.1)
#    # p <- p + geom_rect(aes(xmin = start_pos, xmax = end_pos, ymin = Pct_lo, ymax = Pct_hi, fill = celltype))
#   # p <- p + geom_segment(aes(x = start_pos, xend = end_pos, y = Pct, yend = Pct), linewidth = 2)
#   p <- p + geom_vline(col='lightgrey', xintercept = chrom_length$chrom_end+5e5, linetype=3)
#   p <- p + scale_x_continuous(label = chrom_length$chromosome , breaks = chrom_length$chrom_mid)
#   p <- p + ylab('emVar enrichment at the locus [log10 (FDR)]') + xlab('Chromosome')
#   p <- p + facet_grid(rows = vars(celline))
#   p <- p + theme(legend.box = "vertical")

#   pdf(sprintf("%s/Pct_emVar_per_locus/03a_EmvarEnrichent_of_introgressed_loci_FDR5_gene_annot__byCondition_suggestive.pdf", FIGURE_DIR), height = 6, width = 9)
#     print(p)
#   dev.off()

# FigData <- Pct_emVars_per_gene_condition[grepl('Expression',measure)]
#  p <- ggplot(FigData[FDR>0.05],aes(x= mid_pos,y=-log10(FDR)))+geom_point(col='lightgrey',size=.2)
#  p <- p + geom_point(data=FigData[FDR<=0.05,], aes(x= mid_pos,y=-log10(FDR),fill=COND_ID,shape=measure))
#  p <- p + scale_shape_manual(values=c('ExpressionIncreasing'=24,'ExpressionDecreasing'=25))
#  p <- p + scale_fill_manual(values=color_setup_simplified_norep)
#  p <- p + scale_color_manual(values=color_setup_simplified_norep)+theme_plot(fontsize=11)
#  fig_annotations <- FigData[FDR<=0.05,][order(Pvalue)][!duplicated(gene_name)]
#  p <- p + geom_text_repel(data = fig_annotations, aes(x = mid_pos, y = -log10(FDR), label = gene_name, col=COND_ID), size = 2, segment.size = 0.1)
#    # p <- p + geom_rect(aes(xmin = start_pos, xmax = end_pos, ymin = Pct_lo, ymax = Pct_hi, fill = celltype))
#   # p <- p + geom_segment(aes(x = start_pos, xend = end_pos, y = Pct, yend = Pct), linewidth = 2)
#   p <- p + geom_vline(col='lightgrey', xintercept = chrom_length$chrom_end+5e5, linetype=3)
#   p <- p + scale_x_continuous(label = chrom_length$chromosome , breaks = chrom_length$chrom_mid)
#   p <- p + ylab('emVar enrichment at the locus\n [log10 (FDR)]') + xlab('Chromosome')
#  # p <- p + facet_grid(rows = vars(COND_ID))
#   p <- p + theme(legend.box = "vertical")

#   pdf(sprintf("%s/Pct_emVar_per_locus/03a_EmvarEnrichent_of_gene10kb__byCondition_suggestive.pdf", FIGURE_DIR), height = 4, width = 9)
#     print(p)
#   dev.off()

maxDist=1E4
locus_th=0.01
gene_th=0.05
minSNP=5
minSNP_gene=3
IMMUNE_ONLY=TRUE

emVARs_obs_ctc_gene <- merge(emVARs_obs_ctc, unique(NearbyGenes[distance_to_gene<maxDist,.(crsID,gene_name,gene_type)]), by='crsID', allow.cartesian=TRUE)


Pct_emVars_per_gene_condition <- get_Pct_emVars(emVARs_obs_ctc_gene, split_by = c("gene_name",'gene_type',"introgression_locus","COND_ID"), total_by = "COND_ID", emVarlist='suggestive')
Pct_emVars_per_gene_condition <- make_it_long(Pct_emVars_per_gene_condition, split_by = c("gene_name", 'gene_type',"introgression_locus","COND_ID"))
Pct_emVars_per_gene_condition[N_test >= minSNP_gene, FDR := p.adjust(Pvalue, "fdr"), by = measure]
Pct_emVars_per_gene_condition[N_test >= minSNP_gene][order(Pvalue), head(.SD, 5), by = measure]
Pct_emVars_per_gene_condition <- add_manhattan_coord(Pct_emVars_per_gene_condition)
fwrite(Pct_emVars_per_gene_condition[order(measure,COND_ID,Pvalue),],file=sprintf("%s/Pct_emVar_per_locus/04a_Pct_emVar_per_introgressed_gene100kb__byCondition_suggestive.tsv", FIGURE_DIR),sep='\t')

TableS5d <- merge(Pct_emVars_per_gene_condition, total_condition[,.(N_obs_tot=N_obs,N_test_tot=N_test, COND_ID,measure)], by=c('measure','COND_ID'), allow.cartesian=TRUE)
TableS5d[,c('OR','Fisher_P'):=fisher.test(matrix(c(N_obs,N_obs_tot-N_obs,N_test-N_obs,N_test_tot-N_obs_tot-N_test+N_obs),2,2), alternative = "greater")[c('estimate','p.value')],by=.(measure,COND_ID,introgression_locus,gene_name)]
TableS5d_print <- TableS5d[N_test >= minSNP_gene & paste(COND_ID,introgression_locus, measure) %in% TableS5c[Fisher_P<locus_th, paste(COND_ID,introgression_locus, measure)],]
TableS5d_print[N_test >= minSNP_gene,FDR:=p.adjust(Fisher_P,"fdr"),by=.(measure,COND_ID)]
TableS5d_print[,Dist_kb_from_gene:=maxDist/1000]

TableS5d_print <- TableS5d_print[N_test >= minSNP_gene,.(introgression_locus,COND_ID,gene_name,
      Effect_on_expression=gsub('Expression','',measure),
      N_regtulatory_allele=N_obs,
      N_tested_allele=N_test,
      Pct_regulatory=Pct*100,
      Pct_regulatory_95ci=sprintf('[%s-%s]',round(Pct_lo*100,1),round(Pct_hi*100,1)),
      Pct_regulatory_genome_wide=N_obs_tot/N_test_tot,
      OR, Fisher_P, FDR)]
  fwrite(TableS5d_print,file=sprintf("%s/SupTables/SupTable5d_Pct_emVar_per_gene_10kb__byCondition_suggestive.tsv", FIGURE_DIR),sep='\t')
		
#FigData <- merge(TableS5c, TableS5d, by=c('introgression_locus','measure','COND_ID'),suffix=c('.locus','.gene'))
#FigData <- TableS5d

# FigData <- TableS5d
# #FigData <- Pct_emVars_per_gene_condition[grepl('Expression',measure)]
#  p <- ggplot(FigData[Fisher_P>0.05],aes(x= mid_pos,y=-log10(Fisher_P)))+geom_point(col='lightgrey',size=.2)
#  p <- p + geom_point(data=FigData[Fisher_P<=0.05,], aes(x= mid_pos,y=-log10(FDR),fill=COND_ID,shape=measure))
#  p <- p + scale_shape_manual(values=c('ExpressionIncreasing'=24,'ExpressionDecreasing'=25))
#  p <- p + scale_fill_manual(values=color_setup_simplified_norep)
#  p <- p + scale_color_manual(values=color_setup_simplified_norep)+theme_plot(fontsize=11)
#  fig_annotations <- FigData[Fisher_P<=0.05 & N_test >= minSNP,][order(Pvalue)][!duplicated(gene_name)]
#  p <- p + geom_text_repel(data = fig_annotations, aes(x = mid_pos, y = -log10(Fisher_P), label = gene_name, col=COND_ID), size = 2, segment.size = 0.1)
#    # p <- p + geom_rect(aes(xmin = start_pos, xmax = end_pos, ymin = Pct_lo, ymax = Pct_hi, fill = celltype))
#   # p <- p + geom_segment(aes(x = start_pos, xend = end_pos, y = Pct, yend = Pct), linewidth = 2)
#   p <- p + geom_vline(col='lightgrey', xintercept = chrom_length$chrom_end+5e5, linetype=3)
#   p <- p + scale_x_continuous(label = chrom_length$chromosome , breaks = chrom_length$chrom_mid)
#   p <- p + ylab('emVar enrichment at the locus\n [log10 (FDR)]') + xlab('Chromosome')
#  # p <- p + facet_grid(rows = vars(COND_ID))
#   p <- p + theme(legend.box = "vertical")

#   pdf(sprintf("%s/Pct_emVar_per_locus/03b_EmvarEnrichent_of_gene100kb__byCondition_suggestive_v2.pdf", FIGURE_DIR), height = 4, width = 9)
#     print(p)
#   dev.off()



# emVARs_obs_ctc_gene <- merge(emVARs_obs_ctc, unique(NearbyGenes[distance_to_gene<1e6,.(crsID,gene_name,gene_start,gene_end,chrom)]), by='crsID', allow.cartesian=TRUE)

# minSNP=5
# Pct_emVars_per_gene_condition <- get_Pct_emVars(emVARs_obs_ctc_gene, split_by = c("gene_name","introgression_locus","COND_ID"), total_by = "COND_ID", emVarlist='suggestive')
# Pct_emVars_per_gene_condition <- make_it_long(Pct_emVars_per_gene_condition, split_by = c("gene_name", "introgression_locus","COND_ID"))
# Pct_emVars_per_gene_condition[N_test >= minSNP, FDR := p.adjust(Pvalue, "fdr"), by = measure]
# Pct_emVars_per_gene_condition[N_test >= minSNP][order(Pvalue), head(.SD, 5), by = measure]
# Pct_emVars_per_gene_condition <- add_manhattan_coord(Pct_emVars_per_gene_condition)
# fwrite(Pct_emVars_per_gene_condition[order(measure,COND_ID,Pvalue),],file=sprintf("%s/Pct_emVar_per_locus/04a_Pct_emVar_per_introgressed_gene1Mb__byCondition_suggestive.tsv", FIGURE_DIR),sep='\t')


# FigData <- Pct_emVars_per_gene_condition[grepl('Expression',measure)]
#  p <- ggplot(FigData[FDR>0.05],aes(x= mid_pos,y=-log10(FDR)))+geom_point(col='lightgrey',size=.2)
#  p <- p + geom_point(data=FigData[FDR<=0.05,], aes(x= mid_pos,y=-log10(FDR),fill=COND_ID,shape=measure))
#  p <- p + scale_shape_manual(values=c('ExpressionIncreasing'=24,'ExpressionDecreasing'=25))
#  p <- p + scale_fill_manual(values=color_setup_simplified_norep)
#  p <- p + scale_color_manual(values=color_setup_simplified_norep)+theme_plot(fontsize=11)
#  fig_annotations <- FigData[FDR<=0.05,][order(Pvalue)][!duplicated(gene_name)]
#  p <- p + geom_text_repel(data = fig_annotations, aes(x = mid_pos, y = -log10(FDR), label = gene_name, col=COND_ID), size = 2, segment.size = 0.1)
#    # p <- p + geom_rect(aes(xmin = start_pos, xmax = end_pos, ymin = Pct_lo, ymax = Pct_hi, fill = celltype))
#   # p <- p + geom_segment(aes(x = start_pos, xend = end_pos, y = Pct, yend = Pct), linewidth = 2)
#   p <- p + geom_vline(col='lightgrey', xintercept = chrom_length$chrom_end+5e5, linetype=3)
#   p <- p + scale_x_continuous(label = chrom_length$chromosome , breaks = chrom_length$chrom_mid)
#   p <- p + ylab('emVar enrichment at the locus\n [log10 (FDR)]') + xlab('Chromosome')
#  # p <- p + facet_grid(rows = vars(COND_ID))
#   p <- p + theme(legend.box = "vertical")

#   pdf(sprintf("%s/Pct_emVar_per_locus/03c_EmvarEnrichent_of_gene1Mb__byCondition_suggestive.pdf", FIGURE_DIR), height = 4, width = 9)
#     print(p)
#   dev.off()
# weights vs dist -> setNames(round(exp(-dist_kb/20)*100),dist_kb)
# weightingFun=function(dist){round(exp(-dist_kb/20)*100)}

  # TableS5d[Fisher_P<=0.01 & N_test >= minSNP & paste(COND_ID,introgression_locus, measure) %in% FigData[Fisher_P<.01,paste(COND_ID,introgression_locus, measure)],]

  FigData <- TableS5c
  FigData <-merge(FigData, condition_summary[,.(COND_ID,celline)], by='COND_ID')
  p <- ggplot(FigData[Fisher_P>locus_th & N_test >= minSNP],aes(x= mid_pos,y=-log10(Fisher_P)))
	p <- p + geom_point(col='lightgrey',size=.2)
  p <- p + geom_vline(col='lightgrey', xintercept = chrom_length$chrom_end+5e5, linetype=3)
  p <- p + geom_point(data=FigData[Fisher_P<=locus_th,], aes(x= mid_pos,y=-log10(Fisher_P),fill=COND_ID,shape=measure, size=FDR<.1))
  p <- p + scale_shape_manual(values=c('ExpressionIncreasing'=24,'ExpressionDecreasing'=25))
  p <- p + scale_size_manual(values=c('TRUE'=2,'FALSE'=.7))
  p <- p + scale_fill_manual(values=color_setup_simplified_norep)
  p <- p + scale_color_manual(values=color_setup_simplified_norep)+theme_plot(fontsize=11)
  fig_annotations <- TableS5d[Fisher_P<=gene_th & N_test >= minSNP_gene  & paste(COND_ID,introgression_locus, measure) %in% FigData[Fisher_P<locus_th,paste(COND_ID,introgression_locus, measure)]  ,.(introgression_locus,  COND_ID, measure, gene_name, gene_type, Fisher_P)]
  if(IMMUNE_ONLY){
    fig_annotations <- fig_annotations[grepl('i|a|p',gene_type),]
  }
  #fig_annotations <- TableS5d[Fisher_P<=0.05 & N_test >= minSNP & paste(COND_ID,introgression_locus, measure) %in% FigData[Fisher_P<.01,paste(COND_ID,introgression_locus, measure)]  ,.(introgression_locus,  COND_ID, measure, gene_name, Fisher_P)]
  fig_annotations <- fig_annotations[order(introgression_locus,COND_ID,Fisher_P), .(locus_genes=paste(unique(gene_name),collapse='/')), by=.(introgression_locus,COND_ID,measure)]
  fig_annotations=merge(fig_annotations,FigData[Fisher_P<=locus_th,.(mid_pos,Fisher_P,introgression_locus,COND_ID,celline,measure)],by=c('introgression_locus','COND_ID','measure'),allow.cartesian=TRUE)
    p <- p + geom_text_repel(data = fig_annotations, aes(x = mid_pos, y = -log10(Fisher_P), label = locus_genes, col=COND_ID), size = 2, segment.size = 0.1)
    # p <- p + geom_rect(aes(xmin = start_pos, xmax = end_pos, ymin = Pct_lo, ymax = Pct_hi, fill = celltype))
    # p <- p + geom_segment(aes(x = start_pos, xend = end_pos, y = Pct, yend = Pct), linewidth = 2)
    p <- p + scale_x_continuous(label = chrom_length$chromosome , breaks = chrom_length$chrom_mid)
    p <- p + ylab('emVar enrichment at the locus [-log10(P)]') + xlab('Chromosome')
    p <- p + facet_grid(rows = vars(celline))
    p <- p + theme(legend.box = "vertical")
		p <- p + ylim(c(0,max(4,max(FigData[N_test >= minSNP,-log10(Fisher_P)]))))
    pdf(sprintf("%s/Pct_emVar_per_locus/03a_EmvarEnrichent_of_introgressed_loci_FDR5_gene_annot__byCondition_suggestive_v2.pdf", FIGURE_DIR), height = 6, width = 9)
      print(p)
    dev.off()

cat('All done !\n')
q('no')