
MPRA_DIR <- "/pasteur/helix/projects/evo_immuno_pop/MPRA"
ANALYSIS_DIR <- "MPRA_count_exp6_analysis3"

source(sprintf("%s/scripts/%s/00_starter.R", MPRA_DIR, ANALYSIS_DIR))
source(sprintf("%s/scripts/misc_plots.R", MPRA_DIR))
source(sprintf("%s/scripts/%s/00_load_associations.R", MPRA_DIR, ANALYSIS_DIR))

# Associations_Filtered
# UMI_counts=fread(sprintf('%s/data/%s/01c_MPRA_results_BC_activity/01c_MPRA_results_BC_activity__HepG2_NS_R1__v1.txt.gz',MPRA_DIR,ANALYSIS_DIR))
# UMI_counts=fread(sprintf('%s/data/%s/01c_MPRA_results_BC_activity/01c_MPRA_results_BC_activity__HepG2_TNFa_R1__v1.txt.gz',MPRA_DIR,ANALYSIS_DIR))
UMI_counts <- list()
for (SAMPLE_ID in names(color_setup)) {
  # load UMI per BC
  UMI_counts[[SAMPLE_ID]] <- fread(sprintf("%s/data/%s/01c_MPRA_results_BC_activity/01c_MPRA_results_BC_activity__%s__v1.txt.gz", MPRA_DIR, ANALYSIS_DIR, SAMPLE_ID))
}
UMI_counts <- rbindlist(UMI_counts, idcol = "library")
UMI_counts <- merge(UMI_counts, condition_summary_reps, by.x='library', by.y = "analysis_name")
# load annotations of oligos
oligo_source <- fread(sprintf("%s/data/%s/00_oligo_annot_v2.txt", MPRA_DIR, ANALYSIS_DIR))
UMI_counts=merge(UMI_counts,oligo_source[,.(oligo,GC)])

UMI_counts[,cor.test(DNA,GC,method='s')]
#        rho 
# 0.00599073 & p-value < 2.2e-16

UMI_counts[,cor.test(RNA,GC,method='s')]
#       rho 
# 0.1556474 & p-value < 2.2e-16


UMI_counts[grepl('scrambled',oligo) & DNA>0,cor.test(DNA,GC,method='s')]
UMI_counts[grepl('scrambled',oligo),cor.test(DNA,GC,method='s')]
UMI_counts[!grepl('scrambled',oligo) & DNA>0,cor.test(DNA,GC,method='s')]
UMI_counts[!grepl('scrambled',oligo),cor(DNA,GC,method='s')]

UMI_counts[,.(rho_test_DNA=cor(DNA[!grepl('scrambled|Promoter',oligo)],GC[!grepl('scrambled|Promoter',oligo)],method='s'),
							rho_scrambled_DNA=cor(DNA[grepl('scrambled',oligo)],GC[grepl('scrambled',oligo)],method='s'),
							rho_test_RNA=cor(RNA[!grepl('scrambled|Promoter',oligo)],GC[!grepl('scrambled|Promoter',oligo)],method='s'),
							rho_scrambled_RNA=cor(RNA[grepl('scrambled',oligo)],GC[grepl('scrambled',oligo)],method='s'),
							rho_test_DNAnon0=cor(DNA[DNA>0 & !grepl('scrambled|Promoter',oligo)],GC[DNA>0 & !grepl('scrambled|Promoter',oligo)],method='s'),
							rho_scrambled_DNAnon0=cor(DNA[DNA>0 & grepl('scrambled',oligo)],GC[DNA>0 & grepl('scrambled',oligo)],method='s'),
							rho_test_RNAnon0=cor(RNA[RNA>0 & !grepl('scrambled|Promoter',oligo)],GC[RNA>0 & !grepl('scrambled|Promoter',oligo)],method='s'),
							rho_scrambled_RNAnon0=cor(RNA[RNA>0 & grepl('scrambled',oligo)],GC[RNA>0 & grepl('scrambled',oligo)],method='s')),
							by=COND_ID.x]


UMI_counts[grepl('scrambled',oligo),.(DNA=mean(DNA),RNA=mean(RNA),RNA_overDNA=mean(RNA)/mean(DNA),PctNullDNA=mean(DNA==0),PctNullRNA=mean(RNA==0),DNA_non0=mean(DNA[DNA>0])),keyby=.(GC=cut(GC,10))]
UMI_counts[!grepl('scrambled',oligo),.(DNA=mean(DNA),RNA=mean(RNA),RNA_overDNA=mean(RNA)/mean(DNA),PctNullDNA=mean(DNA==0),PctNullRNA=mean(RNA==0),DNA_non0=mean(DNA[DNA>0])),keyby=.(GC=cut(GC,10))]
#. correlation only at the RNA level, suggest higher activity of GC rich sequyences on average 


by_GC_UMImean <- UMI_counts[grepl('scrambled',oligo),.(.N, DNA=mean(DNA),se_DNA=sd(DNA)/sqrt(.N),RNA=mean(RNA),se_RNA=sd(RNA)/sqrt(.N),RNA_overDNA=mean(RNA)/mean(DNA),PctNullDNA=mean(DNA==0),PctNullRNA=mean(RNA==0),DNA_non0=mean(DNA[DNA>0]),meanGC=mean(GC)),keyby=.(COND_ID.x,GC=cut(GC,10))]
by_GC_UMImean_test <- UMI_counts[grepl('!scrambled',oligo),.(.N, DNA=mean(DNA),se_DNA=sd(DNA)/sqrt(.N),RNA=mean(RNA),se_RNA=sd(RNA)/sqrt(.N),RNA_overDNA=mean(RNA)/mean(DNA),PctNullDNA=mean(DNA==0),PctNullRNA=mean(RNA==0),DNA_non0=mean(DNA[DNA>0]),meanGC=mean(GC)),keyby=.(COND_ID.x,GC=cut(GC,10))]


UMI_counts_scrambled=UMI_counts[grepl('scrambled',oligo),]
by_GC_UMImean <- UMI_counts_scrambled[,.(.N, DNA=mean(Winsorize(as.numeric(DNA))),se_DNA=sd(Winsorize(as.numeric(DNA)))/sqrt(.N),RNA=mean(as.numeric(Winsorize(as.numeric(RNA)))),se_RNA=sd(Winsorize(as.numeric(RNA)))/sqrt(.N),meanGC=mean(as.numeric(GC))),keyby=.(COND_ID.x,GCbin=cut(GC,10))]

p <- ggplot(by_GC_UMImean,aes(x=meanGC,y=DNA,ymin=DNA-2*se_DNA, ymax=DNA+2*se_DNA,col=COND_ID.x))+geom_line()+geom_pointrange()
p <- p + theme_plot(lpos='right')+scale_color_manual(values=color_setup_simplified_norep)+ylab('mean DNA UMI count per barcode ')

pdf(sprintf('%s/figures/%s/GCbias_UMI_DNA_counts_scrambled.pdf',MPRA_DIR,ANALYSIS_DIR))
print(p)
dev.off()


p <- ggplot(by_GC_UMImean,aes(x=meanGC,y=RNA,ymin=RNA-2*se_RNA, ymax=RNA+2*se_RNA,col=COND_ID.x))+geom_line()+geom_pointrange()
p <- p + theme_plot(lpos='right')+scale_color_manual(values=color_setup_simplified_norep)+ylab('mean RNA UMI count per barcode ')

pdf(sprintf('%s/figures/%s/GCbias_UMI_RNA_counts_scrambled.pdf',MPRA_DIR,ANALYSIS_DIR))
print(p)
dev.off()


by_GC_UMImean <- UMI_counts_scrambled[,.(.N, DNA=mean(Winsorize(as.numeric(DNA))),se_DNA=sd(Winsorize(as.numeric(DNA)))/sqrt(.N),RNA=mean(as.numeric(Winsorize(as.numeric(RNA)))),se_RNA=sd(Winsorize(as.numeric(RNA)))/sqrt(.N),meanGC=mean(as.numeric(GC))),keyby=.(COND_ID.x,GCbin=cut(GC,breaks=quantile(GC,seq(0,1,l=6)),include.lowest=T))]

p <- ggplot(by_GC_UMImean,aes(x=meanGC,y=DNA,ymin=DNA-2*se_DNA, ymax=DNA+2*se_DNA,col=COND_ID.x))+geom_line()+geom_pointrange()
p <- p + theme_plot(lpos='right')+scale_color_manual(values=color_setup_simplified_norep)+ylab('mean DNA UMI count per barcode ')

pdf(sprintf('%s/figures/%s/GCbias_UMI_DNA_counts_scrambled_quantiles.pdf',MPRA_DIR,ANALYSIS_DIR))
print(p)
dev.off()


p <- ggplot(by_GC_UMImean,aes(x=meanGC,y=RNA,ymin=RNA-2*se_RNA, ymax=RNA+2*se_RNA,col=COND_ID.x))+geom_line()+geom_pointrange()
p <- p + theme_plot(lpos='right')+scale_color_manual(values=color_setup_simplified_norep)+ylab('mean RNA UMI count per barcode ')

pdf(sprintf('%s/figures/%s/GCbias_UMI_RNA_counts_scrambled_quantiles.pdf',MPRA_DIR,ANALYSIS_DIR))
print(p)
dev.off()


p <- ggplot(UMI_counts_scrambled,aes(x=GC,y=DNA))+geom_point()+geom_smooth(method='lm')
p <- p + theme_plot(lpos='right')+ylab('mean RNA UMI count per barcode ')
pdf(sprintf('%s/figures/%s/GCbias_UMI_DNA_counts_scrambled_scatter.pdf',MPRA_DIR,ANALYSIS_DIR))
print(p)
dev.off()


p <- ggplot(UMI_counts_scrambled,aes(x=GC,y=RNA))+geom_point()+geom_smooth(method='lm')
p <- p + theme_plot(lpos='right')+ylab('mean RNA UMI count per barcode ')

pdf(sprintf('%s/figures/%s/GCbias_UMI_RNA_counts_scrambled_scatter.pdf',MPRA_DIR,ANALYSIS_DIR))
print(p)
dev.off()

by_GC_UMImean <- UMI_counts[,.(.N, DNA=mean(Winsorize(as.numeric(DNA))),se_DNA=sd(Winsorize(as.numeric(DNA)))/sqrt(.N),RNA=mean(as.numeric(Winsorize(as.numeric(RNA)))),se_RNA=sd(Winsorize(as.numeric(RNA)))/sqrt(.N),meanGC=mean(as.numeric(GC))),keyby=.(COND_ID.x,GCbin=cut(GC,10))]


p <- ggplot(by_GC_UMImean,aes(x=meanGC,y=DNA,ymin=DNA-2*se_DNA, ymax=DNA+2*se_DNA,col=COND_ID.x))+geom_line()+geom_pointrange()
p <- p + theme_plot(lpos='right')+scale_color_manual(values=color_setup_simplified_norep)+ylab('mean DNA UMI count per barcode ')

pdf(sprintf('%s/figures/%s/GCbias_UMI_DNA_counts_all.pdf',MPRA_DIR,ANALYSIS_DIR))
print(p)
dev.off()


p <- ggplot(by_GC_UMImean,aes(x=meanGC,y=RNA,ymin=RNA-2*se_RNA, ymax=RNA+2*se_RNA,col=COND_ID.x))+geom_line()+geom_pointrange()
p <- p + theme_plot(lpos='right')+scale_color_manual(values=color_setup_simplified_norep)+ylab('mean RNA UMI count per barcode ')

pdf(sprintf('%s/figures/%s/GCbias_UMI_RNA_counts_all.pdf',MPRA_DIR,ANALYSIS_DIR))
print(p)
dev.off()


p <- ggplot(by_GC_UMImean,aes(x=meanGC,y=RNA/DNA,col=COND_ID.x))+geom_line()+geom_point()
p <- p + theme_plot(lpos='right')+scale_color_manual(values=color_setup_simplified_norep)+ylab('mean RNA/DNA UMI count per barcode ')

pdf(sprintf('%s/figures/%s/GCbias_UMI_RNA_DNA_counts_all.pdf',MPRA_DIR,ANALYSIS_DIR))
print(p)
dev.off()


############## numbers, suppelmentary Note on GC



GC_by_type=merge(tested_and_ctrl_oligos_final,oligo_source[,.(oligo,GC)],by='oligo')[,.(GC=mean(GC)),by=.(type=ifelse(type=='ctrl','scrambled','tested'),posID)]
GC_by_type[,mean(GC),by=type]
# 1: tested 0.4423730
# 2:   ctrl 0.5887387

GC_by_type[,wilcox.test(GC[type=='tested'],GC[type=='ctrl'])$p.value]
#  1.621327e-39

	p <- ggplot(GC_by_type,aes(x=factor(type,c('tested','scrambled')),y=GC))+geom_violin(fill='grey',col=NA)+geom_boxplot(width=0.1,fill='black')
	p <- p + theme_plot(lpos='right')+ xlab('')

	pdf(sprintf('%s/figures/%s/GCbias_scrambled_vs_test.pdf',MPRA_DIR,ANALYSIS_DIR),width=2.5,heigh=4)
	print(p)
	dev.off()


RUN_ID <- "RUN3_Z2_nBC10"
CRITERIA_ACTIVE <- "FDR5_FC1_FC0.2"
ACTIVE_DIR <- sprintf("%s/data/%s/03b_computeFDR/%s/Activity/%s/", MPRA_DIR, ANALYSIS_DIR, RUN_ID, CRITERIA_ACTIVE)


oligo_activity_obs <- fread(sprintf("%s/all_oligos_annotated__%s.tsv.gz", ACTIVE_DIR, CRITERIA_ACTIVE))
oligo_activity_obs_ctc <- oligo_activity_obs[power == 0 & boot == 0 & ANALYSIS_SUBTYPE == "celltype_cond", ]

oligo_activity_perm <- fread(sprintf("%s/all_oligos_annotated_perm__%s.tsv.gz", ACTIVE_DIR, CRITERIA_ACTIVE))
oligo_activity_perm_ctc <- oligo_activity_perm[power == 0 & boot == 0 & ANALYSIS_SUBTYPE == "celltype_cond", ]
oligo_activity_obs_ctc_archaic <- oligo_activity_obs_ctc[oligo %in% tested_and_ctrl_oligos_final[type == "tested", oligo], ]

oligo_activity_obs_ctc[,median(log2(alpha_CRITERIA)[allele.label=="Scrambled"]),by=COND_ID]
compared_activvity_scrambled_tested <- oligo_activity_obs_ctc[,.(scrambled_med=median(log2(alpha_CRITERIA)[allele.label=="scrambled"],na.rm=T),
																					tested_med=median(log2(alpha_CRITERIA)[grepl('INTROGRESSED',allele.label)],na.rm=T),
																					diff_P=wilcox.test(log2(alpha_CRITERIA)[grepl('INTROGRESSED',allele.label)], log2(alpha_CRITERIA)[allele.label=="scrambled"])$p.value
																				),keyby=COND_ID]

compared_activvity_scrambled_tested[,min(scrambled_med-tested_med)]
# [1] 0.143681
compared_activvity_scrambled_tested[,max(diff_P)]
# [1] 5.559766e-06

oligo_activity_obs_ctc[allele.label=="scrambled",.(beta=lm(log2(alpha_CRITERIA)~I(GC*270))$coeff[2]),by=COND_ID][,range(beta)]
oligo_activity_obs_ctc[allele.label=="scrambled" ,.(P=summary(lm(log2(alpha_CRITERIA)~I(GC*270)))$coeff[2,4]),by=COND_ID][,range(P)]

oligo_activity_obs_ctc[allele.label=="scrambled" ,.(P=summary(lm(log2(alpha_CRITERIA)~I(GC*100)))$coeff[2,c(1,4)])]

#               P
#           <num>
# 1: 1.724567e-02
# 2: 2.708120e-19
CRITERIA_ACTIVE <- "FDR5_FC1_FC0.2"
ACTIVE_DIR <- sprintf("%s/data/%s/03b_computeFDR/%s/Activity/%s/", MPRA_DIR, ANALYSIS_DIR, RUN_ID, CRITERIA_ACTIVE)
oligo_activity_obs <- fread(sprintf("%s/all_oligos_annotated__%s.tsv.gz", ACTIVE_DIR, CRITERIA_ACTIVE))
oligo_activity_obs_ctc <- oligo_activity_obs[power == 0 & boot == 0 & ANALYSIS_SUBTYPE == "celltype_cond", ]
raw_active_CRE <- oligo_activity_obs_ctc[oligo_class_CRITERIA!='inactive',.(raw=length(unique(crsID))),by=COND_ID]

CRITERIA_ACTIVE <- "FDR5_FC1_FC0.2_GC"
ACTIVE_DIR <- sprintf("%s/data/%s/03b_computeFDR/%s/Activity/%s/", MPRA_DIR, ANALYSIS_DIR, RUN_ID, CRITERIA_ACTIVE)
oligo_activity_obs <- fread(sprintf("%s/all_oligos_annotated__%s.tsv.gz", ACTIVE_DIR, CRITERIA_ACTIVE))
oligo_activity_obs_ctc <- oligo_activity_obs[power == 0 & boot == 0 & ANALYSIS_SUBTYPE == "celltype_cond", ]
GC_active_CRE <- oligo_activity_obs_ctc[oligo_class_CRITERIA!='inactive',.(GCadj=length(unique(crsID))),by=COND_ID]

CRITERIA_ACTIVE <- "FDR5_FC1_FC0.2_norm"
ACTIVE_DIR <- sprintf("%s/data/%s/03b_computeFDR/%s/Activity/%s/", MPRA_DIR, ANALYSIS_DIR, RUN_ID, CRITERIA_ACTIVE)
oligo_activity_obs <- fread(sprintf("%s/all_oligos_annotated__%s.tsv.gz", ACTIVE_DIR, CRITERIA_ACTIVE))
oligo_activity_obs_ctc <- oligo_activity_obs[power == 0 & boot == 0 & ANALYSIS_SUBTYPE == "celltype_cond", ]
norm_active_CRE <- oligo_activity_obs_ctc[oligo_class_CRITERIA!='inactive',.(scaled=length(unique(crsID))),by=COND_ID]


CRITERIA_ACTIVE <- "FDR5_FC1_FC0.2_GCnorm"
ACTIVE_DIR <- sprintf("%s/data/%s/03b_computeFDR/%s/Activity/%s/", MPRA_DIR, ANALYSIS_DIR, RUN_ID, CRITERIA_ACTIVE)
oligo_activity_obs <- fread(sprintf("%s/all_oligos_annotated__%s.tsv.gz", ACTIVE_DIR, CRITERIA_ACTIVE))
oligo_activity_obs_ctc <- oligo_activity_obs[power == 0 & boot == 0 & ANALYSIS_SUBTYPE == "celltype_cond", ]
GCnorm_active_CRE <- oligo_activity_obs_ctc[oligo_class_CRITERIA!='inactive',.(GCadj_scaled=length(unique(crsID))),by=COND_ID]

active_CREs <- merge(merge(merge(raw_active_CRE,GC_active_CRE),norm_active_CRE),GCnorm_active_CRE)
active_CREs[,increase_GCadj:=GCadj/raw*100]
active_CREs
#        COND_ID   raw GCadj scaled GCadj_scaled increase_GCadj
#           <char> <int> <int>  <int>        <int>          <num>
#  1:     A549_IAV  1006   949   1001          948       94.33400
#  2:      A549_NS  1241  1216   1234         1214       97.98550
#  3:    A549_SARS   864   817    863          817       94.56019
#  4:    A549_TNFa  1302  1265   1291         1261       97.15822
#  5:    HepG2_DEX  2330  2126   2223         2060       91.24464
#  6: HepG2_IFNA2b  2217  2096   2123         2029       94.54217
#  7:     HepG2_NS  2332  2147   2237         2106       92.06690
#  8:   HepG2_TNFa  1836  1683   1807         1668       91.66667
#  9:     K562_DEX   678  1075    706         1116      158.55457
# 10:  K562_IFNA2b   671  1033    700         1075      153.94933
# 11:      K562_NS   778  1060    873         1157      136.24679
# 12:    K562_TNFa   733  1262    756         1289      172.16917

activity_distrib=list()
CRITERIA_ACTIVE <- "FDR5_FC1_FC0.2"
ACTIVE_DIR <- sprintf("%s/data/%s/03b_computeFDR/%s/Activity/%s/", MPRA_DIR, ANALYSIS_DIR, RUN_ID, CRITERIA_ACTIVE)
oligo_activity_obs <- fread(sprintf("%s/all_oligos_annotated__%s.tsv.gz", ACTIVE_DIR, CRITERIA_ACTIVE))
activity_distrib[['obs_raw']] <- oligo_activity_obs[power == 0 & boot == 0 & ANALYSIS_SUBTYPE == "celltype_cond", .(alpha_CRITERIA, celline,  COND_ID, data='obs', normalization='raw')]
oligo_activity_perm <- fread(sprintf("%s/all_oligos_annotated_perm__%s.tsv.gz", ACTIVE_DIR, CRITERIA_ACTIVE))
activity_distrib[['perm_raw']] <- oligo_activity_perm[power == 0 & boot == 0 & ANALYSIS_SUBTYPE == "celltype_cond", .(alpha_CRITERIA, celline,  COND_ID, data='perm', normalization='raw')]


CRITERIA_ACTIVE <- "FDR5_FC1_FC0.2_GCnorm"
ACTIVE_DIR <- sprintf("%s/data/%s/03b_computeFDR/%s/Activity/%s/", MPRA_DIR, ANALYSIS_DIR, RUN_ID, CRITERIA_ACTIVE)
oligo_activity_obs <- fread(sprintf("%s/all_oligos_annotated__%s.tsv.gz", ACTIVE_DIR, CRITERIA_ACTIVE))
activity_distrib[['obs_norm']] <- oligo_activity_obs[power == 0 & boot == 0 & ANALYSIS_SUBTYPE == "celltype_cond", .(alpha_CRITERIA, celline, COND_ID, data='obs', normalization='normalized')]
oligo_activity_perm <- fread(sprintf("%s/all_oligos_annotated_perm__%s.tsv.gz", ACTIVE_DIR, CRITERIA_ACTIVE))
activity_distrib[['perm_norm']] <- oligo_activity_perm[power == 0 & boot == 0 & ANALYSIS_SUBTYPE == "celltype_cond", .(alpha_CRITERIA, celline,  COND_ID, data='perm', normalization='normalized')]
activity_distrib=rbindlist(activity_distrib)

p <- ggplot(activity_distrib, aes(x = pmax(-1,pmin(1, log2(alpha_CRITERIA))), fill = celline, alpha = data)) +
  geom_histogram()
p <- p + facet_grid(rows = vars(celline), cols=vars(factor(normalization,c('raw','normalized')))) + scale_fill_manual(values = color_celline) + scale_alpha_manual(values = c("perm" = 1, "obs" = 0.5))
p <- p + theme_plot(lpos = "right") + xlab("activity (capped at 2)")

pdf(sprintf('%s/figures/%s/activity_distribution_normalization.pdf',MPRA_DIR,ANALYSIS_DIR), height = 2.5, width = 4)
print(p)
dev.off()
########################################################################
########################################################################
########################################################################

######################## simulate data. ################################

nCRS=100
GC=rnorm(nCRS,.5,0.1)

log_activity=rnorm(nCRS)

#log_activity=rep(1,nCRS)

aGC_BC=10
uBC=100
nBC_CRS=rpois(nCRS,uBC+aGC_BC*(GC-0.5))

DT_BC=data.table(crsID=rep(1:nCRS,nBC_CRS),GC_content=rep(GC,nBC_CRS),log_activity=rep(log_activity,nBC_CRS))
DT_BC[,bcID:=1:.N,by=crsID]

uIntegrationBC=300
DT_BC[,nInt_BC:=rpois(.N,uIntegrationBC)]

DT_Frag=DT_BC[,.(fragID=1:nInt_BC),by=.(crsID, GC_content, log_activity, bcID,nInt_BC)]

DNA_complexity=1e5
QCR_cylceDNA=3
RNA_complexity=3e5
QCR_cylceRNA=3


DNA_complexity=1e5
QCR_cylceDNA=3
RNA_complexity=3e5
QCR_cylceRNA=3


DT_Frag[,DNA_prob:=(GC_content**QCR_cylceDNA)]
DT_Frag[,DNA_prob:=DNA_prob/sum(DNA_prob)]
DT_Frag[,DNA_count:=rpois(.N,DNA_complexity*DNA_prob)]
DT_Frag[,RNA_prob:=2^log_activity*(GC_content**QCR_cylceRNA)]
DT_Frag[,RNA_prob:=RNA_prob/sum(RNA_prob)]
DT_Frag[,RNA_count:=rpois(.N,RNA_complexity*RNA_prob)]

DT_BC_count=DT_Frag[,.(RNA=sum(RNA_count>0),DNA=sum(DNA_count>0)),by=.(crsID, GC_content, log_activity, bcID)]

DT_BC_count_true=DT_Frag[,.(RNA=sum(RNA_count>0),DNA=sum(DNA_count>0)),by=.(crsID, GC_content, log_activity, bcID)]

prob_sample=GC**QCR_cylceDNA
QCR_cylceRNA=15

uRNA=15
DNA <- rpois(nCRS,uDNA)
RNA <- rpois(nCRS,uRNA)