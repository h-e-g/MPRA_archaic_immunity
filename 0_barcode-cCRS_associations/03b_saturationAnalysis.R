
.libPaths(c('/pasteur/helix/projects/evo_immuno_pop/single_cell/resources/R_libs/4.1.0/','/pasteur/sonic/homes/mrotival/R/x86_64-pc-linux-gnu-library/4.0'))
shhh <- suppressPackageStartupMessages

shhh(library(dplyr))
shhh(library(ShortRead))
shhh(library(data.table))
shhh(library(Biostrings))
shhh(library(ggplot2))
shhh(library(ggrastr))
library(GGally)
######################################################
######             Saturation analysis           #####
######################################################
MPRA_DIR="/pasteur/helix/projects/evo_immuno_pop/MPRA/"
MPRA_DIR_KERNER=sprintf("%s/outs/RESULTS_gakerner",MPRA_DIR)
FIGURE_DIR=sprintf('%s/figures/02_barcode_lib_2',MPRA_DIR)

OUT_LIB_0_1_2_RAW=sprintf('%s/MPRAflow/outs/LIB012_coords_to_barcodes_raw.tsv.gz',MPRA_DIR)


library(data.table)
library(stringr)
source(sprintf('%s/scripts/misc_plots.R',MPRA_DIR))
### load the output from MPRAflow association
prefix="run2_3lanes_OneIdentifier"
Associations=list()
Associations[['LIB0']]=fread(sprintf('%s/MPRAflow/outs/%s/%s_coords_to_barcodes.tsv',MPRA_DIR_KERNER,prefix,prefix),header=F)
setnames(Associations[['LIB0']],c('V1','V2'),c('CRS','barcode'))

prefix='assoc_B7505-4_LIB1_noMapQ_270M_exactMatch'
Associations[['LIB1']]=fread(sprintf('%s/MPRAflow/outs/%s/%s_coords_to_barcodes.tsv.gz',MPRA_DIR,prefix,prefix),header=F)
setnames(Associations[['LIB1']],c('V1','V2'),c('CRS','barcode'))
#
prefix='assoc_B7505-4_LIB2_noMapQ_270M_exactMatch'
Associations[['LIB2']]=fread(sprintf('%s/MPRAflow/outs/%s/%s_coords_to_barcodes.tsv.gz',MPRA_DIR,prefix,prefix),header=F)
setnames(Associations[['LIB2']],c('V1','V2'),c('CRS','barcode'))

Associations=rbindlist(Associations,idcol='LIB')
setnames(Associations,'barcode','BC')

Associations_LIB1_LIB2=Associations[LIB!="LIB0",]
Associations_LIB1_LIB2[,LIB:='LIB1+LIB2']

Associations=rbind(Associations,Associations_LIB1_LIB2)

Associations_raw=Associations[,.N,by=.(BC,CRS,LIB)]
Associations_raw[,total_count_BC:=sum(N),by=.(BC,LIB)]
Associations_raw[,Pct:=N/sum(N),by=.(BC,LIB)]
Associations_raw[,type:=case_when(N>=3 & Pct==1 ~ '1a_assoc_unique',
                                  N>=3 & Pct>.5 ~ '1b_assoc_multiple',
                                  N>=3 & Pct<=.5 ~ '2_ambiguous',
                                  N<3~ '3_low_coverage')]

CRS_stats=Associations_raw[,.(UniqueBC_perCRS=length(unique(BC[substr(type,1,2)=='1a'])),
                    associatedBC_perCRS=length(unique(BC[substr(type,1,1)=='1']))),by=.(CRS,LIB)]

CRS_stats[,.(UniqueBC=sum(UniqueBC_perCRS),
              assocBC=sum(associatedBC_perCRS),
              median_UniqueBC_perCRS=as.double(median(UniqueBC_perCRS)),
              median_assocBC_perCRS=as.double(median(associatedBC_perCRS)),
              mean_UniqueBC_perCRS=as.double(mean(UniqueBC_perCRS)),
              mean_assocBC_perCRS=as.double(mean(associatedBC_perCRS)),
              Pct_CRS_with_Over10UniqueBC=mean(UniqueBC_perCRS>10),
              Pct_CRS_with_Over10assocBC=mean(associatedBC_perCRS>10),
              Pct_CRS_with_Over30UniqueBC=mean(UniqueBC_perCRS>30),
              Pct_CRS_with_Over30assocBC=mean(associatedBC_perCRS>30)),by=LIB]

# LIB           UniqueBC assocBC median_UniqueBC_perCRS median_assocBC_perCRS mean_UniqueBC_perCRS mean_assocBC_perCRS Pct_CRS_with_Over10UniqueBC Pct_CRS_with_Over10assocBC Pct_CRS_with_Over30UniqueBC Pct_CRS_with_Over30assocBC
# 1:      LIB0  2448645 2936332                    107                   129             211.1993            253.2631                   0.8882180                  0.9026220                   0.7649646                  0.7928239
# 2:      LIB1  1339206 1600586                     80                    96             115.1213            137.5901                   0.9656151                  0.9730938                   0.8281613                  0.8667584
# 3:      LIB2  1394146 1687713                     84                   102             119.8235            145.0548                   0.9666523                  0.9742157                   0.8369575                  0.8770090
# 4: LIB1+LIB2  2197059 2854994                    132                   172             188.6859            245.1901                   0.9823944                  0.9870319                   0.9179835                  0.9470972

Associations_raw[,nT_BC:=str_count(BC,'T')]
Associations_raw[,nA_BC:=str_count(BC,'A')]
Associations_raw[,nG_BC:=str_count(BC,'G')]
Associations_raw[,nC_BC:=str_count(BC,'C')]
Associations_raw[,shannon:=-(nT_BC*log((1+nT_BC)/16)+nC_BC*log((1+nC_BC)/16)+nG_BC*log((1+nG_BC)/16)+nA_BC*log((1+nA_BC)/16))]

fwrite(Associations_raw,file=OUT_LIB_0_1_2_RAW,sep='\t')

####
p <- ggplot(Associations_raw[,.N,by=.(type,LIB)],aes(x=LIB,y=N,fill=type))+geom_bar(stat='Identity',position='stack')
p <- p + theme_plot(rotate.x=90,lpos='right') + scale_fill_manual(values=setNames(quality_5levels[-1],Associations_raw[,sort(unique(type))]))
pdf(sprintf('%s/01_BC_counts_byLIB.pdf',FIGURE_DIR),width=4,height=4)
print(p)
dev.off()

p <- ggplot(Associations_raw,aes(x=type,y=N,fill=LIB,col=LIB))+geom_violin(scale='width',alpha=0.5,col='black')+rasterize(geom_boxplot(notch=T,fill=alpha('white',0.5)),dpi=200)
p <- p + scale_y_continuous(trans='log10') + theme_plot(rotate.x=90)
pdf(sprintf('%s/02_reads_perBC_byLIB.pdf',FIGURE_DIR),width=4,height=4)
print(p)
dev.off()

####
p <- ggplot(Associations_raw,aes(x=type,y=shannon,fill=LIB,col=LIB))+geom_violin(scale='width',alpha=0.5,col='black')+rasterize(geom_boxplot(notch=T,fill=alpha('white',0.5)),dpi=200)
p <- p + theme_plot(rotate.x=90)
pdf(sprintf('%s/03_entropy_perBC_byLIB.pdf',FIGURE_DIR),width=4,height=4)
print(p)
dev.off()


####
p <- ggplot(Associations_raw[,.SD[sample(1:.N,10000)],by=LIB],aes(x=N,y=shannon,col=LIB))+rasterize(geom_point(alpha=.1),dpi=200)
p <- p + scale_x_continuous(trans='log10') + theme_plot(rotate.x=90)
pdf(sprintf('%s/04_entropy_vs_reads_byLIB.pdf',FIGURE_DIR),width=4,height=4)
print(p)
dev.off()

p <- ggplot(BC_perCRS,aes(x=LIB,y=pmin(nBC,500),fill=LIB,col=LIB))+geom_violin(scale='width',alpha=0.5)+rasterize(geom_boxplot(notch=T,fill=alpha('white',0.5),col='black'),dpi=200)
p <- p + scale_y_continuous() + theme_plot(rotate.x=90) + geom_hline(yintercept=30,linetype=2,col='lightgrey')
pdf(sprintf('%s/05a_BCperCRS_byLIB.pdf',FIGURE_DIR),width=4,height=4)
print(p)
dev.off()

BC_perCRS=Associations_raw[type=='1a_assoc_unique',.(nBC=length(unique(BC))),by=.(CRS,LIB)]
p <- ggplot(BC_perCRS,aes(x=LIB,y=nBC,fill=LIB,col=LIB))+geom_violin(scale='width',alpha=0.5)+rasterize(geom_boxplot(notch=T,fill=alpha('white',0.5),col='black'),dpi=200)
p <- p + scale_y_continuous(trans='log10') + theme_plot(rotate.x=90) + geom_hline(yintercept=30,linetype=2,col='lightgrey')
pdf(sprintf('%s/05b_BCperCRS_byLIB_log.pdf',FIGURE_DIR),width=4,height=4)
print(p)
dev.off()

BC_perCRS=merge(BC_perCRS,CRS_annot[,.(CRS,GC,type)],by='CRS',all.x=TRUE)
CRS_annot=fread(sprintf('%s/data/01_barcode_lib_1/02_MPRA_count_exp2-analysis2/CRS_annot_v3.txt',MPRA_DIR))

pdf(sprintf('%s/06_BCperCRS_byGC_log.pdf',FIGURE_DIR),width=4,height=4)
p <- ggplot(BC_perCRS,aes(x=GC,y=nBC,col=LIB))+rasterize(geom_point(alpha=0.1),dpi=200)
p <- p + scale_y_continuous(trans='log10') + theme_plot() + geom_smooth(method='lm',col='black') + facet_wrap(~LIB)
print(p)
dev.off()

OK_LIB1=Associations_raw[LIB=='LIB1' & type=="1a_assoc_unique",unique(BC)]
OK_LIB2=Associations_raw[LIB=='LIB2' & type=="1a_assoc_unique",unique(BC)]
OK_LIB12=Associations_raw[LIB=='LIB1+LIB2' & type=="1a_assoc_unique",unique(BC)]

pdf(sprintf('%s/07a_BCperCRS_LIBpairs_log.pdf',FIGURE_DIR),width=7,height=7)
BC_compared=dcast(BC_perCRS,CRS+cut(GC,quantile(GC,seq(0,1,l=6),na.rm=T))~LIB,value.var='nBC',fill=0,fun.agg=first)
p <- ggpairs(BC_compared,columns=c(3,4,6,5), aes(color = GC, alpha = 0.5))
p <- p + scale_y_continuous(trans='log10') + scale_x_continuous(trans='log10') + theme_plot(rotate.x=90)
print(p)
dev.off()

pdf(sprintf('%s/07b_BCperCRS_LIBpairs.pdf',FIGURE_DIR),width=7,height=7)
BC_compared=dcast(BC_perCRS,CRS+cut(GC,quantile(GC,seq(0,1,l=6),na.rm=T))~LIB,value.var='nBC',fill=0,fun.agg=first)
p <- ggpairs(BC_compared,columns=c(3,4,6,5), aes(color = GC, alpha = 0.5))
p <- p + theme_plot(rotate.x=90)
print(p)
dev.off()

pdf(sprintf('%s/07c_BCperCRS_LIBpairs_log2.pdf',FIGURE_DIR),width=7,height=7)
BC_perCRS[,log_nBC:=log10(nBC)]
BC_compared=dcast(BC_perCRS,CRS+cut(GC,quantile(GC,seq(0,1,l=6),na.rm=T))~LIB,value.var='log_nBC',fill=0,fun.agg=first)
p <- ggpairs(BC_compared,columns=c(3,4,6,5), aes(color = GC, alpha = 0.5))
p <- p + theme_plot(rotate.x=90)
print(p)
dev.off()

sharedBC_LIB12=intersect(Associations_raw[LIB=='LIB1' & type=='1a_assoc_unique',BC],Associations_raw[LIB=='LIB2' & type=='1a_assoc_unique',BC])
sharedASSOC_LIB12=intersect(Associations_raw[LIB=='LIB1' & type=='1a_assoc_unique',paste(BC,CRS)],Associations_raw[LIB=='LIB2' & type=='1a_assoc_unique',paste(BC,CRS)])


Associations_raw[BC %chin% shared_LIB12 & LIB=='LIB1+LIB2' & type!="1a_assoc_unique",]
# BC                              CRS       LIB  N total_count_BC       Pct              type
# 1: AAGCCAAAACCAAAA        2:242799262_G_-_0_archaic LIB1+LIB2  7             17 0.4117647       2_ambiguous
# 2: CACCCAAAAAAGTCG        2:242799262_G_-_0_archaic LIB1+LIB2  6             11 0.5454545 1b_assoc_multiple
# 3: GACACCAAAAACATA        2:242799262_G_-_0_archaic LIB1+LIB2 45             75 0.6000000 1b_assoc_multiple
# 4: TCCGCCCCCAACGCC         7:43698816_T_-_0_archaic LIB1+LIB2 21             43 0.4883721       2_ambiguous
# 5: CTTTGTTTCATTTAG         4:38833173_T_-_0_archaic LIB1+LIB2  5              8 0.6250000 1b_assoc_multiple
# ---
# 37574: CACCTTACCGCCACC        16:71880814_G_+_0_archaic LIB1+LIB2  6             12 0.5000000       2_ambiguous
# 37575: CCCCTGCCCTTTTCT         1:117551429_A_-_0_Purged LIB1+LIB2 16             21 0.7619048 1b_assoc_multiple
# 37576: GATTACATAATTGCA         13:78411062_A_+_0_Purged LIB1+LIB2  3              8 0.3750000       2_ambiguous
# 37577: GCCGCCGAAATCACA          2:162940017_A_+_0_COVID LIB1+LIB2 23             26 0.8846154 1b_assoc_multiple
# 37578: CCCAACTCCCAAACC 14:35873947-35873965_-_scrambled LIB1+LIB2  3              7 0.4285714       2_ambiguous

Associations_raw[BC=='AAGCCAAAACCAAAA',]
#  BC                       CRS       LIB  N total_count_BC       Pct              type
# 1: AAGCCAAAACCAAAA  6:30163095_T_-_0_archaic      LIB0  1              1 1.0000000    3_low_coverage
# 2: AAGCCAAAACCAAAA 2:242799262_G_-_0_archaic      LIB1  7              7 1.0000000   1a_assoc_unique
# 3: AAGCCAAAACCAAAA 15:79224509_C_-_0_archaic      LIB2 10             10 1.0000000   1a_assoc_unique
# 4: AAGCCAAAACCAAAA 2:242799262_G_-_0_archaic LIB1+LIB2  7             17 0.4117647       2_ambiguous
# 5: AAGCCAAAACCAAAA 15:79224509_C_-_0_archaic LIB1+LIB2 10             17 0.5882353 1b_assoc_multiple

Associations_raw[BC=='GCCGCCGAAATCACA',]
# BC                                                 CRS       LIB  N total_count_BC       Pct              type nT_BC nA_BC nG_BC nC_BC  shannon
# 1: GCCGCCGAAATCACA 6:32369123_C_-_-50_archaic|6:32369173_T_-_0_archaic      LIB1  3              3 1.0000000   1a_assoc_unique     1     5     3     6 16.10254
# 2: GCCGCCGAAATCACA                             2:162940017_A_+_0_COVID      LIB2 23             23 1.0000000   1a_assoc_unique     1     5     3     6 16.10254
# 3: GCCGCCGAAATCACA 6:32369123_C_-_-50_archaic|6:32369173_T_-_0_archaic LIB1+LIB2  3             26 0.1153846       2_ambiguous     1     5     3     6 16.10254
# 4: GCCGCCGAAATCACA                             2:162940017_A_+_0_COVID LIB1+LIB2 23             26 0.8846154 1b_assoc_multiple     1     5     3     6 16.10254


# How many barcode per CRS ?
# compute & display the # of reads ber barcode-CRS
# compute & display the # of barcodes per CRS (pre filtering & post filtering, with varying thresholds)
### compute & display the # of CRS per barcode for varying thresholds
# subsample the association dictionary & plot the number of filtered barcode as we increase the sample size (saturation analysis)

filter_barcode=function(Associations){
  BC_count=Associations[,.(count=.N),by=.(CRS,BC,LIB)]
  BC_count[,total_count:=sum(count),by=.(BC,LIB)]
  BC_count[count>=3 & count/total_count==1,]
}
percentage_sample=seq(0.1,1,by=.1)
DT=list()
for (i in 1:length(percentage_sample)){
  sub=Associations[sample(1:.N,.N*percentage_sample[i],replace=F),by=LIB]
  nReads=sub[,.(nReads=.N),by=LIB]
  filtered=filter_barcode(sub)
  DT[[i]]=filtered[,.(N_assoc=.N,
              N_CRS=length(unique(CRS)),
              median_N_barcodes_perCRS=.SD[,.N,by=CRS][,as.double(median(N))],
              mean_N_barcodes_perCRS=.SD[,.N,by=CRS][,as.double(mean(N))],
              N_excessBC_over200perCRS=.SD[,.N,by=CRS][,sum(N[N>=200]-200)],
              N_CRS_over10BC=.SD[,.N,by=CRS][,sum(N>=10)],
              N_CRS_over30BC=.SD[,.N,by=CRS][,sum(N>=30)],
              N_CRS_over50BC=.SD[,.N,by=CRS][,sum(N>=50)],
              N_CRS_over100BC=.SD[,.N,by=CRS][,sum(N>=100)],
              N_CRS_over200BC=.SD[,.N,by=CRS][,sum(N>=200)]),
              by=LIB]
  DT[[i]]=merge(nReads,DT[[i]],keyby=LIB)
  cat(i)
}

DT=rbindlist(DT,idcol="subsample")
DT[,percentage_sample:=percentage_sample[subsample]]
DT[,percentage_excessBC:=N_excessBC_over200perCRS/N_assoc]


fwrite(DT,file=sprintf('%s/SaturationAnalysis.txt',FIGURE_DIR),sep='\t')

DT_long=melt(DT,id.var=c('nReads','percentage_sample','subsample','LIB'))

DT_long[,group:=case_when(
  grepl('N_CRS',variable)~ 'CRS',
  grepl('N_barcode',variable) ~ "barcode per CRS",
    grepl('N_assoc',variable) ~ "Association",
    grepl('percentage_excessBC',variable) ~ "Percentage excess barcodes",
  TRUE ~ 'other')]

pdf(sprintf('%s/08a_SaturationAnalysis_barcodes.pdf',FIGURE_DIR),height=4)
p <- ggplot(DT_long[group!='CRS' & group!='other'],aes(x=nReads,y=value,col=LIB,shape=variable))+facet_wrap(~variable,scales='free')+theme_plot(lpos='right')+geom_point()+geom_line()
print(p)
dev.off()

pdf(sprintf('%s/08b_SaturationAnalysis_CRS.pdf',FIGURE_DIR),height=4)
p <- ggplot(DT_long[group=='CRS'],aes(x=nReads,y=value,col=LIB,shape=variable))+facet_grid(cols=vars(variable),scales='free')+theme_plot(lpos='bottom')+geom_point()+geom_line()
print(p)
dev.off()


count_assoc=Associations_raw[,.(nBC_found=length(unique(BC)),
                    n_uniqBC=length(unique(BC[type=='1a_assoc_unique'])),
                    nReads_association=sum(N),
                    nReads_uniqueBC=sum(N[type=='1a_assoc_unique']))
                    ,by=.(barcode_library)]
count_assoc[,Pct_uniq_BC:=n_uniqBC/nBC_found]
count_assoc[,Pct_read_uniq_BC:=nReads_uniqueBC/nReads_association]