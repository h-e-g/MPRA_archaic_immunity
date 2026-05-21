tic("loading oligo & SNP annotations")
# load annotation of oligos (beforz: CRS)
oligo_source <- fread(sprintf("%s/data/%s/00_oligo_annot_v2.txt", MPRA_DIR, ANALYSIS_DIR))
#oligo_source_v1 <- fread(sprintf("%s/data/%s/00_oligo_anno_v1.txt", MPRA_DIR, ANALYSIS_DIR))
# load annotations of SNPs
SNP_annot <- fread(sprintf("%s/data/%s/00_SNP_annot_v4.txt", MPRA_DIR, ANALYSIS_DIR))
toc()

nBC_oligo=10
nBC_SNP=10

pass_selection_snps <- unlist(fread("/pasteur/helix/projects/evo_immuno_pop/MPRA/data/vcfs/all_common_introgressed/pass_selection_snps.txt"))

length(pass_selection_snps) # number of selected SNPs
# 4349

length(pass_selection_snps)+ 113 # number of tested CREs
# 4462

nTested_oligos=(length(pass_selection_snps)+ 113)*2 # number of tested oligos
nTested_oligos
# 8,924

# What are the CRS with >10 associated barcodes (here,we also remove those that turned out amibiguous on introgressed allele)
tested_snps_posID=SNP_annot[excluded==FALSE & ctrl==FALSE & no_barcodes==FALSE & (nBC_LIB0>nBC_SNP | nBC_LIB2>nBC_SNP) & (INTROGRESSED.allele!="" & is_introgressed) & ID %in% pass_selection_snps,unique(posID)] 
length(tested_snps_posID) # 4015

# how many were removed due to ambiguous introgressed allele ? 
removed_ambiguous=SNP_annot[excluded==FALSE & ctrl==FALSE & no_barcodes==FALSE & (nBC_LIB0>nBC_SNP | nBC_LIB2>nBC_SNP) & (INTROGRESSED.allele=="" & is_introgressed) & ID %in% pass_selection_snps,unique(posID)] 
length(removed_ambiguous) # 184


ctrl_snps_posID=SNP_annot[excluded==FALSE & ctrl==TRUE & no_barcodes==FALSE & (nBC_LIB0>nBC_SNP | nBC_LIB2>nBC_SNP) ,unique(posID)] 
length(ctrl_snps_posID) # 403

ctrl_oligos=oligo_source[excluded==FALSE  & ctrl==TRUE & (nBC_LIB0>nBC_oligo | nBC_LIB2>nBC_oligo) & type_simple!="other" ,.N,by=.(oligo,posID)]
nrow(ctrl_oligos) # 201

ctrl_oligos[,.N,by=grepl('scrambled',oligo)]
#     scrambled     N
# 1:  FALSE   		106
# 2:   TRUE    		95

ctrl_oligos[,.N,by=posID][N>=2,][,.N] # how many controls have paired promoter & scrambled
# 90 
ctrl_oligos[,length(unique(posID))] # how many controls CREs have at least one oligo OK 
# 111

########### get oligos and match posID
# ctrl_snps_oligos = oligo_source[ crsID %in% ctrl_snps_posID & (nBC_LIB0>nBC_oligo | nBC_LIB2>nBC_oligo),.(oligo, posID)]
# nrow(ctrl_snps_oligos) # 890

tested_snps_oligos = oligo_source[ crsID %in% tested_snps_posID & (nBC_LIB0>nBC_oligo | nBC_LIB2>nBC_oligo) & posID %in% tested_snps_posID,.(oligo, posID)]
nrow(tested_snps_oligos) # 7834


########### get lists of tested SNPs with sufficient barcodes per CRS
# ctrl_snps_posID_enoughBC = ctrl_snps_oligos[,.N,by=posID][N>=2,posID]
# length(ctrl_snps_posID_enoughBC) # 403

tested_snps_posID_enoughBC = tested_snps_oligos[,.N,by=posID][N>=2,posID]
length(tested_snps_posID_enoughBC) # 3838

ctrl_posID_enoughBC = ctrl_oligos[,unique(posID)]
length(ctrl_posID_enoughBC) # 111

length(unique(c(ctrl_posID_enoughBC,tested_snps_posID_enoughBC)))
# 3948

tested_and_ctrl_oligos=rbindlist(list(tested=tested_snps_oligos[,.(oligo,posID)],ctrl=ctrl_oligos[,.(oligo,posID)]),idcol='type')
nrow(tested_and_ctrl_oligos) # 8035

tested_and_ctrl_oligos_enoughBC_posID=tested_and_ctrl_oligos[,.N,by=.(posID,type)][N>=2 | type=='ctrl',.(posID)]
nrow(tested_and_ctrl_oligos_enoughBC_posID) # 3949

tested_and_ctrl_oligos_final <- tested_and_ctrl_oligos[posID%in%tested_and_ctrl_oligos_enoughBC_posID$posID,]
nrow(tested_and_ctrl_oligos_final) # 7876

SNP_annot[ posID %in% tested_and_ctrl_oligos_final$posID, length(unique(introgression_locus))]
# 198 loci


DT_oligo_BC=Associations_Filtered[oligo%chin%tested_and_ctrl_oligos$oligo & barcode_library%chin%c('LIB0','LIB2'),.N,by=.(barcode_library,oligo)]
DT_oligo_BC=dcast(DT_oligo_BC,oligo~barcode_library,value.var='N')
DT_oligo=merge(tested_and_ctrl_oligos,DT_oligo_BC,all.x=TRUE)
DT_oligo[is.na(LIB0),LIB0:=0]
DT_oligo[is.na(LIB2),LIB2:=0]

# After association sequencing, each oligo was associated with a median of  > 81 barcodes per plasmids libraries with >90% of oligos being associated with >10 barcodes across plasmids libraries 
length(tested_and_ctrl_oligos)/nTested_oligos
# 0.900381
DT_oligo[,median(LIB0)]# 99
DT_oligo[,median(LIB2)]# 81

DT_oligo[LIB0>10 & LIB2>10,.N,by=posID][N>=2,][,.N] 
# 3324 CREs with >10 BC per CRS in BOTH plasmid libraries

tested_and_ctrl_oligos_final_annot <- merge(tested_and_ctrl_oligos_final,oligo_source[,.(posID, oligo, allele, allele.label, strand, shift, type_selected=type, sequence, nBC_LIB1=nBC_LIB0, nBC_LIB2)],by=c('posID','oligo'))
SNPs_with_more_than_one_introgressed_allele=tested_and_ctrl_oligos_final_annot[type=='tested',length(unique(allele.label)),by=posID][V1==2,posID]
length(unique(SNPs_with_more_than_one_introgressed_allele))
# 3838

dir.create(sprintf("%s/figures/%s/Z_figures/SupTable_1", MPRA_DIR, ANALYSIS_DIR))
fwrite(tested_and_ctrl_oligos_final_annot, file = sprintf("%s/figures/%s/Z_figures/SupTable_1/SupTable1b_tested_oligos.tsv", MPRA_DIR, ANALYSIS_DIR),sep='\t')

