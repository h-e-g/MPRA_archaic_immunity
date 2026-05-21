################################################################################
################################################################################
# File name: 00_load_associations.R
# Author: M. Rotival
################################################################################
################################################################################

# TODO: harmonize how low complexity barcodes are defined and how they are used
# TODO: check where we use Associations_raw and Associations_Filtered

# MPRA_DIR='/pasteur/helix/projects/evo_immuno_pop/users/mrotival/MPRA_Lisa'
# COUNT_DIR='MPRA_count_exp2_v2'
# source(sprintf('%s/scripts/%s/00_starter.R',COUNT_DIR,MPRA_DIR))


##### read association results (raw)
Associations_raw <- fread(OUT_LIB_0_1_2_RAW, header = TRUE)
# Associations_raw=Associations_raw[,.N,by=.(BC,CRS,LIB)]
# Associations_raw[,Pct:=N/sum(N),by=.(BC,LIB)]
# Associations_raw[,type:=case_when(N>=3 & Pct>.5 ~ '1_assoc',
#                                   N>=3 & Pct<=.5 ~ '2_ambiguous',
#                                   N<3~ '3_low_coverage')]
#
BC_group <- Associations_raw[order(LIB, BC, type), ][!duplicated(paste(LIB, BC)), ]
BC_group[, .N, keyby = .(type, LIB)][, .(type, N, Pct = N / sum(N)), by = LIB]
# type       N        Pct
# 1:        1_assoc 2936332 0.44322842
# 2:    2_ambiguous  145606 0.02197869
# 3: 3_low_coverage 3542935 0.53479289

# Associations_raw[,nT_BC:=str_count(BC,'T')]
# Associations_raw[,nA_BC:=str_count(BC,'A')]
# Associations_raw[,nG_BC:=str_count(BC,'G')]
# Associations_raw[,nC_BC:=str_count(BC,'C')]
# Associations_raw[,shannon:=-(nT_BC*log((1+nT_BC)/16)+nC_BC*log((1+nC_BC)/16)+nG_BC*log((1+nG_BC)/16)+nA_BC*log((1+nA_BC)/16))]
setnames(Associations_raw, "LIB", "barcode_library")
setnames(Associations_raw, "CRS", "oligo")
Associations_Filtered <- Associations_raw[N >= 3 & Pct > .5, ] # keep  barcode with > 3 reads and a dominant CRS (>50% of association reads)
#fwrite(Associations_Filtered, file = OUT_LIB_0_1_2_FILTERED, sep = "\t")

Associations_Filtered <- Associations_raw[N >= 3 & Pct == 1 & shannon > 11, ] # keep only barcode with > 3 reads, a unique CRS in the association reads. 
# In addition, to reduce the risk of missagned barcodes, we exclude barcodes with a sequence complexity <11 
# where sequence complexity is defined as the shannon entropy of the nucleotide frequencies within the barcode, 
# given by SHANNON = -(nT_BC*log((1+nT_BC)/16)+nC_BC*log((1+nC_BC)/16)+nG_BC*log((1+nG_BC)/16)+nA_BC*log((1+nA_BC)/16))
# The threshold of 11 was selected to exlcude less than 0.1% of truly random 15bp sequences while removing sequences with a strong overrepresentation of one nucleotide.

# fwrite(Associations_Filtered, file = OUT_LIB_0_1_2_FILTERED_STRICT, sep = "\t")
Associations_Filtered[,BC_num:=paste0('b',seq_len(.N)),by=.(oligo,barcode_library)]


##### read association results (filtered)
# Associations_Filtered=fread(OUT_THREELANES_ONEIDENTIFIER_filtered,header=T)

Associations_Filtered[,length(unique(BC)),by=.(oligo,barcode_library)]

################################################################################
#####     tests :                                                            ###
##### how many barcodes are lost if we are more stringent on barcode sharing ###
################################################################################

# for (Pct_target in c(0.5,0.6,0.7,0.8,0.9,0.95,0.99,1)){
#   cat(Pct_target,
#       as.numeric(Associations_raw[N>=3 & Pct>=Pct_target & shannon>11,.(.N,length(unique(BC)),length(unique(CRS)))]),
#       '\n',sep='\t')
# }
# Pct     N_assoc N_BC    N_CRS
# 0.5     2811378 2803854 11452
# 0.6     2756913 2756913 11448
# 0.7     2701572 2701572 11444
# 0.8     2620472 2620472 11439
# 0.9     2475126 2475126 11432
# 0.95    2410660 2410660 11428
# 0.99    2390440 2390440 11427
# 1       2389138 2389138 11427

################################################################################
#####     tests :                                                            ###
##### What complexity should we expect for truly random sequences            ###
################################################################################

# shannon_rand=c()
# for (i in 1:100000){
#   tab=table(cut(runif(15),seq(0,1,l=5),labels=LETTERS[1:4]))
#   shannon_rand[i]=-sum(tab*log((1+tab)/16))
# }
# range(shannon_rand)
# # [1]  6.65 18.11
# mean(shannon_rand<10)
# # [1] 4e-04
# mean(shannon_rand<11)
# # [1] 0.00093
# mean(shannon_rand<12)
# # [1] 0.00374
# mean(shannon_rand<13)
# # [1] 0.01114
# quantile(shannon_rand,c(1e-4,1e-3))
# #    0.01%      0.1%
# #    8.65      11.21



# shannon_rand=c()
# for (i in 1:100000){
#   tab=table(cut(runif(16),seq(0,1,l=5),labels=LETTERS[1:4]))
#   shannon_rand[i]=-sum(tab*log((1+tab)/17))
# }
# range(shannon_rand)
# # [1]  6.65 18.11
# mean(shannon_rand<10)
# # [1] 4e-04
# mean(shannon_rand<11)
# # [1] 0.00093
# mean(shannon_rand<12)
# # [1] 0.00374
# mean(shannon_rand<13)
# # [1] 0.01114
# quantile(shannon_rand,c(1e-4,1e-3))
#    0.01%      0.1%
#    8.65      11.21
