MPRA_DIR <- "/pasteur/helix/projects/evo_immuno_pop/MPRA"
ANALYSIS_DIR <- "MPRA_count_exp6_analysis3"

source(sprintf("%s/scripts/%s/00_starter.R", MPRA_DIR, ANALYSIS_DIR))
source(sprintf("%s/scripts/%s/00_load_associations.R", MPRA_DIR, ANALYSIS_DIR))


RUN_ID <- "RUN3_Z2_nBC10"
#CRITERIA_ACTIVE <- "FDR5_scrambled5pct_FC0.2_GCnorm"
CRITERIA_ACTIVE <- "FDR5_FC1_FC0.2_GCnorm"

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

source(sprintf("%s/scripts/%s/02z__define_includedCRS.R", MPRA_DIR, ANALYSIS_DIR))

IN_DIR <- sprintf("%s/data/%s/03a_aggMPRA_analyse/%s/", MPRA_DIR, ANALYSIS_DIR, RUN_ID)
ACTIVE_DIR <- sprintf("%s/data/%s/03b_computeFDR/%s/Activity/%s/", MPRA_DIR, ANALYSIS_DIR, RUN_ID, CRITERIA_ACTIVE)
dir.create(ACTIVE_DIR, recursive = TRUE)
FIGURE_DIR <- sprintf("%s/figures/%s/%s", MPRA_DIR, ANALYSIS_DIR, RUN_ID)

ACTIVITY_DIR <- sprintf("%s/02_oligo_activity/%s", FIGURE_DIR, CRITERIA_ACTIVE)
dir.create(sprintf("%s/SupTables", ACTIVITY_DIR), showWarnings = FALSE, recursive = TRUE)
# dir.create(sprintf("%s/figures/%s/Z_figures/SupTable_1", ACTIVITY_DIR))

cat("Numbers , oligo activity")
# oligo_activity_perm <- fread(sprintf("%s/all_oligos_annotated_perm__%s.tsv.gz", ACTIVE_DIR, CRITERIA_ACTIVE))
# oligo_activity_perm_ctc <- oligo_activity_perm[power == 0 & boot == 0 & ANALYSIS_SUBTYPE == "celltype_cond", ]

oligo_activity_obs <- fread(sprintf("%s/all_oligos_annotated__%s.tsv.gz", ACTIVE_DIR, CRITERIA_ACTIVE))
oligo_activity_obs_ctc <- oligo_activity_obs[power == 0 & boot == 0 & ANALYSIS_SUBTYPE == "celltype_cond", ]
oligo_activity_obs_ctc_archaic <- oligo_activity_obs_ctc[oligo %in% tested_and_ctrl_oligos_final[type == "tested", oligo], ]

oligo_activity_obs_sample <- oligo_activity_obs[power == 0 & boot == 0 & ANALYSIS_SUBTYPE == "sample", ]
oligo_activity_obs_sample[, sample := gsub("-ACE2", "", gsub("-P3", "", ANALYSIS_NAME))]
oligo_activity_obs_sample <- merge(oligo_activity_obs_sample, technical_covariates[, .(SETUP_ID, plasmid_library, lentivirus_preparation)], by.x = "ANALYSIS_NAME", by.y = "SETUP_ID")

# ## test SNPS + ctrls
oligo_activity_obs_ctc[, length(unique(crsID))] # 3948
# oligo_activity_obs_ctc[ FDR <= .01, length(unique(crsID))] # 2028
# oligo_activity_obs_ctc[ FDR <= .05, length(unique(crsID))] # 2825
# oligo_activity_obs_ctc[ P_scrambled <= .05, length(unique(crsID))] # 891
# oligo_activity_obs_ctc[ P_scrambled <= .01, length(unique(crsID))] # 397
# oligo_activity_obs_ctc[ FDR_scrambled <= .05, length(unique(crsID))] # 255
# oligo_activity_obs_ctc[ lfdr <= .05, length(unique(crsID))] # 2197
# oligo_activity_obs_ctc[ lfdr <= .2, length(unique(crsID))] # 2906
# oligo_activity_obs_ctc[ lfdr <= .2 & P_scrambled <= .01, length(unique(crsID))] # 386
oligo_activity_obs_ctc[oligo_active_CRITERIA == TRUE, length(unique(crsID))]
oligo_activity_obs_ctc[oligo_strong_CRITERIA == TRUE, length(unique(crsID))]

# ## test SNPs only
oligo_activity_obs_ctc_archaic[, length(unique(crsID))] # 3837
# oligo_activity_obs_ctc_archaic[ FDR <= .01, length(unique(crsID))] # 1923
# oligo_activity_obs_ctc_archaic[ FDR <= .05, length(unique(crsID))] # 2718
# oligo_activity_obs_ctc_archaic[ P_scrambled <= .05, length(unique(crsID))] # 793
# oligo_activity_obs_ctc_archaic[ P_scrambled <= .01, length(unique(crsID))] # 306
# oligo_activity_obs_ctc_archaic[ FDR_scrambled <= .05, length(unique(crsID))] # 171
# oligo_activity_obs_ctc_archaic[ lfdr <= .05, length(unique(crsID))] # 2091
# oligo_activity_obs_ctc_archaic[ lfdr <= .2, length(unique(crsID))] # 2797
# oligo_activity_obs_ctc_archaic[ lfdr <= .2 & P_scrambled <= .01, length(unique(crsID))] # 295
oligo_activity_obs_ctc_archaic[oligo_active_CRITERIA == TRUE, length(unique(crsID))]
oligo_activity_obs_ctc_archaic[oligo_strong_CRITERIA == TRUE, length(unique(crsID))]


########## number of active CREs
active_oligo_counts <- oligo_activity_obs_ctc_archaic[oligo_active_CRITERIA == TRUE,
  .(
    nActive_CRE = length(unique(crsID)),
    nActive_oligo = length(unique(oligo)),
    nEnhancer_oligo = length(unique(oligo[alpha_CRITERIA > 1])),
    nSilencer_oligo = length(unique(oligo[alpha_CRITERIA < 1])),
    nEnhancer_CRE = length(unique(crsID[alpha_CRITERIA > 1])),
    nEnhancerOnly_CRE = length(setdiff(crsID[alpha_CRITERIA > 1], crsID[alpha_CRITERIA < 1])),
    nboth_CRE = length(intersect(crsID[alpha_CRITERIA > 1], crsID[alpha_CRITERIA < 1])),
    nSilencerOnly_CRE = length(setdiff(crsID[alpha_CRITERIA < 1], crsID[alpha_CRITERIA > 1]))
  ),
  by = COND_ID
]
active_oligo_counts
#          COND_ID nActive_CRE nActive_oligo nEnhancer_oligo nSilencer_oligo nEnhancer_CRE nEnhancerOnly_CRE nboth_CRE nSilencerOnly_CRE
#           <char>       <int>         <int>           <int>           <int>         <int>             <int>     <int>             <int>
#  1:     A549_IAV         861          1394             806             588           486               486         0               375
#  2:      A549_NS        1125          1873             996             877           585               584         1               540
#  3:    A549_SARS         728          1160             697             463           423               423         0               305
#  4:    A549_TNFa        1162          1927            1069             858           628               626         2               534
#  5:    HepG2_DEX        1960          3370            1749            1621          1002               999         3               958
#  6: HepG2_IFNA2b        1929          3357            1736            1621           980               978         2               949
#  7:     HepG2_NS        2005          3469            1809            1660          1032              1028         4               973
#  8:   HepG2_TNFa        1571          2648            1431            1217           831               829         2               740
#  9:     K562_DEX        1026          1690             864             826           515               515         0               511
# 10:  K562_IFNA2b         984          1647             873             774           508               508         0               476
# 11:      K562_NS        1063          1802             921             881           533               532         1               530
# 12:    K562_TNFa        1197          2012             960            1052           570               568         2               627

active_oligo_counts[, range(nActive_CRE)]
# 728 2005


# dir.create(sprintf("%s/figures/%s/Z_figures/SupTable_1", MPRA_DIR, ANALYSIS_DIR))
########### Supplementary table 1d
supTable1d_long <- oligo_activity_obs_sample[, .(chrom, start, type_simple, allele.label,sample, crsID, oligo, log2_activity = log2(alpha_CRITERIA), log2_activity.se = log_alpha_se_CRITERIA / log(2))]
cols <- c("log2_activity", "log2_activity.se")

supTable1d <- dcast(supTable1d_long, chrom + start + type_simple + allele.label + crsID + oligo ~ sample, value.var = cols)
my_levels=gsub('-ACE2|-P3', '', condition_summary_reps$analysis_name)
new_col_order <- CJ(my_levels, cols)
oo=order(factor(new_col_order$my_levels, levels = my_levels))
new_col_order <- new_col_order[oo,paste(cols, my_levels, sep = "_")]
setcolorder(supTable1d, c(setdiff(names(supTable1d), new_col_order), new_col_order))
fwrite(supTable1d[order(ifelse(type_simple%in%c('Promoter', 'scrambled'),'ctrl','tested'),chrom,start,allele.label),], file = sprintf("%s/SupTables/SupTable1d_oligos_activity_Sample__%s.tsv", ACTIVITY_DIR, CRITERIA_ACTIVE), sep = "\t", na = "NA")

########### Supplementary table 1d
active_CRS_list <- oligo_activity_obs_ctc_archaic[oligo_active_CRITERIA == TRUE, .(celline, condition, crsID, oligo, log2_activity = log2(alpha_CRITERIA), log2_activity.se = log_alpha_se_CRITERIA / log(2), activity_Z_score = Zscore_CRITERIA, pval_emp, FDR, P_scrambled, ChromHMM_cell_line = ROADMAP_current_celline, oligo_class = oligo_class_CRITERIA)]
fwrite(active_CRS_list, file = sprintf("%s/SupTables/SupTable1f_active_oligos__%s.tsv", ACTIVITY_DIR, CRITERIA_ACTIVE), sep = "\t")

CRS_activity_list <- oligo_activity_obs_ctc_archaic[, .(celline, condition, crsID, oligo, log2_activity = log2(alpha_CRITERIA), log2_activity.se = log_alpha_se_CRITERIA, activity_Z_score = Zscore_CRITERIA, lfdr, FDR, P_scrambled, ChromHMM_cell_line = ROADMAP_current_celline, oligo_class = oligo_class_CRITERIA)]
CRS_activity_list <- merge(CRS_activity_list, oligo_source[, .(crsID = posID, oligo, allele.label)], by = c("crsID", "oligo"))
CRS_activity_list_wide <- dcast(CRS_activity_list, crsID + celline + ChromHMM_cell_line + condition ~ allele.label, value.var = c("log2_activity", "log2_activity.se", "activity_Z_score", "lfdr", "P_scrambled", "oligo_class"))

active_CRE_wide <- merge(CRS_activity_list_wide, active_CRS_list[, .(celline, condition, crsID)], by = c("celline", "condition", "crsID"))

active_CRE_wide[, P_non_introgressed := 2 * pnorm(abs(`activity_Z_score_NON-INTROGRESSED`), low = F)]
active_CRE_wide[, P_introgressed := 2 * pnorm(abs(`activity_Z_score_INTROGRESSED`), low = F)]
active_CRE_wide[, mean(sign(activity_Z_score_INTROGRESSED) == sign(`activity_Z_score_NON-INTROGRESSED`) & P_non_introgressed < 0.01 & P_introgressed < 0.01, na.rm = T)]
# 0.9156448
# At a 20% local FDR, we identified 940-2,034 active CRE per condition of which ~90% displayed consistent effects on expression for both alleles (same direction of effect, replication P-value < 0.01).
# At a 5% local FDR, we identified 728-2,005 active CRE per condition of which >90% displayed consistent effects on expression for both alleles (same direction of effect, replication P-value < 0.01).

active_CRS_list[, .(Pct_enhancer_among_active = mean(grepl("enhancer", oligo_class))), by = .(celline, condition)]
#     celline condition Pct_enhancer_among_active
#      <char>    <char>                     <num>
#  1:    A549       IAV                 0.5781923
#  2:    A549        NS                 0.5317672
#  3:    A549      SARS                 0.6008621
#  4:    A549      TNFa                 0.5547483
#  5:   HepG2       DEX                 0.5189911
#  6:   HepG2    IFNA2b                 0.5171284
#  7:   HepG2        NS                 0.5214759
#  8:   HepG2      TNFa                 0.5404079
#  9:    K562       DEX                 0.5112426
# 10:    K562    IFNA2b                 0.5300546
# 11:    K562        NS                 0.5110988
# 12:    K562      TNFa                 0.4771372

context_specificity_enhancer <- active_CRS_list[grepl("enhancer", oligo_class), .(nContext = length(unique(paste(celline, condition)))), by = crsID]
context_specificity_enhancer[, .N, keyby = nContext]
#     nContext     N
#        <int> <int>
#  1:        1   183
#  2:        2   120
#  3:        3   128
#  4:        4   239
#  5:        5    70
#  6:        6    59
#  7:        7    67
#  8:        8   123
#  9:        9    68
# 10:       10    37
# 11:       11    49
# 12:       12   221

context_specificity_enhancer[, .N, keyby = nContext == 1][, .(context_specific = nContext, N, Pct = N / sum(N))]
#    context_specific     N       Pct
#              <lgcl> <int>     <num>
# 1:            FALSE  1181 0.8658358
# 2:             TRUE   183 0.1341642


context_specificity_silencer <- active_CRS_list[grepl("silencer", oligo_class), .(nContext = length(unique(paste(celline, condition)))), by = crsID]
context_specificity_silencer[, .N, keyby = nContext]

context_specificity_silencer[, .N, keyby = nContext == 1][, .(context_specific = nContext, N, Pct = N / sum(N))]
#    context_specific     N       Pct
# 1:            FALSE  1204 0.8449123
# 2:             TRUE   221 0.1550877

active_CRS_list[, .(is_enhancer = any(grepl("enhancer", oligo_class)), is_silencer = any(grepl("silencer", oligo_class))), by = crsID][, .N, by = .(is_enhancer, is_silencer)][, .(N, Pct = N / sum(N), is_enhancer, is_silencer)]
#        N        Pct is_enhancer is_silencer
#    <int>      <num>      <lgcl>      <lgcl>
# 1:  1250 0.46728972        TRUE       FALSE
# 2:  1311 0.49009346       FALSE        TRUE
# 3:   114 0.04261682        TRUE        TRUE

active_scrambled <- oligo_activity_obs_ctc[grepl("scrambled", oligo),
  .(
    Pct_active = mean(oligo_active_CRITERIA),
    Pct_enhancer = mean(grepl("enhancer", oligo_class_CRITERIA)),
    Pct_silencer = mean(grepl("silencer", oligo_class_CRITERIA)),
    Pct_enhancer_strong = mean(grepl("strong enhancer", oligo_class_CRITERIA)),
    Pct_silencer_strong = mean(grepl("strong silencer", oligo_class_CRITERIA))
  ),
  by = .(celline, condition)
]


active_scrambled
#     celline condition Pct_active Pct_enhancer Pct_silencer Pct_enhancer_strong Pct_silencer_strong
#      <char>    <char>      <num>        <num>        <num>               <num>               <num>
#  1:    A549       IAV  0.2105263    0.1473684   0.06315789          0.09473684          0.02105263
#  2:    A549        NS  0.2631579    0.1684211   0.09473684          0.09473684          0.02105263
#  3:    A549      SARS  0.2210526    0.1368421   0.08421053          0.09473684          0.03157895
#  4:    A549      TNFa  0.3052632    0.2000000   0.10526316          0.09473684          0.02105263
#  5:   HepG2       DEX  0.5263158    0.3263158   0.20000000          0.08421053          0.00000000
#  6:   HepG2    IFNA2b  0.5052632    0.2631579   0.24210526          0.10526316          0.02105263
#  7:   HepG2        NS  0.5578947    0.3157895   0.24210526          0.10526316          0.02105263
#  8:   HepG2      TNFa  0.4315789    0.2526316   0.17894737          0.13684211          0.02105263
#  9:    K562       DEX  0.3578947    0.1684211   0.18947368          0.14736842          0.05263158
# 10:    K562    IFNA2b  0.3263158    0.1684211   0.15789474          0.13684211          0.04210526
# 11:    K562        NS  0.5052632    0.2315789   0.27368421          0.14736842          0.02105263
# 12:    K562      TNFa  0.4526316    0.1894737   0.26315789          0.13684211          0.10526316

oligo_activity_obs_ctc[grepl("scrambled", oligo), .(active_any = any(oligo_active_CRITERIA)), by = crsID][, mean(active_any)]
# 0.7684211
active_scrambled[, mean(Pct_active)]
# [1] 0.3885965


active_count_archaic <- oligo_activity_obs_ctc_archaic[,
  .(
    is_active = any(oligo_active_CRITERIA),
    is_enhancer = any(grepl("enhancer", oligo_class_CRITERIA)),
    is_silencer = any(grepl("silencer", oligo_class_CRITERIA)),
    is_enhancer_strong = any(grepl("strong enhancer", oligo_class_CRITERIA)),
    is_silencer_strong = any(grepl("strong silencer", oligo_class_CRITERIA))
  ),
  by = .(crsID, celline, condition)
]

active_Pct_archaic <- active_count_archaic[,
  .(
    Pct_active = mean(is_active),
    Pct_enhancer = mean(is_enhancer),
    Pct_silencer = mean(is_silencer),
    Pct_enhancer_strong = mean(is_enhancer_strong),
    Pct_silencer_strong = mean(is_silencer_strong)
  ),
  by = .(celline, condition)
]
active_Pct_archaic
## scrambled 5%
#     celline condition Pct_active Pct_enhancer Pct_silencer Pct_enhancer_strong Pct_silencer_strong
#  1:    A549       IAV  0.2253927    0.1272251   0.09816754          0.04554974         0.007853403
#  2:    A549        NS  0.2944255    0.1531013   0.14158597          0.05260403         0.012300445
#  3:    A549      SARS  0.1905260    0.1107040   0.07982204          0.02538602         0.001308558
#  4:    A549      TNFa  0.3041885    0.1643979   0.14031414          0.04345550         0.005235602
#  5:   HepG2       DEX  0.5106826    0.2610735   0.25039083          0.03413236         0.001563314
#  6:   HepG2    IFNA2b  0.5026055    0.2553413   0.24778530          0.05054716         0.005992705
#  7:   HepG2        NS  0.5224075    0.2688900   0.25455967          0.03986451         0.002866076
#  8:   HepG2      TNFa  0.4112565    0.2175393   0.19424084          0.04973822         0.006544503
#  9:    K562       DEX  0.2685161    0.1347815   0.13373462          0.04815493         0.008636483
# 10:    K562    IFNA2b  0.2573222    0.1328452   0.12447699          0.04811715         0.006537657
# 11:    K562        NS  0.2769672    0.1388744   0.13835331          0.04663887         0.004950495
# 12:    K562      TNFa  0.3133508    0.1492147   0.16465969          0.08193717         0.040052356

## log2FC>1
#     celline condition Pct_active Pct_enhancer Pct_silencer Pct_enhancer_strong Pct_silencer_strong
#      <char>    <char>      <num>        <num>        <num>               <num>               <num>
#  1:    A549       IAV  0.2253927    0.1272251   0.09816754          0.01675393        0.0005235602
#  2:    A549        NS  0.2944255    0.1531013   0.14158597          0.01831981        0.0010468464
#  3:    A549      SARS  0.1905260    0.1107040   0.07982204          0.01648783        0.0000000000
#  4:    A549      TNFa  0.3041885    0.1643979   0.14031414          0.02172775        0.0013089005
#  5:   HepG2       DEX  0.5106826    0.2610735   0.25039083          0.02892131        0.0010422095
#  6:   HepG2    IFNA2b  0.5026055    0.2553413   0.24778530          0.03361126        0.0013027619
#  7:   HepG2        NS  0.5224075    0.2688900   0.25455967          0.03178739        0.0013027619
#  8:   HepG2      TNFa  0.4112565    0.2175393   0.19424084          0.02958115        0.0010471204
#  9:    K562       DEX  0.2685161    0.1347815   0.13373462          0.03978016        0.0049725203
# 10:    K562    IFNA2b  0.2573222    0.1328452   0.12447699          0.04053347        0.0033995816
# 11:    K562        NS  0.2769672    0.1388744   0.13835331          0.04689943        0.0049504950
# 12:    K562      TNFa  0.3133508    0.1492147   0.16465969          0.03664921        0.0028795812

######## outdated
# active_Pct_archaic <- oligo_activity_obs_ctc_archaic[,
#   .(
#     Pct_active = mean(oligo_active_CRITERIA),
#     Pct_enhancer = mean(grepl("enhancer", oligo_class_CRITERIA)),
#     Pct_silencer = mean(grepl("silencer", oligo_class_CRITERIA)),
#     Pct_enhancer_strong = mean(grepl("strong enhancer", oligo_class_CRITERIA)),
#     Pct_silencer_strong = mean(grepl("strong silencer", oligo_class_CRITERIA))
#   ),
#   by = .(celline, condition)
# ]
#     celline condition Pct_active Pct_enhancer Pct_silencer Pct_enhancer_strong Pct_silencer_strong
#      <char>    <char>      <num>        <num>        <num>               <num>               <num>
#  1:    A549       IAV  0.1827717   0.10567720   0.07709453          0.03684280        0.0045889603
#  2:    A549        NS  0.2456071   0.13060582   0.11500131          0.04235510        0.0081300813
#  3:    A549      SARS  0.1520913   0.09138587   0.06070539          0.02019143        0.0006555658
#  4:    A549      TNFa  0.2527213   0.14019672   0.11252459          0.03645902        0.0032786885
#  5:   HepG2       DEX  0.4391452   0.22791243   0.21123273          0.03036226        0.0010424811
#  6:   HepG2    IFNA2b  0.4380791   0.22654313   0.21153595          0.04332507        0.0045674018
#  7:   HepG2        NS  0.4519281   0.23566962   0.21625847          0.03621678        0.0022146952
#  8:   HepG2      TNFa  0.3473242   0.18769675   0.15962749          0.04197272        0.0048530955
#  9:    K562       DEX  0.2215522   0.11326691   0.10828526          0.04116413        0.0058993183
# 10:    K562    IFNA2b  0.2156041   0.11428197   0.10132216          0.04149758        0.0044508444
# 11:    K562        NS  0.2347883   0.12000000   0.11478827          0.03921824        0.0035179153
# 12:    K562      TNFa  0.2639035   0.12591815   0.13798531          0.06899265        0.0301678909

######## end outdated
active_Pct_archaic[, range(Pct_enhancer_strong + Pct_silencer_strong)]
## scrambled 5%
# [1]  0.02669458 0.12198953
# log2FC>1
# [1] 0.01648783 0.05184992



	active_Pct_archaic[, mean(Pct_enhancer_strong / (Pct_enhancer_strong + Pct_silencer_strong))]
# 0.8731469
active_Pct_archaic[, range(Pct_enhancer_strong / (Pct_enhancer_strong + Pct_silencer_strong))]

active_CRS_list[, .(is_strong_enhancer = any(grepl("strong enhancer", oligo_class)), is_strong_silencer = any(grepl("strong silencer", oligo_class))), by = crsID][, .N, by = .(is_strong_enhancer, is_strong_silencer)][, .(N, Pct = N / sum(N), is_strong_enhancer, is_strong_silencer)]

# 1:   509 0.1902803738               TRUE              FALSE
# 2:  1977 0.7390654206              FALSE              FALSE
# 3:   187 0.0699065421              FALSE               TRUE
# 4:     2 0.0007476636               TRUE               TRUE

######## @ Venn diagrams
active_CRS_list[condition == "NS" & grepl("strong enhancer", oligo_class), .(cellines = paste(sort(unique(celline)), collapse = " ")), by = .(crsID)][, .N, by = cellines]
#         cellines     N
#             <char> <int>
# 1: A549 HepG2 K562    56
# 2:       A549 K562    26
# 3:      A549 HepG2    35
# 4:            A549    84
# 5:      HepG2 K562    22
# 6:           HepG2    40
# 7:            K562    75

active_CRS_list[condition == "NS" & grepl("strong silencer", oligo_class), .(cellines = paste(sort(unique(celline)), collapse = " ")), by = .(crsID)][, .N, by = cellines]
#           cellines     N
#             <char> <int>
# 1:       A549 K562     5
# 2:            A549    33
# 3: A549 HepG2 K562     3
# 4:      A549 HepG2     6
# 5:           HepG2     2
# 6:            K562    11


active_CRS_list[condition == "NS" & grepl("enhancer", oligo_class), .(cellines = paste(sort(unique(celline)), collapse = " ")), by = .(crsID)][, .N, by = cellines]
#           cellines     N
#             <char> <int>
# 1: A549 HepG2 K562   323
# 2:      A549 HepG2   221
# 3:            A549    35
# 4:       A549 K562     6
# 5:      HepG2 K562   111
# 6:           HepG2   377
# 7:            K562    93
active_CRS_list <- oligo_activity_obs_ctc_archaic[oligo_active_CRITERIA == TRUE, .(celline, condition, crsID, oligo, log2_activity = log2(alpha_CRITERIA), log2_activity.se = log_alpha_se_CRITERIA / log(2), activity_Z_score = Zscore_CRITERIA, lfdr, FDR, P_scrambled, ChromHMM_cell_line = ROADMAP_current_celline, oligo_class = oligo_class_CRITERIA)]

active_CRS_list[condition == "NS" & grepl("silencer", oligo_class), .(cellines = paste(sort(unique(celline)), collapse = " ")), by = .(crsID)][, .N, by = cellines]
#           cellines     N
#             <char> <int>
# 1: A549 HepG2 K562   270
# 2:       A549 K562    29
# 3:            A549    33
# 4:      A549 HepG2   209
# 5:      HepG2 K562    93
# 6:           HepG2   405
# 7:            K562   139

# numbers for Fig 2B to redo separately for enhancers and silencers.
enhancer_count <- oligo_activity_obs_ctc_archaic[condition == "NS" & oligo_active_CRITERIA == TRUE & grepl("enhancer", oligo_class_CRITERIA), .(celline = paste(sort(unique(celline)), collapse = "/"), N_celline = length(unique(celline))), by = posID][, .(nEnhancer = length(unique(posID))), keyby = .(N_celline, celline)]
silencer_count <- oligo_activity_obs_ctc_archaic[condition == "NS" & oligo_active_CRITERIA == TRUE & grepl("silencer", oligo_class_CRITERIA), .(celline = paste(sort(unique(celline)), collapse = "/"), N_celline = length(unique(celline))), by = posID][, .(nSilencer = length(unique(posID))), keyby = .(N_celline, celline)]
strong_enhancer_count <- oligo_activity_obs_ctc_archaic[condition == "NS" & oligo_active_CRITERIA == TRUE & grepl("strong enhancer", oligo_class_CRITERIA), .(celline = paste(sort(unique(celline)), collapse = "/"), N_celline = length(unique(celline))), by = posID][, .(nStrongEnhancer = length(unique(posID))), keyby = .(N_celline, celline)]
strong_silencer_count <- oligo_activity_obs_ctc_archaic[condition == "NS" & oligo_active_CRITERIA == TRUE & grepl("strong silencer", oligo_class_CRITERIA), .(celline = paste(sort(unique(celline)), collapse = "/"), N_celline = length(unique(celline))), by = posID][, .(nStrongSilencer = length(unique(posID))), keyby = .(N_celline, celline)]
all_counts <- merge(enhancer_count, silencer_count, by = c("N_celline", "celline"), all = T)
all_counts_strong <- merge(strong_enhancer_count, strong_silencer_count, by = c("N_celline", "celline"), all = T)
all_counts <- merge(all_counts, all_counts_strong, by = c("N_celline", "celline"), all = T)
all_counts

#    N_celline         celline nEnhancer nSilencer nStrongEnhancer nStrongSilencer
#        <int>          <char>     <int>     <int>           <int>           <int>
# 1:         1            A549        35        33              84              33
# 2:         1           HepG2       377       405              40               2
# 3:         1            K562        93       139              75              11
# 4:         2      A549/HepG2       221       209              35               6
# 5:         2       A549/K562         6        29              26               5
# 6:         2      HepG2/K562       111        93              22              NA
# 7:         3 A549/HepG2/K562       323       270              56               3


melt(all_counts, id.vars = c("N_celline", "celline"))[, .(
  Pct_celline_specific = sum(value[N_celline == 1], na.rm = T) / sum(value, na.rm = T),
  Pct_shared2ormore = sum(value[N_celline > 1], na.rm = T) / sum(value, na.rm = T),
  Pct_shared3 = sum(value[N_celline == 3], na.rm = T) / sum(value, na.rm = T)
), by = variable]
#           variable Pct_celline_specific Pct_shared2ormore Pct_shared3
#             <fctr>                <num>             <num>       <num>
# 1:       nEnhancer            0.4331046         0.5668954   0.2770154
# 2:       nSilencer            0.4898132         0.5101868   0.2292020
# 3: nStrongEnhancer            0.5887574         0.4112426   0.1656805
# 4: nStrongSilencer            0.7666667         0.2333333   0.0500000

cat('##### is CRE activity shared ?: ')
CRE_sharing=oligo_activity_obs_ctc_archaic[condition == "NS" & oligo_active_CRITERIA == TRUE, .(celline_enhancer = paste(sort(unique(celline[grepl("enhancer", oligo_class_CRITERIA)])), collapse = "/"), 
																																										celline_silencer=paste(sort(unique(celline[grepl("silencer", oligo_class_CRITERIA)])), collapse = "/"), 
																																										N_celline_enhancer= length(unique(celline[grepl("enhancer", oligo_class_CRITERIA)])), 
																																										N_celline_silencer= length(unique(celline[grepl("silencer", oligo_class_CRITERIA)])),
																																										N_celline_strong_enhancer= length(unique(celline[grepl("strong enhancer", oligo_class_CRITERIA)])), 
																																										N_celline_strong_silencer= length(unique(celline[grepl("strong silencer", oligo_class_CRITERIA)])) 

																																										), by = posID]
CRE_sharing[,mean(N_celline_enhancer>=2 | N_celline_silencer>=2)]
# 0.5522976
CRE_sharing[(N_celline_strong_enhancer+N_celline_strong_silencer)>0,mean(N_celline_strong_enhancer>=2 | N_celline_strong_silencer>=2)]
## scambled 
# 0.3844221
## log2FC>1 
# 0.3296296
# Although CRE activity was often shared (55%; Fig. 2b), strong CREs showed greater specificity, with only 33 % active in two or more cell lines (Supplementary Fig. 4). 

cat("\ncount enhancers, discovery\n")
all_counts[, .(A549 = sum(nEnhancer * grepl("A549", celline)), HepG2 = sum(nEnhancer * grepl("HepG2", celline)), K562 = sum(nEnhancer * grepl("K562", celline)))]
#    A549 HepG2  K562
#    <int> <int> <int>
# 1:   585  1032   533
cat("\ncount silencers, discovery\n")
all_counts[, .(A549 = sum(nSilencer * grepl("A549", celline)), HepG2 = sum(nSilencer * grepl("HepG2", celline)), K562 = sum(nSilencer * grepl("K562", celline)))]
#     A549 HepG2  K562
#    <int> <int> <int>
# 1:   541   977   531

# numbers for Fig 2B (with marginal Pvalue<0.05 for replication).
enhancer_NS_dicovery <- oligo_activity_obs_ctc_archaic[condition == "NS" & oligo_active_CRITERIA == TRUE & grepl("enhancer", oligo_class_CRITERIA), unique(crsID)]
silencer_NS_dicovery <- oligo_activity_obs_ctc_archaic[condition == "NS" & oligo_active_CRITERIA == TRUE & grepl("silencer", oligo_class_CRITERIA), unique(crsID)]

enhancer_count <- oligo_activity_obs_ctc_archaic[condition == "NS" & crsID %chin% enhancer_NS_dicovery & pval_emp < 0.05 & log2(alpha_CRITERIA) > LOG2FC_TH, .(celline = paste(sort(unique(celline)), collapse = "/"), N_celline = length(unique(celline))), by = posID][, .(nEnhancer = length(unique(posID))), keyby = .(N_celline, celline)]
silencer_count <- oligo_activity_obs_ctc_archaic[condition == "NS" & crsID %chin% silencer_NS_dicovery & pval_emp < 0.05 & log2(alpha_CRITERIA) < -LOG2FC_TH, .(celline = paste(sort(unique(celline)), collapse = "/"), N_celline = length(unique(celline))), by = posID][, .(nSilencer = length(unique(posID))), keyby = .(N_celline, celline)]
strong_enhancer_count <- oligo_activity_obs_ctc_archaic[condition == "NS" & crsID %chin% enhancer_NS_dicovery & pval_emp < 0.05 & log2(alpha_CRITERIA) > LOG2FC_TH & P_scrambled < SCRAMBLED_TH, .(celline = paste(sort(unique(celline)), collapse = "/"), N_celline = length(unique(celline))), by = posID][, .(nStrongEnhancer = length(unique(posID))), keyby = .(N_celline, celline)]
strong_silencer_count <- oligo_activity_obs_ctc_archaic[condition == "NS" & crsID %chin% silencer_NS_dicovery & pval_emp < 0.05 & log2(alpha_CRITERIA) < -LOG2FC_TH & P_scrambled < SCRAMBLED_TH, .(celline = paste(sort(unique(celline)), collapse = "/"), N_celline = length(unique(celline))), by = posID][, .(nStrongSilencer = length(unique(posID))), keyby = .(N_celline, celline)]
all_counts <- merge(enhancer_count, silencer_count, by = c("N_celline", "celline"), all = T)
all_counts_strong <- merge(strong_enhancer_count, strong_silencer_count, by = c("N_celline", "celline"), all = T)
all_counts <- merge(all_counts, all_counts_strong, by = c("N_celline", "celline"), all = T)
all_counts

#    N_celline         celline nEnhancer nSilencer nStrongEnhancer nStrongSilencer
#        <int>          <char>     <int>     <int>           <int>           <int>
# 1:         1            A549        32        24              93              33
# 2:         1           HepG2       221       203              40               2
# 3:         1            K562        77       122              75              11
# 4:         2      A549/HepG2       241       270              34               6
# 5:         2       A549/K562        14        44              28               5
# 6:         2      HepG2/K562        96        89              20              NA
# 7:         3 A549/HepG2/K562       485       426              59               3

melt(all_counts, id.vars = c("N_celline", "celline"))[, .(
  Pct_celline_specific = sum(value[N_celline == 1], na.rm = T) / sum(value, na.rm = T),
  Pct_shared2ormore = sum(value[N_celline > 1], na.rm = T) / sum(value, na.rm = T),
  Pct_shared3 = sum(value[N_celline == 3], na.rm = T) / sum(value, na.rm = T)
), by = variable]

#             <fctr>                <num>             <num>       <num>
# 1:       nEnhancer            0.2830189         0.7169811   0.4159520
# 2:       nSilencer            0.2962649         0.7037351   0.3616299
# 3: nStrongEnhancer            0.5959885         0.4040115   0.1690544
# 4: nStrongSilencer            0.7666667         0.2333333   0.0500000


cat("\ncount enhancers, with lower P for replication\n")
all_counts[, .(A549 = sum(nEnhancer * grepl("A549", celline)), HepG2 = sum(nEnhancer * grepl("HepG2", celline)), K562 = sum(nEnhancer * grepl("K562", celline)))]
#     A549 HepG2  K562
#    <int> <int> <int>
# 1:   772  1043   672
cat("\ncount silencer, with lower P for replication\n")
all_counts[, .(A549 = sum(nSilencer * grepl("A549", celline)), HepG2 = sum(nSilencer * grepl("HepG2", celline)), K562 = sum(nSilencer * grepl("K562", celline)))]
#     A549 HepG2  K562
#    <int> <int> <int>
# 1:   764   988   681
cat("All done !")
q("no")


cre_A549_up_sugg= oligo_activity_obs_ctc_archaic[oligo%in%oligo_A549_up_sugg,unique(crsID)]
cre_HepG2_up_sugg= oligo_activity_obs_ctc_archaic[oligo%in%oligo_HepG2_up_sugg,unique(crsID)]
cre_K562_up_sugg= oligo_activity_obs_ctc_archaic[oligo%in%oligo_K562_up_sugg,unique(crsID)]

cre_A549_up= oligo_activity_obs_ctc_archaic[oligo%in%oligo_A549_up,unique(crsID)]
cre_HepG2_up= oligo_activity_obs_ctc_archaic[oligo%in%oligo_HepG2_up,unique(crsID)]
cre_K562_up= oligo_activity_obs_ctc_archaic[oligo%in%oligo_K562_up,unique(crsID)]

cre_A549_HepG2_diff = oligo_celltype_comp_NS_annot_obs_archaic[ANALYSIS_NAME == "Celltype_A549_HepG2_NS" & pval_emp < 0.05 , unique(crsID)]
cre_A549_K562_diff = oligo_celltype_comp_NS_annot_obs_archaic[ANALYSIS_NAME == "Celltype_A549_K562_NS" & pval_emp < 0.05 , unique(crsID)]
cre_HepG2_K562_diff = oligo_celltype_comp_NS_annot_obs_archaic[ANALYSIS_NAME == "Celltype_HepG2_K562_NS" & pval_emp < 0.05 , unique(crsID)]

CRE_sharing[,diff_sugg_A549_HepG2:=posID %in% cre_A549_HepG2_diff]
CRE_sharing[,diff_sugg_A549_K562:=posID %in% cre_A549_K562_diff]
CRE_sharing[,diff_sugg_HepG2_K562:=posID %in% cre_HepG2_K562_diff]

CRE_sharing[,cre_A549_up_sugg:=posID %in% cre_A549_up_sugg]
CRE_sharing[,cre_HepG2_up_sugg:=posID %in% cre_HepG2_up_sugg]
CRE_sharing[,cre_K562_up_sugg:=posID %in% cre_K562_up_sugg]

CRE_sharing[,cre_A549_up:=posID %in% cre_A549_up]
CRE_sharing[,cre_HepG2_up:=posID %in% cre_HepG2_up]
CRE_sharing[,cre_K562_up:=posID %in% cre_K562_up]



CRE_sharing[N_celline_enhancer>0 |N_celline_enhancer,fisher.test(table(ifelse(N_celline_strong_enhancer>0,'strong','other'), cre_A549_up_sugg | cre_HepG2_up_sugg | cre_K562_up_sugg))]
#### scrambled P<5%
# p-value = 1.08e-07
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  1.776133 3.830925
# sample estimates:
# odds ratio 
#   2.588444 

#### log2FC>1 
# p-value = 1.785e-06
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  1.509825 2.802375
# sample estimates:
# odds ratio 
#   2.050685 

CRE_sharing[N_celline_enhancer>0,fisher.test(table(ifelse(N_celline_strong_enhancer>0,'strong','other'), cre_A549_up | cre_HepG2_up | cre_K562_up))]
#### scrambled P<5%
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  2.91558 6.03342
# sample estimates:
# odds ratio 
#   4.174212 

#### log2FC>1 
# data:  table(ifelse(N_celline_strong_enhancer > 0, "strong", "other"), cre_A549_up | cre_HepG2_up | cre_K562_up)
# p-value = 1.11e-12
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  2.093114 3.818405
# sample estimates:
# odds ratio 
#   2.821544 
