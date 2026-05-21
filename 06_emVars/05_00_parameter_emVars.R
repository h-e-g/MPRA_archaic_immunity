# 05_00_parameter_emvar.R

# current default
# CRITERIA_ACTIVE='lfdr20_GCnorm'
# CRITERIA_EMVARS='Diff_lfdr20_FC.2'

if (CRITERIA_EMVARS == "EmVar_lfdr20_FC.2") {
  ### threshold choice
  LFDR_TH_EMVAR <- 0.2
  FDR_TH_EMVAR <- Inf
  LOG2FC_TH_EMVAR <- 0.2
}

if (CRITERIA_EMVARS == "EmVar_lfdr10_FC.2") {
  ### threshold choice
  LFDR_TH_EMVAR <- 0.1
  FDR_TH_EMVAR <- Inf
  LOG2FC_TH_EMVAR <- 0.2
}

if (CRITERIA_EMVARS == "EmVar_FDR10_FC.2") {
  ### threshold choice
  LFDR_TH_EMVAR <- Inf
  FDR_TH_EMVAR <- 0.1
  LOG2FC_TH_EMVAR <- 0.2
}

if (CRITERIA_EMVARS == "EmVar_FDR5_FC.2") {
  ### threshold choice
  LFDR_TH_EMVAR <- Inf
  FDR_TH_EMVAR <- 0.05
  LOG2FC_TH_EMVAR <- -Inf
}


if (CRITERIA_EMVARS == "EmVar_FDR1_FC.2") {
  ### threshold choice
  LFDR_TH_EMVAR <- Inf
  FDR_TH_EMVAR <- 0.01
  LOG2FC_TH_EMVAR <- 0.2
}



###############################################################################################
########################################## functions ##########################################
###############################################################################################

########## compute percentage of emVars & related metrics within each group ##############
get_Pct_emVars <- function(emVar_table, split_by = NULL, total_by = NULL, total_method = c("recompute", "sum"),emVarlist=c('significant','suggestive')) {
  total_method <- match.arg(total_method)
  emVarlist <- match.arg(emVarlist)
  emVar_table <- emVar_table[seq_len(.N), ]

  if(emVarlist=='significant'){  
    emVar_table[, is_enhancer_CRITERIA := grepl("enhancer", cre_class_CRITERIA)]
    emVar_table[, is_silencer_CRITERIA := grepl("silencer", cre_class_CRITERIA)]
    # count alleles, avective CREs, emVars, emVars increasing activity of active CREs, emVars increasing expression.
    emVar_count <- emVar_table[, .(
      N_introgressed_allele = length(unique(ID)),
      N_CRE_overlaping_allele = length(unique(ID[is_tested_cre == TRUE])),
      N_activity_altering_allele = length(unique(ID[is_emVar_CRITERIA == TRUE])),
      N_activity_increasing_allele = length(unique(ID[ is_tested_cre == TRUE & is_emVar_CRITERIA == TRUE & ((is_enhancer_CRITERIA == TRUE & log2FC_archaic_vs_modern > 0) | (is_silencer_CRITERIA == TRUE & log2FC_archaic_vs_modern < 0))])),
      N_activity_decreasing_allele = length(unique(ID[ is_tested_cre == TRUE & is_emVar_CRITERIA == TRUE & ((is_enhancer_CRITERIA == TRUE & log2FC_archaic_vs_modern < 0) | (is_silencer_CRITERIA == TRUE & log2FC_archaic_vs_modern > 0))])),
      N_expression_increasing_allele = length(unique(ID[ is_tested_cre == TRUE & is_emVar_CRITERIA == TRUE & log2FC_archaic_vs_modern > 0])),
      N_expression_decreasing_allele = length(unique(ID[ is_tested_cre == TRUE & is_emVar_CRITERIA == TRUE & log2FC_archaic_vs_modern < 0]))
    ), keyby = split_by]
  }else{
    emVar_table[, is_enhancer_CRITERIA := (a1+a2)>median(a1+a2)]
    emVar_table[, is_silencer_CRITERIA := (a1+a2)<median(a1+a2)]
    # count alleles, avective CREs, emVars, emVars increasing activity of active CREs, emVars increasing expression.
    emVar_count <- emVar_table[, .(
      N_introgressed_allele = length(unique(ID)),
      N_CRE_overlaping_allele = length(unique(ID[is_tested_cre == TRUE])),
      N_activity_altering_allele = length(unique(ID[ is_tested_cre == TRUE & pval_emp<=0.05])),
      N_activity_increasing_allele = length(unique(ID[ is_tested_cre == TRUE & pval_emp<=0.05 & ((is_enhancer_CRITERIA == TRUE & log2FC_archaic_vs_modern > 0) | (is_silencer_CRITERIA == TRUE & log2FC_archaic_vs_modern < 0))])),
      N_activity_decreasing_allele = length(unique(ID[ is_tested_cre == TRUE & pval_emp<=0.05 & ((is_enhancer_CRITERIA == TRUE & log2FC_archaic_vs_modern < 0) | (is_silencer_CRITERIA == TRUE & log2FC_archaic_vs_modern > 0))])),
      N_expression_increasing_allele = length(unique(ID[ is_tested_cre == TRUE & pval_emp<=0.05 & log2FC_archaic_vs_modern > 0])),
      N_expression_decreasing_allele = length(unique(ID[ is_tested_cre == TRUE & pval_emp<=0.05 & log2FC_archaic_vs_modern < 0]))
    ), keyby = split_by]
  }
  
  if (total_method == "sum") {
    # sum count by specified variants (WARNING if there is overlap, the count will be double counted)
    emVar_count[, N_introgressed_allele_tot := sum(N_introgressed_allele), keyby = total_by]
    emVar_count[, N_CRE_overlaping_allele_tot := sum(N_CRE_overlaping_allele), keyby = total_by]
    emVar_count[, N_activity_altering_allele_tot := sum(N_activity_altering_allele), keyby = total_by]
    emVar_count[, N_activity_increasing_allele_tot := sum(N_activity_increasing_allele), keyby = total_by]
    emVar_count[, N_activity_decreasing_allele_tot := sum(N_activity_decreasing_allele), keyby = total_by]
    emVar_count[, N_expression_increasing_allele_tot := sum(N_expression_increasing_allele), keyby = total_by]
    emVar_count[, N_expression_decreasing_allele_tot := sum(N_expression_decreasing_allele), keyby = total_by]
  } else {
      if(emVarlist=='significant'){  
    # recompute count by specified variables
    emVar_count_tot <- emVar_table[, .(
      N_introgressed_allele_tot = length(unique(ID)),
      N_CRE_overlaping_allele_tot = length(unique(ID[is_tested_cre == TRUE])),
      N_activity_altering_allele_tot = length(unique(ID[is_emVar_CRITERIA == TRUE])),
      N_activity_increasing_allele_tot = length(unique(ID[is_emVar_CRITERIA == TRUE & is_tested_cre == TRUE & ((is_enhancer_CRITERIA == TRUE & log2FC_archaic_vs_modern > 0) | (is_silencer_CRITERIA == TRUE & log2FC_archaic_vs_modern < 0))])),
      N_activity_decreasing_allele_tot = length(unique(ID[is_emVar_CRITERIA == TRUE & is_tested_cre == TRUE & ((is_enhancer_CRITERIA == TRUE & log2FC_archaic_vs_modern < 0) | (is_silencer_CRITERIA == TRUE & log2FC_archaic_vs_modern > 0))])),
      N_expression_increasing_allele_tot = length(unique(ID[is_emVar_CRITERIA == TRUE & log2FC_archaic_vs_modern > 0])),
      N_expression_decreasing_allele_tot = length(unique(ID[is_emVar_CRITERIA == TRUE & log2FC_archaic_vs_modern < 0]))
    ), keyby = total_by]
      }else{
    emVar_count_tot <- emVar_table[, .(
      N_introgressed_allele_tot = length(unique(ID)),
      N_CRE_overlaping_allele_tot = length(unique(ID[is_tested_cre == TRUE])),
      N_activity_altering_allele_tot = length(unique(ID[ is_tested_cre == TRUE & pval_emp<=0.05])),
      N_activity_increasing_allele_tot = length(unique(ID[ is_tested_cre == TRUE & pval_emp<=0.05 & ((is_enhancer_CRITERIA == TRUE & log2FC_archaic_vs_modern > 0) | (is_silencer_CRITERIA == TRUE & log2FC_archaic_vs_modern < 0))])),
      N_activity_decreasing_allele_tot = length(unique(ID[ is_tested_cre == TRUE & pval_emp<=0.05 & ((is_enhancer_CRITERIA == TRUE & log2FC_archaic_vs_modern < 0) | (is_silencer_CRITERIA == TRUE & log2FC_archaic_vs_modern > 0))])),
      N_expression_increasing_allele_tot = length(unique(ID[ is_tested_cre == TRUE & pval_emp<=0.05 & log2FC_archaic_vs_modern > 0])),
      N_expression_decreasing_allele_tot = length(unique(ID[ is_tested_cre == TRUE & pval_emp<=0.05 & log2FC_archaic_vs_modern < 0]))
    ), keyby = total_by]
      }
    }
    shared_by <- intersect(total_by, split_by)
    if (length(shared_by) > 0) {
      emVar_count <- merge(emVar_count, emVar_count_tot, by = shared_by)
    } else {
      emVar_count <- cbind(emVar_count, emVar_count_tot)
    }

  get_Binom_stats <- function(x, N, p) {
    bin.result <- suppressWarnings(try(binom.test(x, N, p = p)))
    if (class(bin.result) == "try-error") {
      list(x, N, estimate = NA, NA, NA, 1)
    } else {
      list(x, N, bin.result$estimate, bin.result$conf.int[1], bin.result$conf.int[2], bin.result$p.value)
    }
  }
  emVar_count[, c("N_CREoverlap", "Ntot_CREoverlap", "Pct_CREoverlap", "Pct_CREoverlap_lo", "Pct_CREoverlap_hi", "Pval_CREoverlap") := get_Binom_stats(N_CRE_overlaping_allele, N_introgressed_allele, p = N_CRE_overlaping_allele_tot / N_introgressed_allele_tot), by = split_by] # N_CRE_overlaping_allele/N_introgressed_allele, by = split_by]
  emVar_count[, c("N_emVar", "Ntot_emVar", "Pct_emVar", "Pct_emVar_lo", "Pct_emVar_hi", "Pval_emVar") := get_Binom_stats(N_activity_altering_allele, N_introgressed_allele, p = N_activity_altering_allele_tot / N_introgressed_allele_tot), by = split_by] # N_CRE_overlaping_allele/N_introgressed_allele, by = split_by]
  emVar_count[, c("N_ActivityIncreasing", "Ntot_ActivityIncreasing", "Pct_ActivityIncreasing", "Pct_ActivityIncreasing_lo", "Pct_ActivityIncreasing_hi", "Pval_ActivityIncreasing") := get_Binom_stats(N_activity_increasing_allele, N_introgressed_allele, p = N_activity_increasing_allele_tot / N_introgressed_allele_tot), by = split_by] # N_CRE_overlaping_allele/N_introgressed_allele, by = split_by]
  emVar_count[, c("N_ActivityDecreasing", "Ntot_ActivityDecreasing", "Pct_ActivityDecreasing", "Pct_ActivityDecreasing_lo", "Pct_ActivityDecreasing_hi", "Pval_ActivityDecreasing") := get_Binom_stats(N_activity_decreasing_allele, N_introgressed_allele, p = N_activity_increasing_allele_tot / N_introgressed_allele_tot), by = split_by] # N_CRE_overlaping_allele/N_introgressed_allele, by = split_by]
  emVar_count[, c("N_ExpressionIncreasing", "Ntot_ExpressionIncreasing", "Pct_ExpressionIncreasing", "Pct_ExpressionIncreasing_lo", "Pct_ExpressionIncreasing_hi", "Pval_ExpressionIncreasing") := get_Binom_stats(N_expression_increasing_allele, N_introgressed_allele, p = N_expression_increasing_allele_tot / N_introgressed_allele_tot), by = split_by]
  emVar_count[, c("N_ExpressionDecreasing", "Ntot_ExpressionDecreasing", "Pct_ExpressionDecreasing", "Pct_ExpressionDecreasing_lo", "Pct_ExpressionDecreasing_hi", "Pval_ExpressionDecreasing") := get_Binom_stats(N_expression_decreasing_allele, N_introgressed_allele, p = N_expression_decreasing_allele_tot / N_introgressed_allele_tot), by = split_by]
  emVar_count[N_activity_altering_allele>0, c("N_sign", "Ntot_sign", "Pct_sign", "Pct_sign_lo", "Pct_sign_hi", "Pval_sign") := get_Binom_stats(N_activity_increasing_allele, N_activity_increasing_allele+N_activity_decreasing_allele, p = 0.5), by = split_by] # N_CRE_overlaping_allele/N_introgressed_allele, by = split_by]
  emVar_count[N_activity_altering_allele>0, c("N_sign_expression", "Ntot_sign_expression", "Pct_sign_expression", "Pct_sign_expression_lo", "Pct_sign_expression_hi", "Pval_sign_expression") := get_Binom_stats(N_expression_increasing_allele, N_activity_altering_allele, p = 0.5), by = split_by] # N_CRE_overlaping_allele/N_introgressed_allele, by = split_by]

  return(emVar_count[seq_len(.N), ])
}

########## do a pivot long on percentage of emVars & related metrics ##############
make_it_long <- function(Pct_emVars, split_by = NULL, measures = c("CREoverlap", "emVar", "ActivityIncreasing", "ActivityDecreasing", "ExpressionIncreasing", "ExpressionDecreasing", "sign", "sign_expression")) {
  Measures_vars <- unlist(lapply(measures, function(x) {
    sprintf(c("N_%s","Ntot_%s","Pct_%s", "Pct_%s_lo", "Pct_%s_hi", "Pval_%s"), x)
  }))
  Pct_emVars_long <- melt(Pct_emVars, id.vars = c(split_by))
  Pct_emVars_long <- Pct_emVars_long[as.character(variable) %chin% Measures_vars]
  Pct_emVars_long[, measure := gsub(sprintf("(Pct|Pval|N|Ntot)_(%s)(_lo|_hi)?", paste(measures, collapse = "|")), "\\2", as.character(variable))]
  Pct_emVars_long[, type := case_when(
    variable == sprintf("Pct_%s", measure) ~ "Pct",
    variable == sprintf("Pct_%s_lo", measure) ~ "Pct_lo",
    variable == sprintf("Pct_%s_hi", measure) ~ "Pct_hi",
    variable == sprintf("Pval_%s", measure) ~ "Pvalue",
	  variable == sprintf("N_%s", measure) ~ "N_obs",
    variable == sprintf("Ntot_%s", measure) ~ "N_test",
  	TRUE ~ "OTHER"
  ), by = .(measure)]
  Pct_emVars_wide <- tidyr::pivot_wider(Pct_emVars_long, id_cols = c(split_by, "measure"), names_from = type, values_from = value)
  Pct_emVars_wide <- as.data.table(Pct_emVars_wide)
  return(Pct_emVars_wide)
}

  chrom_length <- data.table(chromosome = c(1:22, "X", "Y"), length_GRch37 = c(249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431, 135534747, 135006516, 133851895, 115169878, 107349540, 102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566, 155270560, 59373566))
  chrom_length[, cumPos := cumsum(c(0, length_GRch37[-length(length_GRch37)] + 1e6))]
  chrom_length[, chrom_start := cumPos]
  chrom_length[, chrom_end := cumPos + length_GRch37]
  chrom_length[, chrom_mid := cumPos + length_GRch37/2]

########## add x-axis coordinates for manhattan plots ##############
add_manhattan_coord <- function(Pct_table) {
  locus_regexpr <- "([0-9XY]+):([0-9]+)-([0-9]+)"
  Pct_table[, chromosome := gsub(locus_regexpr, "\\1", introgression_locus)]
  Pct_table[, start := as.numeric(gsub(locus_regexpr, "\\2", introgression_locus))]
  Pct_table[, end := as.numeric(gsub(locus_regexpr, "\\3", introgression_locus))]
  Pct_table <- merge(Pct_table, chrom_length, by = "chromosome")
  Pct_table[, start_pos := start + cumPos]
  Pct_table[, end_pos := end + cumPos]
  Pct_table[, mid_pos := (start+end)/2 + cumPos]
  Pct_table[1:.N, ]
}

########## list measures to be analyzed ##############
Measure_table <- data.table(
  measure = c("CREoverlap", "emVar", "ActivityIncreasing",  "ActivityDecreasing", "ExpressionIncreasing",  "ExpressionDecreasing", "sign", "sign_expression"),
  measure_label = c("Pct of CRE overlapping SNPs", "Pct of emVars", "Pct of activity increasing alleles", "Pct of activity decreasing alleles", "Pct of expression increasing alleles","Pct of expression decreasing alleles", "Pct of activity increasing alleles\n(among Emvars)", "Pct of expression increasing alleles\n(among Emvars)"),
  null = c(NA, NA, NA, NA, NA, NA, 0.5, 0.5)
)

model_Pct_emVars <- function(emVar_table, testedVar = NULL, adjustVar = NULL) {
  
emVar_table <- emVar_table[seq_len(.N), ]
  emVar_table[, is_enhancer_CRITERIA := grepl("enhancer", cre_class_CRITERIA)]
  emVar_table[, is_silencer_CRITERIA := grepl("silencer", cre_class_CRITERIA)]
  mod <- glm(is_emVar_CRITERIA ~ split_by, family = "binomial", data = emVar_table[is_tested_CRE == TRUE,])
  library(lme4)
  mod <- glmer(is_emVar_CRITERIA ~ I((nBCs_g1a1+nBCs_g1a2)/100) + (1|ID) + COND_ID + POP_adaptive + Introgression_scenario_top + Introgression_source_top, data=emVARs_obs_ctc[is_tested_cre==TRUE & Introgression_scenario_top!='',],family=binomial(link="logit"))
  summary(mod)
  mod <- glmer(is_emVar_CRITERIA ~ I((nBCs_g1a1+nBCs_g1a2)/100) + (1|ID) + COND_ID + POP_adaptive + Introgression_scenario_top + Introgression_source_top, data=emVARs_obs_ctc[is_tested_cre==TRUE & Introgression_scenario_top!='',],family=binomial(link="logit"))
  summary(mod)
}
