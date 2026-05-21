# 03_00_parameter_activity.R

# current default
# CRITERIA_ACTIVE='lfdr20_GCnorm'
# CRITERIA_DIFF='Diff_lfdr20_FC.2'

if (CRITERIA_DIFF == "Diff_lfdr20_FC.2") {
  ### threshold choice
  LFDR_TH_DIFF <- 0.2
  FDR_TH_DIFF <- Inf
  LOG2FC_TH_DIFF <- 0.2
}

if (CRITERIA_DIFF == "Diff_lfdr10_FC.2") {
  ### threshold choice
  LFDR_TH_DIFF <- 0.1
  FDR_TH_DIFF <- Inf
  LOG2FC_TH_DIFF <- 0.2
}

if (CRITERIA_DIFF == "Diff_FDR10_FC.2") {
  ### threshold choice
  LFDR_TH_DIFF <- Inf
  FDR_TH_DIFF <- 0.1
  LOG2FC_TH_DIFF <- 0.2
}


if (CRITERIA_DIFF == "Diff_FDR5_FC.2") {
  ### threshold choice
  LFDR_TH_DIFF <- Inf
  FDR_TH_DIFF <- 0.05
  LOG2FC_TH_DIFF <- -Inf
}


if (CRITERIA_DIFF == "Diff_FDR1_FC.2") {
  ### threshold choice
  LFDR_TH_DIFF <- Inf
  FDR_TH_DIFF <- 0.01
  LOG2FC_TH_DIFF <- 0.2
}


### threshold & alpha definition choice
# CRITERIA='FDR5_FCthresholds_GCnorm'
# ALPHA_CRITERIA='alpha.GC_norm'
# ALPHA_SE_CRITERIA='alpha.se_norm'
# LFDR_TH_DIFF=Inf
# FDR_TH_DIFF=0.05
# LOG2FC_TH_DIFF=.2
