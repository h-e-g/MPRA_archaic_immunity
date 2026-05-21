# 06_00_parameter_emvar_Diff.R

# current default
# CRITERIA_ACTIVE='lfdr20_GCnorm'
# CRITERIA_EMVARS='Diff_lfdr20_FC.2'

if (CRITERIA_EMVARS_DIFF == "EmVarDiff_lfdr20_FC.2") {
  ### threshold choice
  LFDR_TH_EMVARDIFF <- 0.2
  FDR_TH_EMVARDIFF <- Inf
  LOG2FC_TH_EMVARDIFF <- 0.2
}

if (CRITERIA_EMVARS_DIFF == "EmVarDiff_lfdr10_FC.2") {
  ### threshold choice
  LFDR_TH_EMVARDIFF <- 0.1
  FDR_TH_EMVARDIFF <- Inf
  LOG2FC_TH_EMVARDIFF <- 0.2
}

if (CRITERIA_EMVARS_DIFF == "EmVarDiff_FDR10_FC.2") {
  ### threshold choice
  LFDR_TH_EMVARDIFF <- Inf
  FDR_TH_EMVARDIFF <- 0.1
  LOG2FC_TH_EMVARDIFF <- 0.2
}

if (CRITERIA_EMVARS_DIFF == "EmVarDiff_FDR5_FC.2") {
  ### threshold choice
  LFDR_TH_EMVARDIFF <- Inf
  FDR_TH_EMVARDIFF <- 0.05
  LOG2FC_TH_EMVARDIFF <- -Inf
}

if (CRITERIA_EMVARS_DIFF == "EmVarDiff_FDR1_FC.2") {
  ### threshold choice
  LFDR_TH_EMVARDIFF <- Inf
  FDR_TH_EMVARDIFF <- 0.01
  LOG2FC_TH_EMVARDIFF <- 0.2
}
