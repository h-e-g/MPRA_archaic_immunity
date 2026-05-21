# 03_00_parameter_activity.R

# NORMALIZATION=c('none', 'center','scale_perm','scale_all','scale_nBC_perm')
# GC_CORRECT=c('none','observed','scrambled')

########## ########## ########## ########## ########## ########## ########## ########## 
########## LFDR based ########## ########## ########## ########## ########## ########## 
########## ########## ########## ########## ########## ########## ########## ########## 

if (CRITERIA_ACTIVE == "lfdr20_scrambled1pct_GCnorm") {
  #### threshold & alpha definition choice
  ALPHA_CRITERIA <- "alpha.GC_norm"
  ALPHA_SE_CRITERIA <- "alpha.se_norm"
  LFDR_TH <- 0.2
  FDR_TH <- Inf
  LOG2FC_TH <- -Inf
  SCRAMBLED_TH <- 0.01
  LOG2FC_TH_STRONG <- -Inf
  GC_CORRECT <- "scrambled"
  NORMALIZATION <- "scale_all"
}



if (CRITERIA_ACTIVE == "lfdr10_scrambled1pct_GCnorm") {
  #### threshold & alpha definition choice
  ALPHA_CRITERIA <- "alpha.GC_norm"
  ALPHA_SE_CRITERIA <- "alpha.se_norm"
  LFDR_TH <- 0.1
  FDR_TH <- Inf
  LOG2FC_TH <- -Inf
  SCRAMBLED_TH <- 0.01
  LOG2FC_TH_STRONG <- -Inf
  GC_CORRECT <- "scrambled"
  NORMALIZATION <- "scale_all"
}


if (CRITERIA_ACTIVE == "lfdr20_scrambled1pct") {
  #### threshold & alpha definition choice
  ALPHA_CRITERIA <- "alpha"
  ALPHA_SE_CRITERIA <- "alpha.se"
  LFDR_TH <- 0.2
  FDR_TH <- Inf
  LOG2FC_TH <- -Inf
  SCRAMBLED_TH <- 0.01
  LOG2FC_TH_STRONG <- -Inf
  GC_CORRECT <- "none"
  NORMALIZATION <- "center"
}



if (CRITERIA_ACTIVE == "lfdr20_scrambled1pct_GC") {
  #### threshold & alpha definition choice
  ALPHA_CRITERIA <- "alpha.GC"
  ALPHA_SE_CRITERIA <- "alpha.se"
  LFDR_TH <- 0.2
  FDR_TH <- Inf
  LOG2FC_TH <- -Inf
  SCRAMBLED_TH <- 0.01
  LOG2FC_TH_STRONG <- -Inf
  GC_CORRECT <- "scrambled"
  NORMALIZATION <- "center"
}


if (CRITERIA_ACTIVE == "lfdr20_scrambled1pct_GC_stabVar") {
  #### threshold & alpha definition choice
  ALPHA_CRITERIA <- "alpha.GC"
  ALPHA_SE_CRITERIA <- "alpha.se"
  LFDR_TH <- 0.2
  FDR_TH <- Inf
  LOG2FC_TH <- -Inf
  SCRAMBLED_TH <- 0.01
  LOG2FC_TH_STRONG <- -Inf
  GC_CORRECT <- "scrambled"
  NORMALIZATION <- "scale_nBC_perm"
}


########## ########## ########## ########## ########## ########## ########## ########## 
########## FDR based ########## ########## ########## ########## ########## ########## 
########## ########## ########## ########## ########## ########## ########## ########## 

##########-------------------------------------------------------------------########## 
########## strong defined from scrambled                                         ######
##########-------------------------------------------------------------------########## 

# if (CRITERIA_ACTIVE == "FDR5_scrambled5PctFC1_FC0.3_GCnorm") {
#   #### threshold & alpha definition choice
#   ALPHA_CRITERIA <- "alpha.GC_norm"
#   ALPHA_SE_CRITERIA <- "alpha.se_norm"
#   LFDR_TH <- Inf
#   FDR_TH <- 0.05
#   LOG2FC_TH <- 0.3
#   SCRAMBLED_TH <- 0.05
#   LOG2FC_TH_STRONG <- 1
#   GC_CORRECT <- "scrambled"
#   NORMALIZATION <- "scale_all"
# }

# # if (CRITERIA_ACTIVE == "FDR5_FC1_FC0.2_GCnorm") {
# #   ### threshold & alpha definition choice
# #   ALPHA_CRITERIA <- "alpha.GC_norm"
# #   ALPHA_SE_CRITERIA <- "alpha.se_norm"
# #   LFDR_TH <- Inf
# #   FDR_TH <- 0.05
# #   LOG2FC_TH <- 0.3
# #   SCRAMBLED_TH <- Inf
# #   LOG2FC_TH_STRONG <- 1
# #   GC_CORRECT <- "scrambled"
# #   NORMALIZATION <- "scale_all"
# # }


# if (CRITERIA_ACTIVE == "FDR1_scrambled5pct_GCnorm") {
#   ### threshold & alpha definition choice
#   ALPHA_CRITERIA <- "alpha.GC_norm"
#   ALPHA_SE_CRITERIA <- "alpha.se_norm"
#   LFDR_TH <- Inf
#   FDR_TH <- 0.01
#   LOG2FC_TH <- -Inf
#   SCRAMBLED_TH <- 0.05
#   LOG2FC_TH_STRONG <- -Inf
#   GC_CORRECT <- "scrambled"
#   NORMALIZATION <- "scale_all"
# }


# if (CRITERIA_ACTIVE == "FDR1_scrambled5pct_FC0.2_GCnorm") {
#   ### threshold & alpha definition choice
#   ALPHA_CRITERIA <- "alpha.GC_norm"
#   ALPHA_SE_CRITERIA <- "alpha.se_norm"
#   LFDR_TH <- Inf
#   FDR_TH <- 0.01
#   LOG2FC_TH <- 0.2
#   SCRAMBLED_TH <- 0.05
#   LOG2FC_TH_STRONG <- 0.3
#   GC_CORRECT <- "scrambled"
#   NORMALIZATION <- "scale_all"
# }



# if (CRITERIA_ACTIVE == "FDR1_scrambled5pct_FC0.3_GCnorm") {
#   ### threshold & alpha definition choice
#   ALPHA_CRITERIA <- "alpha.GC_norm"
#   ALPHA_SE_CRITERIA <- "alpha.se_norm"
#   LFDR_TH <- Inf
#   FDR_TH <- 0.01
#   LOG2FC_TH <- 0.3
#   SCRAMBLED_TH <- 0.05
#   LOG2FC_TH_STRONG <- 0.3
#   GC_CORRECT <- "scrambled"
#   NORMALIZATION <- "scale_all"
# }

# if (CRITERIA_ACTIVE == "FDR1_scrambled5pct_FC0.2_GC_stabVar") {
#   ### threshold & alpha definition choice
#   ALPHA_CRITERIA <- "alpha.GC_norm"
#   ALPHA_SE_CRITERIA <- "alpha.se_norm"
#   LFDR_TH <- Inf
#   FDR_TH <- 0.01
#   LOG2FC_TH <- 0.2
#   SCRAMBLED_TH <- 0.05
#   LOG2FC_TH_STRONG <- 0.1
#   GC_CORRECT <- "scrambled"
#   NORMALIZATION <- "scale_all"
# }


# if (CRITERIA_ACTIVE == "FDR1_scrambled1pct_FC0.3_GCnorm") {
#   ### threshold & alpha definition choice
#   ALPHA_CRITERIA <- "alpha.GC_norm"
#   ALPHA_SE_CRITERIA <- "alpha.se_norm"
#   LFDR_TH <- Inf
#   FDR_TH <- 0.01
#   LOG2FC_TH <- 0.3
#   SCRAMBLED_TH <- 0.01
#   LO2FC_TH_STRONG <- 1
#   GC_CORRECT <- "scrambled"
#   NORMALIZATION <- "scale_all"
# }


# if (CRITERIA_ACTIVE == "FDR1_scrambled1pct_FC0.2_GCnorm") {
#   ### threshold & alpha definition choice
#   ALPHA_CRITERIA <- "alpha.GC_norm"
#   ALPHA_SE_CRITERIA <- "alpha.se_norm"
#   LFDR_TH <- Inf
#   FDR_TH <- 0.01
#   LOG2FC_TH <- 0.2
#   SCRAMBLED_TH <- 0.01
#   LO2FC_TH_STRONG <- 1
#   GC_CORRECT <- "scrambled"
#   NORMALIZATION <- "scale_all"
# }


################################
############# tests ############
################################

if (CRITERIA_ACTIVE == "FDR5_scrambled5pct_FC0.2") {
  ### threshold & alpha definition choice
  ALPHA_CRITERIA <- "alpha.GC_norm"
  ALPHA_SE_CRITERIA <- "alpha.se_norm"
  LFDR_TH <- Inf
  FDR_TH <- 0.05
  LOG2FC_TH <- 0.2
  SCRAMBLED_TH <- 0.05
  LOG2FC_TH_STRONG <- -Inf
  GC_CORRECT <- "none"
  NORMALIZATION <- "center"
}

if (CRITERIA_ACTIVE == "FDR5_scrambled5pct_FC0.2_norm") {
  ### threshold & alpha definition choice
  ALPHA_CRITERIA <- "alpha.GC_norm"
  ALPHA_SE_CRITERIA <- "alpha.se_norm"
  LFDR_TH <- Inf
  FDR_TH <- 0.05
  LOG2FC_TH <- 0.2
  SCRAMBLED_TH <- 0.05
  LOG2FC_TH_STRONG <- -Inf
  GC_CORRECT <- "none"
  NORMALIZATION <- "scale_all"
}


if (CRITERIA_ACTIVE == "FDR5_scrambled5pct_FC0.2_GC") {
  ### threshold & alpha definition choice
  ALPHA_CRITERIA <- "alpha.GC_norm"
  ALPHA_SE_CRITERIA <- "alpha.se_norm"
  LFDR_TH <- Inf
  FDR_TH <- 0.05
  LOG2FC_TH <- 0.2
  SCRAMBLED_TH <- 0.05
  LOG2FC_TH_STRONG <- -Inf
  GC_CORRECT <- "scrambled"
  NORMALIZATION <- "center"
}

################################
############# MAIN #############
################################

if (CRITERIA_ACTIVE == "FDR5_scrambled5pct_FC0.2_GCnorm") {
  ### threshold & alpha definition choice
  ALPHA_CRITERIA <- "alpha.GC_norm"
  ALPHA_SE_CRITERIA <- "alpha.se_norm"
  LFDR_TH <- Inf
  FDR_TH <- 0.05
  LOG2FC_TH <- 0.2
  SCRAMBLED_TH <- 0.05
  LOG2FC_TH_STRONG <- -Inf
  GC_CORRECT <- "scrambled"
  NORMALIZATION <- "scale_all"
}

##########-------------------------------------------------------------------########## 
########## strong defined from Fold change                                       ######
##########-------------------------------------------------------------------########## 


################################
############# tests ############
################################

if (CRITERIA_ACTIVE == "FDR5_FC1_FC0.2") {
  ### threshold & alpha definition choice
  ALPHA_CRITERIA <- "alpha.GC_norm"
  ALPHA_SE_CRITERIA <- "alpha.se_norm"
  LFDR_TH <- Inf
  FDR_TH <- 0.05
  LOG2FC_TH <- 0.2
  SCRAMBLED_TH <- Inf
  LOG2FC_TH_STRONG <- 1
  GC_CORRECT <- "none"
  NORMALIZATION <- "center"
}

if (CRITERIA_ACTIVE == "FDR5_FC1_FC0.2_norm") {
  ### threshold & alpha definition choice
  ALPHA_CRITERIA <- "alpha.GC_norm"
  ALPHA_SE_CRITERIA <- "alpha.se_norm"
  LFDR_TH <- Inf
  FDR_TH <- 0.05
  LOG2FC_TH <- 0.2
  SCRAMBLED_TH <- Inf
  LOG2FC_TH_STRONG <- 1
  GC_CORRECT <- "none"
  NORMALIZATION <- "scale_all"
}


if (CRITERIA_ACTIVE == "FDR5_FC1_FC0.2_GC") {
  ### threshold & alpha definition choice
  ALPHA_CRITERIA <- "alpha.GC_norm"
  ALPHA_SE_CRITERIA <- "alpha.se_norm"
  LFDR_TH <- Inf
  FDR_TH <- 0.05
  LOG2FC_TH <- 0.2
  SCRAMBLED_TH <- Inf
  LOG2FC_TH_STRONG <- 1
  GC_CORRECT <- "scrambled"
  NORMALIZATION <- "center"
}

################################
############# MAIN #############
################################

if (CRITERIA_ACTIVE == "FDR5_FC1_FC0.2_GCnorm") {
  ### threshold & alpha definition choice
  ALPHA_CRITERIA <- "alpha.GC_norm"
  ALPHA_SE_CRITERIA <- "alpha.se_norm"
  LFDR_TH <- Inf
  FDR_TH <- 0.05
  LOG2FC_TH <- 0.2
  SCRAMBLED_TH <- Inf
  LOG2FC_TH_STRONG <- 1
  GC_CORRECT <- "scrambled"
  NORMALIZATION <- "scale_all"
}
