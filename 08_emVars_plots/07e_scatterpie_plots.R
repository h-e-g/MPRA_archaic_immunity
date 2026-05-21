# running: sbatch --array=1-96 -p geh,common --mem=40G  $SCRIPT_DIR/00_Rscript.sh MPRA_count_exp6_analysis3/07b_scatterpie_plots.R --batch

# DONE: make sure that posID is include in file name
# DONE: make sure that rsID is included in file name
# DONE: add info of Introgressed_from & SNP_type
# DONE: differentiate between frequency of introgressed allele and intogressed haplotype (to be done in SNP annotation).
# DONE: add annotation of putative targets



MPRA_DIR <- "/pasteur/helix/projects/evo_immuno_pop/MPRA"
ANALYSIS_DIR <- "MPRA_count_exp6_analysis3"
source(sprintf("%s/scripts/%s/00_starter.R", MPRA_DIR, ANALYSIS_DIR))


# Charger les packages nécessaires
library(ggplot2)
library(ggmap)
library(maps)
library(scatterpie)
library(janitor)
library(ggrepel)
library(cowplot)
library(sf)
library(rnaturalearth)
library(ggrepel)


RUN_ID <- "RUN3_Z2_nBC10"
CRITERIA_ACTIVE_OUT <- "FDR5_scrambled5pct_FC0.2_GCnorm"
CRITERIA_DIFF_OUT <- "noDiffFilter_Diff_FDR5_FC.2"
CRITERIA_EMVARS_OUT <- "EmVar_FDR5_FC.2"
CRITERIA_EMVARS_DIFF <- "EmVarDiff_FDR5_FC.2"
STIM_ONLY <- FALSE

cmd <- commandArgs(trailingOnly = TRUE)
print(cmd)

for (i in seq_along(cmd)) {
  if (cmd[i] == "--run_dir" || cmd[i] == "-r") {
    RUN_ID <- cmd[i + 1]
  }
  if (cmd[i] == "--criteria_active" || cmd[i] == "-a") {
    CRITERIA_ACTIVE_OUT <- cmd[i + 1]
  }
  if (cmd[i] == "--criteria_diff" || cmd[i] == "-d") {
    CRITERIA_DIFF_OUT <- cmd[i + 1]
  }
  if (cmd[i] == "--criteria_emvar" || cmd[i] == "-e") {
    CRITERIA_EMVARS_OUT <- cmd[i + 1]
  }
  if (cmd[i] == "--criteria_emvar_diff" || cmd[i] == "-ed") {
    CRITERIA_EMVARS_DIFF <- cmd[i + 1]
  }
  if (cmd[i] == "--stim_only" || cmd[i] == "-s") {
    STIM_ONLY <- as.logical(cmd[i + 1])
  }
  if (cmd[i] == "--batch" || cmd[i] == "-b") {
    BATCH <- as.numeric(cmd[i + 1])
  } # ID of the set of SNPs to plot
}

n_per_batch <- 50

CRITERIA_ACTIVE <- gsub("noActivityFilter_", "", CRITERIA_ACTIVE_OUT)
CRITERIA_DIFF <- gsub("noDiffFilter_", "", CRITERIA_DIFF_OUT)
CRITERIA_EMVARS <- gsub("noEmVarFilter_", "", CRITERIA_EMVARS_OUT)
STIM_ONLY_OUT <- ifelse(STIM_ONLY, "stimOnly", "all")

tic("loading oligo & SNP annotations")
# load annotation of oligos (beforz: CRS)
oligo_source <- fread(sprintf("%s/data/%s/00_oligo_annot_v2.txt", MPRA_DIR, ANALYSIS_DIR))
# load annotations of SNPs
SNP_annot_v4 <- fread(sprintf("%s/data/%s/00_SNP_annot_v4.txt", MPRA_DIR, ANALYSIS_DIR))
# load annotation per POP
selected_annot <- fread(sprintf("%s/data/%s/selected_snp_annotation.tsv.gz", MPRA_DIR, ANALYSIS_DIR))
# load clean annotation of introgression
selected_annot_wide <- fread(sprintf("%s/data/%s/selected_snp_annotation_wide.tsv.gz", MPRA_DIR, ANALYSIS_DIR))
SNP_annot_v5 <- merge(SNP_annot_v4[, -"Introgression_scenario"], selected_annot_wide, by = c("ID", "posID"))
toc()
# load annotations of SNPs
SNP_annot <- SNP_annot_v5

##### active parameters
source(sprintf("%s/scripts/%s/03_00_parameter_activity.R", MPRA_DIR, ANALYSIS_DIR))
##### diff parameters
source(sprintf("%s/scripts/%s/04_00_parameter_diff_activity.R", MPRA_DIR, ANALYSIS_DIR))
##### emVars parameters
source(sprintf("%s/scripts/%s/05_00_parameter_emVars.R", MPRA_DIR, ANALYSIS_DIR))
##### diff parameters
source(sprintf("%s/scripts/%s/06_00_parameter_emVars_Diff.R", MPRA_DIR, ANALYSIS_DIR))

IN_DIR <- sprintf("%s/data/%s/03a_aggMPRA_analyse/%s/", MPRA_DIR, ANALYSIS_DIR, RUN_ID)
ACTIVE_DIR <- sprintf("%s/data/%s/03b_computeFDR/%s/Activity/%s/", MPRA_DIR, ANALYSIS_DIR, RUN_ID, CRITERIA_ACTIVE)
DIFF_DIR <- sprintf("%s/data/%s/03b_computeFDR/%s/Diff/%s/%s/", MPRA_DIR, ANALYSIS_DIR, RUN_ID, CRITERIA_DIFF, CRITERIA_ACTIVE_OUT)
EMVAR_DIR <- sprintf("%s/data/%s/03b_computeFDR/%s/emVars/%s/%s", MPRA_DIR, ANALYSIS_DIR, RUN_ID, CRITERIA_EMVARS, CRITERIA_ACTIVE_OUT)
EMVARDIFF_DIR <- sprintf("%s/data/%s/03b_computeFDR/%s/EmVar_Diff/%s/%s/%s/%s/%s", MPRA_DIR, ANALYSIS_DIR, RUN_ID, CRITERIA_EMVARS_DIFF, CRITERIA_EMVARS_OUT, CRITERIA_DIFF_OUT, CRITERIA_ACTIVE_OUT, STIM_ONLY_OUT)

################### load emVar lists ###################
emVARs_obs <- fread(file = sprintf("%s/all_emVars_annotated_celltype__%s.tsv", EMVAR_DIR, CRITERIA_ACTIVE_OUT))
emVARs_obs[, posID := crsID]
emVARs_obs <- emVARs_obs[is_emVar_CRITERIA == TRUE, ]

response_emVars <- fread(sprintf("%s/all_emVars_diff_obs_response.tsv", EMVARDIFF_DIR))
response_emVars[, posID := crsID]
response_emVars <- response_emVars[is_emVar_Diff_CRITERIA == TRUE]
# myPosID <- response_emVars[1, posID]

celltype_emVars <- fread(sprintf("%s/all_emVars_diff_obs_celltype.tsv", EMVARDIFF_DIR))
celltype_emVars[, posID := crsID]
celltype_emVars <- celltype_emVars[is_emVar_Diff_CRITERIA == TRUE]

###############  define output directories ##############

FIGURE_DIR <- sprintf("%s/figures/%s/%s", MPRA_DIR, ANALYSIS_DIR, RUN_ID)
#FIGURE_DIR <- sprintf("%s/figures/%s/", MPRA_DIR, ANALYSIS_DIR)
FIGURE_EMVAR_DIR <- sprintf("%s/04_emVars/%s/%s", FIGURE_DIR, CRITERIA_EMVARS, CRITERIA_ACTIVE_OUT)
FIGURE_EMVARDIFF_DIR <- sprintf("%s/05_emVarDiff/%s/%s/%s/%s/%s", FIGURE_DIR, CRITERIA_EMVARS_DIFF, CRITERIA_EMVARS_OUT, CRITERIA_DIFF_OUT, CRITERIA_ACTIVE_OUT, STIM_ONLY_OUT)

# emVars
dir.create(sprintf("%s/emVar_plots/all_conditions", FIGURE_EMVAR_DIR), showWarnings = FALSE, recursive = TRUE)
dir.create(sprintf("%s/emVar_plots/NS_only", FIGURE_EMVAR_DIR), showWarnings = FALSE, recursive = TRUE)
dir.create(sprintf("%s/emVar_plots/HepG2_only", FIGURE_EMVAR_DIR), showWarnings = FALSE, recursive = TRUE)
dir.create(sprintf("%s/emVar_plots/A549_only", FIGURE_EMVAR_DIR), showWarnings = FALSE, recursive = TRUE)
dir.create(sprintf("%s/emVar_plots/K562_only", FIGURE_EMVAR_DIR), showWarnings = FALSE, recursive = TRUE)

# response
dir.create(sprintf("%s/response_emVar_plots/all_conditions", FIGURE_EMVARDIFF_DIR), showWarnings = FALSE, recursive = TRUE)
dir.create(sprintf("%s/response_emVar_plots/HepG2_only", FIGURE_EMVARDIFF_DIR), showWarnings = FALSE, recursive = TRUE)
dir.create(sprintf("%s/response_emVar_plots/A549_only", FIGURE_EMVARDIFF_DIR), showWarnings = FALSE, recursive = TRUE)
dir.create(sprintf("%s/response_emVar_plots/K562_only", FIGURE_EMVARDIFF_DIR), showWarnings = FALSE, recursive = TRUE)
# celltypes
dir.create(sprintf("%s/celltype_emVar_plots/all_conditions", FIGURE_EMVARDIFF_DIR), showWarnings = FALSE, recursive = TRUE)
dir.create(sprintf("%s/celltype_emVar_plots/NS_only", FIGURE_EMVARDIFF_DIR), showWarnings = FALSE, recursive = TRUE)

############### define target gene ################

oligo_target <- fread(sprintf("%s/data/%s/00_oligo_targets_v4.txt", MPRA_DIR, ANALYSIS_DIR), sep = "\t")
oligo_target_collapsed <- fread(sprintf("%s/data/%s/00_oligo_targets_collapsed_v4.txt", MPRA_DIR, ANALYSIS_DIR), sep = "\t")
crs_targets <- unique(oligo_target_collapsed[, -"oligo"])
NearestGene <- unique(oligo_target[method == "NearestGene", .(posID = crsID, TargetGene = ifelse(gene_type != "", TargetGene, ""))])
Target_Gene <- oligo_target[gene_type != "", ][order(-grepl("i|a|p", gene_type), -grepl("Nearest", method), -Score1, -Score2), ][!duplicated(crsID)]
Target_Gene[, targetName := paste0(TargetGene, " (", gene_type, "|", method, "|", rank, ")")]


############### load population data ##############

pop_annotation <- clean_names(fread(sprintf("%s/data/igsr_populations.tsv", MPRA_DIR)))
# pop_annotation[,.N,by=population_code]
# SGDP_pos=pop_annotation[grep('SGDP',superpopulation_name),.(latitude=mean(population_latitude),longitude=mean(population_longitude)),by=.(population=superpopulation_name)]
SGDP_pos <- data.table(
  population = c("SGDP_African", "SGDP_WestEurasian", "SGDP_Eastasian", "SGDP_SouthAsian", "SGDP_Agta", "SGDP_Papuan"),
  latitude = c(1.436059, 43.864109, 30.456519, 25.309581, 17.95, -4),
  longitude = c(16.44148, 26.37037, 109.86735, 76.63742, 121.26, 143)
)
OKG_pos <- pop_annotation[population_code != "", .(
  latitude = mean(population_latitude),
  longitude = mean(population_longitude)
),
by = .(population = population_code)
]
OKG_pos[population == "STU", latitude := 7.87] # Wikipedia coordinates for Sri Lanka
OKG_pos[population == "STU", longitude := 80.76] # Wikipedia coordinates for Sri Lanka

population_pos <- rbind(SGDP_pos, OKG_pos)

archaic_pos <- data.table(
  population = c("Vindija", "Altai", "Chagyrskaya", "Denisova"),
  latitude = c(-85, -85, -85, -85), # Position arbitraire en bas de carte
  longitude = c(30, 100, 160, 270)
) # Position arbitraire en bas de carte)

# SNP_freq_introgressed <- fread(file = sprintf("%s/full_annotation/haplotypes_aSNP_with_introgressed_and_frequency_%s.tsv.gz", OUT_DIR, CHR), sep = "\t")
SNP_freq_introgressed <- fread(sprintf("%s/data/%s/00_oligo_annot_v2_freqs.txt", MPRA_DIR, ANALYSIS_DIR))

# SNP_annot[,.(posID,ID,CHROM,POS,REF,ALT,rsID,selection_criteria,Introgressed_in,Introgressed_from.perPop,INTROGRESSED_FINAL=)]
SNP_annot_freqs <- melt(SNP_annot,
  id.vars = c("crsID", "ID", "CHROM", "POS_b37", "REF", "ALT", "ANCESTRAL", "DERIVED", "INTROGRESSED.allele"),
  measure.vars = c("YRI.der", "SGDP_African.der", "GBR.der", "IBS.der", "SGDP_WestEurasian.der", "CHB.der", "JPT.der", "SGDP_Eastasian.der", "PJL.der", "STU.der", "SGDP_Agta.der", "SGDP_Papuan.der"),
  value.name = "derived.freq", variable.name = "population"
)

SNP_annot_freqs[, population := gsub(".der", "", population)]
SNP_annot_freqs <- merge(SNP_annot_freqs, population_pos, by = "population")
SNP_annot_freqs <- merge(SNP_annot_freqs, SNP_freq_introgressed[, .(ID, population = POP, freq_introgressed_allele, freq_introgressed_haplotype)], by = c("population", "ID"), all.x = TRUE)

SNP_annot_archaics <- melt(SNP_annot,
  id.vars = c("crsID", "ID", "CHROM", "POS_b37", "REF", "ALT", "ANCESTRAL", "DERIVED", "INTROGRESSED.allele"),
  measure.vars = c("Vindija.der", "Chagyrskaya.der", "Altai.der", "Denisova.der"), variable.name = "population"
)

SNP_annot_archaics[, population := gsub(".der", "", population)]
SNP_annot_archaics[, derived.count := case_when(value == "0/0" ~ 0, value == "0/1" ~ 1, value == "1/1" ~ 2, TRUE ~ NA)]
SNP_annot_archaics[, derived.freq := derived.count / 2]
SNP_annot_archaics <- merge(SNP_annot_archaics, archaic_pos, by = "population")

# SNP_annot_freqs[, freq_introgressed_allele := ifelse(INTROGRESSED_FINAL == DERIVED, derived.freq, 1 - derived.freq)]
# SNP_annot_freqs[, freq_introgressed_haplotype := ifelse(INTROGRESSED_FINAL == DERIVED, derived.freq, 1 - derived.freq)]
SNP_annot_archaics[, freq_introgressed_allele := ifelse(INTROGRESSED.allele == DERIVED, derived.freq, 1 - derived.freq)]
SNP_annot_archaics[, freq_introgressed_haplotype := freq_introgressed_allele]

# # Exemple de données
# # Les fréquences alléliques par population (remplacer par vos données)
# pop_data <- data.table(
#   population = c("Pop1", "Pop2", "Pop3", "Pop4"),
#   latitude = c(45.0, 30.0, 60.0, 40.0),
#   longitude = c(10.0, -20.0, 100.0, 0),
#   derived.freq = c(0.3, 0.6, 0.2, NA)
# )

# # Les fréquences alléliques des génomes archaïques (remplacer par vos données)
# archaic_data <- data.table(
#   population = c("Altai", "Vindija", "Chagyrskaya", "Denisova"),
#   latitude = c(-85, -85, -85, -85), # Position arbitraire dans le coin supérieur droit
#   longitude = seq(30,270,l=4), # Position arbitraire dans le coin supérieur droit
#   derived.freq = c(0, 0, 0, 0.5)
# )


# Télécharger une carte du monde
# API key is needed for using stadia maps
# world_map <- get_stadiamap(bbox = c(left = -180, bottom = -60, right = 180, top = 90), zoom = 2, maptype = "stamen_watercolor")
#  world_map <- map("world", col = "grey", fill = TRUE, bg = "white", lwd = 0.05, mar = rep(0, 4), border = 0, ylim = c(-80, 80))
# world_map <- spTransform(world_map, crs)

# dev.off()


# Dessiner la carte du monde
world <- ne_countries(scale = "medium", returnclass = "sf")

# Robinson projection creates artifacts on the world map due to polygons that cross the meridian along which the world is split
# I found a trick here on how to solve this:
# (https://stackoverflow.com/questions/56146735/visual-bug-when-changing-robinson-projections-central-meridian-with-ggplot2)
# polygon <- st_polygon(x = list(rbind(c(-0.0001, 90),
#                                      c(0, 90),
#                                      c(0, -90),
#                                      c(-0.0001, -90),
#                                      c(-0.0001, 90)))) %>%
#   st_sfc() %>%
#   st_set_crs(4326)
# world2 <- world %>% st_difference(polygon)
# crs <- "+proj=robin +lon_0=180"
# world <- st_transform(world2, crs)
# p <- ggplot() + theme_void() +
#   geom_sf(data = world, fill = "lightgray", color = "white")


# since we -use long0 150 - we add to adapt the polygon
polygon30 <- st_polygon(x = list(rbind(
  c(-30 - 0.0001, 90),
  c(-30, 90),
  c(-30, -90),
  c(-30 - 0.0001, -90),
  c(-30 - 0.0001, 90)
))) %>%
  st_sfc() %>%
  st_set_crs(4326)
# modify world dataset to remove overlapping portions with world's polygons
world2 <- world %>% st_difference(polygon30)
# perform transformation on modified version of world dataset
CRS_robinson_150 <- "+proj=robin +lon_0=150 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
world_robinson <- st_transform(world2, crs = CRS_robinson_150)

# testing:

RUN_ID <- "RUN2_Z2_nBC10"
# IN_DIR <- sprintf("%s/data/%s/03a_aggMPRA_analyse/%s/", MPRA_DIR, ANALYSIS_DIR, RUN_ID)
IN_DIR <- sprintf("%s/data/%s/03b_computeFDR/%s/", MPRA_DIR, ANALYSIS_DIR, RUN_ID)

# myPosID <- "10:122299806"

# emVars_signif_active <- fread(sprintf("%s/all_emVars_activeCRS_ok_annotated_celltype.tsv", IN_DIR))
# emVars_signif_active <- emVars_signif_active[FDR < 0.05, ]
# emVars_signif_active <- merge(emVars_signif_active, condition_summary, by.x = "ANALYSIS_NAME", by.y = "analysis_name")
# emVars_signif_active_ambiguous <- fread(sprintf("%s/all_emVars_activeCRS_annotated_celltype.tsv", IN_DIR))
# emVars_signif_active_ambiguous <- emVars_signif_active_ambiguous[FDR < 0.05 & !(posID %in% emVars_signif_active[, posID])]
# emVars_signif_active_ambiguous <- merge(emVars_signif_active_ambiguous, condition_summary, by.x = "ANALYSIS_NAME", by.y = "analysis_name")

# response_emVars_active_signif <- fread(sprintf("%s/all_response_emVars_active_signif.tsv", IN_DIR))
# response_emVars_active_signif[, posID := crsID]

posID_list <- intersect(c(emVARs_obs$posID, response_emVars$posID, celltype_emVars$posID), SNP_annot[!is.na(YRI.der) | !is.na(SGDP_African.der), unique(posID)])


min_i <- (BATCH - 1) * n_per_batch + 1
max_i <- min((BATCH - 1) * n_per_batch + n_per_batch, length(posID_list))

if (min_i > length(posID_list)) {
  cat("no more variants to plot\n")
  q("no")
}

# posID_list <- SNP_annot_freqs[, unique(posID)]
for (myPosID in posID_list[min_i:max_i]) {
  cat("\n", myPosID, "\n")
  myRSID <- try(SNP_annot[posID == myPosID, get("rsID")])
  # Introgressed_from <- try(SNP_annot[posID == myPosID, get("Adaptive_from")])
  Introgressed_from <- try(SNP_annot[posID == myPosID, get("Introgression_source_top")])
  allele_match <- try(SNP_annot[posID == myPosID, get("allele_match")])
  Introgression_scenario <- try(SNP_annot[posID == myPosID, get("Introgression_scenario_top")])
  predicted_Target <- try(Target_Gene[crsID == myPosID, targetName])
  if (class(myRSID) == "try-error" || is.na(myRSID) || myRSID == "") {
    myVariantID <- myPosID
    rsID_print <- gsub(":", "-", myPosID)
  } else {
    myVariantID <- paste0(unique(myRSID), " (", myPosID, ")")
    rsID_print <- paste(gsub(":", "-", myPosID), unique(myRSID), sep = "_")
  }

  # criteria_list <- unlist(str_split(SNP_annot[posID == myPosID, selection_criteria], ","))
  # for (type.CRS in criteria_list) {
  # Fusionner les deux jeux de données
  pop_data <- SNP_annot_freqs[crsID == myPosID, .(crsID, ANCESTRAL, DERIVED, INTROGRESSED.allele, population, latitude, longitude, derived.freq, freq_introgressed_allele, freq_introgressed_haplotype, type = "modern")]
  archaic_data <- SNP_annot_archaics[crsID == myPosID, .(crsID, ANCESTRAL, DERIVED, INTROGRESSED.allele, population, latitude, longitude, derived.freq, freq_introgressed_allele, freq_introgressed_haplotype, type = "archaic")]
  combined_data <- rbind(pop_data, archaic_data)
  combined_data[, ancestral.freq := 1 - derived.freq]
  combined_data[, non_introgressed.freq := 1 - freq_introgressed_haplotype]
  combined_data[, non_introgressed.allele := ifelse(INTROGRESSED.allele == DERIVED, ANCESTRAL, DERIVED)]
  combined_data[, freq_introgressed_allele_different_haplotype := freq_introgressed_allele - freq_introgressed_haplotype]

  combined_data[, introgressed_haplotype.freq := freq_introgressed_haplotype]
  combined_data[, introgressed_allele.freq := freq_introgressed_allele_different_haplotype]
  combined_data[, non_introgressed_allele.freq := 1 - freq_introgressed_allele]
  # Transformer les coordonnées des pie charts
  coordinates <- st_as_sf(combined_data, coords = c("longitude", "latitude"), crs = 4326)
  coordinates <- st_transform(coordinates, CRS_robinson_150)
  combined_data$longitude <- st_coordinates(coordinates)[, 1]
  combined_data$latitude <- st_coordinates(coordinates)[, 2]
  combined_data[type == "archaic", latitude := latitude * 1.6]
  combined_data[type == "archaic", longitude := longitude * 1.5]
  combined_data[population == "SGDP_Agta", population := "Agta"]

  color_alleles <- c("derived.freq" = "#1a70cc", "ancestral.freq" = "#f9c14a")
  color_labels <- paste(c(combined_data[, unique(DERIVED)], combined_data[, unique(ANCESTRAL)]), c("(derived)", "(ancestral)"))

  # plot frequency of derived and ancestral alleles
  p <- ggplot() +
    geom_sf(data = world_robinson, fill = "lightgray", color = "white")
  p <- p + geom_scatterpie(
    aes(x = longitude, y = latitude, group = population),
    data = combined_data,
    cols = c("derived.freq", "ancestral.freq"),
    color = NA
  )
  p <- p + theme_light() + scale_fill_manual(values = color_alleles, labels = color_labels, name = "")
  p <- p + scale_x_continuous(breaks = seq(-180, 180, by = 30))
  p <- p + scale_y_continuous(breaks = seq(-180, 180, by = 30))
  p <- p + theme_minimal_grid() + coord_sf(crs = CRS_robinson_150)
  p <- p + geom_text_repel(data = combined_data[type == "modern", ], aes(x = longitude, y = latitude, label = gsub("SGDP_(.*)", "\\1\n(SGDP)", population)))
  p <- p + geom_text(data = combined_data[type == "archaic", ], aes(x = longitude + 1e6, y = latitude, label = population), hjust = 0)
  p <- p + theme(panel.grid.major = element_line(colour = "lightgrey", size = 0.25, linetype = 3))
  p <- p + theme(legend.position = "bottom") + labs(title = paste0("Allele frequency: ", myVariantID), subtitle = paste0("source: ", Introgressed_from, " - ", Introgression_scenario, "\n", "allele match: ", allele_match, "\n", "target: ", predicted_Target)) + xlab("") + ylab("")

  # emVars
  file_name <- sprintf("%s/emVar_plots/all_conditions/variant_%s_01c_scatterpie_plots_B37.pdf", FIGURE_EMVAR_DIR, rsID_print)
  pdf(file_name, height = 6, width = 7)
  print(p)
  dev.off()
  for (SUBDIR in c("NS_only", "HepG2_only", "A549_only", "K562_only")) {
    target_file <- sprintf("%s/emVar_plots/%s/variant_%s_02c_scatterpie_plots_B37.pdf", FIGURE_EMVAR_DIR, SUBDIR, rsID_print)
    system(sprintf("cp %s %s ", file_name, target_file))
  }

  if (myPosID %chin% response_emVars[, unique(posID)]) {
    # response emVars
    file_name <- sprintf("%s/response_emVar_plots/all_conditions/variant_%s_01c_scatterpie_plots_B37.pdf", FIGURE_EMVARDIFF_DIR, rsID_print)
    pdf(file_name, height = 6, width = 7)
    print(p)
    dev.off()

    for (CELLINE in c("HepG2", "A549", "K562")) {
      SUBDIR <- sprintf("%s_only", CELLINE)
      if (CELLINE %in% response_emVars[posID == myPosID, celline]) {
        target_file <- sprintf("%s/response_emVar_plots/%s/variant_%s_02c_scatterpie_plots_B37.pdf", FIGURE_EMVARDIFF_DIR, SUBDIR, rsID_print)
        system(sprintf("cp %s %s ", file_name, target_file))
      }
    }
  }

  if (myPosID %chin% celltype_emVars[, unique(posID)]) {
    # celltype emVars
    file_name <- sprintf("%s/celltype_emVar_plots/all_conditions/variant_%s_01c_scatterpie_plots_B37.pdf", FIGURE_EMVARDIFF_DIR, rsID_print)
    pdf(file_name, height = 6, width = 7)
    print(p)
    dev.off()

    if (any(celltype_emVars[posID == myPosID, grepl("NS", group1_labels) & grepl("NS", group2_labels)])) {
      target_file <- sprintf("%s/celltype_emVar_plots/NS_only/variant_%s_02c_scatterpie_plots_B37.pdf", FIGURE_EMVARDIFF_DIR, rsID_print)
      system(sprintf("cp %s %s ", file_name, target_file))
    }
  }


  # plot frequency of derived and ancestral alleles (separating introgressed and non introgressed haplotypes)
  color_alleles <- c("introgressed_haplotype.freq" = "#1a70cc", "introgressed_allele.freq" = "#cc5b1a", "non_introgressed_allele.freq" = "#f9c14a")
  allele_labels <- c(combined_data[, rep(unique(INTROGRESSED.allele), 2)], combined_data[, unique(non_introgressed.allele)])
  class_labels <- combined_data[, ifelse(INTROGRESSED.allele == DERIVED, c("(derived, introgressed)", "(derived, nonintrogressed)", "(ancestral)"),
    c("(ancestral, introgressed)", "(ancestral, nonintrogressed)", "(derived)")
  )]
  color_labels <- paste(allele_labels, class_labels)

  p <- ggplot() +
    geom_sf(data = world_robinson, fill = "lightgray", color = "white")
  p <- p + geom_scatterpie(
    aes(x = longitude, y = latitude, group = population),
    data = combined_data,
    cols = c("introgressed_haplotype.freq", "introgressed_allele.freq", "non_introgressed_allele.freq"),
    color = NA
  )
  p <- p + theme_light() + scale_fill_manual(values = color_alleles, labels = color_labels, name = "")
  p <- p + scale_x_continuous(breaks = seq(-180, 180, by = 30))
  p <- p + scale_y_continuous(breaks = seq(-180, 180, by = 30))
  p <- p + theme_minimal_grid() + coord_sf(crs = CRS_robinson_150)
  p <- p + geom_text_repel(data = combined_data[type == "modern", ], aes(x = longitude, y = latitude, label = gsub("SGDP_(.*)", "\\1\n(SGDP)", population)))
  p <- p + geom_text(data = combined_data[type == "archaic", ], aes(x = longitude + 1e6, y = latitude, label = population), hjust = 0)
  p <- p + theme(panel.grid.major = element_line(colour = "lightgrey", size = 0.25, linetype = 3))
  p <- p + theme(legend.position = "bottom") + labs(title = paste0("Allele frequency: ", myVariantID), subtitle = paste0("source: ", Introgressed_from, " - ", Introgression_scenario, "\n", "allele match: ", allele_match, "\n", "target: ", predicted_Target)) + xlab("") + ylab("")


  # emVars
  file_name <- sprintf("%s/emVar_plots/all_conditions/variant_%s_01d_scatterpie_plots_detail_B37.pdf", FIGURE_EMVAR_DIR, rsID_print)
  pdf(file_name, height = 6, width = 7)
  print(p)
  dev.off()
  for (SUBDIR in c("NS_only", "HepG2_only", "A549_only", "K562_only")) {
    target_file <- sprintf("%s/emVar_plots/%s/variant_%s_02d_scatterpie_plots_detail_B37.pdf", FIGURE_EMVAR_DIR, SUBDIR, rsID_print)
    system(sprintf("cp %s %s ", file_name, target_file))
  }

  if (myPosID %chin% response_emVars[, unique(posID)]) {
    # response emVars
    file_name <- sprintf("%s/response_emVar_plots/all_conditions/variant_%s_02d_scatterpie_plots_detail_B37.pdf", FIGURE_EMVARDIFF_DIR, rsID_print)
    pdf(file_name, height = 6, width = 7)
    print(p)
    dev.off()

    for (CELLINE in c("HepG2", "A549", "K562")) {
      SUBDIR <- sprintf("%s_only", CELLINE)
      if (CELLINE %in% response_emVars[posID == myPosID, celline]) {
        target_file <- sprintf("%s/response_emVar_plots/%s/variant_%s_02d_scatterpie_plots_detail_B37.pdf", FIGURE_EMVARDIFF_DIR, SUBDIR, rsID_print)
        system(sprintf("cp %s %s ", file_name, target_file))
      }
    }
  }

  if (myPosID %chin% celltype_emVars[, unique(posID)]) {
    # celltype emVars
    file_name <- sprintf("%s/celltype_emVar_plots/all_conditions/variant_%s_02d_scatterpie_plots_detail_B37.pdf", FIGURE_EMVARDIFF_DIR, rsID_print)
    pdf(file_name, height = 6, width = 7)
    print(p)
    dev.off()

    if (any(celltype_emVars[posID == myPosID, grepl("NS", group1_labels) & grepl("NS", group2_labels)])) {
      target_file <- sprintf("%s/celltype_emVar_plots/NS_only/variant_%s_02d_scatterpie_plots_detail_B37.pdf", FIGURE_EMVARDIFF_DIR, rsID_print)
      system(sprintf("cp %s %s ", file_name, target_file))
    }
  }
}

cat("All done\n")
q("no")



# # world_map <- as.data.table(map_data('world'))[region != "Antarctica",]

# # # Dessiner la carte du monde
#  p <- ggplot(world_map, aes(x = long, y = lat, group = group, map_id=region))
#  p <- p + geom_map(fill='lightgrey', colour = "white", linewidth=0.1, map=world_map)
#  p <- p + geom_scatterpie(
#      aes(x = longitude, y = latitude, group = population, r = 5),
#      data = combined_data,
#      cols = c("derived.freq", "ancestral.freq"),
#      color = NA
#    )
# # p <- p + theme_void()
# # Afficher la carte avec les pie charts

# pdf(sprintf("%s/test_scatterpie_plot.pdf", FIGURE_DIR_MAPS), height = 5, width = 5.4)
# print(p)
# dev.off()
