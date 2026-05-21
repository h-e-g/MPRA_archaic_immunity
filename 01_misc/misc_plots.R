################################################################################
# Useful functions
.libPaths(c("/pasteur/helix/projects/evo_immuno_pop/single_cell/resources/R_libs/4.1.0/", "/pasteur/appa/homes/mrotival/R/x86_64-pc-linux-gnu-library/4.0"))
shhh <- suppressPackageStartupMessages



cat('loading libraries\n')
options(width=200)
shhh(library(ggplot2))
shhh(library(ggrastr))
shhh(library(GGally))
shhh(library(viridis))
shhh(library(scales))
shhh(library(data.table))
shhh(library(dplyr))
shhh(library(tictoc))
shhh(library(stringr))
shhh(library(tidyr))
shhh(library(DescTools))
shhh(library(ggrepel))
shhh(library(corrplot))

# library(extrafont)
# # Import system fonts (do this only once)
# font_import(pattern = "Helvetica", prompt = FALSE)
# loadfonts(device = "pdf") 

# library(showtext)
# font_add("Helvetica", regular = sprintf('%s/data/fonts/Helvetica.ttf',MPRA_DIR))
# font_add("Comic", regular = sprintf('%s/data/fonts/Comic Sans MS.ttf',MPRA_DIR))

# showtext_auto()
# ggplot() <- theme_minimal(base_family = "Helvetica")

cat('libraries loaded\n')

# custom plot theme
theme_plot <- function(lpos = "bottom", rotate.x = 0, fontsize=7) {
  ptheme <- theme_bw() +
    theme(
      panel.grid = element_blank(),
      legend.title = element_blank(),
      #text = element_text(size = fontsize, family = "Comic"),
			text = element_text(size = fontsize, family = "sans"),
      legend.position = lpos,
      legend.background = element_rect(fill = "transparent"),
      panel.background = element_rect(fill = "transparent", colour = NA),
      plot.background = element_rect(fill = "transparent", colour = NA),
      strip.background = element_rect(fill = "transparent", colour = NA),
      panel.spacing = unit(0, "pt"),
      strip.text = element_text(color = "black")
    )
  if (rotate.x > 0) {
    ptheme <- ptheme + theme(axis.text.x = element_text(angle = rotate.x, hjust = 1, vjust = 0.5))
  }
  ptheme
}
# https://stackoverflow.com/questions/41631806/change-facet-label-text-and-background-colour

# function to add names to a vector on the fly.
add_names <- function(vector, name) {
  names(vector) <- name
  return(vector)
}

# function to shift legend to empty panel

shift_legend <- function(p) {
  # check if p is a valid object
  if (!"gtable" %in% class(p)) {
    if ("ggplot" %in% class(p)) {
      gp <- ggplotGrob(p) # convert to grob
    } else {
      message("This is neither a ggplot object nor a grob generated from ggplotGrob. Returning original plot.")
      return(p)
    }
  } else {
    gp <- p
  }

  # check for unfilled facet panels
  facet.panels <- grep("^panel", gp[["layout"]][["name"]])
  empty.facet.panels <- sapply(facet.panels, function(i) "zeroGrob" %in% class(gp[["grobs"]][[i]]))
  empty.facet.panels <- facet.panels[empty.facet.panels]
  if (length(empty.facet.panels) == 0) {
    message("There are no unfilled facet panels to shift legend into. Returning original plot.")
    return(p)
  }

  # establish extent of unfilled facet panels (including any axis cells in between)
  empty.facet.panels <- gp[["layout"]][empty.facet.panels, ]
  empty.facet.panels <- list(
    min(empty.facet.panels[["t"]]), min(empty.facet.panels[["l"]]),
    max(empty.facet.panels[["b"]]), max(empty.facet.panels[["r"]])
  )
  names(empty.facet.panels) <- c("t", "l", "b", "r")

  # extract legend & copy over to location of unfilled facet panels
  guide.grob <- which(gp[["layout"]][["name"]] == "guide-box")
  if (length(guide.grob) == 0) {
    message("There is no legend present. Returning original plot.")
    return(p)
  }
  gp <- gtable_add_grob(
    x = gp,
    grobs = gp[["grobs"]][[guide.grob]],
    t = empty.facet.panels[["t"]],
    l = empty.facet.panels[["l"]],
    b = empty.facet.panels[["b"]],
    r = empty.facet.panels[["r"]],
    name = "new-guide-box"
  )

  # squash the original guide box's row / column (whichever applicable)
  # & empty its cell
  guide.grob <- gp[["layout"]][guide.grob, ]
  if (guide.grob[["l"]] == guide.grob[["r"]]) {
    gp <- gtable_squash_cols(gp, cols = guide.grob[["l"]])
  }
  if (guide.grob[["t"]] == guide.grob[["b"]]) {
    gp <- gtable_squash_rows(gp, rows = guide.grob[["t"]])
  }
  gp <- gtable_remove_grobs(gp, "guide-box")

  return(gp)
}

################################################################################
# Color definitions

TE_Class_colors <- c("LINE" = "#B5DBCD", "SINE" = "#DB979D", "LTR" = "#97BBD6", "DNA" = "#EDB09B")
# condition colors, order
condition_color <- c("NS" = "#888888", "IAV" = "#6699CC", "COV" = "#DD5555")
condition_order <- c("NS", "COV", "IAV")

# condition and time point colors, order (time course)
cond_tp_color <- c("T0_T0" = "#888888", "NS_T6" = "#888888", "NS_T24" = "#888888", "IAV_T6" = "#6699CC", "IAV_T24" = "#254d74", "COV_T6" = "#DD5555", "COV_T24" = "#8f0a0a")
cond_tp_order <- c("T0_T0", "NS_T6", "NS_T24", "COV_T6", "COV_T24", "IAV_T6", "IAV_T24")

# lineage colors, order
lineage_order <- c("MONO", "B", "T.CD4", "T.CD8", "NK")
lineage_color <- c("B" = "#9281DE", "MONO" = "#B7637E", "T.CD4" = "#005274", "T.CD8" = "#048c41", "NK" = "#7D9E00")

# celltype colors, order
celltype_color <- c(
  "MONO.CD14" = "#B7637E", "MONO.CD16" = "#8F568A", "pDC" = "#8E6B00", "cDC" = "#E7b315", "MONO.CD14.INFECTED" = "#B9364B",
  "B.N.K" = "#9281DE", "B.N.L" = "#a681de", "B.M.K" = "#6045d9", "B.M.L" = "#8045d9", "Plasmablast" = "#6801A4",
  "T.CD4.N" = "#169FD8", "T.CD4.E" = "#005274", "T.Reg" = "#03C2C0", "T.gd" = "#ffce73",
  "T.CD8.N" = "#00BB54", "T.CD8.CM.EM" = "#0EAD20", "T.CD8.EMRA" = "#048c41", "ILC" = "#3F4F00", "MAIT" = "#004A21",
  "NK.CD56dim" = "#B7DC2A", "NK.CD56brt" = "#7D9E00", "NK.M.LIKE" = "#63C425"
)
celltype_order <- c(
  "MONO.CD14", "MONO.CD16", "MONO.CD14.INFECTED", "cDC", "pDC",
  "B.N.K", "B.N.L", "B.M.K", "B.M.L", "Plasmablast",
  "T.CD4.N", "T.CD4.E", "T.Reg", "T.gd",
  "T.CD8.N", "T.CD8.CM.EM", "T.CD8.EMRA", "ILC", "MAIT",
  "NK.CD56dim", "NK.CD56brt", "NK.M.LIKE"
)

# population colors, order
population_color <- c("YRI" = "#008000", "CEU" = "#eba206", "CHS" = "#71458d")
population_order <- c("YRI", "CEU", "CHS")

population_1kg_color <- c("YRI" = "#008000", "CEU" = "#eba206", "CHS" = "#71458d")
population_1kg_order <- c("YRI", "CEU", "CHS")

# 1kG population colors, order
color_1kG <- c(
  'FIN' = "#54278F", 'GBR' = "#756BB1", 'CEU' = "#9E9AC8", 'TSI' = "#BCBDDC", 'IBS' = "#DADAEB",
  'MXL' = "#BDD7E7", 'PUR' = "#6BAED6", 'CLM' = "#3182BD", 'PEL' = "#08519C", 
  'KHV' = "#A63603", 'JPT' = "#E6550D", 'CHB' = "#FD8D3C", 'CHS' = "#FDAE6B", 'CDX' = "#FDD0A2",
  'GIH' = "#FCBBA1", 'PJL' = "#FC9272", 'BEB' = "#FB6A4A", 'STU' = "#DE2D26", 'ITU' = "#A50F15", 
  'ASW' = "#E5F5E0", 'LWK' = "#C7E9C0", 'GWD' = "#A1D99B", 'MSL' = "#74C476", 'ESN' = "#41AB5D", 'YRI' = "#238B45", 'ACB' = "#005A32",
  'panTro2' = "yellow"
)

color_populations_MPRA=c('YRI' = "#41AB5D", 'AFRICANS' = "#005A32",
  'GBR' = "#08519C", 'IBS' = "#6BAED6", 'WESTEURASIANS' = "#3182BD",
  'STU' = "#54278F", 'PJL' = "#9E9AC8",
  'JPT' = "#A50F15", 'CHB' = "#F87E7E",'EASTEURASIANS' = "#FB4A4A",
  'AGTA' = "#FD8D3C",'PAPUAN' = "#A63603")

# COVID-19 susceptibility and severity colors, order
COVID_color <- c("critical" = "#DD5555", "hospitalized" = "#E68686", "reported" = "#F1BABA")
COVID_order <- c("reported", "hospitalized", "critical")

# archaic colors, order
archaic_color <- c("AMH" = grey(0.8), "DENI" = "#80CDC1", "NEAND" = "#FDAE61", "ARCHAIC" = "#FEE08B", "UNKNOWN_ARCHAIC" = grey(0.5))
color_archaic_MPRA <- c("Denisova" = "#80CDC1", "Vindija" = "#FDAE61", "Vindija/Denisova" = "#FEE08B", "Archaic (unknown)" = grey(0.5))

# add 12 entry divergent colorblind_palette # (taken from the rcartocolor package: safe palette)
safe_colorblind_palette <- c(
  "#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499",
  "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888"
)


quality_5levels <- c("GOOD" = "#264653", "OK" = "#2a9d8f", "AVERAGE" = "#e9c46a", "POOR" = "#f4a261", "BAD" = "#e76f51")
reddish <- "#C59C9D"
blueish <- "#529DD9"
yellowish <- "#E0BA3E"
brownish <- "#A79987"
whitish <- "#E2EDEE"
# brewer.pal(n,'BrBg')

ChromHMM_colors <- fread(sprintf("%s/data/ChromHMM_colors.txt", MPRA_DIR))
ChromHMM_colors[, R := as.numeric(str_split(V5, ",", simplify = T)[, 1])]
ChromHMM_colors[, G := as.numeric(str_split(V5, ",", simplify = T)[, 2])]
ChromHMM_colors[, B := as.numeric(str_split(V5, ",", simplify = T)[, 3])]
ChromHMM_colors[, HEX := rgb(R / 255, G / 255, B / 255)]
ChromHMM_colors[V2=='Enh',HEX:='#E6E600']


mergeCols <- function(col1, col2, pct2 = 0.5) {
  colorRampPalette(c(col1, col2))(100)[round(pct2 * 100)]
}

color_type_simplified=c(other='#AAAAAA',Promoter='#DD5555',scrambled="#08519C")


color_setup_simplified_norep <- c( "A549_NS" = "#36c54c",
                  "A549_IAV" = "#a7e563",
                   "A549_SARS" = "#00a877",
                   "A549_TNFa" = "#007f5c",
                  "HepG2_NS" = "#20AEDD",
                  "HepG2_IFNA2b"= "#9aa8e8",
                  "HepG2_DEX"= "#9acae8",
                  "HepG2_TNFa"= "#4d8cea",
                  "K562_NS"= "#f2a465",
                  "K562_IFNA2b"= "#ff8c69",
                  "K562_DEX"= "#ffd212",
                  "K562_TNFa"= "#e25822")


color_setup_Lisa <- c( "A549_NS" = "#36c54c",
                  "A549_IAV" = "#a7e563",
                   "A549_SARS" = "#00a877",
                   "A549_TNFa" = "#007f5c",
                  "HepG2_NS" = "#20AEDD",
                  "HepG2_IFNA2b"= "#9aa8e8",
                  "HepG2_DEX"= "#9AE7E8",
                  "HepG2_TNFa"= "#4d8cea",
                  "K562_NS"= "#f2a465",
                  "K562_IFNA2b"= "#ff8c69",
                  "K562_DEX"= "#ffa812",
                  "K562_TNFa"= "#e25822")


color_setup_MAX <- c( "A549_NS" = "#8dd296", 
                  "A549_IAV" = "#c7da87",
                   "A549_SARS" = "#88dbbb", 
                   "A549_TNFa" = "#58bf83",
                  "HepG2_NS" = "#7cacd6", 
                  "HepG2_IFNA2b"= "#ab9fe6",
                  "HepG2_DEX"= "#a3ddeb", 
                  "HepG2_TNFa"= "#5e7ebd",
                  "K562_NS"= "#f2a465",
                  "K562_IFNA2b"= "#f2a7ac",
                  "K562_DEX"= "#f6d360",
                  "K562_TNFa"= "#f57c6e"
                  )



# color_setup <- c("Calu3_NS_R1" = "#00cc66", 
#                   "HepG2_NS_R1" = "#6600cc", 
#                   "HepG2_NS_R2" = "#7610dc", 
#                   "HepG2_NS_R3" = "#8620ec",
#                   "HepG2_NS_R4" = "#9632ff",
#                   "HepG2_NS_R5" = "#a543ff",
#                   "HepG2_IFNA2b_R1"= "#a71ed9", 
#                   "HepG2_IFNA2b_R2"= "#b72ee9",
#                   "HepG2_IFNA2b_R3"= "#c73ef9",
#                   "HepG2_IFNA2b_R4"= "#d73fff",
#                   "HepG2_DEX_R1"= "#9796DE", 
#                   "HepG2_DEX_R2"= "#A7A6EE", 
#                   "HepG2_DEX_R3"= "#B7B6FE", 
#                   "HepG2_TNFa_R1"= "#1d1aee", 
#                   "HepG2_TNFa_R2"= "#2d2bee",
#                   "HepG2_TNFa_R3"= "#3d3cfe",
#                   "A549_NS_R1" = "#62cf2f", 
#                   "A549-ACE2_NS_R1" = "#81C16E", 
#                   "A549-ACE2_NS_R2-P3" = "#91D17E", 
#                   "A549-ACE2_NS_R3-P3" = "#A1E18E", 
#                   "A549-ACE2_IAV_R1-P3" = "#9BDAB3", 
#                   "A549-ACE2_IAV_R2" = "#ABEAC3", 
#                   "A549-ACE2_IAV_R3" = "#BCFAD3", 
#                   "A549-ACE2_SARS_R1-P3" = "#9B9D2C", 
#                   "A549-ACE2_SARS_R2-P3" = "#ABAD3C", 
#                   "A549-ACE2_SARS_R3-P3" = "#BBBD4C", 
#                   "A549-ACE2_TNFa_R1" = "#2b9a2b", 
#                   "A549-ACE2_TNFa_R2" = "#3b9b3b",
#                   "A549-ACE2_TNFa_R3" = "#4b9c4b",
#                   "K562_NS_R1"= "#b77718",
#                   "K562_NS_R2"= "#c78728",
#                   "K562_NS_R3"= "#d79738",
#                   "K562_NS_R4"= "#e8a648",
#                   "K562_NS_R5"= "#f9b658",
#                   "K562_IFNA2b_R1"= "#CC5410",
#                   "K562_IFNA2b_R2"= "#DC6420",
#                   "K562_IFNA2b_R3"= "#EC7430",
#                   "K562_IFNA2b_R4"= "#FC8530",
#                   "K562_DEX_R1"= "#81545B",
#                   "K562_DEX_R2"= "#91646B",
#                   "K562_DEX_R3"= "#A1747B",
#                   "K562_TNFa_R1"= "#613f0d",
#                   "K562_TNFa_R2"= "#723f1d",
#                   "K562_TNFa_R3"= "#833f2d"
#                   )

# color_celline_desaturated <- c("HepG2" = "#AE89C4", "A549" = "#81C16E", "K562" = "#BF9454")
# #color_celline <- c("HepG2" = "#6600cc66", "Calu3" = "#00cc6666", "A549" = "#62cf2f", "K562" = "#b77718")
# color_celline_saturated <- c("HepG2" = "#6600cc66", "A549" = "#62cf2f", "K562" = "#b77718")
# color_celline <- color_celline_desaturated

# color_celltype_saturated  <- c("HepG2" = "#6600cc66", "A549-ACE2" = "#62cf2f", "A549" = "#62cf2f", "K562" = "#b77718")
# color_celltype_desaturated  <- c("HepG2" = "#AE89C4", "A549-ACE2" = "#81C16E", "A549" = "#62cf2f", "K562" = "#BF9454")
# color_celltype <- color_celltype_desaturated

color_celline <- color_celltype <- color_setup_simplified_norep[grep('_NS',names(color_setup_simplified_norep))]
names(color_celline)=gsub('_NS','',names(color_celline))
names(color_celltype)=gsub('A549','A549-ACE2',names(color_celline))


color_celline_2levels <- c(color_celline[c("HepG2", "A549", "K562")], sapply(color_celline[c("HepG2", "A549", "K562")], mergeCols, "black"))
names(color_celline_2levels) <- paste(rep(c("HepG2", "A549", "K562"), 2), rep(c("weak", "strong"), e = 3))

color_celline_2levels_light <- c(sapply(color_celline[c("HepG2", "A549", "K562")], mergeCols, "white"), color_celline[c("HepG2", "A549", "K562")] )
names(color_celline_2levels_light) <- paste(rep(c("HepG2", "A549", "K562"), 2), rep(c("weak", "strong"), e = 3))


if(exists('condition_summary_reps')){
	color_setup <- color_setup_simplified_norep[condition_summary_reps$COND_ID]
	names(color_setup) <- condition_summary_reps$analysis_name

	color_setup_norep=color_setup
	names(color_setup_norep)=gsub('^(.*_.*)_(.*)$','\\1',names(color_setup))
	color_setup_norep=color_setup_norep[!duplicated(names(color_setup_norep))]

	color_setup_simplified=color_setup[!names(color_setup) %in% c('Calu3_NS_R1','A549_NS_R1')]
	names(color_setup_simplified)=gsub("-ACE2","",names(color_setup_simplified))

	# color_setup_simplified_norep=color_setup_norep[!names(color_setup_norep) %in% c('Calu3_NS','A549_NS')]
	# names(color_setup_simplified_norep)=gsub("-ACE2","",names(color_setup_simplified_norep))
}
condition_order <- c("NS", "IAV", "SARS","IFNA2b", "DEX", "TNFa")
condition_summary=condition_summary[order(celline,factor(condition,levels=condition_order))]
#color_library=setNames(color_setup[LIB_summary$SETUP_ID], LIB_summary$library)
color_readclass <- c("discordant barcodes, unknown" = "#DD5555", "discordant barcodes, partially known" = "#eba206", "unknown BC, low complexity" = "#AA4499", "unknown BC" = "#332288", "known BC, low complexity UMI" = "#117733", "known BC" = "#44AA99")


CRE_class_colors <- c("strong enhancer" = "#0f50b8", "weak enhancer" = "#4287f5", "inactive"="grey","weak silencer" = "#0fb88e", "strong silencer" = "#047055")

color_TRUEFALSE <- c("FALSE" = "#000000", "TRUE" = "#DD5555")
alpha_TRUEFALSE <- c("FALSE" = 0.5, "TRUE" = 1)


# color_setup <- c( "A549-ACE2_NS_R1" = "#91f7c9", 
#                   "A549-ACE2_IAV_R1-P3" = "#84C3B7",
#                    "A549-ACE2_SARS_R1-P3" = "#88d8db", 
#                    "A549-ACE2_TNFa_R1" = "#58bf83",
#                   "HepG2_NS_R1" = "#20AEDD", 
#                   "HepG2_IFNA2b_R1"= "#4C9BE6",
#                   "HepG2_DEX_R1"= "#9ACAE8", 
#                   "HepG2_TNFa_R1"= "#b8aeeb",
#                   "K562_NS_R1"= "#f2a465",
#                   "K562_IFNA2b_R1"= "#f2a7da",
#                   "K562_DEX_R1"= "#f8e49a",
#                   "K562_TNFa_R1"= "#f57c6e"
#                     )



color_setup_MAX <- c( "A549-ACE2_NS_R1" = "#8dd296", 
                  "A549-ACE2_IAV_R1-P3" = "#c7da87",
                   "A549-ACE2_SARS_R1-P3" = "#88dbbb", 
                   "A549-ACE2_TNFa_R1" = "#58bf83",
                  "HepG2_NS_R1" = "#7cacd6", 
                  "HepG2_IFNA2b_R1"= "#ab9fe6",
                  "HepG2_DEX_R1"= "#a3ddeb", 
                  "HepG2_TNFa_R1"= "#5e7ebd",
                  "K562_NS_R1"= "#f2a465",
                  "K562_IFNA2b_R1"= "#f2a7ac",
                  "K562_DEX_R1"= "#f6d360",
                  "K562_TNFa_R1"= "#f57c6e"
                  )

                  # "A549-ACE2_NS_R1" = "#50d664ff",
                  # "A549-ACE2_IAV_R1-P3" = "#a7e563",
                  # "A549-ACE2_SARS_R1-P3" = "#6be1beff",

# color_setup <- c("Calu3_NS_R1" = "#00cc66", 
#                   "HepG2_NS_R1" = "#20AEDD", 
#                   "HepG2_NS_R2" = "#30BEED", 
#                   "HepG2_NS_R3" = "#40CEFD",
#                   "HepG2_NS_R4" = "#50DEFF",
#                   "HepG2_NS_R5" = "#60EEFF",
#                   "HepG2_IFNA2b_R1"= "#4C9BE6", 
#                   "HepG2_IFNA2b_R2"= "#5CABF6",
#                   "HepG2_IFNA2b_R3"= "#6CBBFF",
#                   "HepG2_IFNA2b_R4"= "#7CCBFF",
#                   "HepG2_DEX_R1"= "#9ACAE8", 
#                   "HepG2_DEX_R2"= "#9ACAE8", 
#                   "HepG2_DEX_R3"= "#9ACAE8", 
#                   "HepG2_TNFa_R1"= "#b8aeeb", 
#                   "HepG2_TNFa_R2"= "#2d2bee",
#                   "HepG2_TNFa_R3"= "#3d3cfe",
#                   "A549_NS_R1" = "#91f7c9", 
#                   "A549-ACE2_NS_R1" = "#91f7c9", 
#                   "A549-ACE2_NS_R2-P3" = "#91D17E", 
#                   "A549-ACE2_NS_R3-P3" = "#A1E18E", 
#                   "A549-ACE2_IAV_R1-P3" = "#84c3b7", 
#                   "A549-ACE2_IAV_R2" = "#ABEAC3", 
#                   "A549-ACE2_IAV_R3" = "#BCFAD3", 
#                   "A549-ACE2_SARS_R1-P3" = "#88d8db", 
#                   "A549-ACE2_SARS_R2-P3" = "#ABAD3C", 
#                   "A549-ACE2_SARS_R3-P3" = "#BBBD4C", 
#                   "A549-ACE2_TNFa_R1" = "#58bf83", 
#                   "A549-ACE2_TNFa_R2" = "#3b9b3b",
#                   "A549-ACE2_TNFa_R3" = "#4b9c4b",
#                   "K562_NS_R1"= "#f2a465",
#                   "K562_NS_R2"= "#c78728",
#                   "K562_NS_R3"= "#d79738",
#                   "K562_NS_R4"= "#e8a648",
#                   "K562_NS_R5"= "#f9b658",
#                   "K562_IFNA2b_R1"= "#f2a7da",
#                   "K562_IFNA2b_R2"= "#DC6420",
#                   "K562_IFNA2b_R3"= "#EC7430",
#                   "K562_IFNA2b_R4"= "#FC8530",
#                   "K562_DEX_R1"= "#f8e49a",
#                   "K562_DEX_R2"= "#fff6ae",
#                   "K562_DEX_R3"= "#ffffbe",
#                   "K562_TNFa_R1"= "#f57c6e",
#                   "K562_TNFa_R2"= "#723f1d",
#                   "K562_TNFa_R3"= "#833f2d"
#                   )


# pdf('users/mrotival/Maxime_generated_text_maestro.pdf')
# ggplot(data=data.frame(x=rnorm(100),y=rnorm(100)),aes(x=x,y=y))+geom_point()+ylab('Helvetica or arial ?')+theme(text = element_text(family = "sans"));
# dev.off()
