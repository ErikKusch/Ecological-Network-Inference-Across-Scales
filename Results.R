#' ####################################################################### #
#' PROJECT: [Data Simplification for Network Inference across Scales] 
#' CONTENTS: 
#'  - Network comparison and Plotting
#'  DEPENDENCIES:
#'  - 0 - Preamble.R
#'  - X - Functions_Data.R
#'  - X - Functions_Plotting.R
#'  - Scripts that need to be run prior:
#'    + YFDP.R 
#'    + FIA.R
#' AUTHOR: [Erik Kusch]
#' ####################################################################### #

# PREAMBLE =================================================================
rm(list=ls())
set.seed(42)

## Sourcing ----------------------------------------------------------------
source("0 - Preamble.R")
source("X - Functions_Data.R")
source("X - Functions_Plotting.R")

# DATA =====================================================================
message("############ LOADING DATA")

## LOADING -----------------------------------------------------------------
YFDP_ls <- Load.Results(Dir = DirEx.YFDP)
Region_ls <- Load.Results(Dir = DirEx.Region)

## MANIPULATION ------------------------------------------------------------
### SPECIES IDENTITIES ####
YFDP_spec <- ReadOut.Species(YFDP_ls) # species which have been analysed
Region_spec <- ReadOut.Species(Region_ls) # species which have been analysed
Shared_spec <- Reduce(intersect, list(YFDP_spec, Region_spec)) # shared species between all analyses

### FLATTENING LISTS ####
YFDP_ls <- Flatten.Lists(YFDP_ls)
Region_ls <- Flatten.Lists(Region_ls)
Region_ls <- BiomeNames.List(Region_ls)
Macro_ls <- Region_ls[grep(names(Region_ls), pattern = "Temperate Conifer")]
Region_ls <- Region_ls[grep(names(Region_ls), pattern = "Yosemite")]

### FULL NETWORK DATA FRAMES ###
message("############ LEGEND ###")
Species_vec <- unique(c(YFDP_spec, 
                 Region_ls[[1]]$Partner1,
                 Macro_ls[[1]]$Partner1))
Species_ls <- base::strsplit(x = as.character(Species_vec), split = " ")
SpeciesLegend <- data.frame(Abbr = paste0(substr(unlist(lapply(Species_ls, "[[", 1)), start = 1, stop = 2),
                                          substr(unlist(lapply(Species_ls, "[[", 2)), start = 1, stop = 2)),
                            Full = Species_vec
)
sink(file = file.path(Dir.Exports, "SpeciesLegend.txt"))
print(SpeciesLegend)
sink()

YFDP_full <- Eff.Data.Frame(List_ls =Limit.Lists(YFDP_ls, YFDP_spec))
YFDP_full$Partner1 <- SpeciesLegend$Abbr[match(YFDP_full$Partner1, SpeciesLegend$Full)]
YFDP_full$Partner2 <- SpeciesLegend$Abbr[match(YFDP_full$Partner2, SpeciesLegend$Full)]

Region_full <- Eff.Data.Frame(List_ls = Limit.Lists(Region_ls, 
                                                    unique(Region_ls[[1]]$Partner1)))
Region_full$Partner1 <- SpeciesLegend$Abbr[match(Region_full$Partner1, SpeciesLegend$Full)]
Region_full$Partner2 <- SpeciesLegend$Abbr[match(Region_full$Partner2, SpeciesLegend$Full)]

Macro_full <- Eff.Data.Frame(List_ls = Limit.Lists(Macro_ls,                                  
                                                   unique(Macro_ls[[1]]$Partner1)))
Macro_full$Partner1 <- SpeciesLegend$Abbr[match(Macro_full$Partner1, SpeciesLegend$Full)]
Macro_full$Partner2 <- SpeciesLegend$Abbr[match(Macro_full$Partner2, SpeciesLegend$Full)]

Shared_abbr <- SpeciesLegend$Abbr[match(Shared_spec, SpeciesLegend$Full)]

# FIGURE 1, S3 & S4 - MODEL SPECIFICATION DRIVES NETWORK INFERENCE =========
message("############ FIGURE 1, S3 & S4")
## DATA --------------------------------------------------------------------
Fig1S3S4_ls <- list(Figure1 = YFDP_full,
                  S3 = Region_full,
                  S4 = Macro_full)
MATplots_ls <- as.list(rep(NA, 3))
names(MATplots_ls) <- c("Plot", "Region", "Macro")
## PLOT PRODUCTION ---------------------------------------------------------
for(i in 1:length(Fig1S3S4_ls)){
  plot_df <- Fig1S3S4_ls[[i]]
  print(names(Fig1S3S4_ls)[i])
  
  if(i == 1){Treatment <- "Pre-Fire"
  Pos1 <- which(plot_df$Partner1 == "Conu" & plot_df$Partner2 == "Coco")
  first_row <- plot_df[Pos1, ]
  Pos2 <- which(plot_df$Partner2 == "Conu" & plot_df$Partner1 == "Coco")
  second_row <- plot_df[Pos2, ]
  plot_df[Pos1, ] <- second_row
  plot_df[Pos2, ] <- first_row
  Pos1 <- which(plot_df$Partner1 == "Cose" & plot_df$Partner2 == "Coco")
  first_row <- plot_df[Pos1, ]
  Pos2 <- which(plot_df$Partner2 == "Cose" & plot_df$Partner1 == "Coco")
  second_row <- plot_df[Pos2, ]
  plot_df[Pos1, ] <- second_row
  plot_df[Pos2, ] <- first_row
  }
  if(i == 2){Treatment <- "Yosemite"
  # Check_sp <- Reg_sp
  }
  if(i == 3){Treatment <- "Temperate Conifer Forests"
  Pos1 <- which(plot_df$Partner1 == "Pien" & plot_df$Partner2 == "Pico")
  first_row <- plot_df[Pos1, ]
  Pos2 <- which(plot_df$Partner2 == "Pien" & plot_df$Partner1 == "Pico")
  second_row <- plot_df[Pos2, ]
  plot_df[Pos1, ] <- second_row
  plot_df[Pos2, ] <- first_row
  # Check_sp <- Macro_sp
  }
  
  Fig1S3S4_ls[[i]] <- plot_df
  bordercol <- function(testx){
    testx$Shared <- 1
    testx$Shared[which(testx$Partner1 %in% Shared_abbr & testx$Partner2 %in% Shared_abbr)] <- 2
    testx$Shared[is.na(testx$Value)] <- 1
    testx$Shared <- as.factor(testx$Shared)
    # testx$Partner1 <- factor(testx$Partner1, levels = c(Shared_abbr, 
    #                                           unique(testx$Partner1[testx$Partner1 %nin% Shared_abbr])))
    # testx$Partner2 <- factor(testx$Partner2, levels = c(Shared_abbr, 
    #                                                     unique(testx$Partner2[testx$Partner2 %nin% Shared_abbr])))
    testx
  }
  
  ### COCCUR ----
  COCC_dat <- Plot.NetMat(data = plot_df, method = "COCCUR", datareturn = TRUE) 
  COCC_dat <- COCC_dat[COCC_dat$Condition == Treatment , ]
  COCC_dat <- bordercol(COCC_dat)
  
  COCC_plot <- ggplot(COCC_dat, aes(x = Partner2, y = Partner1, fill = Value, col = Shared)) +
    geom_tile(lwd = 1,
              linetype = 1,
              width=0.92, height=0.92) + 
    geom_point(aes(shape = Sig), col = "black", size = 2) +
    scale_shape_manual(values=c(32, 15), na.translate = FALSE, name = "Significance") +  
    scale_color_manual(values=c("white", "black"), na.translate = TRUE, name = "", guide = "none") +  
    guides(shape = FALSE) + 
    coord_fixed() + 
    guides(fill = guide_colourbar(barwidth = 2,
                                  barheight = 15,
                                  title = "Association")) + 
    ## axes
    theme_bw() + 
    theme(axis.text.x=element_text(angle = -40, hjust = 0),
          text = element_text(size=15)) +
    scale_fill_gradient2(low = "#ba00bd", high = "#d9db3b") + 
    labs(x = "", y = "")
  
  Empty_plot <- ggplot(COCC_dat, aes(x = Partner2, y = Partner1, fill = Value)) + theme_void()
  
  ### NETASSOC ----
  NETA_dat <- Plot.NetMat(data = plot_df, method = "NETASSOC", datareturn = TRUE) 
  NETA_dat <- NETA_dat[NETA_dat$Condition == Treatment, ]
  NETA_dat <- bordercol(NETA_dat)
  
  NETA_plot <- ggplot(NETA_dat, aes(x = Partner2, y = Partner1, fill = Value, col = Shared)) +
    geom_tile(lwd = 1,
              linetype = 1,
              width=0.92, height=0.92) + 
    geom_point(aes(shape = Sig), col = "black", size = 2) +
    scale_shape_manual(values=c(32, 15), na.translate = FALSE, name = "Significance") +  
    scale_color_manual(values=c("white", "black"), na.translate = TRUE, name = "", guide = "none") +  
    guides(shape = FALSE) + 
    coord_fixed() + 
    guides(fill = guide_colourbar(barwidth = 2,
                                  barheight = 15,
                                  title = "Association")) + 
    ## axes
    theme_bw() + 
    theme(axis.text.x=element_text(angle = -40, hjust = 0),
          text = element_text(size=15)) +
    scale_fill_gradient2(low = "#ba00bd", high = "#d9db3b") + 
    labs(x = "", y = "")
  
  ### HMSC ----
  HMSC_dat <- Plot.NetMat(data = plot_df, method = "HMSC", datareturn = TRUE) 
  HMSC_dat <- HMSC_dat[HMSC_dat$Condition == Treatment, ]
  HMSC_dat <- bordercol(HMSC_dat)
  
  HMSC_dat$Model <- gsub(x = HMSC_dat$Model, pattern = "presence_absence", replacement = "Presence & Absence")
  HMSC_dat$Model <- gsub(x = HMSC_dat$Model, pattern = "abundance", replacement = "Abundance")
  HMSC_dat$Model <- gsub(x = HMSC_dat$Model, pattern = "diametre", replacement = "Performance")
  HMSC_dat$Model <- gsub(x = HMSC_dat$Model, pattern = "biomass", replacement = "Performance")
  
  HMSC_plot1 <- ggplot(HMSC_dat[HMSC_dat$Model == "Presence & Absence",], aes(x = Partner2, y = Partner1, fill = Value, col = Shared)) +
    geom_tile(lwd = 1,
              linetype = 1,
              width=0.92, height=0.92) + 
    geom_point(aes(shape = Sig), col = "black", size = 2) +
    scale_shape_manual(values=c(32, 15), na.translate = FALSE, name = "Significance") +  
    scale_color_manual(values=c("white", "black"), na.translate = TRUE, name = "", guide = "none") +  
    guides(shape = FALSE) + 
    coord_fixed() + 
    guides(fill = guide_colourbar(barwidth = 2,
                                  barheight = 15,
                                  title = "Association")) + 
    ## axes
    theme_bw() + 
    theme(axis.text.x=element_text(angle = -40, hjust = 0),
          text = element_text(size=15)) +
    scale_fill_gradient2(low = "#ba00bd", high = "#d9db3b") + 
    labs(x = "", y = "")
  
  HMSC_plot2 <- ggplot(HMSC_dat[HMSC_dat$Model == "Abundance",], aes(x = Partner2, y = Partner1, fill = Value, col = Shared)) +
    geom_tile(lwd = 1,
              linetype = 1,
              width=0.92, height=0.92) + 
    geom_point(aes(shape = Sig), col = "black", size = 2) +
    scale_shape_manual(values=c(32, 15), na.translate = FALSE, name = "Significance") +  
    scale_color_manual(values=c("white", "black"), na.translate = TRUE, name = "", guide = "none") +  
    guides(shape = FALSE) + 
    coord_fixed() + 
    guides(fill = guide_colourbar(barwidth = 2,
                                  barheight = 15,
                                  title = "Association")) + 
    ## axes
    theme_bw() + 
    theme(axis.text.x=element_text(angle = -40, hjust = 0),
          text = element_text(size=15)) +
    scale_fill_gradient2(low = "#ba00bd", high = "#d9db3b") + 
    labs(x = "", y = "")
  
  HMSC_plot3 <- ggplot(HMSC_dat[HMSC_dat$Model == "Performance",], aes(x = Partner2, y = Partner1, fill = Value, col = Shared)) +
    geom_tile(lwd = 1,
              linetype = 1,
              width=0.92, height=0.92) + 
    geom_point(aes(shape = Sig), col = "black", size = 2) +
    scale_shape_manual(values=c(32, 15), na.translate = FALSE, name = "Significance") +  
    scale_color_manual(values=c("white", "black"), na.translate = TRUE, name = "", guide = "none") +  
    guides(shape = FALSE) + 
    coord_fixed() + 
    guides(fill = guide_colourbar(barwidth = 2,
                                  barheight = 15,
                                  title = "Association")) + 
    ## axes
    theme_bw() + 
    theme(axis.text.x=element_text(angle = -40, hjust = 0),
          text = element_text(size=15)) + 
    scale_fill_gradient2(low = "#ba00bd", high = "#d9db3b") + 
    labs(x = "", y = "")
  
  ### IF_REM ----
  IFREM_dat <- Plot.NetMat(data = plot_df, method = "IF_REM", datareturn = TRUE) 
  IFREM_dat <- IFREM_dat[IFREM_dat$Condition == Treatment, ]
  IFREM_dat <- bordercol(IFREM_dat)
  IFREM_dat$Sig <- factor(IFREM_dat$Sig)
  
  IFREM_plot <- ggplot(IFREM_dat, aes(x = Partner2, y = Partner1, fill = Value, col = Shared)) +
    geom_tile(lwd = 1,
              linetype = 1,
              width=0.92, height=0.92) + 
    geom_point(aes(shape = Sig), col = "black", size = 2) +
    scale_shape_manual(values=c(32, 15), na.translate = FALSE, name = "Significance") +  
    scale_color_manual(values=c("white", "black"), na.translate = TRUE, name = "", guide = "none") +  
    guides(shape = FALSE) + 
    coord_fixed() + 
    guides(fill = guide_colourbar(barwidth = 2,
                                  barheight = 15,
                                  title = "Interaction")) + 
    ## axes
    theme_bw() + 
    theme(axis.text.x=element_text(angle = -40, hjust = 0),
          text = element_text(size=15)) + 
    scale_fill_gradient2(low = "#ba00bd", high = "#d9db3b") + 
    labs(x = "Subject", y = "Actor")
  
  ### IF_REM UNDIRECTED ----
  IFREM_ASSOC_dat <- Plot.NetMat(data = plot_df, method = "IF_REM_Assoc", datareturn = TRUE) 
  IFREM_ASSOC_dat <- IFREM_ASSOC_dat[IFREM_ASSOC_dat$Condition == Treatment, ]
  IFREM_ASSOC_dat <- bordercol(IFREM_ASSOC_dat)
  IFREM_ASSOC_dat$Sig[is.na(IFREM_ASSOC_dat$Sig)] <- 0
  IFREM_ASSOC_dat$Sig <- factor(IFREM_ASSOC_dat$Sig)
  
  IFREM_ASSOC_plot <- ggplot(IFREM_ASSOC_dat, aes(x = Partner2, y = Partner1, fill = Value, col = Shared)) +
    geom_tile(lwd = 1,
              linetype = 1,
              width=0.92, height=0.92) + 
    geom_point(aes(shape = Sig), col = "black", size = 2) +
    scale_shape_manual(values=c(32, 0, 15), na.translate = FALSE, name = "Significance") +  
    scale_color_manual(values=c("white", "black"), na.translate = TRUE, name = "", guide = "none") +  
    guides(shape = FALSE) + 
    coord_fixed() + 
    guides(fill = guide_colourbar(barwidth = 2,
                                  barheight = 15,
                                  title = "Association")) + 
    ## axes
    theme_bw() + 
    theme(axis.text.x=element_text(angle = -40, hjust = 0),
          text = element_text(size=15)) + 
    scale_fill_gradient2(low = "#ba00bd", high = "#d9db3b") + 
    labs(x = "", y = "")
  
  grid2 <-  plot_grid(plotlist = 
              list(
                plot_grid(plotlist = list(COCC_plot, NETA_plot, IFREM_plot), nrow = 1, 
                          labels = c("COOCUR", "NETASSOC", "NDD-RIM")),
                plot_grid(plotlist = list(HMSC_plot1, HMSC_plot2, HMSC_plot3), nrow = 1, 
                          labels = c("HMSC (Presence & Absence)", "HMSC (Abundance)", "HMSC (Performance)"))
              ),
            ncol = 1,
            labels = c("")
  )
  ggsave(grid2, filename = file.path(Dir.Exports, paste0(names(Fig1S3S4_ls)[i], ".jpg")), height = 30, width = 50, units = "cm")
  
  MATplots_ls[[i]] <- list(COCCUR = list(plot = COCC_plot, 
                                         data = COCC_dat),
                           NETASSOC = list(plot = NETA_plot,
                                           data = NETA_dat),
                           HMSC = list(plot = HMSC_plot3,
                                       data = HMSC_dat[HMSC_dat$Model == "Performance",]),
                           NDD_RIM = list(plot = IFREM_plot,
                                           data = IFREM_dat),
                           NDD_RIM_ASSOC = list(plot = IFREM_ASSOC_plot,
                                                data = IFREM_ASSOC_dat)
                           )
}

# FIGURE 2 - NETWORK INFERENCE IS NOT CONSISTENT ============================
message("############ FIGURE 2")
### DATA --------------------------------------------------------------------
# data_ls <- Fig1S3S4_ls
# 
# # YFDP_limited <- YFDP_ls[grep(names(YFDP_ls), pattern = "Pre-Fire")]
# # YFDP_spec <- unique(unlist(lapply(YFDP_limited, FUN = function(x){x[,1]})))
# # Region_limited <- Region_ls[c(grep(names(Region_ls), pattern = "Yosemite"), grep(names(Region_ls), pattern = "Temperate Conifer Forests"))]
# # Region_spec <- unique(unlist(lapply(Region_limited, FUN = function(x){x[,1]})))
# # Shared_spec <- Reduce(intersect, list(YFDP_spec, Region_spec))
# # YFDP_limited <- Limit.Lists(YFDP_limited, Shared_spec)
# # YFDP_df <- Eff.Data.Frame(List_ls = YFDP_limited)
# # Region_limited <- Limit.Lists(Region_limited, Shared_spec)
# # Region_df <- Eff.Data.Frame(List_ls = Region_limited)
# # HMSC_col <- grep(colnames(YFDP_df), pattern = "HMSC")
# # HMSC_col <- -HMSC_col[HMSC_col %nin% grep(colnames(YFDP_df), pattern = "diametre")]
# # YFDP_df <- YFDP_df[,HMSC_col]
# # HMSC_col <- grep(colnames(Region_df), pattern = "HMSC")
# # HMSC_col <- -HMSC_col[HMSC_col %nin% grep(colnames(Region_df), pattern = "biomass")]
# # Region_df <- Region_df[,HMSC_col]
# # data <- cbind(Region_df, YFDP_df)
# 
# ## select only pre-fire data for YFDP
# data_ls[[1]] <- data_ls[[1]][,-grep(colnames(data_ls[[1]]), pattern = "Post-Fire")]
# ## select relevant HMSC models
# HMSC_col <- grep(colnames(data_ls[[1]]), pattern = "HMSC")
# data_ls[[1]] <- data_ls[[1]][, -HMSC_col[HMSC_col %nin% grep(colnames(data_ls[[1]]), pattern = "diametre")]]
# HMSC_col <- grep(colnames(data_ls[[2]]), pattern = "HMSC")
# data_ls[[2]] <- data_ls[[2]][, -HMSC_col[HMSC_col %nin% grep(colnames(data_ls[[2]]), pattern = "biomass")]]
# HMSC_col <- grep(colnames(data_ls[[3]]), pattern = "HMSC")
# data_ls[[3]] <- data_ls[[3]][, -HMSC_col[HMSC_col %nin% grep(colnames(data_ls[[3]]), pattern = "biomass")]]
# 
# ## change names of list items and their data frame columns
# data_ls <- lapply(1:length(data_ls), FUN = function(x){
#   x <- data_ls[[x]]
#   colnames(x) <- gsub(x = colnames(x), pattern = "COCCUR", replacement = "COOCCUR")
#   colnames(x) <- gsub(x = colnames(x), pattern = "IF_REM", replacement = "NDD-RIM")
#   colnames(x) <- gsub(x = colnames(x), pattern = "Yosemite", replacement = "Region")
#   colnames(x) <- gsub(x = colnames(x), pattern = "Pre-Fire", replacement = "Plot")
#   colnames(x) <- gsub(x = colnames(x), pattern = "Temperate Conifer Forests", replacement = "Macro")
#   x
# })
# names(data_ls) <- c("Plot", "Region", "Macro")
# data_ls$Region$COOCCUR.Region.Sig <- FALSE

# Species_ls <- base::strsplit(x = as.character(unique(data$Partner1)), split = " ")
# SpeciesLegend <- data.frame(Abbr = paste0(substr(unlist(lapply(Species_ls, "[[", 1)), start = 1, stop = 2),
#                                           substr(unlist(lapply(Species_ls, "[[", 2)), start = 1, stop = 2)),
#                             Full = unique(data$Partner1)
# )
# print(SpeciesLegend)
# data$Partner1 <- SpeciesLegend$Abbr[match(data$Partner1, SpeciesLegend$Full)]
# data$Partner2 <- SpeciesLegend$Abbr[match(data$Partner2, SpeciesLegend$Full)]
## PLOT PRODUCTION ---------------------------------------------------------
# gplot <- Plot.NetMat.Cross(data = data_ls,
#                            ALLRegionOrder = c("Plot", "Region", "Macro"),
#                            ALLMethodOrder = c("COOCCUR", "NETASSOC", "HMSC", "NDD-RIM"))
# ggsave(gplot, file = file.path(Dir.Exports, "Figure2.jpg"), height = 54, width = 32, units = "cm") 


grid2 <-  plot_grid(plotlist = 
                      list(
                        plot_grid(plotlist = list(MATplots_ls$Plot$COCCUR$plot,
                                                  MATplots_ls$Region$COCCUR$plot,
                                                  MATplots_ls$Macro$COCCUR$plot), nrow = 1, 
                                  labels = c("Plot", "Region", "Macro")),
                        plot_grid(plotlist = list(MATplots_ls$Plot$NETASSOC$plot,
                                                  MATplots_ls$Region$NETASSOC$plot,
                                                  MATplots_ls$Macro$NETASSOC$plot), nrow = 1, 
                                  labels = c("Plot", "Region", "Macro")),
                        plot_grid(plotlist = list(MATplots_ls$Plot$HMSC$plot,
                                                  MATplots_ls$Region$HMSC$plot,
                                                  MATplots_ls$Macro$HMSC$plot), nrow = 1, 
                                  labels = c("Plot", "Region", "Macro")),
                        plot_grid(plotlist = list(MATplots_ls$Plot$NDD_RIM_ASSOC$plot,
                                                  MATplots_ls$Region$NDD_RIM_ASSOC$plot,
                                                  MATplots_ls$Macro$NDD_RIM_ASSOC$plot), nrow = 1, 
                                  labels = c("Plot", "Region", "Macro"))
                      ),
                    ncol = 1,
                    labels = c("")
)
ggsave(grid2, filename = file.path(Dir.Exports, "Figure2.jpg"), 
       height = 60, width = 50, units = "cm")

# TABLES – BETADIVERSITY COMPARISON =========================================
message("############ BETADIVERSITY TABLES")
makeigraph <- function(data1){
  if(0 %in% data1$Sig){
    data1$Sig <- as.character(data1$Sig)
    data1$Sig[which(data1$Sig == 0 | data1$Sig == 0.5)] <- FALSE
    data1$Sig[which(data1$Sig == 1)] <- TRUE
  } # then its NDD-RIM data
  data1$Sig[is.na(data1$Sig)] <- FALSE
  ig1 <- graph_from_data_frame(data1[data1$Sig == TRUE,], directed = FALSE)
  if(length(ig1) != 0){
    E(ig1)$weight <-E(ig1)$Value
    E(ig1)$weight[E(ig1)$weight > 0] <- 1
    E(ig1)$weight[E(ig1)$weight < 0] <- -1
    
    origvert <- names(V(ig1))
    addvert <- unique(data1$Partner1)[unique(data1$Partner1) %nin% names(V(ig1))
    ]
    ig1 <- add_vertices(ig1, 
                        nv = length(addvert))
    ig1 <- set.vertex.attribute(ig1, "name", value = c(origvert, addvert))
  }else{
    ig1 <- make_empty_graph(n = length(unique(data1$Partner1)))
    ig1 <- set.vertex.attribute(ig1, "name", value = unique(data1$Partner1))
  }
  ig1
}

betadiv_calc <- function(Compare = "HMSC"){
  Full_ls <- lapply(MATplots_ls, FUN = function(x){
    x[names(x) != "NDD_RIM_ASSOC"]
  })
  
  firstlevl <- names(Full_ls)
  secondlevl <- names(Full_ls[[1]])
  
  if(Compare %in% firstlevl){
    Compare_ls <- lapply(Full_ls[[Compare]], FUN = function(x){
      makeigraph(x$data)
    })
  }
  
  if(Compare %in% secondlevl){
    Compare_ls <- lapply(Full_ls, FUN = function(x){
      makeigraph(x[[Compare]]$data)
    })
  }
  
  betadiv_df <- betalink::network_betadiversity(N = Compare_ls)
  # ratio of dissimilarity of interactions due to species turnover compared to dissimilarity of interactions in shared species-pairs
  betadiv_df$STvOS <- betadiv_df$ST/betadiv_df$OS
  betadiv_df
}

sink(file = file.path(Dir.Exports, "Betadiversity.txt"))
print("############### COOCCUR")
print(betadiv_calc("COCCUR"))
print("############### NETASSOC")
print(betadiv_calc("NETASSOC"))
print("############### HMSC")
print(betadiv_calc("HMSC"))
print("############### NDD-RIM")
print(betadiv_calc("NDD_RIM"))

print("############### Plot")
print(betadiv_calc("Plot"))
print("############### Region")
print(betadiv_calc("Region"))
print("############### Macro")
print(betadiv_calc("Macro"))
sink()

# S5 - NETWORK VISUALISATION ACROSS SCALES ==================================
message("############ S5")
FName <- "S5"

png(filename = file.path(Dir.Exports, paste0(FName, ".png")), width = 24, height = 40, units = "cm", res = 1000)
Plot.Network(Plot_ls = MATplots_ls)
dev.off()


# FIGURE 3 – MODULARITY OF NETWORKS =========================================
message("############ FIGURE 3")
Comparable_ls <- lapply(MATplots_ls, FUN = function(x){
  lapply(x, FUN = function(y){
    y$data
  })
})
Comparable_ls <- lapply(Comparable_ls, FUN = function(x){
  x[-4]
})

Coni_topo <- Calc.Topology(data_ls = Comparable_ls, 
                           Shared_abbr = Shared_abbr,
                           Sig = TRUE)
ggsave(Coni_topo$plots$Modularity, file = file.path(Dir.Exports, "Figure3.jpg"), height = 16, width = 20, units = "cm") 

# FIGURE 4 – CENTRALITY OF NODES CONTAINED IN NETWORKS INFERRED =============
message("############ FIGURE 4")
ggsave(Coni_topo$plots$Strength, file = file.path(Dir.Exports, "Figure4.jpg"), height = 16, width = 24, units = "cm") 

# # S6 - UNDIRECTED NDD-RIM NETWORKS ==========================================
# message("############ S6")
# NDDRIM_undirected <- Coni_topo$numbers$adj[grep(names(Coni_topo$numbers$adj), pattern = "NDD-RIM")]
# NDDRIM_undirected <- lapply(NDDRIM_undirected, FUN = function(x){
#   x[upper.tri(x)] <- NA
#   x
#   })
# plot_df <- data.frame(
#   Partner1 = data$Partner1,
#   Partner2 = data$Partner2,
#   `NDD-RIM.Plot.effects` = as.vector(NDDRIM_undirected$`NDD-RIM.Plot.effects`),
#   `NDD-RIM.Plot.Sig` = TRUE,
#   `NDD-RIM.Region.effects` = as.vector(NDDRIM_undirected$`NDD-RIM.Region.effects`),
#   `NDD-RIM.Region.Sig` = TRUE,
#   `NDD-RIM.Macro.effects` = as.vector(NDDRIM_undirected$`NDD-RIM.Macro.effects`),
#   `NDD-RIM.Macro.Sig` = TRUE
# )
# gplot <- Plot.NetMat.Cross(data = plot_df,
#                            ALLRegionOrder = c("Plot", "Region", "Macro"),
#                            ALLMethodOrder = c("NDD.RIM")
#                     )
# ggsave(gplot, file = file.path(Dir.Exports, "S4a.jpg"), height = 18, width = 32, units = "cm")
# 
# colnames(plot_df) <- gsub(colnames(plot_df), pattern = "NDD.RIM", replacement = "NDD-RIM")
# 
# png(filename = file.path(Dir.Exports, "S4b.png"), width = 32, height = 14, units = "cm", res = 1000)
# netplot <- Plot.Network(Plot_df = plot_df, ModelOrder = "NDD-RIM",
#                         TreatmentOrder = c("Plot", "Region", "Macro"),
#                         FName = "S7", Directed = FALSE)
# dev.off()
