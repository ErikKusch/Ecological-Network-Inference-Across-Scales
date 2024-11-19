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

YFDP_full <- Eff.Data.Frame(List_ls = Limit.Lists(YFDP_ls, YFDP_spec))
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

# NETWORK INFERENCE IS NOT CONSISTENT ACROSS APPROACHES =========
message("############ NETWORK INFERENCE IS NOT CONSISTENT ACROSS APPROACHES")
## DATA --------------------------------------------------------------------
Fig1S3S4_ls <- list(Plot = YFDP_full,
                  Local = Region_full,
                  Continental = Macro_full)
MATplots_ls <- as.list(rep(NA, 3))
names(MATplots_ls) <- c("Plot", "Region", "Macro")
MATplots_ls <- list(
  SignOnly = MATplots_ls,
  RawData = MATplots_ls
)
MATs_ls <- MATplots_ls

## PLOT PRODUCTION ---------------------------------------------------------
for(k in 1:length(MATplots_ls)){
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
    
    factorcols <- c("#ba00bd", "white", "#d9db3b", "darkgrey")
    names(factorcols) <- c("-1", "0", "1", "NA")
    
    ### COCCUR ----
    COCC_dat <- Plot.NetMat(data = plot_df, method = "COCCUR", datareturn = TRUE) 
    COCC_dat <- COCC_dat[COCC_dat$Condition == Treatment , ]
    COCC_dat <- bordercol(COCC_dat)
    COCC_dat$Sig[is.na(COCC_dat$Sig)] <- FALSE
    COCC_dat$Value[is.na(COCC_dat$Value) & !COCC_dat$Sig] <- 0
    COCC_dat$Value[!COCC_dat$Sig] <- 0
    COCC_dat$Value[is.na(COCC_dat$Value)] <- 0

    if(k == 1){
      COCC_dat$Value <- factor(sign(COCC_dat$Value))
      COCC_plot <- ggplot(COCC_dat, aes(x = Partner2, y = Partner1, fill = Value, col = Shared)) +
        geom_tile(lwd = 1,
                  linetype = 1,
                  width=0.92, height=0.92) + 
        geom_point(aes(shape = Sig), col = "black", size = 2) +
        scale_shape_manual(values=c(32, 15), na.translate = FALSE, name = "Significance") +  
        scale_fill_manual(values = 
                            factorcols[names(factorcols) %in% as.character(COCC_dat$Value)],
                          name = "Association") + 
        guides(shape = FALSE) + 
        coord_fixed() + 
        scale_color_manual(values=c("white", "black"), na.translate = TRUE, name = "", guide = "none") +  
        ## axes
        theme_bw() + 
        theme(axis.text.x=element_text(angle = -40, hjust = 0),
              text = element_text(size=15),
              legend.position="bottom") +
        labs(x = "", y = "") 
    }else{
      COCC_plot <- ggplot(COCC_dat, aes(x = Partner2, y = Partner1, fill = Value, col = Shared)) +
        geom_tile(lwd = 1,
                  linetype = 1,
                  width=0.92, height=0.92) + 
        geom_point(aes(shape = Sig), col = "black", size = 2) +
        scale_shape_manual(values=c(32, 15), na.translate = FALSE, name = "Significance") +  
        guides(shape = FALSE) + 
        coord_fixed() + 
        scale_color_manual(values=c("white", "black"), na.translate = TRUE, name = "", guide = "none") +  
        guides(fill = guide_colourbar(barwidth = 2,
                                      barheight = 15,
                                      title = "Association")) + 
        ## axes
        theme_bw() + 
        theme(axis.text.x=element_text(angle = -40, hjust = 0),
              text = element_text(size=15)) +
        scale_fill_gradient2(low = "#ba00bd", high = "#d9db3b") + 
        labs(x = "", y = "")
    }
    
    
    Empty_plot <- ggplot(COCC_dat, aes(x = Partner2, y = Partner1, fill = Value)) + theme_void()
    
    ### NETASSOC ----
    NETA_dat <- Plot.NetMat(data = plot_df, method = "NETASSOC", datareturn = TRUE) 
    NETA_dat <- NETA_dat[NETA_dat$Condition == Treatment, ]
    NETA_dat <- bordercol(NETA_dat)
    NETA_dat$Value[is.na(NETA_dat$Value) & !NETA_dat$Sig] <- 0
    NETA_dat$Value[!NETA_dat$Sig] <- 0
    NETA_dat$Value[is.na(NETA_dat$Value)] <- 0
    
    if(k == 1){
      NETA_dat$Value <- factor(sign(NETA_dat$Value))
      NETA_plot <- ggplot(NETA_dat, aes(x = Partner2, y = Partner1, fill = Value, col = Shared)) +
        geom_tile(lwd = 1,
                  linetype = 1,
                  width=0.92, height=0.92) + 
        geom_point(aes(shape = Sig), col = "black", size = 2) +
        scale_shape_manual(values=c(32, 15), na.translate = FALSE, name = "Significance") +  
        scale_color_manual(values=c("white", "black"), na.translate = TRUE, name = "", guide = "none") +  
        guides(shape = FALSE) + 
        coord_fixed() + 
        ## axes
        theme_bw() + 
        scale_fill_manual(values = 
                            factorcols[names(factorcols) %in% as.character(NETA_dat$Value)], 
                          name = "Association") +
        theme(axis.text.x=element_text(angle = -40, hjust = 0),
              text = element_text(size=15),
              legend.position="bottom") +
        labs(x = "", y = "")
  
    }else{
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
    }
    
    
    
    ### HMSC ----
    HMSC_dat <- Plot.NetMat(data = plot_df, method = "HMSC", datareturn = TRUE) 
    HMSC_dat <- HMSC_dat[HMSC_dat$Condition == Treatment, ]
    HMSC_dat <- bordercol(HMSC_dat)
    HMSC_dat$Value[is.na(HMSC_dat$Value) & !HMSC_dat$Sig] <- 0
    HMSC_dat$Value[!HMSC_dat$Sig] <- 0
    
    HMSC_dat$Model <- gsub(x = HMSC_dat$Model, pattern = "presence_absence", replacement = "Presence & Absence")
    HMSC_dat$Model <- gsub(x = HMSC_dat$Model, pattern = "abundance", replacement = "Abundance")
    HMSC_dat$Model <- gsub(x = HMSC_dat$Model, pattern = "diametre", replacement = "Performance")
    HMSC_dat$Model <- gsub(x = HMSC_dat$Model, pattern = "biomass", replacement = "Performance")
    
    if(k == 1){
      HMSC_dat$Value <- factor(sign(HMSC_dat$Value))
      
      HMSC_plot1 <- ggplot(HMSC_dat[HMSC_dat$Model == "Presence & Absence",], aes(x = Partner2, y = Partner1, fill = Value, col = Shared)) +
        geom_tile(lwd = 1,
                  linetype = 1,
                  width=0.92, height=0.92) + 
        geom_point(aes(shape = Sig), col = "black", size = 2) +
        scale_shape_manual(values=c(32, 15), na.translate = FALSE, name = "Significance") +  
        scale_color_manual(values=c("white", "black"), na.translate = TRUE, name = "", guide = "none") +  
        guides(shape = FALSE) + 
        coord_fixed() + 
        ## axes
        theme_bw() + 
        scale_fill_manual(values = 
                            factorcols[names(factorcols) %in% 
                                         as.character(HMSC_dat[HMSC_dat$Model == "Presence & Absence",]$Value)],
                          name = "Association") + 
        theme(axis.text.x=element_text(angle = -40, hjust = 0),
              text = element_text(size=15),
              legend.position="bottom") +
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
        ## axes
        theme_bw() + 
        scale_fill_manual(values = 
                            factorcols[names(factorcols) %in% 
                                         as.character(HMSC_dat[HMSC_dat$Model == "Abundance",]$Value)],
                          name = "Association") + 
        theme(axis.text.x=element_text(angle = -40, hjust = 0),
              text = element_text(size=15),
              legend.position="bottom") +
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
        ## axes
        theme_bw() + 
        scale_fill_manual(values = 
                            factorcols[names(factorcols) %in% 
                                         as.character(HMSC_dat[HMSC_dat$Model == "Performance",]$Value)],
                          name = "Association") + 
        theme(axis.text.x=element_text(angle = -40, hjust = 0),
              text = element_text(size=15),
              legend.position="bottom") +
        labs(x = "", y = "")
      
    }else{
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
    }
    
    
    ### IF_REM ----
    IFREM_dat <- Plot.NetMat(data = plot_df, method = "IF_REM", datareturn = TRUE) 
    IFREM_dat <- IFREM_dat[IFREM_dat$Condition == Treatment, ]
    IFREM_dat <- bordercol(IFREM_dat)
    
    # IFREM_dat$Value[is.na(IFREM_dat$Value) & IFREM_dat$Sig != 1] <- 0
    IFREM_dat$Value[!(IFREM_dat$Sig) & !is.na(IFREM_dat$Value)] <- 0
    IFREM_dat$Sig <- factor(IFREM_dat$Sig)
    # IFREM_ASSOC_dat$Value[is.na(IFREM_ASSOC_dat$Value)] <- 0
    
    if(k == 1){
      IFREM_dat$Value <- factor(sign(IFREM_dat$Value))
      
      IFREM_plot <- ggplot(IFREM_dat, aes(x = Partner2, y = Partner1, fill = Value, col = Shared)) +
        geom_tile(lwd = 1,
                  linetype = 1,
                  width=0.92, height=0.92) + 
        geom_point(aes(shape = Sig), col = "black", size = 2) +
        scale_shape_manual(values=c(32, 15), na.translate = FALSE, name = "Significance") +  
        scale_color_manual(values=c("white", "black"), na.translate = TRUE, name = "", guide = "none") +  
        guides(shape = FALSE) + 
        coord_fixed() + 
        ## axes
        theme_bw() + 
        scale_fill_manual(values = 
                            factorcols[names(factorcols) %in% as.character(IFREM_dat$Value)], 
                          name = "Interaction") +
        theme(axis.text.x=element_text(angle = -40, hjust = 0),
              text = element_text(size=15),
              legend.position="bottom") +
        labs(x = "Subject", y = "Actor")
      
    }else{
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
    }
    
    ### IF_REM UNDIRECTED ----
    IFREM_ASSOC_dat <- Plot.NetMat(data = plot_df, method = "IF_REM_Assoc", datareturn = TRUE) 
    IFREM_ASSOC_dat <- IFREM_ASSOC_dat[IFREM_ASSOC_dat$Condition == Treatment, ]
    IFREM_ASSOC_dat <- bordercol(IFREM_ASSOC_dat)
    IFREM_ASSOC_dat$Sig[is.na(IFREM_ASSOC_dat$Sig)] <- 0
    IFREM_ASSOC_dat$Sig <- factor(IFREM_ASSOC_dat$Sig)
    
    IFREM_ASSOC_dat$Value[is.na(IFREM_ASSOC_dat$Value) & IFREM_ASSOC_dat$Sig != 1] <- 0
    IFREM_ASSOC_dat$Value[IFREM_ASSOC_dat$Sig != 1] <- 0
    IFREM_ASSOC_dat$Value[is.na(IFREM_ASSOC_dat$Value)] <- 0
    
    if(k == 1){
      IFREM_ASSOC_dat$Value <- factor(sign(IFREM_ASSOC_dat$Value))
      
      IFREM_ASSOC_plot <- ggplot(IFREM_ASSOC_dat, aes(x = Partner2, y = Partner1, fill = Value, col = Shared)) +
        geom_tile(lwd = 1,
                  linetype = 1,
                  width=0.92, height=0.92) + 
        geom_point(aes(shape = Sig), col = "black", size = 2) +
        scale_shape_manual(values=c(32, 0, 15), na.translate = FALSE, name = "Significance") +  
        scale_color_manual(values=c("white", "black"), na.translate = TRUE, name = "", guide = "none") +  
        guides(shape = FALSE) + 
        coord_fixed() + 
        ## axes
        theme_bw() + 
        scale_fill_manual(values = 
                            factorcols[names(factorcols) %in% as.character(IFREM_ASSOC_dat$Value)], 
                          name = "Association") +
        theme(axis.text.x=element_text(angle = -40, hjust = 0),
              text = element_text(size=15),
              legend.position="bottom") +
        labs(x = "", y = "")
    
    }else{
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
    }
    
    if(k == 1){
      grid2 <-  plot_grid(plotlist = 
                            list(
                              plot_grid(plotlist = lapply(list(COCC_plot, NETA_plot, IFREM_plot), FUN = function(x){
                                x + theme(legend.position = "none")
                              }), 
                              nrow = 1, 
                                        labels = c("COOCUR", "NETASSOC", "NDD-RIM")),
                              plot_grid(plotlist = lapply(list(HMSC_plot1, HMSC_plot2, HMSC_plot3), FUN = function(x){
                                x + theme(legend.position = "none")
                              }),
                              nrow = 1, 
                                        labels = c("HMSC (Presence & Absence)", "HMSC (Abundance)", "HMSC (Performance)"))
                              ,
                              as_ggplot(get_legend(HMSC_plot1))
                            ),
                          rel_heights = c(1, 1, 0.1), 
                          ncol = 1,
                          labels = c("")
      )
      ggsave(grid2, filename = file.path(Dir.Exports, paste0("Approaches_", names(Fig1S3S4_ls)[i], "_SIGNONLY.jpg")), height = 37, width = 50, units = "cm")
    }else{
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
      ggsave(grid2, filename = file.path(Dir.Exports, paste0("Approaches_", names(Fig1S3S4_ls)[i], ".jpg")), height = 30, width = 50, units = "cm") 
    }
    
    MATplots_ls[[k]][[i]] <- list(COCCUR = list(plot = COCC_plot, 
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
}

# NETWORK INFERENCE IS NOT CONSISTENT ACROSS SCALES ============================
message("############ NETWORK INFERENCE IS NOT CONSISTENT ACROSS SCALES")
grid2 <-  plot_grid(plotlist = 
                      list(
                        plot_grid(plotlist = 
                                    lapply(
                                      list(MATplots_ls$SignOnly$Plot$COCCUR$plot,
                                           MATplots_ls$SignOnly$Region$COCCUR$plot,
                                           MATplots_ls$SignOnly$Macro$COCCUR$plot), 
                                      FUN = function(x){
                                        x + theme(legend.position = "none")
                                      }), 
                                  nrow = 1, 
                                  labels = c("Plot", "Region", "Macro")),
                        plot_grid(plotlist = 
                                    lapply(
                                      list(MATplots_ls$SignOnly$Plot$NETASSOC$plot,
                                           MATplots_ls$SignOnly$Region$NETASSOC$plot,
                                           MATplots_ls$SignOnly$Macro$NETASSOC$plot), 
                                      FUN = function(x){
                                        x + theme(legend.position = "none")
                                      }), 
                                  nrow = 1, 
                                  labels = c("Plot", "Region", "Macro")),
                        plot_grid(plotlist = 
                                    lapply(
                                      list(MATplots_ls$SignOnly$Plot$HMSC$plot,
                                           MATplots_ls$SignOnly$Region$HMSC$plot,
                                           MATplots_ls$SignOnly$Macro$HMSC$plot), 
                                      FUN = function(x){
                                        x + theme(legend.position = "none")
                                      }), 
                                  nrow = 1, 
                                  labels = c("Plot", "Region", "Macro")),
                        plot_grid(plotlist = 
                                    lapply(
                                      list(MATplots_ls$SignOnly$Plot$NDD_RIM_ASSOC$plot,
                                           MATplots_ls$SignOnly$Region$NDD_RIM_ASSOC$plot,
                                           MATplots_ls$SignOnly$Macro$NDD_RIM_ASSOC$plot), 
                                      FUN = function(x){
                                        x + theme(legend.position = "none")
                                      }), 
                                      nrow = 1, 
                                  labels = c("Plot", "Region", "Macro")),
                        as_ggplot(get_legend(MATplots_ls$SignOnly$Plot$HMSC$plot))
                      ),
                    rel_heights = c(1, 1, 1, 1, 0.1), 
                    ncol = 1,
                    labels = c("")
)
ggsave(grid2, filename = file.path(Dir.Exports, "Scales_SIGNONLY.jpg"), 
       height = 60, width = 38, units = "cm")


grid2 <-  plot_grid(plotlist = 
                      list(
                        plot_grid(plotlist = list(MATplots_ls$RawData$Plot$COCCUR$plot,
                                                  MATplots_ls$RawData$Region$COCCUR$plot,
                                                  MATplots_ls$RawData$Macro$COCCUR$plot), nrow = 1, 
                                  labels = c("Plot", "Region", "Macro")),
                        plot_grid(plotlist = list(MATplots_ls$RawData$Plot$NETASSOC$plot,
                                                  MATplots_ls$RawData$Region$NETASSOC$plot,
                                                  MATplots_ls$RawData$Macro$NETASSOC$plot), nrow = 1, 
                                  labels = c("Plot", "Region", "Macro")),
                        plot_grid(plotlist = list(MATplots_ls$RawData$Plot$HMSC$plot,
                                                  MATplots_ls$RawData$Region$HMSC$plot,
                                                  MATplots_ls$RawData$Macro$HMSC$plot), nrow = 1, 
                                  labels = c("Plot", "Region", "Macro")),
                        plot_grid(plotlist = list(MATplots_ls$RawData$Plot$NDD_RIM_ASSOC$plot,
                                                  MATplots_ls$RawData$Region$NDD_RIM_ASSOC$plot,
                                                  MATplots_ls$RawData$Macro$NDD_RIM_ASSOC$plot), nrow = 1, 
                                  labels = c("Plot", "Region", "Macro"))
                      ),
                    ncol = 1,
                    labels = c("")
)
ggsave(grid2, filename = file.path(Dir.Exports, "Scales.jpg"), 
       height = 60, width = 50, units = "cm")

# TABLES – BETADIVERSITY COMPARISON =========================================
message("############ BETADIVERSITY TABLES")
makeigraph <- function(data1){
  # data1$Sig <- as.numeric(data1$Sig)
  if(class(data1$Sig) == "factor"){
    data1$Sig <- as.character(data1$Sig)
    data1$Sig[which(data1$Sig == 0 | data1$Sig == 0.5)] <- FALSE
    data1$Sig[which(data1$Sig == 1)] <- TRUE
  } # then its NDD-RIM data
  data1$Sig[is.na(data1$Sig)] <- FALSE
  ig1 <- graph_from_data_frame(data1[data1$Sig == TRUE,], directed = FALSE)
  if(length(igraph::E(ig1)) != 0){
    igraph::E(ig1)$weight <-igraph::E(ig1)$Value
    igraph::E(ig1)$weight[igraph::E(ig1)$weight > 0] <- 1
    igraph::E(ig1)$weight[igraph::E(ig1)$weight < 0] <- -1
    
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

betadiv_calc <- function(Compare = "HMSC", Approach = "Matrix"){
  Full_ls <- lapply(MATplots_ls$RawData, FUN = function(x){
    lapply(x[names(x) != "NDD_RIM"], "[[", "data")
  })
  
  firstlevl <- names(Full_ls)
  secondlevl <- names(Full_ls[[1]])
  
  if(Compare %in% firstlevl){
    Compare_ls <- lapply(Full_ls[[Compare]], FUN = function(x){
      makeigraph(x)
    })
  }
  
  if(Compare %in% secondlevl){
    Compare_ls <- lapply(Full_ls, FUN = function(x){
      makeigraph(x[[Compare]])
    })
  }
  
  if(Approach == "Matrix"){
    Cnames <- names(Compare_ls)
    specs <- sort(Reduce(intersect, lapply(Compare_ls, FUN = function(y){V(y)$name})))
    Compare_ls <- lapply(Compare_ls, FUN = function(x){
      if(length(igraph::E(x)) == 0){
        z <- as.matrix(as_adjacency_matrix(x))
      }else{
        z <- as.matrix(as_adjacency_matrix(x, attr = "weight"))
      }
      keep <- na.omit(match(
        specs,
        colnames(z))
      )
      # which(colnames(z) %in% c("Abco", "Abma", "Cade", "Conu", "Pila", "Pipo", "Psme", "Quke"))
      z <- z[keep,keep]
      z[lower.tri(z)] <- NA
      diag(z) <- NA
      z
    })
    names(Compare_ls) <- Cnames
    
    FUN_Matcomparison <- function(mat1, mat2){
      eq <- mat1==mat2 # avoid to later compute this twice
      # eq <- ifelse(eq, 0, mat2) # get the desired matrix
      sum(!eq, na.rm = TRUE)/sum(!is.na(eq)) # get the percentage of non equal values
      #[1] 33.33
    }
    
    comp_df <- as.data.frame(t(combn(names(Compare_ls), 2)))
    colnames(comp_df) <- c("i", "j")
    comp_df$OS <- apply(comp_df, MARGIN = 1, FUN = function(z){
      FUN_Matcomparison(Compare_ls[[z[1]]], Compare_ls[[z[2]]])
    })
    betadiv <- comp_df
  }else{
    betadiv_df <- betalink::network_betadiversity(N = Compare_ls)
    # ratio of dissimilarity of interactions due to species turnover compared to dissimilarity of interactions in shared species-pairs
    # betadiv_df$STvOS <- betadiv_df$ST/betadiv_df$OS
    betadiv <- betadiv_df[, c("i", "j", "OS")]
  }
  betadiv
}

for(ComApp_iter in c("Matrix", "Betadiv")){
  sink(file = file.path(Dir.Exports, paste0("Betadiversity_", ComApp_iter,".txt")))
  print("############### COOCCUR")
  print(betadiv_calc("COCCUR", Approach = ComApp_iter))
  print("############### NETASSOC")
  print(betadiv_calc("NETASSOC", Approach = ComApp_iter))
  print("############### HMSC")
  print(betadiv_calc("HMSC", Approach = ComApp_iter))
  print("############### NDD-RIM")
  print(betadiv_calc("NDD_RIM_ASSOC", Approach = ComApp_iter))
  
  print("############### Plot")
  print(betadiv_calc("Plot", Approach = ComApp_iter))
  print("############### Region")
  print(betadiv_calc("Region", Approach = ComApp_iter))
  print("############### Macro")
  print(betadiv_calc("Macro", Approach = ComApp_iter))
  sink()
  
  betadivWithin_gg <- rbind(
    cbind(betadiv_calc("COCCUR", Approach = ComApp_iter), Method = "COCCUR")#[,c(1,2,4,8)]
    ,
    cbind(betadiv_calc("NETASSOC", Approach = ComApp_iter), Method = "NETASSOC")#[,c(1,2,4,8)]
    ,
    cbind(betadiv_calc("HMSC", Approach = ComApp_iter), Method = "HMSC")#[,c(1,2,4,8)]
    ,
    cbind(betadiv_calc("NDD_RIM_ASSOC", Approach = ComApp_iter), Method = "NDD_RIM_ASSOC")#[,c(1,2,4,8)]
  )
  colnames(betadivWithin_gg)[3] <- "Dissimilarity"
  betadivWithin_gg$Dissimilarity <- round(betadivWithin_gg$Dissimilarity*100, 2)
  betadivWithin_gg$Method <- gsub(pattern = "COCCUR", replacement = "COOCCUR", x = betadivWithin_gg$Method)
  betadivWithin_gg$Method <- gsub(pattern = "NDD_RIM_ASSOC", replacement = "NDD-RIM", x = betadivWithin_gg$Method)
  
  betadivWithin_gg <- ggplot(betadivWithin_gg, aes(x = factor(j, level = c("Plot", "Region", "Macro")), 
                                                   y = factor(i, level = c("Plot", "Region", "Macro"))
  )) +
    geom_tile(aes(fill = Dissimilarity)) +
    geom_label(aes(label= paste0(Dissimilarity, "%"))) + 
    facet_wrap(~ factor(Method, level = c("COOCCUR", "NETASSOC", "HMSC", "NDD-RIM"))) + 
    scale_fill_viridis_c(option = "F", direction = -1, begin = 0.3) + 
    theme_bw() + labs(x = "Scale", y = "Scale")
  ggsave(betadivWithin_gg, file = file.path(Dir.Exports, paste0("FigureBetaDivWithin", ComApp_iter,".jpg")), height = 15, width = 21, units = "cm") 
  
  betadivAcross_gg <- rbind(
    cbind(betadiv_calc("Macro", Approach = ComApp_iter), Scale = "Macro")#[,c(1,2,4,8)]
    ,
    cbind(betadiv_calc("Region", Approach = ComApp_iter), Scale = "Region")#[,c(1,2,4,8)]
    ,
    cbind(betadiv_calc("Plot", Approach = ComApp_iter), Scale = "Plot")#[,c(1,2,4,8)]
  )
  colnames(betadivAcross_gg)[3] <- "Dissimilarity"
  betadivAcross_gg$Dissimilarity <- round(betadivAcross_gg$Dissimilarity*100, 2)
  betadivAcross_gg$i <- gsub(pattern = "COCCUR", replacement = "COOCCUR", x = betadivAcross_gg$i)
  betadivAcross_gg$j <- gsub(pattern = "NDD_RIM_ASSOC", replacement = "NDD-RIM", x = betadivAcross_gg$j)
  
  betadivAcross_gg <- ggplot(betadivAcross_gg, aes(x = factor(j, level = c("COOCCUR", "NETASSOC", "HMSC", "NDD-RIM")), 
                                                   y = factor(i, level = c("COOCCUR", "NETASSOC", "HMSC", "NDD-RIM"))
  )) +
    geom_tile(aes(fill = Dissimilarity)) +
    geom_label(aes(label= paste0(Dissimilarity, "%"))) + 
    facet_wrap(~ factor(Scale, level = c("Plot", "Region", "Macro"))) + 
    scale_fill_viridis_c(option = "F", direction = -1, begin = 0.3) + 
    theme_bw() + labs(x = "Approach", y = "Approach")  
  betadivAcross_gg
  ggsave(betadivAcross_gg, file = file.path(Dir.Exports, paste0("FigureBetaDivAcross", ComApp_iter,".jpg")), height = 12, width = 30, units = "cm") 
}

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

# S7&8 - INPUT DATA TYPES AND VARIATION =======================================
load(file.path(Dir.YFDP,"Pre-Fire_ModelFrames.RData"))
LocalFrames_ls <- ModelFrames_ls
load(file.path(Dir.FIA, "FIABiome1003.RData"))
RegionalFrames_ls <- ModelFrames_ls
load(file.path(Dir.FIA, "FIABiome9.RData"))
MacroFrames_ls <- ModelFrames_ls

## Within-Scale Data Correspondence ----
WithinCorres_ls <- list(Local = LocalFrames_ls,
                        Region = RegionalFrames_ls,
                        Macro = MacroFrames_ls)

WithinCorres_ls <- lapply(WithinCorres_ls, FUN = function(x){
  cor_df <- data.frame(Performance = unlist(c(as.matrix(x$FitCom[,-1]))),
                       Abundance =  unlist(c(as.matrix(x$Community[,-1]))),
                       Occurrence = unlist(c(sign(as.matrix(x$Community[,-1])))))
  cor_df <- na.omit(cor_df)
  # cor(cor_df$Performance, cor_df$Abundance, use = "complete.obs")
  
  # cowplot::plot_grid(
    ggplot(data = cor_df, aes(x = Performance, y = Abundance)) +
      geom_point() +
      theme_bw()
  #   ,
  #   ggplot(data = data.frame(Value = c(cor_df$Performance, cor_df$Abundance),
  #                            Input = rep(c("Performance", "Abundance"), each = nrow(cor_df))), 
  #          aes(x = Input, y = Value)) +
  #     geom_boxplot() +
  #     # stat_compare_means(paired = TRUE) + 
  #     theme_bw()
  # )
})
ggsave(plot_grid(plotlist = WithinCorres_ls, ncol = 1, labels = "AUTO"),
       file = file.path(Dir.Exports, "S7.jpg"), height = 18*2, width = 16*2, units = "cm"
)

## Across-Scale Data Correspondence ----
AcrossCorres_ls <- list(Local = LocalFrames_ls,
                        Region = RegionalFrames_ls,
                        Macro = MacroFrames_ls)


AcrossCorres_ls <- lapply(AcrossCorres_ls, FUN = function(x){
  rbind(
    data.frame(x$FitCom[,Shared_spec], Input = "Performance"),
    data.frame(x$Community[,Shared_spec], Input = "Abundance")
  )
})

AcrossCorres_ls <- lapply(gsub(Shared_spec, pattern = " ", replacement = "."), FUN = function(x){
  plot_df <- do.call(rbind, lapply(names(AcrossCorres_ls), FUN = function(y){
    ret_df <- AcrossCorres_ls[[y]][, c(x, "Input")]
    ret_df$Scale <- y
    colnames(ret_df)[1] <- "Value"
    ret_df
  }))
  ggplot(plot_df, aes(y = Value, x = factor(Scale, levels = c("Local", "Region", "Macro")))) + 
    geom_boxplot() + 
    facet_wrap(~Input) + 
    # stat_compare_means() +
    theme_bw() + labs(title = gsub(x, pattern = "\\.", replacement = " "), x = "Scale")
})
ggsave(plot_grid(plotlist = AcrossCorres_ls, ncol = 1, labels = "AUTO"),
       file = file.path(Dir.Exports, "S8.jpg"), height = 32*2, width = 16*2, units = "cm")
