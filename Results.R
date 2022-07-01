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

# WITHIN-YFDP, WITHIN METHOD ===============================================
message("############ WITHIN-YFDP; WITHIN-METHOD")
## DATA --------------------------------------------------------------------
YFDP_limited <- Limit.Lists(YFDP_ls, YFDP_spec)
YFDP_df <- Eff.Data.Frame(List_ls = YFDP_limited)

## MATRICES VISUALISATIONS -------------------------------------------------
gplot <- Plot.NetMat(data = YFDP_df, method = "HMSC") 
ggsave(gplot, file = file.path(DirEx.YFDP, "HMSC-Matrix.png"), height = 32, width = 32, units = "cm") 
gplot <- Plot.NetMat(data = YFDP_df, method = "IF_REM") 
ggsave(gplot, file = file.path(DirEx.YFDP, "IF_REM-Matrix.png"), height = 16, width = 32, units = "cm") 
gplot <- Plot.NetMat(data = YFDP_df, method = "NETASSOC") 
ggsave(gplot, file = file.path(DirEx.YFDP, "NETASSOC-Matrix.png"), height = 16, width = 32, units = "cm") 
gplot <- Plot.NetMat(data = YFDP_df, method = "COCCUR") 
ggsave(gplot, file = file.path(DirEx.YFDP, "COCCUR-Matrix.png"), height = 16, width = 32, units = "cm") 
discard <- grep(colnames(YFDP_df), pattern = "HMSC")[grep(colnames(YFDP_df), pattern = "HMSC") %nin% grep(colnames(YFDP_df), pattern = "diametre")]
gplot <- Plot.NetMat.Cross(data = YFDP_df[,-discard],
                           ALLRegionOrder = c("Pre-Fire", "Post-Fire", "Difference"),
                           ALLMethodOrder = c("COCCUR", "NETASSOC", "HMSC", "IF_REM"))
ggsave(gplot, file = file.path(DirEx.YFDP, "Cross-YFDP.png"), height = 54, width = 32, units = "cm")


## PUBLICATION PLOT --------------------------------------------------------

### COCCUR ----
COCC_dat <- Plot.NetMat(data = YFDP_df, method = "COCCUR", datareturn = TRUE) 
COCC_dat <- COCC_dat[COCC_dat$Condition == "Pre-Fire", ]

COCC_plot <- ggplot(COCC_dat, aes(x = Partner2, y = Partner1, fill = Value)) +
  geom_tile(color = "white",
            lwd = 1.5,
            linetype = 1) + 
  geom_point(aes(shape = Sig), size = 2) +
  scale_shape_manual(values=c(32, 15), na.translate = FALSE, name = "Significance") +  
  guides(shape = FALSE) + 
  coord_fixed() + 
  guides(fill = guide_colourbar(barwidth = 2,
                                barheight = 15,
                                title = "Association")) + 
  ## axes
  theme_bw() + 
  theme(axis.text.x=element_text(angle = -40, hjust = 0)) + 
  scale_fill_gradient2(low = "#ba00bd", high = "#d9db3b") + 
  labs(x = "", y = "")

Empty_plot <- ggplot(COCC_dat, aes(x = Partner2, y = Partner1, fill = Value)) + theme_void()

### NETASSOC ----
NETA_dat <- Plot.NetMat(data = YFDP_df, method = "NETASSOC", datareturn = TRUE) 
NETA_dat <- NETA_dat[NETA_dat$Condition == "Pre-Fire", ]

NETA_plot <- ggplot(NETA_dat, aes(x = Partner2, y = Partner1, fill = Value)) +
  geom_tile(color = "white",
            lwd = 1.5,
            linetype = 1) + 
  geom_point(aes(shape = Sig), size = 2) +
  scale_shape_manual(values=c(32, 15), na.translate = FALSE, name = "Significance") +  
  guides(shape = FALSE) + 
  coord_fixed() + 
  guides(fill = guide_colourbar(barwidth = 2,
                                barheight = 15,
                                title = "Association")) + 
  ## axes
  theme_bw() + 
  theme(axis.text.x=element_text(angle = -40, hjust = 0)) + 
  scale_fill_gradient2(low = "#ba00bd", high = "#d9db3b") + 
  labs(x = "", y = "")

### HMSC ----
HMSC_dat <- Plot.NetMat(data = YFDP_df, method = "HMSC", datareturn = TRUE) 
HMSC_dat <- HMSC_dat[HMSC_dat$Condition == "Pre-Fire", ]

HMSC_dat$Model <- gsub(x = HMSC_dat$Model, pattern = "presence_absence", replacement = "Presence & Absence")
HMSC_dat$Model <- gsub(x = HMSC_dat$Model, pattern = "abundance", replacement = "Abundance")
HMSC_dat$Model <- gsub(x = HMSC_dat$Model, pattern = "diametre", replacement = "Performance")

HMSC_plot1 <- ggplot(HMSC_dat[HMSC_dat$Model == "Presence & Absence",], aes(x = Partner2, y = Partner1, fill = Value)) +
  geom_tile(color = "white",
            lwd = 1.5,
            linetype = 1) + 
  geom_point(aes(shape = Sig), size = 2) +
  scale_shape_manual(values=c(32, 15), na.translate = FALSE, name = "Significance") +  
  guides(shape = FALSE) + 
  coord_fixed() + 
  guides(fill = guide_colourbar(barwidth = 2,
                                barheight = 15,
                                title = "Association")) + 
  ## axes
  theme_bw() + 
  theme(axis.text.x=element_text(angle = -40, hjust = 0)) + 
  scale_fill_gradient2(low = "#ba00bd", high = "#d9db3b") + 
  labs(x = "", y = "")

HMSC_plot2 <- ggplot(HMSC_dat[HMSC_dat$Model == "Abundance",], aes(x = Partner2, y = Partner1, fill = Value)) +
  geom_tile(color = "white",
            lwd = 1.5,
            linetype = 1) + 
  geom_point(aes(shape = Sig), size = 2) +
  scale_shape_manual(values=c(32, 15), na.translate = FALSE, name = "Significance") +  
  guides(shape = FALSE) + 
  coord_fixed() + 
  guides(fill = guide_colourbar(barwidth = 2,
                                barheight = 15,
                                title = "Association")) + 
  ## axes
  theme_bw() + 
  theme(axis.text.x=element_text(angle = -40, hjust = 0)) + 
  scale_fill_gradient2(low = "#ba00bd", high = "#d9db3b") + 
  labs(x = "", y = "")

HMSC_plot3 <- ggplot(HMSC_dat[HMSC_dat$Model == "Performance",], aes(x = Partner2, y = Partner1, fill = Value)) +
  geom_tile(color = "white",
            lwd = 1.5,
            linetype = 1) + 
  geom_point(aes(shape = Sig), size = 2) +
  scale_shape_manual(values=c(32, 15), na.translate = FALSE, name = "Significance") +  
  guides(shape = FALSE) + 
  coord_fixed() + 
  guides(fill = guide_colourbar(barwidth = 2,
                                barheight = 15,
                                title = "Association")) + 
  ## axes
  theme_bw() + 
  theme(axis.text.x=element_text(angle = -40, hjust = 0)) + 
  scale_fill_gradient2(low = "#ba00bd", high = "#d9db3b") + 
  labs(x = "", y = "")

### IF_REM ----
IFREM_dat <- Plot.NetMat(data = YFDP_df, method = "IF_REM", datareturn = TRUE) 
IFREM_dat <- IFREM_dat[IFREM_dat$Condition == "Pre-Fire", ]

IFREM_plot <- ggplot(IFREM_dat, aes(x = Partner2, y = Partner1, fill = Value)) +
  geom_tile(color = "white",
            lwd = 1.5,
            linetype = 1) + 
  geom_point(aes(shape = Sig), size = 2) +
  scale_shape_manual(values=c(32, 15), na.translate = FALSE, name = "Significance") +  
  guides(shape = FALSE) + 
  coord_fixed() + 
  guides(fill = guide_colourbar(barwidth = 2,
                                barheight = 15,
                                title = "Association")) + 
  ## axes
  theme_bw() + 
  theme(axis.text.x=element_text(angle = -40, hjust = 0)) + 
  scale_fill_gradient2(low = "#ba00bd", high = "#d9db3b") + 
  labs(x = "Subject", y = "Actor")

### FINAL OUTPUT ----
plot_grid(plotlist = 
list(
plot_grid(plotlist = list(COCC_plot, Empty_plot, Empty_plot), nrow = 1),
plot_grid(plotlist = list(Empty_plot, NETA_plot, Empty_plot), nrow = 1),
plot_grid(plotlist = list(HMSC_plot1, HMSC_plot2, HMSC_plot3), nrow = 1),
plot_grid(plotlist = list(Empty_plot, Empty_plot, IFREM_plot), nrow = 1)
),
ncol = 1,
labels = c("")
)

plot_grid(plotlist = 
            list(
              plot_grid(plotlist = list(COCC_plot, NETA_plot, IFREM_plot), nrow = 1, 
                        labels = c("COOCUR", "NETASSOC", "IF-REM")),
              plot_grid(plotlist = list(HMSC_plot1, HMSC_plot2, HMSC_plot3), nrow = 1, 
                        labels = c("HMSC (Presence & Absence)", "HMSC (Abundance)", "HMSC (Performance)"))
              ),
          ncol = 1,
          labels = c("")
)
ggsave(filename = file.path(Dir.Exports, "ModelVis.png"), height = 32, width = 54, units = "cm")


## YFDP: Pre-Fire & FIA: Yosemite, Temperate Conifer -----------------------
### DATA --------------------------------------------------------------------
YFDP_limited <- YFDP_ls[grep(names(YFDP_ls), pattern = "Pre-Fire")]
YFDP_spec <- unique(unlist(lapply(YFDP_limited, FUN = function(x){x[,1]})))
Region_limited <- Region_ls[c(grep(names(Region_ls), pattern = "Yosemite"), grep(names(Region_ls), pattern = "Temperate Conifer Forests"))]
Region_spec <- unique(unlist(lapply(Region_limited, FUN = function(x){x[,1]})))
Shared_spec <- Reduce(intersect, list(YFDP_spec, Region_spec))
YFDP_limited <- Limit.Lists(YFDP_limited, Shared_spec)
YFDP_df <- Eff.Data.Frame(List_ls = YFDP_limited)
Region_limited <- Limit.Lists(Region_limited, Shared_spec)
Region_df <- Eff.Data.Frame(List_ls = Region_limited)
HMSC_col <- grep(colnames(YFDP_df), pattern = "HMSC")
HMSC_col <- -HMSC_col[HMSC_col %nin% grep(colnames(YFDP_df), pattern = "diametre")]
YFDP_df <- YFDP_df[,HMSC_col]
HMSC_col <- grep(colnames(Region_df), pattern = "HMSC")
HMSC_col <- -HMSC_col[HMSC_col %nin% grep(colnames(Region_df), pattern = "biomass")]
Region_df <- Region_df[,HMSC_col]
data <- cbind(Region_df, YFDP_df)
colnames(data) <- gsub(x = colnames(data), pattern = "COCCUR", replacement = "COOCCUR")
colnames(data) <- gsub(x = colnames(data), pattern = "IF_REM", replacement = "IF-REM")
colnames(data) <- gsub(x = colnames(data), pattern = "Yosemite", replacement = "Region")
colnames(data) <- gsub(x = colnames(data), pattern = "Pre-Fire", replacement = "Plot")
colnames(data) <- gsub(x = colnames(data), pattern = "Temperate Conifer Forests", replacement = "Macro")
data$COOCCUR.Region.Sig <- FALSE
### TOPOLOGY ----------------------------------------------------------------
Coni_topo <- Calc.Topology(data = data, Sig = TRUE, Model = "ALL", TreatmentOrder = c("Plot", "Region", "Macro"))
ggsave(Coni_topo$plots$Strength, file = file.path(Dir.Exports, "Coni-Topo_Strength.png"), height = 32, width = 32, units = "cm") 
ggsave(Coni_topo$plots$Eigenvector, file = file.path(Dir.Exports, "Coni-Topo_Eigenvector.png"), height = 32, width = 32, units = "cm") 
ggsave(Coni_topo$plots$Modularity, file = file.path(Dir.Exports, "Coni-Topo_Modularity.png"), height = 24, width = 32, units = "cm") 
ggsave(Coni_topo$plots$Nestedness, file = file.path(Dir.Exports, "Coni-Topo_Nestedness.png"), height = 24, width = 32, units = "cm")  
### MATRICES VISUALISATIONS -------------------------------------------------
gplot <- Plot.NetMat.Cross(data = data,
                           ALLRegionOrder = c("Plot", "Region", "Macro"),
                           ALLMethodOrder = c("COOCCUR", "NETASSOC", "HMSC", "IF-REM"))
ggsave(gplot, file = file.path(Dir.Exports, "Cross-Coni.png"), height = 54, width = 32, units = "cm") 
# ### NETWORKS VISUALISED -----------------------------------------------------
FName <- "Cross-Coni2"
png(filename = file.path(Dir.Exports, paste0(FName, ".png")), width = 24, height = 32, units = "cm", res = 1000)
Plot.Network(Plot_df = data,
             ModelOrder = c("COOCCUR", "NETASSOC", "HMSC", "IF-REM"),
             TreatmentOrder = c("Plot", "Region", "Macro"),
             FName = FName,
             Dir = Dir.Exports)
dev.off()
