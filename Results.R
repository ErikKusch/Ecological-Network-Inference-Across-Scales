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



## TOPOLOGY ----------------------------------------------------------------
YFDP_df2 <- YFDP_df[,-grep(pattern = "NETASSOC", colnames(YFDP_df))]
YFDP_topo <- Calc.Topology(data = YFDP_df2, Sig = TRUE, Model = "ALL", TreatmentOrder = c("Pre-Fire", "Post-Fire"))
ggsave(YFDP_topo$plots$Strength, file = file.path(DirEx.YFDP, "Topo_Strength.png"), height = 32, width = 32, units = "cm") 
ggsave(YFDP_topo$plots$Eigenvector, file = file.path(DirEx.YFDP, "Topo_Eigenvector.png"), height = 32, width = 32, units = "cm") 
ggsave(YFDP_topo$plots$Modularity, file = file.path(DirEx.YFDP, "Topo_Modularity.png"), height = 24, width = 32, units = "cm") 
ggsave(YFDP_topo$plots$Nestedness, file = file.path(DirEx.YFDP, "Topo_Nestedness.png"), height = 24, width = 32, units = "cm") 
# YFDP_topo$numbers
## NETWORKS VISUALISED -----------------------------------------------------
### all methods ----
FName <- "Cross-YFDP2"
png(filename = file.path(DirEx.YFDP, paste0(FName, ".png")), width = 16, height = 32, units = "cm", res = 1000)
discard <- grep(colnames(YFDP_df), pattern = "HMSC")[grep(colnames(YFDP_df), pattern = "HMSC") %nin% grep(colnames(YFDP_df), pattern = "diametre")]
Plot.Network(Plot_df = YFDP_df[,-discard],  
             ModelOrder = c("COCCUR", "NETASSOC", "HMSC", "IF_REM"),
             TreatmentOrder = c("Pre-Fire", "Post-Fire"), 
             FName = FName,
             Dir = DirEx.YFDP)
dev.off()
### HMSC only ----
FName <- "HMSC-YFDP"
png(filename = file.path(DirEx.YFDP, paste0(FName, ".png")), width = 16, height = 24, units = "cm", res = 1000)
keep <- grep(colnames(YFDP_df), pattern = "HMSC")
Part_df <- YFDP_df[,1:2]
YFDP_df <- YFDP_df[,keep]
names <- strsplit(colnames(YFDP_df), split = "[.]")
colnames(YFDP_df) <- paste(lapply(names, "[[", 3), 
                           lapply(names, "[[", 2), 
                           lapply(names, "[[", 4),
                           sep = ".")
YFDP_df <- cbind(Part_df, YFDP_df)
Plot.Network(Plot_df = YFDP_df,  
             ModelOrder = c("presence_absence", "abundance", "diametre"),
             TreatmentOrder = c("Pre-Fire", "Post-Fire"), 
             FName = FName,
             Dir = DirEx.YFDP)
dev.off()

# WITHIN-FIA, WITHIN METHOD ================================================
message("############ WITHIN-FIA; WITHIN-METHOD")
## FIA1: Vermont, Maine, Temperate Conifer ---------------------------------
### DATA --------------------------------------------------------------------
Region_ls2 <- Region_ls
Region_ls2 <- Region_ls2[c(grep(names(Region_ls), pattern = "Vermont"), 
                           grep(names(Region_ls), pattern = "Maine"),
                           grep(names(Region_ls), pattern = "Temperate Broadleaf & Mixed Forests")
)]
Shared_spec <- Reduce(intersect, list(
  Region_ls2$COCCUR.Vermont$Partner1,
  Region_ls2$COCCUR.Maine$Partner1,
  Region_ls2$`COCCUR.Temperate Broadleaf & Mixed Forests`$Partner1
))
Region_limited <- Limit.Lists(Region_ls2, Shared_spec)
Region_df <- Eff.Data.Frame(List_ls = Region_limited)
### TOPOLOGY ----------------------------------------------------------------
Region_df2 <- Region_df[,-grep(pattern = "COCCUR", colnames(Region_df))]
Region_topo <- Calc.Topology(data = Region_df2, Sig = TRUE, Model = "ALL", TreatmentOrder = c("Vermont", "Maine", "Temperate Broadleaf & Mixed Forests"))
ggsave(Region_topo$plots$Strength, file = file.path(DirEx.Region, "Broad-Topo_Strength.png"), height = 32, width = 32, units = "cm") 
ggsave(Region_topo$plots$Eigenvector, file = file.path(DirEx.Region, "Broad-Topo_Eigenvector.png"), height = 32, width = 32, units = "cm") 
ggsave(Region_topo$plots$Modularity, file = file.path(DirEx.Region, "Broad-Topo_Modularity.png"), height = 16, width = 32, units = "cm") 
ggsave(Region_topo$plots$Nestedness, file = file.path(DirEx.Region, "Broad-Topo_Nestedness.png"), height = 16, width = 32, units = "cm")  
### MATRICES VISUALISATIONS -------------------------------------------------
gplot <- Plot.NetMat(data = Region_df, method = "HMSC", 
                     TreatmentOrder = c("Vermont", "Maine", "Temperate Broadleaf & Mixed Forests")) 
ggsave(gplot, file = file.path(DirEx.Region, "Broad-HMSC-Matrix.png"), height = 32, width = 32, units = "cm") 
gplot <- Plot.NetMat(data = Region_df, method = "IF_REM", 
                     TreatmentOrder = c("Vermont", "Maine", "Temperate Broadleaf & Mixed Forests")) 
ggsave(gplot, file = file.path(DirEx.Region, "Broad-IF_REM-Matrix.png"), height = 16, width = 32, units = "cm") 
gplot <- Plot.NetMat(data = Region_df, method = "NETASSOC", 
                     TreatmentOrder = c("Vermont", "Maine", "Temperate Broadleaf & Mixed Forests")) 
ggsave(gplot, file = file.path(DirEx.Region, "Broad-NETASSOC-Matrix.png"), height = 16, width = 32, units = "cm") 
gplot <- Plot.NetMat(data = Region_df, method = "COCCUR", 
                     TreatmentOrder = c("Vermont", "Maine", "Temperate Broadleaf & Mixed Forests")) 
ggsave(gplot, file = file.path(DirEx.Region, "Broad-COCCUR-Matrix.png"), height = 16, width = 32, units = "cm") 
### NETWORKS VISUALISED -----------------------------------------------------
#### all methods ----
FName <- "Cross-Broad_FIA"
png(filename = file.path(DirEx.Region, paste0(FName, ".png")), width = 24, height = 32, units = "cm", res = 1000)
discard <- grep(colnames(Region_df), pattern = "HMSC")[grep(colnames(Region_df), pattern = "HMSC") %nin% grep(colnames(Region_df), pattern = "biomass")]
Plot.Network(Plot_df = Region_df[,-discard],  
             ModelOrder = c("COCCUR", "NETASSOC", "HMSC", "IF_REM"),
             TreatmentOrder = c("Vermont", "Maine", "Temperate Broadleaf & Mixed Forests"), 
             FName = FName,
             Dir = DirEx.YFDP)
dev.off()
#### HMSC only ----
FName <- "HMSC-Broad_FIA"
png(filename = file.path(DirEx.Region, paste0(FName, ".png")), width = 24, height = 24, units = "cm", res = 1000)
keep <- grep(colnames(Region_df), pattern = "HMSC")
Part_df <- Region_df[,1:2]
Region_df <- Region_df[,keep]
names <- strsplit(colnames(Region_df), split = "[.]")
colnames(Region_df) <- paste(lapply(names, "[[", 3), 
                           lapply(names, "[[", 2), 
                           lapply(names, "[[", 4),
                           sep = ".")
Region_df <- cbind(Part_df, Region_df)
Plot.Network(Plot_df = Region_df,  
             ModelOrder = c("presence_absence", "abundance", "biomass"),
             TreatmentOrder = c("Vermont", "Maine", "Temperate Broadleaf & Mixed Forests"), 
             FName = FName,
             Dir = DirEx.YFDP)
dev.off()

## FIA2: Yosemite, Temperate Broadleaf & Mixed Forests ---------------------
### DATA --------------------------------------------------------------------
Region_ls2 <- Region_ls
Region_ls2 <- Region_ls2[c(grep(names(Region_ls), pattern = "Yosemite"), 
                           grep(names(Region_ls), pattern = "Temperate Conifer Forests")
)]
Shared_spec <- Reduce(intersect, list(
  Region_ls2$COCCUR.Yosemite$Partner1,
  Region_ls2$`COCCUR.Temperate Conifer Forests`$Partner1
))
Region_limited <- Limit.Lists(Region_ls2, Shared_spec)
Region_df <- Eff.Data.Frame(List_ls = Region_limited)
### TOPOLOGY ----------------------------------------------------------------
Region_topo <- Calc.Topology(data = Region_df, Sig = TRUE, Model = "ALL", TreatmentOrder = c("Yosemite", "Temperate Conifer Forests"))
ggsave(Region_topo$plots$Strength, file = file.path(DirEx.Region, "Coni-Topo_Strength.png"), height = 32, width = 32, units = "cm") 
ggsave(Region_topo$plots$Eigenvector, file = file.path(DirEx.Region, "Coni-Topo_Eigenvector.png"), height = 32, width = 32, units = "cm") 
ggsave(Region_topo$plots$Modularity, file = file.path(DirEx.Region, "Coni-Topo_Modularity.png"), height = 24, width = 32, units = "cm") 
ggsave(Region_topo$plots$Nestedness, file = file.path(DirEx.Region, "Coni-Topo_Nestedness.png"), height = 24, width = 32, units = "cm")  
### MATRICES VISUALISATIONS -------------------------------------------------
gplot <- Plot.NetMat(data = Region_df, method = "HMSC", 
                     TreatmentOrder = c("Yosemite", "Temperate Conifer Forests", "Difference")) 
ggsave(gplot, file = file.path(DirEx.Region, "Coni-HMSC-Matrix.png"), height = 32, width = 32, units = "cm") 
gplot <- Plot.NetMat(data = Region_df, method = "IF_REM", 
                     TreatmentOrder = c("Yosemite", "Temperate Conifer Forests", "Difference"))
ggsave(gplot, file = file.path(DirEx.Region, "Coni-IF_REM-Matrix.png"), height = 16, width = 32, units = "cm") 
gplot <- Plot.NetMat(data = Region_df, method = "NETASSOC", 
                     TreatmentOrder = c("Yosemite", "Temperate Conifer Forests", "Difference"))
ggsave(gplot, file = file.path(DirEx.Region, "Coni-NETASSOC-Matrix.png"), height = 16, width = 32, units = "cm") 
gplot <- Plot.NetMat(data = Region_df, method = "COCCUR", 
                     TreatmentOrder = c("Yosemite", "Temperate Conifer Forests", "Difference"))
ggsave(gplot, file = file.path(DirEx.Region, "Coni-COCCUR-Matrix.png"), height = 16, width = 32, units = "cm") 
### NETWORKS VISUALISED -----------------------------------------------------
#### all methods ----
FName <- "Cross-Coni_FIA"
png(filename = file.path(DirEx.Region, paste0(FName, ".png")), width = 16, height = 32, units = "cm", res = 1000)
discard <- grep(colnames(Region_df), pattern = "HMSC")[grep(colnames(Region_df), pattern = "HMSC") %nin% grep(colnames(Region_df), pattern = "biomass")]
Plot.Network(Plot_df = Region_df[,-discard],  
             ModelOrder = c("COCCUR", "NETASSOC", "HMSC", "IF_REM"),
             TreatmentOrder = c("Yosemite", "Temperate Conifer Forests"), 
             FName = FName,
             Dir = DirEx.Region)
dev.off()
#### HMSC only ----
FName <- "HMSC-Coni_FIA"
png(filename = file.path(DirEx.Region, paste0(FName, ".png")), width = 16, height = 24, units = "cm", res = 1000)
keep <- grep(colnames(Region_df), pattern = "HMSC")
Part_df <- Region_df[,1:2]
Region_df <- Region_df[,keep]
names <- strsplit(colnames(Region_df), split = "[.]")
colnames(Region_df) <- paste(lapply(names, "[[", 3), 
                             lapply(names, "[[", 2), 
                             lapply(names, "[[", 4),
                             sep = ".")
Region_df <- cbind(Part_df, Region_df)
Plot.Network(Plot_df = Region_df,  
             ModelOrder = c("presence_absence", "abundance", "biomass"),
             TreatmentOrder = c("Yosemite", "Temperate Conifer Forests"), 
             FName = FName,
             Dir = DirEx.Region)
dev.off()

# CROSS-SCALE, WITHIN METHOD ===============================================
## for this, we need to show the same species across scales
message("############ CROSS-SCALE; WITHIN-METHOD")

## FIA: Vermont, Maine, Temperate Broadleaf --------------------------------
### DATA --------------------------------------------------------------------
Region_ls2 <- Region_ls
Region_ls2 <- Region_ls2[c(grep(names(Region_ls), pattern = "Vermont"), 
                           grep(names(Region_ls), pattern = "Maine"),
                           grep(names(Region_ls), pattern = "Temperate Broadleaf & Mixed Forests")
)]
Shared_spec <- Reduce(intersect, list(
  Region_ls2$COCCUR.Vermont$Partner1,
  Region_ls2$COCCUR.Maine$Partner1,
  Region_ls2$`COCCUR.Temperate Broadleaf & Mixed Forests`$Partner1
))
Region_limited <- Limit.Lists(Region_ls2, Shared_spec)
Region_df <- Eff.Data.Frame(List_ls = Region_limited)

HMSC_cols <- grep(colnames(Region_df), pattern = "HMSC")
TargetHMSC <- c(grep(colnames(Region_df), pattern = "biomass"), grep(colnames(Region_df), pattern = "diametre"))
Region_df <- Region_df[,-HMSC_cols[HMSC_cols %nin% TargetHMSC]]

### MATRICES VISUALISATIONS -------------------------------------------------
gplot <- Plot.NetMat.Cross(data = Region_df,
                           ALLRegionOrder = c("Vermont", "Maine", "Temperate Broadleaf & Mixed Forests"),
                           ALLMethodOrder = c("COCCUR", "NETASSOC", "HMSC", "IF_REM"))
ggsave(gplot, file = file.path(DirEx.Region, "Cross-Broad.png"), height = 54, width = 32, units = "cm") 
### NETWORKS VISUALISED -----------------------------------------------------
### this is already done above

### TOPOLOGY ----------------------------------------------------------------
### this is already done above

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
### TOPOLOGY ----------------------------------------------------------------
Coni_topo <- Calc.Topology(data = data, Sig = TRUE, Model = "ALL", TreatmentOrder = c("Pre-Fire", "Yosemite", "Temperate Conifer Forests"))
ggsave(Coni_topo$plots$Strength, file = file.path(Dir.Exports, "Coni-Topo_Strength.png"), height = 32, width = 32, units = "cm") 
ggsave(Coni_topo$plots$Eigenvector, file = file.path(Dir.Exports, "Coni-Topo_Eigenvector.png"), height = 32, width = 32, units = "cm") 
ggsave(Coni_topo$plots$Modularity, file = file.path(Dir.Exports, "Coni-Topo_Modularity.png"), height = 24, width = 32, units = "cm") 
ggsave(Coni_topo$plots$Nestedness, file = file.path(Dir.Exports, "Coni-Topo_Nestedness.png"), height = 24, width = 32, units = "cm")  
### MATRICES VISUALISATIONS -------------------------------------------------
gplot <- Plot.NetMat.Cross(data = data,
                           ALLRegionOrder = c("Pre-Fire", "Yosemite", "Temperate Conifer Forests"),
                           ALLMethodOrder = c("COCCUR", "NETASSOC", "HMSC", "IF_REM"))
ggsave(gplot, file = file.path(Dir.Exports, "Cross-Coni.png"), height = 54, width = 32, units = "cm") 
### NETWORKS VISUALISED -----------------------------------------------------
FName <- "Cross-Coni2"
png(filename = file.path(Dir.Exports, paste0(FName, ".png")), width = 24, height = 32, units = "cm", res = 1000)
Plot.Network(Plot_df = data,  
             ModelOrder = c("COCCUR", "NETASSOC", "HMSC", "IF_REM"),
             TreatmentOrder = c("Pre-Fire", "Yosemite", "Temperate Conifer Forests"), 
             FName = FName,
             Dir = Dir.Exports)
dev.off()

# CROSS-SCALE, ACROSS METHOD ================================================
message("############ CROSS-SCALE, ACROSS METHOD")
## for this, we need to show the same species per scale
## we only use signs here as magnitudes of effects aren't comparable

## FIA: Vermont, Maine, Temperate Conifer ----------------------------------
### DATA --------------------------------------------------------------------
### MATRICES VISUALISATIONS -------------------------------------------------
### TOPOLOGY ----------------------------------------------------------------

## YFDP: Pre-Fire & FIA: Yosemite, Temperate Broadleaf & Mixed Forests -----
### DATA --------------------------------------------------------------------
### MATRICES VISUALISATIONS -------------------------------------------------
### TOPOLOGY ----------------------------------------------------------------
