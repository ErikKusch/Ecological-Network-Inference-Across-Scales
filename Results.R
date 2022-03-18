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
## VISUALISATIONS ----------------------------------------------------------
gplot <- Plot.NetMat(data = YFDP_df, method = "HMSC") 
ggsave(gplot, file = file.path(DirEx.YFDP, "HMSC-Matrix.png"), height = 32, width = 32, units = "cm") 
gplot <- Plot.NetMat(data = YFDP_df, method = "IF_REM") 
ggsave(gplot, file = file.path(DirEx.YFDP, "IF_REM-Matrix.png"), height = 16, width = 32, units = "cm") 
gplot <- Plot.NetMat(data = YFDP_df, method = "NETASSOC") 
ggsave(gplot, file = file.path(DirEx.YFDP, "NETASSOC-Matrix.png"), height = 16, width = 32, units = "cm") 
gplot <- Plot.NetMat(data = YFDP_df, method = "COCCUR") 
ggsave(gplot, file = file.path(DirEx.YFDP, "COCCUR-Matrix.png"), height = 16, width = 32, units = "cm") 
## TOPOLOGY ----------------------------------------------------------------
YFDP_df <- YFDP_df[,-grep(pattern = "NETASSOC", colnames(YFDP_df))]
YFDP_topo <- Calc.Topology(data = YFDP_df, Sig = TRUE, Model = "ALL", TreatmentOrder = c("Pre-Fire", "Post-Fire"))
ggsave(YFDP_topo$plots$Strength, file = file.path(DirEx.YFDP, "Topo_Strength.png"), height = 32, width = 32, units = "cm") 
ggsave(YFDP_topo$plots$Eigenvector, file = file.path(DirEx.YFDP, "Topo_Eigenvector.png"), height = 32, width = 32, units = "cm") 
ggsave(YFDP_topo$plots$Modularity, file = file.path(DirEx.YFDP, "Topo_Modularity.png"), height = 24, width = 32, units = "cm") 
ggsave(YFDP_topo$plots$Nestedness, file = file.path(DirEx.YFDP, "Topo_Nestedness.png"), height = 24, width = 32, units = "cm") 
# YFDP_topo$numbers

# # WITHIN-FIA, WITHIN METHOD ================================================
# message("############ WITHIN-FIA; WITHIN-METHOD")
# ## FIA1: Vermont, Maine, Temperate Conifer ---------------------------------
# ### DATA ----
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
# ### VISUALISATIONS ----
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
# ### TOPOLOGY ----
Region_df <- Region_df[,-grep(pattern = "COCCUR", colnames(Region_df))]
Region_topo <- Calc.Topology(data = Region_df, Sig = TRUE, Model = "ALL", TreatmentOrder = c("Vermont", "Maine", "Temperate Broadleaf & Mixed Forests"))
ggsave(Region_topo$plots$Strength, file = file.path(DirEx.Region, "Broad-Topo_Strength.png"), height = 32, width = 32, units = "cm") 
ggsave(Region_topo$plots$Eigenvector, file = file.path(DirEx.Region, "Broad-Topo_Eigenvector.png"), height = 32, width = 32, units = "cm") 
ggsave(Region_topo$plots$Modularity, file = file.path(DirEx.Region, "Broad-Topo_Modularity.png"), height = 16, width = 32, units = "cm") 
ggsave(Region_topo$plots$Nestedness, file = file.path(DirEx.Region, "Broad-Topo_Nestedness.png"), height = 16, width = 32, units = "cm")  

# ## FIA2: Yosemite, Temperate Broadleaf & Mixed Forests ---------------------
# ### DATA ----
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
# ### VISUALISATIONS ----
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
# ### TOPOLOGY ----
Region_topo <- Calc.Topology(data = Region_df, Sig = TRUE, Model = "ALL", TreatmentOrder = c("Yosemite", "Temperate Conifer Forests"))
ggsave(Region_topo$plots$Strength, file = file.path(DirEx.Region, "Coni-Topo_Strength.png"), height = 32, width = 32, units = "cm") 
ggsave(Region_topo$plots$Eigenvector, file = file.path(DirEx.Region, "Coni-Topo_Eigenvector.png"), height = 32, width = 32, units = "cm") 
ggsave(Region_topo$plots$Modularity, file = file.path(DirEx.Region, "Coni-Topo_Modularity.png"), height = 24, width = 32, units = "cm") 
ggsave(Region_topo$plots$Nestedness, file = file.path(DirEx.Region, "Coni-Topo_Nestedness.png"), height = 24, width = 32, units = "cm")  

# CROSS-SCALE, WITHIN METHOD ===============================================
## for this, we need to show the same species across scales
message("############ CROSS-SCALE; WITHIN-METHOD")

## FIA: Vermont, Maine, Temperate Broadleaf --------------------------------
### DATA ----
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

### VISUALISATIONS ----
gplot <- Plot.NetMat.Cross(data = Region_df,
                           ALLRegionOrder = c("Vermont", "Maine", "Temperate Broadleaf & Mixed Forests"),
                           ALLMethodOrder = c("COCCUR", "NETASSOC", "HMSC", "IF_REM"))
ggsave(gplot, file = file.path(Dir.Exports, "Cross-Broad.png"), height = 54, width = 32, units = "cm") 
### TOPOLOGY ----
### this is already done above

## YFDP: Pre-Fire & FIA: Yosemite, Temperate Conifer -----------------------
### DATA ----
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
### VISUALISATIONS ----
gplot <- Plot.NetMat.Cross(data = data,
                           ALLRegionOrder = c("Pre-Fire", "Yosemite", "Temperate Conifer Forests"),
                           ALLMethodOrder = c("COCCUR", "NETASSOC", "HMSC", "IF_REM"))
ggsave(gplot, file = file.path(Dir.Exports, "Cross-Coni.png"), height = 54, width = 32, units = "cm") 

### TOPOLOGY ----
Coni_topo <- Calc.Topology(data = data, Sig = TRUE, Model = "ALL", TreatmentOrder = c("Pre-Fire", "Yosemite", "Temperate Conifer Forests"))
ggsave(Coni_topo$plots$Strength, file = file.path(Dir.Exports, "Coni-Topo_Strength.png"), height = 32, width = 32, units = "cm") 
ggsave(Coni_topo$plots$Eigenvector, file = file.path(Dir.Exports, "Coni-Topo_Eigenvector.png"), height = 32, width = 32, units = "cm") 
ggsave(Coni_topo$plots$Modularity, file = file.path(Dir.Exports, "Coni-Topo_Modularity.png"), height = 24, width = 32, units = "cm") 
ggsave(Coni_topo$plots$Nestedness, file = file.path(Dir.Exports, "Coni-Topo_Nestedness.png"), height = 24, width = 32, units = "cm")  

# CROSS-SCALE, ACROSS METHOD ================================================
message("############ CROSS-SCALE, ACROSS METHOD")
## for this, we need to show the same species per scale
## we only use signs here as magnitudes of effects aren't comparable

## FIA: Vermont, Maine, Temperate Conifer ----------------------------------
### DATA ----
### VISUALISATIONS ----
### TOPOLOGY ----

## YFDP: Pre-Fire & FIA: Yosemite, Temperate Broadleaf & Mixed Forests -----
### DATA ----
### VISUALISATIONS ----
### TOPOLOGY ----


# NETWORK PLOTTING =========================================================
message("############ NETWORK PLOTTING")
Dir.WithinScale <- file.path(Dir.Exports, "Networks")

YFDP_limited <- Limit.Lists(YFDP_ls, Shared_spec)
YFDP_df <- Eff.Data.Frame(List_ls = YFDP_limited)
YFDP_plots <- Plot.Network(Plot_df = YFDP_df)
Region_limited <- Limit.Lists(Region_ls, Shared_spec)
Region_df <- Eff.Data.Frame(List_ls = Region_limited)
Region_plots <- Plot.Network(Plot_df = Region_df)

#### HMSC model specification within scales ####
Header_vec <- c("YFDP")
if(!dir.exists(Dir.WithinScale)){dir.create(Dir.WithinScale)}
for(Plots_Iter in 1:length(Header_vec)){
  Plots.Compare(Compare_df = list(YFDP = YFDP_df,
                                  Region = Region_df),
                Header = Header_vec[Plots_Iter], SubsetHeader = "HMSC",
                X = "Model", SubsetX = NULL,
                Y = "Treatment", SubsetY = NULL,
                HMSCFilter = "ALL",
                FName = paste0("HMSC-Specifications_", Header_vec[Plots_Iter]),
                Dir = Dir.WithinScale)
}
# --> most information in diametre-driven HMSC models. Using those going forward

#### Methods within scales  ####
Header_vec <- c("YFDP", "Region")
if(!dir.exists(Dir.WithinScale)){dir.create(Dir.WithinScale)}
for(Plots_Iter in 1:length(Header_vec)){
  SubsetY <- ifelse(Plots_Iter == 1, "Pre-Fire", "Temperate Broadleaf & Mixed Forests")
  Plots.Compare(Compare_df = list(YFDP = YFDP_df,
                                  Region = Region_df),
                Header = Header_vec[Plots_Iter],
                SubsetHeader = NULL,
                X = "Method", SubsetX = NULL,
                Y = "Treatment", SubsetY = SubsetY,
                HMSCFilter = "diametre",
                FName = paste0("Within-Scale_", Header_vec[Plots_Iter]),
                Dir = Dir.WithinScale)
}

#### Methods across Treatments  ####
Header_vec <- c("YFDP", "Region")
if(!dir.exists(Dir.WithinScale)){dir.create(Dir.WithinScale)}
for(Plots_Iter in 1:length(Header_vec)){
  Plots.Compare(Compare_df = list(YFDP = YFDP_df,
                                  Region = Region_df),
                Header = Header_vec[Plots_Iter],
                SubsetHeader = NULL,
                X = "Treatment", SubsetX = NULL,
                Y = "Method", SubsetY = NULL,
                HMSCFilter = "diametre",
                FName = paste0("Across-Treatment_", Header_vec[Plots_Iter]),
                Dir = Dir.WithinScale)
}

#### Methods across Scales  ####
Header_vec <- c("COCCUR")
if(!dir.exists(Dir.CrossScale)){dir.create(Dir.CrossScale)}
for(Plots_Iter in 1:length(Header_vec)){
  Plots.Compare(Compare_df = list(YFDP = YFDP_df,
                                  Region = Region_df),
                Header = Header_vec[Plots_Iter],
                X = "Data Set", SubsetX = NULL,
                Y = "Treatment", SubsetY = c("Pre-Fire", "Temperate Broadleaf & Mixed Forests"),
                HMSCFilter = "diametre",
                FName = paste0("Cross-Scale_", Header_vec[Plots_Iter]),
                Dir = Dir.CrossScale)
}