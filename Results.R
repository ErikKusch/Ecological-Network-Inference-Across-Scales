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
#### SPECIES IDENTITIES ####
YFDP_spec <- ReadOut.Species(YFDP_ls) # species which have been analysed

Shared_spec <- YFDP_spec

# Region_spec <- ReadOut.Species(Region_ls) # species which have been analysed
# Shared_spec <- Reduce(intersect, list(YFDP_spec, Region_spec)) # shared species between all analyses

#### FLATTENING LISTS ####
YFDP_ls <- Flatten.Lists(YFDP_ls)
# Region_ls <- Flatten.Lists(Region_ls)
# Region_ls <- BiomeNames.List(Region_ls)

# VISUALISATIONS ===========================================================

## SHARED SPECIES ----------------------------------------------------------
# this is where we will compare interactions/associations 1-1 between shared species of all analyses
### LIMITING RESULTS TO SHARED SPECIES BETWEEN DATA SETS #####
YFDP_limited <- Limit.Lists(YFDP_ls, Shared_spec)
# Region_limited <- Limit.Lists(Region_ls, Shared_spec)

### INTERACTION EXTRACTION ####
YFDP_df <- Eff.Data.Frame(List_ls = YFDP_limited)
# Region_df <- Eff.Data.Frame(List_ls = Region_limited)

### WITHIN-SCALE COMPARISON --------------------------
#### YFDP ####
gplot <- Plot.NetMat(data = YFDP_df, method = "HMSC") 
ggsave(gplot, file = file.path(DirEx.YFDP, "HMSC-Matrix.png"), height = 32, width = 32, units = "cm") 
gplot <- Plot.NetMat(data = YFDP_df, method = "IF_REM") 
ggsave(gplot, file = file.path(DirEx.YFDP, "IF_REM-Matrix.png"), height = 16, width = 32, units = "cm") 
gplot <- Plot.NetMat(data = YFDP_df, method = "NETASSOC") 
ggsave(gplot, file = file.path(DirEx.YFDP, "NETASSOC-Matrix.png"), height = 16, width = 32, units = "cm") 
gplot <- Plot.NetMat(data = YFDP_df, method = "COCCUR") 
ggsave(gplot, file = file.path(DirEx.YFDP, "COCCUR-Matrix.png"), height = 16, width = 32, units = "cm") 

#### FIA ####

### ACROSS-SCALE COMPARISON --------------------------

#### YFDP vs. FIA - Yosemite ####

#### YFDP vs. FIA - Temperate Mixed ####

#### FIA - VT - FIA - Whatever covers VT ####
#### FIA - Maine - FIA - Whatever covers Maine ####

#### FIA - Yosemite - FIA - Whatever covers Yosemite ####

### NETWORK PLOTTING ---------------------------------
YFDP_plots <- Plot.Network(Plot_df = YFDP_df)
Region_plots <- Plot.Network(Plot_df = Region_df)

#### HMSC model specification within scales ####
Header_vec <- c("YFDP")
Dir.WithinScale <- file.path(Dir.Exports, "Within-Scale")
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
Dir.WithinScale <- file.path(Dir.Exports, "Within-Scale")
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
Dir.WithinScale <- file.path(Dir.Exports, "Within-Scale")
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
Dir.CrossScale <- file.path(Dir.Exports, "Cross-Scale")
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


# TOPOLOGY =================================================================
# # this is where we will calculate whole-network topology metrics for comparison of entire networks

## WITHIN-SCALE COMPARISON --------------------------
## SHARED SPECIES ----------------------------------------------------------
# this is where we will compare interactions/associations 1-1 between shared species of all analyses

### WITHIN-SCALE COMPARISON --------------------------
#### YFDP ####

#### FIA ####

### ACROSS-SCALE COMPARISON --------------------------

#### YFDP vs. FIA - Yosemite ####

#### YFDP vs. FIA - Temperate Mixed ####

## INDIVIDUAL NETWORKS -----------------------------------------------------
YFDP_df <- Eff.Data.Frame(List_ls = YFDP_ls)
Region_df <- Eff.Data.Frame(List_ls = Region_ls)

### WITHIN-SCALE COMPARISON --------------------------
#### YFDP ####

#### FIA ####

### ACROSS-SCALE COMPARISON --------------------------

#### YFDP vs. FIA - Yosemite ####

#### YFDP vs. FIA - Temperate Mixed ####
