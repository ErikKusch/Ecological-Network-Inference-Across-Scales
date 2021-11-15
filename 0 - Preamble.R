#' ####################################################################### #
#' PROJECT: [Data Simplification for Network Inference across Scales] 
#' CONTENTS: 
#'  - Directory Establishment
#'  - Package Loading
#'  - Helper Functions
#'  DEPENDENCIES:
#'  - 
#' AUTHOR: [Erik Kusch]
#' ####################################################################### #

# DIRECTORIES ===============================================================
Dir.Base <- getwd() # read out the project directory
## DATA ---------------------------------------------------------------------
Dir.Data <- file.path(Dir.Base, "Data")
Dir.YFDP <- file.path(Dir.Data, "YFDP")
Dir.Region <- file.path(Dir.Data, "Regional")
Dir.FIA <- file.path(Dir.Data, "RAW_FIA")
Dir.Observations <- file.path(Dir.Data, "Observations")
Dir.Shapes <- file.path(Dir.Data, "Shapes")
DataDirs <- c(Dir.Data, Dir.Shapes, Dir.YFDP, Dir.Region, Dir.FIA, Dir.Observations)
CreateDir <- sapply(DataDirs, function(x) if(!dir.exists(x)) dir.create(x))
## EXPORTS ------------------------------------------------------------------
Dir.Exports <- file.path(Dir.Base, "Exports")
DirEx.YFDP <- file.path(Dir.Exports, "YFDP")
DirEx.Region <- file.path(Dir.Exports, "Region")
DirEx.Observations <- file.path(Dir.Exports, "Observations")
ExportDirs <- c(Dir.Exports, DirEx.YFDP, DirEx.Region, DirEx.Observations)
CreateDir <- sapply(ExportDirs, function(x) if(!dir.exists(x)) dir.create(x))
rm(list = c("CreateDir", "ExportDirs", "DataDirs"))

# PACKAGES ================================================================
## CRAN -------------------------------------------------------------------
install.load.package <- function(x) {
  if (!require(x, character.only = TRUE))
    install.packages(x, repos='http://cran.us.r-project.org')
  require(x, character.only = TRUE)
}

package_vec <- c(
  "devtools", # needed for non-cran packages further down
  "rgeos",
  "readxl", # for reading xlsx files
  "sp", # for handling spatialpolygondataframes
  "rgdal", # for loading shapefiles of species ranges
  "raster", # for storing data as rasters
  "ncdf4", # for ncdf namespace when loading nertcdf files
  "fasterize", # for establishing richness maps in a timely manner
  "sf", # for use of SpatialpolygonsDataFrame objects in fasterize
  "gimms", # for downloading the reference raster for NDVI data
  "ggplot2", # for plotting various things
  "dplyr", # for data manipulation
  "stringr", # for padding numbers
  "rFIA", # for downloading and using Forest Inventory Analysis (FIA) data
  "rstan", # for access to stan
  "pheatmap", # for heatmaps of interactions
  "ape", # for calculating phylogenetic distances
  "qgraph", # for network visualisation
  "gawdis", # for balanced gower distance
  "reshape2", # for community matrix generation from abundance list
  "phytools", # for averaging of phylogenies
  "Hmsc", # for HMSC models
  "colorspace", # for HMSC evaluation
  "vioplot", # for HMSC evaluation
  "snow", # for HMSC evaluation
  "corrplot", # for HMSC evaluation
  "writexl",
  "cooccur",
  "netassoc",
  "BIEN"
)
sapply(package_vec, install.load.package)

## KrigR ------------------------------------------------------------------
if("KrigR" %in% rownames(installed.packages()) == FALSE){ # KrigR check
  Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true")
  devtools::install_github("https://github.com/ErikKusch/KrigR", force = TRUE)
}
library(KrigR) 
try(source("X - PersonalSettings.R")) # I do this here to specify number of cores and API credentials and am thus not sharing this file
# CDS API (needed for ERA5-Land downloads)
if(!exists("API_Key") | !exists("API_User")){ # CS API check: if CDS API credentials have not been specified elsewhere
  API_User <- readline(prompt = "Please enter your Climate Data Store API user number and hit ENTER.")
  API_Key <- readline(prompt = "Please enter your Climate Data Store API key number and hit ENTER.")
} # end of CDS API check
# NUMBER OF CORES
if(!exists("numberOfCores")){ # Core check: if number of cores for parallel processing has not been set yet
  numberOfCores <- as.numeric(readline(prompt = paste("How many cores do you want to allocate to these processes? Your machine has", parallel::detectCores())))
} # end of Core check

## PhyloMaker -------------------------------------------------------------
if("V.PhyloMaker" %in% rownames(installed.packages()) == FALSE){ # PhyloMaker check
  Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true")
  devtools::install_github("jinyizju/V.PhyloMaker")
}
library(V.PhyloMaker) 

## Rethinking -------------------------------------------------------------
if("rethinking" %in% rownames(installed.packages()) == FALSE){ # PhyloMaker check
  Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true")
  devtools::install_github("stan-dev/cmdstanr")
  install.packages(c("coda","mvtnorm","devtools","loo","dagitty"))
  devtools::install_github("rmcelreath/rethinking")
}
library(rethinking)

# FUNCTIONALITY =============================================================
`%nin%` <- Negate(`%in%`)

hush <- function(code){
  sink("NUL") # use /dev/null in UNIX
  tmp = code
  sink()
  return(tmp)
}

Sort.DF <- function(Data = NULL, Column = NULL){
  Data[order(Data[ , Column] ), ]
}
