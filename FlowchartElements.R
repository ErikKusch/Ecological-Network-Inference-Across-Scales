#' ####################################################################### #
#' PROJECT: [Data Simplification for Network Inference across Scales] 
#' CONTENTS: 
#'  - Plotting
#'  DEPENDENCIES:
#'  - 0 - Preamble.R
#'  - 0 - ShapeFiles.R
#'  - X - Functions_Plotting.R
#' AUTHOR: [Erik Kusch]
#' ####################################################################### #

# PREAMBLE =================================================================
rm(list=ls())
set.seed(42)

## Sourcing ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
source("0 - Preamble.R")
source("0 - ShapeFiles.R")


# DATA MANIPULATION =========================================================
US_shp <- CountryMask[CountryMask$NAME == "United States of America",]
US_shp <- crop(US_shp, extent(-124.7,-66.9,25.1,49))


# PLOTTING ==================================================================
plot(US_shp)
plot(FIA_shp[FIA_shp$BIOME == 4,], col = "forestgreen", lty = 0, add = TRUE)

plot(US_shp)
plot(FIA_shp[FIA_shp$BIOME == 5,], col = "forestgreen", lty = 0, add = TRUE)

plot(crop(US_shp, extent(-85,-66.9,25.1,49)))
plot(States_shp[States_shp$BIOME == 1000 | States_shp$BIOME == 1001,], col = "darkred", lty = 0, add = TRUE)

plot(StateMask[StateMask$name == "California",])
plot(States_shp[States_shp$BIOME == 1003,], col = "darkred", lty = 0, add = TRUE)

plot(US_shp)
