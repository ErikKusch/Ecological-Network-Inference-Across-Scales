#' ####################################################################### #
#' PROJECT: [Data Simplification for Network Inference across Scales] 
#' CONTENTS: Reference raster and Shapefiles
#'  DEPENDENCIES:
#'  - 0 - Preamble.R must be run prior
#' AUTHOR: [Erik Kusch]
#' ####################################################################### #

#### LAND MASK (for masking species in the sea which are terrestrial and marine) -----------------------------------------------------------
if(!file.exists(file.path(Dir.Shapes, "LandMask.zip"))){ # if land mask has not been downloaded yet
  download.file("https://www.naturalearthdata.com/http//www.naturalearthdata.com/download/10m/physical/ne_10m_land.zip", destfile = paste(Dir.Shapes, "LandMask.zip", sep="/")) # download cultural vector
  unzip(paste(Dir.Shapes, "LandMask.zip", sep="/"), exdir = Dir.Shapes) # unzip the data
}
LandMask <- readOGR(Dir.Shapes, "ne_10m_land", verbose = FALSE) # read land mask in

#### ECOREGIONS -----------------------------------------------------------
if(!file.exists(file.path(Dir.Shapes, "WWF_ecoregions"))){
  download.file("http://assets.worldwildlife.org/publications/15/files/original/official_teow.zip", destfile = file.path(Dir.Shapes, "wwf_ecoregions.zip"))
  unzip(file.path(Dir.Shapes, "wwf_ecoregions.zip"), exdir = file.path(Dir.Shapes, "WWF_ecoregions"))
}
EcoregionsMask <- readOGR(file.path(Dir.Shapes, "WWF_ecoregions", "official", "wwf_terr_ecos.shp"), verbose = FALSE) # loading shapefile for wwf ecoregions
### vectors for WWF region naming
Abbr_Realms <- levels((EcoregionsMask$REALM))
Full_Realms <- c("Australasia", "Antarctic", "Afrotropics", "IndoMalay", "Nearctic", "Neotropics", "Oceania", "Palearctic")
Abbr_Biomes <- c(1:14, 99, 98)
Full_Biomes <- c("Tropical & Subtropical Moist Broadleaf Forests",
                 "Tropical & Subtropical Dry Broadleaf Forests",
                 "Tropical & Subtropical Coniferous Forests",
                 "Temperate Broadleaf & Mixed Forests",
                 "Temperate Conifer Forests",
                 "Boreal Forests/Taiga",
                 "Tropical & Subtropical Grasslands, Savannas & Shrublands",
                 "Temperate Grasslands, Savannas & Shrublands",
                 "Flooded Grasslands & Savannas",
                 "Montane Grasslands & Shrublands",
                 "Tundra",
                 "Mediterranean Forests, Woodlands & Scrub",
                 "Deserts & Xeric Shrublands",
                 "Mangroves",
                 "Snow-Covered/Barren",
                 "Limnic Bodies of Water")

#### FIA REGIONS -----------------------------------------------------------
if(!file.exists(file.path(Dir.Shapes, "FIAMask.rds"))){
  US_shp <- CountryMask[which(CountryMask$NAME == "United States of America"),]
  FIA_shp <- crop(EcoregionsMask, US_shp) # cropping to
  saveRDS(FIA_shp, file.path(Dir.Shapes, "FIAMask.rds"))
}else{
  FIA_shp <- readRDS(file.path(Dir.Shapes, "FIAMask.rds"))
}

#### ERA5-LAND MASK --------------------------------------------------------
if(!file.exists(file.path(Dir.Shapes, "LandMask.nc"))) { # if not downloaded yet
  LandMaskERA5 <- download_ERA(
    Variable = "2m_temperature",
    DataSet = "era5-land",
    DateStart = "1995-01-01",
    DateStop = "1995-01-01",
    TResolution = "day",
    TStep = 1,
    FileName = "LandMask",
    Dir = Dir.Base,
    API_User = API_User,
    API_Key = API_Key
  )
}else{
  LandMaskERA5 <- raster(file.path(Dir.Shapes, "LandMask.nc"))
}
