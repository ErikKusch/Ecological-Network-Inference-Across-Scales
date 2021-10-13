# ####################################################################### #
# PROJECT: [PhD; 2D - SDM NETWORKS] 
# CONTENTS: Obtaining global SDM outputs and generate species-association networks from them
# AUTHOR: Erik Kusch
# EDIT: 17/08/2020
# ####################################################################### #
############## RANGE LOADING -------------------------------------------------
if(!file.exists(file.path(Dir.Ranges, "Ranges_ls.RData")) & !file.exists(file.path(Dir.Ranges, "Rasters_ls.RData"))){
  ## IUCN ----
  Dirs_vec <- c(Dir.Ranges.Amphibians, Dir.Ranges.Reptiles, Dir.Ranges.Mammals) # IUCN directories
  Ranges_ls <- as.list(rep(NA, length(Dirs_vec))) # List for range shapefiles in R
  Rasters_ls <- as.list(rep(NA, length(Dirs_vec))) # List for richness rasters in R
  for(Ranges_Iter in 1:length(Dirs_vec)){ # loop over all IUCN directories
    Name <- strsplit(x = Dirs_vec[Ranges_Iter], split = paste0(Dir.Ranges, "/"))[[1]][2] # isolate name of directory
    Shp <- readOGR(dsn = Dirs_vec[Ranges_Iter], verbose = TRUE) # load range shape data
    Shp <- Shp[which(Shp$terrestial == "true"), ] # mask for terrestrial class
    Ranges_ls[[Ranges_Iter]] <- Shp # save shapefiles to list
    if(!file.exists(file.path(Dir.Ranges, paste0("Richness_", Name, ".nc")))){ # if richness raster has not been established yet
      RichRas <- fasterize(st_as_sf(Shp), raster = Reference_ras, field = NULL, fun = "count") # establish richness raster
      RichRas <- mask(RichRas, LandMask) # maks for land mask 
      RichRas <- mask(RichRas, LakeMask, inverse = TRUE) # mask for lakes
      writeRaster(x = RichRas, filename = file.path(Dir.Ranges, paste0("Richness_", Name)), format="CDF") # write raster with directory name to range directory
    }else{ # if richness raster already exists
      RichRas <-  raster(file.path(Dir.Ranges, paste0("Richness_", Name, ".nc"))) # load richness raster
    }
    Rasters_ls[[Ranges_Iter]] <- RichRas # save richness raster to list
  }
  names(Ranges_ls) <- gsub(pattern = paste0(Dir.Ranges, "/"), replacement = "", x = Dirs_vec) # set names of list
  names(Rasters_ls) <- gsub(pattern = paste0(Dir.Ranges, "/"), replacement = "", x = Dirs_vec) # set names of list
  rm(list = c("Dirs_vec", "Name", "Shp", "RichRas", "Ranges_Iter")) # clean environment
  
  ## BIRDLIFE ----
  Birds_shp <- sf::st_read(dsn = file.path(Dir.Ranges.Birds, "BOTW.gdb")) # load BirdLife data
  classes <- class(Birds_shp$Shape[1])[[1]] # read class of first shape object
  for(Classes_Iter in 2:length(Birds_shp$Shape)){ # loop over all shapes in the shapefile
    classes <- c(classes, class(Birds_shp$Shape[Classes_Iter])[[1]]) # append the class of current object
  }
  Birds_shp <- Birds_shp[which(classes == "sfc_MULTIPOLYGON"), ] # retain MULTIPOLYGON objects, this gets rid of MULTISURFACE objects which fasterize can't handle
  if(!file.exists(file.path(Dir.Ranges, "Richness_BIRDS.nc"))){ # if richness raster has not been established yet
    Birdsrich_ras <- fasterize(Birds_shp, raster = Reference_ras, field = NULL, fun = "count") # establish richness raster
    Birdsrich_ras <- mask(Birdsrich_ras, LandMask) # maks for land mask 
    Birdsrich_ras <- mask(Birdsrich_ras, LakeMask, inverse = TRUE) # mask for lakes
    writeRaster(x = Birdsrich_ras, filename = file.path(Dir.Ranges, "Richness_BIRDS"), format="CDF") # write raster with
  }else{ # if richness raster already exists
    Birdsrich_ras <-  raster(file.path(Dir.Ranges, paste0("Richness_BIRDS.nc"))) # load richness raster
  }
  
  ## BIEN ----
  
  
  ## SAVING TO DISK ----
  Ranges_ls[["BIRDS"]] <- Birds_shp # save to list
  save(Ranges_ls, file = file.path(Dir.Ranges, "Ranges_ls.RData"))
  Rasters_ls[["Total"]] <- sum(stack(Rasters_ls), na.rm = TRUE) # calculate total species richness
  Rasters_ls[["Total"]] <- mask(Rasters_ls[["Total"]], LandMask) # maks for land mask 
  Rasters_ls[["Total"]] <- mask(Rasters_ls[["Total"]], LakeMask, inverse = TRUE) # mask for lakes
  Rasters_ls[["BIRDS"]] <- Birdsrich_ras # save to list
  save(Rasters_ls, file = file.path(Dir.Ranges, "Rasters_ls.RData"))
  rm(list = c("Birds_shp", "Birdsrich_ras", "classes", "Classes_Iter")) # cleaning up environment
}else{
  load(file.path(Dir.Ranges, "Ranges_ls.RData"))
  load(file.path(Dir.Ranges, "Rasters_ls.RData"))
}

############## PLOTTING RICHNESS MAPS ----------------------------------------
## MAPVIEW HTML ----

### PICK IT UP HERE!! - SHOW BIOMES AND REALMS AND PROTECTED AREAS ON THE MAP AS WELL AS RIVERS

Ecoregions_sf <- st_as_sf(EcoregionsMask) 

Realms_sf <- Ecoregions_sf %>% 
  group_by(REALM) %>%
  summarise(geometry = sf::st_union(geometry)) %>%
  ungroup()
Realms_sf <- Realms_sf[which(Realms_sf$REALM %in% Abbr_Realms), ]
Realms_sf$REALM <- Full_Realms

Biomes_sf <- Ecoregions_sf %>% 
  group_by(BIOME) %>%
  summarise(geometry = sf::st_union(geometry)) %>%
  ungroup()
Biomes_sf <- Biomes_sf[which(Biomes_sf$BIOME %in% Abbr_Biomes), ]
Biomes_sf$BIOME <- Full_Biomes

# Union_sf <- st_union(Realms_sf, Biomes_sf)
RealmBiomes_sf <- st_intersection(Realms_sf, Biomes_sf)
RealmBiomes_sf$REALMBIOMES <- paste(RealmBiomes_sf$REALM, RealmBiomes_sf$BIOME, sep=" - ")
RealmBiomes_sf <- RealmBiomes_sf[, c(-1,-2)]

# ProtectedAreas_mv <- mapview(ProtectedAreasMask, col.regions = "green", color = "black", alpha.regions = 0.3)
# Countries_mv <- mapview(CountryMask, color = "black", alpha.regions = 0)
Rivers_mv <- mapview(RiversMask, color = "blue", alpha.regions = 0)
# RealmBiomes_mv <- mapview(RealmBiomes_sf, zcol = "REALM", labels = paste(Intersec$REALM, Intersec$BIOME, sep=" - "))
RealmBiomes_mv <- mapview(RealmBiomes_sf, color = "black")
# Biomes_mv <- mapview(Biomes_sf, color = "red", legend = FALSE, alpha.regions = 0)
Amphibians_mv <- mapview(layer.name = "Amphibian Species Richness", Rasters_ls$AMPHIBIANS, legend = TRUE, 
                         maxpixels =  ncell(Reference_ras), na.color = "#FFFFFF00") 
Reptiles_mv <- mapview(layer.name = "Reptilian Species Richness", Rasters_ls$REPTILES, legend = TRUE, 
                       maxpixels =  ncell(Reference_ras), na.color = "#FFFFFF00")
Mammals_mv <- mapview(layer.name = "Mammalian Species Richness", Rasters_ls$TERRESTRIAL_MAMMALS, legend = TRUE, 
                      maxpixels =  ncell(Reference_ras), na.color = "#FFFFFF00")
Birds_mv <- mapview(layer.name = "Avian Species Richness", Rasters_ls$BIRDS, legend = TRUE, 
                    maxpixels =  ncell(Reference_ras), na.color = "#FFFFFF00")
# Plants_mv <- mapview(layer.name = "Plant Species Richness", Rasters_ls$PLANTS, legend = TRUE, 
#                     maxpixels =  ncell(Reference_ras), na.color = "#FFFFFF00")
Total_mv <- mapview(layer.name = "Total Species Richness", Rasters_ls$Total, legend = TRUE, 
                    maxpixels =  ncell(Reference_ras), na.color = "#FFFFFF00")
Richness_mv <- Countries_mv + Amphibians_mv + Reptiles_mv + Mammals_mv + Birds_mv + Total_mv # Combine all maps into 1 map
mapshot(Richness_mv, url = paste0(Dir.Richness.Current, "/RICHNESSMaps.html")) # export to disk
rm(list = c("Countries_mv", "Amphibians_mv", "Reptiles_mv", "Mammals_mv", "Birds_mv", "Total_mv", "Richness_mv"))

## JPEGS ----
col.mapview <- got(n = max(unlist(lapply(X = Rasters_ls, FUN = maxValue))), alpha = 1, begin = 0, end = 1, direction = 1, option = "daenerys")
Title_vec <- c("Amphibian", "Reptilian", "Mammalian", "Avian", 
               # "Plant", 
               "Total")
for(Plot_Iter in 1:length(Rasters_ls)){
  jpeg(file=paste(Dir.Richness.Current, "/", Title_vec[[Plot_Iter]], "_Richness.jpeg", sep = ""), width = 16, height = 12, units = "cm", quality = 100, res = 1000)
  plot(LandMask, col = "grey", bg = "black", main = paste(Title_vec[[Plot_Iter]], "Species Richness"))
  plot(Rasters_ls[[Plot_Iter]], col = col.mapview, add = TRUE, horizontal = TRUE, 
       legend.args = list(text='Number of Species'))
  plot(LakeMask, col = "black", add = TRUE)
  dev.off()
}
rm(list = c("Plot_Iter", "Title_vec", "col.mapview"))

####### SPECIES LIST -----------------------------------------------------
if(!file.exists(file.path(Dir.Ranges, "Species_df.RData"))){
  SpeciesNames_vec <- NA
  DataSets_vec <- NA
  for(Names_Iter in 1:length(Ranges_ls)){ # loop over all range data sets
    if(names(Ranges_ls)[[Names_Iter]] != "BIRDS"){ # BIRD data is stored differently than IUCN data
      SpeciesNamesAdd_vec <- as.character(Ranges_ls[[Names_Iter]]$binomial) # extract species names
    }else{
      SpeciesNamesAdd_vec <- as.character(Ranges_ls[[Names_Iter]]$SCINAME) # extract species names
    }
    SpeciesNames_vec <- c(SpeciesNames_vec, SpeciesNamesAdd_vec) # combine species names
    DataSets_vec <- c(DataSets_vec, rep(Names_Iter, length(SpeciesNamesAdd_vec))) # note range data set (i.e. place in Ranges_ls)
  }
  # Make Species data frame listing names, IDs, and data sets
  Species_df <- data.frame(ID = rep(NA, length(DataSets_vec)),
                           Name = SpeciesNames_vec,
                           Group = DataSets_vec)
  Species_df <- Species_df[-1, ] # remove initial NA row
  Species_df <- Species_df[!duplicated(Species_df), ] # remove multiple mentions of the same species
  Species_df <- Species_df[order(Species_df$Name),] # sort by alphabet
  Species_df$ID <- 1:dim(Species_df)[1] # Species ID for later analysis
  save(Species_df, file = file.path(Dir.Ranges, "Species_df.RData"))
  rm(list = c("DataSets_vec", "SpeciesNames_vec", "Names_Iter"))
}else{
  load(file.path(Dir.Ranges, "Species_df.RData"))
}

####### ESTABLISH CO-OCCURRENCES -----------------------------------------------------
### OCCURRENCE DATA FRAME BIOMES IN REALMS ----
if(!file.exists(file.path(Dir.Ranges, "SpeciesCells_df.RData"))){
  if(!file.exists(file.path(Dir.Ranges, "SpeciesCells_ls.RData"))){
    Realms_ras <- fasterize(sf = st_as_sf(EcoregionsMask), raster = Reference_ras, field = "REALM", fun = "first")
    Biomes_ras <- fasterize(sf = st_as_sf(EcoregionsMask), raster = Reference_ras, field = "BIOME", fun = "first")
    SpeciesCells_ls <- as.list(rep(NA, dim(Species_df)[1])) # a list which will hold sparse data frames for each species and it's CellIDs
    names(SpeciesCells_ls) <- Species_df$Name # set names of the list positions
    Cells_pb <- txtProgressBar(min = 0, max = dim(Species_df)[1], style = 3) # make progress bar
    for(Cells_Iter in Species_df$ID){ # loop over all species IDs
      DataSet <- Ranges_ls[[Species_df$Group[[Cells_Iter]]]] # isolate the data set of the current species
      if(names(Ranges_ls)[[Species_df$Group[Cells_Iter]]] != "BIRDS"){ # BIRD data is stored differently than IUCN data
        Polys <- which(as.character(DataSet$binomial) == as.character(Species_df$Name[[Cells_Iter]])) # identify the polygon position(s)
      }else{
        Polys <- which(as.character(DataSet$SCINAME) == as.character(Species_df$Name[[Cells_Iter]])) # identify the polygon position(s)
      }
      Presence_ras <- fasterize(st_as_sf(DataSet[Polys,]), raster = Reference_ras, field = NULL, fun = "max") # raserize species range
      Cells <- which(!is.na(values(Presence_ras))) # identify cells with data
      # create data frame with species ID and CellIDs
      SpeciesCells_df <- data.frame(SpeciesID = rep(Species_df$ID[[Cells_Iter]], length(Cells)), 
                                    CellID = Cells,
                                    RealmID = values(Realms_ras)[Cells],
                                    BiomeID = values(Biomes_ras)[Cells])
      SpeciesCells_ls[[Cells_Iter]] <- SpeciesCells_df # save dta frame to list
      print(as.character(Species_df$Name[[Cells_Iter]]))
      setTxtProgressBar(Cells_pb, Cells_Iter) # update progress bar
    }
    rm(list = c("SpeciesCells_df", "Cells_pb", "Names_Iter", "Presence_ras", "Polys", "SpeciesNamesAdd_vec", "Cells_Iter", "Cells", "DataSet"))
    save(SpeciesCells_ls, file = file.path(Dir.Ranges, "SpeciesCells_ls.RData"))
  }else{
    load(file.path(Dir.Ranges, "SpeciesCells_ls.RData"))
  }
  SpeciesCells_df <- bind_rows(SpeciesCells_ls, .id = 'column_label')
  rm(SpeciesCells_ls)
  SpeciesCells_df <- na.omit(SpeciesCells_df)
  SpeciesCells_df$RealmID <- str_pad(as.character(SpeciesCells_df$RealmID), 2, "left","0")
  SpeciesCells_df$BiomeID <- str_pad(as.character(SpeciesCells_df$BiomeID), 2, "left","0")
  SpeciesCells_df$BiomesInRealms<- with(SpeciesCells_df, paste(RealmID, BiomeID, sep ="_"))
  save(SpeciesCells_df, file = file.path(Dir.Ranges, "SpeciesCells_df.RData"))
}else{
  load(file.path(Dir.Ranges, "SpeciesCells_df.RData"))
}

### COMMUNITY DATA FRAME BY BIOMES IN REALMS ----
if(!file.exists(file.path(Dir.Ranges, "SpeciesCells_df.RData"))){
  BiomeRealms <- sort(unique(SpeciesCells_df$BiomesInRealms))
  SiteCommunities_ls <- as.list(rep(NA, length(BiomeRealms)))
  names(SiteCommunities_ls) <- BiomeRealms
  for(Realm_Iter in 1:length(BiomeRealms)){
    print(BiomeRealms[Realm_Iter])
    CurrRealm <- BiomeRealms[Realm_Iter]
    CurrRows <- which(SpeciesCells_df$BiomesInRealms == CurrRealm)
    CurrCells <- SpeciesCells_df[CurrRows,]
    SiteCommunities_ls[[Realm_Iter]] <- table(CurrCells[,2:3])
  }
  save(SiteCommunities_ls, file = file.path(Dir.Ranges, "SiteCommunities_ls.RData"))
}else{
  load(file.path(Dir.Ranges, "SiteCommunities_ls.RData"))
}
rm(SpeciesCells_df)

### SPECIES ASSOCIATION NETWORKS ----
# install.packages("netassoc")
library(netassoc)

assoc_net <- make_netassoc_network(obs = SiteCommunities_ls[[2]], verbose = TRUE)
plot_netassoc_network(assoc_net$network_all)



### OUTPUT ----


