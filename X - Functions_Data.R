#' ####################################################################### #
#' PROJECT: [PhD; X - DATA FUNCTIONALITY] 
#' CONTENTS: 
#'  - Functionality for phylogenetic distance calculation for differing levels of taxonomic resolution
#'  DEPENDENCIES:
#'  - 
#' AUTHOR: [Erik Kusch]
#' ####################################################################### #

# PHYLOGENETIC DISTANCE ==================================================
FUN.PhyloDist <- function(SpeciesNames = NULL, Boot = 1e3){
  verbatim <- TRUE
  ## Object Creation ----
  Phylo_ls <- as.list(rep(NA, Boot)) # empty list to put distance matrices in
  PhyloDist_ls <- as.list(rep(NA, Boot)) # empty list to put distance matrices in
  pb <- txtProgressBar(min = 0, max = Boot, style = 3)
  
  ## Species Name Identification ----
  SpeciesNames <- unique(gsub(pattern = " " , replacement = "_", x = SpeciesNames)) # make sure that binary nomenclature is separated by underscores
  Status <- rep(3, length(SpeciesNames)) # this tracks at which level data was matched with phylogeny, 3 indexes failure
  ### Fully recognised species ----
  Recognised_Spec <- SpeciesNames[SpeciesNames %in% V.PhyloMaker::tips.info$species]
  Status[SpeciesNames %in% V.PhyloMaker::tips.info$species] <- 1 # match at species-level
  if(isTRUE(verbatim)){message(paste(length(Recognised_Spec), "are recognised at species level by V.Phylomaker"))}
  
  ### Genus-Level recognised species ----
  Unrecognised_Spec <- SpeciesNames[SpeciesNames %nin% Recognised_Spec] ## these are the species we were given, but aren't resolved to species-level
  Reported_Gen <- unlist(lapply(strsplit(x = Unrecognised_Spec, split = "_"), '[[', 1)) # genuses of non-fully recongised species
  Known_Gen <- unlist(lapply(strsplit(x = V.PhyloMaker::tips.info$species, split = "_"), '[[', 1)) # all genuses contained within phylogeny megatree
  Recognised_Gen <- Reported_Gen[Reported_Gen %in% Known_Gen]
  Status[SpeciesNames %in% Unrecognised_Spec[Reported_Gen %in% Known_Gen]] <- 2 # genus-level match
  if(isTRUE(verbatim)){message(paste(length(Recognised_Gen), "are recognised at genus level by V.Phylomaker"))}
  
  ### Non-recognised species ----
  Failed_Spec <- Unrecognised_Spec[Reported_Gen %nin% Known_Gen]
  if(isTRUE(verbatim)){message(paste(length(Failed_Spec), "are recognised neither at species or genus level by V.Phylomaker"))}
  if(length(Failed_Spec > 0) & isTRUE(verbatim)){
    message(paste("These species are:", paste(Failed_Spec, collapse = ", ")))
  }
  
  ## Phylogeny Creation Loop ----
  for(Phylo_Boot in 1:Boot){
    verbatim <- ifelse(Phylo_Boot == 1, TRUE, FALSE) # output of phylo function only on first iteration
    ### Random Sampling of Genus-level species ----
    Phylo_Spec <- c()
    if(length(Recognised_Gen) > 0){
      Sample_Spec <- V.PhyloMaker::tips.info$species ## all species we can sample from
      Sample_Gen <- Known_Gen ## all genuses for matching
      Sample_Gen <- Known_Gen[Sample_Spec %nin% Recognised_Spec] ## remove already fully resolved species from sample spectrum
      Sample_Spec <- Sample_Spec[Sample_Spec %nin% Recognised_Spec] ## remove already fully resolved species from sample spectrum
      ## loop over all sampling genuses
      Phylo_Spec <- NA
      for(Samp_Iter in 1:length(Recognised_Gen)){ 
        if(sum(Sample_Gen %in% Recognised_Gen[Samp_Iter]) == 0){## there are not enough species left over for us to randomly sample one for such a record
          Phylo_Spec[Samp_Iter] <- NA
        }else{
          Phylo_Spec[Samp_Iter] <- sample(Sample_Spec[Sample_Gen %in% Recognised_Gen[Samp_Iter]], 1)
        }
        Sample_Gen <- Sample_Gen[Sample_Spec %nin% Phylo_Spec[Samp_Iter]]
        Sample_Spec <- Sample_Spec[Sample_Spec %nin% Phylo_Spec[Samp_Iter]]
      }
      Phylo_Spec <- c(Recognised_Spec, Phylo_Spec)
      names(Phylo_Spec) <- c(SpeciesNames[Status == 1], SpeciesNames[Status == 2])
      Phylo_Spec <- na.omit(Phylo_Spec)
    }else{
      Phylo_Spec <- Recognised_Spec
    }
    ### Phylogeny Building ----
    ## prepare phylogeny building
    phylo_df <- data.frame(species = Phylo_Spec,
                           genus =  unlist(lapply(strsplit(Phylo_Spec, split = "_"), `[[`, 1)),
                           family =  V.PhyloMaker::tips.info$family[na.omit(base:: match(Phylo_Spec, V.PhyloMaker::tips.info$species))]
    )
    ## build phylogeny
    phylo_phylo <- hush(V.PhyloMaker::phylo.maker(sp.list = phylo_df, scenarios = "S3"))
    tree <- phylo_phylo$scenario.3
    tree$edge.length <-  tree$edge.length + 0.001
    tree$tip.label <- names(Phylo_Spec[match(tree$tip.label, Phylo_Spec)])
    
    ### Phylogenetic Distance ----
    phylo_dist <- ape::cophenetic.phylo(tree)
    
    ### Saving to objects ----
    Phylo_ls[[Phylo_Boot]] <- tree
    PhyloDist_ls[[Phylo_Boot]] <- phylo_dist
    setTxtProgressBar(pb, Phylo_Boot)
  }
  # Average out the phylogeny ----
  class(Phylo_ls) <- "multiPhylo"
  # Avg_phylo <- phytools::averageTree(Phylo_ls)
  ## calculate distance aggregates
  PhyloDist_mean <- apply(simplify2array(PhyloDist_ls), 1:2, mean) # build mean of all distances
  PhyloDist_SD <- apply(simplify2array(PhyloDist_ls), 1:2, sd) # build sd of all distances
  Avg_phylo <- ape::nj(PhyloDist_mean) # build average tree
  ## combine it into one object
  Phylo_ls <- list(Avg_Phylo = Avg_phylo,
                   Dist_Mean = PhyloDist_mean,
                   Dist_SD = PhyloDist_SD)
  ## Export of Results -----
  return(Phylo_ls) 
}

# FOREST INVENTORY ANALYSIS DATA =========================================
FUN.FIA <- function(states = c("DE","MD"), nCores = parallel::detectCores()/2){
  ### EXISTENCE CHECK
  Check_vec <- states %nin% substring(list.files(Dir.Plots.FIA), 1, 2)
  if(length(substring(list.files(Dir.Plots.FIA), 1, 2)) != 0){
    if(unique(substring(list.files(Dir.Plots.FIA), 1, 2) %in% states) != TRUE){stop("Your FIA directory contains data for more than the states you asked to analyse here. Please remove these or change the state argument here to include these files already present.")}
  }
  
  ## BIOME SHAPE PREPARATION
  FIAMerged_shp <- aggregate(FIA_shp, by = "BIOME") # Merge shapes according to biome type
  FIAMerged_shp@data$Names <- Full_Biomes[match(FIAMerged_shp@data$BIOME, Abbr_Biomes)] # Assign corresponding full text biome names
  
  ## CALCULATION OF FITNESS AS APPROXIMATED BY BIOMASS
  if(!file.exists(file.path(Dir.Plots, "FIABiomes_df.rds"))){
    # might need to run devtools::install_github('hunter-stanke/rFIA') to circumvent "Error in rbindlist(inTables..." as per https://github.com/hunter-stanke/rFIA/issues/7
    if(sum(Check_vec) != 0){
      FIA_df <- rFIA::getFIA(states = states[Check_vec], dir = Dir.Plots.FIA, nCores = nCores) # download FIA state data and save it to the FIA directory
    }else{
      FIA_df <- rFIA::readFIA(dir = Dir.Plots.FIA) # load all of the data in the FIA directory
    }
    FIABiomass_df <- biomass(db = FIA_df, # which data base to use
                             polys = FIAMerged_shp,
                             bySpecies = TRUE, # group by Species
                             byPlot = TRUE, # group by plot
                             nCores = nCores,
                             treeType = "live",
                             returnSpatial = TRUE
    )
    FIABiomass_df<- FIABiomass_df[which(FIABiomass_df$YEAR >= 1986 & FIABiomass_df$YEAR < 2020), ] # subsetting for year range we can cover with climate data with a five-year buffer on the front
    
    ## CLIMATE DATA EXTRACTION 
    Layer_seq <- seq.Date(as.Date("1981-01-01"), as.Date("2020-12-31"), by = "month")
    for(Clim_Iter in 1:length(ECV_vec)){
      print(ECV_vec[Clim_Iter])
      FIABiomass_df$XYZ <- NA
      colnames(FIABiomass_df)[ncol(FIABiomass_df)] <- ECV_vec[Clim_Iter]
      FIABiomass_df$XYZ <- NA
      colnames(FIABiomass_df)[ncol(FIABiomass_df)] <- paste0(ECV_vec[Clim_Iter], "_SD")
      FIABiomass_df$XYZ <- NA
      colnames(FIABiomass_df)[ncol(FIABiomass_df)] <- paste0(ECV_vec[Clim_Iter], "_UC")
      Extrac_temp <- raster::extract(stack(file.path(Dir.Plots, paste0("FIA_", ECV_vec[Clim_Iter], ".nc"))), 
                                     FIABiomass_df)
      Uncert_temp <- raster::extract(stack(file.path(Dir.Plots, paste0("FIA_UC", ECV_vec[Clim_Iter], ".nc"))), 
                                     FIABiomass_df)
      pb <- txtProgressBar(min = 0, max = nrow(Extrac_temp), style = 3) 
      for(Plot_Iter in 1:nrow(Extrac_temp)){
        ## only retain the last ten years leading up to data collection
        Need_seq <- seq.Date(as.Date(paste0(FIABiomass_df[Plot_Iter, ]$YEAR-10, "-01-01")), 
                             as.Date(paste0(FIABiomass_df[Plot_Iter, ]$YEAR, "-01-01")), 
                             by = "month")
        Time_seq <- Extrac_temp[Plot_Iter, which(Layer_seq %in% Need_seq)]
        Uncert_seq <- Uncert_temp[Plot_Iter, which(Layer_seq %in% Need_seq)]
        FIABiomass_df[Plot_Iter, ECV_vec[Clim_Iter]] <- mean(Time_seq, na.rm = TRUE)
        FIABiomass_df[Plot_Iter, paste0(ECV_vec[Clim_Iter], "_SD")] <- sd(Time_seq, na.rm = TRUE)
        FIABiomass_df[Plot_Iter, paste0(ECV_vec[Clim_Iter], "_UC")] <- mean(Uncert_seq, na.rm = TRUE)
        setTxtProgressBar(pb, Plot_Iter)
      }
    }
    FIABiomass_df <- na.omit(FIABiomass_df) # remove NAs (produced by era5-land not having data for certain plots)
    saveRDS(FIABiomass_df, file.path(Dir.Plots, "FIABiomes_df.rds"))
  }else{
    FIABiomass_df <- readRDS(file.path(Dir.Plots, "FIABiomes_df.rds"))
  }
  
  ## PHYLOGENY 
  if(!file.exists(file.path(Dir.Plots, "Phylogeny.RData"))){
    Phylo_ls <- FUN.PhyloDist(FIABiomass_df$SCIENTIFIC_NAME)
    save(Phylo_ls, file = file.path(Dir.Plots, "Phylogeny.RData"))
  }else{
    load(file.path(Dir.Plots, "Phylogeny.RData"))
  }
  Phylo_specs <- Phylo_ls$Avg_Phylo$tip.label
  Phylo_specs <- gsub(pattern = "_", replacement = " ", x =  Phylo_specs)
  Phylo_ls$Avg_Phylo$tip.label <- Phylo_specs
  ## limiting to recognised species
  FIABiomass_df <- FIABiomass_df[FIABiomass_df$SCIENTIFIC_NAME %in% Phylo_specs, ]
  
  ## SPLITTING INTO BIOMES
  FIASplit_ls <- split(FIABiomass_df, FIABiomass_df$polyID) # extract for each biome to separate data frame
  FIASplit_ls <- c(list(FIABiomass_df), FIASplit_ls)
  names(FIASplit_ls) <- c("ALL", FIAMerged_shp@data$Names) # apply correct names of biomes
  FIASplit_ls <- FIASplit_ls[-c((length(FIASplit_ls)-1):length(FIASplit_ls))] # remove 98 and 99 biome which is barren and limnic
  for(Biome_Iter in 1:length(FIASplit_ls)){
    BiomeName <- names(FIASplit_ls)[Biome_Iter]
    print(BiomeName)
    if(file.exists(file.path(Dir.Plots, paste0("FIABiome", Biome_Iter, ".RData")))){next()}
    FIAIter_df <- FIASplit_ls[[Biome_Iter]]
    
    ### Metadata_df ----
    Metadata_df <- as.data.frame(FIAIter_df[,c("pltID", "YEAR", paste0(rep(ECV_vec, each = 3), c("", "_SD", "_UC")))])[-15]
    colnames(Metadata_df)[1] <- "SiteID"
    Metadata_df <- Metadata_df[!duplicated(Metadata_df$SiteID), ]
    ### Fitness ----
    Weight_df <- as.data.frame(FIAIter_df[,c("pltID", "SCIENTIFIC_NAME", "BIO_ACRE")])[,-4]
    colnames(Weight_df) <- c("SiteID", "taxon", "value")
    ### FitCom ----
    Weight_cast <- reshape2::dcast(Weight_df, SiteID ~ taxon, value.var = 'value', fun.aggregate = mean)
    Weight_cast[Weight_cast == -Inf] <- 0
    Weight_cast[is.na(Weight_cast)] <- 0
    ### Community ----
    Abund_df <- as.data.frame(FIAIter_df[,c("pltID", "SCIENTIFIC_NAME", "nStems")])[,-4]
    colnames(Abund_df) <- c("SiteID", "taxon", "individual_nr")
    Abund_cast <- reshape2::dcast(Abund_df, SiteID ~ taxon, value.var = 'individual_nr', fun.aggregate = max)
    Abund_cast[Abund_cast == -Inf] <- 0
    Abund_cast <- Abund_cast[ , colnames(Abund_cast) %in% colnames(Weight_cast)]
    ## LIST MERGING ----
    ModelFrames_ls <- list(
      Fitness = Weight_df,
      FitCom = Weight_cast,
      Community = Abund_cast
    )
    ### Phylo_ls ----
    Phylo_Iter <- Phylo_ls
    Phylo_Iter$Avg_Phylo <- drop.tip(phy = Phylo_Iter$Avg_Phylo, 
                                     tip =  Phylo_Iter$Avg_Phylo$tip.label[Phylo_Iter$Avg_Phylo$tip.label %nin% colnames(ModelFrames_ls$FitCom)[-1]])
    Pos_Safe <- colnames(Phylo_Iter$Dist_Mean) %in% gsub(x = Phylo_Iter$Avg_Phylo$tip.label, pattern = " ", replacement = "_")
    Phylo_Iter$Dist_Mean <- Phylo_Iter$Dist_Mean[Pos_Safe, Pos_Safe]
    Phylo_Iter$Dist_SD <- Phylo_Iter$Dist_SD[Pos_Safe, Pos_Safe]
    save(BiomeName, Metadata_df, ModelFrames_ls, Phylo_Iter, 
         file = file.path(Dir.Plots, paste0("FIABiome", Biome_Iter, ".RData")))
  }
}

# FOREST INVENTORY ANALYSIS CLIMATE DATA =================================
FUN.FIACLIM <- function(ECV = NULL, Shp = NULL, Dir = NULL){
  PrecipFix <- ifelse(startsWith(ECV_vec[ECV], "total"), TRUE, FALSE)
  if(!file.exists(file.path(Dir, paste0("FIA_",  ECV_vec[ECV], ".nc")))){
    Temp_ras <- download_ERA(Variable = ECV_vec[ECV],
                             DateStart = "1981-01-01",
                             DateStop = "2020-12-31",
                             TResolution = "month",
                             TStep = 1,
                             API_Key = API_Key,
                             API_User = API_User,
                             Extent = Shp,
                             Dir = Dir,
                             FileName = paste0("FIA_",  ECV_vec[ECV]),
                             Cores = Cores,
                             TryDown = 42,
                             SingularDL = TRUE,
                             verbose = TRUE,
                             PrecipFix = PrecipFix
    )
  }
  if(!file.exists(file.path(Dir, paste0("FIA_UC",  ECV_vec[ECV], ".nc")))){
    Temp_ras <- download_ERA(Variable = ECV_vec[ECV],
                             DataSet = "era5",
                             Type = "ensemble_members",
                             DateStart = "1981-01-01",
                             DateStop = "2020-12-31",
                             TResolution = "month",
                             TStep = 1,
                             API_Key = API_Key,
                             API_User = API_User,
                             Extent = Shp,
                             Dir = Dir,
                             FileName = paste0("FIA_UC",  ECV_vec[ECV]),
                             Cores = parallel::detectCores(),
                             TryDown = 42,
                             SingularDL = TRUE,
                             verbose = TRUE,
                             PrecipFix = PrecipFix
    ) 
    print("Aggregating Ensembles")
    Temp_ras <- stackApply(Temp_ras, rep(1:(nlayers(Temp_ras)/10), each = 10), sd, progress = "text")
    writeRaster(x = Temp_ras, filename = file.path(Dir, paste0("FIA_UC",  ECV_vec[ECV], ".nc")), overwrite = TRUE)
  }
} 