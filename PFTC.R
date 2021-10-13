#' ######################################################################## #
#' PROJECT: [Data Simplification for Network Inference across Scales] 
#' CONTENTS: 
#'  - Generate PFTC species-association/-interaction networks
#'  DEPENDENCIES:
#'  - 0 - Preamble.R
#'  - X - Functions_Bayes.R
#'  - X - Functions_Plotting.R
#'  - PFTC3-Puna-PFTC5_Peru_2018-2020_LeafTraits_clean.csv (available from PFTC group efforts; https://osf.io/hjpwt/)
#'  - PU.10_PFTC3.10_2020_Peru_Coordinates.xlsx (available from PFTC group efforts; https://osf.io/uk85w/)
#' AUTHOR: [Erik Kusch]
#' ######################################################################## #

# PREAMBLE =================================================================
rm(list=ls())
set.seed(42)

## Sourcing ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
source("0 - Preamble.R")
source("0 - ShapeFiles.R")
source("X - Functions_Plotting.R")
source("X - Functions_Bayes.R")
source("X - Functions_Data.R")

## Bayes Settings ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
nSamples <- 7000
nWarmup <- 700
nChains <- 8

# DATA =====================================================================

## PHYLOGENY ---------------------------------------------------------------
#### Bootstrapped Phylogenetic Distance Creation 
if(!file.exists(file.path(Dir.PFTC, "Phylogeny.RData"))){ # check if phylogenetic distance has already been established
  download.file(url = "https://osf.io/hjpwt/download", 
                destfile = file.path(Dir.PFTC, "PFTC3-Puna-PFTC5_Peru_2018-2020_LeafTraits_clean.csv"), mode = "wb")
  Raw_df <- read.csv(file.path(Dir.PFTC, "PFTC3-Puna-PFTC5_Peru_2018-2020_LeafTraits_clean.csv")) # load pftc data
  
  Phylo_ls <- FUN.PhyloDist(SpeciesNames = Raw_df$taxon) # get phylogenetic distance with mixed species-/genus-pool
  save(Phylo_ls, file = file.path(Dir.PFTC, "Phylogeny.RData")) # save phylo dist and sd
  unlink(file.path(Dir.PFTC, "PFTC3-Puna-PFTC5_Peru_2018-2020_LeafTraits_clean.csv"))
}else{ # if phylo has already been established prior
  load(file.path(Dir.PFTC, "Phylogeny.RData")) # load object containing phylo distance
}
Phylo_Specs <- gsub(pattern = "_", replacement = " ", x =  Phylo_ls$Avg_Phylo$tip.label)

## META & CLIMATE DATA ----------------------------------------------------
ECV_vec <- c("2m_temperature", "volumetric_soil_water_layer_1", "total_precipitation", "potential_evaporation")
if(!file.exists(file.path(Dir.PFTC, "Metadata_df.csv"))){ # bioclimatic data not established yet
  download.file(url = "https://osf.io/uk85w/download", 
                destfile = file.path(Dir.PFTC, "PU.10_PFTC3.10_2020_Peru_Coordinates.xlsx"), mode = "wb")
  Metadata_df <- as.data.frame(readxl::read_xlsx(file.path(Dir.PFTC, "PU.10_PFTC3.10_2020_Peru_Coordinates.xlsx"))) # read metadata file, available here https://osf.io/uk85w/
  Metadata_df$SiteID <- with(Metadata_df, paste(Site, Treatment, PlotID, sep = "_")) # create ID column which combines site, treatment, and plotid
  Metadata_df$Year <- 2018 # set observation years to 2018 for all sites (that was the field season most of the data was collected throughout)
  colnames(Metadata_df)[5] <- "Lat"
  colnames(Metadata_df)[6] <- "Lon"
  
  ### THIS IS WHERE NEW AND DIFFERENT CLIMATE DATA NEEDS TO BE DOWNLOADED
  for(Variable_Iter in ECV_vec){
    if(Variable_Iter == "total_precipitation"){
      PrecipFix <- TRUE
    }else{
      PrecipFix <- FALSE
    }
    print(Variable_Iter)
    ## Download Anylsis Data at Highest Resolution
    Era5Data_ras <- download_ERA(Variable = Variable_Iter,
                                 DateStart = "1981-01-01",
                                 DateStop = paste0(unique(Metadata_df$Year),"-12-31"),
                                 TResolution = "month",
                                 TStep = 1,
                                 API_Key = API_Key,
                                 API_User = API_User,
                                 Extent = Metadata_df,
                                 ID = "SiteID",
                                 Buffer = 0.1,
                                 Dir = Dir.Data,
                                 FileName = paste0(Variable_Iter, "_Data"),
                                 Cores = Cores,
                                 TryDown = 42,
                                 SingularDL = TRUE,
                                 verbose = TRUE,
                                 PrecipFix = PrecipFix
    )
    ## Download Uncertainty Data
    Era5Uncert_ras <- download_ERA(Variable = Variable_Iter,
                                   DataSet = "era5",
                                   Type = "ensemble_members",
                                   DateStart = "1981-01-01",
                                   DateStop = paste0(unique(Metadata_df$Year),"-12-31"),
                                   TResolution = "month",
                                   TStep = 1,
                                   API_Key = API_Key,
                                   API_User = API_User,
                                   Extent = Metadata_df,
                                   Buffer = 0.5,
                                   ID = "SiteID",
                                   Dir = Dir.Data,
                                   FileName = paste0(Variable_Iter, "_Uncert.nc"),
                                   Cores = parallel::detectCores(),
                                   TryDown = 42,
                                   SingularDL = TRUE,
                                   verbose = TRUE,
                                   PrecipFix = PrecipFix
    )
    Era5Uncert_ras <- stackApply(Era5Uncert_ras, rep(1:(nlayers(Era5Uncert_ras)/10), each = 10), sd)
    
    ## Extract Data to Locations
    Extract_df <- cbind(Metadata_df$Lon, Metadata_df$Lat)
    colnames(Extract_df) <- c("x", "y")
    Extract_sp <- SpatialPoints(Extract_df)
    Data_mat <- raster::extract(Era5Data_ras, Extract_sp)
    Uncert_mat <- raster::extract(Era5Uncert_ras, Extract_sp)
    
    ## Save Mean Values to Original Data Source
    Metadata_df$XYZ <- NA
    colnames(Metadata_df)[ncol(Metadata_df)] <- Variable_Iter
    Metadata_df[, ncol(Metadata_df)] <- rowMeans(Data_mat)
    Metadata_df$XYZ <- NA
    colnames(Metadata_df)[ncol(Metadata_df)] <- paste0(Variable_Iter, "_SD")
    Metadata_df[, ncol(Metadata_df)] <- apply(Data_mat, 1, sd)
    Metadata_df$XYZ <- NA
    colnames(Metadata_df)[ncol(Metadata_df)] <- paste0(Variable_Iter, "_UC")
    Metadata_df[, ncol(Metadata_df)] <- rowMeans(Uncert_mat)

    ## writing result
    write.csv(Metadata_df, file.path(Dir.PFTC, "Metadata_df.csv")) 
  }
  unlink(file.path(Dir.PFTC, "PU.10_PFTC3.10_2020_Peru_Coordinates.xlsx"))
  unlink(file.path(Dir.Data, list.files(Dir.Data, pattern = ".nc")))
}else{ # bioclimatic data has been established
  Metadata_df <- read.csv(file.path(Dir.PFTC, "Metadata_df.csv"))[-1] # load the metadata with attached bioclimatic data
}

## MODEL DATA FRAMES -------------------------------------------------------
if(!file.exists(file.path(Dir.PFTC, "ModelFrames.RData"))){
  download.file(url = "https://osf.io/hjpwt/download", 
                destfile = file.path(Dir.PFTC, "PFTC3-Puna-PFTC5_Peru_2018-2020_LeafTraits_clean.csv"), mode = "wb")
  Raw_df <- read.csv(file.path(Dir.PFTC, "PFTC3-Puna-PFTC5_Peru_2018-2020_LeafTraits_clean.csv")) # load pftc data
  Raw_df <- Raw_df[Raw_df$taxon %in% Phylo_Specs, ] # limitting to phylogeny-recognised species
  ### FITNESS ####
  Raw_df <- Raw_df[Raw_df$trait == "dry_mass_g", c("site", "treatment", "plot_id", "taxon", "trait", "value")] # only dry biomass rows and select only relevant columns for speedier data handling 
  Raw_df$SiteID <- with(Raw_df, paste(site, treatment, plot_id, sep = "_")) # create SiteID index
  Weights_df <- Raw_df[,c("SiteID", "taxon", "value")] # limit data to necessary parts
  ### ABUNDANCE ####
  Raw_df <- read.csv(file.path(Dir.PFTC, "PFTC3-Puna-PFTC5_Peru_2018-2020_LeafTraits_clean.csv")) # load pftc data
  Raw_df <- Raw_df[Raw_df$taxon %in% Phylo_Specs, ] # limitting to phylogeny-recognised species
  Raw_df <- Raw_df[, c("site", "treatment", "plot_id", "taxon", "individual_nr")] # only dry biomass rows and select only relevant columns for speedier data handling 
  Raw_df$SiteID <- with(Raw_df, paste(site, treatment, plot_id, sep = "_")) # create SiteID index
  Abund_df <- Raw_df[,c("SiteID", "taxon", "individual_nr")] # limit data to necessary parts
  Abund_cast <- reshape2::dcast(Abund_df, SiteID ~ taxon, value.var = 'individual_nr', fun.aggregate = max)
  Abund_cast[Abund_cast == -Inf] <- 0
  ### LIST MERGING ####
  ModelFrames_ls <- list(
    Fitness = Weights_df,
    Community = Abund_cast
  )
  save(ModelFrames_ls, file = file.path(Dir.PFTC, "ModelFrames.RData"))
  unlink(file.path(Dir.PFTC, "PFTC3-Puna-PFTC5_Peru_2018-2020_LeafTraits_clean.csv"))
}else{
  load(file.path(Dir.PFTC, "ModelFrames.RData"))
}

## FUNCTIONAL TRAITS -------------------------------------------------------
## download PFTC data if needed
if(!file.exists(file.path(Dir.PFTC, "Traits_df.csv"))){
  download.file(url = "https://osf.io/hjpwt/download", 
                destfile = file.path(Dir.PFTC, "PFTC3-Puna-PFTC5_Peru_2018-2020_LeafTraits_clean.csv"), mode = "wb")
  Raw_df <- read.csv(file.path(Dir.PFTC, "PFTC3-Puna-PFTC5_Peru_2018-2020_LeafTraits_clean.csv")) #
  Raw_df$SiteID <- with(Raw_df, paste(site, treatment, plot_id, sep = "_")) # create SiteID index
  Raw_df <- Raw_df[Raw_df$taxon %in% Phylo_Specs, ] # limitting to phylogeny-recognised species
  Traits_df <- data.frame(
    Observation = Raw_df$id,
                        # Species = Raw_df$taxon,
                        Traits = Raw_df$trait,
                        Values = Raw_df$value
                        # SiteID = Raw_df$SiteID
  )
  ## make data frame wide for clearer representation of data
  Traits_df <- reshape(Traits_df, direction = "wide", 
                     idvar = "Observation", timevar = "Traits")
  colnames(Traits_df) <- gsub(pattern = "Values.", replacement = "", x = colnames(Traits_df))
  Traits_df$Species <- Raw_df$taxon[match(Traits_df$Observation, Raw_df$id)]
  Traits_df$SiteID <- Raw_df$SiteID[match(Traits_df$Observation, Raw_df$id)]
  ## Remove species not recognized by taxonomy here
  Traits_df <- Traits_df[Traits_df$Species %in% colnames(Phylo_ls[[1]]), ]
  write.csv(Traits_df, file = file.path(Dir.PFTC, "Traits_df.csv"))
  unlink(file.path(Dir.PFTC, "PFTC3-Puna-PFTC5_Peru_2018-2020_LeafTraits_clean.csv"))
}else{
  Traits_df <- read.csv(file.path(Dir.PFTC, "Traits_df.csv"))[ , -1]
}


# ANALYSIS =================================================================
## HMSC --------------------------------------------------------------------

## IF-REM ------------------------------------------------------------------
### INDEX ####
Index_df <- Weights_df
Index_df <- na.omit(Index_df) # remove NA rows
colnames(Index_df) <- c("SiteID", "Species", "Outcome")
## Remove species not recognized by taxonomy here
Index_df <- Index_df[Index_df$Species %in% colnames(Phylo_ls[[1]]), ]

### NEIGHBOURS  ####
## aggregating predictors of individual plants at plots to species-level means
Raw_df <- read.csv(file.path(Dir.PFTC, "PFTC3-Puna-PFTC5_Peru_2018-2020_LeafTraits_clean.csv")) # load pftc data
Raw_df <- Raw_df[ , c("site", "treatment", "plot_id", "taxon", "trait", "value")] # only keep individual numbers
Raw_df$SiteID <- with(Raw_df, paste(site, treatment, plot_id, sep = "_")) # create SiteID index
Raw_df <- Raw_df[Raw_df$trait == "dry_mass_g", ]
Neigh_df <- Raw_df[,c("SiteID", "taxon", "value")] # limit data to necessary parts
## Remove species not recognized by taxonomy here
Neigh_df <- Neigh_df[Neigh_df$taxon %in% colnames(Phylo_ls[[1]]), ]
## continue data handling
Neigh_Data_df <- aggregate.data.frame(x = Neigh_df$value, 
                                      by = list(Neigh_df$SiteID, Neigh_df$taxon), 
                                      FUN = "mean")
colnames(Neigh_Data_df) <- c("SiteID", "Species", "Predictor") # assign more readable column names
Neigh_SD_df <- aggregate.data.frame(x = Neigh_df$value, 
                                    by = list(Neigh_df$SiteID, Neigh_df$taxon), 
                                    FUN = "sd")
colnames(Neigh_SD_df) <- c("SiteID", "Species", "Predictor") # assign more readable column names
## make data frame wide for clearer representation of data
Neigh_Data_df <- reshape(data = Neigh_Data_df, direction = "wide", 
                         idvar = "SiteID", timevar = "Species")
Neigh_Data_df[is.na(Neigh_Data_df)] <- 0 # set all NAs to zero, DO NOT DO THIS WITH ANYTHING BUT ABUNDANCES
colnames(Neigh_Data_df) <- gsub(colnames(Neigh_Data_df), pattern = "Predictor.", replacement = "") # fix column names which reshape() messed up
Neigh_SD_df <- reshape(data = Neigh_SD_df, direction = "wide", 
                       idvar = "SiteID", timevar = "Species")
Neigh_SD_df[is.na(Neigh_SD_df)] <- 0 # set all NAs to zero, DO NOT DO THIS WITH ANYTHING BUT ABUNDANCES
colnames(Neigh_SD_df) <- gsub(colnames(Neigh_SD_df), pattern = "Predictor.", replacement = "") # fix column names which reshape() messed up
## combine data into final Neighbour list
Neigh_ls <- list(Mean = Neigh_Data_df,
                 SD = Neigh_SD_df) 

### ENVIRONMENT ---------------------------------------
## isolating only climate information
Envir_ls <- list(Mean = Metadata_df[ , c(9, seq(11, ncol(Metadata_df), 2))],
                 SD = Metadata_df[ , c(9, seq(12, ncol(Metadata_df), 2))] 
)
colnames(Envir_ls[[1]]) <- c("SiteID", ECV_vec)
colnames(Envir_ls[[2]]) <- c("SiteID", ECV_vec)

### TRAITS --------------------------------------------
## preparing data frame of traits and observations
Raw_df <- read.csv(file.path(Dir.PFTC, "PFTC3-Puna-PFTC5_Peru_2018-2020_LeafTraits_clean.csv")) # load pftc data
Raw_df$SiteID <- with(Raw_df, paste(site, treatment, plot_id, sep = "_")) # create SiteID index
PFTC_df <- data.frame(Observation = Raw_df$id,
                      # Species = Raw_df$taxon,
                      Traits = Raw_df$trait,
                      Values = Raw_df$value
                      # SiteID = Raw_df$SiteID
)
## make data frame wide for clearer representation of data
PFTC_df <- reshape(PFTC_df, direction = "wide", 
                   idvar = "Observation", timevar = "Traits")
colnames(PFTC_df) <- gsub(pattern = "Values.", replacement = "", x = colnames(PFTC_df))
PFTC_df$Species <- Raw_df$taxon[match(PFTC_df$Observation, Raw_df$id)]
PFTC_df$SiteID <- Raw_df$SiteID[match(PFTC_df$Observation, Raw_df$id)]
## Remove species not recognized by taxonomy here
PFTC_df <- PFTC_df[PFTC_df$Species %in% colnames(Phylo_ls[[1]]), ]


### TREATMENT IDENIFICATION  ####
## split SiteID so we can identify observations by Site and Treatment separately 
Treatments_df <- data.frame(Treatment = sapply(str_split(Index_df$SiteID, "_"), "[[", 2),
                            Plot = sapply(str_split(Index_df$SiteID, "_"), "[[", 1)
)
Treatments_vec <- c("ALL", as.character(unique(Treatments_df$Treatment)), as.character(unique(Treatments_df$Plot))) # find treatmend and site identifiers
Treatments_vec <- Treatments_vec[Treatments_vec != "NA"] # remove NA identifier

### TREATMENT LOOP  ####
## loop over all treatments/sites and run network methodology
for(Treatment_Iter in Treatments_vec){ 
  #### Data Preparation and Subsetting ####
  ## Create iteration versions of the base data
  Run_Index_df <- Index_df
  Run_Neigh_df <- Neigh_ls
  Run_Envir_df <- Envir_ls
  Run_PFTC_df <- PFTC_df
  ## subset if needed according to iteration
  if(Treatment_Iter != "ALL"){
    ## Identify Treatments in Envir
    NeighTreatment_df <- data.frame(Treatment = sapply(str_split(Run_Envir_df$Mean$SiteID, "_"), "[[", 2),
                                    Site = sapply(str_split(Run_Envir_df$Mean$SiteID, "_"), "[[", 1)
    )
    Run_Envir_df$Mean <- Run_Envir_df$Mean[with(NeighTreatment_df, Site == Treatment_Iter | Treatment == Treatment_Iter), ]
    Run_Envir_df$SD <- Run_Envir_df$SD[with(NeighTreatment_df, Site == Treatment_Iter | Treatment == Treatment_Iter), ]
    
    ## Identify Treatments in Neigh
    NeighTreatment_df <- data.frame(Treatment = sapply(str_split(Run_Neigh_df$Mean$SiteID, "_"), "[[", 2),
                                    Site = sapply(str_split(Run_Neigh_df$Mean$SiteID, "_"), "[[", 1)
    )
    Run_Neigh_df$Mean <- Run_Neigh_df$Mean[with(NeighTreatment_df, Site == Treatment_Iter | Treatment == Treatment_Iter), ]
    Run_Neigh_df$SD <- Run_Neigh_df$SD[with(NeighTreatment_df, Site == Treatment_Iter | Treatment == Treatment_Iter), ]
    
    ## Identify Treatments in Index
    NeighTreatment_df <- data.frame(Treatment = sapply(str_split(Run_Index_df$SiteID, "_"), "[[", 2),
                                    Site = sapply(str_split(Run_Index_df$SiteID, "_"), "[[", 1)
    )
    Run_Index_df <- Run_Index_df[with(NeighTreatment_df, Site == Treatment_Iter | Treatment == Treatment_Iter), ]
    
    ## Identify Treatments in Traits
    PFTCTreatment_df <- data.frame(Treatment = sapply(str_split(Run_PFTC_df$SiteID, "_"), "[[", 2),
                                   Site = sapply(str_split(Run_PFTC_df$SiteID, "_"), "[[", 1)
    )
    Run_PFTC_df <- Run_PFTC_df[with(PFTCTreatment_df, Site == Treatment_Iter | Treatment == Treatment_Iter), -1]
  }
  
  ### Data Sanity Checks ####
  #### Assuring only sites with information are retained
  Sites_vec <- intersect(unique(Run_Index_df$SiteID), unique(Run_Neigh_df$Mean$SiteID))
  Sites_vec <- intersect(Sites_vec, unique(Run_Envir_df$Mean$SiteID))
  Run_Index_df <- Run_Index_df[Run_Index_df$SiteID %in% Sites_vec, ]
  Run_Neigh_df$Mean <- Run_Neigh_df$Mean[Run_Neigh_df$Mean$SiteID %in% Sites_vec, ]
  Run_Neigh_df$SD <- Run_Neigh_df$SD[Run_Neigh_df$SD$SiteID %in% Sites_vec, ]
  Run_Envir_df$Mean <- Run_Envir_df$Mean[Run_Envir_df$Mean$SiteID %in% Sites_vec, ]
  Run_Envir_df$SD <- Run_Envir_df$SD[Run_Envir_df$SD$SiteID %in% Sites_vec, ]
  Run_PFTC_df <- Run_PFTC_df[Run_PFTC_df$SiteID %in% Sites_vec, ]
  #### Assuring only observed species are retained
  Run_Index_df <- Run_Index_df[Run_Index_df$Outcome != 0, ]
  Run_Neigh_df$Mean <- Run_Neigh_df$Mean[, c(1,which(colSums(Run_Neigh_df$Mean[,-1], na.rm = TRUE) != 0)+1)] # retain only neighbours which are observed at least once at each site/treatment
  Run_Neigh_df$SD <- Run_Neigh_df$SD[, c(1,which(colSums(Run_Neigh_df$Mean[,-1], na.rm = TRUE) != 0)+1)]
  OmitCols <- -(which(colnames(Run_Neigh_df$Mean)[-1] %nin% Run_Index_df$Species)+1)
  if(length(OmitCols) != 0){
    Run_Neigh_df$Mean <- Run_Neigh_df$Mean[, OmitCols]
    Run_Neigh_df$SD <- Run_Neigh_df$SD[, OmitCols]
  }
  #### Limiting Phylogeny to observed species
  Phylo_specs <- which(colnames(Phylo_ls[[1]]) %in% unique(Run_Index_df$Species))
  Phylo_ls$Mean <- Phylo_ls$Mean[Phylo_specs, Phylo_specs]
  Phylo_ls$SD <- Phylo_ls$SD[Phylo_specs, Phylo_specs]
  
  ### Trait Distances ####
  Traits_ls <- list(Mean = aggregate(.~Species, data = Run_PFTC_df[,-ncol(Run_PFTC_df)], mean, na.rm = TRUE),
                    SD = aggregate(.~Species, data = Run_PFTC_df[,-ncol(Run_PFTC_df)], sd, na.rm = TRUE))
  rownames(Traits_ls$Mean) <- Traits_ls$Mean$Species
  Traits_ls$Mean <- Traits_ls$Mean[,-1]
  TraitDistance <- gawdis(x = Traits_ls$Mean
                          # , w.type ="optimized"
  )
  Traits_ls$Mean <- as.matrix(TraitDistance)
  
  ### Model Data List ####
  Stan_list <- FUN.StanList(
    Outcome = "Outcome", # name of the outcome variable column in Index_df
    Index_df = Run_Index_df, # data frame containing columns for ID, Fitness, Site, and Species
    Neigh_df = Run_Neigh_df$Mean, # data frame containing column for site and columns for all species; cells containing predictor variable of each species at each site
    Neigh_UC = Run_Neigh_df$SD,
    Envir_df = Run_Envir_df$Mean, # data frame containing column for site and columns for climate variables; cells containing climate variable values at each site
    Envir_UC = Run_Envir_df$SD,
    Phylo_df = Phylo_ls$Mean, # phylogenetic distances
    Phylo_UC = Phylo_ls$SD,
    Traits_df = Traits_ls$Mean # trait similarities
  )
  
  #### Preferences
  rstan_options(auto_write = TRUE)
  rstan_options(javascript = FALSE)
  options(mc.cores = nChains) 
  
  ## CHECKING DATA
  FUN.DataDims(data = Stan_list)
  if(Stan_list$N < 100){ # a pretty arbitrary cut-off for now!!!
    message(paste("Not enough samples available for treatment/site", Treatment_Iter))
    next()
  }
  
  ### Directory Creation ####
  Dir.PlotNets.PFTC.Iter <- file.path(Dir.PlotNets.PFTC, Treatment_Iter)
  if(file.exists(file.path(Dir.PlotNets.PFTC, paste0(Treatment_Iter, ".jpeg")))){
    print(paste("PFTC models already run for Treatment/Plot, ", Treatment_Iter,". You can find them here:", Dir.PlotNets.PFTC.Iter))
    next()
  }
  dir.create(Dir.PlotNets.PFTC.Iter)  
  
  ### Running Model ####
  Stan_model <- stan(file = 'StanModel_Development.stan',
                     data =  Stan_list,               # named list of data
                     chains = nChains,
                     warmup = nWarmup,          # number of warmup iterations per chain
                     iter = nSamples,            # total number of iterations per chain
                     refresh = 100,         # show progress every 'refresh' iterations
                     control = list(max_treedepth = 10)
  )
  save(Stan_model, file = file.path(Dir.PlotNets.PFTC.Iter, "Model.RData"))
  # load(file.path(Dir.PlotNets.PFTC.Iter, "Model.RData")) # for testing
  
  ### Model Diagnostics ####
  # Get the full posteriors 
  joint.post.draws <- extract.samples(Stan_model)
  # Select parameters of interest
  param.vec <- c('a', 'beta_ij', 'effect', 'response', 're', 'inter_mat', 'mu', 'sigma_alph')
  # Draw 1000 samples from the 80% posterior interval for each parameter of interest
  p.samples <- list()
  p.samples <- sapply(param.vec[param.vec != 'sigma_alph' & param.vec != 'inter_mat'], function(p) {
    p.samples[[p]] <- apply(joint.post.draws[[p]], 2, function(x){
      sample(x[x > quantile(x, 0.1) & x < quantile(x, 0.9)], size = 1000)
    })  # this only works for parameters which are vectors
  })
  # there is only one sigma_alph parameter so we must sample differently:
  p.samples[['sigma_alph']] <- sample(joint.post.draws$sigma[
    joint.post.draws$sigma > quantile(joint.post.draws$sigma, 0.1) & 
      joint.post.draws$sigma < quantile(joint.post.draws$sigma, 0.9)], size = 1000)
  # WARNING: in the STAN model, parameter 'a' lies within a logarithmic, and must thus be logarithmitised to return estimates of intrinsic performance
  intrinsic.perf <- log(p.samples$a)
  colnames(intrinsic.perf) <- levels(Stan_list$Index$Species)
  inter_mat <- return_inter_array(joint.post.draws, 
                                  response = p.samples$response,
                                  effect = p.samples$effect,
                                  levels(Stan_list$Index$Species),
                                  colnames(Run_Neigh_df)[-1])
  # inter_mat is now a 3 dimensional array, where rows = focals, columns = neighbours and 3rd dim = samples from the posterior; inter_mat[ , , 1] should return a matrix consisting of one sample for every interaction 
  param.vec <- c('a', 'beta_ij', 'effect', 'response', 're', 'inter_mat', 'mu', 'sigma')
  try(stan_model_check(fit = Stan_model,
                       results_folder = Dir.PlotNets.PFTC.Iter,
                       params = param.vec))
  
  ### Plotting ####
  FUN.PlotNetUncert(Model = inter_mat, Dir = Dir.PlotNets.PFTC, Name = Treatment_Iter)
}

## NETASSOC ----------------------------------------------------------------

## COOCCUR -----------------------------------------------------------------
