#' ######################################################################## #
#' PROJECT: [Data Simplification for Network Inference across Scales] 
#' CONTENTS: 
#'  - Generate PFTC species-association/-interaction networks
#'  DEPENDENCIES:
#'  - 0 - Preamble.R
#'  - X - Functions_Bayes.R
#'  - X - Functions_Data.R
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
source("X - Functions_Bayes.R")
source("X - Functions_Data.R")

## Bayes Settings ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
nSamples <- 9000
nWarmup <- 1000
nChains <- 4
thin <- 1

# DATA =====================================================================
message("############ PREPARING DATA")
## PHYLOGENY ---------------------------------------------------------------
#### Bootstrapped Phylogenetic Distance Creation 
if(!file.exists(file.path(Dir.PFTC, "Phylogeny.RData"))){ # check if phylogenetic distance has already been established
  message("Generating phylogeny")
  download.file(url = "https://osf.io/hjpwt/download", 
                destfile = file.path(Dir.PFTC, "PFTC3-Puna-PFTC5_Peru_2018-2020_LeafTraits_clean.csv"), mode = "wb")
  Raw_df <- read.csv(file.path(Dir.PFTC, "PFTC3-Puna-PFTC5_Peru_2018-2020_LeafTraits_clean.csv")) # load pftc data
  
  Phylo_ls <- FUN.PhyloDist(SpeciesNames = Raw_df$taxon) # get phylogenetic distance with mixed species-/genus-pool
  save(Phylo_ls, file = file.path(Dir.PFTC, "Phylogeny.RData")) # save phylo dist and sd
  unlink(file.path(Dir.PFTC, "PFTC3-Puna-PFTC5_Peru_2018-2020_LeafTraits_clean.csv"))
}else{ # if phylo has already been established prior
  message("Phylogeny already generated")
  load(file.path(Dir.PFTC, "Phylogeny.RData")) # load object containing phylo distance
}
Phylo_Specs <- gsub(pattern = "_", replacement = " ", x =  Phylo_ls$Avg_Phylo$tip.label)
Phylo_ls$Avg_Phylo$tip.label <- Phylo_Specs

## META & CLIMATE DATA ----------------------------------------------------
ECV_vec <- c("2m_temperature", "volumetric_soil_water_layer_1", "total_precipitation", "potential_evaporation")
if(!file.exists(file.path(Dir.PFTC, "Metadata_df.csv"))){ # bioclimatic data not established yet
  message("Obtaining covariate data")
  download.file(url = "https://osf.io/uk85w/download", 
                destfile = file.path(Dir.PFTC, "PU.10_PFTC3.10_2020_Peru_Coordinates.xlsx"), mode = "wb")
  Metadata_df <- as.data.frame(readxl::read_xlsx(file.path(Dir.PFTC, "PU.10_PFTC3.10_2020_Peru_Coordinates.xlsx"))) # read metadata file, available here https://osf.io/uk85w/
  Metadata_df <- Metadata_df[-which(is.na(Metadata_df$PlotID)), ]
  Metadata_df$SiteID <- with(Metadata_df, paste(Site, Treatment, PlotID, sep = "_")) # create ID column which combines site, treatment, and plotid
  Metadata_df$Year <- 2018 # set observation years to 2018 for all sites (that was the field season most of the data was collected throughout)
  colnames(Metadata_df)[5] <- "Lat"
  colnames(Metadata_df)[6] <- "Lon"
  Shp <- KrigR::buffer_Points(Points = Metadata_df, Buffer = 0.1, ID = "SiteID")
  for(Clim_Iter in 1:length(ECV_vec)){
    ## climate data download
    FUN.CLIM(ECV = Clim_Iter, Shp = Shp, Dir = Dir.PFTC)
    Era5Data_ras <- stack(file.path(Dir.PFTC, paste0(ECV_vec[Clim_Iter], ".nc")))
    Era5Uncert_ras <- stack(file.path(Dir.PFTC, paste0("UC",  ECV_vec[Clim_Iter], ".nc")))
    ## Extract Data to Locations
    Extract_df <- cbind(Metadata_df$Lon, Metadata_df$Lat)
    colnames(Extract_df) <- c("x", "y")
    Extract_sp <- SpatialPoints(Extract_df)
    Data_mat <- raster::extract(Era5Data_ras, Extract_sp)
    Uncert_mat <- raster::extract(Era5Uncert_ras, Extract_sp)
    ## Save Mean Values to Original Data Source
    Metadata_df$XYZ <- NA
    colnames(Metadata_df)[ncol(Metadata_df)] <- ECV_vec[Clim_Iter]
    Metadata_df[, ncol(Metadata_df)] <- rowMeans(Data_mat)
    Metadata_df$XYZ <- NA
    colnames(Metadata_df)[ncol(Metadata_df)] <- paste0(ECV_vec[Clim_Iter], "_SD")
    Metadata_df[, ncol(Metadata_df)] <- apply(Data_mat, 1, sd)
    Metadata_df$XYZ <- NA
    colnames(Metadata_df)[ncol(Metadata_df)] <- paste0(ECV_vec[Clim_Iter], "_UC")
    Metadata_df[, ncol(Metadata_df)] <- rowMeans(Uncert_mat)
    ## writing result
    write.csv(Metadata_df, file.path(Dir.PFTC, "Metadata_df.csv")) 
  }
  unlink(file.path(Dir.PFTC, "PU.10_PFTC3.10_2020_Peru_Coordinates.xlsx"))
  unlink(file.path(Dir.Data, list.files(Dir.Data, pattern = ".nc")))
}else{ # bioclimatic data has been established
  message("Covariate data already retrieved")
  Metadata_df <- read.csv(file.path(Dir.PFTC, "Metadata_df.csv"))[-1] # load the metadata with attached bioclimatic data
}
# set this because 2m_temperature is stored as x2m_temperature in Metadata_df
ECV_vec[ECV_vec == "2m_temperature"] <- "X2m_temperature"
Metadata_df <- Sort.DF(Metadata_df, "SiteID")

## MODEL DATA FRAMES -------------------------------------------------------
if(!file.exists(file.path(Dir.PFTC, "ModelFrames.RData"))){
  message("Readying model data frames")
  download.file(url = "https://osf.io/hjpwt/download", 
                destfile = file.path(Dir.PFTC, "PFTC3-Puna-PFTC5_Peru_2018-2020_LeafTraits_clean.csv"), mode = "wb")
  Raw_df <- read.csv(file.path(Dir.PFTC, "PFTC3-Puna-PFTC5_Peru_2018-2020_LeafTraits_clean.csv")) # load pftc data
  Raw_df <- Raw_df[-which(is.na(Raw_df$plot_id)), ]
  Raw_df$SiteID <- with(Raw_df, paste(site, treatment, plot_id, sep = "_")) # create SiteID index
  Raw_df <- Raw_df[-which(Raw_df$SiteID %nin% Metadata_df$SiteID), ]
  Raw_df <- Raw_df[Raw_df$taxon %in% Phylo_Specs, ] # limiting to phylogeny-recognised species
  ### FITNESS ####
  Weight_df <- Raw_df[Raw_df$trait == "dry_mass_g", c("SiteID", "taxon", "value")] # only dry biomass rows and select only relevant columns for speedier data handling 
  Weight_cast <- reshape2::dcast(Weight_df, SiteID ~ taxon, value.var = 'value', fun.aggregate = mean)
  Weight_cast[Weight_cast == -Inf] <- 0
  ## ABUNDANCE ####
  Raw_df <- read.csv(file.path(Dir.PFTC, "PFTC3-Puna-PFTC5_Peru_2018-2020_LeafTraits_clean.csv")) # load pftc data
  Raw_df <- Raw_df[-which(is.na(Raw_df$plot_id)), ]
  Raw_df$SiteID <- with(Raw_df, paste(site, treatment, plot_id, sep = "_")) # create SiteID index
  Raw_df <- Raw_df[-which(Raw_df$SiteID %nin% Metadata_df$SiteID), ]
  Raw_df <- Raw_df[Raw_df$taxon %in% Phylo_Specs, ] # limitting to phylogeny-recognised species
  Abund_df <- Raw_df[,c("SiteID", "taxon", "individual_nr")] # limit data to necessary parts
  Abund_cast <- reshape2::dcast(Abund_df, SiteID ~ taxon, value.var = 'individual_nr', fun.aggregate = max)
  Abund_cast[Abund_cast == -Inf] <- 0
  Abund_cast <- Abund_cast[ , colnames(Abund_cast) %in% colnames(Weight_cast)]
  ## LIST MERGING ####
  ModelFrames_ls <- list(
    Fitness = Weight_df,
    FitCom = Weight_cast,
    Community = Abund_cast
  )
  save(ModelFrames_ls, file = file.path(Dir.PFTC, "ModelFrames.RData"))
  unlink(file.path(Dir.PFTC, "PFTC3-Puna-PFTC5_Peru_2018-2020_LeafTraits_clean.csv"))
}else{
  message("Model data frames already compiled")
  load(file.path(Dir.PFTC, "ModelFrames.RData"))
}

## FUNCTIONAL TRAITS -------------------------------------------------------
## download PFTC data if needed
if(!file.exists(file.path(Dir.PFTC, "Traits_df.csv"))){
  message("Obtaining trait data")
  download.file(url = "https://osf.io/hjpwt/download", 
                destfile = file.path(Dir.PFTC, "PFTC3-Puna-PFTC5_Peru_2018-2020_LeafTraits_clean.csv"), mode = "wb")
  Raw_df <- read.csv(file.path(Dir.PFTC, "PFTC3-Puna-PFTC5_Peru_2018-2020_LeafTraits_clean.csv")) #
  Raw_df <- Raw_df[-which(is.na(Raw_df$plot_id)), ]
  Raw_df$SiteID <- with(Raw_df, paste(site, treatment, plot_id, sep = "_")) # create SiteID index
  Raw_df <- Raw_df[-which(Raw_df$SiteID %nin% Metadata_df$SiteID), ]
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
  Traits_df <- Traits_df[Traits_df$Species %in% Phylo_Specs, ]
  write.csv(Traits_df, file = file.path(Dir.PFTC, "Traits_df.csv"))
  unlink(file.path(Dir.PFTC, "PFTC3-Puna-PFTC5_Peru_2018-2020_LeafTraits_clean.csv"))
}else{
  message("Trait data already obtained")
  Traits_df <- read.csv(file.path(Dir.PFTC, "Traits_df.csv"))[ , -1]
}

## SITEID MATCHING ---------------------------------------------------------
# we don't have measures at some of the sites. Here, we remove these sites from the metadata frame
Metadata_df <- Metadata_df[Metadata_df$SiteID %in% ModelFrames_ls$FitCom$SiteID, ]
if(sum(ModelFrames_ls$FitCom$SiteID != Metadata_df$SiteID)){stop("ModelFrames and Metadata Frame are not in the same order")}
Phylo_ls$Avg_Phylo <- drop.tip(phy = Phylo_ls$Avg_Phylo, tip = Phylo_ls$Avg_Phylo$tip.label[Phylo_ls$Avg_Phylo$tip.label %nin% colnames(ModelFrames_ls$FitCom)[-1]])

## TREATMENT IDENTIFICATION ------------------------------------------------
## split SiteID so we can identify observations by Site and Treatment separately 
Treatments_vec <- unique(c(Metadata_df$Site, Metadata_df$Treatment))
Treatments_vec <- c("ALL", Treatments_vec)

# ANALYSIS =================================================================
## IF-REM ------------------------------------------------------------------
message("############ STARTING IF-REM ANALYSES")
Dir.IFREM <- file.path(DirEx.PFTC, "IF_REM")
if(!dir.exists(Dir.IFREM)){dir.create(Dir.IFREM)}

## combine SiteID, focal fitness, and neighbour counts
Index_df <- cbind(ModelFrames_ls$Fitness, 
                  ModelFrames_ls$Community[match(ModelFrames_ls$Fitness$SiteID,
                                                 ModelFrames_ls$Community$SiteID),  -1])

for(Treatment_Iter in Treatments_vec){ # HMSC treatment loop
  message(paste("### Treatment:", Treatment_Iter))
  Dir.TreatmentIter <- file.path(Dir.IFREM, Treatment_Iter)
  if(!dir.exists(Dir.TreatmentIter)){dir.create(Dir.TreatmentIter)}
  
  ### DATA SUBSETTING ####
  Index_Iter <- Index_df
  if(Treatment_Iter != "ALL"){
    Treat_df <- data.frame(Treatment = sapply(str_split(Index_Iter$SiteID, "_"), "[[", 2),
                           Site = sapply(str_split(Index_Iter$SiteID, "_"), "[[", 1)
    )
    Index_Iter <- Index_Iter[with(Treat_df, Site == Treatment_Iter | Treatment == Treatment_Iter), ]
  }
  Index_Iter <- Index_Iter[, c(1:3, which(colSums(Index_Iter[, -1:-3]) != 0) + 3)] # ensuring only present species are retained
  
  ### DATA PREPRATION ####
  StanList_Iter <- FUN.StanList(Fitness = "value", data = Index_Iter)
  
  ### DATA CHECKS ####
  FUN.DataDims(data = StanList_Iter)
  
  #### Preferences
  rstan_options(auto_write = TRUE)
  rstan_options(javascript = FALSE)
  options(mc.cores = nChains) 
  
  ### Model Execution ----
  Stan_model <- stan(file = 'joint_model.stan',
                     data =  StanList_Iter,
                     chains = 1,
                     warmup = nWarmup*nChains/2,
                     iter = nSamples*nChains/2,
                     refresh = 100,
                     control = list(max_treedepth = 10)
  )
  save(Stan_model, file = file.path(Dir.TreatmentIter, "Model.RData"))
  
  ### Model Diagnostics ----
  # Get the full posteriors 
  joint.post.draws <- extract.samples(Stan_model)
  # Select parameters of interest
  param.vec <- c('beta_i0', 'beta_ij', 'effect', 'response', 're', 'inter_mat', 'mu')
  # Draw 1000 samples from the 80% posterior interval for each parameter of interest
  p.samples <- list()
  p.samples <- sapply(param.vec[param.vec != 'inter_mat'], function(p) {
    p.samples[[p]] <- apply(joint.post.draws[[p]], 2, function(x){
      sample(x[x > quantile(x, 0.1) & x < quantile(x, 0.9)], size = nSamples)
    })  # this only works for parameters which are vectors
  })
  # there is only one sigma_alph parameter so we must sample differently:
  p.samples[['sigma_alph']] <- sample(joint.post.draws$sigma[
    joint.post.draws$sigma > quantile(joint.post.draws$sigma, 0.1) & 
      joint.post.draws$sigma < quantile(joint.post.draws$sigma, 0.9)], size = nSamples)
  # WARNING: in the STAN model, parameter 'a' lies within a logarithmic, and must thus be logarithmitised to return estimates of intrinsic performance
  intrinsic.perf <- log(p.samples$beta_i0)
  colnames(intrinsic.perf) <- levels(factor(Index_Iter$taxon))
  inter_mat <- return_inter_array(joint.post.draws = joint.post.draws, 
                                  response = p.samples$response,
                                  effect = p.samples$effect,
                                  focalID = levels(factor(Index_Iter$taxon)),
                                  neighbourID = colnames(Index_Iter[, -1:-3]))
  # inter_mat is now a 3 dimensional array, where rows = focals, columns = neighbours and 3rd dim = samples from the posterior; inter_mat[ , , 1] should return a matrix consisting of one sample for every interaction 
  try(stan_model_check(fit = Stan_model,
                       results_folder = Dir.TreatmentIter,
                       params = param.vec))
  
  ### Interaction/Association Matrix ----
  Interaction_mean <- apply(inter_mat, c(1, 2), mean) # will return the mean estimate for every interaction (NB: this is the mean of the 80% posterior interval, so will be slightly different to the mean value returned from summary(fit), which is calculated from the full posterior distribution)  
  Interaction_mean <- Interaction_mean*-1 # need to switch sign of results
  diag(Interaction_mean) <- NA
  Interaction_hpdi <- apply(inter_mat, c(1, 2), HPDI, prob = 0.89)
  Interaction_min <- -Interaction_hpdi[1,,]
  diag(Interaction_min) <- NA
  Interaction_max <- -Interaction_hpdi[2,,]
  diag(Interaction_max) <- NA
  Interactions_igraph <- data.frame(Actor = rep(dimnames(Interaction_mean)$neighbour, length(dimnames(Interaction_mean)$species)),
                                    Subject = rep(dimnames(Interaction_mean)$species, each = length(dimnames(Interaction_mean)$neighbour)),
                                    Inter_mean = as.vector(t(Interaction_mean)),
                                    Inter_min = as.vector(t(Interaction_min)),
                                    Inter_max = as.vector(t(Interaction_max))
  )
  Interactions_IFREM <- Interactions_igraph[order(abs(Interactions_igraph$Inter_mean), decreasing = TRUE), ]
  Interactions_IFREM <- na.omit(Interactions_IFREM)
  save(Interactions_IFREM, file = file.path(Dir.TreatmentIter, "Interac.RData"))
  
  # ### Plotting ####
  # FUN.PlotNetUncert(Model = inter_mat, Dir = Dir.PlotNets.PFTC, Name = Treatment_Iter)
} 

## HMSC --------------------------------------------------------------------
message("############ STARTING HMSC ANALYSES")
Dir.HMSC <- file.path(DirEx.PFTC, "HMSC")
if(!dir.exists(Dir.HMSC)){dir.create(Dir.HMSC)}

for(Treatment_Iter in Treatments_vec){ # HMSC treatment loop
  message(paste("### Treatment:", Treatment_Iter))
  Dir.TreatmentIter <- file.path(Dir.HMSC, Treatment_Iter)
  if(!dir.exists(Dir.TreatmentIter)){dir.create(Dir.TreatmentIter)}
  
  ### DATA SUBSETTING ####
  ModelFrames_Iter <- ModelFrames_ls
  Metadata_Iter <- Metadata_df
  Phylo_Iter <- Phylo_ls$Avg_Phylo
  Traits_Iter <- Traits_df
  
  if(Treatment_Iter != "ALL"){
    ## Fitness List
    Treat_df <- data.frame(Treatment = sapply(str_split(ModelFrames_Iter$Fitness$SiteID, "_"), "[[", 2),
                           Site = sapply(str_split(ModelFrames_Iter$Fitness$SiteID, "_"), "[[", 1)
    )
    ModelFrames_Iter$Fitness <- ModelFrames_Iter$Fitness[with(Treat_df, Site == Treatment_Iter | Treatment == Treatment_Iter), ]
    ## Community Matrices
    Treat_df <- data.frame(Treatment = sapply(str_split(ModelFrames_Iter$FitCom$SiteID, "_"), "[[", 2),
                           Site = sapply(str_split(ModelFrames_Iter$FitCom$SiteID, "_"), "[[", 1)
    )
    ModelFrames_Iter$FitCom <- ModelFrames_Iter$FitCom[with(Treat_df, Site == Treatment_Iter | Treatment == Treatment_Iter), ]
    ModelFrames_Iter$Community <- ModelFrames_Iter$Community[with(Treat_df, Site == Treatment_Iter | Treatment == Treatment_Iter), ]
    ## Metadata 
    Treat_df <- data.frame(Treatment = sapply(str_split(Metadata_Iter$SiteID, "_"), "[[", 2),
                           Site = sapply(str_split(Metadata_Iter$SiteID, "_"), "[[", 1)
    )
    Metadata_Iter <- Metadata_Iter[with(Treat_df, Site == Treatment_Iter | Treatment == Treatment_Iter), ]
  }
  Traits_Iter <- aggregate(. ~ Species, data = Traits_Iter[,c(-1,-10)], FUN = mean)
  rownames(Traits_Iter) <- Traits_Iter$Species
  
  ### DATA PREPRATION ####
  S <- Metadata_Iter[ , c("SiteID", "Site", "Treatment", "PlotID", "Lat", "Lon")] # S: study design, including units of study and their possible coordinates, If you don't have variables that define the study design, indicate this by S=NULL
  X <- Metadata_Iter[ , c("SiteID", "Site", "Treatment", "PlotID", "Elevation", ECV_vec)] # X: covariates to be used as predictors, If you don't have covariate data, indicate this by X=NULL
  Y_BM <- ModelFrames_Iter$FitCom[ , -1] # Y: species data
  Y_AB <- ModelFrames_Iter$Community[ , -1] # Y: species data
  P <- Phylo_Iter # P: phylogenetic information given by taxonomical levels, e.g. order, family, genus, species; If TP does not have phylogenetic data (because you don't have such data at all, or because, it is given in tree-format, like is the case in this example), indicate this with P=NULL 
  Tr <- Traits_Iter # Tr: species traits (note that T is a reserved word in R and that's why we use Tr); If you don't have trait data, indicate this by Tr=NULL. 
  
  ### DATA CHECKS ####
  if(all(dim(Y_AB) == dim(Y_BM))){print("Community matrices are the same dimensions")}else{stop("Community matrices have unequal dimensions")}
  if(is.numeric(as.matrix(Y_AB)) || is.logical(as.matrix(Y_AB)) && is.finite(sum(Y_AB, na.rm=TRUE))){print("Species data looks ok")
  }else{stop("Species data should be numeric and have finite values")}
  if(any(is.na(S))){stop("study design has NA values - not allowed for")
  }else{print("study design looks ok")}
  if(any(is.na(X))){stop("Covariate data has NA values - not allowed for")
  }else{print("Covariate data looks ok")}
  if(all(colnames(Y_AB) %in% rownames(Tr))){print("species names in TP and SXY match")
  }else{stop("species names in TP and SXY do not match")}
  if(any(is.na(Tr))){stop("Tr has NA values - not allowed for")
  }else{print("Tr looks ok")}
  if(any(is.na(P))){stop("P has NA values - not allowed for")
  }else{print("P looks ok")}
  if(all(sort(P$tip.label) == sort(colnames(Y_AB)))){print("species names in P and SXY match")}else{stop("species names in P and SXY do not match")}
  
  ### Model Specification ----
  ## removing rare species
  Ypa <- 1*(Y_AB>0) # identify presences (1) and absences (0) in species data
  # rarespecies defined as less than three occurrences at all sites
  # rarespecies <- which(colSums(Ypa)<3)
  # message(paste("Removing", length(rarespecies), "rare species."))
  # Y_AB <- Y_AB[,-rarespecies]
  # Y_BM <- Y_BM[,-rarespecies]
  # Ypa <- Ypa[,-rarespecies]
  # Tr <- droplevels(Tr[-rarespecies,])
  # P <- drop.tip(P, rarespecies) 
  ## recoding treatments
  TreatN <- vector(mode="numeric", length=nrow(X))
  TreatN[X$Treatment=="C"] <- 0
  TreatN[X$Treatment=="B"] <- 1
  TreatN[X$Treatment=="BB"] <- 2
  TreatN[X$Treatment=="NB"] <- 3
  X$TreatN = factor(TreatN, ordered = TRUE)
  X <- X[, c("TreatN", "Elevation", ECV_vec)]
  ## Model Formulae
  XFormula <- as.formula(paste("~ TreatN + Elevation", paste(ECV_vec, collapse = " + "), sep = " + "))
  TrFormula <- ~log(plant_height_cm) + leaf_area_cm2 + leaf_thickness_mm
  ## StudyDesign
  unique_plot <- paste(S$Site, S$Plot, sep="_") # Spatial coordinates are included, but are clustered at the five sites, and by treatment within sites. Other S data are Site (n=5) and within-site and treatment plots (1-5). Elevation is already in the X table.
  studyDesign <- data.frame(site = as.factor(S$Site), unique_plot = as.factor(unique_plot))
  St <- studyDesign$site
  rL.site <- HmscRandomLevel(units = levels(St))
  Pl <- studyDesign$unique_plot
  rL.plot <- HmscRandomLevel(units = levels(Pl))
  ## Model Objects
  Ypa <- 1*(Y_AB>0)
  Yabu <- Y_AB
  Yabu[Y_AB==0] <- NA
  Yabu <- log(Yabu)
  Ybiom <- Y_BM
  Ybiom[Y_BM==0] <- NA
  Ybiom <- log(Ybiom)
  ## Models
  m1 <- Hmsc(Y=Ypa, XData = X,  XFormula = XFormula,
             TrData = Tr,TrFormula = TrFormula,
             distr="probit",
             studyDesign=studyDesign,
             ranLevels={list("site" = rL.site, "unique_plot" = rL.plot)})
  m2 <- Hmsc(Y=Yabu, YScale = TRUE,
             XData = X,  XFormula = XFormula,
             TrData = Tr,TrFormula = TrFormula,
             distr="normal",
             studyDesign=studyDesign,
             ranLevels={list("site" = rL.site, "unique_plot" = rL.plot)})
  m3 <- Hmsc(Y=Ybiom, YScale = TRUE,
             XData = X,  XFormula = XFormula,
             TrData = Tr,TrFormula = TrFormula,
             distr="normal",
             studyDesign=studyDesign,
             ranLevels={list("site" = rL.site, "unique_plot" = rL.plot)})
  models <- list(m1,m2, m3)
  modelnames <- c("presence_absence","abundance", "biomass")
  
  ### Modelling ----
  samples_list <- nSamples
  message(paste0("thin = ",as.character(thin),"; samples = ",as.character(nSamples)))
  for(HMSC_Iter in 1:length(modelnames)){ # HMSC model loop
    #### Model Execution ----
    hmsc_modelname <- modelnames[HMSC_Iter]
    hmsc_model <- models[[HMSC_Iter]]
    message(paste0("model = ",hmsc_modelname))
    filename <- file.path(Dir.TreatmentIter, paste0(hmsc_modelname, ".RData"))
    if(!file.exists(filename)){
      hmsc_model <- sampleMcmc(hmsc_model, samples = nSamples, thin = thin,
                               transient = nWarmup,
                               nChains = nChains,
                               nParallel = nChains) 
      save(hmsc_model,hmsc_modelname,file=filename)
    }else{
      load(filename)
    }
    
    ### Model Evaluation ----
    message("Evaluation")
    vals <- HMSC.Eval(Model = hmsc_model, Dir = Dir.TreatmentIter, Name = hmsc_modelname, thin = thin, nSamples = nSamples, nChains = nChains)
    ### Interaction/Association Matrix ----
    Interaction_mean <- vals$`Posterior mean`[,-1]
    Interaction_ProbPos <- vals$`Pr(x>0)`[,-1]
    Interaction_ProbNeg <- vals$`Pr(x<0)`[,-1]
    Partner2 <- c()
    for(i in 1:(length(colnames(Interaction_mean))-1)){
      Partner2 <- c(Partner2, colnames(Interaction_mean)[-c(1:i)])
    }
    Interactions_igraph <- data.frame(Partner1 = rep(rownames(Interaction_mean), times = (length(colnames(Interaction_mean))-1):0),
                                      Partner2 = Partner2,
                                      Inter_mean = t(Interaction_mean)[lower.tri(t(Interaction_mean), diag = FALSE)],
                                      Inter_ProbPos = t(Interaction_ProbPos)[lower.tri(t(Interaction_ProbPos), diag = FALSE)],
                                      Inter_ProbNeg = t(Interaction_ProbNeg)[lower.tri(t(Interaction_ProbNeg), diag = FALSE)]
    )
    Interactions_HMSC <- Interactions_igraph[order(abs(Interactions_igraph$Inter_mean), decreasing = TRUE), ] 
    save(Interactions_HMSC, file = file.path(Dir.TreatmentIter, paste0(Name, "_Interac.RData")))
  } # end of HMSC model loop
} # end of HMSC treatment loop

## NETASSOC ----------------------------------------------------------------
message("############ STARTING NETASSOC ANALYSES")
Dir.NETASSOC <- file.path(DirEx.PFTC, "NETASSOC")
if(!dir.exists(Dir.NETASSOC)){dir.create(Dir.NETASSOC)}

## COOCCUR -----------------------------------------------------------------
message("############ STARTING COCCUR ANALYSES")
Dir.COOCCUR <- file.path(DirEx.PFTC, "COCCUR")
if(!dir.exists(Dir.COOCCUR)){dir.create(Dir.COOCCUR)}