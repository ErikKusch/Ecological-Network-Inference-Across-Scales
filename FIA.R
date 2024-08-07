#' ####################################################################### #
#' PROJECT: [Data Simplification for Network Inference across Scales] 
#' CONTENTS: 
#'  - Generate plot data base species-association/-interaction networks
#'  DEPENDENCIES:
#'  - 0 - Preamble.R
#'  - 0 - ShapeFiles.R
#'  - X - Functions_Data.R
#'  - X - Functions_Bayes.R
#'  - X - Functions_Plotting.R
#'  - No data dependencies
#' AUTHOR: [Erik Kusch]
#' ####################################################################### #

# PREAMBLE =================================================================
rm(list=ls())
set.seed(42)

## Sourcing ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
source("0 - Preamble.R")
source("0 - ShapeFiles.R")
source("X - Functions_Bayes.R")
source("X - Functions_Data.R")

## Bayes Settings ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
nSamples <- 1000
nWarmup <- 500
nChains <- 4
thin <- 1

# DATA =====================================================================
Shared_spec <- c("Abies concolor", "Abies magnifica", "Calocedrus decurrens","Cornus nuttallii", "Pinus lambertiana", "Pinus ponderosa", "Pseudotsuga menziesii", "Quercus kelloggii")

## CLIMATE DATA RETRIEVAL --------------------------------------------------
ECV_vec <- c("2m_temperature", "volumetric_soil_water_layer_1", "total_precipitation", "potential_evaporation")
FIA_shp <- crop(FIA_shp, extent(extent(FIA_shp)[1], -59.5, extent(FIA_shp)[3], extent(FIA_shp)[4]))
if(!file.exists(file.path(Dir.FIA, "FIABiomes_df.rds"))){
  for(Clim_Iter in 1:length(ECV_vec)){
    FUN.CLIM(ECV = Clim_Iter, Shp = FIA_shp, Dir = Dir.FIA)
  }
}

## FIA DATA RETRIEVAL ------------------------------------------------------
if(sum(file.exists(file.path(Dir.FIA, paste0("FIABiome", 9, ".RData")))) != 1){ #1:13
  FUN.FIA(states = c("AL", "AK", "AZ", "AR", "CA", "CO", "CT", "DE", "FL", "GA", "HI", "ID", "IL", "IN", "IA", "KS", "KY", "LA", "ME", "MD", "MA", "MI", "MN", "MS", "MO", "MT", "NE", "NV", "NH", "NJ", "NM", "NY", "NC", "ND", "OH", "OK", "OR", "PA", "RI", "SC", "SD", "TN", "TX", "UT", "VT", "VA", "WA", "WV", "WI", "WY"), 
          nCores = parallel::detectCores(), 
          Dir.FIA = Dir.FIA,
          Shape = FIA_shp,
          Name = "BIOMES")
}
if(sum(file.exists(file.path(Dir.FIA, paste0("FIABiome", c(1003), ".RData")))) != 1){ #1000:1001, 
  FUN.FIA(states = c("AL", "AK", "AZ", "AR", "CA", "CO", "CT", "DE", "FL", "GA", "HI", "ID", "IL", "IN", "IA", "KS", "KY", "LA", "ME", "MD", "MA", "MI", "MN", "MS", "MO", "MT", "NE", "NV", "NH", "NJ", "NM", "NY", "NC", "ND", "OH", "OK", "OR", "PA", "RI", "SC", "SD", "TN", "TX", "UT", "VT", "VA", "WA", "WV", "WI", "WY"), 
          nCores = parallel::detectCores(), 
          Dir.FIA = Dir.FIA,
          Shape = States_shp,
          Name = "STATES")
}
# FIABiomes_fs <- list.files(path = Dir.FIA, pattern = "FIABiome")

# ANALYSIS =================================================================
Treatment_Vec <- c(1003, 9 # this is the temperate confier comparison
                   #, 1000, 1001, 8 # this would be a temperate broadleaf comparison
                   )

## HMSC --------------------------------------------------------------------
message("############ STARTING HMSC ANALYSES")
Dir.HMSC <- file.path(DirEx.Region, "HMSC")
if(!dir.exists(Dir.HMSC)){dir.create(Dir.HMSC)}

for(Treatment_Iter in Treatment_Vec){ # HMSC treatment loop
  load(file.path(Dir.FIA, paste0("FIABiome",Treatment_Iter,".RData")))
  ECV_vec[1] <- "X2m_temperature"
  colnames(Metadata_df) <- gsub(colnames(Metadata_df), pattern = "2m_temperature", replacement = ECV_vec[1])
  message(paste("### Biome:", BiomeName, "(", nrow(ModelFrames_ls$Fitness), "Observations )"))
  Dir.TreatmentIter <- file.path(Dir.HMSC, Treatment_Iter)
  if(!dir.exists(Dir.TreatmentIter)){dir.create(Dir.TreatmentIter)}
  
  ### REMOVAL OF NON-FOCAL SPECIES FOR COMPARISON ###
  if(Treatment_Iter == 8){
    Shared_spec <- c("Acer rubrum", "Fagus grandifolia", "Acer saccharum", "Quercus rubra", "Fraxinus americana", "Pinus strobus", "Populus tremuloides", "Betula alleghaniensis", "Abies balsamea", "Betula papyrifera" )
    Phylo_Iter$Avg_Phylo <- drop.tip(Phylo_Iter$Avg_Phylo, Phylo_Iter$Avg_Phylo$tip.label[Phylo_Iter$Avg_Phylo$tip.label %nin% Shared_spec])
    if(is.null(Phylo_Iter$Avg_Phylo)){
      print("All shared species absent")
      next()
    }
    ModelFrames_ls$FitCom <- ModelFrames_ls$FitCom[, colnames(ModelFrames_ls$FitCom) %in% Shared_spec]
    ModelFrames_ls$Community <- ModelFrames_ls$Community[, colnames(ModelFrames_ls$Community) %in% Shared_spec]
  }
  
  # table(ModelFrames_ls$Fitness$taxon)
  
  sink(file.path(Dir.TreatmentIter, "Biome.txt"))
  print("BIOME")
  print(BiomeName)
  print("OBSERVATIONS")
  print(sum(ModelFrames_ls$Fitness$taxon %in% Phylo_Iter$Avg_Phylo$tip.label))
  print("SPECIES")
  print(length(Phylo_Iter$Avg_Phylo$tip.label))
  print("SITES")
  print(nrow(Metadata_df))
  sink()
  
  ### DATA PREPRATION ####
  Phylo_Iter <- Phylo_Iter$Avg_Phylo
  S <- Metadata_df[, c("SiteID", "YEAR")] # S: study design, including units of study and their possible coordinates, If you don't have variables that define the study design, indicate this by S=NULL
  X <- Metadata_df[,-1:-2] # X: covariates to be used as predictors, If you don't have covariate data, indicate this by X=NULL
  if(Treatment_Iter == 8){
    Y_BM <- ModelFrames_ls$FitCom
    Y_AB <- ModelFrames_ls$Community
  }else{
    Y_BM <- ModelFrames_ls$FitCom[ , -1] # Y: species data
    Y_AB <- ModelFrames_ls$Community[ , -1] # Y: species data
  }
  
  P <- Phylo_Iter # P: phylogenetic information given by taxonomical levels, e.g. order, family, genus, species; If TP does not have phylogenetic data (because you don't have such data at all, or because, it is given in tree-format, like is the case in this example), indicate this with P=NULL
  Tr <- NULL # Tr: species traits (note that T is a reserved word in R and that's why we use Tr); If you don't have trait data, indicate this by Tr=NULL.
  
  P <- drop.tip(P, P$tip.label[P$tip.label %nin% colnames(Y_AB)])
  
  ### DATA CHECKS ####
  if(all(dim(Y_AB) == dim(Y_BM))){print("Community matrices are the same dimensions")}else{stop("Community matrices have unequal dimensions")}
  if(is.numeric(as.matrix(Y_AB)) || is.logical(as.matrix(Y_AB)) && is.finite(sum(Y_AB, na.rm=TRUE))){print("Species data looks ok")
  }else{stop("Species data should be numeric and have finite values")}
  if(any(is.na(S))){stop("study design has NA values - not allowed for")
  }else{print("study design looks ok")}
  if(any(is.na(X))){stop("Covariate data has NA values - not allowed for")
  }else{print("Covariate data looks ok")}
  if(any(is.na(P))){stop("P has NA values - not allowed for")
  }else{print("P looks ok")}
  if(all(sort(P$tip.label) == sort(colnames(Y_AB)))){print("species names in P and SXY match")}else{stop("species names in P and SXY do not match")}
  
  ### Model Specification ----
  ## removing rare species
  Ypa <- 1*(Y_AB>0) # identify presences (1) and absences (0) in species data
  ## Model Formulae
  XFormula <- as.formula(paste0("~", paste(ECV_vec[1:3], collapse = " + ")))
  TrFormula <- NULL
  ## StudyDesign
  unique_plot <- paste(S$SiteID, sep="_")
  studyDesign <- data.frame(site = as.factor(S$SiteID))
  St <- studyDesign$site
  rL.site <- HmscRandomLevel(units = levels(St))
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
             ranLevels={list("site" = rL.site)})
  m2 <- Hmsc(Y=Yabu, YScale = TRUE,
             XData = X,  XFormula = XFormula,
             TrData = Tr,TrFormula = TrFormula,
             distr="poisson",
             studyDesign=studyDesign,
             ranLevels={list("site" = rL.site)})
  m3 <- Hmsc(Y=Ybiom, YScale = TRUE,
             XData = X,  XFormula = XFormula,
             TrData = Tr,TrFormula = TrFormula,
             distr="normal",
             studyDesign=studyDesign,
             ranLevels={list("site" = rL.site)})
  models <- list(m1,m2, m3)
  modelnames <- c("presence_absence","abundance", "biomass")
  
  ### Modelling ----
  samples_list <- nSamples
  message(paste0("thin = ",as.character(thin),"; samples = ",as.character(nSamples)))
  for(HMSC_Iter in 1:length(modelnames)){ # HMSC model loop
    #### Model Execution ----
    hmsc_modelname <- modelnames[HMSC_Iter]
    if(file.exists(file.path(Dir.TreatmentIter, paste0(hmsc_modelname, "_Interac.RData")))){
      message("Already sampled & evaluated")
      next()
    }
    
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
      message("Already sampled")
      load(filename)
    }
    
    ### Model Evaluation ----
    message("Evaluation")
    vals <- HMSC.Eval(hmsc_model = hmsc_model, Dir = Dir.TreatmentIter, Name = hmsc_modelname,
                      thin = thin, nSamples = nSamples, nChains = nChains)
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
    save(Interactions_HMSC, file = file.path(Dir.TreatmentIter, paste0(hmsc_modelname, "_Interac.RData")))
  } # end of HMSC model loop
} # end of HMSC treatment loop

# IF-REM ------------------------------------------------------------------
message("############ STARTING IF-REM ANALYSES")
Dir.IFREM <- file.path(DirEx.Region, "IF_REM")
if(!dir.exists(Dir.IFREM)){dir.create(Dir.IFREM)}

for(Treatment_Iter in Treatment_Vec){ # "Vermont, Maine, Yosemite, Temperate Broadleaf, Coniferous
  closeAllConnections()
  load(file.path(Dir.FIA, paste0("FIABiome",Treatment_Iter,".RData")))
  message(paste("### Biome:", BiomeName, "(", nrow(ModelFrames_ls$Fitness), "Observations )"))
  Dir.TreatmentIter <- file.path(Dir.IFREM, Treatment_Iter)
  if(!dir.exists(Dir.TreatmentIter)){dir.create(Dir.TreatmentIter)}
  if(file.exists(file.path(Dir.TreatmentIter, "Interac.RData"))){
    message("Model already compiled and evaluated")
    next()
  }
  sink(file.path(Dir.TreatmentIter, "Biome.txt"))
  print("BIOME")
  print(BiomeName)
  print("OBSERVATIONS")
  print(nrow(ModelFrames_ls$Fitness))
  print("SPECIES")
  print(nrow(Phylo_Iter$Dist_Mean))
  print("SITES")
  print(nrow(Metadata_df))
  sink()
  
  ## Simplification by genus
  # GenusFitness <- ModelFrames_ls$Fitness
  # GenusFitness$taxon <- sapply(strsplit(GenusFitness$taxon, split = " "), "[[", 1)
  # GenusFitness <- aggregate(x = GenusFitness$value,
  #           by = list(GenusFitness$SiteID, GenusFitness$taxon),
  #           FUN = "mean")
  # colnames(GenusFitness) <- c("SiteID", "taxon", "value")
  #
  # GenusCommunity <- ModelFrames_ls$Community
  # colnames(GenusCommunity) <- c("SiteID", sapply(strsplit(colnames(GenusCommunity[,-1]), split = " "), "[[", 1))
  # GenusComm <- as.data.frame(t(rowsum(t(GenusCommunity[, -1]), group = colnames(GenusCommunity)[-1], na.rm = TRUE)))
  # GenusComm$SiteID <- GenusCommunity$SiteID
  # ## combine SiteID, focal fitness, and neighbour counts
  # Index_df <- cbind(GenusFitness,
  #                   GenusComm[match(GenusFitness$SiteID, GenusComm$SiteID),  -ncol(GenusComm)])
  # # Index_df <- Index_df[, -which(colSums(Index_df[,-1:-2]) == 0)-2]
  # Index_df <- Index_df[which(Index_df$value > 0), ]
  
  # SPecies-Level analysis
  Index_df <- cbind(ModelFrames_ls$Fitness,
                    ModelFrames_ls$Community[match(ModelFrames_ls$Fitness$SiteID,
                                                   ModelFrames_ls$Community$SiteID),  -1])
  Index_df <- Index_df[, c(1:3, which(colSums(Index_df[, -1:-3]) != 0) + 3)] # ensuring only present species are retained
  Index_df <- Index_df[which(Index_df$value > 0), ]
  # Index_df <- Index_df[which(Index_df$value > 2.5), ]
  
  ### DATA PREPRATION ####
  StanList_Iter <- data_prep(perform = "value",
                             focal = "taxon",
                             plot = "SiteID",
                             df = Index_df,
                             nonNcols = 3)
    # FUN.StanList(Fitness = "value", data = Index_df)
  
  ### DATA CHECKS ####
  df <- Index_df
  neighbourID <- colnames(df)[-1:-3]
  
  N_all <- df[ , neighbourID]
  N_all <- apply(N_all, c(1,2), as.numeric)
  X_all <- cbind(model.matrix(~as.factor(df$taxon)), N_all)
  R_all <- pracma::rref(X_all)
  Z_all <- t(R_all) %*% R_all
  indep <- sapply(seq(1, dim(Z_all)[1], 1), function(k){ 
    ifelse(Z_all[k, k] == 1 & sum(Z_all[k, -k]) == 0, 1, 0)
  }) #
  all(indep == 1) # if TRUE then neighbours are linearly independent and we can continue
  if(!all(indep == 1)) warning('WARNING neighbours are not linearly independent') 
  FUN.DataDims(data = StanList_Iter)
  
  #### Preferences
  rstan_options(auto_write = TRUE)
  rstan_options(javascript = FALSE)
  options(mc.cores = 4)
  
  ### Model Execution ----
  if(file.exists(file.path(Dir.TreatmentIter, "Model.RData"))){
    message("Model already compiled")
    load(file.path(Dir.TreatmentIter, "Model.RData"))
  }else{
    unlink(file.path(Dir.Base, "joint_model_MEE_FIA.exe"))
    unlink(file.path(Dir.Base, "joint_model_MEE_FIA.rds"))
    Stan_model <- rstan::stan(file = 'joint_model_MEE_FIA.stan',
                              data =  StanList_Iter,
                              chains = nChains,
                              warmup = round(nWarmup/1, 0),
                              iter = round(nSamples/1, 0),
                              cores = nChains,
                              include = TRUE,
                              pars = "ndd_betaij",
                              refresh = 10,
                              control = list(max_treedepth = 10,
                                             adapt_delta = 0.9),
                              model_name = BiomeName,
                              seed = 42
    )
    save(Stan_model, file = file.path(Dir.TreatmentIter, "Model.RData"))
  }
  # }
  # ### Model Diagnostics ----
  # # Get the full posteriors
  joint.post.draws <- extract.samples(Stan_model)
  # # Select parameters of interest
  # param.vec <- Stan_model@model_pars[Stan_model@model_pars %nin% c('response1', 'responseSm1', 'lp__')]
  # # Draw 1000 samples from the 80% posterior interval for each parameter of interest
  # p.samples <- list()
  # p.samples <- sapply(param.vec[param.vec %nin% c('ri_betaij', 'ndd_betaij',
  #                                                 "weight", "mu", "mu2"
  # )], function(p) {
  #   print(p)
  #   p.samples[[p]] <- apply(joint.post.draws[[p]], 2, function(x){
  #     sample(x[x > quantile(x, 0.1) & x < quantile(x, 0.9)], size = 100)
  #   })  # this only works for parameters which are vectors
  # })
  # # there is only one sigma_alph parameter so we must sample differently:
  # intrinsic.perf <- p.samples$gamma_i
  # colnames(intrinsic.perf) <- levels(factor(Index_df$taxon))
  # # inter_mat <- return_inter_array(joint.post.draws = joint.post.draws,
  # #                                 response = p.samples$response,
  # #                                 effect = p.samples$effect,
  # #                                 focalID = levels(factor(Index_Iter$taxon)),
  # #                                 neighbourID = colnames(Index_Iter[, -1:-3]))
  inter_mat <- aperm(joint.post.draws$ndd_betaij, c(2, 3, 1))
  rownames(inter_mat) <- levels(factor(Index_df$taxon))
  colnames(inter_mat) <- levels(factor(Index_df$taxon))
  # inter_mat is now a 3 dimensional array, where rows = focals, columns = neighbours and 3rd dim = samples from the posterior; inter_mat[ , , 1] should return a matrix consisting of one sample for every interaction 
  # try(stan_model_check(fit = Stan_model,
  #                      results_folder = Dir.TreatmentIter,
  #                      params = param.vec))
  
  ### Interaction/Association Matrix ----
  Interaction_mean <- apply(inter_mat, c(1, 2), mean) # will return the mean estimate for every interaction (NB: this is the mean of the 80% posterior interval, so will be slightly different to the mean value returned from summary(fit), which is calculated from the full posterior distribution)  
  Interaction_mean <- Interaction_mean*-1 # need to switch sign of results
  diag(Interaction_mean) <- NA
  Interaction_hpdi <- apply(inter_mat, c(1, 2), HPDI, prob = 0.89)
  Interaction_min <- -Interaction_hpdi[1,,]
  diag(Interaction_min) <- NA
  Interaction_max <- -Interaction_hpdi[2,,]
  diag(Interaction_max) <- NA
  Interactions_igraph <- data.frame(Actor = rep(dimnames(Interaction_mean)[[2]], length(dimnames(Interaction_mean)[[1]])),
                                    Subject = rep(dimnames(Interaction_mean)[[1]], each = length(dimnames(Interaction_mean)[[2]])),
                                    Inter_mean = as.vector(t(Interaction_mean)),
                                    Inter_min = as.vector(t(Interaction_min)),
                                    Inter_max = as.vector(t(Interaction_max))
  )
  Interactions_IFREM <- Interactions_igraph[order(abs(Interactions_igraph$Inter_mean), decreasing = TRUE), ]
  Interactions_IFREM <- na.omit(Interactions_IFREM)
  save(Interactions_IFREM, file = file.path(Dir.TreatmentIter, "Interac.RData"))
}

## NETASSOC ----------------------------------------------------------------
message("############ STARTING NETASSOC ANALYSES")
Dir.NETASSOC <- file.path(DirEx.Region, "NETASSOC")
if(!dir.exists(Dir.NETASSOC)){dir.create(Dir.NETASSOC)}

### Null model ranges ----
FIABiomass_df <- rbind(readRDS(file.path(Dir.FIA, "FIABiomes_df_BIOMES.rds")),
      readRDS(file.path(Dir.FIA, "FIABiomes_df_STATES.rds")))

if(file.exists(file.path(Dir.FIA, "Ranges_poly.RData"))){
  load(file.path(Dir.FIA, "Ranges_poly.RData"))
}else{
  Range_specs <- unique(FIABiomass_df$SCIENTIFIC_NAME)
  Ranges_spoly <- BIEN_ranges_load_species(species = Range_specs)
  Ranges_spoly <- st_as_sf(Ranges_spoly)
  Ranges_spoly <- sf::st_make_valid(Ranges_spoly)
  RangesFIA_spoly <- st_crop(Ranges_spoly, FIA_shp)
  save(RangesFIA_spoly, file = file.path(Dir.FIA, "Ranges_poly.RData"))
}

### Analysis loop ----
for(Treatment_Iter in Treatment_Vec){ # HMSC treatment loop
  load(file.path(Dir.FIA, paste0("FIABiome",Treatment_Iter,".RData")))
  message(paste("### Biome:", BiomeName, "(", nrow(ModelFrames_ls$FitCom), "Sites )"))
  
  Dir.TreatmentIter <- file.path(Dir.NETASSOC, Treatment_Iter)
  if(!dir.exists(Dir.TreatmentIter)){dir.create(Dir.TreatmentIter)}
  sink(file.path(Dir.TreatmentIter, "Biome.txt"))
  print("BIOME")
  print(BiomeName)
  print("OBSERVATIONS")
  print(nrow(ModelFrames_ls$Fitness))
  print("SPECIES")
  print(nrow(Phylo_Iter$Dist_Mean))
  print("SITES")
  print(nrow(Metadata_df))
  sink()
  
  
  mat_Iter <- ModelFrames_ls$Community[ , -1]
  
  
  Sites_sf <- FIABiomass_df[FIABiomass_df$pltID %in% ModelFrames_ls$Community$SiteID , ]
  Sites_sf <- Sites_sf[!duplicated(Sites_sf$pltID), ]
  RangesIntersec_sf <- st_intersects(Sites_sf, RangesFIA_spoly)
  
  Null_mat <- mat_Iter 
  Null_mat[] <- 0
  
  for(Null_Iter in 1:nrow(Null_mat)){
    Pres_spec <- RangesFIA_spoly$species[RangesIntersec_sf[Null_Iter][[1]]]
    Pres_spec <- gsub(Pres_spec, pattern = "_", replacement = " ")
    Null_mat[Null_Iter, which(colnames(Null_mat) %in% Pres_spec)] <- 1
  }
  
  rownames(mat_Iter) <- ModelFrames_ls$Community[ , 1]
  mat_Iter[is.na(mat_Iter)] <- 0
  # Removing species without observations
  Null_mat <- Null_mat[,colSums(mat_Iter) != 0]
  mat_Iter <- mat_Iter[,colSums(mat_Iter) != 0]
  # Removing species without range presence
  mat_Iter <- mat_Iter[,colSums(Null_mat) != 0]
  Null_mat <- Null_mat[,colSums(Null_mat) != 0]
  # Removing sites without species
  Null_mat <- Null_mat[rowSums(mat_Iter)!=0,]
  mat_Iter <- mat_Iter[rowSums(mat_Iter)!=0,]
  # Removing sites without ranges
  mat_Iter <- mat_Iter[rowSums(Null_mat)!=0,]
  Null_mat <- Null_mat[rowSums(Null_mat)!=0,]
  
  # data check
  if(nrow(mat_Iter) < 2){
    sink(file.path(Dir.TreatmentIter, "DataIssue.txt"))
    print("Not enough data")
    sink()
    next()
  }
  ### Model Execution ####
  if(file.exists(file.path(Dir.TreatmentIter, "Model.RData"))){
    print("Model already compiled")
    load(file.path(Dir.TreatmentIter, "Model.RData"))
  }else{
    model_netassoc <- make_netassoc_network(obs = t(mat_Iter), 
                                            plot = FALSE, verbose = FALSE)
    save(model_netassoc, file = file.path(Dir.TreatmentIter, "Model.RData"))
  }
  
  ### Interaction/Association Matrix ----
  # Interac_df <- model_netassoc$matrix_spsp_ses_thresholded
  Interac_df <- list(Effects = model_netassoc$matrix_spsp_ses_all,
                     p = model_netassoc$matrix_spsp_pvalue)
  Interac_df <- data.frame(Partner1 = rep(rownames(Interac_df$Effects), each = ncol(Interac_df$Effects)),
                           Partner2 = rep(colnames(Interac_df$Effects), nrow(Interac_df$Effects)),
                           effects = as.numeric(Interac_df$Effects), 
                           p = as.numeric(Interac_df$p),
                           Sig = as.numeric(Interac_df$p) < 0.05
  )
  
  Interac_df <- Interac_df[as.vector(upper.tri(model_netassoc$matrix_spsp_ses_all)), ]
  save(Interac_df, file = file.path(Dir.TreatmentIter, "Interac.RData")) 
}

## COOCCUR -----------------------------------------------------------------
message("############ STARTING COCCUR ANALYSES")
Dir.COOCCUR <- file.path(DirEx.Region, "COCCUR")
if(!dir.exists(Dir.COOCCUR)){dir.create(Dir.COOCCUR)}

for(Treatment_Iter in Treatment_Vec){ # HMSC treatment loop
  load(file.path(Dir.FIA, paste0("FIABiome",Treatment_Iter,".RData")))
  message(paste("### Biome:", BiomeName, "(", nrow(ModelFrames_ls$FitCom), "Sites )"))
  Dir.TreatmentIter <- file.path(Dir.COOCCUR, Treatment_Iter)
  if(!dir.exists(Dir.TreatmentIter)){dir.create(Dir.TreatmentIter)}
  sink(file.path(Dir.TreatmentIter, "Biome.txt"))
  print("BIOME")
  print(BiomeName)
  print("OBSERVATIONS")
  print(nrow(ModelFrames_ls$Fitness))
  print("SPECIES")
  print(nrow(Phylo_Iter$Dist_Mean))
  print("SITES")
  print(nrow(Metadata_df))
  sink()
  
  mat_Iter <- ModelFrames_ls$Community[,-1]
  rownames(mat_Iter) <- ModelFrames_ls$Community[ , 1]
  mat_Iter[is.na(mat_Iter)] <- 0
  mat_Iter <- mat_Iter[colSums(mat_Iter) != 0]
  ### Model Execution ####
  if(file.exists(file.path(Dir.TreatmentIter, "Model.RData"))){
    print("Model already compiled")
    load(file.path(Dir.TreatmentIter, "Model.RData"))
  }else{
    model_coccurr <- cooccur(mat = t(mat_Iter), 
                             type = "spp_site", thresh = FALSE, spp_names = TRUE)
    save(model_coccurr, file = file.path(Dir.TreatmentIter, "Model.RData"))
  }
  
  ### Interaction/Association Matrix ----
  # test <- tryCatch(plot(model_coccurr))
  # if(nrow(test$data) != 0){
  #   jpeg(file=file.path(Dir.TreatmentIter, "Cooccur_Assocs.jpeg"), width = 32, height = 32, units = "cm", res = 100)
  #   print(test)
  #   dev.off() 
  # }
  Interac_df <- effect.sizes(model_coccurr, standardized = TRUE)
  Interac_df$pLT <- prob.table(model_coccurr)$p_lt
  Interac_df$pGT <- prob.table(model_coccurr)$p_gt
  Interac_df$Sig <- Interac_df[,4] < 0.05 | Interac_df[,5] < 0.05
  colnames(Interac_df)[1:2] <- c("Partner1", "Partner2")
  save(Interac_df, file = file.path(Dir.TreatmentIter, "Interac.RData")) # In standardized form, these values are bounded from -1 to 1, with positive values indicating positive associations and negative values indication negative associations; see 10.18637/jss.v069.c02
}

