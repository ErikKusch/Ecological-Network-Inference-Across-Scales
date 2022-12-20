#' ######################################################################## #
#' PROJECT: [Data Simplification for Network Inference across Scales] 
#' CONTENTS: 
#'  - Generate YFDP species-association/-interaction networks
#'  DEPENDENCIES:
#'  - 0 - Preamble.R
#'  - X - Functions_Bayes.R
#'  - X - Functions_Data.R
#'  - X - Functions_Plotting.R
#'  - yfdp_CTFS_Census_1_20211109_for_Erik.csv (available through http://yfdp.org/index.html upon request)
#'  - yfdp_CTFS_Census_2_20211109_for_Erik.csv (available through http://yfdp.org/index.html upon request)
#'  - YFDP_Tagged_Species_Traits.csv (available here: yfdp.org/Summary/YFDP_Tagged_Species_Traits.csv)
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

## HOUSEKEEPING ------------------------------------------------------------
## lookup objects for species acronyms and treatment identifiers
Species_df <- data.frame(Short = c("PILA", "ABCO", "CADE", 
                                   "CONU", "PRVI", "QUKE", 
                                   "FRCA", "PIPO", "SASC", 
                                   "ABMA", "PREM", "COSE", 
                                   "UNKN", "PSME", "COCOC"),
                         Long = c("Pinus lambertiana", "Abies concolor", "Calocedrus decurrens",
                                  "Cornus nuttallii", "Prunus virginiana", "Quercus kelloggii",
                                  "Frangula californica", "Pinus ponderosa", "Salix scouleriana",
                                  "Abies magnifica", "Prunus emarginata", "Cornus sericea", 
                                  "Unknown", "Pseudotsuga menziesii", "Corylus cornuta")
)
Treatments_ls <- list(Name = c("Pre-Fire", "Post-Fire"),
                      File = c("yfdp_CTFS_Census_1_20211109_for_Erik.csv",
                               "yfdp_CTFS_Census_2_20211109_for_Erik.csv")
)

## FUNCTIONAL TRAITS -------------------------------------------------------
## download YFDP data if needed
if(!file.exists(file.path(Dir.YFDP, "Traits_df.csv"))){
  message("Obtaining trait data")
  Raw_df <- read.csv(file.path(Dir.YFDP, "YFDP_Tagged_Species_Traits.csv"))
  Raw_df$taxon <- Species_df$Long[match(Raw_df$CODE, Species_df$Short)]
  Raw_df[!is.na(Raw_df$taxon), ]
  Raw_df <- Raw_df[Raw_df$taxon %in% Phylo_Specs, ] # limitting to phylogeny-recognised species
  Traits_df <- data.frame(
    taxon = Raw_df$taxon,
    Shape = Raw_df$SHAPE,
    SLA = as.numeric(Raw_df$SLA),
    LEAF_C = as.numeric(Raw_df$LEAF_C),
    LEAF_N = as.numeric(Raw_df$LEAF_N)
  )
  Traits_df <- Traits_df[!is.na(Traits_df$LEAF_C), ]
  Traits_df <- Traits_df[order(Traits_df$taxon),]
  write.csv(Traits_df, file = file.path(Dir.YFDP, "Traits_df.csv"))
}else{
  message("Trait data already obtained")
  Traits_df <- read.csv(file.path(Dir.YFDP, "Traits_df.csv"))[,-1]
}

## REFOMRATTING RAW DATA ---------------------------------------------------
for(Treatment_Iter in 1:length(Treatments_ls$Name)){
  message(paste("Preparing data for", Treatments_ls$Name[Treatment_Iter], "treatment"))
  Raw_df <- read.csv(file.path(Dir.YFDP, Treatments_ls$File[Treatment_Iter]))
  Raw_df <- Raw_df[Raw_df$SPECIES != "UNKN", ]
  Raw_df <- Raw_df[Raw_df$DBH > 1, ]
  if(Treatment_Iter == 2){ # remove DBHs that have been altered by fire in post-fire data
    Raw_df <- Raw_df[Raw_df$BURN_AT_DBH != 3, ]
  }
  YFDP_df <- data.frame(SiteID = Raw_df$QUADRAT,
                        taxon = Raw_df$SPECIES,
                        value = Raw_df$DBH,
                        Lon = Raw_df$UTM_X,
                        Lat = Raw_df$UTM_Y,
                        Year = year(Raw_df$DATE))
  YFDP_df$taxon <- Species_df$Long[match(YFDP_df$taxon, Species_df$Short)]
  YFDP_df[, c("Lon", "Lat")] <- ConvertCoordinates(YFDP_df$Lon, YFDP_df$Lat)
  YFDP_df <- YFDP_df[YFDP_df$taxon %in% Traits_df$taxon, ]
  write.csv(YFDP_df, file = file.path(Dir.YFDP, paste0(Treatments_ls$Name[Treatment_Iter], ".csv")))
}

## GENERATING PHYLOGENY ---------------------------------------------------
if(!file.exists(file.path(Dir.YFDP, "Phylogeny.RData"))){
  Metadata1_df <- read.csv(file.path(Dir.YFDP, paste0(Treatments_ls$Name[1], ".csv")))
  Metadata2_df <- read.csv(file.path(Dir.YFDP, paste0(Treatments_ls$Name[2], ".csv")))
  Metadata_df <- rbind(Metadata1_df, Metadata2_df)
  message("Generating phylogeny")
  Phylo_ls <- FUN.PhyloDist(SpeciesNames = Metadata_df$taxon)
  save(Phylo_ls, file = file.path(Dir.YFDP, "Phylogeny.RData"))
}else{
  message("Phylogeny already generated")
  load(file.path(Dir.YFDP, "Phylogeny.RData"))
}
Phylo_Specs <- gsub(pattern = "_", replacement = " ", x =  Phylo_ls$Avg_Phylo$tip.label)
Phylo_ls$Avg_Phylo$tip.label <- Phylo_Specs

## META & CLIMATE DATA ----------------------------------------------------
ECV_vec <- c("2m_temperature", "volumetric_soil_water_layer_1", "total_precipitation", "potential_evaporation")
if(sum(file.exists(file.path(Dir.YFDP, paste0(Treatments_ls$Name, "_Metadata.csv")))) < length(Treatments_ls$Name)){ # bioclimatic data not established yet
  message("Obtaining covariate data")
  Metadata1_df <- read.csv(file.path(Dir.YFDP, paste0(Treatments_ls$Name[1], ".csv")))
  Metadata1_df$Census <- "Pre"
  Metadata2_df <- read.csv(file.path(Dir.YFDP, paste0(Treatments_ls$Name[2], ".csv")))
  Metadata2_df$Census <- "Post"
  Metadata_df <- rbind(Metadata1_df, Metadata2_df)
  Metadata_df$X <- 1:nrow(Metadata_df)
  
  ## CLIMATE DATA DOWNLOAD
  Shp <- KrigR::buffer_Points(Points = Metadata_df, Buffer = 0.1, ID = "X")
  for(Clim_Iter in 1:length(ECV_vec)){
    FUN.CLIM(ECV = Clim_Iter, Shp = Shp, Dir = Dir.YFDP)
  }
  coordinates(Metadata_df) <- ~Lon+Lat
  
  ## CLIMATE DATA EXTRACTION 
  Layer_seq <- seq.Date(as.Date("1981-01-01"), as.Date("2020-12-31"), by = "month")
  for(Clim_Iter in 1:length(ECV_vec)){
    print(ECV_vec[Clim_Iter])
    Metadata_df$XYZ <- NA
    names(Metadata_df)[ncol(Metadata_df)] <- ECV_vec[Clim_Iter]
    Metadata_df$XYZ <- NA
    names(Metadata_df)[ncol(Metadata_df)] <- paste0(ECV_vec[Clim_Iter], "_SD")
    Metadata_df$XYZ <- NA
    names(Metadata_df)[ncol(Metadata_df)] <- paste0(ECV_vec[Clim_Iter], "_UC")
    Extrac_temp <- raster::extract(stack(file.path(Dir.YFDP, paste0(ECV_vec[Clim_Iter], ".nc"))), 
                                   Metadata_df)
    Uncert_temp <- raster::extract(stack(file.path(Dir.YFDP, paste0("UC", ECV_vec[Clim_Iter], ".nc"))),
                                   Metadata_df)
    pb <- txtProgressBar(min = 0, max = nrow(Extrac_temp), style = 3) 
    for(Plot_Iter in 1:nrow(Extrac_temp)){
      ## only retain the last ten years leading up to data collection
      Need_seq <- seq.Date(as.Date(paste0(Metadata_df[Plot_Iter, ]$Year-10, "-01-01")), 
                           as.Date(paste0(Metadata_df[Plot_Iter, ]$Year, "-01-01")), 
                           by = "month")
      Time_seq <- Extrac_temp[Plot_Iter, which(Layer_seq %in% Need_seq)]
      Uncert_seq <- Uncert_temp[Plot_Iter, which(Layer_seq %in% Need_seq)]
      Metadata_df[Plot_Iter, ECV_vec[Clim_Iter]] <- mean(Time_seq, na.rm = TRUE)
      Metadata_df[Plot_Iter, paste0(ECV_vec[Clim_Iter], "_SD")] <- sd(Time_seq, na.rm = TRUE)
      Metadata_df[Plot_Iter, paste0(ECV_vec[Clim_Iter], "_UC")] <- mean(Uncert_seq, na.rm = TRUE)
      setTxtProgressBar(pb, Plot_Iter)
    }
  }
  
  ## average by plot and census
  Metadata_df <- aggregate(. ~ SiteID + Census, 
                           data = as.data.frame(Metadata_df[,c(-1,-3:-5)]), 
                           FUN = mean)
  
  ## break apart by Census column
  Metadata_ls <- split(Metadata_df, Metadata_df$Census)
  names(Metadata_ls) <- sort(Treatments_ls$Name)
  
  ## writing results
  sapply(names(Metadata_ls), 
    function (x) write.csv(Metadata_ls[[x]], file=file.path(Dir.YFDP, paste0(x, "_Metadata.csv")))
    )
}else{ # bioclimatic data has been established
  message("Covariate data already retrieved")
}
# set this because 2m_temperature is stored as x2m_temperature in Metadata_df
ECV_vec[ECV_vec == "2m_temperature"] <- "X2m_temperature"

## MODEL DATA FRAMES -------------------------------------------------------
if(sum(file.exists(file.path(Dir.YFDP, paste0(Treatments_ls$Name, "_ModelFrames.RData")))) < length(Treatments_ls$Name)){
  message("Readying model data frames")
  for(Treatment_Iter in 1:length(Treatments_ls$Name)){
    message(paste("Preparing data for", Treatments_ls$Name[Treatment_Iter], "treatment"))
    Raw_df <- read.csv(file.path(Dir.YFDP, paste0(Treatments_ls$Name[Treatment_Iter], ".csv")))
    Raw_df <- Raw_df[Raw_df$taxon %in% Phylo_Specs, ] # limiting to phylogeny-recognised species
    ### FITNESS ####
    Weight_df <- Raw_df[ , c("SiteID", "taxon", "value")] # only dry biomass rows and select only relevant columns for speedier data handling 
    Weight_cast <- reshape2::dcast(Weight_df, SiteID ~ taxon, value.var = 'value', fun.aggregate = mean)
    ## ABUNDANCE ####
    Abund_cast <- dcast(Raw_df, SiteID ~ taxon, fun.aggregate = length)
    ## LIST MERGING ####
    ModelFrames_ls <- list(
      Fitness = Weight_df,
      FitCom = Weight_cast,
      Community = Abund_cast
    )
    save(ModelFrames_ls, file = file.path(Dir.YFDP, paste0(Treatments_ls$Name[Treatment_Iter], "_ModelFrames.RData")))
  }
}else{
  message("Model data frames already compiled")
}

# ANALYSIS =================================================================
## HMSC --------------------------------------------------------------------
message("############ STARTING HMSC ANALYSES")
Dir.HMSC <- file.path(DirEx.YFDP, "HMSC")
if(!dir.exists(Dir.HMSC)){dir.create(Dir.HMSC)}

for(Treatment_Iter in Treatments_ls$Name[1]){ # HMSC treatment loop
  message(paste("### Treatment:", Treatment_Iter))
  Dir.TreatmentIter <- file.path(Dir.HMSC, Treatment_Iter)
  if(!dir.exists(Dir.TreatmentIter)){dir.create(Dir.TreatmentIter)}
  ### DATA LOADING ####
  load(file.path(Dir.YFDP, paste0(Treatment_Iter, "_ModelFrames.RData")))
  ModelFrames_Iter <- ModelFrames_ls
  Metadata_Iter <- read.csv(file.path(Dir.YFDP, paste0(Treatment_Iter, "_Metadata.csv")))
  Phylo_Iter <- Phylo_ls$Avg_Phylo
  Traits_Iter <- read.csv(file.path(Dir.YFDP, "Traits_df.csv"))[,-1]
  rownames(Traits_Iter) <- Traits_Iter$taxon
  Traits_Iter <- Traits_Iter[,-1]

  ### DATA PREPRATION ####
  S <- Metadata_Iter[ , c("SiteID", "Lat", "Lon")] # S: study design, including units of study and their possible coordinates, If you don't have variables that define the study design, indicate this by S=NULL
  X <- Metadata_Iter[ , c("SiteID", ECV_vec)] # X: covariates to be used as predictors, If you don't have covariate data, indicate this by X=NULL
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
  X <- X[, c(ECV_vec)]
  ## Model Formulae
  XFormula <- as.formula(paste("~ ", paste(ECV_vec[1:3], collapse = " + "), sep = " + "))
  TrFormula <- ~SLA + LEAF_C + LEAF_N
  ## StudyDesign
  unique_plot <- paste(S$SiteID, sep="_")
  studyDesign <- data.frame(site = as.factor(S$SiteID))
  St <- studyDesign$site
  rL.site <- HmscRandomLevel(units = levels(St))
  ## Model Objects
  Ypa <- 1*(Y_AB>0)
  Yabu <- Y_AB
  # Yabu[Y_AB==0] <- NA
  # Yabu <- log(Yabu)
  Ybiom <- Y_BM
  Ybiom[is.na(Y_BM)] <- 0
  # Ybiom <- log(Ybiom)
  ## Models
  m1 <- Hmsc(Y=Ypa, XData = X,  XFormula = XFormula,
             TrData = Tr,TrFormula = TrFormula,
             distr="probit",
             studyDesign=studyDesign,
             ranLevels={list("site" = rL.site)})
  m2 <- Hmsc(Y=Yabu, YScale = TRUE,
             XData = X,  XFormula = XFormula,
             TrData = Tr,TrFormula = TrFormula,
             distr="normal",
             studyDesign=studyDesign,
             ranLevels={list("site" = rL.site)})
  m3 <- Hmsc(Y=Ybiom, YScale = TRUE,
             XData = X,  XFormula = XFormula,
             TrData = Tr,TrFormula = TrFormula,
             distr="normal",
             studyDesign=studyDesign,
             ranLevels={list("site" = rL.site)})
  models <- list(m1,m2, m3)
  modelnames <- c("presence_absence","abundance", "diametre")

  ### Modelling ----
  samples_list <- nSamples
  message(paste0("thin = ",as.character(thin),"; samples = ",as.character(nSamples)))
  for(HMSC_Iter in 1:length(modelnames)){ # HMSC model loop
    hmsc_modelname <- modelnames[HMSC_Iter]
    if(file.exists(file.path(Dir.TreatmentIter, paste0(hmsc_modelname, "_Interac.RData")))){
      message("Already computed and evaluated")
      next()
    }
    #### Model Execution ----
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
      message("Model already compiled")
      load(filename)
    }

    ### Model Evaluation ----
    message("Evaluation")
    vals <- HMSC.Eval(hmsc_model = hmsc_model, Dir = Dir.TreatmentIter, Name = hmsc_modelname, thin = thin, nSamples = nSamples, nChains = nChains)
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

## IF-REM ------------------------------------------------------------------
message("############ STARTING IF-REM ANALYSES")
Dir.IFREM <- file.path(DirEx.YFDP, "IF_REM")
if(!dir.exists(Dir.IFREM)){dir.create(Dir.IFREM)}

for(Treatment_Iter in Treatments_ls$Name[1]){ # HMSC treatment loop
  message(paste("### Treatment:", Treatment_Iter))
  Dir.TreatmentIter <- file.path(Dir.IFREM, Treatment_Iter)
  if(file.exists(file.path(Dir.TreatmentIter, "Interac.RData"))){
    message("Model already compiled and evaluated")
    next()
  }
  if(!dir.exists(Dir.TreatmentIter)){dir.create(Dir.TreatmentIter)}
  if(file.exists(file.path(Dir.TreatmentIter, "Interac.RData"))){
    message("Already computed and evaluated")
    next()
  }
  ### DATA LOADING ####
  load(file.path(Dir.YFDP, paste0(Treatment_Iter, "_ModelFrames.RData")))
  
  Index_df <- cbind(ModelFrames_ls$Fitness, 
                    ModelFrames_ls$Community[match(ModelFrames_ls$Fitness$SiteID,
                                                   ModelFrames_ls$Community$SiteID),  -1])
  Index_df <- Index_df[, c(1:3, which(colSums(Index_df[, -1:-3]) != 0) + 3)] # ensuring only present species are retained
  Index_df <- Index_df[which(Index_df$value > 0), ]
  
  ### DATA PREPRATION ####
  StanList_Iter <- data_prep(perform = "value",
                             focal = "taxon",
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
    unlink(file.path(Dir.Base, "joint_model_MEE.exe"))
    unlink(file.path(Dir.Base, "joint_model_MEE.rds"))
    Stan_model <- rstan::stan(file = 'joint_model_MEE.stan',
                       data =  StanList_Iter,
                       chains = nChains,
                       warmup = round(nWarmup/1, 0),
                       iter = round(nSamples/1, 0),
                       cores = nChains,
                       include = TRUE,
                       pars = "ndd_betaij",
                       refresh = 100,
                       control = list(max_treedepth = 10,
                                      adapt_delta = 0.9),
                       # model_name = BiomeName,
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
Dir.NETASSOC <- file.path(DirEx.YFDP, "NETASSOC")
if(!dir.exists(Dir.NETASSOC)){dir.create(Dir.NETASSOC)}

for(Treatment_Iter in Treatments_ls$Name[1]){ # HMSC treatment loop
  message(paste("### Treatment:", Treatment_Iter))
  Dir.TreatmentIter <- file.path(Dir.NETASSOC, Treatment_Iter)
  if(!dir.exists(Dir.TreatmentIter)){dir.create(Dir.TreatmentIter)}
  ### DATA LOADING ####
  load(file.path(Dir.YFDP, paste0(Treatment_Iter, "_ModelFrames.RData")))
  ModelFrames_Iter <- ModelFrames_ls
  mat_Iter <- ModelFrames_Iter$Community[ , -1]
    # ModelFrames_Iter$FitCom[, -1]
  rownames(mat_Iter) <- ModelFrames_Iter$Community[ , 1]
  mat_Iter[is.na(mat_Iter)] <- 0
  mat_Iter <- mat_Iter[colSums(mat_Iter) != 0]
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
Dir.COOCCUR <- file.path(DirEx.YFDP, "COCCUR")
if(!dir.exists(Dir.COOCCUR)){dir.create(Dir.COOCCUR)}

for(Treatment_Iter in Treatments_ls$Name[1]){ # HMSC treatment loop
  message(paste("### Treatment:", Treatment_Iter))
  Dir.TreatmentIter <- file.path(Dir.COOCCUR, Treatment_Iter)
  if(!dir.exists(Dir.TreatmentIter)){dir.create(Dir.TreatmentIter)}
  ### DATA LOADING ####
  load(file.path(Dir.YFDP, paste0(Treatment_Iter, "_ModelFrames.RData")))
  ModelFrames_Iter <- ModelFrames_ls
  mat_Iter <- ModelFrames_Iter$Community[ , -1]
  rownames(mat_Iter) <- ModelFrames_Iter$Community[ , 1]
  mat_Iter[is.na(mat_Iter)] <- 0
  mat_Iter[mat_Iter > 1] <- 1
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
