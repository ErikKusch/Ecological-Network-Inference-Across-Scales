#' ####################################################################### #
#' PROJECT: [PhD; X - BAYESIAN FUNCTIONALITY] 
#' CONTENTS: 
#'  - Functionality for execution and inspection of intrinsic fitness method
#'  DEPENDENCIES:
#'  - 
#' AUTHOR: [Malyon Bimler, Erik Kusch]
#' ####################################################################### #

# PREPARING DATA ===========================================================
## KUSCH -------------------------------------------------------------------
# takes character and data frame arguments and returns a list object in preparation for network model execution
FUN.StanList <- function(
  Outcome = "Outcome", # name of the outcome variable column in Index_df
  Index_df = NULL, # data frame containing columns for ID, Fitness, Site, and Species
  Neigh_df = NULL, # data frame containing column for site and columns for all species; cells containing containing predictor variable of each species at each site
  Neigh_UC = NULL,
  Envir_df = NULL, # data frame containing column for site and columns for climate variables; cells containing climate variable values at each site
  Envir_UC = NULL,
  Phylo_df = NULL, # phylogenetic distances
  Phylo_UC = NULL,
  Traits_df = NULL
  ){
  # Basic List and Data Cleaning ----
  stan.data <- list() # set up the data in list format as preferred by STAN
  Index_df[,Outcome] <- as.numeric(Index_df[,Outcome]) # ensure that outcome is numeric
  Index_df$Species <- as.character(Index_df$Species) # ensure that species are stored as characters
  Index_df$SiteID <- as.character(Index_df$SiteID) # ensure that sites are stored as characters
  Index_df$ObsID <- 1:nrow(Index_df) # establish observation ID column
  Neigh_df$SiteID <- as.character(Neigh_df$SiteID) # ensure SiteIDs aren't factors
  Neigh_UC$SiteID <- as.character(Neigh_UC$SiteID) # ensure SiteIDs aren't factors
  Envir_df$SiteID <- as.character(Envir_df$SiteID) # ensure SiteIDs aren't factors
  Envir_UC$SiteID <- as.character(Envir_UC$SiteID) # ensure SiteIDs aren't factors
  
  # Data frame sorting ----
  Neigh_df <- Neigh_df[order(Neigh_df$SiteID), ]
  Neigh_UC <- Neigh_UC[order(Neigh_UC$SiteID), ]
  Envir_df <- Envir_df[order(Envir_df$SiteID), ]
  Envir_UC <- Envir_UC[order(Envir_UC$SiteID), ]
  
  # Number of Interactions ----
  SitesMatch <- base::match(Index_df$SiteID, Neigh_df$SiteID) # find matching sites between index and neighbours
  Interac_df <- cbind(Index_df,Neigh_df[SitesMatch,-1]) # build large data frame containing all neighbour predictors
  Interac_df <- Interac_df[ , -c(1:4)] # remove first four columns (SiteID, Species, Outcome, ObsID)
  Interac_df[Interac_df > 0] <- 1 # set all predictors which indicate presence to 1
  Interac_df[is.na(Interac_df)] <- 0 # set NAs to zero to indicate no observed interaction
  Counts <- split(Interac_df, Index_df$Species) # split into groups by species
  obs <- do.call(rbind, lapply(Counts, colSums)) # count observed placement at same plot(s) for all species pairs
  
  # Data Dimensions ----
  ## Integers ----
  stan.data$N <- nrow(Index_df) # number of observations
  stan.data$N_Species <- length(unique(Index_df$Species)) # number of species
  stan.data$N_Sites <- nrow(Envir_df) # number of sites
  stan.data$N_Neigh <- ncol(Neigh_df)-1 # number of neighbours
  stan.data$N_Interac <- length(obs[obs > 0]) # number of observed interactions
  stan.data$N_Clim <- ncol(Envir_df)-1 # number of climate variables
  
  ## Vectors ----
  stan.data$ID_Species <- as.numeric(as.factor(Index_df$Species)) # numeric indices for focal species
  stan.data$ID_Site <- as.numeric(as.factor(Index_df$SiteID)) # numeric indices for sites
  stan.data$fitness <- Index_df[, Outcome] # fitness proxy
  
  ## Indices for later back-translation of numeric values ----
  ## named factor vector of species names and numeric IDs
  stan.data$Index$Species <- as.factor(Index_df$Species)
  names(stan.data$Index$Species) <- as.numeric(as.factor(Index_df$Species))
  ## named factor vector of site IDs names and numeric IDs
  stan.data$Index$Sites <- as.factor(Envir_df$SiteID)
  names(stan.data$Index$Sites) <- as.numeric(as.factor(Envir_df$SiteID))
  
  ## Indices for Interactions in Alpha-Matrix (matching species to interactions) ----
  stan.data$inter_per_species <- obs # first count the number of interactions observed for each focal species
  stan.data$inter_per_species[stan.data$inter_per_species > 0] <- 1 # set all interactions to one
  stan.data$inter_per_species <- rowSums(stan.data$inter_per_species, na.rm = TRUE) # count number of interactions per species
  stan.data$icol <- unlist(apply(ifelse(obs > 0, TRUE, FALSE), 1, which)) # column index in the alpha matrix for each observed interaction
  names(stan.data$icol) <- NULL # remove names of object
  stan.data$irow <- rep(1, stan.data$inter_per_species[[1]]) # begin the row index
  stan.data$istart <- 1 # begin the start and end indices for the vector of interactions per species 
  stan.data$iend <- stan.data$inter_per_species[[1]] # begin the start and end indices for the vector of interactions per species 
  for(s in 2:stan.data$N_Species){ # populate indices for all the other species
    # starting position of a_ij's for i in the vector of observed interactions (ie the 1st 'j')
    stan.data$istart[s] <- sum(stan.data$inter_per_species[s-1:s])+1
    # end position of a_ij's for i in the vector of observed interactions (the last 'j')
    stan.data$iend[s] <-  sum(stan.data$inter_per_species[1:s])
    # row index in the alpha matrix for each observed interaction
    stan.data$irow <- c(stan.data$irow, rep(s, stan.data$inter_per_species[[s]]))
  }
  
  # ## Number of observations per interaction ----  I did not see this needed in the actual model
  # stan.data$Obs <- as.vector(apply(obs, 1, c)) # vector of the number of observations for each interactions
  # stan.data$Obs <- stan.data$Obs[stan.data$Obs > 0] # remove unobserved interactions
  
  # Model Matrices ----
  ## Neighbourhood Matrix
  stan.data$Mat_Neigh <- as.matrix(Neigh_df[ , -1])
  stan.data$UC_Neigh <- as.matrix(Neigh_UC[ , -1])
  # stan.data$Mat_Neigh <- cbind(as.numeric(as.factor(Neigh_df$SiteID)), stan.data$Mat_Neigh)
  # colnames(stan.data$Mat_Neigh)[1] <- "SiteID"
  ## Environment Matrix
  stan.data$Mat_Clim <- as.matrix(Envir_df[ , -1])  
  stan.data$UC_Clim <- as.matrix(Envir_UC[ , -1])  
  # stan.data$Mat_Clim <- cbind(as.numeric(as.factor(Envir_df$SiteID)), stan.data$Mat_Clim)
  # colnames(stan.data$Mat_Clim)[1] <- "SiteID"
  ## Phylogenetic Distance
  stan.data$Mat_Phylo <- Phylo_df
  stan.data$UC_Phylo <- Phylo_UC
  ## Trait similarity
  stan.data$Mat_Traits <- Traits_df
  
  # Return list ----
  return(stan.data)
}

# ## BIMLER ------------------------------------------------------------------
# # works on a data frame where the first four columns are plotID, fitness proxy (identified with Fitness argument), focal (species membership), focalID (speciesID and plotID)
# Fun_StanList <- function(Fitness = "fit", data = NULL){
#   ## Basic List and Data Cleaning ----
#   stan.data <- list() # set up the data in list format as preferred by STAN
#   data[ , -c(1:4)] <- sapply(data[ , -c(1:4)], as.numeric)
#   data$focal <- as.character(data$focal)
#   
#   ## Number of Observations ----
#   Counts <- data[ , -c(1:4)] # remove first four columns (plot, fitness proxy, focal, focalID)
#   Counts[Counts > 0] <- 1 # set all counts/abundances which indicate presence to 1
#   Counts <- split(Counts, data$focal) # split into groups by species
#   obs <- do.call(rbind, lapply(Counts, colSums))
#   
#   ## Integers ----
#   stan.data$S <- length(unique(data$focal))  # number of species
#   stan.data$N <- nrow(data)                  # number of observations
#   stan.data$K <- ncol(data[ , -c(1:4)])      # number of neighbours
#   stan.data$I <- length(obs[obs > 0])        # number of observed interactions
#   stan.data$Z <- stan.data$S*stan.data$K   # total number of possible interactions         
#   
#   ## Vectors ----
#   stan.data$species_ID <- as.numeric(as.factor(data$focal)) # numeric indices for focal species
#   stan.data$fitness <- data[, Fitness] # fitness proxy
#   
#   ## Indices for Interactions in Alpha-Matrix ----
#   stan.data$inter_per_species <- obs # first count the number of interactions observed for each focal species
#   stan.data$inter_per_species[stan.data$inter_per_species > 0] <- 1 # set all interactions to one
#   stan.data$inter_per_species <- rowSums(stan.data$inter_per_species)
#   stan.data$inter_per_species_perneighbour <- obs
#   stan.data$icol <- unlist(apply(ifelse(obs > 0, T, F), 1, which)) # column index in the alpha matrix for each observed interaction
#   names(stan.data$icol) <- NULL
#   stan.data$icol <- as.vector(stan.data$icol)
#   stan.data$irow <- rep(1, stan.data$inter_per_species[[1]]) # begin the row index
#   stan.data$istart <- 1 # begin the start and end indices for the vector of interactions per species 
#   stan.data$iend <- stan.data$inter_per_species[[1]] # begin the start and end indices for the vector of interactions per species 
#   for(s in 2:stan.data$S){ # populate indices for all the other species
#     # starting position of a_ij's for i in the vector of observed interactions (ie the 1st 'j')
#     stan.data$istart[s] <- sum(stan.data$inter_per_species[s-1:s])+1
#     # end position of a_ij's for i in the vector of observed interactions (the last 'j')
#     stan.data$iend[s] <-  sum(stan.data$inter_per_species[1:s])
#     # row index in the alpha matrix for each observed interaction
#     stan.data$irow <- c(stan.data$irow, rep(s, stan.data$inter_per_species[[s]]))
#   }
#   
#   ## Model Matrix ----
#   stan.data$X <- as.matrix(data[ , -c(1:4)])  
#   
#   ## Number of observations per interaction ----
#   stan.data$Obs <- as.vector(apply(obs, 1, c)) # vector of the number of observations for each interactions
#   stan.data$Obs <- stan.data$Obs[stan.data$Obs > 0] # remove unobserved interactions
#   
#   return(stan.data)
# }

# CHECKING DATA ============================================================
## KUSCH -------------------------------------------------------------------
FUN.DataDims <- function(data = NULL){
  message(paste0('Number of observations = ', data$N))
  message(paste0('Number of focal species = ', data$N_Species))
  message(paste0('Number of neighbouring species = ', data$N_Neigh))
  message(paste0('Number of observed interactions = ', data$N_Interac))
  message(paste0('Number of sites = ', data$N_Sites))
  message(paste0('Number of environmental variables = ', data$N_Clim))
}

# ## BIMLER ------------------------------------------------------------------
# Fun_PreCheck <- function(data = NULL){
#   key_speciesID <- unique(data$focal)
#   key_neighbourID <- colnames(data[ , -c(1:4)])
#   # message(paste0('Dataset selected: ', comm))
#   message(paste0('Fitness data dimensions = ', dim(data)[1], ', ', dim(data)[2]))
#   message(paste0('Number of focals = ', length(key_speciesID)))
#   message(paste0('Number of neighbours = ', length(key_neighbourID)))
# }

################# MODEL SPECIFICATION ###
# model specified in 'Supplement - StanModel.stan'

################# MODEL INSPECTION ###
# This function takes the posterior draws for interaction estimates extracted from 
# the STAN model fit object and returns a focal x neighbour x sample array of all 
# interactions, both observed (IFM) and unrealised (RIM)
return_inter_array <- function(joint.post.draws, # posterior draws extracted using the extract.samples() function
                               response = p.samples$response, # samples for the response parameters
                               effect = p.samples$effect,     # samples for the effect parameters
                               focalID,   # vector of focal identities (names)
                               neighbourID,  # vector of neighbour identities (names)
                               ...){
  
  
  # extract the IFM interaction estimates - though observed interactions are already sampled from 
  # the 'beta_ij' parameters, sampling from inter_mat has the benefit of including columns of 0's for 
  # unrealised interactions (which can then be filled with the corresponding RIM estimates)
  ifm_inters <- joint.post.draws$inter_mat
  ifm_inters <- as.data.frame(aperm(ifm_inters, perm = c(1, 3, 2)))
  # there is now one column per interaction - unrealised interactions will simply be a column of 0's
  # sample from the 80% posterior interval for each interaction
  ifm_inters <- apply(ifm_inters, 2, function(x) {
    inter <- x[x > quantile(x, 0.1) & x < quantile(x, 0.9)]
    if (length(inter > 0)) {sample(inter, size = 1000)} 
    else {rep(0, 1000)} # this is for those unobserved interactions (0)
  })
  
  # calculate unrealised interactions 
  rim_inters <- lapply(c(1:length(focalID)), function(x) {
    # randomly re-order samples from the response parameter for each focal 
    r <- sample(response[ , x], dim(response)[1]) 
    # multiply focal's response by all neighbour effect parameters 
    return(r*effect)})
  rim_inters <- do.call(cbind, rim_inters) # this returns RI model estimates for every possible interaction 
  # every column is an interaction
  
  # replace unobserved interactions (columns of 0 in ifm_inters) with the values predicted by the RIM
  all_inters <- ifm_inters
  all_inters[ , apply(all_inters, 2, function(x) {all(x == 0)})] <- 
    rim_inters[ , apply(all_inters, 2, function(x) {all(x == 0)})]
  length(colSums(all_inters)[colSums(all_inters) == 0]) # this should be equal to 0 
  
  # reconstruct species x neighbour x sample interaction arrays  
  inter_mat <- array(data = unlist(all_inters), 
                     dim = c(nrow(all_inters), length(neighbourID), length(focalID)), 
                     dimnames = list('samples' = seq(1, nrow(all_inters)), 
                                     'neighbour' = neighbourID, 
                                     'species' = focalID))
  inter_mat <- aperm(inter_mat, c(3,2,1))
  # inter_mat is now a 3 dimensional array, where rows = focals, columns = neighbours and 3rd dim = samples from the posterior
  # inter_mat[ , , 1] should return a matrix consisting of one sample for every interaction 
  # apply(inter_mat, c(1, 2), mean) will return the mean estimate for every interaction (NB: this is the 
  # mean of the 80% posterior interval, so will be slightly different to the mean value returned from 
  # summary(fit), which is calculated from the full posterior distribution) 
  
}

# This functions spits out diagnostics for convergence
stan_diagnostic <- function(fit, 
                            results_folder,
                            ...) {
  
  # Print diagnostics
  check_hmc_diagnostics(fit)
  
  check_treedepth(fit)
  # if fit saturates the threshold, rerun with larger maximum tree depth and recheck saturation
  check_energy(fit) # check E-BFMI
  check_divergences(fit) # check divergence (validity concern)
  
  fit_summary <- summary(fit)$summary
  # get the range of rhat values
  print(c('rhat range: ', range(na.omit(fit_summary[ , 'Rhat']))))
  # and n_eff
  print(c('n_eff range: ', range(na.omit(fit_summary[ , 'n_eff']))))
  
  # print some diagnostic plots 
  png(paste0(results_folder, '/diag_rhat.png'), width=800, height=800)
  x <- stan_rhat(fit)  # distribution of Rhat
  print(x)
  dev.off()
  png(paste0(results_folder, '/diag_effSS_totalSS.png'), width=800, height=800)
  x <- stan_ess(fit)   # ratio of effective sample size to total sample size
  print(x)
  dev.off()
  png(paste0(results_folder, '/diag_MCse_postsd.png'), width=800, height=800)
  x <- stan_mcse(fit)  # ratio of Monte Carlo standard error to posterior standard deviation
  print(x)
  dev.off() 
  
  #Histograms of lp__ and accept_stat, and their joint distribution
  png(paste0(results_folder, '/diag_sample.png'), width=800, height=800)
  stan_diag(fit, 'sample')
  dev.off()
  # Violin plots showing the distributions of lp__ and accept_stat at each of the sampled step sizes (one per chain)
  png(paste0(results_folder, '/diag_stepsize.png'), width=800, height=800)
  stan_diag(fit, 'stepsize')
  dev.off()
  # Histogram of treedepth and violin plots showing the distributions of lp__ and accept_stat for each value of treedepth
  png(paste0(results_folder, '/diag_treedepth.png'), width=800, height=800)
  stan_diag(fit, 'treedepth')
  dev.off()
  # Violin plots showing the distributions of lp__ and accept_stat for iterations that encountered divergent transitions (divergent=1) and those that did not (divergent=0)
  png(paste0(results_folder, '/diag_divergence.png'), width=800, height=800)
  stan_diag(fit, 'divergence')
  dev.off()
  
  # NB: for further diagnostics, I can explore with
  # - stan_par(fit, 'specific parameter')
  # - stan_ac(fit, 'param') for auto-correlation
}

# This function checks the validity of the RE model
stan_model_check <- function(fit,
                             results_folder,
                             params = c('disp_dev', 'a', 'effect', 'response'),
                             ...){
  
  fit_summary <- summary(fit)$summary
  
  # remove 'mu' from the parameter vector and do it after because it is too long
  params <- params[params != 'mu']
  
  sapply(params, function(x) { 
    
    # this is just to length the plots so they can show all parameters
    N <- length(grep(x, rownames(fit_summary)))
    # some exceptions: 
    if (x == 'a') {N <- 20}
    if (x == 're') {N <- 800}
    if (x == 'sigma_alph') {N <- 10}
    
    # save traceplot of the posterior draws
    png(paste0(results_folder, '/tplot_', x, '.png'), width = 800, height = (N/5*70))
    tplot <- rstan::traceplot(fit, pars = x, inc_warmup = F, ncol = 5)
    print(tplot)
    dev.off()
    
    # plot posterior uncertainty intervals
    png(paste0(results_folder, '/postui_', x, '.png'), width = 600, height = (N/5*70))
    postint <- plot(fit, pars = x, show_density = T) 
    print(postint)
    dev.off()
  })
}

# This function plots the posterior predicted seed number vs the actual data
stan_post_pred_check <- function(post.draws,
                                 results_folder,
                                 stan.data,
                                 ...) {
  
  # currently using the loo package, can switch to rethinking
  
  # phi is the overdispersion parameter for the neg binom model
  # mu is the mean for predicted seed number 
  
  # extract mu and phi
  mu <- post.draws$mu # matrix with nrow = draws and ncol = observations
  disp_dev <- post.draws$disp_dev
  phi <- (disp_dev^2)^(-1)
  
  # generating posterior predictions
  seed_pred <- matrix(nrow = dim(mu)[1], ncol = dim(mu)[2])
  for (i in 1:dim(mu)[1]) {     # for each posterior draw
    for (j in 1:dim(mu)[2]) {    # for each observation 
      # draw from the predicted distribution
      seed_pred[i, j] <- rnbinom(1, mu = mu[i, j], size = phi[i])  
    }
  }
  # get maximum density for plot limits
  max.density <- max(c(apply(seed_pred, 1, function(x) {max(density(x)$y)}), 
                       max(density(stan.data$seeds)$y)))
  
  # dev.new(noRStudioGD = T)
  png(paste0(results_folder, '/postpredch.png'), width=800, height=800)
  # start a plot with the first draw 
  ppc.plot <- plot(density(seed_pred[1, ]), ylim = c(0, max.density), col = 'darkgrey',
                   ylab = 'Seed density',
                   main = 'Post. pred. check',
                   sub = '(grey = predicted, black = observed)') 
  for (i in 2:dim(seed_pred)[1]) {
    # add a line for each draw
    ppc.plot <- lines(density(seed_pred[i, ]), col = 'darkgrey')
  }
  # add the actual data
  ppc.plot <- lines(density(stan.data$seeds), col = 'black', lwd = 2)  
  print(ppc.plot)
  dev.off()
}
