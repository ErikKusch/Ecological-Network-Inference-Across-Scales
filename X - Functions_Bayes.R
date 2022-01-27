#' ####################################################################### #
#' PROJECT: [PhD; X - BAYESIAN FUNCTIONALITY] 
#' CONTENTS: 
#'  - Functionality for execution and inspection of Bayesian methodology/models
#'  DEPENDENCIES:
#'  - 
#' AUTHOR: [Malyon Bimler, Erik Kusch]
#' ####################################################################### #

# HMSC =====================================================================
HMSC.Eval <- function(Model = NULL, Name = NULL, Dir = getwd(), thin = thin, nSamples = nSamples, nChains = nChains){
  ma <- NULL
  na <- NULL
  ### MCMC convergence ----
  mpost <- convertToCodaObject(Model, spNamesNumbers = c(T,F), covNamesNumbers = c(T,F))
  psrf.beta <- gelman.diag(mpost$Beta,multivariate=FALSE)$psrf # "if these are smaller than 1.05, we are quite happy"
  tmp <- summary(psrf.beta)
  if(is.null(ma)){
    ma=psrf.beta[,1]
    na = paste0(as.character(thin),",",as.character(nSamples))
  } else {
    ma = cbind(ma,psrf.beta[,1])
    if(j==1){
      na = c(na,paste0(as.character(thin),",",as.character(nSamples)))
    } else {
      na = c(na,"")
    }
  }
  jpeg(file=file.path(Dir, paste0(Name, "_MCMC_convergence.jpeg")), width = 32, height = 32, units = "cm", res = 100)
  par(mfrow=c(2,1))
  vioplot(ma,col=rainbow_hcl(1),names=na,ylim=c(0,max(ma)),main="psrf(beta)")
  vioplot(ma,col=rainbow_hcl(1),names=na,ylim=c(0.9,1.1),main="psrf(beta)")
  dev.off()
  
  ### Model Fit ----
  m <- Model
  preds <- tryCatch(computePredictedValues(m)) # this produces "first error: the leading minor of order 2 is not positive definite" in some instances, which has been linked to Inf values on the github issues page for HMSC, but aren't present in my models
  MF <- evaluateModelFit(hM=m, predY=preds)
  # partition <- createPartition(m, nfolds = 2)
  # preds2 <- computePredictedValues(m,partition=partition, nParallel = nChains)  # this produces "first error: the leading minor of order 2 is not positive definite" in some instances, which has been linked to Inf values on the github issues page for HMSC, but aren't present in my models
  # MFCV <- evaluateModelFit(hM=m, predY=preds2)
  # WAIC <- computeWAIC(m)
  # filename_out <- file.path(Dir, paste0(Name, "_MF_thin_", as.character(thin),
  #                                       "_samples_", as.character(nSamples),
  #                                       "_chains_",as.character(nChains),
  #                                       ".Rdata"))
  # save(MF,MFCV,WAIC,Name,file = filename_out)
  # jpeg(file=file.path(Dir, paste0(Name, "_model_fit.jpeg")), width = 32, height = 32, units = "cm", res = 100)
  # if(!is.null(MF$TjurR2)){
  #   plot(MF$TjurR2,MFCV$TjurR2,xlim=c(-1,1),ylim=c(-1,1),
  #        xlab = "explanatory power",
  #        ylab = "predictive power",
  #        main=paste0(Name,", thin = ",as.character(thin),", samples = ",as.character(nSamples),": Tjur R2"))
  #   abline(0,1)
  #   abline(v=0)
  #   abline(h=0)
  # }
  # if(!is.null(MF$R2)){
  #   plot(MF$R2,MFCV$R2,xlim=c(-1,1),ylim=c(-1,1),
  #        xlab = "explanatory power",
  #        ylab = "predictive power",
  #        main=paste0(Name,", thin = ",as.character(thin),", samples = ",as.character(nSamples),": R2"))
  #   abline(0,1)
  #   abline(v=0)
  #   abline(h=0)
  # }
  # if(!is.null(MF$AUC)){
  #   plot(MF$AUC,MFCV$AUC,xlim=c(0,1),ylim=c(0,1),
  #        xlab = "explanatory power",
  #        ylab = "predictive power",
  #        main=paste0(Name,", thin = ",as.character(thin),", samples = ",as.character(nSamples),": AUC"))
  #   abline(0,1)
  #   abline(v=0.5)
  #   abline(h=0.5)
  # }
  # dev.off()
  
  ### PARAMETER ESTIMATES ----
  pdf(file = file.path(Dir, paste0(Name, "_parameter_estimates.pdf")))
  VP <- computeVariancePartitioning(Model)
  vals <- VP$vals
  mycols <- rainbow(nrow(VP$vals))
  plotVariancePartitioning(hM=Model, VP=VP,cols = mycols, args.leg=list(bg="white",cex=0.7),
                           main = paste0("Proportion of explained variance, ",Name),cex.main=0.8)
  R2 <- NULL
  if(!is.null(MF$TjurR2)){
    TjurR2 <- MF$TjurR2
    vals <- rbind(vals,TjurR2)
    R2 <- TjurR2
  }
  if(!is.null(MF$R2)){
    R2 <- MF$R2
    vals <- rbind(vals,R2)
  }
  write.csv(vals,file=file.path(Dir, paste0(Name, "_parameter_estimates_VP.csv")))
  if(!is.null(R2)){
    VPr <- VP
    for(k in 1:m$ns){
      VPr$vals[,k] <- R2[k]*VPr$vals[,k]
    }
    VPr$vals <- VPr$vals[,order(-R2)]
    plotVariancePartitioning(hM=Model, VP=VPr,cols = mycols, args.leg=list(bg="white",cex=0.7),ylim=c(0,1),
                             main=paste0("Proportion of raw variance, ",Name),cex.main=0.8)
  }
  postBeta = getPostEstimate(Model, parName="Beta")
  show.sp.names = (is.null(Model$phyloTree) && Model$ns<=20) 
  plotBeta(Model, post=postBeta, supportLevel = 0.95,param="Sign",
           plotTree = !is.null(Model$phyloTree),
           covNamesNumbers = c(TRUE,FALSE),
           spNamesNumbers=c(show.sp.names,FALSE),
           cex=c(0.6,0.6,0.8))
  mymain = paste0("BetaPlot, ",Name)
  if(!is.null(Model$phyloTree)){
    mpost = convertToCodaObject(Model)
    rhovals = unlist(poolMcmcChains(mpost$Rho))
    mymain = paste0(mymain,", E[rho] = ",round(mean(rhovals),2),", Pr[rho>0] = ",round(mean(rhovals>0),2))
  }
  title(main=mymain, line=2.5, cex.main=0.8)
  me = as.data.frame(t(postBeta$mean))
  me = cbind(Model$spNames,me)
  colnames(me) = c("Species",Model$covNames)
  po = as.data.frame(t(postBeta$support))
  po = cbind(Model$spNames,po)
  colnames(po) = c("Species",Model$covNames)
  ne = as.data.frame(t(postBeta$supportNeg))
  ne = cbind(Model$spNames,ne)
  colnames(ne) = c("Species",Model$covNames)
  vals = list("Posterior mean"=me,"Pr(x>0)"=po,"Pr(x<0)"=ne)
  writexl::write_xlsx(vals,path = file.path(Dir, paste0(Name, "_parameter_estimates_Beta.xlsx")))
  if(Model$nt>1){
    postGamma = getPostEstimate(Model, parName="Gamma")
    plotGamma(Model, post=postGamma, supportLevel = 0.9, param="Sign",
              covNamesNumbers = c(TRUE,FALSE),
              trNamesNumbers=c(m$nt<21,FALSE),
              cex=c(0.6,0.6,0.8))
    title(main=paste0("GammaPlot ",Name), line=2.5,cex.main=0.8)
  }
  OmegaCor = computeAssociations(Model)
  supportLevel = 0.95
  for (r in 1:Model$nr){
    plotOrder = corrMatOrder(OmegaCor[[r]]$mean,order="AOE")
    toPlot = ((OmegaCor[[r]]$support>supportLevel) + (OmegaCor[[r]]$support<(1-supportLevel))>0)*sign(OmegaCor[[r]]$mean)
    if(Model$ns>20){
      colnames(toPlot)=rep("",Model$ns)
      rownames(toPlot)=rep("",Model$ns)
    }
    mymain = paste0("Associations, ",Name, ": ",names(Model$ranLevels)[[r]])
    if(Model$ranLevels[[r]]$sDim>0){
      mpost = convertToCodaObject(Model)
      alphavals = unlist(poolMcmcChains(mpost$Alpha[[1]][,1]))
      mymain = paste0(mymain,", E[alpha1] = ",round(mean(alphavals),2),", Pr[alpha1>0] = ",round(mean(alphavals>0),2))
    }
    corrplot(toPlot[plotOrder,plotOrder], method = "color",
             col=colorRampPalette(c("blue","white","red"))(3),
             mar=c(0,0,1,0),
             main=mymain,cex.main=0.8)
    
    me = as.data.frame(OmegaCor[[r]]$mean)
    me = cbind(Model$spNames,me)
    colnames(me)[1] = ""
    po = as.data.frame(OmegaCor[[r]]$support)
    po = cbind(Model$spNames,po)
    colnames(po)[1] = ""
    ne = as.data.frame(1-OmegaCor[[r]]$support)
    ne = cbind(Model$spNames,ne)
    colnames(ne)[1] = ""
    vals = list("Posterior mean"=me,"Pr(x>0)"=po,"Pr(x<0)"=ne)
    writexl::write_xlsx(vals,path = file.path(Dir, paste0(Name, "_parameter_estimates_Omega_",names(Model$ranLevels)[[r]],".xlsx")))
  }
  dev.off()
  return(vals)
}

# IF-REM ===================================================================
## PREPARING DATA ----------------------------------------------------------
# works on a data frame where the first four columns are plotID, fitness proxy (identified with Fitness argument), focal (species membership), focalID (speciesID and plotID)
FUN.StanList <- function(Fitness = "value", data = NULL){
  ## Basic List and Data Cleaning ----
  stan.data <- list() # set up the data in list format as preferred by STAN
  data[ , -c(1:3)] <- sapply(data[ , -c(1:3)], as.numeric) # ensure that neighbourhood protion of the data frame is numeric, first three columns are SiteID, taxon, and fitness proxy
  data$taxon <- as.character(data$taxon)
  
  ## Number of Observations ----
  Counts <- data[ , -c(1:3)] # remove first four columns (SiteID, taxon, and fitness proxy)
  Counts[Counts > 0] <- 1 # set all counts/abundances which indicate presence to 1
  Counts <- split(Counts, data$taxon) # split into groups by species
  obs <- do.call(rbind, lapply(Counts, colSums))
  
  ## Inferable Interactions ---
  # df is the dataframe, the 'focal' column gives you the focal species name for each observation,
  # and nonNcols (3) is a variable which tells you how many columns in df are NOT neighbour abundances
  Q <- t(sapply(levels(as.factor(data$taxon)), function(f){
    # inferrable interactions are identified for each focal species independently
    N_i <- as.matrix(data[data$taxon == f, -c(1:3)]) # I just need the matrix of neighbour
    # abundances for all observations for that one taxon
    X_i <- cbind(1,N_i)
    R_i <- pracma::rref(X_i)
    Z_i <- t(R_i) %*% R_i
    # param k is inferrable if its corresponding row/column is all 0 except for the k'th element
    # ignore intercept because we always want to include it (beta_i0)
    sapply(seq(2, dim(Z_i)[1], 1), function(k){
      ifelse(Z_i[k, k] == 1 & sum(Z_i[k, -k]) == 0, 1, 0)
    })
  }))
  stan.data <- list() # storing data in the format preffered by stan
  stan.data$Q <- Q # Q is a matrix of taxon x neighbours, if Q[i, j] = 1 then the interaction between i and # j is inferrable
  
  ## Integers ----
  stan.data$S <- length(unique(data$taxon))  # number of species
  stan.data$N <- nrow(data)                  # number of observations
  stan.data$K <- ncol(data[ , -c(1:3)])      # number of neighbours
  stan.data$I <- sum(Q)                      # # number of interactions which can be inferred by a single parameter
  stan.data$Z <- stan.data$S*stan.data$K   # total number of possible interactions
  
  ## Vectors ----
  stan.data$species_ID <- as.numeric(as.factor(data$taxon)) # numeric indices for taxon species
  stan.data$fitness <- data[, Fitness] # fitness proxy
  
  ## Indices for Interactions in Alpha-Matrix ----
  stan.data$inter_per_species <- rowSums(Q)
  # column index in the interactions matrix for each inferrable interaction
  stan.data$icol <- unlist(apply(ifelse(Q > 0, T, F), 1, which))
  names(stan.data$icol) <- NULL
  stan.data$icol <- as.vector(stan.data$icol)
  # begin the row index
  stan.data$irow <- rep(1, stan.data$inter_per_species[[1]])
  # begin the start and end indices for the vector of realised interactions per species
  stan.data$istart <- 1
  stan.data$iend <- stan.data$inter_per_species[[1]] #???
  # populate indices for all the other species
  for(s in 2:stan.data$S) {
    # starting position of beta_ij's for i in the vector of observed interactions (ie the 1st 'j')
    stan.data$istart[s] <- sum(stan.data$inter_per_species[s-1:s])+1
    # end position of beta_ij's for i in the vector of observed interactions (the last 'j')
    stan.data$iend[s] <- sum(stan.data$inter_per_species[1:s])
    # row index in the interaction matrix for each observed interaction
    stan.data$irow <- c(stan.data$irow, rep(s, stan.data$inter_per_species[[s]]))
  }
  
  ## Model Matrix ----
  stan.data$X <- as.matrix(data[ , -c(1:3)])
  
  ## Number of observations per interaction ----
  stan.data$Obs <- as.vector(apply(obs, 1, c)) # vector of the number of observations for each interactions
  stan.data$Obs <- stan.data$Obs[stan.data$Obs > 0] # remove unobserved interactions
  
  return(stan.data)
}

## CHECKING DATA -----------------------------------------------------------
FUN.DataDims <- function(data = NULL){
  message(paste0('Number of observations = ', data$N))
  message(paste0('Number of focal species = ', data$S))
  message(paste0('Number of neighbouring species = ', data$K))
  message(paste0('Number of inferable/realised interactions = ', data$I))
}

## MODEL SPECIFICATION -----------------------------------------------------
# model specified in 'joint_model_updated.stan'

## MODEL EVALUATION --------------------------------------------------------
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
    if (length(inter > 0)) {sample(inter, size = dim(response)[1])} 
    else {rep(0, 1000)} # this is for those unobserved interactions (0)
  })
  ifm_inters <- do.call(cbind, ifm_inters) # this returns IF model estimates for every possible interaction 
  # every column is an interaction
  
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
