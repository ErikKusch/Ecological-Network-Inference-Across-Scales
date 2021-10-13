#' ####################################################################### #
#' PROJECT: [PhD; 2B - PLOT NETWORKS] 
#' CONTENTS: 
#'  - Generate plot data base species-association/-interaction networks
#'  DEPENDENCIES:
#'  - 0 - Preamble.R
#'  - 0 - ShapeFiles.R
#'  - X - Functions_Bayes.R
#'  - X - Functions_Plotting.R
#'  - No data dependencies
#' AUTHOR: [Erik Kusch]
#' ####################################################################### #
rm(list=ls())

# PREAMBLE ================================================================
source("0 - Preamble.R")
source("0 - ShapeFiles.R")
source("X - Functions_Plotting.R")
source("X - Functions_Bayes.R")
nSamples <- 7000
nWarmup <- 700
nChains <- 4

####### FIA DATA RETRIEVAL ---------------------------------------------------------
FUN_PlotData_FIA <- function(states = c("DE","MD"), nCores = parallel::detectCores()/2){
  ### EXISTENCE CHECK
  Check_vec <- states %nin% substring(list.files(Dir.Plots.FIA), 1, 2)
  if(length(substring(list.files(Dir.Plots.FIA), 1, 2)) != 0){
    if(unique(substring(list.files(Dir.Plots.FIA), 1, 2) %in% states) != TRUE){stop("Your FIA directory contains data for more than the states you asked to analyse here. Please remove these or change the state argument here to include these files already present.")}
  }
  # might need to run devtools::install_github('hunter-stanke/rFIA') to circumvent "Error in rbindlist(inTables..." as per https://github.com/hunter-stanke/rFIA/issues/7
  if(sum(Check_vec) != 0){
    FIA_df <- rFIA::getFIA(states = states[Check_vec], dir = Dir.Plots.FIA, nCores = nCores) # download FIA state data and save it to the FIA directory
  }else{
    FIA_df <- rFIA::readFIA(dir = Dir.Plots.FIA) # load all of the data in the FIA directory
  }
  
  ## BIOME SHAPE PREPARATION
  FIAMerged_shp <- aggregate(FIA_shp, by = "BIOME") # Merge shapes according to biome type
  FIAMerged_shp@data$Names <- Full_Biomes[match(FIAMerged_shp@data$BIOME, Abbr_Biomes)] # Assign corresponding full text biome names
  if(!file.exists(paste0(Dir.PlotNets.FIA, "/FIA_SitesBiomes.html"))){
    FIAPlots_df <- unique(FIA_df$PLOT[, c("LON", "LAT")]) # plot extraction with lat and lon
    FIAPlots_df <- na.omit(FIAPlots_df) # remove any NA rows
    ## plot biomes and FIA sites as mapview object
    Biomes_mv <- mapview(FIAMerged_shp, color = "black", layer.name = "Biomes", zcol = "Names", col.regions = heat.colors)
    FIA_mv <- mapview(FIAPlots_df, xcol = "LON", ycol = "LAT", legend = FALSE, cex = 3.5, layer.name = "FIA Sites", color = "Black", grid = FALSE)
    save(Biomes_mv, FIA_mv, FIAPlots_df, FIAMerged_shp, file = file.path(Dir.Plots, "MapData_FIA.RData"))
    ## may need to run remotes::install_github("r-spatial/mapview") to get this to work
    mapshot(Biomes_mv+FIA_mv, url = paste0(Dir.PlotNets.FIA, "/FIA_SitesBiomes.html"), fgb = FALSE)
  }
  
  ## CALCULATION OF FITNESS AS APPROXIMATED BY BIOMASS
  if(!file.exists(file.path(Dir.Plots, "FIABiomes_df.rds"))){
    FIABiomass_df <- biomass(db = FIA_df, # which data base to use
                             polys = FIAMerged_shp, 
                             bySpecies = TRUE, # group by Species
                             byPlot = TRUE, # group by plot
                             nCores = nCores, # use half the machines cores for this
                             treeType = "live" 
    )
    saveRDS(FIABiomass_df, file.path(Dir.Plots, "FIABiomes_df.rds"))
  }else{
    load(file.path(Dir.Plots, "FIABiomes_df.rds"))
  }
  
  ## SPLITTING INTO BIOMES
  FIASplit_ls <- split(FIABiomass_df, FIABiomass_df$polyID) # extract for each biome to separate data frame
  names(FIASplit_ls) <- FIAMerged_shp@data$Names # apply correct names of biomes
  FIASplit_ls <- FIASplit_ls[-c((length(FIASplit_ls)-1):length(FIASplit_ls))] # remove 98 and 99 biome which is barren and limnic
  for(Biome_Iter in 1:length(FIASplit_ls)){
    BiomeName <- names(FIASplit_ls)[Biome_Iter]
    print(BiomeName)
    if(file.exists(file.path(Dir.Plots, paste0("FIABiome", Biome_Iter, ".RData")))){next()}
    FIAIter_df <- FIASplit_ls[[Biome_Iter]]
    FIAIter_df <- FIAIter_df[,c("pltID", "BIO_ACRE", "SCIENTIFIC_NAME", "nStems", "YEAR")] # select columns we need
    colnames(FIAIter_df) <- c("plot", "biomass", "focal", "Number at Plot", "Year") # assign new column names
    Species_vec <- unique(FIAIter_df$focal) # identify all species names
    Plots_vec <- unique(FIAIter_df$plot) # identify all plot IDs
    Interaction_ls <- as.list(rep(NA, length = length(Plots_vec))) # establish empty list with one slot for each plot
    names(Interaction_ls) <- Plots_vec # set name of the list positions to reflect the plotIDs
    counter <- 1 # create counter to index where to put the data frame created in the loop into the list
    for(Iter_plot in Plots_vec){ # plot loop: loop over all plots
      Iter_df <- FIAIter_df[FIAIter_df$plot == Iter_plot, ] # select data for currently focussed plot
      Iter_df <- Iter_df[order(Iter_df$Year, Iter_df$focal), ] # sort by year first and then by focal alphabetically in each year
      # Iter_df <- group_by(.data = Iter_df, .dots=c("focal", "Year")) %>% # group by species
      #   summarise_at(.vars = c("biomass", "Number at Plot"), .funs = median) # summarise biomass and number grouped by species
      
      Counts_mat <- matrix(rep(0, length = dim(Iter_df)[1]*length(Species_vec)), nrow = dim(Iter_df)[1]) # create empty matrix
      colnames(Counts_mat) <- Species_vec # assign column names of species names
      for(k in unique(Iter_df$Year)){
        k_df <- Iter_df[Iter_df$Year == k,] # select data for year at that plot
        Matches_vec <- base::match(x = k_df$focal, table = Species_vec) # identify position of matches of species names
        Counts_k <- rep(Iter_df$`Number at Plot`[Iter_df$Year == k], # repeat counts of each species
                        each = sum(Iter_df$Year == k) # as often as there are observations at that plot in that year
        )
        Counts_mat[Iter_df$Year == k,Matches_vec] <- Counts_k # save counts to right positions in neighbour matrix
      }
      Interaction_ls[[counter]] <- cbind(as.data.frame(Iter_df), as.data.frame(Counts_mat)) # combine matrix with biomass data
      counter <- counter +1 # raise counter
    }
    Interaction_df <- bind_rows(Interaction_ls, .id = "plot") # combine data frames in list elements into one big data frame
    Interaction_df <- cbind(Interaction_df$plot, Interaction_df$biomass, Interaction_df$focal, Interaction_df$Year, Interaction_df[,6:dim(Interaction_df)[2]])
    colnames(Interaction_df)[1:4] <- c("plot", "biomass", "focal", "year")
    # Interaction_df <- Interaction_df[-c(which(Interaction_df$biomass == 0)), ] # remove 0 biomass entries
    print(paste("Data dimensions:", paste(dim(Interaction_df), collapse = " & ")))
    save(BiomeName, FIAIter_df, Interaction_df, file = file.path(Dir.Plots, paste0("FIABiome", Biome_Iter, ".RData")))
  }
}

if(sum(file.exists(file.path(Dir.Plots, paste0("FIABiome", 1:12, ".RData")))) != 12){
  FUN_PlotData_FIA(states = c("AL", "AK", "AZ", "AR", "CA", "CO", "CT", "DE", "FL", "GA", "HI", "ID", "IL", "IN", "IA", "KS", "KY", "LA", "ME", "MD", "MA", "MI", "MN", "MS", "MO", "MT", "NE", "NV", "NH", "NJ", "NM", "NY", "NC", "ND", "OH", "OK", "OR", "PA", "RI", "SC", "SD", "TN", "TX", "UT", "VT", "VA", "WA", "WV", "WI", "WY"), nCores = parallel::detectCores()/2)
}

################# RUN NETWORK METHOD FOR EACH FIA SUBSET ########################################### --------------------------

# probably only want to run this for 3, 6, 7, 10, 11 ==> c(3,6,7,11)
# Model_Iter = 11 # or 7&10&11 for high data
FIABiomes_fs <- list.files(path = Dir.Plots, pattern = "FIABiome")


for(Model_Iter in 1:length(FIABiomes_fs)){ # 10 for biggest data set
  load(file.path(Dir.Plots, FIABiomes_fs[[Model_Iter]]))
  print(BiomeName)
  print(dim(Interaction_df))
  
  Dir.FIABiome <- file.path(Dir.PlotNets.FIA, FIABiomes_fs[[Model_Iter]])
  if(dir.exists(Dir.FIABiome)){
    warning(paste(FIABiomes_fs[[Model_Iter]], "directory already exists. Skipping."))
    next()
  }
  dir.create(Dir.FIABiome)
  
  Test_df <- Interaction_df[Interaction_df$biomass != 0, ]
  
  tmp <- Test_df[, 1:3]
  tmp$focalID <- with(Test_df, paste(plot, focal, year, sep ="_"))
  Test_df <- cbind(tmp, Test_df[, 4:ncol(Test_df)])
  colnames(Test_df)[1:5] <- c("plotID", "fit", "focal", "focalID", "time")
  
  
  Mal_df <- Test_df[,-5] # data frame without time
  
  OmitCols <- -(which(colnames(Mal_df)[-c(1:4)] %nin% Mal_df$focal)+4)
  if(length(OmitCols) != 0){
    Mal_df <- Mal_df[, OmitCols]
  }else{
    Mal_df <- Mal_df
  }
  
  FIA_StanList <- Fun_StanList(Fitness = "fit", data = Mal_df)
  
  #### Preferences
  rstan_options(auto_write = TRUE)
  rstan_options(javascript = FALSE)
  options(mc.cores = nChains) 
  
  
  ## CHECKING DATA
  focalID <- unique(Mal_df$focal)
  neighbourID <- colnames(Mal_df[ , -c(1:4)])
  Fun_PreCheck(data = Mal_df)
  ## RUNNING MODEL
  fit <- stan(file = 'Supplement - StanModel.stan',
              data =  FIA_StanList,               # named list of data
              chains = nChains,
              warmup = nWarmup,          # number of warmup iterations per chain
              iter = nSamples,            # total number of iterations per chain
              refresh = 100,         # show progress every 'refresh' iterations
              control = list(max_treedepth = 10)
  )
  save(fit, file = file.path(Dir.FIABiome, "Model.RData"))
  
  ## MODEL DIAGNOSTICS
  # Get the full posteriors 
  joint.post.draws <- extract.samples(fit)
  # Select parameters of interest
  param.vec <- c('a', 'beta_ij', 'effect', 'response', 're', 'inter_mat', 'mu', 'sigma_alph'
                 # 'disp_dev' # including this one leads to "Error in apply(joint.post.draws[[p]], 2, function(x) { : dim(X) must have a positive length"
  )
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
  colnames(intrinsic.perf) <- focalID
  inter_mat <- return_inter_array(joint.post.draws, 
                                  response = p.samples$response,
                                  effect = p.samples$effect,
                                  sort(focalID),
                                  neighbourID)
  # inter_mat is now a 3 dimensional array, where rows = focals, columns = neighbours and 3rd dim = samples from the posterior; inter_mat[ , , 1] should return a matrix consisting of one sample for every interaction 
  param.vec <- c('a', 'beta_ij', 'effect', 'response', 're', 'inter_mat', 'mu', 'sigma')
  try(stan_model_check(fit = fit,
                       results_folder = Dir.FIABiome,
                       params = param.vec))
  # stan_post_pred_check(post.draws = joint.post.draws, # error here on missing data!!!
  #                       results_folder = Dir.Malyon,
  #                       stan.data = StanList_ls)
  # Interactions can now be divided by the appropriate scaling (intrinsic performance, and demographic rates
  # if a population dynamics model is used) in order to return per capita interaction strengths. 
  Interaction_mean <- apply(inter_mat, c(1, 2), mean) # will return the mean estimate for every interaction (NB: this is the mean of the 80% posterior interval, so will be slightly different to the mean value returned from summary(fit), which is calculated from the full posterior distribution)  
  Interaction_mean <- Interaction_mean[, order(colnames(Interaction_mean))] # sort columns
  Interaction_mean <- Interaction_mean*-1 # need to switch sign of results
  diag(Interaction_mean) <- NA
  Interaction_hpdi <- apply(inter_mat, c(1, 2), HPDI, prob = 0.89) *-1 # need to switch sign of results
  Interaction_min <- Interaction_hpdi[1,,]
  diag(Interaction_min) <- NA
  Interaction_max <- Interaction_hpdi[2,,]
  diag(Interaction_max) <- NA
  Interactions_igraph <- data.frame(Actor = rep(dimnames(Interaction_mean)$neighbour, length(dimnames(Interaction_mean)$species)),
                                    Subject = rep(dimnames(Interaction_mean)$species, each = length(dimnames(Interaction_mean)$neighbour)),
                                    Inter_mean = as.vector(t(Interaction_mean)),
                                    Inter_min = as.vector(t(Interaction_min)),
                                    Inter_max = as.vector(t(Interaction_max))
  )
  # Interactions_Malyon <- Interactions_igraph[order(abs(Interactions_igraph$Inter_mean), decreasing = TRUE), ]
  # Interactions_igraph <- na.omit(Interactions_Malyon)
  
  Interactions_DAG <- Interactions_igraph[order(Interactions_igraph$Actor), ]
  Interactions_DAG <- na.omit(Interactions_DAG)
  
  Signs_DAG <- with(Interactions_DAG, data.frame(
    Intrer_mean = sign(Inter_mean),
    Intrer_max = sign(Inter_max),
    Intrer_min = sign(Inter_min)
  ))
  
  Interactions_DAG <- Interactions_DAG[abs(apply(Signs_DAG, 1, sum)) == 3, ]
  
  Interactions_DAG$Actor <- gsub(x = Interactions_DAG$Actor, pattern = " ", replacement = "_")
  Interactions_DAG$Subject <- gsub(x = Interactions_DAG$Subject, pattern = " ", replacement = "_")
  
  
  ## Graph Plot
  Dag_Paths <- paste(paste(Interactions_DAG$Actor, "->", Interactions_DAG$Subject), collapse = " ")
  dag <- dagitty(x = paste0("dag {", Dag_Paths, "}"))
  tidy_dag <- tidy_dagitty(dag, layout = "fr")
  tidy_dag$data$weight <- c(abs(Interactions_DAG$Inter_mean), rep(NA, sum(is.na(tidy_dag$data$direction))))
  tidy_dag$data$label <- c(round(Interactions_DAG$Inter_mean, 2), rep(NA, sum(is.na(tidy_dag$data$direction))))
  tidy_dag$data$colour <- c(ifelse(c(Interactions_DAG$Inter_mean) > 0, "green", "red"), rep("white", sum(is.na(tidy_dag$data$direction))))
  WeightMod <- 1/max(tidy_dag$data$weight, na.rm = TRUE)
  
  DAG_plot <- ggplot(tidy_dag, aes(x = x, y = y, xend = xend, yend = yend), layout = "linear") +
    geom_dag_point(colour = "grey", shape = 18, size = 10) +
    geom_dag_text(colour = "black", size = 2) +
    geom_dag_edges_arc(aes(edge_colour = colour, edge_width = exp(weight)*WeightMod)) +
    # geom_dag_edges_diagonal(aes(edge_colour = colour, edge_width = weight*WeightMod)) + 
    labs(title = BiomeName, subtitle = "Realised Interaction (Bimler et al. Method") + 
    theme_dag()
  ggsave(DAG_plot, units = "cm", width = 32, height = 18, filename = file.path(Dir.FIABiome, "DAG.png"))
  
  dag_data <- na.omit(tidy_dag$data)
  df <- data.frame(
    A = dag_data$name, 
    B = dag_data$to)
  df.g <- graph.data.frame(d = df, directed = TRUE)
  E(df.g)$weight <- dag_data$weight
  E(df.g)$label <- dag_data$label
  E(df.g)$colour <- dag_data$colour
  
  igraph_plot <- ggraph(df.g, layout = 'linear', circular = TRUE) + 
    # geom_edge_arc(aes(col = colour, width = weight*WeightMod, alpha = ..index..)) + 
    # scale_edge_alpha('Edge direction', guide = 'edge_direction') + 
    geom_edge_arc(aes(col = colour, width = exp(dag_data$weight)*WeightMod), show.legend = FALSE) + 
    scale_edge_color_manual(values = rev(unique(dag_data$colour))) +
    geom_node_point(color = "black", size = 2) + 
    geom_node_text(aes(label = name),  repel = TRUE)+
    coord_fixed() + 
    labs(title = BiomeName, subtitle = "Realised Interaction (Bimler et al. Method)") + 
    theme_graph()
  ggsave(igraph_plot, units = "cm", width = 32, height = 32, filename = file.path(Dir.FIABiome, "IGraph.png"))
  
  diag(Interaction_mean) <- 0
  Interaction_mean[abs(sign(Interaction_mean) + sign(Interaction_min) + sign(Interaction_max)) != 3] <- 0
  
  
  
  fv.colors = colorRampPalette(c("red","white","green")) ## define the color ramp
  colorlut = fv.colors(100)[c(1,seq(50,100,length.out=99))] ## select colors to use
  
  seq(sum(Interaction_mean>0))
  
  jpeg(filename = file.path(Dir.FIABiome, "Interactions.jpeg"), units = "cm", width = 32, height = 18, res = 1000)
  pheatmap(Interaction_mean, display_numbers = T, 
           color = c("red", "white","green"), 
           breaks = c(min(Interaction_mean, na.rm = TRUE), -0.01, 0.01, max(Interaction_mean, na.rm = TRUE)), main = "Actors (Columns) x Subject (Rows)",
           fontsize = 5)
  dev.off()
}
