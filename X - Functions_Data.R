#' ####################################################################### #
#' PROJECT: [PhD; X - DATA FUNCTIONALITY] 
#' CONTENTS: 
#'  - Functionality for retrieval of environmental data for observations
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
  Avg_phylo <- phytools::averageTree(Phylo_ls)
  ## calculate distance aggregates
  PhyloDist_mean <- apply(simplify2array(PhyloDist_ls), 1:2, mean) # build mean of all distances
  PhyloDist_SD <- apply(simplify2array(PhyloDist_ls), 1:2, sd) # build sd of all distances
  ## combine it into one object
  Phylo_ls <- list(Avg_Phylo = Avg_phylo,
                   Dist_Mean = PhyloDist_mean,
                   Dist_SD = PhyloDist_SD)
  ## Export of Results -----
  return(Phylo_ls) 
}