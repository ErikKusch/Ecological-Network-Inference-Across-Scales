#' ####################################################################### #
#' PROJECT: [PhD; X - NETWORK PLOTTING] 
#' CONTENTS: 
#'  - Functionality for plotting of networks
#'  DEPENDENCIES:
#'  - 
#' AUTHOR: [Erik Kusch]
#' ####################################################################### #

# NETWORK PLOTTING =======================================================
# plotting a network from a three-column data frame
Plot.Network <- function(Plot_ls = NULL){
  ## column identification
  ## plot generation
  Plots_ls <- as.list(rep(NA, length(Plot_ls)))
  names(Plots_ls) <- names(Plot_ls)
  for(Plot_Iter in 1:length(Plots_ls)){
    
    plots_iter <- Plot_ls[[Plot_Iter]]
    
    regplots_ls <- as.list(rep(NA, length(plots_iter)))
    names(regplots_ls) <- names(plots_iter)
    
    for(i in 1:length(plots_iter)){
      Plot_df <- plots_iter[[i]]$data
      name_plot <- names(plots_iter)[[i]]
      if(startsWith(name_plot, "NDD")){
        Plot_df$Sig <- ifelse(as.numeric(as.character(Plot_df$Sig)) < 1, FALSE, TRUE)
      }
      Plot_df[which(is.na(Plot_df[,4])),4] <- FALSE
      Plot_df[,3][Plot_df[,3] == 0] <- NA
      
      directed <- ifelse(names(plots_iter)[[i]] == "NDD_RIM", TRUE, FALSE)
      regplots_ls[[i]] <- qgraph::qgraph(Plot_df[,1:3], # the three-column data frame of "Partner1", "Partner2", and "Effect"
                                         bg = "white",                       
                                         color = "white", # colour of nodes
                                         layout = 'groups', # layout (I have found this one to be consistent across al plots)
                                         esize = 7, # line width of strongest connection/edge
                                         negCol = "#ba00bd", # colour of negative edges
                                         posCol = "#d9db3b", # colour of positive edges
                                         palette = "colorblind", # colour palette
                                         fade = FALSE, # fading effect on edges relative to strength
                                         directed = directed, # directed network or not
                                         DoNotPlot = TRUE,
                                         lty = Plot_df[,4]*-1+2, # significant interactions are solid, non-significant ones are dashed
                                         title = name_plot
      )
    }
    
    Plots_ls[[Plot_Iter]] <- regplots_ls
  }
  
  
  Plotting_df <- expand.grid(names(Plots_ls), names(regplots_ls))
  
  layout(matrix(1:nrow(Plotting_df), 
                nrow = length(names(regplots_ls)), 
                ncol = length(names(Plots_ls)), 
                byrow = TRUE))
  
  for(Plot_Iter in 1:nrow(Plotting_df)){
    Plot_pos1 <- which(names(Plots_ls) == Plotting_df$Var1[Plot_Iter])
    Plot_pos2 <- which(names(regplots_ls) == Plotting_df$Var2[Plot_Iter])
    plot(Plots_ls[[Plot_pos1]][[Plot_pos2]])
  } 
}

# NETWORKS AS MATRICES ===================================================
# plotting of networks belonging to the same method as matrices
Plot.NetMat <- function(method = "COCCUR", data = YFDP_df, TreatmentOrder = c("Pre-Fire", "Post-Fire"), datareturn = FALSE){ ## currently only works for YFDP
  plot_df <- data[, c(1:2, grep(colnames(data), pattern = method))]
  if(method == "IF_REM"){
    plot_df <- plot_df[, -grep(colnames(plot_df), pattern = "IF_REM_Assoc")]
  }
  
  if("Difference" %in% TreatmentOrder){
    if(method == "HMSC"){
      
      EffectCol <- c(grep(colnames(plot_df), pattern = "effects"), ncol(plot_df)+c(1,3,5) )
      SigCol <- c(grep(colnames(plot_df), pattern = "Sig"), ncol(plot_df)+c(2,4,6) )
      
      plot_df$Difference.abundance <- plot_df[, grep(x = colnames(plot_df), pattern = "abundance")[1]] - plot_df[, grep(x = colnames(plot_df), pattern = "abundance")[3]]
      plot_df$Difference.abundance.Sig <- ifelse(rowSums(plot_df[, grep(x = colnames(plot_df), pattern = "abundance")[c(2,4)]]) == 2, TRUE, FALSE)
      
      if(sum(grepl(pattern = "diametre", x = colnames(data))) == 0){
        plot_df$Difference.diametre <- plot_df[, grep(x = colnames(plot_df), pattern = "biomass")[1]] - plot_df[, grep(x = colnames(plot_df), pattern = "biomass")[3]]
        plot_df$Difference.diametre.Sig <- ifelse(rowSums(plot_df[, grep(x = colnames(plot_df), pattern = "biomass")[c(2,4)]]) == 2, TRUE, FALSE)
      }else{
        plot_df$Difference.diametre <- plot_df[, grep(x = colnames(plot_df), pattern = "diametre")[1]] - plot_df[, grep(x = colnames(plot_df), pattern = "diametre")[3]]
        plot_df$Difference.diametre.Sig <- ifelse(rowSums(plot_df[, grep(x = colnames(plot_df), pattern = "diametre")[c(2,4)]]) == 2, TRUE, FALSE)
      }
      
      plot_df$Difference.presence_absence <- plot_df[, grep(x = colnames(plot_df), pattern = "presence_absence")[1]] - plot_df[, grep(x = colnames(plot_df), pattern = "presence_absence")[3]]
      plot_df$Difference.presence_absence.Sig <- ifelse(rowSums(plot_df[, grep(x = colnames(plot_df), pattern = "presence_absence")[c(2,4)]]) == 2, TRUE, FALSE)
      
    }else{
      EffectCol <- c(grep(colnames(plot_df), pattern = "effects"), ncol(plot_df)+1)
      SigCol <- c(grep(colnames(plot_df), pattern = "Sig"), ncol(plot_df)+2)
      plot_df$Difference <- plot_df[,EffectCol[1]] - plot_df[,EffectCol[2]]
      plot_df$Difference.Sig <- ifelse(rowSums(plot_df[,SigCol[1:2]]) == 2, TRUE, FALSE)
    }
  }else{
    EffectCol <- c(grep(colnames(plot_df), pattern = "effects"))
    SigCol <- c(grep(colnames(plot_df), pattern = "Sig"))
  }
  
  effect_df <- reshape(plot_df, 
                     direction = "long",
                     varying = list(names(plot_df)[EffectCol]),
                     v.names = "Value",
                     idvar = c("Partner1", "Partner2"),
                     timevar = "Condition",
                     times = names(plot_df)[EffectCol])
  
  sig_df <- reshape(plot_df, 
                       direction = "long",
                       varying = list(names(plot_df)[SigCol]),
                       v.names = "Sig",
                       idvar = c("Partner1", "Partner2"),
                       timevar = "Condition",
                       times = names(plot_df)[SigCol])
  
  plot_df <- data.frame(Partner1 = effect_df$Partner1,
                        Partner2 = effect_df$Partner2,
                        Value = effect_df$Value,
                        Sig = sig_df$Sig,
                        Condition = effect_df$Condition)
  plot_df$Condition <-  gsub(plot_df$Condition, pattern = ".effects", replacement = "")
  plot_df$Condition <- gsub(plot_df$Condition, pattern = paste0(method, "."), replacement = "")
  
  if(method == "HMSC"){
    plot_df$Model <- unlist(lapply(strsplit(plot_df$Condition, ".", fixed = TRUE), "[[", 2))
    plot_df$Condition <- unlist(lapply(strsplit(plot_df$Condition, ".", fixed = TRUE), "[[", 1))
  }
  
  ## reshuffling data so all values fall into upper diagonal
  if(method != "IF_REM"){
    for(Row_i in 1:nrow(plot_df)){
      if(method == "HMSC"){
        Op_pos <- which(plot_df$Partner2 == plot_df[Row_i, ]$Partner1 &
                          plot_df$Partner1 == plot_df[Row_i, ]$Partner2 &
                          plot_df$Condition == plot_df[Row_i, ]$Condition &
                          plot_df$Model == plot_df[Row_i, ]$Model)
      }else{
        Op_pos <- which(plot_df$Partner2 == plot_df[Row_i, ]$Partner1 &
                          plot_df$Partner1 == plot_df[Row_i, ]$Partner2 &
                          plot_df$Condition == plot_df[Row_i, ]$Condition)
      }
      if(!is.na(plot_df$Value[Op_pos])){
        plot_df$Value[Row_i] <- plot_df$Value[Op_pos]
        plot_df$Value[Op_pos] <- NA
        plot_df$Sig[Row_i] <- plot_df$Sig[Op_pos]
        plot_df$Sig[Op_pos] <- FALSE
      }
    }
  }
  
  levels(plot_df$Partner1) <- sort(levels(plot_df$Partner1))
  levels(plot_df$Partner2) <- sort(levels(plot_df$Partner2))
  
  gplot <- ggplot(plot_df, aes(x = Partner2, y = Partner1, fill = Value))
  gplot <- gplot +
    geom_tile(color = "white",
              lwd = 1.5,
              linetype = 1) + 
    ## colours, look, and legend
    geom_point(aes(shape = Sig), size = 2) +
    scale_shape_manual(values=c(32, 15), na.translate = FALSE, name = "Significance") +  
    guides(shape = FALSE) + 
    coord_fixed() + 
    guides(fill = guide_colourbar(barwidth = 2,
                                  barheight = 15,
                                  title = "Association")) + 
    ## axes
    theme_bw() + 
    theme(axis.text.x=element_text(angle = -40, hjust = 0))
  
  if(all(plot_df$Value[!is.na(plot_df$Value)] > 0)){
    gplot <- gplot + scale_fill_gradient(high = "#5ab4ac")
  }else{
    if(all(plot_df$Value[!is.na(plot_df$Value)] < 0)){
      gplot <- gplot + scale_fill_gradient(high = "#d8b365")
    }else{
      gplot <- gplot + scale_fill_gradient2(low = "#5ab4ac", high = "#d8b365")
    }
  }
  
  if(method != "IF_REM"){
    gplot <- gplot + labs(x = "", y = "")
  }else{
    gplot <- gplot + labs(x = "Subject", y = "Actor")
  }
  
  if(method != "HMSC"){
    gplot <- gplot + facet_wrap(~factor(Condition, levels = TreatmentOrder))
  }else{
    if(sum(grepl(pattern = "diametre", x = colnames(data))) == 0){
      HMSCOrder <- c("presence_absence", "abundance", "biomass")
    }else{
      HMSCOrder <- c("presence_absence", "abundance", "diametre")
    }
    
    gplot <- gplot + facet_wrap(~ factor(Model, levels = HMSCOrder) +
                                  factor(Condition, levels = TreatmentOrder)
    )
  }
  if(datareturn){
    return(plot_df)
  }else{
  return(gplot)
  }
}

# # NETWORKS AS MATRICES ACROSS SCALES =====================================
# # plotting of networks belonging to the same method as matrices across scales
# Plot.NetMat.Cross <- function(data = NULL,
#                               ALLRegionOrder = c("Pre-Fire", "Yosemite", "Temperate Broadleaf & Mixed Forests"),
#                               ALLMethodOrder = c("COCCUR", "NETASSOC", "HMSC", "NDD-RIM")){
#   
#   plot_ls <- as.list(rep(NA, length(ALLRegionOrder)))
#   names(plot_ls) <- ALLRegionOrder
#   
#   for(region in ALLRegionOrder){
#     print(region)
#     plot_df <- data_ls[[region]]
#     plot_df <- plot_df[, -grep(colnames(plot_df), pattern = "_Assoc")]
#     
#     EffectCol <- c(grep(colnames(plot_df), pattern = "effects"))
#     SigCol <- c(grep(colnames(plot_df), pattern = "Sig"))
#     # plot_df$Difference <- plot_df[,EffectCol[1]] - plot_df[,EffectCol[2]]
#     # plot_df$Difference.Sig <- ifelse(rowSums(plot_df[,SigCol[1:2]]) == 2, TRUE, FALSE)
#     
#     effect_df <- reshape(plot_df, 
#                          direction = "long",
#                          varying = list(names(plot_df)[EffectCol]),
#                          v.names = "Value",
#                          idvar = c("Partner1", "Partner2"),
#                          timevar = "Condition",
#                          times = names(plot_df)[EffectCol])
#     
#     sig_df <- reshape(plot_df, 
#                       direction = "long",
#                       varying = list(names(plot_df)[SigCol]),
#                       v.names = "Sig",
#                       idvar = c("Partner1", "Partner2"),
#                       timevar = "Condition",
#                       times = names(plot_df)[SigCol])
#     
#     plot_df <- data.frame(Partner1 = effect_df$Partner1,
#                           Partner2 = effect_df$Partner2,
#                           Value = effect_df$Value,
#                           Sig = sig_df$Sig,
#                           Condition = effect_df$Condition)
#     plot_df$Condition <-  unlist(lapply(strsplit(plot_df$Condition, split = "\\.") , "[[", 1))
#     
#     
#     
#     
#   }
#   
#   
#   
#   
#   
#   
#   
#   
#   
#   
#   
#   
#   
#   
#   
#   
#   
#   
#   
#   
#   
#   
#   
#   
#   
#   
#   
#   
#   
#   
#   
#   
#   
#   
#   for(method in ALLMethodOrder){
#     print(method)
#     
#     
#     method_ls <- lapply(names(data), FUN = function(x){
#       x <- data[[x]]
#       x[, c(1:2, grep(colnames(x), pattern = method))]
#     })
#     
#     do.call(cbind, method_ls)
#     
#     plot_df <- plot_ls[[]]
#     
#     plot_df <- data[, c(1:2, grep(colnames(data), pattern = method))]
#     
#     if("Difference" %in% ALLRegionOrder){
#       if(method == "HMSC"){
#         
#         if(sum(grepl(pattern = "diametre", x = colnames(data))) == 0){
#           plot_df$Difference.diametre <- plot_df[, grep(x = colnames(plot_df), pattern = "biomass")[1]] - plot_df[, grep(x = colnames(plot_df), pattern = "biomass")[3]]
#           plot_df$Difference.diametre.Sig <- ifelse(rowSums(plot_df[, grep(x = colnames(plot_df), pattern = "biomass")[c(2,4)]]) == 2, TRUE, FALSE)
#         }else{
#           plot_df$Difference.diametre <- plot_df[, grep(x = colnames(plot_df), pattern = "diametre")[1]] - plot_df[, grep(x = colnames(plot_df), pattern = "diametre")[3]]
#           plot_df$Difference.diametre.Sig <- ifelse(rowSums(plot_df[, grep(x = colnames(plot_df), pattern = "diametre")[c(2,4)]]) == 2, TRUE, FALSE)
#         }
#         
#       }else{
#         EffectCol <- c(grep(colnames(plot_df), pattern = "effects"), ncol(plot_df)+1)
#         SigCol <- c(grep(colnames(plot_df), pattern = "Sig"), ncol(plot_df)+2)
#         plot_df$Difference <- plot_df[,EffectCol[1]] - plot_df[,EffectCol[2]]
#         plot_df$Difference.Sig <- ifelse(rowSums(plot_df[,SigCol[1:2]]) == 2, TRUE, FALSE)
#       }
#     }else{
#       EffectCol <- c(grep(colnames(plot_df), pattern = "effects"))
#       SigCol <- c(grep(colnames(plot_df), pattern = "Sig"))
#     }
#     
#     effect_df <- reshape(plot_df, 
#                          direction = "long",
#                          varying = list(names(plot_df)[EffectCol]),
#                          v.names = "Value",
#                          idvar = c("Partner1", "Partner2"),
#                          timevar = "Condition",
#                          times = names(plot_df)[EffectCol])
#     
#     sig_df <- reshape(plot_df, 
#                       direction = "long",
#                       varying = list(names(plot_df)[SigCol]),
#                       v.names = "Sig",
#                       idvar = c("Partner1", "Partner2"),
#                       timevar = "Condition",
#                       times = names(plot_df)[SigCol])
#     
#     plot_df <- data.frame(Partner1 = effect_df$Partner1,
#                           Partner2 = effect_df$Partner2,
#                           Value = effect_df$Value,
#                           Sig = sig_df$Sig,
#                           Condition = effect_df$Condition)
#     plot_df$Condition <-  gsub(plot_df$Condition, pattern = ".effects", replacement = "")
#     plot_df$Condition <- gsub(plot_df$Condition, pattern = paste0(method, "."), replacement = "")
#     if(method == "HMSC"){
#       plot_df$Condition <- gsub(plot_df$Condition, pattern = ".biomass", replacement = "")
#       plot_df$Condition <- gsub(plot_df$Condition, pattern = ".diametre", replacement = "")
#     } 
#   
#     if(method == "NETASSOC" | method == "NDD-RIM"){
#       gplot <- ggplot(plot_df, aes(x = Partner2, y = Partner1, fill = Value))
#     }else{
#       gplot <- ggplot(plot_df, aes(x = Partner1, y = Partner2, fill = Value))
#     }
#     gplot <- gplot +
#       geom_tile(color = "white",
#                 lwd = 1.5,
#                 linetype = 1) + 
#       ## colours, look, and legend
#       geom_point(aes(shape = as.factor(Sig)), size = 4) +
#       scale_shape_manual(values=c(32, 15), na.translate = FALSE, name = "Significance") +  
#       guides(shape = FALSE) + 
#       coord_fixed() + 
#       ## axes
#       theme_bw() + 
#       theme(axis.text.x=element_text(angle = -20, hjust = 0),
#             text = element_text(size=15))
#     
#     if(all(plot_df$Value[!is.na(plot_df$Value)] > 0)){
#       gplot <- gplot + scale_fill_gradient(high = "#d9db3b")
#     }else{
#       if(all(plot_df$Value[!is.na(plot_df$Value)] < 0)){
#         gplot <- gplot + scale_fill_gradient(high = "#ba00bd")
#       }else{
#         gplot <- gplot + scale_fill_gradient2(low = "#ba00bd", high = "#d9db3b")
#       }
#     }
#     
#     if(method != "NDD-RIM"){
#       gplot <- gplot + labs(x = "", y = "") + 
#         guides(fill = guide_colourbar(barwidth = 2,
#                                       barheight = 15,
#                                       title = "Association"))
#     }else{
#       gplot <- gplot + labs(y = "Actor", x = "Subject") +
#         guides(fill = guide_colourbar(barwidth = 2,
#                                       barheight = 15,
#                                       title = "Interaction \n Strength")) 
#     }
#     gplot <- gplot + facet_wrap(~factor(Condition, levels = ALLRegionOrder))
#     plot_ls[[method]] <- gplot + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
#     
#   }
#     return(plot_grid(plotlist=plot_ls, ncol = 1, labels = ALLMethodOrder))
#   }

# NETWORK TOPOLOGIES =====================================================
## for calculation and visualisation of topology metrics
Calc.Topology <- function(data_ls = Comparable_ls,
                          Shared_abbr = Shared_abbr,
                          Sig = TRUE){
  
  Plotting_df <- expand.grid(names(data_ls), names(data_ls[[1]]))
  
  # output objects
  Base_df <- data.frame(Species = NA,
                        Value = NA,
                        Model = NA,
                        Treatment = NA)
  output_ls <- list(Centrality = list(Strength = Base_df,
                                      Eigenvector = Base_df),
                    Modularity = Base_df[,-1],
                    adj = c())
  
  # loop over models here
  for(i in 1:nrow(Plotting_df)){
    print(paste("#####", i))
    
    name_Region <- as.character(Plotting_df$Var1[[i]])
    name_Model <- as.character(Plotting_df$Var2[[i]])
    
    data_i <- data_ls[[name_Region]][[name_Model]]
    data_i <- data_i[data_i$Shared == 2, ]
    
    ## network building
    ## only keep singificant interactions when Sig == TRUE
    if(Sig){
      data_i$Sig[is.na(data_i$Sig)] <- FALSE
      if(startsWith(name_Model, "NDD")){
        SigPos <- as.numeric(as.character(data_i$Sig)) > 0.5
      }else{
        SigPos <- data_i$Sig
      }
      
      if(sum(SigPos)>0){
        data_i <- data_i[SigPos,]
      }else{
        data_i <- NA
        
        strength_df <- data.frame(Species = NA,
                                  Value = NA,
                                  Model = NA,
                                  Treatment = NA)
        rownames(strength_df) <- c()
        colnames(strength_df) <- colnames(Base_df)
        output_ls$Centrality$Strength <- rbind(output_ls$Centrality$Strength, strength_df)
        output_ls$Centrality$Eigenvector <- rbind(output_ls$Centrality$Eigenvector, strength_df)
        output_ls$Modularity <- rbind(output_ls$Modularity, c(NA, name_Model, name_Region))
        
        next()
      }
    }else{
      data_i <- data_i
    }
    
    data_i <- na.omit(data_i[,c("Partner1", "Partner2", "Value")])
    colnames(data_i)[3] <- "weight"
    
    
    if(name_Model == "NDD_RIM"){directed=TRUE}else{directed=FALSE}
    
    g <- igraph::graph_from_data_frame(d = data_i, directed = directed)
    
    adj <- igraph::as_adjacency_matrix(g, attr = "weight")
    
    g_mod <- g
    if(!is.null(E(g)$weight)){
      E(g)$weight <- abs(E(g)$weight)
    }
    
    ## Node-Level Topology ----
    ### Strength
    print("Node Strength")
    strength <- range01(sort(strength(g), decreasing = TRUE))
    strength_df <- data.frame(Species = names(strength),
                              Value = strength,
                              Model = name_Model,
                              Treatment = name_Region)
    rownames(strength_df) <- c()
    colnames(strength_df) <- colnames(Base_df)
    output_ls$Centrality$Strength <- rbind(output_ls$Centrality$Strength, strength_df)
    
    ### Eigenvector Centrality
    print("Eigenvectors")
    ecvent <- sort(evcent(g, scale = TRUE)$vector, decreasing = TRUE)
    ecvent_df <- data.frame(Species = names(ecvent),
                            Value = ecvent,
                            Model = name_Model,
                            Treatment = name_Region)
    rownames(ecvent_df) <- c()
    colnames(ecvent_df) <- colnames(Base_df)
    output_ls$Centrality$Eigenvector <- rbind(output_ls$Centrality$Eigenvector, ecvent_df)
    
    ## Network-Level Topology ----
    ### Modularity
    print("Modularity")
    Fun.Modu <- function(g_mod){
      E(g_mod)$weight[is.na(E(g_mod)$weight)] <- 0
      e <- simplify(g_mod)
    
        fc <- try(cluster_spinglass(g_mod, implementation = "neg", spins = 1e2), silent = TRUE)
        if(class(fc) == "try-error"){
          # stop()
          modularity <- NA
          rm(fc)
          # fc <- cluster_optimal(g)
          # modularity <- modularity(e, membership = fc$membership)
        }else{
          modularity <- fc$modularity
        }
      
      return(modularity)
    }
    Modu_infer <- Fun.Modu(g_mod)
    
    output_ls$Modularity <- rbind(output_ls$Modularity, c(Modu_infer, name_Model, name_Region))
    
  }  
  ## cleaning output object
  output_ls$Centrality$Strength <- output_ls$Centrality$Strength[-1,]
  output_ls$Centrality$Strength$Value[which(output_ls$Centrality$Strength$Value == "NaN")] <- 1
  output_ls$Centrality$Eigenvector <- output_ls$Centrality$Eigenvector[-1,]
  output_ls$Centrality$Eigenvector$Value[which(output_ls$Centrality$Eigenvector$Value == "NaN")] <- 1
  output_ls$Modularity <- output_ls$Modularity[-1,]
  output_ls$Modularity$Value <- as.numeric(output_ls$Modularity$Value)
  
  ## Modularity ----
  gplot_df <- output_ls$Modularity
  Mod_gg <- ggplot(gplot_df, aes(x = factor(Model, 
                                            levels = names(data_ls[[1]])), 
                                 y = factor(Treatment, 
                                            levels = names(data_ls)), 
                                 fill = Value)
  ) +
    geom_tile(color = "black",
              lwd = 0.5,
              linetype = 1) + 
    coord_fixed() +
    guides(fill = guide_colourbar(barwidth = 2,
                                  barheight = 15,
                                  title = "Modularity")) + 
    theme_bw() + 
    theme(axis.text.x=element_text(angle = -20, hjust = 0)) + 
    scale_fill_gradient2(low = "#5ab4ac", high = "#d8b365") + 
    labs(y = "Treatment") + 
    scale_size_continuous(range = c(2, 10), breaks=pretty_breaks(8), name = "Effect Size") +
    guides(size = guide_legend(override.aes = list(shape = 25), order = 3))
  
  ## Strength ----
  # difference calculation
  output_ls$Centrality$Strength <- na.omit(output_ls$Centrality$Strength)
  NodeTopo_df <- expand.grid(Shared_abbr, names(data_ls[[1]]), names(data_ls))
  colnames(NodeTopo_df) <- colnames(output_ls$Centrality$Strength)[-2]
  NodeTopo_df$Value <- NA
  
  match_vec <- match(with(output_ls$Centrality$Strength, paste(Species, Model, Treatment)),
                     with(NodeTopo_df, paste(Species, Model, Treatment)))
  
  NodeTopo_df$Value[match_vec] <- output_ls$Centrality$Strength$Value
  
  Strength_gg <- ggplot(NodeTopo_df, 
                        aes(x = factor(Model, levels = names(data_ls[[1]])), 
                            y = Species, fill = Value)) + 
    geom_tile(color = "black",
              lwd = 0.5,
              linetype = 1) + 
    coord_fixed() + 
    geom_label(label= round(NodeTopo_df$Value, 2)) + 
    guides(fill = guide_colourbar(barwidth = 2,
                                  barheight = 15,
                                  title = "Node Strength")) + 
    ## axes
    theme_bw() + 
    theme(axis.text.x=element_text(angle = -20, hjust = 0)) + 
    facet_wrap(~factor(Treatment, levels = names(data_ls))) +
    scale_fill_gradient2(low = "#5ab4ac", high = "#d8b365")
  
  # ## Eigenvector ----
  # # difference calculation
  # data2 <- output_ls$Centrality$Eigenvector
  # if(length(levels(factor(data2$Treatment))) == 2){
  #   for(k in unique(data2$Model)){
  #     diff_ls <- split(data2[data2$Model == k,], data2[data2$Model == k,"Treatment"])
  #     diff_ls <- lapply(diff_ls, FUN = Sort.DF, Column = "Species")
  #     Diffs <- diff_ls[[TreatmentOrder[1]]]$Value - diff_ls[[TreatmentOrder[2]]]$Value
  # 
  # 
  #     data3 <- data.frame(Species = diff_ls[[TreatmentOrder[1]]]$Species,
  #                         Value = Diffs,
  #                         Model = rep(k, length(Diffs)),
  #                         Treatment = rep("Difference", length(Diffs)))
  #     data2 <- rbind(data2, data3)
  #   }
  # }
  # output_ls$Centrality$Eigenvector <- data2
  # Eigenvector_gg <- ggplot(output_ls$Centrality$Eigenvector, aes(x = factor(Model, levels = ModelOrder), y = Species, fill = Value)) +
  #   geom_tile(color = "black",
  #             lwd = 0.5,
  #             linetype = 1) +
  #   coord_fixed() +
  #   geom_label(label= round(output_ls$Centrality$Eigenvector$Value, 2)) +
  #   guides(fill = guide_colourbar(barwidth = 2,
  #                                 barheight = 15,
  #                                 title = "Eigenvector Centrality")) +
  #   ## axes
  #   theme_bw() +
  #   theme(axis.text.x=element_text(angle = -20, hjust = 0)) +
  #   facet_wrap(~factor(Treatment, levels = TreatmentOrder)) +
  #   scale_fill_gradient2(low = "#5ab4ac", high = "#d8b365")
  # 
  
  plot_ls <- list(Strength = Strength_gg,
                  # Eigenvector = Eigenvector_gg,
                  Modularity = Mod_gg
                  #, Nestedness = Nest_gg
  )
  return(list(plots = plot_ls,
              numbers = output_ls)
  )
}
