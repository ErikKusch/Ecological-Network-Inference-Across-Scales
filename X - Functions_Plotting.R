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
Plot.Network <- function(Plot_df = NULL, 
                         ModelOrder = c("COCCUR", "NETASSOC", "HMSC", "IF_REM"),
                         TreatmentOrder = c("Pre-Fire", "Post-Fire"),
                         FName = "Cross-YFDP", 
                         Dir){
  ## column identification
  Plot_cols <- grep(colnames(Plot_df), pattern = "effects")
  
  ## species indetification
  Species_ls <- base::strsplit(x = as.character(unique(Plot_df$Partner1)), split = " ")
  SpeciesLegend <- data.frame(Abbr = paste0(substr(unlist(lapply(Species_ls, "[[", 1)), start = 1, stop = 2),
                                            substr(unlist(lapply(Species_ls, "[[", 2)), start = 1, stop = 2)),
                              Full = unique(Plot_df$Partner1)
  )
  ## plot generation
  Plots_ls <- as.list(rep(NA, length(Plot_cols)+1))
  names(Plots_ls) <- c(colnames(Plot_df)[Plot_cols], "Legend")
  for(Plot_Iter in 1:(length(Plots_ls)-1)){
    Plot_Counter <- Plot_cols[Plot_Iter]
    directed <- ifelse(length(grep(colnames(Plot_df)[Plot_Counter], pattern = "IF_REM")) == 1, TRUE, FALSE)
    Iter_df <- Plot_df[,c(1:2, Plot_Counter, Plot_Counter+1)]
    Iter_df[which(is.na(Iter_df[,4])),4] <- FALSE
    Iter_df$Partner1 <- SpeciesLegend$Abbr[match(Iter_df$Partner1 , SpeciesLegend$Full)]
    Iter_df$Partner2 <- SpeciesLegend$Abbr[match(Iter_df$Partner2 , SpeciesLegend$Full)]
    Iter_df[,3][Iter_df[,3] == 0] <- NA
    name_plot <- names(Plots_ls)[Plot_Iter]
    name_plot <- strsplit(x = name_plot, split = "[.]")
    name_plot <- paste(unlist(lapply(name_plot, "[[", 1)), unlist(lapply(name_plot, "[[", 2)), sep = " - ")
    Plots_ls[[Plot_Iter]] <- qgraph::qgraph(Iter_df[,1:3], # the three-column data frame of "Partner1", "Partner2", and "Effect"
                                            bg = "white",                       
                                            color = "white", # colour of nodes
                                            layout = 'groups', # layout (I have found this one to be consistent across al plots)
                                            esize = 7, # line width of strongest connection/edge
                                            negCol = "#ba00bd", # colour of negative edges
                                            posCol = "#d9db3b", # colour of positive edges
                                            palette = "colorblind", # colour palette
                                            fade = TRUE, # fading effect on edges relative to strength
                                            directed = directed, # directed network or not
                                            DoNotPlot = TRUE,
                                            lty = Iter_df[,4]*-1+2, # significant interactions are solid, non-significant ones are dashed
                                            title = name_plot
    )
  }
  Plots_ls[[(length(Plot_cols)+1)]] <- SpeciesLegend
  ggsave(grid.table(Plots_ls[[length(Plots_ls)]]), filename = file.path(Dir, paste0(FName, "_Legend.png")))
  Identifiers <- strsplit(x = names(Plots_ls)[-length(Plots_ls)], split = "[.]")
  Id_df <- data.frame(Model = unlist(lapply(Identifiers, "[[", 1)),
                      Treatment = unlist(lapply(Identifiers, "[[", 2)))
  Models <- unique(Id_df$Model)
  Treatments <- unique(Id_df$Treatment)
  Plotting_df <- expand.grid(TreatmentOrder, ModelOrder)
  layout(matrix(1:(length(Plots_ls)-1), 
                nrow = length(Models), 
                ncol = length(Treatments), 
                byrow = TRUE))
  for(Plot_Iter in 1:nrow(Plotting_df)){
    Plot_pos <- intersect(which(Plotting_df$Var1[Plot_Iter] == Id_df$Treatment),
                          which(Plotting_df$Var2[Plot_Iter] == Id_df$Model))[1]
    plot(Plots_ls[[Plot_pos]])
  } 
}

# NETWORKS AS MATRICES ===================================================
# plotting of networks belonging to the same method as matrices
Plot.NetMat <- function(method = "COCCUR", data = YFDP_df, TreatmentOrder = c("Pre-Fire", "Post-Fire", "Difference"), datareturn = FALSE){ ## currently only works for YFDP
  plot_df <- data[, c(1:2, grep(colnames(data), pattern = method))]
  
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

# NETWORKS AS MATRICES ACROSS SCALES =====================================
# plotting of networks belonging to the same method as matrices across scales
Plot.NetMat.Cross <- function(data = NULL,
                              ALLRegionOrder = c("Pre-Fire", "Yosemite", "Temperate Broadleaf & Mixed Forests"),
                              ALLMethodOrder = c("COCCUR", "NETASSOC", "HMSC", "IF_REM")){
  
  plot_ls <- as.list(rep(NA, length(ALLMethodOrder)))
  names(plot_ls) <- ALLMethodOrder
  
  for(method in ALLMethodOrder){
    
    plot_df <- data[, c(1:2, grep(colnames(data), pattern = method))]
    
    if("Difference" %in% ALLRegionOrder){
      if(method == "HMSC"){
        
        if(sum(grepl(pattern = "diametre", x = colnames(data))) == 0){
          plot_df$Difference.diametre <- plot_df[, grep(x = colnames(plot_df), pattern = "biomass")[1]] - plot_df[, grep(x = colnames(plot_df), pattern = "biomass")[3]]
          plot_df$Difference.diametre.Sig <- ifelse(rowSums(plot_df[, grep(x = colnames(plot_df), pattern = "biomass")[c(2,4)]]) == 2, TRUE, FALSE)
        }else{
          plot_df$Difference.diametre <- plot_df[, grep(x = colnames(plot_df), pattern = "diametre")[1]] - plot_df[, grep(x = colnames(plot_df), pattern = "diametre")[3]]
          plot_df$Difference.diametre.Sig <- ifelse(rowSums(plot_df[, grep(x = colnames(plot_df), pattern = "diametre")[c(2,4)]]) == 2, TRUE, FALSE)
        }
        
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
      plot_df$Condition <- gsub(plot_df$Condition, pattern = ".biomass", replacement = "")
      plot_df$Condition <- gsub(plot_df$Condition, pattern = ".diametre", replacement = "")
    } 
    
    if(method == "NETASSOC" | method == "IF_REM"){
      gplot <- ggplot(plot_df, aes(x = Partner2, y = Partner1, fill = Value))
    }else{
      gplot <- ggplot(plot_df, aes(x = Partner1, y = Partner2, fill = Value))
    }
    gplot <- gplot +
      geom_tile(color = "white",
                lwd = 1.5,
                linetype = 1) + 
      ## colours, look, and legend
      geom_point(aes(shape = Sig), size = 4) +
      scale_shape_manual(values=c(32, 15), na.translate = FALSE, name = "Significance") +  
      guides(shape = FALSE) + 
      coord_fixed() + 
      guides(fill = guide_colourbar(barwidth = 2,
                                    barheight = 15,
                                    title = "Association")) + 
      ## axes
      theme_bw() + 
      theme(axis.text.x=element_text(angle = -20, hjust = 0))
    
    if(all(plot_df$Value[!is.na(plot_df$Value)] > 0)){
      gplot <- gplot + scale_fill_gradient(high = "#d9db3b")
    }else{
      if(all(plot_df$Value[!is.na(plot_df$Value)] < 0)){
        gplot <- gplot + scale_fill_gradient(high = "#ba00bd")
      }else{
        gplot <- gplot + scale_fill_gradient2(low = "#ba00bd", high = "#d9db3b")
      }
    }
    
    if(method != "IF_REM"){
      gplot <- gplot + labs(x = "", y = "")
    }else{
      gplot <- gplot + labs(y = "Actor", x = "Subject")
    }
    gplot <- gplot + facet_wrap(~factor(Condition, levels = ALLRegionOrder))
    plot_ls[[method]] <- gplot + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
  }
    return(plot_grid(plotlist=plot_ls, ncol = 1, labels = ALLMethodOrder))
  }

# NETWORK TOPOLOGIES =====================================================
## for calculation and visualisation of topology metrics
Calc.Topology <- function(data = YFDP_df, Sig = TRUE, Model = "HMSC", TreatmentOrder = c("Pre-Fire", "Post-Fire"), EffSize = FALSE){
  if(length(TreatmentOrder) == 2){TreatmentOrder <- c(TreatmentOrder, "Difference")}
  
  if(Model != "ALL"){
    data <- data[,c(1:2,grep(colnames(data), pattern = Model))]
  }else{
    HMSC_col <- grep(colnames(data), pattern = "HMSC")
    Targetcol <- c(grep(colnames(data), pattern = "diametre"), grep(colnames(data), pattern = "biomass"))
    HMSC_discard <- HMSC_col[HMSC_col %nin% Targetcol]
    if(length(HMSC_discard)>0){data <- data[,-HMSC_discard]}
  }
  
  # model identification
  Models_vec <- colnames(data)[grep(colnames(data), pattern = "effects")]
  
  # output objects
  Base_df <- data.frame(Species = NA,
                        Value = NA,
                        Model = NA,
                        Treatment = NA)
  output_ls <- list(Centrality = list(Strength = Base_df,
                                      Eigenvector = Base_df),
                    Modularity = Base_df,
                    Nestedness = Base_df)
  output_ls$Nestedness$Warning = FALSE
  
  # loop over models here
  for(i in Models_vec){
    print(paste("#####", i))
    ## model characteristics
    Char_vec <- unlist(strsplit(x = i, split = "[.]"))
    if(Model == "HMSC"){
      Char_vec <- c(Char_vec[3], Char_vec[2])
    }else{
      Char_vec <- c(Char_vec[1], Char_vec[2])
    }
    
    ## network building
    ## only keep singificant interactions when Sig == TRUE
    if(Sig){
      NonSigPos <- which(!data[,gsub(i, pattern = ".effects", replacement = ".Sig")])
      if(length(NonSigPos)>0){
        data_i <- data[-NonSigPos,]
      }else{
        data_i <- data
      }
    }else{
      data_i <- data
    }
    data_i <- na.omit(data_i[,c("Partner1", "Partner2", i)])
    colnames(data_i)[3] <- "weight"
    g <- igraph::graph_from_data_frame(d = data_i, directed = TRUE)
    g <- as.undirected(g, mode = "collapse")
    Vnames <- V(g)$name
    addv <- as.character(unique(data$Partner1)[unique(data$Partner1) %nin% Vnames])
    g <- add_vertices(g, length(addv))
    g <- set.vertex.attribute(g, "name", value = c(Vnames, addv))
    if(!is.null(E(g)$weight)){
      E(g)$weight <- abs(E(g)$weight)
    }
    
    ## Node-Level Topology ----
    ### Strength
    print("Node Strength")
    strength <- range01(sort(strength(g), decreasing = TRUE))
    strength_df <- data.frame(Species = names(strength),
                              Value = strength,
                              matrix(rep(Char_vec, each = length(strength)), ncol = 2))
    rownames(strength_df) <- c()
    colnames(strength_df) <- colnames(Base_df)
    output_ls$Centrality$Strength <- rbind(output_ls$Centrality$Strength, strength_df)
    
    ### Eigenvector Centrality
    print("Eigenvectors")
    ecvent <- sort(evcent(g, scale = TRUE)$vector, decreasing = TRUE)
    ecvent_df <- data.frame(Species = names(ecvent),
                            Value = ecvent,
                            matrix(rep(Char_vec, each = length(ecvent)), ncol = 2))
    rownames(ecvent_df) <- c()
    colnames(ecvent_df) <- colnames(Base_df)
    output_ls$Centrality$Eigenvector <- rbind(output_ls$Centrality$Eigenvector, ecvent_df)
    
    ## Network-Level Topology ----
    if(length(E(g)) > 0){
      ### randomisation of links for effect size calculation
      graph_ls <- as.list(rep(NA, 1e3))
      adjmat <- as.matrix(as_adjacency_matrix(g, attr = "weight"))
      set.seed(42)
      for(boot in 1:1e3){
        adjmati <- adjmat
        adjmati[1:length(adjmat)] <- adjmat[sample(1:length(adjmat), length(adjmat), replace = FALSE)]
        graph_ls[[boot]] <- graph_from_adjacency_matrix(adjmati, mode = "undirected", weighted = TRUE)
      }
      ### Modularity
      print("Modularity")
      Fun.Modu <- function(g){
        E(g)$weight[is.na(E(g)$weight)] <- 0
        e <- simplify(g)
        fc <- cluster_optimal(g)
        modularity <- modularity(e, membership = fc$membership, weights = E(e)$weight)
        return(modularity)
      }
      Modu_infer <- Fun.Modu(g)
      if(EffSize){
        Modu_rand <- unlist(lapply(graph_ls, FUN = Fun.Modu))
        Modu_effsize <- Modu_infer-mean(Modu_rand)/sd(Modu_rand)
      }else{
        Modu_effsize <- NA
      }
      output_ls$Modularity <- rbind(output_ls$Modularity, c(Modu_infer, Char_vec, Modu_effsize))
      
      ### Nestedness
      print("Nestedness")
      Fun.Nest <- function(g){
        maxnodf_mat <- as.matrix(as_adjacency_matrix(g, attr = "weight"))
        maxnodf_mat <- maxnodf_mat[rowSums(maxnodf_mat) > 0, colSums(maxnodf_mat) > 0]
        maxnodf_mat[maxnodf_mat > 0] <- 1
        nestedness <- tryCatch(maxnodf(maxnodf_mat), warning = function(w) { })
        if(is.null(nestedness)){
          Nest_vec <- c(NA, Char_vec, TRUE)
        }else{
          test <- tryCatch(maxnodf(maxnodf_mat), warning = function(w)w)
          if(!is.null(test$message)){
            if(startsWith(test$message, "Number of links does not satisfy Number of links")){
              Nest_vec <- c(nestedness$max_nodf, Char_vec, TRUE)
            }else{
              Nest_vec <- c(nestedness$max_nodf, Char_vec, FALSE)
            }
          }else{
            Nest_vec <- c(nestedness$max_nodf, Char_vec, FALSE)
          }
        }
        return(Nest_vec)
      }
      Nest_infer <- hush(Fun.Nest(g))
      if(EffSize){
        Nest_rand <- hush(lapply(graph_ls, FUN = Fun.Nest))
        Nest_rand <- data.frame(Value = as.numeric(unlist(lapply(Nest_rand, "[[", 1))),
                                Warning = as.logical(unlist(lapply(Nest_rand, "[[", 4))))
        Nest_rand <- Nest_rand$Value[!Nest_rand$Warning]
        if(!is.na(Nest_infer[1])){
          Nest_effsize <- as.numeric(Nest_infer[1])-mean(Nest_rand)/sd(Nest_rand)
        }else{
          Nest_effsize <- NA
        }
      }else{
        Nest_effsize <- NA
      }
      output_ls$Nestedness <- rbind(output_ls$Nestedness, c(Nest_infer, Nest_effsize))
    }else{
      output_ls$Modularity <- rbind(output_ls$Modularity, c(NA, Char_vec, NA))
      output_ls$Nestedness <- rbind(output_ls$Nestedness, c(NA, Char_vec, NA))
    }
  }
  ## cleaning output object
  colnames(output_ls$Nestedness) <- c(colnames(output_ls$Nestedness)[-1], "EffectSize")
  colnames(output_ls$Modularity) <- c(colnames(output_ls$Modularity)[-1], "EffectSize")
  output_ls$Centrality$Strength <- output_ls$Centrality$Strength[-1,]
  output_ls$Centrality$Eigenvector <- output_ls$Centrality$Eigenvector[-1,]
  output_ls$Modularity <- output_ls$Modularity[-1,]
  output_ls$Modularity$Value <- as.numeric(output_ls$Modularity$Value)
  output_ls$Nestedness <- output_ls$Nestedness[-1,]
  output_ls$Nestedness$Value <- as.numeric(output_ls$Nestedness$Value)
  
  ## Nestedness ----
  # difference calculation
  data2 <- output_ls$Nestedness
  if(length(levels(factor(data2$Treatment))) == 2){
    for(k in unique(data2$Model)){
      diff_ls <- split(data2[data2$Model == k,], data2[data2$Model == k,"Treatment"])
      Diffs <- diff_ls[[TreatmentOrder[1]]]$Value - diff_ls[[TreatmentOrder[2]]]$Value
      data3 <- data.frame(Value = Diffs,
                          Model = k,
                          Treatment = "Difference",
                          Warning = !any(c(diff_ls[[1]]$Warning, diff_ls[[2]]$Warning) == "FALSE"),
                          "EffectSize" = NA
      )
      data2 <- rbind(data2, data3)
    }
  }
  output_ls$Nestedness <- data2
  output_ls$Nestedness$Warning <- !as.logical(output_ls$Nestedness$Warning) # convert to significance
  
  gplot_df <- output_ls$Nestedness
  gplot_df$EffectSize <- as.numeric(gplot_df$EffectSize)
  gplot_df$Sign <- sign(gplot_df$EffectSize)
  
  Nest_gg <- ggplot(gplot_df, aes(x = Model, y = factor(Treatment, levels = TreatmentOrder), fill = Value, shape = as.factor(Sign), size = abs(EffectSize))) +
    geom_tile(color = "black",
              lwd = 0.5,
              linetype = 1) + 
    coord_fixed() +
    guides(fill = guide_colourbar(barwidth = 2,
                                  barheight = 15,
                                  title = "Nestedness")) + 
    theme_bw() + 
    theme(axis.text.x=element_text(angle = -20, hjust = 0)) + 
    scale_fill_gradient2(low = "#5ab4ac", high = "#d8b365") + 
    labs(y = "Treatment") + 
    geom_point(fill = "black") + 
    scale_size_continuous(range = c(2, 10), breaks=pretty_breaks(8), name = "Effect Size") + 
    scale_shape_manual(values=c(25, 24), na.translate = FALSE, name = "Effect Sign") + 
    guides(shape = guide_legend(override.aes = list(size = 4), order = 2),
           size = guide_legend(override.aes = list(shape = 25), order = 3))
  
  ## Modularity ----
  # difference calculation
  data2 <- output_ls$Modularity
  if(length(levels(factor(data2$Treatment))) == 2){
    for(k in unique(data2$Model)){
      diff_ls <- split(data2[data2$Model == k,], data2[data2$Model == k,"Treatment"])
      Diffs <- diff_ls[[TreatmentOrder[1]]]$Value - diff_ls[[TreatmentOrder[2]]]$Value
      data3 <- data.frame(Value = Diffs,
                          Model = k,
                          Treatment = "Difference",
                          EffectSize = NA
      )
      data2 <- rbind(data2, data3)
    }
  }
  output_ls$Modularity <- data2
  
  gplot_df <- output_ls$Modularity
  gplot_df$EffectSize <- as.numeric(gplot_df$EffectSize)
  gplot_df$Sign <- sign(gplot_df$EffectSize)
  
  Mod_gg <- ggplot(gplot_df, aes(x = Model, y = factor(Treatment, levels = TreatmentOrder), fill = Value, shape = as.factor(Sign), size = abs(EffectSize))) +
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
    geom_point(fill = "black") + 
    scale_size_continuous(range = c(2, 10), breaks=pretty_breaks(8), name = "Effect Size") + 
    scale_shape_manual(values=c(25, 24), na.translate = FALSE, name = "Effect Sign") + 
    guides(shape = guide_legend(override.aes = list(size = 4), order = 2),
           size = guide_legend(override.aes = list(shape = 25), order = 3))

  ## Strength ----
  # difference calculation
  data2 <- output_ls$Centrality$Strength
  if(length(levels(factor(data2$Treatment))) == 2){
    # TreatmentOrder <- c(TreatmentOrder, "Difference")
    for(k in unique(data2$Model)){
      diff_ls <- split(data2[data2$Model == k,], data2[data2$Model == k,"Treatment"])
      diff_ls <- lapply(diff_ls, FUN = Sort.DF, Column = "Species")
      Diffs <- diff_ls[[TreatmentOrder[1]]]$Value - diff_ls[[TreatmentOrder[2]]]$Value
      
      
      data3 <- data.frame(Species = diff_ls[[TreatmentOrder[1]]]$Species,
                          Value = Diffs,
                          Model = rep(k, length(Diffs)),
                          Treatment = rep("Difference", length(Diffs)))
      data2 <- rbind(data2, data3)
    }
  }
  output_ls$Centrality$Strength <- data2
  
  Strength_gg <- ggplot(output_ls$Centrality$Strength, aes(x = Model, y = Species, fill = Value)) + 
    geom_tile(color = "black",
              lwd = 0.5,
              linetype = 1) + 
    coord_fixed() + 
    geom_label(label= round(output_ls$Centrality$Strength$Value, 2)) + 
    guides(fill = guide_colourbar(barwidth = 2,
                                  barheight = 15,
                                  title = "Node Strength")) + 
    ## axes
    theme_bw() + 
    theme(axis.text.x=element_text(angle = -20, hjust = 0)) + 
    facet_wrap(~factor(Treatment, levels = TreatmentOrder)) +
    scale_fill_gradient2(low = "#5ab4ac", high = "#d8b365")
  
  ## Eigenvector ----
  # difference calculation
  data2 <- output_ls$Centrality$Eigenvector
  if(length(levels(factor(data2$Treatment))) == 2){
    for(k in unique(data2$Model)){
      diff_ls <- split(data2[data2$Model == k,], data2[data2$Model == k,"Treatment"])
      diff_ls <- lapply(diff_ls, FUN = Sort.DF, Column = "Species")
      Diffs <- diff_ls[[TreatmentOrder[1]]]$Value - diff_ls[[TreatmentOrder[2]]]$Value
      
      
      data3 <- data.frame(Species = diff_ls[[TreatmentOrder[1]]]$Species,
                          Value = Diffs,
                          Model = rep(k, length(Diffs)),
                          Treatment = rep("Difference", length(Diffs)))
      data2 <- rbind(data2, data3)
    }
  }
  output_ls$Centrality$Eigenvector <- data2
  Eigenvector_gg <- ggplot(output_ls$Centrality$Eigenvector, aes(x = Model, y = Species, fill = Value)) + 
    geom_tile(color = "black",
              lwd = 0.5,
              linetype = 1) + 
    coord_fixed() + 
    geom_label(label= round(output_ls$Centrality$Eigenvector$Value, 2)) + 
    guides(fill = guide_colourbar(barwidth = 2,
                                  barheight = 15,
                                  title = "Eigenvector Centrality")) + 
    ## axes
    theme_bw() + 
    theme(axis.text.x=element_text(angle = -20, hjust = 0)) + 
    facet_wrap(~factor(Treatment, levels = TreatmentOrder)) +
    scale_fill_gradient2(low = "#5ab4ac", high = "#d8b365")
  
  
  plot_ls <- list(Strength = Strength_gg,
                  Eigenvector = Eigenvector_gg,
                  Modularity = Mod_gg,
                  Nestedness = Nest_gg)
  return(list(plots = plot_ls,
              numbers = output_ls)
  )
}