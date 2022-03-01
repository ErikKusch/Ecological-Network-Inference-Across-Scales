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
Plot.Network <- function(Plot_df = NULL){
  Plots_ls <- as.list(rep(NA, ncol(Plot_df)-2))
  names(Plots_ls) <- colnames(Plot_df)[-1:-2]
  for(Plot_Iter in 1:length(Plots_ls)){
    Plot_Counter <- Plot_Iter+2
    directed <- ifelse(length(grep(colnames(Plot_df)[Plot_Counter], pattern = "IF_REM")) == 1, TRUE, FALSE)
    Iter_df <- Plot_df[,c(1:2, Plot_Counter)]
    Iter_df[Iter_df == 0] <- NA
    Plots_ls[[Plot_Iter]] <- qgraph::qgraph(Iter_df, # the three-column data frame of "Partner1", "Partner2", and "Effect"
                                  color = "grey", # colour of nodes
                                  layout = 'groups', # layout (I have found this one to be consistent across al plots)
                                  esize = 7, # line width of strongest connection/edge
                                  # negCol = 'darkred', # colour of negative edges
                                  # posCol = 'green', # colour of positive edges
                                  palette = "colorblind", # colour palette
                                  fade = TRUE, # fading effect on edges relative to strength
                                  directed = directed, # directed network or not
                                  DoNotPlot = TRUE
    )
  }
  return(Plots_ls)
}

# PLOT COMPARISON ========================================================
# plotting of networks belonging to the same method/treatment/data set next to each other
Plots.Compare <- function(Compare_df = NULL,
                          Header = "COCCUR", # can be anything occurring in the Method, Treatment, or Model column that is created within this function from Compare_df
                          SubsetHeader = NULL, 
                          X = "Data Set", # X axis of resulting panel plot
                          Y = "Treatment", # y axis of resulting panel plot
                          SubsetX = NULL, # further subsetting of X
                          SubsetY = c("Pre-Fire", "ALL"), # further subsetting of Y
                          FName = "Test",
                          Dir = Dir.Exports, 
                          HMSCFilter = "diameter") 
  {
  # Data Reformatting ####
  Compare_long <- melt(Compare_df, id.vars = c("Partner1", "Partner2"), variable.name = "Identifier")
  colnames(Compare_long)[4] <- "Inference"
  colnames(Compare_long)[5] <- "Data Set"
  Compare_long$Partner1 <- as.character(Compare_long$Partner1)
  Compare_long$Partner2 <- as.character(Compare_long$Partner2)
  Compare_long$`Data Set` <- as.character(Compare_long$`Data Set`)
  Identifier_ls <- strsplit(x = as.character(Compare_long$Identifier), split = "[.]")
  Compare_long$Method <- unlist(lapply(Identifier_ls, `[`, 1))
  Compare_long$Treatment <- unlist(lapply(Identifier_ls, `[`, 2))
  Compare_long$Model <- unlist(lapply(Identifier_ls, `[`, 3))
  Compare_long <- Compare_long[,-3]
  Compare_long$Sp1 <- Compare_long$Partner1
  Compare_long$Partner1 <- unlist(lapply(lapply(strsplit(Compare_long$Partner1, " "), FUN = substr, start = 1, stop = 2), paste, collapse=""))
  Compare_long$Sp2 <- Compare_long$Partner2
  Compare_long$Partner2 <- unlist(lapply(lapply(strsplit(Compare_long$Partner2, " "), FUN = substr, start = 1, stop = 2), paste, collapse=""))
  
  # Subsetting ####
  ## Subsetting according to Header
  Iter_pos <- which(Compare_long$Method == Header | 
                      Compare_long$Treatment == Header | 
                      Compare_long$Model == Header |
                      Compare_long$`Data Set` == Header)
  Iter_df <- Compare_long[Iter_pos, ]
  
  ## Further subsetting 
  if(!is.null(SubsetHeader)){
    Iter_pos <- which(Compare_long$Method == SubsetHeader | 
                        Compare_long$Treatment == SubsetHeader | 
                        Compare_long$Model == SubsetHeader |
                        Compare_long$`Data Set` == SubsetHeader)
    Iter_df <- Compare_long[Iter_pos, ]
  }
  if(!is.null(SubsetX)){
    Iter_df <- Iter_df[which(Iter_df[ , X] %in% SubsetX), ]
  }
  if(!is.null(SubsetY)){
    Iter_df <- Iter_df[which(Iter_df[ , Y] %in% SubsetY), ]
  }
  if(HMSCFilter != "ALL"){
    HMSC_pos <- Iter_df$Model == HMSCFilter
    HMSC_pos[is.na(HMSC_pos)] <- TRUE
    Iter_df <- Iter_df[HMSC_pos, ]
  }
  
  # SPECIES LEGEND -----
  SpeciesLegend <- data.frame(Code = unique(Iter_df$Partner1),
                              Species = unique(Iter_df$Sp1))
  
  # NETWORKS -----
  ## Plot Layout and Combinations of X and Y ####
  X_pos <- unique(Iter_df[, X])
  Y_pos <- unique(Iter_df[, Y])
  Combines_df <- expand.grid(X_pos, Y_pos)
  colnames(Combines_df) <- c("X", "Y")
  if(all(is.na(Iter_df$Inference))){
    stop("No inferred interactions present for current selection")
  }else{
    Max_df <- aggregate(abs(Inference) ~ Method, data = Iter_df, max)
  }
  if(!is.null(SubsetHeader)){
    if(SubsetHeader == "HMSC" | Header == "HMSC"){
      Max_df <- aggregate(abs(Inference) ~ Method+Model, data = Iter_df, max)
    }
  }
  
  ## Plot Creation ####
  for(Max_Iter in 1:2){
    Name <- ifelse(Max_Iter == 2, "Shared", "Indvidual")
    Plots_ls <- as.list(rep(NA, nrow(Combines_df)))
    names(Plots_ls) <- paste0(paste(Combines_df$Y, Combines_df$X, sep =" ("), ")")
    for(Counter in 1:nrow(Combines_df)){
      Plot_df <- Iter_df[Iter_df[, X] == Combines_df$X[Counter] & Iter_df[, Y] == Combines_df$Y[Counter], 1:3]
      Curr_Method <- Combines_df[Counter, c(sum(Max_df$Method %in% Combines_df$X), sum(Max_df$Method %in% Combines_df$Y)) > 0]
      if(class(Curr_Method) == "data.frame"){
        Curr_Method <- Combines_df[Counter, c(sum(Max_df$Model %in% Combines_df$X), sum(Max_df$Model %in% Combines_df$Y)) > 0]
      }
      if(class(Curr_Method) == "data.frame"){Curr_Method <- Header}
      Maximum <- ifelse(Max_Iter == 2, Max_df[which(Max_df$Method == Curr_Method), 2], 
                        max(abs(Plot_df$Inference), na.rm = TRUE))
      if(is.na(Maximum)){
        if(all(is.na(Plot_df$Inference))){
          Maximum = NULL
        }else{
          Maximum <- ifelse(Max_Iter == 2, Max_df[which(Max_df$Model == Curr_Method), 3], 
                            max(abs(Plot_df$Inference), na.rm = TRUE))
        }
      }
      if(nrow(Plot_df) == 0){next()}
      directed <- ifelse(Combines_df$Y[Counter] == "IF_REM" | 
                           Combines_df$X[Counter] == "IF_REM" | 
                           Header == "IF_REM", 
                         TRUE, FALSE)
      Plots_ls[[Counter]] <- qgraph::qgraph(Plot_df, # the three-column data frame of "Partner1", "Partner2", and "Effect"
                                            bg = "darkgrey", # background colour
                                            color = "white", # colour of nodes
                                            layout = 'groups', # layout (I have found this one to be consistent across al plots)
                                            # esize = 8, # line width of strongest connection/edge
                                            palette = "colorblind", # colour palette
                                            fade = FALSE, # fading effect on edges relative to strength
                                            directed = directed, # directed network or not
                                            DoNotPlot = TRUE,
                                            title = names(Plots_ls)[Counter],
                                            title.cex = 1, 
                                            maximum = Maximum
      )
    }
    
    ## Plot Combination into single output ####
    Combines2_df <- Combines_df[which(!is.na(Plots_ls)), ]
    Plots_ls <- Plots_ls[which(!is.na(Plots_ls))]
    png(file.path(Dir, paste0(FName, "_", Name, ".png")),
        height = 18*length(unique(Combines2_df$Y)),
        width = 12*length(unique(Combines2_df$X)),
        unit = "cm", res = 100)
    layout(matrix(1:length(Plots_ls), 
                  nrow = length(unique(Combines2_df$Y)), 
                  ncol = length(unique(Combines2_df$X)), 
                  byrow = TRUE))
    for(Plot_Iter in 1:nrow(Combines2_df)){
      plot(Plots_ls[[Plot_Iter]])
    } 
    dev.off()
  }
  ggsave(gridExtra::grid.table(SpeciesLegend), filename = file.path(Dir, paste0(FName, "_XLegend.png")))
}




# NETWORKS AS MATRICES ===================================================
# plotting of networks belonging to the same method as matrices
Plot.NetMat <- function(method = "COCCUR", data = YFDP_df){ ## currently only works for YFDP
  plot_df <- data[, c(1:2, grep(colnames(data), pattern = method))]
  if(method == "HMSC"){
    plot_df$Difference.abundance <- plot_df[, grep(x = colnames(plot_df), pattern = "abundance")[1]] - plot_df[, grep(x = colnames(plot_df), pattern = "abundance")[2]]
    plot_df$Difference.diametre <- plot_df[, grep(x = colnames(plot_df), pattern = "diametre")[1]] - plot_df[, grep(x = colnames(plot_df), pattern = "diametre")[2]]
    plot_df$Difference.presence_absence <- plot_df[, grep(x = colnames(plot_df), pattern = "presence_absence")[1]] - plot_df[, grep(x = colnames(plot_df), pattern = "presence_absence")[2]]
  }else{
    plot_df$Difference <- plot_df[,3] - plot_df[,4]
  }
  
  
  plot_df <- reshape(plot_df, 
                     direction = "long",
                     varying = list(names(plot_df)[3:ncol(plot_df)]),
                     v.names = "Value",
                     idvar = c("Partner1", "Partner2"),
                     timevar = "Condition",
                     times = names(plot_df)[3:ncol(plot_df)])
  plot_df$Condition <- gsub(plot_df$Condition, pattern = paste0(method, "."), replacement = "")
  
  if(method == "HMSC"){
    plot_df$Model <- unlist(lapply(strsplit(plot_df$Condition, ".", fixed = TRUE), "[[", 2))
    plot_df$Condition <- unlist(lapply(strsplit(plot_df$Condition, ".", fixed = TRUE), "[[", 1))
  }
  
  
  gplot <- ggplot(plot_df, aes(x = Partner1, y = Partner2, fill = Value)) +
    geom_tile(color = "white",
              lwd = 1.5,
              linetype = 1) + 
    ## colours, look, and legend
    
    coord_fixed() + 
    guides(fill = guide_colourbar(barwidth = 2,
                                  barheight = 20,
                                  title = "Association")) + 
    ## axes
    theme_bw() + 
    theme(axis.text.x=element_text(angle = -20, hjust = 0))
  
  if(all(plot_df$Value[!is.na(plot_df$Value)] > 0)){
    gplot <- gplot + scale_fill_gradient(high = "forestgreen")
  }else{
    if(all(plot_df$Value[!is.na(plot_df$Value)] < 0)){
      gplot <- gplot + scale_fill_gradient(high = "darkred")
  }else{
      gplot <- gplot + scale_fill_gradient2(low = "darkred", high = "forestgreen")
    }
    }
  
  if(method != "IF_REM"){
    gplot <- gplot + labs(x = "", y = "")
  }else{
    gplot <- gplot + labs(x = "Actor", y = "Subject")
  }
  
  if(method != "HMSC"){
    gplot <- gplot + facet_wrap(~factor(Condition, levels = c("Pre-Fire", "Post-Fire", "Difference")))
  }else{
    gplot <- gplot + facet_wrap(~ factor(Model, levels = c("presence_absence", "abundance", "diametre")) +
                                  factor(Condition, levels = c("Pre-Fire", "Post-Fire", "Difference"))
                                )
  }
  
  gplot
  return(gplot)
}










# FUN.PlotNetUncert <- function(Model = inter_mat, Dir = getwd(), Name = "Plot"){
#   ####### DATA LOADING ---------------------------------------------------------
#   mean_df <- apply(Model, c(1, 2), mean) # will return the mean estimate for every interaction (NB: this is the mean of the 80% posterior interval, so will be slightly different to the mean value returned from summary(fit), which is calculated from the full posterior distribution)  
#   sd_df <- apply(Model, c(1, 2), sd)
#   
#   png(filename = file.path(Dir, paste0(Name,".png")), units = "cm", width = 32, height = 32, res = 1000)
#   GGDAG_graph <- qgraph::qgraph(mean_df, 
#                                 edge.width = 1-sd_df, 
#                                 layout = 'circle',
#                                 negCol = 'darkred',
#                                 posCol = 'green',
#                                 labels = dimnames(mean_df)$species,
#                                 fade = TRUE,
#                                 directed = TRUE,
#                                 label.cex = 1.2
#   )
#   dev.off()
#   
#   
#   # colnames(mean_df) <- gsub(pattern = " ", replacement = "_", x = colnames(mean_df))
#   # rownames(mean_df) <- gsub(pattern = " ", replacement = "_", x = rownames(mean_df))
#   # 
#   # mean_df <- mean_df[, order(colnames(mean_df))] # sort columns
#   # mean_df <- mean_df*-1 # need to switch sign of results
#   # diag(mean_df) <- NA
#   # Interaction_hpdi <- apply(inter_mat, c(1, 2), HPDI, prob = 0.89) *-1 # need to switch sign of results
#   # min_df <- Interaction_hpdi[1,,]
#   # diag(min_df) <- NA
#   # max_df <- Interaction_hpdi[2,,]
#   # diag(max_df) <- NA
#   # 
#   # ####### DATA MANIPULATION ---------------------------------------------------------
#   # uncertainty_df <- max_df - min_df
#   # rownames(uncertainty_df) <- rownames(mean_df)
#   # colnames(uncertainty_df) <- colnames(mean_df)
#   # Interactions_igraph <- data.frame(Actor = rep(colnames(mean_df), each = nrow(mean_df)),
#   #                                   Subject = rep(rownames(mean_df), ncol(mean_df)),
#   #                                   Estimate = as.vector(unlist(mean_df)),
#   #                                   Uncertainty = as.vector(unlist(uncertainty_df))
#   # )
#   # Interactions_DAG <- Interactions_igraph
#   # Interactions_DAG <- na.omit(Interactions_DAG)
#   # 
#   # ####### PLOTTING ---------------------------------------------------------
#   # opacity <- 1-Interactions_DAG$Uncertainty/abs(Interactions_DAG$Estimate) # alpha = 1-percentage of uncertainty to corresponding estimate; this is questionable because uncertainties now get treated differently depending on estimate of the effect they are associated with
#   # opacity[opacity<.1] <- 0.05 # make all edges for which uncertainty is bigger than the effect (opacity < 0) or very large slightly visible
#   # opacity[opacity>1] <- 1
#   # 
#   # width <- abs(Interactions_DAG$Estimate)
#   # width <- rescale(width, to = c(0.01, 1))
#   # ###. GGDAG ----
#   # Dag_Paths <- paste(paste(Interactions_DAG$Actor, "->", Interactions_DAG$Subject), collapse = " ")
#   # dag <- dagitty(x = paste0("dag {", Dag_Paths, "}"))
#   # tidy_dag <- tidy_dagitty(.dagitty = dag, layout = "circle")
#   # tidy_dag$data$weight <- width
#   # tidy_dag$data$label <- round(Interactions_DAG$Estimate, 2)
#   # tidy_dag$data$colour <- ifelse(Interactions_DAG$Estimate > 0, "#009933", "#cc0000")
#   # tidy_dag$data$alpha <- opacity
#   # # set.seed(42)
#   # # tidy_dag$data <- tidy_dag$data[sample(1:3e3, size = 1e2),] # for testing only
#   # GGDAG_graph <- ggplot(tidy_dag, aes(x = x, y = y, xend = xend, yend = yend)) +
#   #   geom_dag_point(colour = "purple", size = 4) +
#   #   geom_dag_text(colour = "black", size = 3) +
#   #   theme_dag() + 
#   #   geom_dag_edges_arc(aes(edge_colour = colour, 
#   #                          # edge_width = weight,
#   #                          # label = label,
#   #                          edge_alpha = weight #alpha
#   #   ), 
#   #   angle_calc = 'along', label_dodge = grid::unit(8, "points"), show.legend = FALSE) + 
#   #   ggtitle("Network + Uncertainty")
#   # ggsave(filename = file.path(Dir, paste0(Name,".png")), plot = GGDAG_graph, units = "cm", width = 24, height = 24)
#   # 
#   # ###. VISNETWORTK ---- 
#   # node_list <- tibble(id = unique(colnames(mean_df)))
#   # node_list <- node_list %>%
#   #   mutate(label = id)
#   # 
#   # edge_list <- tibble(from = Interactions_DAG$Actor, to = Interactions_DAG$Subject)
#   # edge_list <- edge_list %>%
#   #   add_column(width = width) %>%
#   #   add_column(color.color = ifelse(Interactions_DAG$Estimate > 0, "#009933", "#cc0000")) %>%
#   #   add_column(color.highlight = ifelse(Interactions_DAG$Estimate > 0, "#009933", "#cc0000")) %>%
#   #   add_column(color.opacity = opacity)
#   # 
#   # VisNet_graph <- visNetwork(nodes = node_list, edges = edge_list, height = 1920) %>% 
#   #   visIgraphLayout(layout = "layout_with_fr", randomSeed = 42) %>%
#   #   visNodes(shape = "circle") %>%
#   #   visEdges(arrows = "to", smooth = list(roundness = 0.3), color = list(opacity = opacity)) 
#   # visSave(VisNet_graph, file = file.path(Dir, paste0(Name, ".html")), selfcontained = TRUE, background = "darkgrey")  
#   
#   ###. HEATMAP ---- 
#   fv.colors = colorRampPalette(c("red","white","green")) ## define the color ramp
#   colorlut = fv.colors(100)[c(1,seq(50,100,length.out=99))] ## select colors to use
#   
#   diag(mean_df) <- 0
#   jpeg(filename = file.path(Dir.PlotNets.PFTC, paste0(Name, ".jpeg")), units = "cm", width = 32, height = 32, res = 1000)
#   pheatmap(mean_df, display_numbers = TRUE, number_format =  "%.3f",
#            color = c("red", "white","green"), 
#            breaks = c(min(mean_df, na.rm = TRUE), -0.00001, 0.00001, max(mean_df, na.rm = TRUE)), main = "Actors (Columns) x Subject (Rows)",
#            fontsize = 10, cluster_rows = FALSE, cluster_cols=FALSE)
#   dev.off()
# }
