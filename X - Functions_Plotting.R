#' ####################################################################### #
#' PROJECT: [PhD; X - NETWORK PLOTTING] 
#' CONTENTS: 
#'  - Functionality for plotting of networks
#'  DEPENDENCIES:
#'  - 
#' AUTHOR: [Erik Kusch]
#' ####################################################################### #

FUN.PlotNetUncert <- function(Model = inter_mat, Dir = getwd(), Name = "Plot"){
  ####### DATA LOADING ---------------------------------------------------------
  mean_df <- apply(Model, c(1, 2), mean) # will return the mean estimate for every interaction (NB: this is the mean of the 80% posterior interval, so will be slightly different to the mean value returned from summary(fit), which is calculated from the full posterior distribution)  
  sd_df <- apply(Model, c(1, 2), sd)
  
  png(filename = file.path(Dir, paste0(Name,".png")), units = "cm", width = 32, height = 32, res = 1000)
  GGDAG_graph <- qgraph::qgraph(mean_df, 
                                edge.width = 1-sd_df, 
                                layout = 'circle',
                                negCol = 'darkred',
                                posCol = 'green',
                                labels = dimnames(mean_df)$species,
                                fade = TRUE,
                                directed = TRUE,
                                label.cex = 1.2
  )
  dev.off()
  
  
  # colnames(mean_df) <- gsub(pattern = " ", replacement = "_", x = colnames(mean_df))
  # rownames(mean_df) <- gsub(pattern = " ", replacement = "_", x = rownames(mean_df))
  # 
  # mean_df <- mean_df[, order(colnames(mean_df))] # sort columns
  # mean_df <- mean_df*-1 # need to switch sign of results
  # diag(mean_df) <- NA
  # Interaction_hpdi <- apply(inter_mat, c(1, 2), HPDI, prob = 0.89) *-1 # need to switch sign of results
  # min_df <- Interaction_hpdi[1,,]
  # diag(min_df) <- NA
  # max_df <- Interaction_hpdi[2,,]
  # diag(max_df) <- NA
  # 
  # ####### DATA MANIPULATION ---------------------------------------------------------
  # uncertainty_df <- max_df - min_df
  # rownames(uncertainty_df) <- rownames(mean_df)
  # colnames(uncertainty_df) <- colnames(mean_df)
  # Interactions_igraph <- data.frame(Actor = rep(colnames(mean_df), each = nrow(mean_df)),
  #                                   Subject = rep(rownames(mean_df), ncol(mean_df)),
  #                                   Estimate = as.vector(unlist(mean_df)),
  #                                   Uncertainty = as.vector(unlist(uncertainty_df))
  # )
  # Interactions_DAG <- Interactions_igraph
  # Interactions_DAG <- na.omit(Interactions_DAG)
  # 
  # ####### PLOTTING ---------------------------------------------------------
  # opacity <- 1-Interactions_DAG$Uncertainty/abs(Interactions_DAG$Estimate) # alpha = 1-percentage of uncertainty to corresponding estimate; this is questionable because uncertainties now get treated differently depending on estimate of the effect they are associated with
  # opacity[opacity<.1] <- 0.05 # make all edges for which uncertainty is bigger than the effect (opacity < 0) or very large slightly visible
  # opacity[opacity>1] <- 1
  # 
  # width <- abs(Interactions_DAG$Estimate)
  # width <- rescale(width, to = c(0.01, 1))
  # ###. GGDAG ----
  # Dag_Paths <- paste(paste(Interactions_DAG$Actor, "->", Interactions_DAG$Subject), collapse = " ")
  # dag <- dagitty(x = paste0("dag {", Dag_Paths, "}"))
  # tidy_dag <- tidy_dagitty(.dagitty = dag, layout = "circle")
  # tidy_dag$data$weight <- width
  # tidy_dag$data$label <- round(Interactions_DAG$Estimate, 2)
  # tidy_dag$data$colour <- ifelse(Interactions_DAG$Estimate > 0, "#009933", "#cc0000")
  # tidy_dag$data$alpha <- opacity
  # # set.seed(42)
  # # tidy_dag$data <- tidy_dag$data[sample(1:3e3, size = 1e2),] # for testing only
  # GGDAG_graph <- ggplot(tidy_dag, aes(x = x, y = y, xend = xend, yend = yend)) +
  #   geom_dag_point(colour = "purple", size = 4) +
  #   geom_dag_text(colour = "black", size = 3) +
  #   theme_dag() + 
  #   geom_dag_edges_arc(aes(edge_colour = colour, 
  #                          # edge_width = weight,
  #                          # label = label,
  #                          edge_alpha = weight #alpha
  #   ), 
  #   angle_calc = 'along', label_dodge = grid::unit(8, "points"), show.legend = FALSE) + 
  #   ggtitle("Network + Uncertainty")
  # ggsave(filename = file.path(Dir, paste0(Name,".png")), plot = GGDAG_graph, units = "cm", width = 24, height = 24)
  # 
  # ###. VISNETWORTK ---- 
  # node_list <- tibble(id = unique(colnames(mean_df)))
  # node_list <- node_list %>%
  #   mutate(label = id)
  # 
  # edge_list <- tibble(from = Interactions_DAG$Actor, to = Interactions_DAG$Subject)
  # edge_list <- edge_list %>%
  #   add_column(width = width) %>%
  #   add_column(color.color = ifelse(Interactions_DAG$Estimate > 0, "#009933", "#cc0000")) %>%
  #   add_column(color.highlight = ifelse(Interactions_DAG$Estimate > 0, "#009933", "#cc0000")) %>%
  #   add_column(color.opacity = opacity)
  # 
  # VisNet_graph <- visNetwork(nodes = node_list, edges = edge_list, height = 1920) %>% 
  #   visIgraphLayout(layout = "layout_with_fr", randomSeed = 42) %>%
  #   visNodes(shape = "circle") %>%
  #   visEdges(arrows = "to", smooth = list(roundness = 0.3), color = list(opacity = opacity)) 
  # visSave(VisNet_graph, file = file.path(Dir, paste0(Name, ".html")), selfcontained = TRUE, background = "darkgrey")  
  
  ###. HEATMAP ---- 
  fv.colors = colorRampPalette(c("red","white","green")) ## define the color ramp
  colorlut = fv.colors(100)[c(1,seq(50,100,length.out=99))] ## select colors to use
  
  diag(mean_df) <- 0
  jpeg(filename = file.path(Dir.PlotNets.PFTC, paste0(Name, ".jpeg")), units = "cm", width = 32, height = 32, res = 1000)
  pheatmap(mean_df, display_numbers = TRUE, number_format =  "%.3f",
           color = c("red", "white","green"), 
           breaks = c(min(mean_df, na.rm = TRUE), -0.00001, 0.00001, max(mean_df, na.rm = TRUE)), main = "Actors (Columns) x Subject (Rows)",
           fontsize = 10, cluster_rows = FALSE, cluster_cols=FALSE)
  dev.off()
}
