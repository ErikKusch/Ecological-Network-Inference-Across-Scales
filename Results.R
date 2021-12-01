#' ####################################################################### #
#' PROJECT: [Data Simplification for Network Inference across Scales] 
#' CONTENTS: 
#'  - Network comparison and Plotting
#'  DEPENDENCIES:
#'  - 0 - Preamble.R
#'  - X - Functions_Data.R
#'  - X - Functions_Plotting.R
#'  - Scripts that need to be run prior:
#'    + YFDP.R 
#'    + FIA.R
#' AUTHOR: [Erik Kusch]
#' ####################################################################### #

# PREAMBLE =================================================================
rm(list=ls())
set.seed(42)

## Sourcing ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
source("0 - Preamble.R")
source("X - Functions_Data.R")
source("X - Functions_Plotting.R")

# DATA =====================================================================
message("############ LOADING DATA")

## LOADING +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
YFDP_ls <- Load.Results(Dir = DirEx.YFDP)
Region_ls <- Load.Results(Dir = DirEx.Region)

## MANIPULATION ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#### SPECIES IDENTITIES ####
YFDP_spec <- ReadOut.Species(YFDP_ls) # species which have been analysed
Region_spec <- ReadOut.Species(Region_ls) # species which have been analysed
Shared_spec <- Reduce(intersect, list(YFDP_spec, Region_spec)) # shared species between all analyses

#### FLATTENING LISTS ####
YFDP_ls <- Flatten.Lists(YFDP_ls)
Region_ls <- Flatten.Lists(Region_ls)
Region_ls <- BiomeNames.List(Region_ls)

# RESULTS ==================================================================

## SHARED SPECIES ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# this is where we will compare interactions/associations 1-1 between shared species of all analyses
#### LIMITING RESULTS TO SHARED SPECIES BETWEEN DATA SETS #####
YFDP_limited <- Limit.Lists(YFDP_ls, Shared_spec)
Region_limited <- Limit.Lists(Region_ls[-grep(names(Region_ls), pattern = "HMSC")], # currently omitting all HMSC results because they aren't ready yet
                              Shared_spec)

#### INTERACTION EXTRACTION ####
YFDP_df <- Eff.Data.Frame(List_ls = YFDP_limited)
Region_df <- Eff.Data.Frame(List_ls = Region_limited)

#### PLOTTING NETWORKS ####
YFDP_plots <- Plot.Network(Plot_df = YFDP_df)
Region_plots <- Plot.Network(Plot_df = Region_df)

#### COMPARING NETWORKS ####

#### HMSC model specification within scales --------------
Header_vec <- c("YFDP")
Dir.WithinScale <- file.path(Dir.Exports, "Within-Scale")
if(!dir.exists(Dir.WithinScale)){dir.create(Dir.WithinScale)}
for(Plots_Iter in 1:length(Header_vec)){
  Plots.Compare(Compare_df = list(YFDP = YFDP_df, 
                                  Region = Region_df), 
                Header = Header_vec[Plots_Iter], SubsetHeader = "HMSC",
                X = "Model", SubsetX = NULL, 
                Y = "Treatment", SubsetY = NULL, 
                HMSCFilter = "ALL",
                FName = paste0("HMSC-Specifications_", Header_vec[Plots_Iter]),
                Dir = Dir.WithinScale)
}   
# --> most information in diametre-driven HMSC models. Using those going forward

#### Methods within scales --------------
Header_vec <- c("YFDP", "Region")
Dir.WithinScale <- file.path(Dir.Exports, "Within-Scale")
if(!dir.exists(Dir.WithinScale)){dir.create(Dir.WithinScale)}
for(Plots_Iter in 1:length(Header_vec)){
  SubsetY <- ifelse(Plots_Iter == 1, "Pre-Fire", "Temperate Broadleaf & Mixed Forests")
  Plots.Compare(Compare_df = list(YFDP = YFDP_df, 
                                  Region = Region_df), 
                Header = Header_vec[Plots_Iter],
                SubsetHeader = NULL,
                X = "Method", SubsetX = NULL, 
                Y = "Treatment", SubsetY = SubsetY, 
                HMSCFilter = "diametre",
                FName = paste0("Within-Scale_", Header_vec[Plots_Iter]),
                Dir = Dir.WithinScale)
}

#### Methods across Treatments --------------
Header_vec <- c("YFDP", "Region")
Dir.WithinScale <- file.path(Dir.Exports, "Within-Scale")
if(!dir.exists(Dir.WithinScale)){dir.create(Dir.WithinScale)}
for(Plots_Iter in 1:length(Header_vec)){
  Plots.Compare(Compare_df = list(YFDP = YFDP_df, 
                                  Region = Region_df), 
                Header = Header_vec[Plots_Iter],
                SubsetHeader = NULL,
                X = "Treatment", SubsetX = NULL, 
                Y = "Method", SubsetY = NULL, 
                HMSCFilter = "diametre",
                FName = paste0("Across-Treatment_", Header_vec[Plots_Iter]),
                Dir = Dir.WithinScale)
} 

#### Methods across Scales --------------
Header_vec <- c("COCCUR")
Dir.CrossScale <- file.path(Dir.Exports, "Cross-Scale")
if(!dir.exists(Dir.CrossScale)){dir.create(Dir.CrossScale)}
for(Plots_Iter in 1:length(Header_vec)){
  Plots.Compare(Compare_df = list(YFDP = YFDP_df, 
                                  Region = Region_df), 
                Header = Header_vec[Plots_Iter],
                X = "Data Set", SubsetX = NULL, 
                Y = "Treatment", SubsetY = c("Pre-Fire", "Temperate Broadleaf & Mixed Forests"), 
                HMSCFilter = "diametre",
                FName = paste0("Cross-Scale_", Header_vec[Plots_Iter]),
                Dir = Dir.CrossScale)
}

#### COMPARING HEATMAPS ####

# HEATMAP -----


for(Counter in 1:nrow(Combines_df)){
  Plot_df <- Iter_df[Iter_df[, X] == Combines_df$X[Counter] & Iter_df[, Y] == Combines_df$Y[Counter], 1:3]
  if(nrow(Plot_df) == 0){next()}
  directed <- ifelse(Combines_df$Y[Counter] == "IF_REM" | 
                       Combines_df$X[Counter] == "IF_REM" | 
                       Header == "IF_REM", 
                     TRUE, FALSE)
  
}






col.pos <- colorRampPalette(c("palegreen", "lightgreen", "green", "forestgreen"))(1e4)
col.neg <- colorRampPalette(c("lightsalmon", "orangered", "red", "darkred"))(1e4)
seq.pos <- seq(0.0001, max(abs(Iter_df$Inference), na.rm = TRUE), length.out = 1e4)
seq.neg <- seq(-0.0001, -max(abs(Iter_df$Inference), na.rm = TRUE), length.out = 1e4)





mean_df <- reshape(Plot_df, idvar = "Partner1", timevar = "Partner2", direction = "wide")
rownames(mean_df) <- mean_df$Partner1
mean_df <- mean_df[ , -1]
colnames(mean_df) <- unlist(lapply(strsplit(colnames(mean_df), split = "[.]"), `[`, 2))
mean_df <- as.matrix(mean_df)


## Make vector of colors for values below threshold
rc1 <- col.neg[seq.neg >= min(Plot_df$Inference, na.rm = TRUE)]
## Make vector of colors for values above threshold
rc2 <- col.pos[seq.pos <= max(Plot_df$Inference, na.rm = TRUE)]
rampcols <- c(rc1, rc2)
rampcols[c(length(rc1), length(rc1)+1)] <- rgb(t(col2rgb("white")), maxColorValue=256) 

pheatmap(mean_df, display_numbers = TRUE, 
         number_format =  "%.3f",
         color = rampcols,
         main = "Actors (Columns) x Subject (Rows)",
         fontsize = 10, 
         cluster_rows = FALSE, 
         cluster_cols=FALSE,
         scale = "none"
)

ggplot() + 
  geom_bar(data=df2[which(df2$KK==10),], aes(k, value, fill = variable),stat = 'identity',position="dodge") +
  geom_bar(data=df2[which(df2$KK==30),], aes(k, value, fill = variable),stat = 'identity',position="dodge",alpha=0.5) +
  theme(legend.position = 'top')


rampcols <- c(rev(col.neg), col.pos)
rampcols[c(length(col.neg), length(col.pos)+1)] <- rgb(t(col2rgb("white")), maxColorValue=256) 

plot(c(rev(seq.neg), seq.pos), col = rampcols)


# some toy data
d <- data.frame(x = 1:length(rampcols), y = 1:length(rampcols))

# interpolate values from zero to y and create corresponding number of x values
vals <- lapply(d$y, function(y) seq(0, y, by = 0.01))
y <- unlist(vals)
mid <- rep(d$x, lengths(vals))
d2 <- data.frame(x = mid - 0.0001,
                 xend = mid + 0.0001,
                 y = y,
                 yend = y)

ggplot(data = d2, aes(x = x, xend = xend, y = y, yend = yend, color = y)) +
  geom_segment(size = 2) +
  scale_color_gradient2(low = "red", mid = "yellow", high = "green", 
                        midpoint = max(d2$y)/2) 

















































## INDIVIDUAL NETWORKS +++++++++++++++++++++++++++++++++++++++++++++++++++++
# this is where we will calculate whole-network topology metrics for comparison of entire networks









