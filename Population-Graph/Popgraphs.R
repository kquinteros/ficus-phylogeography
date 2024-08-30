# ---------------------------------------------------------------------------- #
# Author: Kevin Quinteros
# Date: October 12, 2023
# ---------------------------------------------------------------------------- #
# Purpose: This R code aims to analyze and visualize the distribution of 
# Pegoscapus population graphs as described by (Dyer & Nason 2004).
#
#
# 1. Plot Pegoscapus populations graph using the popgraph package
#
# 2. Include geographical coordinates along with elevation data
#
# 3. Evaluate population connectivity of Pegoscapus population graphs using 
#.   permutation
#
# 4. To compute and visualize a histogram of permutation and observed
#
# ----------------------------------------------------------------------------#

# ------------------------------LOAD LIBRARIES--------------------------------
#load library
library(popgraph)
library(igraph)
library(gstudio)
library(tidyverse)
library(data.table)
library(raster)
library(rgdal)
library(elevatr)

# ------------------------------SET ENVIRONMENT--------------------------------
# set workding directory
# setwd("~/Projects/Fig_wasp_phylogeography/") 

# ------------------------------READ IN DATA-----------------------------------
#Read in Genetic Data
pego.df <- read.csv("hf_pego_sep.csv") 
pego <- read_population("hf_pego_sep.csv", type="separated",
                        locus.columns = 2:1416,
                        phased= TRUE, sep = ",",
                        header = TRUE)

#fix population colum
pego$Pop <- pego.df$Pop

#Read in geospatial Data
geo <- read_csv("locality_pego.csv")

# --------------------------DATA MANIPULATION-----------------------------------
#group some site together,since they're in the same cluster based off PCA
#and STRUCTURE or they're representative of a general geographical area
which(pego$Pop == "214") #group with 209
which(pego$Pop == "229") #group with 226
which(pego$Pop == "228") #group with 226
which(pego$Pop == "100T") #group with 100
which(pego$Pop == "220") #group with 219
which(pego$Pop == "201") #group with 95
which(pego$Pop == "96") #group with 70
which(pego$Pop == "39") #group with 70
which(pego$Pop == "215") #group with 217
which(pego$Pop == "231") #group with 230
which(pego$Pop == "172") #group with 158

#replacing values
pego$Pop[31] <- "209"
pego$Pop[c(83,84,85)] <- "226"
pego$Pop[c(15,32)] <- "100"
pego$Pop[60] <- "219" 
pego$Pop[c(22,23)] <- "95"
pego$Pop[c(2,9,10)] <- "70"
pego$Pop[c(24,25)] <- "217"
pego$Pop[c(91,92,93)] <- "230"
pego$Pop[c(20,21)] <- "158"

# Match rows in data_frame1 with data_frame2 based on the 'ID' column
match_indices <- match(pego.df$Ind, geo$Sample)
# Print the matched indices
cat("Matched indices:\n",match_indices )

#Add Geographical data to GENEPOP object
pego <- cbind( geo[match_indices,3:6], pego)

#Remove column
pego$Ind <- NULL

#convert data
pego.pg <- to_mv(pego) 

#refine geo dataframe
geo <- cbind(geo[match_indices,3:6],pego$Pop,geo$fine_region[match_indices])
names(geo) <- c("Latitude", "Longitude","Elevation", "State","Population","Region")

# Convert to data.table
setDT(geo)

remove(avg_elevation)
#find avearge of elevt
avg_elevation <- geo %>%
  group_by(Population) %>%
  summarise(Average_Elevation = mean(Elevation))


# Average Latitude, Longitude and Elevation
geo <- geo %>%
  group_by(Population) %>%
  mutate(avg_elevation = mean(Elevation)) 
geo <- geo %>%
  group_by(Population) %>%
  mutate(avg_latitude = mean(Latitude)) 
geo <- geo %>%
  group_by(Population) %>%
  mutate(avg_longitude = mean(Longitude)) 

#Code above doesn't word for every group, try by just taking average
avg_latitude <- tapply(geo$Latitude, geo$Population, mean)
avg_longitude <- tapply(geo$Longitude, geo$Population, mean)
geo.short <- cbind(avg_elevation,avg_latitude,avg_longitude)
geo.short$Population <- factor(geo.short$Population)

#Add factors
geo$Region <- factor(geo$Region)
geo$Population <- factor(geo$Population)
geo$State <- factor(geo$State)

#remove objects
remove(avg_latitude,avg_longitude, match_indices)
# --------------CODING POPULATION GRAPH-----------------------------------------
graph <- popgraph(x=pego.pg, groups = geo$Population)
graph <- decorate_graph(graph, geo.short, stratum = "Population")


#Adding Geo and elevation data
V(graph)$Latitude <- geo.short$avg_latitude
V(graph)$Longitude <- geo.short$avg_longitude
V(graph)$Elevation <- geo.short$Average_Elevation

# Download elevation data
mexico.alt <- raster::getData("alt", country = "mexico", path = tempdir(), mask = TRUE)

#Convert to Information
graph.nodes <- to_SpatialPoints(graph)
graph.nodes
graph.edges <- to_SpatialLines(graph)
graph.edges

#plot on top of raster map 
plot(mexico.alt, terrain.colors(255))
plot( graph.edges, add=TRUE, col="#555555" )
plot( graph.nodes, add=TRUE, col="black", cex=1.5 )
plot( graph.nodes, add=TRUE, col="red", pch=16, cex=1.5 )

# ----------------Extracting Spatial Data Using Population Graphs---------------
df.nodes <- data.frame(Pop=V(graph)$name, Latitude=V(graph)$Latitude, Longitude=V(graph)$Longitude)
df.nodes$Elevation <- geo.short$Average_Elevation
summary(df.nodes)

#Extracting Data Along Popgraph Edges
edges <- to_SpatialLines(graph)
proj4string(edges) <- CRS( proj4string( mexico.alt ))
plot( mexico.alt, legend=FALSE)
plot(edges,add=TRUE)

#determine which edge is the longest
edge_lengths <- SpatialLinesLengths( edges )
longest <- sort( edge_lengths,decreasing = TRUE )[1]
longest

#is there a preference for environment 
allpops <- V(graph)$name

#make an adjacency matrix connecting all pairs of populations
A <- matrix(1,nrow=length(allpops),ncol=length(allpops))
diag(A) <- 0
rownames(A) <- colnames(A) <- allpops
saturated_graph <- graph.adjacency(A, mode = "undirected")
saturated_graph <- as.popgraph( saturated_graph )
plot(saturated_graph)
coords
#pull all the edges as SpatialLines objects
coords <- cbind(geo.short[,1],geo.short[,4:3])
names(coords) <- c("Population","Longitude","Latitude")
saturated_graph <- decorate_graph( saturated_graph, coords, stratum = "Population")
all_edges <- to_SpatialLines( saturated_graph )

#we can extract elevation data from the elevation raster.
edge_values <- raster::extract(mexico.alt, all_edges, fun=max, na.rm=TRUE, df=TRUE)
edge_names <- as_edgelist( saturated_graph )
edge_values$Nodes <- paste( edge_names[,1], edge_names[,2], sep="-")

#extract the edges that we observed in the original Population Graph
e <- as_edgelist( graph )
obs <- edge_values$MEX_msk_alt[ edge_values$Nodes %in% paste( e[,1], e[,2], sep="-") ]
mean(obs) #mean elevation

# ---------------------------------PERMUTATION----------------------------------
# permute the network a moderate number of times and take the values of permuted
#elevation to see if our observed are smaller than all potential elevations for 
#this specific network.
perm_elev <- rep(NA,999)
for( i in 1:length(perm_elev) ) {
  perm_graph <- randomize_graph( graph )
  e <- as_edgelist( perm_graph )
  perm_val <- edge_values$MEX_msk_alt[ edge_values$Nodes %in% paste( e[,1], e[,2], sep="-") ]
  perm_elev[i] <- mean(perm_val)
}

#Plot out results
library(cowplot)
df <- data.frame( Elevation=c(mean(obs),perm_elev), Category=c("Observed",rep("Permuted",999)))
dist <- ggplot( df, aes(x=Elevation,fill=Category)) +
  geom_histogram(stat="bin", bins=40) + 
  xlab("Elevation (m)") +
  ylab("Distribution of Permuted Elevations") +
  theme(
    axis.title.x = element_text(size = 8),   # Set x-axis label text size
    axis.title.y = element_text(size = 8),   # Set y-axis label text size
    axis.text.x = element_text(size = 6),    # Set x-axis tick text size
    axis.text.y = element_text(size = 6)     # Set y-axis tick text size
  ) +
  theme_cowplot()
plot(dist)
sum( mean(obs) >= perm_elev )

# Save to file 
ggsave(filename = "output/Popgraph/Dist-Permutation-Elevations.pdf", plot = dist, width = 6, height = 5)


