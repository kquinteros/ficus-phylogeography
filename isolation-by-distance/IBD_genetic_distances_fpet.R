#Author: Kevin Quinteros
#Date: Feb 9, 2022
#This code preforms a mantel test on genetic distances 
#this will be for the ficus petiolaris 
#Author: Kevin Quinteros
#Date: Feb 9, 2022
#This code preforms a mantel test on genetic distances 
#this will be for the ficus petiolaris 
#load libraries
library(tidyverse)
library(ade4)
library(adegenet)

#setworking directories
setwd("~/Projects/git_repositories/github/Fig_wasp_phylogeography/")

#####################import structure data (unlinked)###########################
data <- read.structure("Data/Fpet_tree/Second_Finn_dataset/output/sub_fpet-s2_trim25-s3_clu_th90-s7_min125_outfiles/sub_fpet-s2_trim25-s3_clu_th90-s7_min125.stru",
                       n.ind = 250,
                       n.loc = 1192,
                       onerowperind = F,
                       NA.char = "-9",
                       row.marknames = 0,
                       col.lab = 1,
                       col.pop = 0, quiet = T)

x <- read.delim("Data/STRUCTURE_FILES/Tree/Fpet_tree_subsetted.str", sep = " ",header = F)
sub <- unique(x$V1)
remove(x)
###################subset the data set according to structure subset############
sites <- read.csv("Data/Fpet_tree/Second_Finn_dataset/fpet_tree_locality.csv")
sites <- apply(sites,2,function(x) gsub("-Tetas", "T", x)) %>% as.data.frame()
pop(data)  <- sites$Site
strata(data) <- sites
rownames(data$tab) <- sites$Sample_ID 
#subset
data.sub <- data[rownames(data$tab) %in% sub,] 
################create the subset datasets####################################
No_south <- sites %>%
  filter(Phylogroups!= "SOUTH") %>%
  filter(Phylogroups!= "Out") %>%
  filter(Phylogroups!= "BAJA")
No_south %<>% distinct(Site, .keep_all = TRUE)
Baja <- sites %>% filter(State == "BAJA")
Baja %<>% distinct(Site, .keep_all = TRUE)

#########################convert to hierfstat object#######################3####
hf <- genind2hierfstat(data, pop = data$pop)
hf.sub <- genind2hierfstat(data.sub, pop = data.sub$pop)

#remove all samples that are not BAJA
hf <- hf %>% filter(pop %in% Baja$Site) #subset hierfasta objects
hf <- droplevels(hf)

#remove all samples that are from Oaxaca 
hf.sub <- hf.sub %>% filter(pop %in% No_south$Site)
hf.sub <- droplevels(hf.sub)



#########################calculate genetic distances###########################
#calculating genetic distance using different methods
# see Takezaki and Nei (1996)
#Cavalli-Sforza and Edwards Chord Distance 
Dgen <- genet.dist(hf, diploid = T, method = "Dch") %>% as.matrix()
Dgen.sub <- genet.dist(hf.sub, diploid = T, method = "Dch") %>% as.matrix()
Cavalli.fpet <- list(Dgen, Dgen.sub)
names(Cavalli.fpet) <- c("Dgen", "Dgen.sub")
save(Cavalli.fpet, file = "Data/Gen_dist/fpet_cavalli-sforza_edwards_genetic_dist.RDA")
#Nei's genetic distance  
Dgen <- genet.dist(hf, diploid = T, method = "Da") %>% as.matrix()
Dgen.sub <- genet.dist(hf.sub, diploid = T, method = "Da") %>% as.matrix()
Nei.fpet <- list(Dgen,Dgen.sub)
names(Nei.fpet) <- c("Dgen", "Dgen.sub")
save(Nei.fpet, file = "Data/Gen_dist/fpet_Nei_83_genetic_dist.RDA")
#######################calculate geographic distance############################
Baja$lat <- as.numeric(Baja$lat)
Baja$lon <- as.numeric(Baja$lon)
No_south$lat <- as.numeric(No_south$lat)
No_south$lon <- as.numeric(No_south$lon)

Dgeo.b <- fields::rdist.earth(as.matrix(Baja[3:2]), miles = F)
Dgeo.ns <- fields::rdist.earth(as.matrix(No_south[3:2]), miles = F)
geo_dist <- list(Dgeo.b, Dgeo.ns)
names(geo_dist) <- c("Dgeo.b", "Dgeo.ns")
save(geo_dist, file = "Data/Gen_dist/fpet_geo_distance.RDA")
