#load required libraries
library(tidyverse)
library(adegenet)
library(hierfstat)
library(magrittr)

#set working directory#######
setwd("~/Projects/git_repositories/github/Fig_wasp_phylogeography/")

#read in locality data and format#####
sites <- read_csv("Data/Locality_Fpet.csv")
sites <- sites[-3,] #remove FW268 since it did not phase
sites <- sites[-1]
sites <- sites[-6]
sites %<>% distinct(Site, .keep_all = TRUE)

#subsetting of data####
moax <- sites %>%
  filter(State%in% c("OAX", "MOR"))
No_moax <- sites %>%
  filter(State!= "OAX") %>%
  filter(State != "MOR")
Baja <- sites %>% filter(State == "BAJA" )
Nor <- sites %>%
  filter(State != "OAX") %>%
  filter(State != "MOR") %>%
  filter(State != "BAJA")

#Calculate geographic distance (kilometers#########
Dgeo <- fields::rdist.earth(as.matrix(sites[3:2]), miles = F)
Dgeo.m <- fields::rdist.earth(as.matrix(moax[3:2]), miles = F)
Dgeo.nm <- fields::rdist.earth(as.matrix(No_moax[3:2]), miles = F)
Dgeo.b <- fields::rdist.earth(as.matrix(Baja[3:2]), miles = F)
Dgeo.n <- fields::rdist.earth(as.matrix(Nor[3:2]), miles = F)


###This section is used to convert genetic data into genetic distance matrix####
#read in genetic data
load("Data/Phased_data/UCE_50p/Method1/Phased_snps/phased_snp_50.Rda")

#remove P342 hybrid individual
phased_snp_50 <- phased_snp_50[rownames(phased_snp_50$tab)!="P342_petiolaris",]

#convert to heifstat object
hf_p50 <- genind2hierfstat(phased_snp_50, pop = phased_snp_50$pop)
#remvoe rows with more than 50% 
#hf_p50 <- hf_p50[which(rowMeans(!is.na(hf_p50)) > 0.5), ] 
#subset dataset
###########Calculate Genetic distance with Oaxaca and Morelos
hf_m <- hf_p50 %>% filter(pop %in% moax$Site) #subset hierfasta objects
hf_m <- droplevels(hf_m)
###########Calculate Genetic distance with No Oaxaca and Morelos
hf_nm <- hf_p50 %>% filter(pop %in% No_moax$Site)
hf_nm <- droplevels(hf_nm)
###########Calculate Genetic distance with baja
hf_b <- hf_p50 %>% filter(pop %in% Baja$Site)
hf_b <- droplevels(hf_b)
###########Calculate Genetic distance with Northern Inland Pops###########
hf_nl <- hf_p50 %>% filter(pop %in% Nor$Site)
hf_nl <- droplevels(hf_nl)


#calculating genetic distance using different methods
# see Takezaki and Nei (1996)
#Cavalli-Sforza and Edwards Chord Distance 
Dgen <- genet.dist(hf_p50, diploid = T, method = "Dch") %>% as.matrix()
Dgen.m <- genet.dist(hf_m, diploid = T, method = "Dch") %>% as.matrix()
Dgen.nm <- genet.dist(hf_nm, diploid = T, method = "Dch") %>% as.matrix()
Dgen.b <- genet.dist(hf_b, diploid = T, method = "Dch") %>% as.matrix()
Dgen.n <- genet.dist(hf_nl, diploid = T, method = "Dch") %>% as.matrix()

Cavalli <- list(Dgen,Dgen.m,Dgen.nm,Dgen.b,Dgen.n)
names(Cavalli) <- c("Dgen", "Dgen.m","Dgen.nm", "Dgen.b","Dgen.n")
save(Cavalli, file = "Data/Gen_dist/cavalli-sforza_and_edwards_genetic_dist.RDA")
