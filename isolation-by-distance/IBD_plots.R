#author: Kevin Quinteros
#date: oct 6th 2020
#purpose: create scatter plots of genetic distance vs geographic distance data

#Set the working directory 
setwd("~/Projects/git_repositories/github/Fig_wasp_phylogeography/")

#Load libraries
library(tidyverse)
library(magrittr)
library(gridExtra)

#small funciton so I can intirate through list of mantel results 
c.ann <- function(x){
  tmp <- paste('R =',
               round(x$obs, digits=3),
               '|', 'p =',
               round(x$pvalue, 
                     digits=3), sep=' ')
  return(tmp)
  remove(tmp)
}

####################section for unsubsetted data ###############################
#Load genetic distance matrices
tmp <- load("Data/Gen_dist/Gen_dist_p50.RDA") %>% as.vector()

# replace diagonals with NA
dmat <- lapply(tmp, get) 
names(dmat) <- tmp
dmat <- lapply(dmat, function(x){diag(x)<- NA; x} )
remove(tmp)

#reformat them in vectors
dmat <- lapply(dmat, c)
#make vectors into tibbles and then bind them into a tibble
dmat <- lapply(dmat, tibble) %>% bind_cols
dmat %<>% filter_all(any_vars(!is.na(.)))
#rename tibble cols
colnames(dmat) <-c("Dgeo", "Dgen.1", "Dgen.2", "Dgen.3", "Dgen.4",
                   "Dgen.5", "Dgen.6", "Dgen.7", "Dgen.8")
#take log of distance data
dmat %<>% add_column(Dgeo.l = log(dmat$Dgeo))
#load in mantel test
load("output/Mantel_test/Mantel_p50.R")
#creating list of mantel results
result <- list(result.1,
               result.2,
               result.3,
               result.4,
               result.5,
               result.6,
               result.7,
               result.8)
#creatiung annotations for scatter plots based off mantel test results
#included R and P values
annot<- lapply(result, c.ann)
#creating basic scatter plot using ggplot2
sp.1 <- ggplot(dmat, aes(Dgeo,Dgen.1) ) + 
  geom_point() + 
  ylab("Nei's distance (1972)") +
  ggtitle(paste("Isolation by Distance", annot[1], sep = " "))  +
  xlab("Geographic Distance") +
  theme_classic() 

#ploting genetic distance vs geographic distance
sp.2 <- ggplot(dmat, aes(x=Dgeo, y=Dgen.2)) + 
  geom_point() + 
  ylab("Cavalli-Sforza and Edwards Chord distance") +
  ggtitle(paste("Isolation by Distance", annot[2], sep = " ")) +
  xlab("Geographic Distance") +
  theme_classic() 
 
#ploting genetic distance vs geographic distance
sp.3 <- ggplot(dmat, aes(x=Dgeo, y=Dgen.3)) +
  geom_point() + 
  ylab("Weir and Cockerham (1984) Fst") +
  ggtitle(paste("Isolation by Distance", annot[3], sep = " "))  +
  xlab("Geographic Distance") +
  theme_classic() 
  
#ploting genetic distance vs geographic distance
sp.4 <- ggplot(dmat, aes(x=Dgeo, y=Dgen.4)) + 
  geom_point() + 
  ylab("Rogers' distance")+
  ggtitle(paste("Isolation by Distance", annot[4], sep = " "))  +
  xlab("Geographic Distance") +
  theme_classic() 

#ploting genetic distance vs geographic distance
sp.5 <- ggplot(dmat, aes(x=Dgeo, y=Dgen.5)) +
  geom_point() + 
  ylab("Prevosti's distance")+
  ggtitle(paste("Isolation by Distance", annot[5], sep = " "))  +
  xlab("Geographic Distance") +
  theme_classic() 

#ploting genetic distance vs geographic distance
sp.6 <- ggplot(dmat, aes(x=Dgeo, y=Dgen.6)) +
  geom_point() + 
  ylab("FSTs according to Nei (1987)")+
  ggtitle(paste("Isolation by Distance", annot[6], sep = " "))  +
  xlab("Geographic Distance (euclidean)") +
  theme_classic() 

#ploting genetic distance vs geographic distance
sp.7 <- ggplot(dmat, aes(x=Dgeo.l, y=Dgen.7)) +
  geom_point() + 
  ylab("Linear Fst, WC Fst")+
  ggtitle(paste("Isolation by Distance", annot[7], sep = " ")) +
  xlab("Log Geographic Distance") +
  theme_classic()  

#ploting genetic distance vs geographic distance
sp.8 <- ggplot(dmat, aes(x=Dgeo.l, y=Dgen.8)) +
  geom_point() + 
  ylab("Linear Fst, Nei (1987) Fst")+
  ggtitle(paste("Isolation by Distance", annot[8], sep = " "))  +
  xlab("Log Geographic Distance") +
  theme_classic()  
#save plots to file in a grid pattern. 
pdf(file = "output/IBD/IBD_plots.pdf",   # The directory you want to save the file in
    width = 16.5, # The width of the plot in inches
    height = 12) # The height of the plot in inches
grid.arrange(sp.1, sp.2, sp.3, sp.4, sp.5, sp.6, sp.7, sp.8, nrow = 2 )
dev.off()

###################Subsetted  Dataset: Oaxaca and Morelos Removed###############
#Load genetic distance matrices
tmp <- load("Data/Gen_dist/Gen_dist_p50_nm.RDA") %>% as.vector()
# replace diagonals with NA
dmat <- lapply(tmp, get) 
names(dmat) <- tmp
dmat <- lapply(dmat, function(x){diag(x)<- NA; x} )
remove(tmp)

#reformat them in vectors
dmat <- lapply(dmat, c)
#make vectors into tibbles and then bind them into a tibble
dmat <- lapply(dmat, tibble) %>% bind_cols
dmat %<>% filter_all(any_vars(!is.na(.)))
#rename tibble cols
colnames(dmat) <-c("Dgeo", "Dgen.1", "Dgen.2", "Dgen.3", "Dgen.4",
                   "Dgen.5", "Dgen.6", "Dgen.7", "Dgen.8")
#take log of distance data
dmat %<>% add_column(Dgeo.l = log(dmat$Dgeo))
#load in mantel test
load("output/Mantel_test/Mantel_p50_nm.R")
#creating list of mantel results
result <- list(result.1nm,
               result.2nm,
               result.3nm,
               result.4nm,
               result.5nm,
               result.6nm,
               result.7nm,
               result.8nm)
#creatiung annotations for scatter plots based off mantel test results
#included R and P values
annot<- lapply(result, c.ann)
#creating basic scatter plot using ggplot2
sp.1 <- ggplot(dmat, aes(Dgeo,Dgen.1) ) + 
  geom_point() + 
  ylab("Nei's distance (1972)") +
  ggtitle(paste("Isolation by Distance", annot[1], sep = " "))  +
  xlab("Geographic Distance") +
  theme_classic() 

#ploting genetic distance vs geographic distance
sp.2 <- ggplot(dmat, aes(x=Dgeo, y=Dgen.2)) + 
  geom_point() + 
  ylab("Cavalli-Sforza and Edwards Chord distance") +
  ggtitle(paste("Isolation by Distance", annot[2], sep = " ")) +
  xlab("Geographic Distance") +
  theme_classic() 

#ploting genetic distance vs geographic distance
sp.3 <- ggplot(dmat, aes(x=Dgeo, y=Dgen.3)) +
  geom_point() + 
  ylab("Weir and Cockerham (1984) Fst") +
  ggtitle(paste("Isolation by Distance", annot[3], sep = " "))  +
  xlab("Geographic Distance") +
  theme_classic() 

#ploting genetic distance vs geographic distance
sp.4 <- ggplot(dmat, aes(x=Dgeo, y=Dgen.4)) + 
  geom_point() + 
  ylab("Rogers' distance")+
  ggtitle(paste("Isolation by Distance", annot[4], sep = " "))  +
  xlab("Geographic Distance") +
  theme_classic() 

#ploting genetic distance vs geographic distance
sp.5 <- ggplot(dmat, aes(x=Dgeo, y=Dgen.5)) +
  geom_point() + 
  ylab("Prevosti's distance")+
  ggtitle(paste("Isolation by Distance", annot[5], sep = " "))  +
  xlab("Geographic Distance") +
  theme_classic() 

#ploting genetic distance vs geographic distance
sp.6 <- ggplot(dmat, aes(x=Dgeo, y=Dgen.6)) +
  geom_point() + 
  ylab("FSTs according to Nei (1987)")+
  ggtitle(paste("Isolation by Distance", annot[6], sep = " "))  +
  xlab("Geographic Distance (euclidean)") +
  theme_classic() 

#ploting genetic distance vs geographic distance
sp.7 <- ggplot(dmat, aes(x=Dgeo,l, y=Dgen.7)) +
  geom_point() + 
  ylab("Linear Fst, WC Fst")+
  ggtitle(paste("Isolation by Distance", annot[7], sep = " ")) +
  xlab("Log Geographic Distance") +
  theme_classic()  

#ploting genetic distance vs geographic distance
sp.8 <- ggplot(dmat, aes(x=Dgeo.l, y=Dgen.8)) +
  geom_point() + 
  ylab("Linear Fst, Nei (1987) Fst")+
  ggtitle(paste("Isolation by Distance", annot[8], sep = " "))  +
  xlab("Log Geographic Distance") +
  theme_classic() 

#save plots to file in a grid pattern. 
pdf(file = "output/IBD/IBD_plots_nm.pdf",   # The directory you want to save the file in
    width = 16.5, # The width of the plot in inches
    height = 12) # The height of the plot in inches
grid.arrange(sp.1, sp.2, sp.3, sp.4, sp.5, sp.6, sp.7, sp.8, nrow = 2 )
dev.off()


####################Subsetted Dataset: Baja only###############################
#Load genetic distance matrices
tmp <- load("Data/Gen_dist/Gen_dist_p50_b.RDA") %>% as.vector()

# replace diagonals with NA
dmat <- lapply(tmp, get) 
names(dmat) <- tmp
dmat <- lapply(dmat, function(x){diag(x)<- NA; x} )
remove(tmp)

#reformat them in vectors
dmat <- lapply(dmat, c)
#make vectors into tibbles and then bind them into a tibble
dmat <- lapply(dmat, tibble) %>% bind_cols
dmat %<>% filter_all(any_vars(!is.na(.)))
#rename tibble cols
colnames(dmat) <-c("Dgeo", "Dgen.1", "Dgen.2", "Dgen.3", "Dgen.4",
                   "Dgen.5", "Dgen.6", "Dgen.7", "Dgen.8")
#take log of distance data
dmat %<>% add_column(Dgeo.l = log(dmat$Dgeo))
#load in mantel test
load("output/Mantel_test/Mantel_p50_b.R")
#creating list of mantel results
result <- list(result.1b,
               result.2b,
               result.3b,
               result.4b,
               result.5b,
               result.6b,
               result.7b,
               result.8b)
#creatiung annotations for scatter plots based off mantel test results
#included R and P values
annot<- lapply(result, c.ann)
#creating basic scatter plot using ggplot2
sp.1 <- ggplot(dmat, aes(Dgeo,Dgen.1) ) + 
  geom_point() + 
  ylab("Nei's distance (1972)") +
  ggtitle(paste("Isolation by Distance", annot[1], sep = " "))  +
  xlab("Geographic Distance") +
  theme_classic() 

#ploting genetic distance vs geographic distance
sp.2 <- ggplot(dmat, aes(x=Dgeo, y=Dgen.2)) + 
  geom_point() + 
  ylab("Cavalli-Sforza and Edwards Chord distance") +
  ggtitle(paste("Isolation by Distance", annot[2], sep = " ")) +
  xlab("Geographic Distance") +
  theme_classic() 

#ploting genetic distance vs geographic distance
sp.3 <- ggplot(dmat, aes(x=Dgeo, y=Dgen.3)) +
  geom_point() + 
  ylab("Weir and Cockerham (1984) Fst") +
  ggtitle(paste("Isolation by Distance", annot[3], sep = " "))  +
  xlab("Geographic Distance") +
  theme_classic() 

#ploting genetic distance vs geographic distance
sp.4 <- ggplot(dmat, aes(x=Dgeo, y=Dgen.4)) + 
  geom_point() + 
  ylab("Rogers' distance")+
  ggtitle(paste("Isolation by Distance", annot[4], sep = " "))  +
  xlab("Geographic Distance") +
  theme_classic() 

#ploting genetic distance vs geographic distance
sp.5 <- ggplot(dmat, aes(x=Dgeo, y=Dgen.5)) +
  geom_point() + 
  ylab("Prevosti's distance")+
  ggtitle(paste("Isolation by Distance", annot[5], sep = " "))  +
  xlab("Geographic Distance") +
  theme_classic() 

#ploting genetic distance vs geographic distance
sp.6 <- ggplot(dmat, aes(x=Dgeo, y=Dgen.6)) +
  geom_point() + 
  ylab("FSTs according to Nei (1987)")+
  ggtitle(paste("Isolation by Distance", annot[6], sep = " "))  +
  xlab("Geographic Distance (euclidean)") +
  theme_classic() 

#ploting genetic distance vs geographic distance
sp.7 <- ggplot(dmat, aes(x=Dgeo.l, y=Dgen.7)) +
  geom_point() + 
  ylab("Linear Fst, WC Fst")+
  ggtitle(paste("Isolation by Distance", annot[7], sep = " ")) +
  xlab("Log Geographic Distance") +
  theme_classic()  

#ploting genetic distance vs geographic distance
sp.8 <- ggplot(dmat, aes(x=Dgeo.l, y=Dgen.8)) +
  geom_point() + 
  ylab("Linear Fst, Nei (1987) Fst")+
  ggtitle(paste("Isolation by Distance", annot[8], sep = " "))  +
  xlab("Log Geographic Distance") +
  theme_classic() 
#save plots to file in a grid pattern. 
pdf(file = "output/IBD/IBD_plots_b.pdf",   # The directory you want to save the file in
    width = 16.5, # The width of the plot in inches
    height = 12) # The height of the plot in inches
grid.arrange(sp.1, sp.2, sp.3, sp.4, sp.5, sp.6, sp.7, sp.8, nrow = 2 )
dev.off()

###############Subsetted Dataset: Northern inland populations###################
#Load genetic distance matrices
tmp <-load("Data/Gen_dist/Gen_dist_p50_nl.RDA") %>% as.vector()
# replace diagonals with NA
dmat <- lapply(tmp, get) 
names(dmat) <- tmp
dmat <- lapply(dmat, function(x){diag(x)<- NA; x} )
remove(tmp)

#reformat them in vectors
dmat <- lapply(dmat, c)
#make vectors into tibbles and then bind them into a tibble
dmat <- lapply(dmat, tibble) %>% bind_cols
dmat %<>% filter_all(any_vars(!is.na(.)))
#rename tibble cols
colnames(dmat) <-c("Dgeo", "Dgen.1", "Dgen.2", "Dgen.3", "Dgen.4",
                   "Dgen.5", "Dgen.6", "Dgen.7", "Dgen.8")
#take log of distance data
dmat %<>% add_column(Dgeo.l = log(dmat$Dgeo))
#load in mantel test
load("output/Mantel_test/Mantel_p50_nl.R")
#creating list of mantel results
result <- list(result.1nl,
               result.2nl,
               result.3nl,
               result.4nl,
               result.5nl,
               result.6nl,
               result.7nl,
               result.8nl)
#creatiung annotations for scatter plots based off mantel test results
#included R and P values
annot<- lapply(result, c.ann)
#creating basic scatter plot using ggplot2
sp.1 <- ggplot(dmat, aes(Dgeo,Dgen.1) ) + 
  geom_point() + 
  ylab("Nei's distance (1972)") +
  ggtitle(paste("Isolation by Distance", annot[1], sep = " "))  +
  xlab("Geographic Distance") +
  theme_classic() 

#ploting genetic distance vs geographic distance
sp.2 <- ggplot(dmat, aes(x=Dgeo, y=Dgen.2)) + 
  geom_point() + 
  ylab("Cavalli-Sforza and Edwards Chord distance") +
  ggtitle(paste("Isolation by Distance", annot[2], sep = " ")) +
  xlab("Geographic Distance") +
  theme_classic() 

#ploting genetic distance vs geographic distance
sp.3 <- ggplot(dmat, aes(x=Dgeo, y=Dgen.3)) +
  geom_point() + 
  ylab("Weir and Cockerham (1984) Fst") +
  ggtitle(paste("Isolation by Distance", annot[3], sep = " "))  +
  xlab("Geographic Distance") +
  theme_classic() 

#ploting genetic distance vs geographic distance
sp.4 <- ggplot(dmat, aes(x=Dgeo, y=Dgen.4)) + 
  geom_point() + 
  ylab("Rogers' distance")+
  ggtitle(paste("Isolation by Distance", annot[4], sep = " "))  +
  xlab("Geographic Distance") +
  theme_classic() 

#ploting genetic distance vs geographic distance
sp.5 <- ggplot(dmat, aes(x=Dgeo, y=Dgen.5)) +
  geom_point() + 
  ylab("Prevosti's distance")+
  ggtitle(paste("Isolation by Distance", annot[5], sep = " "))  +
  xlab("Geographic Distance") +
  theme_classic() 

#ploting genetic distance vs geographic distance
sp.6 <- ggplot(dmat, aes(x=Dgeo, y=Dgen.6)) +
  geom_point() + 
  ylab("FSTs according to Nei (1987)")+
  ggtitle(paste("Isolation by Distance", annot[6], sep = " "))  +
  xlab("Geographic Distance (euclidean)") +
  theme_classic() 

#ploting genetic distance vs geographic distance
sp.7 <- ggplot(dmat, aes(x=Dgeo.l, y=Dgen.7)) +
  geom_point() + 
  ylab("Linear Fst, WC Fst")+
  ggtitle(paste("Isolation by Distance", annot[7], sep = " ")) +
  xlab("Log Geographic Distance") +
  theme_classic()  

#ploting genetic distance vs geographic distance
sp.8 <- ggplot(dmat, aes(x=Dgeo.l, y=Dgen.8)) +
  geom_point() + 
  ylab("Linear Fst, Nei (1987) Fst")+
  ggtitle(paste("Isolation by Distance", annot[8], sep = " "))  +
  xlab("Log Geographic Distance") +
  theme_classic() 
#save plots to file in a grid pattern. 
pdf(file = "output/IBD/IBD_plots_nl.pdf",   # The directory you want to save the file in
    width = 16.5, # The width of the plot in inches
    height = 12) # The height of the plot in inches
grid.arrange(sp.1, sp.2, sp.3, sp.4, sp.5, sp.6, sp.7, sp.8, nrow = 2 )
dev.off()

##############Subsetted Dataset:Oaxaca and Morelos#############################
#Load genetic distance matrices
tmp <- load("Data/Gen_dist/Gen_dist_p50_m.RDA") %>% as.vector 
# replace diagonals with NA
dmat <- lapply(tmp, get) 
names(dmat) <- tmp
dmat <- lapply(dmat, function(x){diag(x)<- NA; x} )
remove(tmp)

#reformat them in vectors
dmat <- lapply(dmat, c)
#make vectors into tibbles and then bind them into a tibble
dmat <- lapply(dmat, tibble) %>% bind_cols
dmat %<>% filter_all(any_vars(!is.na(.)))
#rename tibble cols
colnames(dmat) <-c("Dgeo", "Dgen.1", "Dgen.2", "Dgen.3",
                   "Dgen.6", "Dgen.7", "Dgen.8")
#take log of distance data
dmat %<>% add_column(Dgeo.l = log(dmat$Dgeo))
#load in mantel test
load("output/Mantel_test/Mantel_p50_m.R")
#creating list of mantel results
result <- list(result.1m,
               result.2m,
               result.3m,
               result.6m,
               result.7m,
               result.8m)
#creatiung annotations for scatter plots based off mantel test results
#included R and P values
annot<- lapply(result, c.ann)
#creating basic scatter plot using ggplot2
sp.1 <- ggplot(dmat, aes(Dgeo,Dgen.1) ) + 
  geom_point() + 
  ylab("Nei's distance (1972)") +
  ggtitle(paste("Isolation by Distance", annot[1], sep = " "))  +
  xlab("Geographic Distance") +
  theme_classic() 

#ploting genetic distance vs geographic distance
sp.2 <- ggplot(dmat, aes(x=Dgeo, y=Dgen.2)) + 
  geom_point() + 
  ylab("Cavalli-Sforza and Edwards Chord distance") +
  ggtitle(paste("Isolation by Distance", annot[2], sep = " ")) +
  xlab("Geographic Distance") +
  theme_classic() 

#ploting genetic distance vs geographic distance
sp.3 <- ggplot(dmat, aes(x=Dgeo, y=Dgen.3)) +
  geom_point() + 
  ylab("Weir and Cockerham (1984) Fst") +
  ggtitle(paste("Isolation by Distance", annot[3], sep = " "))  +
  xlab("Geographic Distance") +
  theme_classic() 

#ploting genetic distance vs geographic distance
sp.6 <- ggplot(dmat, aes(x=Dgeo, y=Dgen.6)) +
  geom_point() + 
  ylab("FSTs according to Nei (1987)")+
  ggtitle(paste("Isolation by Distance", annot[4], sep = " "))  +
  xlab("Geographic Distance (euclidean)") +
  theme_classic() 

#ploting genetic distance vs geographic distance
sp.7 <- ggplot(dmat, aes(x=Dgeo.l , y=Dgen.7)) +
  geom_point() + 
  ylab("Linear Fst, WC Fst")+
  ggtitle(paste("Isolation by Distance", annot[5], sep = " ")) +
  xlab("Log Geographic Distance") +
  theme_classic()  

#ploting genetic distance vs geographic distance
sp.8 <- ggplot(dmat, aes(x=Dgeo.l, y=Dgen.8)) +
  geom_point() + 
  ylab("Linear Fst, Nei (1987) Fst")+
  ggtitle(paste("Isolation by Distance", annot[6], sep = " "))  +
  xlab("Log Geographic Distance") +
  theme_classic() 

#save plots to file in a grid pattern. 
pdf(file = "output/IBD/IBD_plots_m.pdf",   # The directory you want to save the file in
    width = 16.5, # The width of the plot in inches
    height = 12) # The height of the plot in inches
grid.arrange(sp.1, sp.2, sp.3, sp.6, sp.7, sp.8, nrow = 2 )
dev.off()
