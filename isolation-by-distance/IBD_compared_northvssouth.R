#Author: Kevin Quinteros
#Date: Feb 8, 2022
#This code preforms a mantel test on genetic distances 
#this will be for the ficus petiolaris 
#Author: Kevin Quinteros
#I'm supposed to just compare North vs South populations 

#load libraries
library(tidyverse)
library(ade4)
library(adegenet)
library(magrittr)

#setworking directories
setwd("~/Projects/git_repositories/github/Fig_wasp_phylogeography/")

#functions
`%nin%` = Negate(`%in%`)
####################load data###############
load("Data/Gen_dist/cavalli-sforza_and_edwards_genetic_dist.RDA") #genetic distance 1
load("Data/Gen_dist/Nei_83_genetic_dist.RDA") #genetic distance 2
load("Data/Gen_dist/geo_distance.RDA") #geographical distance
##################Load site information#########################
sites <- read_csv("Data/Locality_Fpet.csv")
sites <- sites[-3,] #remove FW268 since it did not phase
sites <- sites[-1]
sites <- sites[-6]
sites %<>% distinct(Site, .keep_all = TRUE)

#################remove within group comparisons##################
#give names col and row names to matrices
rownames(Cavalli$Dgen) <- sites$Site
colnames(Cavalli$Dgen) <- sites$Site
#Nei's distances 
rownames(Nei.83$Dgen) <- sites$Site
colnames(Nei.83$Dgen) <- sites$Site
#geographical distances
rownames(geo_distance$Dgeo) <- sites$Site
colnames(geo_distance$Dgeo) <- sites$Site

#sites that we can se
south <- sites[sites$Region == "South",] %>% 
  select(Site) 
north <- sites[sites$Region != "South",] %>%
  select(Site) 
south <- south$Site
north <- north$Site

#compare only the north and south regions 
Cavalli$Dgen <- Cavalli$Dgen[,colnames(Cavalli$Dgen) %nin% south]
Cavalli$Dgen <- Cavalli$Dgen[rownames(Cavalli$Dgen) %nin% north,]

Nei.83$Dgen <- Nei.83$Dgen[,colnames(Nei.83$Dgen) %nin% south]
Nei.83$Dgen <- Nei.83$Dgen[rownames(Nei.83$Dgen) %nin% north,]

geo_distance$Dgeo <- geo_distance$Dgeo[,colnames(geo_distance$Dgeo) %nin% south]
geo_distance$Dgeo <- geo_distance$Dgeo[rownames(geo_distance$Dgeo) %nin% north,]
###################perform mantel test##########################
geo_dist <- geo_distance$Dgeo
nei_dist <- Nei.83$Dgen
cav_dist <- Cavalli$Dgen
save(geo_dist,nei_dist, cav_dist, file = "Data/Gen_dist/pego_editted_data_northvssouth.RDA")

results.cav <- mantel.randtest(as.dist(cav_dist),as.dist(geo_dist), nrepet = 999)
result.nei <- mantel.randtest(as.dist(nei_dist),as.dist(geo_dist), nrepet = 999)
save(result.nei,results.cav, file = "output/Mantel_test/pego_editted_data_northvssouth_results.RDA")

#######################Create annotations for ggplot##############################

tmp.cav <- paste('R =',
                  round(results.cav$obs, digits=3),
                  '|', 'p =',
                  round(results.cav$pvalue, 
                        digits=3), sep=' ')
tmp.nei <- paste('R =',
                 round(result.nei$obs, digits=3),
                 '|', 'p =',
                 round(result.nei$pvalue, 
                       digits=3), sep=' ')
###########################Create dataframes#################################
#these dataframes will be used for plotting
dmat.cav <- cbind(as.vector(cav_dist), as.vector(geo_dist)) %>% as.data.frame()
dmat.cav <- dmat.cav[rowSums(dmat.cav[])>1,]

dmat.nei <- cbind(as.vector(nei_dist), as.vector(geo_dist)) %>% as.data.frame()
dmat.nei <- dmat.nei[rowSums(dmat.nei[])>1,]

###########################Create ggplots####################################
library(cowplot)
library(ggpubr)
p.cav <- ggplot(dmat.cav,aes(x = V2, y = V1)) + 
  geom_point(size = 2) + 
  ylab("Cavalli-Sforza and Edwards Chord distance") +
  xlab("Geographic Distance (km)") +
  ggtitle("Pegoscapus: North-Central Mexico \t
          vs South Mexico ")  +
  theme_cowplot() +
  theme(aspect.ratio = 1,
        plot.title = element_text(size = 11),
        axis.text = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12))

p.cav.lm <- p.cav + geom_smooth(method = lm, formula = y~x, se = F) +
  stat_regline_equation( label.x=1000, label.y=0.08, size = 4)+
  annotate( 'text',y = 0.09, x = 1000, label = tmp.cav, size=4)

plot(p.cav.lm)
ggsave(p.cav.lm, filename = "output/final_figures/Figure_5/IBD_cav_compare_northvssouth_lm.pdf", device = "pdf", 
       width = 4.5, height = 4, units = "in")

ggsave(p.cav, filename = "output/final_figures/Figure_5/IBD_cav_compare_northvssouth.pdf",
       device = "pdf", width = 4, height = 4, units = "in")


p.nei <- ggplot(dmat.nei,aes(x = V2, y = V1)) + 
  geom_point(size = 2) + 
  ylab("Nei's genetic distance (1983)") +
  xlab("Geographic Distance (km)") +
  ggtitle("Pegoscapus: North-Central Mexico \t
          vs South Mexico ")  +
  theme_cowplot() +
  theme(aspect.ratio = 1,
        plot.title = element_text(size = 11),
        axis.text = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12))

p.nei.lm <- p.nei + geom_smooth(method = lm, formula = y~x, se = F) +
  stat_regline_equation( label.x=500, label.y=0.08, size = 4)+
  annotate( 'text',y = 0.105, x = 1000, label = tmp.nei, size=4)

plot(p.nei.lm)
ggsave(p.nei.lm, filename = "output/final_figures/Figure_5/IBD_nei_compare_northvssouth_lm.pdf",
       device = "pdf", width = 4, height = 4, units = "in")

ggsave(p.nei, filename = "output/final_figures/Figure_5/IBD_nei_compare_northvssouth.pdf",
       device = "pdf", width = 4, height = 4, units = "in")
################remove all sites not in central, MOR and NOR OAX################
load("Data/Gen_dist/cavalli-sforza_and_edwards_genetic_dist.RDA") #genetic distance 1
load("Data/Gen_dist/Nei_83_genetic_dist.RDA") #genetic distance 2
load("Data/Gen_dist/geo_distance.RDA") #geographical distance
#select sites I want
sites <- sites[sites$fine_region %in% c("JAL","ZAC","MOR","OAX_NOR"),]
central <- sites$Site[sites$fine_region %in% c("JAL","ZAC")]
south_nor <- sites$Site[sites$fine_region %in% c("MOR","OAX_NOR")]

#remove the sites I don't want
view(Cavalli$Dgen)
Cavalli$Dgen <- Cavalli$Dgen[,colnames(Cavalli$Dgen) %in% central]
Cavalli$Dgen <- Cavalli$Dgen[rownames(Cavalli$Dgen) %in% south_nor,]

Nei.83$Dgen <- Nei.83$Dgen[,colnames(Nei.83$Dgen) %in% central]
Nei.83$Dgen <- Nei.83$Dgen[rownames(Nei.83$Dgen) %in% south_nor,]

geo_distance$Dgeo <- geo_distance$Dgeo[,colnames(geo_distance$Dgeo) %in% central]
geo_distance$Dgeo <- geo_distance$Dgeo[rownames(geo_distance$Dgeo) %in% south_nor,]

###################perform mantel test##########################
geo_dist <- geo_distance$Dgeo
nei_dist <- Nei.83$Dgen
cav_dist <- Cavalli$Dgen
save(geo_dist,nei_dist, cav_dist, file = "Data/Gen_dist/pego_editted_data_centralvssout_norh.RDA")

results.cav <- mantel.randtest(as.dist(cav_dist),as.dist(geo_dist), nrepet = 999)
result.nei <- mantel.randtest(as.dist(nei_dist),as.dist(geo_dist), nrepet = 999)
save(result.nei,results.cav, file = "output/Mantel_test/pego_editted_data_centralvssout_norh.RDA")

#######################Create annotations for ggplot##############################

tmp.cav <- paste('R =',
                 round(results.cav$obs, digits=3),
                 '|', 'p =',
                 round(results.cav$pvalue, 
                       digits=3), sep=' ')
tmp.nei <- paste('R =',
                 round(result.nei$obs, digits=3),
                 '|', 'p =',
                 round(result.nei$pvalue, 
                       digits=3), sep=' ')
###########################Create dataframes#################################
#these dataframes will be used for plotting
dmat.cav <- cbind(as.vector(cav_dist), as.vector(geo_dist)) %>% as.data.frame()
dmat.cav <- dmat.cav[rowSums(dmat.cav[])>1,]

dmat.nei <- cbind(as.vector(nei_dist), as.vector(geo_dist)) %>% as.data.frame()
dmat.nei <- dmat.nei[rowSums(dmat.nei[])>1,]

p.nei <- ggplot(dmat.nei,aes(x = V2, y = V1)) + 
  geom_point(size = 2) + 
  ylab("Nei's genetic distance (1983)") +
  xlab("Geographic Distance (km)") +
  ggtitle("Pegoscapus: Central Mexico \t
          vs South-North Mexico ")  +
  geom_smooth(method = lm, formula = y~x, se = F) +
  stat_regline_equation( label.x=500, label.y=0.08, size = 4)+
  annotate( 'text',y = 0.085, x = 525, label = tmp.nei, size=4)+
  theme_cowplot() +
  theme(aspect.ratio = 1,
        plot.title = element_text(size = 11),
        axis.text = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12))
plot(p.nei)

ggsave(p.cav, filename = "output/final_figures/Figure_5/IBD_cav_compare_centralvssouth-north.pdf",
       device = "pdf", width = 4, height = 4, units = "in")

p.cav <- ggplot(dmat.cav,aes(x = V2, y = V1)) + 
  geom_point(size = 2) + 
  ylab("Cavalli-Sforza and Edwards Chord distance") +
  xlab("Geographic Distance (km)") +
  ggtitle("Pegoscapus: Central Mexico \t
          vs South_north Mexico ")  +
  geom_smooth(method = lm, formula = y~x, se = F) +
  stat_regline_equation( label.x=480, label.y=0.085, size = 4)+
  annotate( 'text',y = 0.09, x = 550, label = tmp.cav, size=4)+
  theme_cowplot() +
  theme(aspect.ratio = 1,
        plot.title = element_text(size = 11),
        axis.text = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12))
plot(p.cav)
ggsave(p.nei, filename = "output/final_figures/Figure_5/IBD_nei_compare_centralvssouth-north.pdf",
       device = "pdf", width = 4, height = 4, units = "in")

ggsave(p.cav, filename = "output/final_figures/Figure_5/IBD_cav_compare_centralvssouth-north.pdf",
       device = "pdf", width = 4, height = 4, units = "in")
