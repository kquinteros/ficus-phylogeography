#Author: Kevin Quinteros
#Date: Feb 9, 2022
#This code preforms a mantel test on genetic distances 

#load libraries
library(tidyverse)
library(ade4)
library(adegenet)
library(ggpubr)

#setworking directories
setwd("~/Projects/git_repositories/github/Fig_wasp_phylogeography/")
########################Load Dataset############################################
load("Data/Gen_dist/cavalli-sforza_and_edwards_genetic_dist.RDA") #genetic distance 1
load("Data/Gen_dist/Nei_83_genetic_dist.RDA") #genetic distance 2
load("Data/Gen_dist/geo_distance.RDA") #geographical distance 

########################preform mantel test#####################################
#mantel for Northern and Southern population using Cavalli-sforza and edwards chord distance
result.1.nm = mantel.randtest(as.dist(Cavalli$Dgen.nm), as.dist(geo_distance$Dgeo.nm), nrepet = 999)
result.1.m = mantel.randtest(as.dist(Cavalli$Dgen.m), as.dist(geo_distance$Dgeo.m), nrepet = 999)
#mantel for Northern and Southern population using Nei's 1983
result.2.nm = mantel.randtest(as.dist(Nei.83$Dgen.nm), as.dist(geo_distance$Dgeo.nm), nrepet = 999)
result.2.m = mantel.randtest(as.dist(Nei.83$Dgen.m), as.dist(geo_distance$Dgeo.m), nrepet = 999)

save(result.1.nm,result.1.m, result.2.nm, result.2.m, file = "output/Mantel_test/pego_mantel_test.RDA")
#######################Create annotations for ggplot##############################
tmp.1.nm <- paste('R =',
             round(result.1.nm$obs, digits=3),
             '|', 'p =',
             round(result.1.nm$pvalue, 
                   digits=3), sep=' ')
tmp.1.m <- paste('R =',
                  round(result.1.m$obs, digits=3),
                  '|', 'p =',
                  round(result.1.m$pvalue, 
                        digits=3), sep=' ')
tmp.2.nm <- paste('R =',
                  round(result.2.nm$obs, digits=3),
                  '|', 'p =',
                  round(result.2.nm$pvalue, 
                        digits=3), sep=' ') 
tmp.2.m <- paste('R =',
                  round(result.2.m$obs, digits=3),
                  '|', 'p =',
                  round(result.2.m$pvalue, 
                        digits=3), sep=' ') 
###########################Create dataframes#################################
#these dataframes will be used for plotting
dmat.1.nm <- cbind(as.vector(Cavalli$Dgen.nm), as.vector(geo_distance$Dgeo.nm)) %>% as.data.frame()
dmat.1.nm <- dmat.1.nm[rowSums(dmat.1.nm[])>1,]
dmat.1.m <- cbind(as.vector(Cavalli$Dgen.m), as.vector(geo_distance$Dgeo.m)) %>% as.data.frame()
dmat.1.m <- dmat.1.m[rowSums(dmat.1.m[])>1,]

dmat.2.nm <- cbind(as.vector(Nei.83$Dgen.nm), as.vector(geo_distance$Dgeo.nm)) %>% as.data.frame()
dmat.2.nm <- dmat.2.nm[rowSums(dmat.2.nm[])>1,]
dmat.2.m <- cbind(as.vector(Nei.83$Dgen.m), as.vector(geo_distance$Dgeo.m)) %>% as.data.frame()
dmat.2.m <- dmat.2.m[rowSums(dmat.2.m[])>1,]

###########################Create ggplots####################################
library(cowplot)
p.1.nm <- ggplot(dmat.1.nm,aes(x = V2, y = V1)) + 
        geom_point(size = 2) + 
        ylab("Cavalli-Sforza and Edwards Chord distance") +
        xlab("Geographic Distance (km)") +
        ggtitle("Pegoscapus: Northern and Central Mexico")  +
        stat_smooth(method = lm, formula = y~x, se = F) + 
        stat_regline_equation( label.x=900, label.y=0.010, size = 3)+
        annotate( 'text',y = 0.011, x = 1300, label = tmp.1.nm, size=3)+
        theme_cowplot() +
        theme(aspect.ratio = 1,
              plot.title = element_text(size = 11),
              axis.text = element_text(size = 12),
              axis.title.x = element_text(size = 12),
              axis.title.y = element_text(size = 12))
plot(p.1.nm)

ggsave(p.1.nm, filename = "output/final_figures/Figure_5/IBD_cavalli_no_moax.pdf",
       device = "pdf", width = 4, height = 4, units = "in")

p.1.m <- ggplot(dmat.1.m,aes(x = V2, y = V1)) + 
        geom_point(size = 2) + 
        ylab("Cavalli-Sforza and Edwards Chord distance") +
        xlab("Geographic Distance (km)") +
        ggtitle("Pegoscapus: Morelos and Oaxaca")  +
        stat_smooth(method = lm, formula = y~x, se = F) + 
        stat_regline_equation( label.x=200, label.y=0.018, size = 4)+
        annotate( 'text',y = 0.021, x = 300, label = tmp.1.m, size=4)+
        theme_cowplot() +
        theme(aspect.ratio = 1,
              plot.title = element_text(size = 11),
              axis.text = element_text(size = 12),
              axis.title.x = element_text(size = 12),
              axis.title.y = element_text(size = 12))
plot(p.1.m)
ggsave(p.1.m, filename = "output/final_figures/Figure_5/IBD_cavalli_moax.pdf",
       device = "pdf", width = 4, height = 4, units = "in")

p.2.nm <- ggplot(dmat.2.nm,aes(x = V2, y = V1)) + 
        geom_point(size = 2) + 
        ylab("Nei's genetic distance (1983)") +
        xlab("Geographic Distance (km)") +
        ggtitle("Pegoscapus: Northern and Central Mexico")  +
        stat_smooth(method = lm, formula = y~x, se = F) + 
        stat_regline_equation( label.x=750, label.y=0.0050, size = 4)+
        annotate( 'text',y = 0.006, x = 900, label = tmp.2.nm, size=4)+
        theme_cowplot() +
        theme(aspect.ratio = 1,
              plot.title = element_text(size = 11),
              axis.text = element_text(size = 12),
              axis.title.x = element_text(size = 12),
              axis.title.y = element_text(size = 12))
plot(p.2.nm)
ggsave(p.2.nm, filename = "output/final_figures/Figure_5/IBD_Nei_83_no_moax.pdf",
       device = "pdf", width = 4, height = 4, units = "in")

p.2.m <- ggplot(dmat.2.m,aes(x = V2, y = V1)) + 
        geom_point(size = 2) + 
        ylab("Nei's genetic distance (1983)") +
        xlab("Geographic Distance (km)") +
        ggtitle("Pegoscapus: Morelos and Oaxaca")  +
        stat_smooth(method = lm, formula = y~x, se = F) + 
        stat_regline_equation( label.x=200, label.y=0.018, size = 4)+
        annotate( 'text',y = 0.021, x = 300, label = tmp.2.m, size=4)+
        theme_cowplot() +
        theme(aspect.ratio = 1,
              plot.title = element_text(size = 11),
              axis.text = element_text(size = 12),
              axis.title.x = element_text(size = 12),
              axis.title.y = element_text(size = 12))
plot(p.2.m)
ggsave(p.2.m, filename = "output/final_figures/Figure_5/IBD_Nei_83_moax.pdf",
       device = "pdf", width = 4, height = 4, units = "in")
