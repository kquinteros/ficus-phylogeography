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
#################################Load Data######################################
load("Data/Gen_dist/fpet_cavalli-sforza_edwards_genetic_dist.RDA")
load("Data/Gen_dist/fpet_geo_distance.RDA")
load("Data/Gen_dist/fpet_Nei_83_genetic_dist.RDA")
########################preform mantel test#####################################
#mantel for Northern and Southern population using Cavalli-sforza and edwards chord distance
result.1.ns = mantel.randtest(as.dist(Cavalli.fpet$Dgen.sub), as.dist(geo_dist$Dgeo.ns), nrepet = 999)
result.1 = mantel.randtest(as.dist(Cavalli.fpet$Dgen), as.dist(geo_dist$Dgeo.b), nrepet = 999)
#mantel for Northern and Southern population using Nei's 1983
result.2.ns = mantel.randtest(as.dist(Nei.fpet$Dgen.sub), as.dist(geo_dist$Dgeo.ns), nrepet = 999)
result.2 = mantel.randtest(as.dist(Nei.fpet$Dgen), as.dist(geo_dist$Dgeo.b), nrepet = 999)

save(result.1, result.1.ns, result.2, result.2.ns, file = "output/Mantel_test/fpet_mantel_test.RDA")


#######################Create annotations for ggplot##############################
tmp.1.ns <- paste('R =',
                  round(result.1.ns$obs, digits=3),
                  '|', 'p =',
                  round(result.1.ns$pvalue, 
                        digits=3), sep=' ')
tmp.1 <- paste('R =',
                 round(result.1$obs, digits=3),
                 '|', 'p =',
                 round(result.1$pvalue, 
                       digits=3), sep=' ')
tmp.2.ns <- paste('R =',
                  round(result.2.ns$obs, digits=3),
                  '|', 'p =',
                  round(result.2.ns$pvalue, 
                        digits=3), sep=' ') 
tmp.2 <- paste('R =',
                 round(result.2$obs, digits=3),
                 '|', 'p =',
                 round(result.2$pvalue, 
                       digits=3), sep=' ') 
###########################Create dataframes#################################
#these dataframes will be used for plotting
dmat.1.ns <- cbind(as.vector(Cavalli.fpet$Dgen.sub), as.vector(geo_dist$Dgeo.ns)) %>% as.data.frame()
dmat.1.ns <- dmat.1.ns[rowSums(dmat.1.ns[])>1,]
dmat.1 <- cbind(as.vector(Cavalli.fpet$Dgen), as.vector(geo_dist$Dgeo.b)) %>% as.data.frame()
dmat.1 <- dmat.1[rowSums(dmat.1[])>1,]

dmat.2.ns <- cbind(as.vector(Nei.fpet$Dgen.sub), as.vector(geo_dist$Dgeo.ns)) %>% as.data.frame()
dmat.2.ns <- dmat.2.ns[rowSums(dmat.2.ns[])>1,]
dmat.2 <- cbind(as.vector(Nei.fpet$Dgen), as.vector(geo_dist$Dgeo.b)) %>% as.data.frame()
dmat.2 <- dmat.2[rowSums(dmat.2[])>1,]

###########################Create ggplots####################################
library(cowplot)
library(ggpubr)
p.1.ns <- ggplot(dmat.1.ns,aes(x = V2, y = V1)) + 
  geom_point(size = 2) + 
  ylab("Cavalli-Sforza and Edwards Chord distance") +
  xlab("Geographic Distance (km)") +
  ggtitle("F. petiolaris: Sonora and Central Mexico")  +
  stat_smooth(method = lm, formula = y~x, se = F) + 
  stat_regline_equation( label.x=200, label.y=0.018, size = 4)+
  annotate( 'text',y = 0.021, x = 300, label = tmp.1.ns, size=4)+
  theme_cowplot() +
  theme(aspect.ratio = 1,
        plot.title = element_text(size = 11),
        axis.text = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12))
plot(p.1.ns)

ggsave(p.1.ns, filename = "output/final_figures/Figure_5/IBD_fpet_cavalli_no_south.pdf", device = "pdf", width = 4, height = 4, units = "in")

p.1 <- ggplot(dmat.1,aes(x = V2, y = V1)) + 
  geom_point(size = 2) + 
  ylab("Cavalli-Sforza and Edwards Chord distance") +
  xlab("Geographic Distance (km)") +
  ggtitle("F. petiolaris: Baja California")  +
  stat_smooth(method = lm, formula = y~x, se = F) + 
  stat_regline_equation( label.x=400, label.y=0.006, size = 4)+
  annotate( 'text',y = 0.005, x = 600, label = tmp.1, size=4)+
  theme_cowplot() +
  theme(aspect.ratio = 1,
        plot.title = element_text(size = 11),
        axis.text = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12))
plot(p.1 )

ggsave(p.1 , filename = "output/final_figures/Figure_5/IBD_fpet_cavalli_baja.pdf",
       device = "pdf", width = 4, height = 4, units = "in")

p.2.ns <- ggplot(dmat.2.ns,aes(x = V2, y = V1)) + 
  geom_point(size = 2) + 
  ylab("Nei's genetic distance (1983)") +
  xlab("Geographic Distance (km)") +
  ggtitle("F. petiolaris: Sonora and Central Mexico")  +
  stat_smooth(method = lm, formula = y~x, se = F) + 
  stat_regline_equation( label.x=400, label.y=0.008, size = 4)+
  annotate( 'text',y = 0.01, x = 600, label = tmp.2.ns, size=4)+
  theme_cowplot() +
  theme(aspect.ratio = 1,
        plot.title = element_text(size = 11),
        axis.text = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12))
plot(p.2.ns)
ggsave(p.2.ns, filename = "output/final_figures/Figure_5/IBD_fpet_Nei_83_no_south.pdf", device = "pdf", width = 4, height = 4, units = "in")

p.2 <- ggplot(dmat.2,aes(x = V2, y = V1)) + 
  geom_point(size = 2) + 
  ylab("Nei's genetic distance (1983)") +
  xlab("Geographic Distance (km)") +
  ggtitle("F. petiolaris: Baja California")  +
  stat_smooth(method = lm, formula = y~x, se = F) + 
  stat_regline_equation( label.x=400, label.y=0.0055, size = 4)+
  annotate( 'text',y = 0.005, x = 600, label = tmp.2, size=4)+
  theme_cowplot() +
  theme(aspect.ratio = 1,
        plot.title = element_text(size = 11),
        axis.text = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12))
plot(p.2)
ggsave(p.2, filename = "output/final_figures/Figure_5/IBD_fpet_Nei_83_baja.pdf",
       device = "pdf", width = 4, height = 4, units = "in")
