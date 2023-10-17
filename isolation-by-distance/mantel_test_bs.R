#Author: Kevin Quinteros
#Date: Feb 15, 2022
#This code preforms a mantel test on genetic distances 
#this will be for the pegoscapus
#Author: Kevin Quinteros
#I'm supposed to just compare North vs South populations 

#load libraries
library(tidyverse)
library(ade4)
library(adegenet)
library(ggpubr)

########################Load Dataset############################################
load("Data/Gen_dist/Nei_83_genetic_dist.RDA") #genetic distance 2
load("Data/Gen_dist/geo_distance.RDA") #geographical distance 
##################Load site information#########################
sites <- read_csv("Data/Locality_Fpet.csv")
sites <- sites[-3,] #remove FW268 since it did not phase
sites <- sites[-1]
sites <- sites[-6]
sites %<>% distinct(Site, .keep_all = TRUE)
########################preform mantel test#####################################
result.2.bs = mantel.randtest(as.dist(Nei.83$Dgen.bs), as.dist(geo_distance$Dgeo.bs), nrepet = 999)

view(Nei.83$Dgen.bs)
view(geo_distance$Dgeo.bs)

tmp.bs <- paste('R =',
                 round(result.2.bs$obs, digits=3),
                 '|', 'p =',
                 round(result.2.bs$pvalue, 
                       digits=3), sep=' ')
save(result.2.bs ,file = "output/Mantel_test/pego_mantel_test_Baja_Sono_sina.RDA")

dmat.2.bs <- cbind(as.vector(Nei.83$Dgen.bs), as.vector(geo_distance$Dgeo.bs)) %>% as.data.frame()
dmat.2.bs <- dmat.2.bs[rowSums(dmat.2.bs[])>9.504539e-05,]

p.2.bs <- ggplot(dmat.2.bs,aes(x = V2, y = V1)) + 
  geom_point(size = 2) + 
  ylab("Nei's genetic distance (1983)") +
  xlab("Geographic Distance (km)") +
  ggtitle("Pegoscapus: Baja, Sonora and Sinaloa")  +
  stat_smooth(method = lm, formula = y~x, se = F) + 
  stat_regline_equation( label.x=60, label.y=0.005, size = 4)+
  annotate( 'text',y = 0.004, x = 300, label = tmp.bs, size=4)+
  theme_cowplot() +
  theme(aspect.ratio = 1,
        plot.title = element_text(size = 11),
        axis.text = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12))
plot(p.2.bs)
ggsave(p.2.bs, filename = "output/final_figures/Figure_5/IBD_Nei_83_bs.pdf",
       device = "pdf", width = 4, height = 4, units = "in")


#########################
#Nei's distances 
rownames(Nei.83$Dgen.bs) <- Baja_Sin$Site
colnames(Nei.83$Dgen.bs) <- Baja_Sin$Site
#geographical distances
rownames(geo_distance$Dgeo.bs) <- Baja_Sin$Site
colnames(geo_distance$Dgeo.bs) <- Baja_Sin$Site

Baja_Sin$Region
baja <- Baja_Sin[Baja_Sin$State == "BAJA",] %>% 
  select(Site) 
no_baja <- Baja_Sin[Baja_Sin$State != "BAJA",] %>%
  select(Site) 

baja<- baja$Site
no_baja <- no_baja$Site

#compare only the baja  and  sono & sina regions
#functions
`%nin%` = Negate(`%in%`)
Nei.83$Dgen.bs <- Nei.83$Dgen.bs[,colnames(Nei.83$Dgen.bs) %nin% baja]
Nei.83$Dgen.bs <- Nei.83$Dgen.bs[rownames(Nei.83$Dgen.bs) %nin% no_baja,]

geo_distance$Dgeo.bs <- geo_distance$Dgeo.bs[,colnames(geo_distance$Dgeo.bs) %nin% baja]
geo_distance$Dgeo.bs <- geo_distance$Dgeo.bs[rownames(geo_distance$Dgeo.bs) %nin%  no_baja,]


###distance mat
geo_dist <- geo_distance$Dgeo.bs
nei_dist <- Nei.83$Dgen.bs


result.bs <- mantel.randtest(as.dist(nei_dist),as.dist(geo_dist), nrepet = 999)

tmp.bs <- paste('R =',
                 round(result.bs$obs, digits=3),
                 '|', 'p =',
                 round(result.bs$pvalue, 
                       digits=3), sep=' ')
dmat.nei <- cbind(as.vector(nei_dist), as.vector(geo_dist)) %>% as.data.frame()
dmat.nei <- dmat.nei[rowSums(dmat.nei[])>1,]

p.nei <- ggplot(dmat.nei,aes(x = V2, y = V1)) + 
  geom_point(size = 2) + 
  ylab("Nei's genetic distance (1983)") +
  xlab("Geographic Distance (km)") +
  ggtitle("Pegoscapus: baja 
          vs sono-sina")  +
  theme_cowplot() +
  theme(aspect.ratio = 1,
        plot.title = element_text(size = 11),
        axis.text = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12))

p.nei.lm <- p.nei + geom_smooth(method = lm, formula = y~x, se = F) +
  stat_regline_equation( label.x=50, label.y=0.01, size = 4)+
  annotate( 'text',y = 0.009, x = 700, label = tmp.bs, size=4)

plot(p.nei.lm)
ggsave(p.nei.lm, filename = "output/final_figures/Figure_5/IBD_nei_compare_baja_sino_sono_lm.pdf",
       device = "pdf", width = 4, height = 4, units = "in")
