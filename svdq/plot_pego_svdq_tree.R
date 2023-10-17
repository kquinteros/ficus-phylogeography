#Author: Kevin Quinteros
#Date: Feb 7th, 2022 
#Purpose: Plot phylogenies of F.petiolaris
#two trees will need to plotted. The lineage tree and the sub (OTUs)

#load library
library(ape)
library(tidyverse)
library(ggtree)

#setwd
setwd("~/Projects/git_repositories/github/Fig_wasp_phylogeography/")

##################Plot lineae tree ##################################
#read in lineage tree
line <- ape::read.tree("output/svdquartets/wasp/pegoscapus_comp_50p_lineage_rooted.tre")
#trim the outgroup species.
line <- drop.tip(line, c("P274_petiolaris", "P341_petiolaris"))
#set up 
loct <- read.csv("Data/Locality_Fpet.csv")
loct$fine_region[loct$fine_region == "BAJA"] <- "Baja California"
loct$fine_region[loct$fine_region == "SON_COS"] <- "Coastal Sonora"
loct$fine_region[loct$fine_region == "SON_IN_LAND"] <- "Inland Sonora"
loct$fine_region[loct$fine_region == "SINA_COS"] <- "Coastal Sinaloa"
loct$fine_region[loct$fine_region == "SINA_IN_LAND"] <- "Inland Sinaloa"
loct$fine_region[loct$fine_region == "ZAC"] <- "Zacatecas"
loct$fine_region[loct$fine_region == "JAL"] <- "Jalisco"
loct$fine_region[loct$fine_region == "MOR"] <- "Morelos"
loct$fine_region[loct$fine_region == "OAX_NOR"] <- "Northern Oaxaca"
loct$fine_region[loct$fine_region == "OAX"] <- "Oaxaca"

loct <- loct[loct$Sample %in% line$tip.label,]
loct <- loct[match(line$tip.label, loct$Sample),] 
loct <- loct[,c(1,12)]

#set order 
reference <- c("Baja California", "Coastal Sonora","Inland Sonora", 
               "Coastal Sinaloa","Inland Sinaloa", "Zacatecas",
               "Jalisco","Morelos", "Northern Oaxaca", "Oaxaca") #site reference
#color palette 
safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", 
                             "#117733", "#332288", "#AA4499","#44AA99",
                             "#999933", "#882255", "#661100", "#6699CC", "#888888")
#preliminary tree
trel <- ggtree(line, branch.length="none")  
fin <- trel %<+% loct + 
  geom_tippoint(size=2.5, shape = 21,
                aes(fill = factor(fine_region)),
                color = "black") +
  scale_fill_manual(values = safe_colorblind_palette,
                    breaks= reference) +
  xlim(0,22)+
  theme(legend.position = c(0.7,0.2), 
        legend.title = element_blank(), # no title
        legend.key = element_blank()) # no keys
plot(fin)
ggsave(fin, filename = "output/final_figures/Figure_4/svdq_wasp_lineage_outgroup_removed.pdf",
       device = "pdf",
       width = 4,
       height = 10,
       units = "in",
       dpi = 300)

##############################Trim additional Morelos samples##################
#these samples had high amounts of missing data. 
line <- drop.tip(line, tip = c("P345_petiolaris", "P346_petiolaris", "P347_petiolaris"))
#
loct <- loct[loct$Sample %in% line$tip.label,]

trel.1 <- ggtree(line, branch.length="none")  
fin.1 <- trel.1 %<+% loct + 
  geom_tippoint(size=2.5, shape = 21,
                aes(fill = factor(fine_region)),
                color = "black") +
  scale_fill_manual(values = safe_colorblind_palette,
                    breaks= reference) +
  xlim(0,22)+
  theme(legend.position = c(0.7,0.2), 
        legend.title = element_blank(), # no title
        legend.key = element_blank()) # no keys
plot(fin.1)
ggsave(fin.1, filename = "output/final_figures/Figure_4/svdq_wasp_lineage_outgroup_removed_morelos_sample_removed.pdf",
       device = "pdf",
       width = 4,
       height = 10,
       units = "in",
       dpi = 300)
