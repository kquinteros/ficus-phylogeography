#Author: Kevin Quinteros
#Date: Feb 6th, 2022 
#Purpose: Plot phylogenies of F.petiolaris
#two trees will need to plotted. The lineage tree and the sub (OTUs)

#load library
library(ape)
library(tidyverse)
library(ggtree)

#setwd
setwd("~/Projects/git_repositories/github/Fig_wasp_phylogeography/")
##################################import tree##################################
sub <- read.tree("output/svdquartets/fpet/ficus_trim25-s3_clu_th90-s7_min125_sub.tre")
##############################rooted trees ###################################
sub <- root(sub,outgroup = "South", resolve.root = T)
sub$tip.label[sub$tip.label == "Baja"] <- "Baja California"
sub$tip.label[sub$tip.label == "South"] <- " Oaxaca"
sub$tip.label[sub$tip.label == "Sonora_coast"] <- "Coast Sonora"
sub$tip.label[sub$tip.label == "Sonora_land"] <- "Inland Sonora"
sub$tip.label[sub$tip.label == "Central"] <- "Central Mexico"
sub <- ape::rotate(sub, node = 7)
sub <- ape::rotate(sub, node = 8)
sub <- ape::rotate(sub, node = 6)
sub <- ape::rotate(sub, node = 9)
plot(sub)
#######################plotting of tree #######################################
tres <- ggtree(sub,branch.length="none") + geom_tiplab(offset = 0.1) 
#used to show the node number of tree
#ggtree(sub, branch.length="none") + geom_text(aes(label=node), hjust=-.3)
tres <- tres + theme_tree2() +  xlim(0, 4.5) + theme_tree() 
plot(tres)
ggsave(tres, filename = "output/final_figures/Figure_2/fpet_sub.pdf",
       width = 4, height = 4, units = "in", dpi =300)

##################Plot lineae tree ##################################
#read in lineage tree
line <- read.tree("output/svdquartets/fpet/ficus_trim25-s3_clu_th90-s7_min125_sub_lineage_rooted.tre")
site <- read.csv("Data/Fpet_tree/Second_Finn_dataset/fpet_tree_locality.csv")
site$Phylogroups[site$Phylogroups == "SOUTH"] <- "Oaxaca"
site$Phylogroups[site$Phylogroups == "BAJA"] <- "Baja California"
site$Phylogroups[site$Phylogroups == "CENTRAL"] <- "Central"
site$Phylogroups[site$Phylogroups == "SON_COS"] <- "Coastal Sonora"
site$Phylogroups[site$Phylogroups == "SON_LAND"] <- "Inland Sonora"
site <- site %>% 
  mutate(Sample_ID = str_replace(Sample_ID, "-Tetas-", "T_"))
site <- site %>% 
  mutate(Sample_ID = str_replace(Sample_ID, "-", "_"))
site <- site[site$Sample_ID %in% line$tip.label,]
site <- site[match(line$tip.label, site$Sample_ID),] 
site <- site[,c(1,9)]

#preliminary tree
trel <- ggtree(line, branch.length="none")  
fin <- trel %<+% site + 
  geom_tippoint(size=3, shape = 21,
                aes(fill = factor(Phylogroups)),
              color = "black") +
  scale_fill_manual(values = c("#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7"),
                    breaks= c("Baja California","Coastal Sonora","Inland Sonora","Central","Oaxaca")) +
  xlim(0,22)+
  theme(legend.position = c(0.7,0.2), 
        legend.title = element_blank(), # no title
        legend.key = element_blank()) # no keys
plot(fin)
ggsave(fin, filename = "output/final_figures/Figure_2/svdq_fpet_lineage.pdf",
       device = "pdf",
       width = 3,
       height = 10,
       units = "in",
       dpi = 300)

