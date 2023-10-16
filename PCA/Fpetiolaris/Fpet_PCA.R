#Author: Kevin Quinteros
#Date: February 2, 2021
#Purpose: Calculate principal components for fpetiolaris data
#additional notes" We replace missing values with average of the site 
#sites with only 1 individual are removed
#color individuals by fine_regions instead of site

#Principal components analysis on Fpetiolaris ddRAD data
library(adegenet)
library(ggplot2)
library(ggrepel)
library(gridExtra)
library(tidyverse)
library(magrittr)
library(gridExtra)
library(cowplot)


###Set working directory--------------
setwd("~/Projects/git_repositories/github/Fig_wasp_phylogeography/")

#########################Loading recurring functions############################
#purpose: combine separate site data frames into one data frame
#and remove any columns which still have missing data
pca <- function(seq_tab){
  MG.sep <- do.call("rbind", seq_tab) #bind rows into data frame
  #subset matrix to exclude missing data
  #for example, if one individual was not sampled at a particular locus
  MG.sep.nomiss <- MG.sep[, !colSums(is.na(MG.sep)), drop = FALSE]
  return(MG.sep.nomiss)
}
###################### import and filter data (always run)######################
#import structure data (unlinked)
#import structure data (linked snps)
#I'm going with the unlinked data
#import the subset data set. Baja california had too many samples, so we 
#reduced individuals, but still kept representative from each sampling site.
data <- read.structure("Data/STRUCTURE_FILES/Tree/Fpet_tree_subsetted.str",
                                      n.ind = 60,
                                      n.loc = 1192,
                                      onerowperind = F,
                                      NA.char = "-9",
                                      row.marknames = 0,
                                      col.lab = 1,
                                      col.pop = 0)

#import structure data (unlinked)
#data <- read.structure("Data/Fpet_tree/Second_Finn_dataset/output/sub_fpet-s2_trim25-s3_clu_th90-s7_min125_outfiles/sub_fpet-s2_trim25-s3_clu_th90-s7_min125.stru",
#               n.ind = 250,
#               n.loc = 1192,
#              onerowperind = F,
#               NA.char = "-9",
#              row.marknames = 0,
#               col.lab = 1,
#               col.pop = 0)
#read in the locality data
loct <- read.csv("Data/Fpet_tree/Second_Finn_dataset/fpet_tree_locality.csv")

loct$Sample_ID <- gsub('-Tetas', 'T', loct$Sample_ID)
loct$Site <- gsub('-Tetas', 'T', loct$Site)
loct <- loct[loct$Sample_ID %in% rownames(data$tab),]

pop(data)  <- loct$Site
strata(data) <- loct
loct <- loct[!loct$Phylogroups == "Out",]
loct$Phylogroups[loct$Phylogroups == "SON_LAND"] <- "Inland Sonora"
loct$Phylogroups[loct$Phylogroups == "SON_COS"] <- "Coastal Sonora"
loct$Phylogroups[loct$Phylogroups == "BAJA"] <- "Baja California"
loct$Phylogroups[loct$Phylogroups == "CENTRAL"] <- "Central Mexico"
loct$Phylogroups[loct$Phylogroups == "SOUTH"] <- "Oaxaca"
#############################All Sites##########################################
########################### split sites (always run) %%%%%%%%%%%%%%%%%%%%%%%%%%%
#phased data set
seq_pop <- seppop(data) #separate genotypes per population
seq_pop$Out <- NULL
seq_tab <- lapply(seq_pop , tab, freq=T, NA.method= c("mean"))
seq_tab <- lapply(seq_tab ,
                    function(x) 
                      if (!is.matrix(x)) {x <- NULL} 
                    else {x}
)
seq_tab <- Filter(Negate(is.null), seq_tab)
####################calculate principal components %%%%%%%%%%%%%%%%%%
#combine individuals of interest
#run custom function on data set 
MG.sep  <- pca(seq_tab)
#calcalate principal components 
pc  <- dudi.pca(MG.sep , scale=F, scannf=F, nf=3)
#phased(calculate axis contribution to variance)
pc.ax1 <- pc$eig[1]/sum(pc$eig) %>% as.numeric()
pc.ax2 <- pc$eig[2]/sum(pc$eig) %>% as.numeric()
pc.ax3 <- pc$eig[3]/sum(pc$eig) %>% as.numeric()

######################Pop factors: All%%%%%%%%%%%%%%%%%%%%%%%%%%%
MG.sep <- as.data.frame(MG.sep)
MG.sep <- cbind(pc$li$Axis1, pc$li$Axis2,pc$li$Axis3, MG.sep)
names(MG.sep)[names(MG.sep) == 'pc$li$Axis1'] <- 'Axis1'
names(MG.sep)[names(MG.sep) == 'pc$li$Axis2'] <- 'Axis2'
names(MG.sep)[names(MG.sep) == 'pc$li$Axis3'] <- 'Axis3'
loct$Phylogroups <- factor(loct$Phylogroups)
#############plot PCA: no moax%%%%%%%%%%%%%%%%%
library(RColorBrewer)

#plot Axis 1 and Axis 2
fpet.all <- ggplot(MG.sep) +
  geom_point(size=3, shape = 21, color = "black", aes(x=Axis1,
                                                      y=Axis2,
                                                      fill=loct$Phylogroups)) +
  xlab(sprintf("PC 1 - %s%%", format(round(pc.ax1 * 100, 2), nsmall=2))) +
  ylab(sprintf("PC 2 - %s%%", format(round(pc.ax2 * 100, 2), nsmall=2))) +
  labs(fill="Regional Sites")+
  coord_fixed() +
  theme_cowplot() +
  theme(aspect.ratio = 1, 
        legend.text=element_text(size=9),
        legend.title=element_text(size=9),
        legend.position = c(0.15,0.75),
        axis.text = element_text(size = 9),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10)) +
  scale_fill_manual(values = c("#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7"),
                    breaks= c("Baja California","Coastal Sonora","Inland Sonora","Central Mexico","Oaxaca")) 
plot(fpet.all)

ggsave(fpet.all, filename = "output/final_figures/Figure_2/fpet_all_pca12.pdf", device = "pdf", width = 4, height = 4, units = "in", dpi=300,)

#plot Axis 1 and Axis 3

fpet.all <- ggplot(MG.sep) +
  geom_point(size=3, shape = 21, color = "black", aes(x=Axis1,
                                                      y=Axis3,
                                                      fill=loct$Phylogroups)) +
  xlab(sprintf("PC 1 - %s%%", format(round(pc.ax1 * 100, 2), nsmall=2))) +
  ylab(sprintf("PC 3 - %s%%", format(round(pc.ax3 * 100, 2), nsmall=2))) +
  labs(fill="Regional Sites")+
  coord_fixed() +
  theme_cowplot() +
  theme(aspect.ratio = 1, 
        legend.text=element_text(size=9),
        legend.title=element_text(size=9),
        legend.position = c(0.05,0.2),
        axis.text = element_text(size = 9),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10)) +
  scale_fill_manual(values = c("#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7"),
                    breaks= c("Baja California","Coastal Sonora","Inland Sonora","Central Mexico","Oaxaca")) 
plot(fpet.all)

ggsave(fpet.all, filename = "output/final_figures/Figure_2/fpet_all_pca13.pdf", device = "pdf", width = 4, height = 4, units = "in", dpi=300)

#plot Axis 2 and Axis 3
fpet.all <- ggplot(MG.sep) +
  geom_point(size=3, shape = 21, color = "black", aes(x=Axis2,
                                                      y=Axis3,
                                                      fill=loct$Phylogroups)) +
  xlab(sprintf("PC 2 - %s%%", format(round(pc.ax2 * 100, 2), nsmall=2))) +
  ylab(sprintf("PC 3 - %s%%", format(round(pc.ax3 * 100, 2), nsmall=2))) +
  labs(fill="Regional Sites")+
  coord_fixed() +
  theme_cowplot() +
  theme(aspect.ratio = 1, 
        legend.text=element_text(size=9),
        legend.title=element_text(size=9),
        legend.position = c(0.65,0.20),
        axis.text = element_text(size = 9),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10)) +
  scale_fill_manual(values = c("#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7"),
                    breaks= c("Baja California","Coastal Sonora","Inland Sonora","Central Mexico","Oaxaca")) 
plot(fpet.all)

#save to file
ggsave(fpet.all, filename = "output/final_figures/Figure_2/fpet_all_pca23.pdf", device = "pdf", width = 4, height = 4, units = "in", dpi=300)
