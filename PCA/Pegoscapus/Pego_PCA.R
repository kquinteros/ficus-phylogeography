#Author: Kevin Quinteros
#Date: February 2, 2021
#Purpose: Calculate principal components for pollinator data
#additional notes" We replace missing values with average of the site 
#sites with only 1 individual are removed
#color individuals by fine_regions instead of site
#final PCA for figure 2 using the phased dataset

#Principal components analysis on Pollinator UCE SNP data
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
#Phased SNPs
load("Data/Phased_data/UCE_50p/Method1/Phased_snps/phased_snp_50.Rda")
#############################All Sites##########################################
########################### split sites (always run) %%%%%%%%%%%%%%%%%%%%%%%%%%%
#phased data set
seq_pop.p <- seppop(phased_snp_50) #separate genotypes per population
seq_tab.p <- lapply(seq_pop.p, tab, freq=T, NA.method= c("mean"))
seq_tab.p <- lapply(seq_tab.p,
                    function(x) 
                      if (!is.matrix(x)) {x <- NULL} 
                    else {x}
)
seq_tab.p <- Filter(Negate(is.null), seq_tab.p )

####################calculate principal components %%%%%%%%%%%%%%%%%%
#combine individuals of interest
#run custom function on data set 
MG.sep.p <- pca(seq_tab.p)

#calcalate principal components 
pca.p <- dudi.pca(MG.sep.p, scale=F, scannf=F, nf=3)

#phased(calculate axis contribution to variance)
pca.p.ax1 <- pca.p$eig[1]/sum(pca.p$eig)
pca.p.ax2 <- pca.p$eig[2]/sum(pca.p$eig)
pca.p.ax3 <- pca.p$eig[3]/sum(pca.p$eig)

######################Pop factors: All%%%%%%%%%%%%%%%%%%%%%%%%%%%
MG.sep.p <- as.data.frame(MG.sep.p)
MG.sep.p <- cbind(pca.p$li$Axis1, pca.p$li$Axis2, MG.sep.p)
names(MG.sep.p)[names(MG.sep.p) == 'pca.p$li$Axis1'] <- 'Axis1'
names(MG.sep.p)[names(MG.sep.p) == 'pca.p$li$Axis2'] <- 'Axis2'
loct <- read_csv("Data/Locality_Fpet.csv")
#change name of fine regions 
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

reference <- c("Baja California", "Coastal Sonora","Inland Sonora", 
               "Coastal Sinaloa","Inland Sinaloa", "Zacatecas",
               "Jalisco","Morelos", "Northern Oaxaca", "Oaxaca") #site reference
loct <- loct[order(factor(loct$fine_region, levels = reference)),] #align reference to
loct <- loct %>% filter(Sample %in% rownames(MG.sep.p))
MG.sep.p <-  MG.sep.p[match(loct$Sample,rownames(MG.sep.p)),]
loct$fine_region <- factor(loct$fine_region)
#set my color palette
safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", 
                             "#117733", "#332288", "#AA4499",
                             "#44AA99", "#999933", "#882255", 
                             "#661100", "#6699CC", "#888888")
####################plot PCA%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
peg.all <- ggplot(MG.sep.p[3:1828], aes(x = MG.sep.p$Axis1,
                                        y = MG.sep.p$Axis2,
                                        fill=loct$fine_region)) +
  geom_point(size=3, shape = 21, color="black") +
  scale_fill_manual(values = safe_colorblind_palette,breaks = reference ) +
  xlab(sprintf("PC 1 - %s%%", format(round(pca.p.ax1 * 100, 2), nsmall=2))) +
  ylab(sprintf("PC 2 - %s%%", format(round(pca.p.ax2 * 100, 2), nsmall=2))) +
  labs(fill= "Phylogroups") +
  coord_fixed () +
  theme_cowplot() +
  theme(aspect.ratio = 1,
        legend.text=element_text(size=9),
        legend.title=element_text(size=9),
        legend.position = c(1,0.5),
        axis.text = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10))
plot(peg.all)
ggsave(peg.all, filename = "output/final_figures/Figure_3/peg_all_pca.pdf",
       device = "pdf", width = 6.5, height = 4, units = "in", dpi=300)

#############################No Morelos and Oaxaca##############################
seq_pop.p <- phased_snp_50
#removing 218 hybrid individual 
seq_pop.p <- seq_pop.p[!row.names(seq_pop.p$tab) %in% c("P342_petiolaris")]
seq_pop.p <- seppop(seq_pop.p) #separate genotypes per population
seq_tab.p <- lapply(seq_pop.p, tab, freq=T, NA.method= c("mean"))
seq_tab.p <- lapply(seq_tab.p,
                    function(x) 
                      if (!is.matrix(x)) {x <- NULL} 
                    else {x}
)
seq_tab.p <- Filter(Negate(is.null), seq_tab.p )
##format seq_pop objects#
sites <- read_csv("Data/Locality_Fpet.csv")
moax <- sites %>%
  filter(State %in% c("OAX", "MOR"))
moax <- unique(moax$Site)
#remove all sites form Morelos and Oaxaca 
toRemove <- which(names(seq_tab.p) %in% moax)
toRemove <- rev(toRemove)
for (i in toRemove){
  seq_tab.p[[i]] <- NULL
}

####################calculate principal components 
#combine individuals of interest
#run custom function on data set 
MG.sep.p <- pca(seq_tab.p)
#calcalate principal components 
pca.p <- dudi.pca(MG.sep.p, scale=F, scannf=F, nf=3)

#phased(calculate axis contribution to variance)
pca.p.ax1 <- pca.p$eig[1]/sum(pca.p$eig)
pca.p.ax2 <- pca.p$eig[2]/sum(pca.p$eig)
pca.p.ax3 <- pca.p$eig[3]/sum(pca.p$eig)

######################Pop factors: All%%%%%%%%%%%%%%%%%%%%%%%%%%%
MG.sep.p <- as.data.frame(MG.sep.p)
MG.sep.p <- cbind(pca.p$li$Axis1, pca.p$li$Axis2, MG.sep.p)
names(MG.sep.p)[names(MG.sep.p) == 'pca.p$li$Axis1'] <- 'Axis1'
names(MG.sep.p)[names(MG.sep.p) == 'pca.p$li$Axis2'] <- 'Axis2'

loct <- read_csv("Data/Locality_Fpet.csv")
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

reference <- c("Baja California", "Coastal Sonora","Inland Sonora",
               "Coastal Sinaloa","Inland Sinaloa", "Zacatecas",
               "Jalisco","Morelos", "Northern Oaxaca", "Oaxaca") #site reference
loct <- loct[order(factor(loct$fine_region, levels = reference)),] #align reference to
loct <- loct %>% filter(Sample %in% rownames(MG.sep.p))
MG.sep.p <-  MG.sep.p[match(loct$Sample,rownames(MG.sep.p)),]
loct$fine_region <- factor(loct$fine_region)
#set color paletter
safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77",
                             "#117733", "#332288", "#AA4499","#44AA99",
                             "#999933", "#882255", "#661100", "#6699CC", "#888888")

#############plot PCA: no moax%%%%%%%%%%%%%%%%%
peg.no_moax <- ggplot(MG.sep.p[3:1998], aes(x = MG.sep.p$Axis1, y = MG.sep.p$Axis2, color=loct$fine_region)) +
  geom_point(size=3) +
  #scale_color_manual(values = wes_palette("Darjeeling1",7,type = c( "continuous")),breaks = reference ) +
  scale_color_manual(values = safe_colorblind_palette[1:7],breaks = reference ) +
  xlab(sprintf("PC 1 - %s%%", format(round(pca.p.ax1 * 100, 2), nsmall=2))) +
  ylab(sprintf("PC 2 - %s%%", format(round(pca.p.ax2 * 100, 2), nsmall=2))) +
  coord_fixed() +
  theme_cowplot() +
  theme(aspect.ratio = 1,
        legend.position = "none",
        axis.text = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10))
plot(peg.no_moax)

ggsave(peg.no_moax,
       filename = "output/final_figures/Figure_3/peg_no_moax_pca.pdf", 
       device = "pdf", width = 4, height = 4, units = "in", dpi=300)
#############################Morelos and Oaxaca#################################
seq_pop.p <- seppop(phased_snp_50) #separate genotypes per population
seq_tab.p <- lapply(seq_pop.p, tab, freq=T, NA.method= c("mean"))
seq_tab.p <- lapply(seq_tab.p,
                    function(x) 
                      if (!is.matrix(x)) {x <- NULL} 
                    else {x}
)

seq_tab.p <- Filter(Negate(is.null), seq_tab.p )
##format seq_pop objects#
sites <- read_csv("Data/Locality_Fpet.csv")
moax <- sites %>%
  filter(State %in% c("OAX", "MOR"))
moax <- unique(moax$Site)
`%notin%` = Negate(`%in%`)

#remove all sites form Morelos and Oaxaca 
toRemove <- which(names(seq_tab.p) %notin% moax)
toRemove <- rev(toRemove)
for (i in toRemove){
  seq_tab.p[[i]] <- NULL
}

####################calculate principal components
#combine individuals of interest
#run custom function on data set 
MG.sep.p <- pca(seq_tab.p)

#calcalate principal components 
pca.p <- dudi.pca(MG.sep.p, scale=F, scannf=F, nf=3)

#phased(calculate axis contribution to variance)
pca.p.ax1 <- pca.p$eig[1]/sum(pca.p$eig)
pca.p.ax2 <- pca.p$eig[2]/sum(pca.p$eig)
pca.p.ax3 <- pca.p$eig[3]/sum(pca.p$eig)

######################Pop factors: All%%%%%%%%%%%%%%%%%%%%%%%%%%%
MG.sep.p <- as.data.frame(MG.sep.p)
MG.sep.p <- cbind(pca.p$li$Axis1, pca.p$li$Axis2, MG.sep.p)
names(MG.sep.p)[names(MG.sep.p) == 'pca.p$li$Axis1'] <- 'Axis1'
names(MG.sep.p)[names(MG.sep.p) == 'pca.p$li$Axis2'] <- 'Axis2'

##########################Get population Factors
loct <- read_csv("Data/Locality_Fpet.csv")
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

reference <- c("Baja California", "Coastal Sonora","Inland Sonora", 
               "Coastal Sinaloa","Inland Sinaloa", "Zacatecas", "Jalisco",
               "Morelos", "Northern Oaxaca", "Oaxaca") #site reference
loct <- loct[order(factor(loct$fine_region, levels = reference)),] #align reference to
loct <- loct %>% filter(Sample %in% rownames(MG.sep.p))
MG.sep.p <-  MG.sep.p[match(loct$Sample,rownames(MG.sep.p)),]
loct$fine_region <- factor(loct$fine_region)
#set color paletter
safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733",
                             "#332288", "#AA4499","#44AA99", "#999933",
                             "#882255", "#661100", "#6699CC", "#888888")


#############plot PCA: no moax%%%%%%%%%%%%%%%%%
peg.moax <- ggplot(MG.sep.p[3:2498],
                   aes(x = MG.sep.p$Axis1, y = MG.sep.p$Axis2, fill=loct$fine_region)) +
  geom_point(size=3, shape = 21, color = "black") +
  scale_fill_manual(values = c( "#999933","#882255", "#661100")) +
  xlab(sprintf("PC 1 - %s%%", format(round(pca.p.ax1 * 100, 2), nsmall=2))) +
  ylab(sprintf("PC 2 - %s%%", format(round(pca.p.ax2 * 100, 2), nsmall=2))) +
  coord_fixed() +
  theme_cowplot() +
  theme(aspect.ratio = 1,
        legend.position = "none",
        axis.text = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10))
plot(peg.moax)

ggsave(peg.moax, filename = "output/final_figures/Figure_3/peg_moax_pca.pdf",
       device = "pdf", width = 4, height = 4, units = "in", dpi=300)
