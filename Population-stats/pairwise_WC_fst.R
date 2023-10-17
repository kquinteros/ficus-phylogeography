#Date: March 11, 2020
#Author: Kevin Quinteros 
#Updated: June 10, 2020
#Purpose: Calculated pairwise Fst between sampling Regions and inferred genetic
#clusters
#Calculateing Fst between cluster for fpet was done in separate 
#file. See Cluster Pairwise Fst.R (fpet)
##########################set working environment###############################
#Load required libraries
library(adegenet)
library(tidyverse)
library(hierfstat)
library(magrittr)

#Set working directory
setwd("~/Projects/git_repositories/github/Fig_wasp_phylogeography/")


###############################read in data#####################################
pego <- read.csv("output/Genetic_Stats/data/Pego_filtered-n≥2-≤50%-NA_cluster data.csv")

pego_pairwise_fst <- pairwise.WCfst(pego[-2]) #calculate WC_FST

#write to file
write.csv(pego_pairwise_fst, file = "output/Genetic_Stats/Fst/pego_pairwise_WC_fst_region.csv")

pego_pairwise_fst_bp <- boot.ppfst(pego[-2], nboot = 999) #bootstrap
#save to file
save(pego_pairwise_fst_bp, file = "output/Genetic_Stats/Fst/pego_pairwise_WC_fst_bp_region.csv")

ll <- pego_pairwise_fst_bp$ll #lower confidence limit
ul <- pego_pairwise_fst_bp$ul #upper confidence limit

write.csv(ll, "output/Genetic_Stats/Fst/pego_pairwise_WC_fst_ll_region.csv")
write.csv(ul, "output/Genetic_Stats/Fst/pego_pairwise_WC_fst_ul_region.csv")

########################format the complete table###############################

#round the table
pego_pairwise_fst %<>% round(digits = 4)
ll %<>% round(digits = 4)
ul %<>% round(digits = 4)

#format the table
tmp <- matrix( paste0("(",ll, " - ", ul, ")"), 
        nrow=nrow(ll), dimnames=dimnames(ll))
pego_Fst_table <- tmp <- matrix( paste0(pego_pairwise_fst, "\n", tmp), 
                      nrow=nrow(tmp), dimnames=dimnames(tmp))
remove(tmp)

#give rownames and column names
rownames(pego_Fst_table) <- c("Baja California", "Coastal Sonora","Inland Sonora",
  "Coastal Sinaloa", "Inland Sinaloa","Zacatecas",
  "Jalisco", "Morelos", "Northern  Oaxaca", "Southern Oaxaca") 
colnames(pego_Fst_table) <- c("Baja California", "Coastal Sonora","Inland Sonora",
                              "Coastal Sinaloa", "Inland Sinaloa","Zacatecas",
                              "Jalisco", "Morelos", "Northern  Oaxaca", "Southern Oaxaca") 
#change the diagonals and the lower triangle of matrix
diag(pego_Fst_table) <- NA
pego_Fst_table[lower.tri(pego_Fst_table)] <- NA

#write to file
write.table(pego_Fst_table,
            file = "output/Genetic_Stats/Fst/pego_Fst-table_region.txt",
            sep = "\t")

############################Fst for phylogroups#################################
#do the same but with inferred phylogroups
pego <- add_column(pego, cluster = NA, .before = "pop")
pego <- pego %>% mutate(cluster = case_when(
  pop==1 ~ 1, # Baja
  pop==1 ~ 1,
  pop==1 ~ 1,
  pop==1 ~ 1,
  pop==1 ~ 1,
  pop==1 ~ 1,
  pop==1 ~ 1,
  pop==2 ~ 1,  #Coastal Sonora
  pop==2 ~ 1,
  pop==2 ~ 1,
  pop==2 ~ 1,
  pop==3 ~ 1, #Inland Sonora
  pop==3 ~ 1, 
  pop==3 ~ 1,
  pop==3 ~ 1, 
  pop==4 ~ 1, #Coastal Sinaloa
  pop==4 ~ 1,
  pop==4 ~ 1, 
  pop==5 ~ 1, #Inland Sinaloa
  pop==5 ~ 1, 
  pop==6 ~ 1, #Zacatecas
  pop==7 ~ 1, #Jalisco
  pop==7 ~ 1,
  pop==8 ~ 2, #Morelos
  pop==8 ~ 2, 
  pop==9 ~ 3, #Northern Oaxaca
  pop==10 ~ 4, #Southern Oaxaca
  pop==10 ~ 4, 
  pop==10 ~ 4
))

pego_pairwise_fst.cluster <- pairwise.WCfst(pego[-c(2:3)]) #calculate WC_FST

#write to file
write.csv(pego_pairwise_fst.cluster, file = "output/Genetic_Stats/Fst/pego_pairwise_WC_fst_cluster.csv")
#bootstrap
pego_pairwise_fst_bp.cluster <- boot.ppfst(pego[-c(2:3)], nboot = 999)
#save to file
save(pego_pairwise_fst_bp.cluster, file = "output/Genetic_Stats/Fst/pego_pairwise_WC_fst_bp_cluster.csv")

ll.cluster <- pego_pairwise_fst_bp.cluster$ll #lower confidence limit
ul.cluster <- pego_pairwise_fst_bp.cluster$ul #upper confidence limit

write.csv(ll.cluster, "output/Genetic_Stats/Fst/pego_pairwise_WC_fst_ll_cluster.csv")
write.csv(ul.cluster, "output/Genetic_Stats/Fst/pego_pairwise_WC_fst_ul_cluster.csv")

######################format the complete cluster###############################

#round the table
pego_pairwise_fst.cluster %<>% round(digits = 4)
ll.cluster %<>% round(digits = 4)
ul.cluster %<>% round(digits = 4)

#format the table
tmp <- matrix( paste0("(",ll.cluster, " - ", ul.cluster, ")"), 
               nrow=nrow(ll.cluster), dimnames=dimnames(ll.cluster))
pego_Fst_table.cluster <- matrix( paste0(pego_pairwise_fst.cluster, "\n", tmp), 
                                 nrow=nrow(tmp), dimnames=dimnames(tmp))
remove(tmp)
view(pego_Fst_table.cluster)
#give rownames and column names
rownames(pego_Fst_table.cluster) <- c("North", "Morelos",
                              "Northern  Oaxaca", "Southern Oaxaca") 
colnames(pego_Fst_table.cluster) <- c("North", "Morelos",
                              "Northern  Oaxaca", "Southern Oaxaca")  
#change the diagonals and the lower triangle of matrix
diag(pego_Fst_table.cluster) <- NA
pego_Fst_table.cluster[lower.tri(pego_Fst_table.cluster)] <- NA

#write to file
write.table(pego_Fst_table.cluster,
            file = "output/Genetic_Stats/Fst/pego_Fst-table_cluster.txt",
            sep = "\t")

######################do the same for the fpetiolaris##########################
#read the data
fpet <- read.csv(file = "output/Genetic_Stats/data/Fpet filtered-n≥3-≤50%-NA_cluster data.csv")
#pairwise WC_Fst
fpet_pairwise_fst <- pairwise.WCfst(fpet[-2])
#write to file
write.csv(fpet_pairwise_fst , file = "output/Genetic_Stats/Fst/fpet_pairwise_WC_fst.csv")
#bootstrap
fpet_pairwise_fst_bp <- boot.ppfst(fpet[-2], nboot = 999)
#save to file
save(fpet_pairwise_fst_bp,
     file ="output/Genetic_Stats/Fst/fpet_pairwise_WC_fst_bp.csv")

ll.fpet <- fpet_pairwise_fst_bp$ll #lower confidence limit
ul.fpet <- fpet_pairwise_fst_bp$ul #upper confidence limit
#write to file
write.csv(ll.fpet, "output/Genetic_Stats/Fst/fpet_pairwise_WC_fst_ll.csv")
write.csv(ul.fpet, "output/Genetic_Stats/Fst/fpet_pairwise_WC_fst_ul.csv")

#####################format the complete table fpet#############################
#round the table
fpet_pairwise_fst %<>% round(digits = 4)
ll.fpet %<>% round(digits = 4)
ul.fpet %<>% round(digits = 4)

#format the table
tmp <- matrix( paste0("(",ll.fpet, " - ", ul.fpet, ")"), 
               nrow=nrow(ll.fpet), dimnames=dimnames(ll.fpet))
fpet_Fst_table <- matrix( paste0(fpet_pairwise_fst, "\n", tmp), 
                                  nrow=nrow(tmp), dimnames=dimnames(tmp))
remove(tmp)
view(fpet_Fst_table)
#give rownames and column names
rownames(fpet_Fst_table) <- c("Baja California","Coastal Sonora",
                                      "Inland Sonora","Central Mexico",
                                      "Southern Oaxaca") 
colnames(fpet_Fst_table) <- c("Baja California","Coastal Sonora",
                                      "Inland Sonora","Central Mexico",
                                      "Southern Oaxaca")  
#change the diagonals and the lower triangle of matrix
diag(fpet_Fst_table) <- NA
fpet_Fst_table[lower.tri(fpet_Fst_table)] <- NA

#write to file
write.table(fpet_Fst_table,
            file = "output/Genetic_Stats/Fst/fpet_Fst-table.txt",
            sep = "\t")

################################test line A##############################
above <- c("112","158","172")
below <- c("113","95","201","179","204")
A <- loct[loct$Site %in% above,]
fpet.A <- fpet[rownames(fpet) %in% A$Sample_ID,]
fpet.A$pop <- fpet.A$pop[fpet.A$pop=="BAJA"] <- "Above"
droplevels()

B <- loct[loct$Site %in% below,]
fpet.B <- fpet[rownames(fpet) %in% B$Sample_ID,]
fpet.B$pop <- fpet.B$pop[fpet.B$pop=="BAJA"] <- "Below"

fpet.AB <- rbind(fpet.A,fpet.B)
remove(above, below, fpet.A, fpet.B,A,B)
fpet_pairwise_fst_A <- pairwise.WCfst(fpet.AB)
fpet.AB$pop <- factor(fpet.AB$pop)
fpet_pairwise_fst_bp.A <- boot.ppfst(fpet.AB, nboot = 999)
write.csv(fpet_pairwise_fst_AB, file = "output/Genetic_Stats/fpet_pairwise_WC_fst_lineA.csv")
save(fpet_pairwise_fst_bp.A, file = "output/Genetic_Stats/fpet_pairwise_WC_fst_bp_lineA.RDA")

################################test line B##############################
above <- c("113","95","201","179","204")
below <- c("39","70","96","205")

A <- loct[loct$Site %in% above,]
fpet.A <- fpet[rownames(fpet) %in% A$Sample_ID,]
fpet.A$pop <- fpet.A$pop[fpet.A$pop=="BAJA"] <- "Above"
droplevels()

B <- loct[loct$Site %in% below,]
fpet.B <- fpet[rownames(fpet) %in% B$Sample_ID,]
fpet.B$pop <- fpet.B$pop[fpet.B$pop=="BAJA"] <- "Below"

fpet.AB <- rbind(fpet.A,fpet.B)
remove(above, below, fpet.A, fpet.B,A,B)
fpet_pairwise_fst_B <- pairwise.WCfst(fpet.AB)
fpet.AB$pop <- factor(fpet.AB$pop)
fpet_pairwise_fst_bp.B <- boot.ppfst(fpet.AB, nboot = 999)
write.csv(fpet_pairwise_fst_B, file = "output/Genetic_Stats/fpet_pairwise_WC_fst_lineB.csv")
save(fpet_pairwise_fst_bp.B, file = "output/Genetic_Stats/fpet_pairwise_WC_fst_bp_lineB.RDA")
