#Author: Kevin Quinteros
#Purpose: Calculate expected Heterozygosity
library(adegenet)
library(tidyverse)
library(hierfstat)
library(cowplot)
library(ggrepel)
###########################set working directory###############################
setwd("~/Projects/Fig_wasp_phylogeography/")
##########################read in Fpetiolaris data#############################
#read in data
#data <- read.table("Data/Fpet_tree/Second_Finn_dataset/output/sub_fpet-s2_trim25-s3_clu_th90-s7_min125_outfiles/sub_fpet-s2_trim25-s3_clu_th90-s7_min125.stru")

#filter data set
#see log for more information: "output/Genetic_Stats/Fpet filtered, n≥3, ≤50% NA log.csv
#At least 3 individual needed per site and removed individuals with more 50% missing data
data <- read.csv(file = "output/Genetic_Stats/data/Fpet filtered, n≥3, ≤50% NA data.csv")

#Add cluster information based off structurre and phylogenetics
data <- add_column(data, cluster = NA, .before = "pop")
data <- data%>%mutate(cluster = case_when(
  pop=="112" ~ 1, # Baja
  pop=="113" ~ 1,
  pop=="158" ~ 1,
  pop=="172" ~ 1,
  pop=="179" ~ 1,
  pop=="201" ~ 1,
  pop=="204" ~ 1,
  pop=="205" ~ 1,
  pop=="39" ~ 1,
  pop=="70" ~ 1,
  pop=="95" ~ 1,
  pop=="96" ~ 1,
  pop=="100" ~ 2, #Coastal Sonora
  pop=="104" ~ 2,
  pop=="103" ~ 3, # Inland Sonora
  pop=="215" ~ 4, # Central Mexico
  pop=="217" ~ 4,
  pop=="210" ~ 5, #Oaxaca
  pop=="214" ~ 5
))
########################convert to other data formats##########################
#genind object
ade <- adegenet::df2genind(data[4:732], 
                           ploidy = 2, ind.names = data$ind,
                           pop = data$pop, sep = "",
                           type = "codom", NA.char = NA,
                           ncode = 2)
#strata information
strata(ade) <- tibble(data[1:3]) 
#genepop
gen <- genind2genpop(ade)
###############Calculate expected and observed heterozygosity##################
#using the basic.stats function
stats <- basic.stats(data = data[c(2,4:732)], diploid = T) #calculate basic stats 
Avg_Hs <- apply(stats$Hs, 2, mean, na.rm = T) # calculates average expected Hs (na removed)
Avg_Ho <- apply(stats$Ho, 2, mean, na.rm = T) #calculates average observed

#calculate the expected heterozygosity with the adegenet function 
Hs.ade <- Hs(ade) 

#calculate  expected heterozygosity with nason method
#use this output to plot for the paper. 
output <- betas(data[c(2,4:732)]) 
Hs.betas <- rowMeans(output$Hi, na.rm = T)

##########################Read in Fis table####################################
Fis <- read.csv(file = "output/Genetic_Stats/Fis/Fis_fpet_nason_approach.csv")
###################Add Expected Heterozygosity to Fis table####################

#change order to match Fis table order of population
Avg_Hs <- Avg_Hs[order(factor(names(Avg_Hs)))] %>% round(digits = 4)
Avg_Ho <- Avg_Ho[order(factor(names(Avg_Ho)))] %>% round(digits = 4)
Hs.ade <- Hs.ade[order(factor(names(Hs.ade)))] %>% round(digits = 4)
Hs.betas <- Hs.betas[order(factor(names(Hs.betas)))] %>% round(digits = 4)

#combine different heterozygosity measurements 
Fis <- cbind(Fis, Avg_Ho, Avg_Hs, Hs.ade, Hs.betas)

####################Read in Locality information###############################
lat <- readxl::read_xlsx("Data/Summary_coord/Fpet_plot_sites.xlsx")
####################Add Latitude and Longitude#################################
lat$`Sample Site`[lat$`Sample Site` == "100T"] <- "100"
lat <- lat[order(factor(lat$`Sample Site`)),]
Fis <- cbind(Fis, lat$Longitude, lat$Latitude, lat$`Samples Trees`)
names(Fis)[11:13] <- c("Longitude","Latitude","Trees_sampled")

output_file_name <- "Fpet_expected_Hs_and_Fis_filtered_n≥3_≤50%_NA_data.csv"
title <- paste0("Ficus petiolaris population-specific Fis estimates and 95% confidence limits.\n",
                " - Generating script:  Fis-w-confidence-limits.R and Exp_Hs.R\n",
                " - Input data files: Fpet filtered, n≥3, ≤50% NA data.csv, Fpet_plot_sites.xlsx, and Fis_fpet_nason_approach.csv \n",
                " - Date: ", Sys.Date(), "\n",
                " - Pop : populations \n",
                " - n : number of samples after filteration \n",
                " - Fis : inbreeding coefficient \n",
                " - n-est: number of loci used in the estimate of Fis \n",
                " - ll: lower confidence interval for Fis \n",
                " - hl: higher confidence interal for Fis \n",
                " - Avg_Ho: average observed heterozygosity using basic.stats() \n",
                " - Avg_Hs: avearge expected heterozygosity using basic.stats() \n",
                " - Hs.ade: expected heterozygosity using the Hs() function \n",
                " - Hs.betas: expected heterozygosity using betas() function (was used in the paper)\n",
                " - Longitude: Longitude \n",
                " - Latitude: Latitude \n",
                " - Samples Tree: the actual numer of Trees sampled and sequenced \n")

writeLines(title, output_file_name)

suppressWarnings(write.table(Fis, sep=",", append=TRUE, quote=FALSE, 
                             col.names=TRUE, row.names = FALSE, 
                             file=output_file_name))

#########################Add Cluster information###############################
#Add cluster information based off structurre and phylogenetics
Fis <- add_column(Fis, cluster = NA, .before = "pop")
Fis <- Fis%>%mutate(cluster = case_when(
  pop=="112" ~ "Baja California", # Baja
  pop=="113" ~ "Baja California",
  pop=="158" ~ "Baja California",
  pop=="172" ~ "Baja California",
  pop=="179" ~ "Baja California",
  pop=="201" ~ "Baja California",
  pop=="204" ~ "Baja California",
  pop=="205" ~ "Baja California",
  pop=="39" ~ "Baja California",
  pop=="70" ~ "Baja California",
  pop=="95" ~ "Baja California",
  pop=="96" ~ "Baja California",
  pop=="100" ~ "Coastal Sonora", #Coastal Sonora
  pop=="104" ~ "Coastal Sonora",
  pop=="103" ~ "Inland Sonora", # Inland Sonora
  pop=="215" ~ "Central Mexico", # Central Mexico
  pop=="217" ~ "Central Mexico",
  pop=="210" ~ "Oaxaca", #Oaxaca
  pop=="214" ~ "Oaxaca"
))

library(cowplot) #plotting theme
#################################Plot the Hs.betas##############################
#This is the figure 
h.plot <- ggplot() +
  geom_point(data = Fis,
             aes(x = Latitude, y = Hs.betas, fill = cluster),
             size=3, shape=21, color = "black") +
  scale_fill_manual(values = c("#56B4E9","#009E73","#F0E442",
                               "#0072B2","#D55E00","#CC79A7"),
                    breaks= c("Baja California","Coastal Sonora",
                              "Inland Sonora","Central Mexico",
                              "Oaxaca")) +
  xlab("Latitude") + 
  ylab("Expected Heterozygosity")+
  labs(fill= "Regional Sites") +
  theme_cowplot()+
  theme(aspect.ratio = 1,
        legend.text=element_text(size=9),
        legend.title=element_text(size=9),
        legend.position = c(1,0.5),
        axis.text = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10)) 
plot(h.plot)

ggsave("output/final_figures/Figure5/Fpet_Exp_Hs_and_Lat_nason_function.pdf",
       plot = h.plot, device = "pdf",
       width = 6, height = 4, units = "in", dpi = 300)
