#Author: Kevin Quinteros
#Purpose: Calculate expected Heterozygosity for Pegoscapus
library(adegenet)
library(tidyverse)
library(hierfstat)
library(cowplot)
library(ggrepel)
library(magrittr)

######################read in filter Adegenet Object##########################
#Site must have 2 or more individuals
#individuals with more than 50% removed
#Combined site 209 and 214 since they are close in proximity
data <- read.csv(file = "Data/Pego filtered, n≥2, ≤50% NA data.csv")

#remove hybrid individual
data <- data[!data$ind== "P342_petiolaris",]

#Formatting dataset
data <- add_column(data, region = NA, .before = "pop")
data <- add_column(data, cluster = NA, .before = "region")
data <- data %>% mutate(region = case_when(
  pop=="39" ~ "Baja California", # Baja
  pop=="70" ~ "Baja California",
  pop=="95" ~ "Baja California",
  pop=="96" ~ "Baja California",
  pop=="158" ~ "Baja California",
  pop=="172" ~ "Baja California",
  pop=="201" ~ "Baja California",
  pop=="36" ~ "Coastal Sonora",  #Coastal Sonora
  pop=="100" ~ "Coastal Sonora",
  pop=="104" ~ "Coastal Sonora",
  pop=="101" ~ "Coastal Sonora",
  pop=="103" ~ "Inland Sonora", #Inland Sonora
  pop=="106" ~ "Inland Sonora", 
  pop=="107" ~ "Inland Sonora",
  pop=="108" ~ "Inland Sonora", 
  pop=="226" ~ "Coastal Sinaloa", #Coastal Sinaloa
  pop=="230" ~ "Coastal Sinaloa",
  pop=="231" ~ "Coastal Sinaloa", 
  pop=="228" ~ "Inland Sonora", #Inland Sinaloa
  pop=="229" ~ "Inland Sonora", 
  pop=="217" ~ "Zacatecas", #Zacatecas
  pop=="215" ~ "Jalisco", #Jalisco
  pop=="218" ~ "Jalisco",
  pop=="219" ~ "Morelos", #Morelos
  pop=="220" ~ "Morelos", 
  pop=="221" ~ "Northern Oaxaca", #Northern Oaxaca
  pop=="209" ~ "Southern Oaxaca", #Southern Oaxaca
  pop=="214" ~ "Southern Oaxaca", 
  pop=="222" ~ "Southern Oaxaca"
))

data <- data %>% mutate(cluster = case_when(
  pop=="39" ~ "Northern Mexico", # Baja
  pop=="70" ~ "Northern Mexico",
  pop=="95" ~ "Northern Mexico",
  pop=="96" ~ "Northern Mexico",
  pop=="158" ~ "Northern Mexico",
  pop=="172" ~ "Northern Mexico",
  pop=="201" ~ "Northern Mexico",
  pop=="36" ~ "Northern Mexico",  #Coastal Sonora
  pop=="100" ~ "Northern Mexico",
  pop=="104" ~ "Northern Mexico",
  pop=="101" ~ "Northern Mexico",
  pop=="103" ~ "Northern Mexico", #Inland Sonora
  pop=="106" ~ "Northern Mexico", 
  pop=="107" ~ "Northern Mexico",
  pop=="108" ~ "Northern Mexico", 
  pop=="226" ~ "Northern Mexico", #Coastal Sinaloa
  pop=="230" ~ "Northern Mexico",
  pop=="231" ~ "Northern Mexico", 
  pop=="228" ~ "Northern Mexico", #Inland Sinaloa
  pop=="229" ~ "Northern Mexico", 
  pop=="217" ~ "Northern Mexico", #Zacatecas
  pop=="215" ~ "Northern Mexico", #Jalisco
  pop=="218" ~ "Northern Mexico",
  pop=="219" ~ "Southern Mexico", #Morelos
  pop=="220" ~ "Southern Mexico", 
  pop=="221" ~ "Southern Mexico", #Northern Oaxaca
  pop=="209" ~ "Southern Mexico", #Southern Oaxaca
  pop=="214" ~ "Southern Mexico", 
  pop=="222" ~ "Southern Mexico"
))


########################convert to other data formats##########################
#genind object
ade <- adegenet::df2genind(data[5:919], 
                           ploidy = 2, ind.names = data$ind,
                           pop = data$pop, sep = "",
                           type = "codom", NA.char = NA,
                           ncode = 2)
#strata information
strata(ade) <- tibble(data[1:4]) 
#genepop
gen <- genind2genpop(ade)


###############Calculate expected and observed heterozygosity##################
#using the basic.stats function
stats <- basic.stats(data = data[c(3,5:919)], diploid = T) #calculate basic stats 
Avg_Hs <- apply(stats$Hs, 2, mean, na.rm = T) # calculates average expected Hs (na removed)
Avg_Ho <- apply(stats$Ho, 2, mean, na.rm = T) #calculates average observed

#calculate the expected heterozygosity with the adegenet function 
Hs.ade <- Hs(ade) 

#calculate  expected heterozygosity with nason method
#use this output to plot for the paper. 
output <- betas(data[c(3,5:919)]) 
Hs.betas <- rowMeans(output$Hi, na.rm = T)

##########################Read in Fis table####################################
Fis <- read.csv(file = "output/Genetic_Stats/Fis/Fis_pego_nason_approach.csv")

###################Add Expected Heterozygosity to Fis table####################

#change order to match Fis table order of population
Avg_Hs %<>% round(digits = 4)
Avg_Ho %<>% round(digits = 4)
Hs.ade %<>% round(digits = 4)
Hs.betas %<>% round(digits = 4)

#combine different heterozygosit measurements 
Fis <- cbind(Fis, Avg_Ho, Avg_Hs, Hs.ade, Hs.betas)

####################Read in Locality information###############################
lat <- readxl::read_xlsx("Data/Pego_pop_cord.xlsx")
####################Add Latitude and Longitude#################################
lat$`Sampling Site`[lat$`Sampling Site` == "100T"] <- "101" #rename 100T
#remove sites that are not in filtered data
lat <- lat[which(lat$`Sampling Site` %in% Fis$pop),] 
#match order of latitude table to that of Fis table
lat <- lat[match(Fis$pop, lat$`Sampling Site`),]
#combine coordinates table with Fis table
Fis <- cbind(Fis, lat$Longitude, lat$Latitude, lat$`Sampled Wasps`)
names(Fis)[11:13] <- c("Longitude","Latitude","Wasps_sampled")


output_file_name <- "Pego_expected_Hs_and_Fis_filtered_n≥3_≤50%_NA_data.csv"
title <- paste0("Pegoscapus population-specific Fis estimates and 95% confidence limits.\n",
                " - Generating script:  Fis-w-confidence-limits.R and Exp_Hs.R\n",
                " - Input data files: Pego filtered, n≥3, ≤50% NA data.csv, Pego_pop_cord.xlsx, and Fis_pego_nason_approach.csv \n",
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
                " - Wasps_sampled: the actual numer of wasps sampled and sequenced \n")

writeLines(title, output_file_name)

suppressWarnings(write.table(Fis, sep=",", append=TRUE, quote=FALSE, 
                             col.names=TRUE, row.names = FALSE, 
                             file=output_file_name))

###############################add cluster information##########################
#Formatting dataset
Fis <- add_column(Fis, region = NA, .before = "pop")
Fis <- add_column(Fis, cluster = NA, .before = "region")
 Fis <- Fis %>% mutate(region = case_when(
  pop=="39" ~ "Baja California", # Baja
  pop=="70" ~ "Baja California",
  pop=="95" ~ "Baja California",
  pop=="96" ~ "Baja California",
  pop=="158" ~ "Baja California",
  pop=="172" ~ "Baja California",
  pop=="201" ~ "Baja California",
  pop=="36" ~ "Coastal Sonora",  #Coastal Sonora
  pop=="100" ~ "Coastal Sonora",
  pop=="104" ~ "Coastal Sonora",
  pop=="101" ~ "Coastal Sonora",
  pop=="103" ~ "Inland Sonora", #Inland Sonora
  pop=="106" ~ "Inland Sonora", 
  pop=="107" ~ "Inland Sonora",
  pop=="108" ~ "Inland Sonora", 
  pop=="226" ~ "Coastal Sinaloa", #Coastal Sinaloa
  pop=="230" ~ "Coastal Sinaloa",
  pop=="231" ~ "Coastal Sinaloa", 
  pop=="228" ~ "Inland Sinaloa", #Inland Sinaloa
  pop=="229" ~ "Inland Sinaloa", 
  pop=="217" ~ "Zacatecas", #Zacatecas
  pop=="215" ~ "Jalisco", #Jalisco
  pop=="218" ~ "Jalisco",
  pop=="219" ~ "Morelos", #Morelos
  pop=="220" ~ "Morelos", 
  pop=="221" ~ "Northern Oaxaca", #Northern Oaxaca
  pop=="209" ~ "Southern Oaxaca", #Southern Oaxaca
  pop=="214" ~ "Southern Oaxaca", 
  pop=="222" ~ "Southern Oaxaca"
))

Fis <- Fis %>% mutate(cluster = case_when(
  pop=="39" ~ "Northern Mexico", # Baja
  pop=="70" ~ "Northern Mexico",
  pop=="95" ~ "Northern Mexico",
  pop=="96" ~ "Northern Mexico",
  pop=="158" ~ "Northern Mexico",
  pop=="172" ~ "Northern Mexico",
  pop=="201" ~ "Northern Mexico",
  pop=="36" ~ "Northern Mexico",  #Coastal Sonora
  pop=="100" ~ "Northern Mexico",
  pop=="104" ~ "Northern Mexico",
  pop=="101" ~ "Northern Mexico",
  pop=="103" ~ "Northern Mexico", #Inland Sonora
  pop=="106" ~ "Northern Mexico", 
  pop=="107" ~ "Northern Mexico",
  pop=="108" ~ "Northern Mexico", 
  pop=="226" ~ "Northern Mexico", #Coastal Sinaloa
  pop=="230" ~ "Northern Mexico",
  pop=="231" ~ "Northern Mexico", 
  pop=="228" ~ "Northern Mexico", #Inland Sinaloa
  pop=="229" ~ "Northern Mexico", 
  pop=="217" ~ "Northern Mexico", #Zacatecas
  pop=="215" ~ "Northern Mexico", #Jalisco
  pop=="218" ~ "Northern Mexico",
  pop=="219" ~ "Southern Mexico", #Morelos
  pop=="220" ~ "Southern Mexico", 
  pop=="221" ~ "Southern Mexico", #Northern Oaxaca
  pop=="209" ~ "Southern Mexico", #Southern Oaxaca
  pop=="214" ~ "Southern Mexico", 
  pop=="222" ~ "Southern Mexico"
))

#############################prepare to plot####################################
#color palette
Fis <- add_column(Fis, palette = NA, .before = "pop")

#safe_colorblind_palette
Fis <- Fis %>% mutate(palette = case_when(
  region== "Baja California" ~ "#88CCEE", 
  region== "Coastal Sonora" ~ "#CC6677",
  region== "Inland Sonora" ~ "#DDCC77",
  region== "Coastal Sinaloa" ~ "#117733",
  region== "Inland Sinaloa" ~ "#332288",
  region== "Zacatecas" ~ "#AA4499",
  region== "Jalisco" ~ "#44AA99",
  region== "Morelos" ~ "#999933",
  region== "Northern Oaxaca" ~ "#882255",  #Coastal Sonora
  region== "Southern Oaxaca" ~ "#661100",  #Coastal Sonora
))

#reorder the dataframe 
reference <- c("Baja California", "Coastal Sonora","Inland Sonora", 
               "Coastal Sinaloa","Inland Sinaloa", "Zacatecas",
               "Jalisco","Morelos", "Northern Oaxaca", "Southern Oaxaca")

Fis <- Fis %>% arrange(factor(region, levels = reference))

library(cowplot) #plotting themes
###############################plot Hs.betas####################################
#this figure for paper.
#calculate the second expected heterozygosity with the adegent function 
h1.plot <- ggplot() +
  geom_point(data = Fis,
             aes(x = Latitude, y = Hs.betas,
                 shape = cluster, fill = region)
             ,size=3, color = "black") +
  scale_fill_manual(values = Fis$palette, breaks = Fis$region ) +
  scale_shape_manual(values = c(21, 24))+
  xlab("Latitude") + 
  ylab("Expected Heterozygosity")+
  labs(color= "Phylogroups", shape = "Main Pops.") +
  theme_cowplot()+
  theme(aspect.ratio = 1,
        legend.text=element_text(size=9),
        legend.title=element_text(size=9),
        legend.position = c(1,0.5),
        axis.text = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10)) 
plot(h1.plot)

ggsave("Pego_Exp_Hs_and_Lat_Hs_betas.pdf",
       plot = h1.plot, device = "pdf",
       width = 7, height = 4, units = "in", dpi = 400)


