

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
#####                              Overview                                #####
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #

# Fis-w-confidence-limits.R
# John Nason
# September 3, 2022

# This script estimates Wright's Fis with confidence limits.

# Wright's F-statistics: Overview:

# Wright's F-statistics have multiple interpretations, including measures of 
# inbreeding.
# Fis: measures the average inbreeding of individuals (I) within subpopulations 
#     (S). 
# Fst: measures the average inbreeding within subpopulations (S) relative to the 
#     total population (T).
# Fit: measures the average inbreeding of individuals (I) relative to the total 
#     population (T).

# In many studies, there is an additional hierarchical level in the data, with
# subpopulations (or populations) nested within regions (or similar). The example 
# we will work with below includes this additional hierarchical level, with 
# individuals nested within patches, which are nested within localities, which 
# are nested within the total population. 

# Here we focus on the estimation of population-specific estimates of Fis. For
# a tutorial estimating Wright's F-statistics in general, see 
# Fstats_w_Confidence_Limits.R.

# Importantly, Fis can only be interpreted as local inbreeding when calculated 
# at the local population level, which in the example below is at the level of 
# Patch. At this level, Fis summarizes deviations from Hardy-Weinberg 
# equilibrium due to non-random mating, with values greater than zero 
# interpretable as inbreeding between relatives. However, if populations are 
# pooled and then Fis calculated, this estimate of Fis will include deviations 
# from HWE attributable to both local inbreeding and genetic differentiation 
# among populations. Consequently, the pooling of populations is not advisable 
# in estimating Fis since then its biological interpretation is no longer #
# straight forward.

# That said, it is possible to properly obtain the average Fis of populations 
# nested within regions using Weir and Cockerham's estimators of Wright's 
# F-statistics. This properly approach accounts for distinct populations within 
# regions as opposed to pooling them.


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
#####                      1. Prepare environment                          #####
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #

# Clean the R environment.
rm(list = ls())

# Load required libraries.
library(rstudioapi)
library(hierfstat)
library(tidyverse)


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
#####           2. Set working environment and source scripts              #####
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #

# Get the path to this script.
main_script_path <- dirname(rstudioapi::getSourceEditorContext()$path)

# Set the current working directory to the folder PopToolsJN and access 
# the following R files located within that directory.
#helper_script_path <- "Fig_wasp_phylogeography/"
#setwd(helper_script_path)
source("get-file-name-path.R")


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
#####                    3. Import input data file                         #####
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #

# Here are file paths to unfiltered and filtered Fpet pop gen data files. Note 
# that this filtered file is in .rds format. There is also a .csv version 
# available.

# Fpet filtered:
input_file_path <- paste0("Fig_wasp_phylogeography/output/Genetic_Stats/data/Fpet filtered, n≥3, ≤50% NA data.csv")

#Pego filtered: Remove P342 since its a hybrid
input_file_path <- paste0("Fig_wasp_phylogeography/output/Genetic_Stats/data/Pego filtered, n≥2, ≤50% NA data.csv")
# get.file.name() and get.file.path() are functions in custom R file 
# get-file-name-path.R.
#input_file_name <- get.file.name(input_file_path)
#input_file_dir <- get.file.dir.path(input_file_path, endchar="/")

# Use this to read the input data file if in .csv format.
input_df <- read.csv(file = input_file_path)

#if running the pegoscapus dataset remove P342 since it's a hybrid
input_df <- input_df[!input_df$ind == "P342_petiolaris",]

# Use this to read the input data file if in .rds format.
#input_df <- readRDS(input_file_path)

# The total number of samples
length(input_df[,1])

# Number of columns.
length(input_df[1,])

# View the data frame to check contents and understand its structure.
View(input_df)

# Check the class of the data in the first 10 columns in the data frame. Confirm
# the pop names, their order, and that the genotypes are integer.
str(input_df[,c(1:10)])

# In column pop, change the name of 100-Tetas to 100.
input_df["pop"][input_df["pop"] == "100-Tetas"] <- "100"

View(input_df)

# Delete column 2, which contains individual sample IDs.
input_df <- input_df[,-2]
View(input_df)
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
#####             4. Estimate Population-specific Fis and CLs              #####
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #

# Population-specific Fis estimates and confidence limits can be obtained using 
# the functions basic.stats and boot.ppfis. Fis estimated by locus and grouping, 
# which for data set gtrunchier are Locality and Patch. Note that values of NaN 
# are returned for loci that are monomorphic or have no data.

# The function basic.stats returns a variety of genetic measures at the per
# locus and overall levels.
bs <- basic.stats(input_df)

# Kevin's approach for obtaining population-specific Fis estimates.
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
# For each population, calculate mean Ho over loci. Ditto for Hs.
Ho <- bs$Ho %>% colMeans(na.rm = T) #average site observed heterozygosity per site
Hs <- bs$Hs %>% colMeans(na.rm = T) #average site expected heterozygosity per site
Fis.mean <- 1 - Ho/Hs #average Fis per site 
Ho
Hs
Fis.mean

#John's approach
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
# From these basic stats, extract estimates of Fis by locus (row) and
# population (column).
locus_in_pop_Fis <- bs$Fis
locus_in_pop_Fis

# We can now obtain population-level Fis estimates as means over loci.
pop_Fis <- colMeans(locus_in_pop_Fis, na.rm = TRUE, dims = 1)
pop_Fis
grand_mean_Fis <- mean(pop_Fis)
pop_mean_Fis <- colMeans(locus_in_pop_Fis, na.rm = T)

#compare 
pop_Fis
pop_mean_Fis
Fis.mean
# Next, we estimate 95% confidence intervals for pop-level Fis estimates. The 
# results are returned in a list with the lower limits in column ll and upper 
# limits in column hl. Note that values of NA are returned for loci that are 
# monomorphic or have no data.
Fis_CIs <- boot.ppfis(dat=input_df, nboot=1000, quant=c(0.025,0.975), diploid=TRUE)
Fis_CIs
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #


# Now organize the Fis estimates and associated information in a data frame.
# First, extract a table of population samples sizes as a data frame.
pop_table_df <- as.data.frame(table(input_df$pop))
colnames(pop_table_df) <- c("pop", "n")
pop_table_df

# Append columns for population sample size and Fis estimates.
r_pop_Fis <- pop_Fis %>% round(digits = 4)
pop_table_df <- cbind(pop_table_df, r_pop_Fis)
colnames(pop_table_df) <- c("pop", "n", "Fis")
row.names(pop_table_df) <- 1:length(pop_table_df[,1])
pop_table_df

# Determine for each population the number of loci contributing to it's Fis
# estimate.
n_pops <- length(unique(input_df[,1]))
n_loci <- ncol(input_df)-1
# Create data frame to hold pop names and number of individual locus Fst 
# estimates per population.
Fis_ests_per_pop_df <- data.frame(matrix(NA,    # Create empty data frame
                          nrow = n_pops,
                          ncol = 2))
colnames(Fis_ests_per_pop_df) <- c("pop", "Fis_ests")
c <- 0
for (pop in colnames(locus_in_pop_Fis)){
  c <- c+1
  Fis_ests <- n_loci - sum(sapply(locus_in_pop_Fis[,c], function(x) sum(is.na(x))))
  Fis_ests_per_pop_df[c,1] <- pop
  Fis_ests_per_pop_df[c,2] <- Fis_ests
}
Fis_ests_per_pop_df

# Append the column of Fis estimates per pop.
pop_table_df <- cbind(pop_table_df, Fis_ests_per_pop_df[,2])
colnames(pop_table_df) <- c("pop", "n", "Fis", "n_ests")
pop_table_df
mean(pop_table_df$Fis)

# Append the lower (ll) and upper (hl) confidence limits on Fis to the table.
# Note that Fis_CIs is a list with Fis_CIs[[1]] being the function call and 
# Fis_CIs[[2]] being a data frame containing the lower and upper confidence 
# limits.
pop_table_df <- cbind(pop_table_df, Fis_CIs[[2]])
pop_table_df

#write to file (depending on which dataset you are running, choice one)
write.csv(pop_table_df, "Fig_wasp_phylogeography/output/Genetic_Stats/Fis/Fis_fpet_nason_approach.csv", row.names = F)
write.csv(pop_table_df, "Fig_wasp_phylogeography/output/Genetic_Stats/Fis/Fis_pego_nason_approach.csv", row.names = F)

# Nei estimator of Fis used by hierfstat is positively biased by amount 
# 1/((2n)-1). Therefore, try adjusting Fis downward by amount 1/((2n)-1).
# This is a crude way of doing it since it assumes that for a population of
# sample size n, n genotypes are available at every locus, which is definitely 
# not the case.
pop_table_df$Fis
mean(pop_table_df$Fis)
pop_table_adj_df <- pop_table_df
pop_table_adj_df$Fis <- pop_table_adj_df$Fis - (1/((2*pop_table_adj_df$n)-1))
pop_table_adj_df$Fis
mean(pop_table_adj_df$Fis)

# Barplot of Fis values across populations.
# main_label <- paste0("Fis across populations (mean = ", format(round(grand_mean_Fis, 5), nsmall=5), ")")
# barCenters <- barplot(pop_table_df$Fis, las=2, main=main_label, xlab="Population",
#                       ylab="Fis", ylim = c(min(pop_table_df$ll), max(pop_table_df$hl)), names.arg = unique(pop_table_df[,1]), cex.names=1) 
# arrows(x0=barCenters, y0=pop_table_df$hl, x1=barCenters, y1=pop_table_df$ll, angle=90, length=0.075, code=3)
# abline(h=grand_mean_Fis, col="black", lwd=1.2, lty=2)


# To visualize the results, create a barplot of Fis values across populations.
library(ggplot2)
library(scales) # provides labels = comma, which converts the 
                # y-axis values from scientific to standard notation.
for (c in 1){
  main_label <- paste0("\nFis across populations (mean = ", 
                       format(round(grand_mean_Fis, 5), nsmall=5), ")\n")
  lower_limit <- floor(min(pop_table_df$ll)*10)/10
  lower_limit
  upper_limit <- ceiling(max(pop_table_df$hl)*10)/10
  upper_limit
}
p <- ggplot(pop_table_df, aes(x=pop, y=Fis)) +
  geom_bar(stat="identity", color="black", fill="lightgray") +
  geom_hline(yintercept=grand_mean_Fis,linetype=2) +
  geom_errorbar(aes(ymin=ll, ymax=hl), width=.2,
                position=position_dodge(.9)) +
  geom_text(aes(label=n_ests), vjust=1.6, color="black", size=3.5) +
  labs(title=main_label, x="\nPopulation\n", y="    Fis   ") +
  scale_y_continuous(breaks = seq(lower_limit, upper_limit, by = 0.1), 
                     limits=c(lower_limit, upper_limit), labels = comma) +
  theme_void() # No plot background, border, or ticks. Will add below.
p + theme(axis.text = element_text(size = 12), 
          axis.text.x = element_text(angle = 90, vjust=0.5, margin = margin(t = 4)), # t is for top
          axis.text.y = element_text(margin = margin(r = 4)), # r is for right
          axis.title = element_text(size = 14),
          plot.title = element_text(size = 16, hjust=0.5),
          panel.border = element_rect(color = "black",
                                      fill = NA,
                                      size = 0.5),
          axis.ticks = element_line(colour = 'black', size = 0.5),
          axis.ticks.length = unit(1, "mm"))


# To get a sense of the mean and variance of a population-specific Fis estimate, 
# plot a histogram of the frequency of locus-specific Fis estimates for that
# population.

# Indicate a specific population.
pop_name = "214"
for (c in 1){
  temp <- locus_in_pop_Fis[,pop_name]
  mean_Fis <- pop_mean_Fis[pop_name]
  title <- paste0("Individual-locus Fis estimates in population ", pop_name, 
                  "\n(mean = ", format(round(mean_Fis, 5), nsmall=5),
                  ", median = ", median(temp, na.rm=TRUE), 
                  ", n = ", length(temp), ")")
  hist(temp, main=title, xlab="Fis estimates", right=FALSE)
  abline(v=mean_Fis, col="black", lwd=1, lty=2)
}


# A plot of Fis versus sample size (n).
plot(x=pop_table_df$n, y=pop_table_df$Fis, main="Fis versus sample size (n)",
     xlab="n", ylab="Fis")
abline(lm(pop_table_df$Fis ~ pop_table_df$n))

#this method was not used in the paper. 
#plot(x=pop_table_df$Fis, y=Fis.mean, main="Population-level Fis_Quinteros versus Fis_Nason",
#     xlab="Fis_Nason", ylab="Fis_Quinteros", xlim=c(-0.1,0.1), ylim=c(-0.1,0.15))
#text_label <- "Fis_Quinteros = 1 - Ho_mean/He_mean\nFis_Nason = mean(Fis_locus)"
#text(x=-0.05, y=0.13, label=text_label)
#abline(coef = c(0, 1))


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
#####              4. Calculate population-level Fis "by hand"             #####
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #

# Here we estimate Fis by hand directly from the data using the estimator on
# page 65 from Weir (1996).
# - The estimator from Weir returns Fis = 0 when the sample data for a locus is
#   in HWE, that is when the observed equals the expected heterozygosity.
# - hierfstat uses the estimator of Fis from equation 7.41 on page 164 of Nei
#   (1987). This estimator returns Fis = 1/(2n-1) when the sample data for a 
#   locus is in HWE. This is because it aims to estimate the Fis for the natural
#   population from which the data are obtained and in natural populations with
#   random mating the expectation for Fis is 1/(2Ne-1), where Ne is the
#   population size. So the Nei estimator assumes 1/(2n-1) is a good 
#   approximation of 1/(2Ne-1), that is, that the sample size n is a good 
#   approximation of the effective size Ne. This seems a stretch to me and may
#   be a fair approximation when n is fairly large, but will be a very poor
#   approximation when a small sample is obtained from a population of large Ne,
#   that is when n << Ne. This is very likely to be the case with data sets
#   generated by next gen sequencing, where the number of individuals sampled 
#   population is often less than 10.

# Example with Ne = 100
# For example, consider a natural population with Ne = 100 that is subject to 
# random mating. In this case, the actual Fis in the source population should 
# be 1/(2Ne-1) = 1/99 = 0.01. For comparison to the actual Fis, below we plot 
# Fis estimates using Nei (1987) and Weir (1996) where the sample data are in 
# HWE. You will see that at Ne = 100, Nei rather dramatically over-estimates the 
# when the sample n is actual Fis is small. For example if the sample size is n,
# then Fis_Nei87 = 1/(2n-1) = 0.11, which is 11x greater than the actual Fis!. 
# In contrast, Weir underestimates the actual Fis but only by a small amount, 
# with a bias of -1/(2Ne-1) = 0.01.

# Example with Ne = 10
# An Ne this small is unlikely to occur in nature for most organisms, but it 
# illustrates the point that if Ne is very small, then the Nei estimator is a
# better choice for n >= 6, but Weir is better for n <= 4. Since this is a wash
# and Weir is better when Ne is of a more typical size, I'd go with Weir.

Ne <- 100
n <- 2:Ne
Fis_Nei87 <- 1/(2*n-1)
Fis_Weir96 <- rep.int(0, Ne-1)
Fis_actual <- 1/(2*Ne-1)
xlim_max <- min(Ne,30)
# Plot Fis_Nei87.
plot(x=n, y=Fis_Nei87, main="Population-level Fis estimates for samples at HWE",
     xlab="Sample size (n)", ylab="Estimated Fis", pch=1, type="b", xlim=c(2,xlim_max), ylim=c(-0.01,max(Fis_Nei87)))
# Add Fis_Weir96.
lines(x=n, y=Fis_Weir96, pch=20, col="black", type="b")
# Add Fis_actual.
lines(x=n, y=rep(Fis_actual,Ne-1), lty=2, col="black")
# Add a legend.
legend(x=round(0.78*xlim_max), y=0.98*max(Fis_Nei87), legend=c("Nei (1987)", "Weir (1996)", "Actual"),
       col=c("black", "black", "black"), pch=c(1,20,NA), lty=c(1,1,2), cex=0.9)

# Kevin's approach for obtaining population-specific Fis estimates from 
# hierfstat, which is equivalent to the Nei (1987) estimator of Fis. For each 
# population, calculate mean Ho and Hs over loci. This estimate of He is the 
# same as in Nei (1987). Then estimate Fis as 1 - (Ho/Hs). In contrast, above
# let hierfstat estimate Fis per locus per population, and then I take the 
# average of the locus-specific Fis estimates to obtain a population-level
# estimate.
Ho <- bs$Ho %>% colMeans(na.rm = T) #average site observed heterozygosity per site
Hs <- bs$Hs %>% colMeans(na.rm = T) #average site expected heterozygosity per site
Fis.mean <- 1 - Ho/Hs #average Fis per site 
Fis.mean
mean(Fis.mean)


# Now we calculate pop-level Fis using equation 2.21 on page 64 of Weir's 
# "Genetic Data Analysis" (1996).

# The minimum and maximum alleles are 1 and 4:
min(input_df[,-1], na.rm = T)
max(input_df[,-1], na.rm = T)

# Here we go!
pop_Fis_Weir_vec <- as.numeric()
pop_Fis_Nei_vec <- as.numeric()
n_actual <- as.numeric()
for (pop in 1:length(pop_table_df[,1])){
  # Obtain the population data
  #pop <- 1
  writeLines(paste0(pop_table_df[pop,1]))
  pop_df <- input_df[input_df[,1] == pop_table_df[pop,1], ]
  #View(pop_df)
  locus_Fis_Weir_vec <- as.numeric()
  locus_Fis_Nei_vec <- as.numeric()
    for (loc in 2:ncol(input_df)){
      #writeLines(paste0("locus = ", loc))
      geno_mat <- matrix(0, nrow = 4, ncol = 4)
      n <- 0
      alleles <- as.numeric()
      for (ind in 1:nrow(pop_df)){
        if (!is.na(pop_df[ind,loc])){
          n <- n + 1
          a1 <- pop_df[ind,loc] %% 10 
          a2 <- pop_df[ind,loc] %/% 10
          alleles <- c(alleles, a1, a2)
          if (a1 == a2){
            geno_mat[a1,a1] <- geno_mat[a1,a1] + 1
          } else{
            geno_mat[a1,a2] <- geno_mat[a1,a2] + 1
            geno_mat[a2,a1] <- geno_mat[a2,a1] + 1
          }
        }
      } # end individual
      n_actual <- c(n_actual, n)
      #print(geno_mat)
      #print(alleles)
      # Calculate Fis for each allele
      allele_Fis_Weir_vec <- as.numeric()
      allele_Fis_Nei_vec <- as.numeric()
      if (length(unique(alleles)) ==1){
        Fis_Weir <- NaN
        Fis_Nei <- NaN
      } else
      for (al in unique(alleles)){
        geno_mat2 <- geno_mat
        x_ua <- unique(alleles)
        x_al <- al
        pop_size <- nrow(pop_df)
        x_n <- n
        x_loc <- loc
        
        # Weir.
        p <- (1/(2*n))*(sum(geno_mat[al,])+geno_mat[al,al]) # Allele frequency: count the homozygote twice
        n_het <- sum(geno_mat[al,]) - geno_mat[al,al] # Number of heterozygotes
        Fis_Weir <- 1 - (n_het)/(2*n*p*(1-p)) 
        Fis_Weir <- Fis_Weir 
        allele_Fis_Weir_vec <- c(allele_Fis_Weir_vec, Fis_Weir)
        
        # Nei.
        homo_freq <- 1 - p^2 - (1-p)^2 # The sum of the expected homozygote frequencies
        Ho <- (n-n_het)/n  # The observed frequency of heterozygotes
        Hs <- (n/(n-1))*(1 - homo_freq - (Ho/(2*n))) # The expected frequency of heterozygotes
        Fis_Nei <- 1 - (Ho/Hs)
        allele_Fis_Nei_vec <- c(allele_Fis_Nei_vec, Fis_Nei)
      } # end allele
      if (length(allele_Fis_Weir_vec) > 0){
        locus_Fis_Weir_vec <- c(locus_Fis_Weir_vec, mean(allele_Fis_Weir_vec, na.rm = T))
      }
      if (length(allele_Fis_Nei_vec) > 0){
        locus_Fis_Nei_vec <- c(locus_Fis_Nei_vec, mean(allele_Fis_Nei_vec, na.rm = T))
      }
    } # end locus
  pop_Fis_Weir_vec <- c(pop_Fis_Weir_vec, mean(locus_Fis_Weir_vec, na.rm = T))
  pop_Fis_Nei_vec <- c(pop_Fis_Nei_vec, mean(locus_Fis_Nei_vec, na.rm = T))
} # end pop

pop_Fis_Weir_vec
mean(pop_Fis_Weir_vec)

pop_Fis_Nei_vec
mean(pop_Fis_Nei_vec)

# The approximate bias in the Nei Fis estimate relative to a random mating 
# source population with large Ne:
m <- mean(n_actual)
Nei_bias <- 1/((2*m)-1)
Nei_bias

# For interest, subtract off the bias from the mean Nei Fis estimate over 
# populations.
mean(pop_Fis_Nei_vec)-Nei_bias

# Consider tallying up number of Fis estimates per population...

# Make a data frame with columns population name, size, and Fis_Weir by 
# population.
pop_table_Weir <- cbind(pop_table_df[,1:2], pop_Fis_Weir_vec)
colnames(pop_table_Weir) <- c("pop", "n", "Fis_Weir")
View(pop_table_Weir)

# Note that the following plotting code is copied from above. Need to put it in 
# a function...
# To visualize the results, create a barplot of Fis values across populations.
library(ggplot2)
library(scales) # provides labels = comma, which converts the 
# y-axis values from scientific to standard notation.
for (c in 1){
  main_label <- paste0("\nFis_Weir across populations (mean = ", 
                       format(round(mean(pop_Fis_Weir_vec), 5), nsmall=5), ")\n")
  # lower_limit <- floor(min(pop_table_df$ll)*10)/10
  # lower_limit
  # upper_limit <- ceiling(max(pop_table_df$hl)*10)/10
  # upper_limit
}
p <- ggplot(pop_table_Weir, aes(x=pop, y=Fis_Weir)) +
  geom_bar(stat="identity", color="black", fill="lightgray") +
  geom_hline(yintercept=mean(pop_Fis_Weir_vec),linetype=2) +
  # geom_errorbar(aes(ymin=ll, ymax=hl), width=.2,
  #               position=position_dodge(.9)) +
  #geom_text(aes(label=n_ests), vjust=1.6, color="black", size=3.5) +
  labs(title=main_label, x="\nPopulation\n", y="    Fis   ") +
  # scale_y_continuous(breaks = seq(lower_limit, upper_limit, by = 0.1), 
  #                    limits=c(lower_limit, upper_limit), labels = comma) +
  theme_void() # No plot background, border, or ticks. Will add below.
p + theme(axis.text = element_text(size = 12), 
          axis.text.x = element_text(angle = 90, vjust=0.5, margin = margin(t = 4)), # t is for top
          axis.text.y = element_text(margin = margin(r = 4)), # r is for right
          axis.title = element_text(size = 14),
          plot.title = element_text(size = 16, hjust=0.5),
          panel.border = element_rect(color = "black",
                                      fill = NA,
                                      size = 0.5),
          axis.ticks = element_line(colour = 'black', size = 0.5),
          axis.ticks.length = unit(1, "mm"))

# A plot of Fis_Weir versus sample size (n).
plot(x=pop_table_Weir$n, y=pop_table_Weir$Fis_Weir, main="Fis_Weir versus sample size (n)",
     xlab="n", ylab="Fis")
abline(lm(pop_table_Weir$Fis_Weir ~ pop_table_Weir$n))
temp
# Plot Fis_Nei versus the Fis_Weir. Here, Fis_Nei is the Fis I obtained from 
# and what I call Fis_Nason, above. 
axis_min <- min(pop_table_Weir$Fis_Weir, pop_table_df$Fis)
axis_max <- max(pop_table_Weir$Fis_Weir, pop_table_df$Fis)
plot(x=pop_table_Weir$Fis_Weir, y=pop_table_df$Fis, main="Fis_Nason versus Fis_Weir",
     xlab="Fis_Weir", ylab="Fis_Nei", xlim=c(axis_min,axis_max), ylim=c(axis_min,axis_max))
abline(coef = c(0, 1))
text(x=0.05, y=-0.075, label=paste0("• Fis_Nason is the Nei (1987)\nestimate from hierfstat\n\n",
                                  "• Fis_Weir is the hand-coded\nWeir (1996) estimate"))


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
#####                     4. Export results to file                        #####
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #

# Save data frame pop_Fis_df to file in the current working directory. We use 
# write.table here instead of write.csv since the latter does not allow 
# appending to an existing file.

# Set the current working directory to the folder containing this script so that  
# the export files will be located there.
setwd(main_script_path)

#produced same results as method above. 
output_file_name <- "Pego Pop Fis filtered.csv" #used in the paper

# First, create a title and write to new file "Fpet Pop Fis filtered.csv".
title <- paste0("Pego population-specific Fis estimates and 95% confidence limits.\n",
                " - Generating script:  Fis-w-confidence-limits.R\n",
                " - Input data file:  Pego filtered, n≥2, ≤50% NA data.csv \n",
                " - Date: ", Sys.Date(), "\n")
writeLines(title)
writeLines(title, output_file_name)

# Now append the data frame to the file. Here write.table() is wrapped in 
# suppressWarnings() because R otherwise generates an error when you try and
# append column names to a file as it normally expects them to be at the top
# of a file. Here though we want a descriptive title at the top.
suppressWarnings(write.table(pop_table_df, sep=",", append=TRUE, quote=FALSE, 
            col.names=TRUE, row.names = FALSE, 
            file=output_file_name))

