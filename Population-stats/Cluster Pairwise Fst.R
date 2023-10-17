
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
#####                             Overview                                 #####
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #

# John Nason
# July 4, 2022

# This code takes a Ficus petiolaris SNP data set and estimates Fst and its
# 95% confidence limits between pairs of "genetic clusters". hierfstat functions
# ppfst() and boot.ppfst() are used to estimate Fst and 95% CIs, respectively.

# The input data is expected to be formatted as follows:
#   Col 1: population names
#   Col 2: individual sample names
#   Remaining cols: genotypes (integer genotypes with missing values coded NA)
# 
# The inferred genetic clusters are assumed to be:
#   Baja: Sites 112, 113, 158, 172, 179, 201, 204, 205, 39, 70, 95, 96
#   Coastal Sonora: Sites 100-Tetas, 104
#   Inland Sonora: Site 103
#   Central Mexico: 215, 217
#   Oaxaca: Sites 210, 214


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
#####                         Prepare environment                          #####
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #

# Clean the R environment.
rm(list = ls())

# Load required libraries.
library(tibble)
library(dplyr)
library(hierfstat)
 
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
#####                       Import input data file                         #####
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #

# Read csv data file into dataframe object.
Fpet.df <- read.csv(file = "~/Projects/Fig_wasp_phylogeography/output/Genetic_Stats/hf_fpet_sub.csv")

# View the data frame to check check contents and understand its structure.
#View(Fpet.df)

# Check the data type of the first 10 columns in the data frame.
str(Fpet.df[,c(1:10)])

# The first column denotes the population name, the second column the individual 
# ID, some of which contain alphabetic characters, and the remaining columns 
# are loci with their genotypes. The genotypes must be integer.

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
#####                       Format input data file                         #####
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #

# We will estimate F-statistics and CIs using the package hierfstat and will 
# need to prep the Fpet.df. for this. Specifically, the first column of the 
# dataframe will need to contain a numeric value indicating for each population 
# the genetic cluster to which it has been assigned based on prior, independent 
# analyses, such as PCA or STRUCTURE.

# Add a new first column to Fpet.df to contain the cluster information. 
# Here add_column() is from tibble.
Fpet.df <- add_column(Fpet.df, cluster = NA, .before = "pop")
#View(Fpet.df)

# Assign cluster information for each population. Here mutate() is from dplyr.
Fpet.df <- Fpet.df%>%mutate(cluster = case_when(
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
  pop=="100-Tetas" ~ 2, #Coastal Sonora
  pop=="104" ~ 2,
  pop=="103" ~ 3, # Inland Sonora
  pop=="215" ~ 4, # Central Mexico
  pop=="217" ~ 4,
  pop=="210" ~ 5, #Oaxaca
  pop=="214" ~ 5
))
#View(Fpet.df)

# For printing purposes below, make a list of the names of the clusters, ordered
# as above..
clustername <- c("Baja", 
                 "Coastal Sonora", 
                 "Inland Sonora", 
                 "Central Mexico", 
                 "Southern Oaxaca")
Fpet.df <- Fpet.df[-which(is.na(Fpet.df$cluster)),]

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
#####                   User-defined filtering settings                    #####
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #

# First, set the maximum proportions of NA values (missing genotypes) allowable   
# for each locus and for individual sample. Loci and individuals exceeding these  
# proportions will be filtered from the dataset. Note that the run time of the
# code below will increase as ind.na and locus.na are increased.
ind.NA <- 0.5
locus.NA <- 0.5


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
#####                       Analyze input data file                        #####
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #

# A for loop is used here to subset pair of genetic clusters from Fpet.df,
# filter them for the number of missing genotypes (NAs) per locus and individual
# sample, and then estimate Fst and its 95% confidence limits. Results are 
# output for each pair of genetic clusters.

results.mat = matrix(data=NA, nrow = 5, ncol = 5)
rownames(results.mat) <- c("Baja",
                           "Coastal Sonora", 
                           "Inland Sonora", 
                           "Central Mexico", 
                           "Southern Oaxaca")
colnames(results.mat) <- c("Baja",
                           "Coastal Sonora", 
                           "Inland Sonora", 
                           "Central Mexico", 
                           "Southern Oaxaca")


for (c1 in unique(Fpet.df$cluster)) {
  for (c2 in unique(Fpet.df$cluster)[(c1+1):5]) {
    
    # Here we are iterating through genetic clusters c1 and c2 with c2 > c1.
    
    # We'll track the run time and when finished print the total elapsed time 
    # in the for loop.
    if (c1==1 & c2==2){
      start_time <- Sys.time()
    }
    
    # Subset the Fpet.df into the current pair of clusters, c1 and c2.
    df <- subset(Fpet.df, cluster==c1 | cluster==c2)
    
    # Print a header for each cluster-wise comparison.
    writeLines("###############################################\n")
    print(paste(clustername[c1], "versus", clustername[c2]), quote = FALSE)
    writeLines("")
    
    # hierfstat expects clusters to be number sequentially beginning with 1, so
    # here, for each pair, we renumber c1 to 1 and c2 to 2.
    df["cluster"][df["cluster"] == c1] <- 1
    df["cluster"][df["cluster"] == c2] <- 2
    
    # Filter the subsetted data, retaining loci at which the proportion of 
    # missing loci (NAs) is no greater than the maximum locus.NA.
    filtered.df <- df
    for (x in unique(filtered.df$pop)) {
      temp <- subset(filtered.df, pop==x)
      good_loci <- names(which(colSums(is.na(temp)) <= round(locus.NA*length(temp[,1]))))
      filtered.df <- filtered.df[ , names(filtered.df) %in% good_loci]
    } 
    
    # Set the max number of NA values (missing genotypes) allowable per 
    # individual.
    filtered.df.ind.NA <- round(ind.NA * length(filtered.df[1,-c(1:3)])) 
    filtered.df.ind.NA
    
    # Filter the subsetted data, retaining loci at which the number of missing
    # loci (NAs) is no greater than the maximum filtered.df.ind.NA.
    filtered.df <- filtered.df[rowSums(is.na(filtered.df)) <= filtered.df.ind.NA, ]
    
    print(paste("Number of individuals post filtering:", length(filtered.df[,1])), quote = FALSE)
    print(paste("Number of loci post filtering:       ", length(filtered.df[1,-c(1:3)])), quote = FALSE)
    writeLines("")
    
    # The final filtering step, check that there aren't any populations having 
    # loci that are all NA. If there are, then the program terminates with an 
    # error indicating the problematic population and loci.
    for (x in unique(filtered.df$pop)) {
      temp <- subset(filtered.df, pop==x)
      allNAcols <- apply(temp, 2, function(x)all(is.na(x))) 
      colswithallNA <-names(allNAcols[allNAcols>0])
      if (length(colswithallNA) > 0) {
        print(paste("Error: In population", x, ", the following columns are all NA:", colswithallNA))
        print("Program terminated.")
        break
      }
    } 
    
    # Estimate Fst for the current pair of clusters.
    pp <- pp.fst(dat=filtered.df[,-c(2:3)], diploid=TRUE)
    writeLines(paste("Pairwise Fst =", round(pp$fst.pp[1,2], digits=5)))
    
    # Estimate 95% CIs on Fst for the current pair of clusters.
    boot.pp <- boot.ppfst(dat=filtered.df[,-c(2:3)], nboot=100, quant=c(0.025,0.975), diploid=TRUE)
    writeLines(paste("Upper and lower 95% CIs:", 
                     round(boot.pp$ul[1,2], digits=5), ",",
                     round(boot.pp$ll[1,2], digits=5)))
    writeLines("")
    
    results.mat[c1,c2] <- as.character(round(pp$fst.pp[1,2], digits=5))
    results.mat[c2,c1] <- paste(as.character(round(boot.pp$ul[1,2], digits=5)), ",", as.character(signif(boot.pp$ll[1,2], digits=5)))
    
    if (c1==4 & c2==5) {
      writeLines("###############################################\n")
      print("Results in matrix form, with Fst estimates in the upper")
      writeLines("triangle and upper and lower 95% CLs in the lower triangle.\n")
      print(results.mat)
      writeLines("")
      end_time <- Sys.time()
      print(paste("For loop runtime:", round(end_time - start_time, digits=5)), quote = FALSE)
    }

  }
} 

# The matrix of results can also be viewed as a data frame, which is easier to
# read.
results.mat.df <- data.frame(results.mat)
View(results.mat.df)

# And the data frame can be saved to file once you set the file path.
write.csv(results.mat.df, "Your file path/Fpet Pairwise Cluster Fst.csv")


