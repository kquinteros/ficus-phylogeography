#Author: Kevin Quinteros
#Date: Jan 26, 2022
#Analysis: SNPCLUSTER
#Purpose: determine if 218 individual is hybrid from two parental populations
################################set environment#################################
library(adegenet)

################################Load dataset####################################
load("phased_snp_50_fixed.Rda")
pop(phased_snp_50__unedit) <- phased_snp_50__unedit$strata$fine_region
phased_snp_50__unedit <- phased_snp_50__unedit[order(as.numeric(as.character(phased_snp_50__unedit$strata$Lon)))]
#############################run with STRUCTURE data############################
data <- read.structure("pha_random_50p.str",
                       n.ind = 100,
                       n.loc = 1880,
                       onerowperind = F,
                       col.lab = 1,
                       col.others = 0,
                       row.marknames = 0, ask = F)
strata(data) <- phased_snp_50__unedit$strata
data$strata$Phrapl_region[data$strata$Phrapl_region == "OCCIDENTAL"] <- "NORTH"
pop(data) <- data$strata$Phrapl_region
data <- data[order(as.numeric(as.character(data$strata$Lon)))]
#We ended up removing the individual P345-P347, since they had a lot of missing data
data <- data[-c(86:88),]
data$strata[76,]


res.no.hyb <- snapclust(data, k = 2, hybrids = FALSE)
pdf("output/snpcluster/structuredata_snapclust_wo_hybrids.pdf", width = 8, height = 5)
compoplot(res.no.hyb, n.col = 2, col.pal = hybridpal(),
          main = "snapclust without hybrids",
          txt.leg = c("North", "South"))
dev.off()

#run model with consideration for hybrids but no backcrosses 
#To be formally identified,hybrids need to be modelled as a separate population, whose allele frequency distribution is somewhere between the parental populations.

res.hyb <- snapclust(data, k = 2, hybrids = TRUE, parent.lab = unique(pop(data)))
pdf("output/snpcluster/structuredata_snapclust_w_hybrids_no_backcrosses.pdf",
    width = 8, height = 5)
compoplot(res.hyb, n.col = 2, col.pal = hybridpal(), 
          main = "snapclust with hybrids but no backcrosses", 
          , show.lab = T, lab = data$strata$Label)
dev.off()

res2.back <- snapclust(data, k=2, hybrids = TRUE, hybrid.coef = c(.25, .5),parent.lab = unique(pop(data)))
pdf("output/snpcluster/structuredata_snapclust_w_hybrids_w_backcrosses.pdf", 
    width = 8, height = 5)
compoplot(res2.back, col.pal = funky,
          main = "snapclust with F1 and backcrosses")
dev.off()

res3.back <- snapclust(data, k=2, hybrids = TRUE, hybrid.coef = c(.125,.25, .5))

compoplot(res3.back, col.pal = funky,
          main = "snapclust with F1 and backcrosses")
#P342 is the hybrid individual 
