#Author: Kevin Quinteros
#Date: Jan 26, 2022
#Analysis: SNPCLUSTER
#Purpose: determine if 218 individual is hybrid from two parental populations
################################set environment#################################
library(adegenet)
setwd("~/Projects/git_repositories/github/Fig_wasp_phylogeography/")
################################Load dataset####################################
#Unphased SNPs
load("Data/SNP_data/SNP_50p/Method1/gen_50_fixed.Rda")
pop(gen_snp_50) <- gen_snp_50$strata$fine_region 
gen_snp_50 <- gen_snp_50[order(as.numeric(as.character(gen_snp_50$strata$Lon)))]

#Phased SNPs
#load("Data/Phased_data/UCE_50p/Method1/Phased_snps/phased_snp_50.Rda")
#pop(phased_snp_50) <- phased_snp_50$strata$fine_region
#phased_snp_50 <- phased_snp_50[order(as.numeric(as.character(phased_snp_50$strata$Lon)))]

load("Data/Phased_data/UCE_50p/Method3/Phased_snps/phased_snp_50_fixed.Rda")
pop(phased_snp_50__unedit) <- phased_snp_50__unedit$strata$fine_region
phased_snp_50__unedit <- phased_snp_50__unedit[order(as.numeric(as.character(phased_snp_50__unedit$strata$Lon)))]
##########################SNPCLUSTER with unphased data#########################
#What's the optimal value of clusters (K) using AIC for model comparision
a.aic <- snapclust.choose.k(20, gen_snp_50)
#plot results 
pdf('output/snpcluster/unphased_snpcluster_optimal_K.pdf', width = 8, height = 5)
plot(a.aic, type = "b", cex = 2, xlab = "k", ylab = "AIC")
points(which.min(a.aic), min(a.aic), col = "blue", pch = 20, cex = 2)
abline(v = 6, lty = 2, col = "red")
dev.off()
#optimal K-cluster was K=2 
#run clsuter analysis with K=2 
a.clust <- snapclust(gen_snp_50, k = 2)
a.clust$group
#assign population information to individual snpcluster results 
a.tab <- table(gen_snp_50$pop, a.clust$group) #assign groups from our pop info
#population assignment 
a.tab
pdf('output/snpcluster/unphased_assignment_of_cluster_to_ind.pdf', width = 8, height = 5)
table.value(a.tab, col.labels = 1:2)
dev.off()
gen_snp_50$strata[87,]
#P340 (218) is an individual from Jalisco but assigned to the southern population
#P345 (219) is an individual from Morelos but is split between North and south 
#create structure plot

pdf("output/snpcluster/unphased_snpclusterK2.pdf", width = 8, height = 5)
compoplot(a.clust) 
dev.off()
######################DAPC############################################
#We ended up removing the individual P345-P347, since they had alot of missing data
gen_snp_50 <- gen_snp_50[-c(87:89),]
View(gen_snp_50$tab)

pdf("output/snpcluster/unphased_sdapc.pdf", width = 8, height = 5)
d.dapc <- dapc(gen_snp_50, n.pca = 20, n.da = 2)
scatter(d.dapc, clab = 0.85, col = funky(24),
             posi.da="topleft", posi.pca = "bottomleft", scree.pca = TRUE)
dev.off()

##########################SNAPCLUSTER Models########################
#run model, with no consideration for hybrids
res.no.hyb <- snapclust(gen_snp_50, k = 2, hybrids = FALSE)
pdf("output/snpcluster/unphased_snapclust_wo_hybrids.pdf", width = 8, height = 5)
compoplot(res.no.hyb, n.col = 2, col.pal = hybridpal(),
          main = "snapclust without hybrids",
          txt.leg = c("North", "South"))
dev.off()
#run model with consideration for hybrids but no backcrosses 
#To be formally identified,hybrids need to be modelled as a separate population, whose allele frequency distribution is somewhere between the parental populations.
res.hyb <- snapclust(gen_snp_50, k = 2, hybrids = TRUE)
pdf("output/snpcluster/unphased_snapclust_w_hybrids_no_backcrosses.pdf",
    width = 8, height = 5)
compoplot(res.hyb, n.col = 2, col.pal = hybridpal(), 
          main = "snapclust with hybrids but no backcrosses",
          txt.leg = c("North", "South", "0.5_A-0.5_B"))
dev.off()

res2.back <- snapclust(gen_snp_50, k=2, hybrids = TRUE, hybrid.coef = c(.25, .5))
pdf("output/snpcluster/unphased_snapclust_w_hybrids_w_backcrosses.pdf", 
    width = 8, height = 5)
compoplot(res2.back, col.pal = funky,
          main = "snapclust with F1 and backcrosses")
dev.off()

##########################run with the phased dataset##########################
#What's the optimal value of clusters (K) using AIC for model comparision
b.aic <- snapclust.choose.k(20, phased_snp_50__unedit)
#plot results 
pdf('output/snpcluster/phased_snpcluster_optimal_K.pdf', width = 8, height = 5)
plot(b.aic, type = "b", cex = 2, xlab = "k", ylab = "AIC")
points(which.min(b.aic), min(b.aic), col = "blue", pch = 20, cex = 2)
abline(v = 6, lty = 2, col = "red")
dev.off()
#optimal K-cluster was K=2 
#run clsuter analysis with K=2 
b.clust <- snapclust(phased_snp_50__unedit, k = 2)
b.clust$group
#assign population information to individual snpcluster results 
b.tab <- table(phased_snp_50__unedit$pop, b.clust$group) #assign groups from our pop info
#population assignment 
b.tab
pdf('output/snpcluster/phased_assignment_of_cluster_to_ind.pdf', width = 8, height = 5)
table.value(b.tab, col.labels = 1:2)
dev.off()
View(phased_snp_50__unedit$tab)
#P342 (218) is an individual from Jalisco but assigned to the southern population
#P342 is the hybrid. In a Tina voice "Dumbass" 
#P345 (219) is an individual from Morelos but is split between North and south 
#create structure plot
pdf("output/snpcluster/phased_snpclusterK2.pdf", width = 8, height = 5)
compoplot(b.clust) 
dev.off()

######################DAPC############################################
#We ended up removing the individual P344-P347, since they had alot of missing data
phased_snp_50__unedit <- phased_snp_50__unedit[-c(85:88),]

pdf("output/snpcluster/phased_sdapc.pdf", width = 8, height = 5)
d.dapc <- dapc(phased_snp_50__unedit, n.pca = 20, n.da = 2)
scatter(d.dapc, clab = 0.85, col = funky(24),
        posi.da="topleft", posi.pca = "bottomleft", scree.pca = TRUE)
dev.off()
############################SNAPCLUSTER##############################
#run model, with no consideration for hybrids
res.no.hyb <- snapclust(phased_snp_50__unedit, k = 2, hybrids = FALSE)
pdf("output/snpcluster/phased_snapclust_wo_hybrids.pdf", width = 20, height = 10)
compoplot(res.no.hyb, n.col = 2, col.pal = hybridpal(),
          main = "snapclust without hybrids",
          txt.leg = c("North", "South"), show.lab = T, cex.names =0.5)

dev.off()
?compoplot
#run model with consideration for hybrids but no backcrosses 
#To be formally identified,hybrids need to be modelled as a separate population, whose allele frequency distribution is somewhere between the parental populations.
res.hyb <- snapclust(phased_snp_50__unedit, k = 2, hybrids = TRUE)
pdf("output/snpcluster/phased_snapclust_w_hybrids_no_backcrosses.pdf",
    width = 8, height = 5)
compoplot(res.hyb, n.col = 2, col.pal = hybridpal(), 
          main = "snapclust with hybrids but no backcrosses",
          txt.leg = c("North", "South", "0.5_A-0.5_B"))
dev.off()

res2.back <- snapclust(phased_snp_50__unedit, k=2, hybrids = TRUE, hybrid.coef = c(.25, .5))
pdf("output/snpcluster/phased_snapclust_w_hybrids_w_backcrosses.pdf", 
    width = 8, height = 5)
compoplot(res2.back, col.pal = funky,
          main = "snapclust with F1 and backcrosses")
dev.off()
#############################run with STRUCTURE data############################
data <- read.structure("Data/STRUCTURE_FILES/WASP/pha_random_50p.str",
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
?compoplot
res2.back <- snapclust(data, k=2, hybrids = TRUE, hybrid.coef = c(.25, .5),parent.lab = unique(pop(data)))
pdf("output/snpcluster/structuredata_snapclust_w_hybrids_w_backcrosses.pdf", 
    width = 8, height = 5)
compoplot(res2.back, col.pal = funky,
          main = "snapclust with F1 and backcrosses")
dev.off()
compoplot

res3.back <- snapclust(data, k=2, hybrids = TRUE, hybrid.coef = c(.125,.25, .5))

compoplot(res3.back, col.pal = funky,
          main = "snapclust with F1 and backcrosses")
#P342 is the hybrid individual 
