#load required libraries
library(tidyverse)
library(adegenet)
library(hierfstat)
library(magrittr)

#set working directory#######
setwd("~/Projects/git_repositories/github/Fig_wasp_phylogeography/")

#read in locality data and format#####
sites <- read_csv("Data/Locality_Fpet.csv")
sites %<>% distinct(Site, .keep_all = TRUE)
`%not_in%` = Negate(`%in%`)
#these sites only have one 
sites %<>% filter(Site %not_in% c("96", "214", "220", "228"))

#subsetting of data####
moax <- sites %>%
  filter(State%in% c("OAX", "MOR"))
No_moax <- sites %>%
  filter(State!= "OAX") %>%
  filter(State != "MOR")
Baja <- sites %>% filter(State == "BAJA" )
Nor <- sites %>%
  filter(State != "OAX") %>%
  filter(State != "MOR") %>%
  filter(State != "BAJA")

#Calculate geographic distance#########
Dgeo <- fields::rdist.earth(as.matrix(sites[4:3]),miles = F)
Dgeo.m <- fields::rdist.earth(as.matrix(moax[4:3]),miles = F)
Dgeo.nm <- fields::rdist.earth(as.matrix(No_moax[4:3]),miles = F)
Dgeo.b <- fields::rdist.earth(as.matrix(Baja[4:3]),miles = F)
Dgeo.n <- fields::rdist.earth(as.matrix(Nor[4:3]),miles = F)

###This section is used to convert genetic data into genetic distance matrix####
#read in genetic data
load("Data/Phased_data/UCE_50p/Method1/Phased_snps/phased_snp_50.Rda")

#convert to heifstat object
hf_p50 <- genind2hierfstat(phased_snp_50, pop = phased_snp_50$pop)
hf_p50 %<>% filter(pop%not_in% c("96", "214", "220", "228"))
#subset dataset
###########Calculate Genetic distance with Oaxaca and Morelos
hf_m <- hf_p50 %>% filter(pop %in% moax$Site)
hf_m <- droplevels(hf_m)
###########Calculate Genetic distance with No Oaxaca and Morelos
hf_nm <- hf_p50 %>% filter(pop %in% No_moax$Site)
hf_nm <- droplevels(hf_nm)
###########Calculate Genetic distance with baja
hf_b <- hf_p50 %>% filter(pop %in% Baja$Site)
hf_b <- droplevels(hf_b)
###########Calculate Genetic distance with Northern Inland Pops###########
hf_nl <- hf_p50 %>% filter(pop %in% Nor$Site)
hf_nl <- droplevels(hf_nl)

#calculating genetic distance using different methods
# see Takezaki and Nei (1996)
#Nei's genetic distance 
Dgen <- genet.dist(hf_p50, diploid = T, method = "Nei87" ) %>% as.matrix()
Dgen.m <- genet.dist(hf_m, diploid = T, method = "Nei87" ) %>% as.matrix()
Dgen.nm <- genet.dist(hf_nm, diploid = T, method = "Nei87" ) %>% as.matrix()
Dgen.b <- genet.dist(hf_b, diploid = T, method = "Nei87" ) %>% as.matrix()
Dgen.n <- genet.dist(hf_nl, diploid = T, method = "Nei87" ) %>% as.matrix()

Nei87 <- list(Dgen,Dgen.n, Dgen.m, Dgen.b, Dgen.nm)
saveRDS(Nei87, file = "Data/Gen_dist/Nei87.RDS")




#mantel test

#results for all data
result = mantel.randtest(as.dist(Dgen),
                         quasieuclid(as.dist(Dgeo)),
                         nrepet = 999)
#results for morelos and oaxaca
result.m = mantel.randtest(as.dist(Dgen.m),
                           quasieuclid(as.dist(Dgeo.m)),
                           nrepet = 999)
#results for all northern samples. Excludes morelos and oaxaca
result.nm = mantel.randtest(as.dist(Dgen.nm), 
                            quasieuclid(as.dist(Dgeo.nm)),
                            nrepet = 999)
#results for baja 
result.b  = mantel.randtest(as.dist(Dgen.b),
                            quasieuclid(as.dist(Dgeo.b)),
                            nrepet = 999)
#results for northen populations without morelos, baja, and oaxaca 
result.n = mantel.randtest(as.dist(Dgen.n),
                           quasieuclid(as.dist(Dgeo.n)), 
                           nrepet = 999)

#create dataframe
dmat <- cbind(as.vector(Dgeo), as.vector(Dgen)) %>% as.data.frame()
dmat <- dmat[rowSums(dmat[])>1,]
dmat.m <- cbind(as.vector(Dgeo.m), as.vector(Dgen.m)) %>% as.data.frame()
dmat.m <- dmat.m[rowSums(dmat.m[])>1,]
dmat.nm <- cbind(as.vector(Dgeo.nm), as.vector(Dgen.nm)) %>% as.data.frame()
dmat.nm <- dmat.nm[rowSums(dmat.nm[])>1,]
dmat.b <- cbind(as.vector(Dgeo.b), as.vector(Dgen.b)) %>% as.data.frame()
dmat.b <- dmat.b[rowSums(dmat.b[])>1,]
dmat.n <- cbind(as.vector(Dgeo.n), as.vector(Dgen.n))%>% as.data.frame()
dmat.n <- dmat.n[rowSums(dmat.n[])>1,]

#createa annotations
tmp <- paste('R =',
             round(result$obs, digits=3),
             '|', 'p =',
             round(result$pvalue, 
                   digits=3), sep=' ')
tmp.m <- paste('R =',
               round(result.m$obs, digits=3),
               '|', 'p =',
               round(result.m$pvalue, 
                     digits=3), sep=' ')
tmp.nm <- paste('R =',
                round(result.nm$obs, digits=3),
                '|', 'p =',
                round(result.nm$pvalue, 
                      digits=3), sep=' ')
tmp.b <- paste('R =',
               round(result.b$obs, digits=3),
               '|', 'p =',
               round(result.b$pvalue, 
                     digits=3), sep=' ')
tmp.n <- paste('R =',
               round(result.n$obs, digits=3),
               '|', 'p =',
               round(result.n$pvalue, 
                     digits=3), sep=' ')


p <- ggplot(dmat, aes(V1,V2) ) + 
  geom_point() + 
  ylab("Nei's distance (1987)") +
  ggtitle(paste("Isolation by Distance", tmp, sep = " "))  +
  xlab("Geographic Distance") +
  theme_classic() 

ggsave( filename = "output/IBD/Nei_1987_p50.pdf",
        plot = p,
        device = "pdf",
        dpi = 300)

p.m <- ggplot(dmat.m, aes(V1,V2) ) + 
  geom_point() + 
  ylab("Nei's distance (1987)") +
  ggtitle(paste("Isolation by Distance", tmp.m, sep = " "))  +
  xlab("Geographic Distance") +
  theme_classic() 
ggsave( filename = "output/IBD/Nei_1987_p50_moax.pdf",
        plot = p.m,
        device = "pdf",
        dpi = 300)

p.nm <- ggplot(dmat.nm, aes(V1,V2) ) + 
  geom_point() + 
  ylab("Nei's distance (1987)") +
  ggtitle(paste("Isolation by Distance", tmp.nm, sep = " "))  +
  xlab("Geographic Distance") +
  theme_classic() 
ggsave( filename = "output/IBD/Nei_1987_p50_no-moax.pdf",
        plot = p.nm,
        device = "pdf",
        dpi = 300)

p.b <- ggplot(dmat.b , aes(V1,V2) ) + 
  geom_point() + 
  ylab("Nei's distance (1987)") +
  ggtitle(paste("Isolation by Distance", tmp.b, sep = " "))  +
  xlab("Geographic Distance") +
  theme_classic() 
ggsave( filename = "output/IBD/Nei_1987_p50_baja.pdf",
        plot = p.b,
        device = "pdf",
        dpi = 300)

p.n <- ggplot(dmat.n , aes(V1,V2) ) + 
  geom_point() + 
  ylab("Nei's distance (1987)") +
  ggtitle(paste("Isolation by Distance", tmp.n, sep = " "))  +
  xlab("Geographic Distance") +
  theme_classic() 
ggsave( filename = "output/IBD/Nei_1987_p50_nor.pdf",
        plot = p.n,
        device = "pdf",
        dpi = 300)

