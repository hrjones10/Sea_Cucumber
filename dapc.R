### LIBRARIES
library(adegenet)
library(ade4)
library(ggplot2)
library(reshape)
library(stringr)
library(stackr) #v.0.4.5 
#stackr download--> older package https://rdrr.io/github/thierrygosselin/stackr/
install.packages("remotes")
remotes::install_github("thierrygosselin/stackr")

strata <- read.table("strata.cucumbers.717ind.txt", header = T)

### Read in the VCF file as genind object
#can't seem to use function "vcfr2genind" ####
CUCUMBER_genind <- vcfR2genind(
  data = "filtered_2719neutral_snps.recode.vcf",
  strata = strata,
  imputation.method = NULL,
  parallel.core = 4)

gen.dat <- CUCUMBER_genind$genind.no.imputation

#####################################
### No prior - use find.clusters
####################################

### Find clusters using K-means
grp <- find.clusters(gen.dat, max.n.clust=25)

# -- when prompted select all PCs and K with the lowest BIC (2)

names(grp) 
individuals <- grp$grp 
data_kmeans <- data.frame(strata[,1],grp$grp)
names(data_kmeans)
colnames(data_kmeans) = c("IND","K")
data_kmeans
individuals

grp$size

table(pop(gen.dat), grp$grp)

# Perform DAPC with K = 2 clusters 
dapc_no_apriori <-dapc(gen.dat, grp$grp)

# -- when prompted select 200 PCs and only 1 eigenvalue

# number of PCs to retain:
temp <- optim.a.score(dapc_no_apriori, n.sim = 20)

# re-run DAPC
dapc_no_apriori <-dapc(gen.dat, grp$grp, n.pca = temp$best)
dapc_no_apriori <-dapc(gen.dat, grp$grp, n.pca = 20)

dapc_no_apriori$eig/sum(dapc_no_apriori$eig)*100

### Plot discriminant function 1
scatter(dapc_no_apriori,scree.da = TRUE, posi.da = "bottomleft", scree.pca = FALSE,  bg="white", legend=TRUE, solid=.4, cstar = 0, clabel = 0)



#######################################
####### With sampling location prior
#######################################

### Determine the alpha-score 
dapc2 <- dapc(gen.dat, n.pca=200) #200 for dataset with 717 inds
temp <- optim.a.score(dapc2, n.sim = 20)

### Run the DAPC with n.pca = temp$best
dapc2 <-dapc(gen.dat, gen.dat$pop, n.pca = temp$best)


### Plot discriminant functions 1 and 2
scatter(dapc2,1,2,  bg="white", scree.da=TRUE, scree.pca = FALSE, posi.da = "bottomright", legend=FALSE, solid=.4, cstar = 1, clabel = 0)
