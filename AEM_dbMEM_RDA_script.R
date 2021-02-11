## LIBRARIES

library(vegan)
library(reshape2)
library(dplyr)
library(SoDA) #no package?
library(adespatial) #no package?
library(ade4)
library(spdep)


#### SITE GEOGRAPHIC COORDINATES

sites <- read.table("24sites_lat_long.txt", header=TRUE)
sites <- sites[-c(21:24), ]  # remove Alaska sites
geo <- sites[, c(2:3)]
geo

geo$CAR <- geoXY(geo$LAT, geo$LONG) # Tranform geographic coordinates into cartesian coordinates; order is Lat, Long
geo

#### GENETIC DATA 

freq <- read.table("2719neutral_snps_noAK.frq.strat", header = T)
data2 <- data[,2:3]
data3 <- cbind(data2, data$MAF)
colnames(data3) <- c("SNP", "POP", "MAF")
head(data3)

### Convert from long to wide
geno <- dcast(data3, POP~SNP)

# Hellinger transformed data 
geno.h <- decostand(geno, "hellinger")

# pca 
geno.pca <- prcomp(geno.h, scale = T)

summary(geno.pca)
plot(geno.pca)
screeplot(geno.pca, npcs = 20, type = "lines", bstick = F)

# select 7 PC axes that explain > 5% of the total genetic variation
geno.pca.axes <- geno.pca$x[,1:7] 
head(geno.pca.axes)


#### (1a) dbMEM VARIABLES (straight-line)

euclidian_distances <- dist(geo$CAR, method="euclidean")
dbMEM <- dbmem(euclidian_distances, MEM.autocor = "non-null", store.listw = TRUE)
names(dbMEM)
dbMEM.vectors <- as.data.frame(dbMEM)
dbMEM.vectors

write.table(dbMEM.vectors, file = "dbMEM_straight_all_vectors.txt", col.names = T, row.names = F, quote = F, sep = "\t")


#### (1b) dbMEM VARIABLES (in-water)

in.sea <- read.table(file = "insea_matrix.txt", header = T)
rownames(in.sea) <- in.sea[,1]
in.sea <- in.sea[,-1]
in.sea <- in.sea[-c(21:24), -c(21:24)] # remove Alaska sites

names <- as.vector(sites$SITE)

in.sea.mat <- as.matrix(in.sea)
in.sea.mat <- in.sea.mat[names,names]
in.sea.dist <- as.dist(in.sea.mat)

dbMEM.inwater <- dbmem(in.sea.dist, MEM.autocor = 'non-null', store.listw = TRUE)
dbMEM.vectors.inwater <- as.data.frame(dbMEM.inwater)
dbMEM.vectors.inwater

write.table(dbMEM.vectors.inwater, file = "dbMEM_inwater_all_vectors.txt", col.names = T, row.names = F, quote = F, sep = "\t")

#### (2) AEM VARIABLES

# connectivity matrix: 
con_mat <- read.table("con_mat_mean_avg10years_pld120.txt", header = TRUE) # connectivity matrix based on dispersal/connectivity probabilities 
head(con_mat)
dim(con_mat)

con_mat <- as.matrix(con_mat)
head(con_mat) 
con_mat <- con_mat[as.vector(sites$Abb), as.vector(sites$Abb)]
#View(con_mat)

# highest probability between two sites: 
upper <- con_mat[upper.tri(con_mat, diag = TRUE)]
lower <- con_mat[lower.tri(con_mat, diag = TRUE)]

highest <- vector()

for (i in 1 : length(upper)) {
  if (upper[i] < lower[i]) {
    highest[i] <- lower[i]
  } else {
    highest[i] <- upper[i]
  }
}


# Replace the uppper matrix of con_mat with the highest values; lower matrix with 0s
con_mat2 <- con_mat
con_mat2[upper.tri(con_mat2, diag = TRUE)] <- highest

con_mat2[upper.tri(con_mat2, diag = TRUE)] == highest
head(con_mat2)

con_mat2[lower.tri(con_mat2, diag = T)] <- 0
diag(con_mat2) <- diag(con_mat)
dim(con_mat2)

## coords object: 
xy <- sites[,4:6]
xy$Nb <- seq(1, 20, 1)
xy <- xy[,c(4,3,2)]  # order is long, lat
xy

## links object: 
con_mat_ex <- expand.grid(con_mat2)
con_mat_ex$POP1 <- rep(colnames(con_mat2), each = 20)
con_mat_ex$POP2 <- rep(colnames(con_mat2), 20)
con_mat_ex$POP1id <- rep(seq(1,20,1), each = 20)
con_mat_ex$POP2id <- rep(seq(1,20,1), 20)

con_mat_ex$POP_id <- paste(con_mat_ex$POP1, con_mat_ex$POP2, sep = "-")
con_mat_ex <- con_mat_ex[,-c(2,3)]
con_mat_ex
colnames(con_mat_ex) <- c("PROB", "POP1", "POP2", "PAIR")

con.mat.highest.prob.no0 <- filter(con_mat_ex, PROB > 0)
nrow(con.mat.highest.prob.no0)
head(con.mat.highest.prob.no0)

edges.highest.prob <- con.mat.highest.prob.no0[,c(2,3)]
edges.highest.prob
class(edges.highest.prob)

edges <- edges.highest.prob
edges
edges <- filter(edges, POP1 != POP2)
nrow(edges)

# create binary matrix:
bin.mat <- aem.build.binary(coords=xy,link=edges, plot.connexions = TRUE)
str(bin.mat)
bin.mat

# edge weights = dispersal probabilities:
weights <- filter(con.mat.highest.prob.no0, POP1 != POP2)
weight.vec <- as.vector(weights[,1])
length(weight.vec)

# calculate AEM variables: 
cuke.aem.wt <- adespatial::aem(aem.build.binary = bin.mat, weight = weight.vec, rm.link0 = TRUE)
cuke.aem.nowt <- adespatial::aem(aem.build.binary = bin.mat, rm.link0 = TRUE) # no weight

AEM.vectors.wt <- as.data.frame(cuke.aem.wt$vectors)
AEM.vectors.nowt <- as.data.frame(cuke.aem.nowt$vectors)

write.table(AEM.vectors.wt, file = "AEM_weight_all_vectors.txt", sep = "\t", quote = F, col.names = T, row.names = F)

###############################
###############################
########### RDAs ##############
###############################
###############################


#### (1a) dbMEM (straight)

dbmem.mod0 <- rda(geno.pca.axes ~1, dbMEM.vectors)
dbmem.mod1 <- rda(geno.pca.axes ~., dbMEM.vectors)
summary(dbmem.mod1)

dbmem.all.ord <- ordistep(dbmem.mod0, scope = formula(dbmem.mod1), direction = "forward", permutations = 999)
dbmem.all.ord

dbmem.sub <- dbMEM.vectors[,c(1,2)]

dbmem.mod.sel <- rda(geno.pca.axes ~., as.data.frame(dbmem.sub))
summary(dbmem.mod.sel)


RsquareAdj(dbmem.mod.sel)
anova(dbmem.mod.sel, permutations = 999)
anova(dbmem.mod.sel, by = "margin", permutations = 999)
anova(dbmem.mod.sel, by = "axis", permutations = 999)


#### (1b) dbMEM (in-water)

dbmem.water.mod0 <- rda(geno.pca.axes ~1, dbMEM.vectors.inwater)
dbmem.water.mod1 <- rda(geno.pca.axes ~., dbMEM.vectors.inwater)
summary(dbmem.water.mod1)

dbmem.water.all.ord <- ordistep(dbmem.water.mod0, scope = formula(dbmem.water.mod1), direction = "forward", permutations = 999)
dbmem.water.all.ord

dbmem.water.sub <- dbMEM.vectors.inwater[,c(2,19)] # the dbMEM vectors selected by ordistep

dbmem.water.mod.sel <- rda(geno.pca.axes ~., as.data.frame(dbmem.water.sub))
summary(dbmem.water.mod.sel)

RsquareAdj(dbmem.water.mod.sel)
anova(dbmem.water.mod.sel, permutations = 999)
anova(dbmem.water.mod.sel, by = "margin", permutations = 999)
anova(dbmem.water.mod.sel, by = "axis", permutations = 999)


#### (2) AEM 

aem.mod0 <- rda(geno.pca.axes ~1, AEM.vectors.wt)
aem.mod1 <- rda(geno.pca.axes ~., AEM.vectors.wt)
summary(aem.mod1)

aem.sel.ord <- ordistep(aem.mod0, scope = formula(aem.mod1), direction = "forward", permutations = 999)
aem.sel.ord


aem.nowt.mod0 <- rda(geno.pca.axes ~1, AEM.vectors.nowt)
aem.nowt.mod1 <- rda(geno.pca.axes ~., AEM.vectors.nowt)
summary(aem.mod1)

aem.nowt.sel.ord <- ordistep(aem.nowt.mod0, scope = formula(aem.nowt.mod1), direction = "forward", permutations = 999)
aem.nowt.sel.ord

aem.sub <- AEM.vectors.wt[,c(1,2,4,11,8,17)]
aem.nowt.sub <- AEM.vectors.nowt[,c(1,3,5)]

mod.sel <- rda(geno.pca.axes ~., as.data.frame(aem.sub))
summary(mod.sel)
anova(mod.sel, permuations = 999)  
anova(mod.sel, by = "margin", permutations = 999) 
anova(mod.sel, by = "axis", permutations = 999)
RsquareAdj(mod.sel) 


mod.nowt.sel <- rda(geno.pca.axes ~., ad.data.frame(aem.nowt.sub))
summary(mod.nowt.sel)
anova(mod.nowt.sel, permutations = 999)
anova(mod.nowt.sel, by = "margin", permutations = 999)
anova(mod.nowt.sel, by = "axis", permutations = 999)
RsquareAdj(mod.nowt.sel)



### (3) FULL MODEL: AEM.SUB AND dbMEM.SUB
full.pred <- cbind(dbmem.sub, aem.sub, dbmem.water.sub)
full.pred
colnames(full.pred) <- c("dbMEM1s", "dbMEM2s", "AEM1", "AEM2", "AEM4", "AEM11", "AEM8", "AEM17", "dbMEM2w", "dbMEM19w")
full.rda <- rda(geno.pca.axes ~., full.pred)

summary(full.rda)
RsquareAdj(full.rda)
anova(full.rda, permutations = 999)
anova(full.rda, by = "margin", permutations = 999)
anova(full.rda, by = "axis", permutations = 999)


### (4) PARTIAL RDA 

prda1 <- rda(geno.pca.axes ~ AEM1 + AEM2 + AEM4 + AEM11 + AEM8 + AEM17 + Condition(dbMEM1s + dbMEM2s + dbMEM2w + dbMEM19w), data = full.pred)

anova(prda1, permutations = 999)
anova(prda1, by = "margin", permutations = 999)
summary(prda1)
RsquareAdj(prda1)



### (5) Variance Partitioning 
vp1 <- varpart(geno.pca.axes, aem.sub, dbmem.water.sub)
vp1
aFrac <- rda(geno.pca.axes, aem.sub, dbmem.sub, dbmem.water.sub)
anova(aFrac)

showvarparts(3, bg=2:4)
plot(vp1, bg = 2:4)

