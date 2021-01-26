## This script performs a RDA and detects canddiate SNPs associated with environmental variables 

library(usdm)
library(vegan)

## SNP DATA 
snps <- read.table("3699snps_freqs.txt", header = T)

# convert to matrix
snp.mat <- as.matrix(snps)

# create vector of site names
popnames <- rownames(snp.mat)
popnames

# Hellinger transformation
snp.hel <- decostand(snp.mat, method = "hellinger")

## ENV DATA

env1 <- read.table("env_predictors.txt", header = T)
rownames(env1) <- env1$site
rownames(env1)

# order rows to match snps
env1 <- env1[match(popnames, env1$site), ]

# separate lat/long and env vars
coords <- env1[,2:3]
env.vars <- env1[,4:ncol(env1)]

# check variance inflation factor (VIF) of environmental variables
vif(env.vars)

# scale predictors
env.scale <- scale(env.vars, scale = T, center = T)

# run RDA between hellinger transformed SNPs and all environmental variables (9 in total)
rda1 <- rda(snp.hel ~., data = as.data.frame(env.scale))
rda0 <- rda(snp.hel ~ 1, data = as.data.frame(env.scale))

summary(rda1)
RsquareAdj(rda1) 
anova(rda1) 
anova(rda1, by = "axis") 


##################################################################

#### DETECT OUTLIERS   

# Function for detecting outliers (from Forester et al. 2018)
outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)                   # find +/- z sd from mean loading     
  x[x < lims[1] | x > lims[2]]   # locus names in these tails
}


# Outliers from the rda with all env vars

## get the loadings
load.rda <- summary(rda1)$species[,1:3]
load.rda[,1]
load.rda

# for axes 1 and 2: 
cand1.3SD <- outliers(load.rda[,1], 3)
cand2.3SD <- outliers(load.rda[,2], 3)

ncand1.3SD <- length(cand1.3SD)
ncand2.3SD <- length(cand2.3SD)
ncand1.3SD
ncand2.3SD

ncand <- ncand1.3SD+ncand2.3SD  
ncand

### 
cand1.3SD.df <- cbind.data.frame(rep(1, times = length(cand1.3SD)), names(cand1.3SD), unname(cand1.3SD)); colnames(cand1.3SD.df) <- c("axis", "snp", "loading")

cand2.3SD.df <- cbind.data.frame(rep(2, times = length(cand2.3SD)), names(cand2.3SD), unname(cand2.3SD)); colnames(cand2.3SD.df) <- c("axis", "snp", "loading")

cand <- rbind(cand1.3SD.df, cand2.3SD.df)
cand$snp <- as.character(cand$snp)

cand.mat <- matrix(nrow=(ncand), ncol=8)  # ncol = number of predictors
colnames(cand.mat) <- c("PCA1", "PCA2", "PCA3", "BO2_curvelmean_bdmean", "BO2_tempmean_bdmean", "BO2_tempmin_bdmean", "BO2_salinitymean_bdmean", "BO2_curvelmean_ss")


for (i in 1:length(cand$snp)) {
  nam <- cand[i,2]
  snp.gen <- snp.mat[,nam]
  cand.mat[i,] <- apply(env.scale,2,function(x) cor(x,snp.gen))
}


full.cand.df <- cbind(cand, cand.mat)
full.cand.df

cand$snp[duplicated(cand$snp)]  # check for duplicates

# To determine which of the predictors each candidate SNP is most strongly correlated with:

for (i in 1:length(full.cand.df$snp)) {
  bar <- full.cand.df[i,]
  full.cand.df[i,12] <- names(which.max(abs(bar[4:11]))) # gives the variable
  full.cand.df[i,13] <- max(abs(bar[4:11]))              # gives the correlation
}

colnames(full.cand.df)[12] <- "predictor"
colnames(full.cand.df)[13] <- "correlation"

full.cand.df

table(full.cand.df$predictor)  
table(full.cand.df$axis) 

