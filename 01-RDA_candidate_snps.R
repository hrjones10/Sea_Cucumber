## This script performs a RDA and detects candidate SNPs associated with environmental variables 

####HJ seeing if she can edit this from Rstudio  

library(usdm)
library(vegan)

## SNP DATA 
snps <- read.table("3699snps_freqs.txt", header = T)
View(snps)# convert to matrix
snp.mat <- as.matrix(snps)
View(snp.mat)

# create vector of Site names
popnames <- rownames(snp.mat)
popnames 
#"OGD" "QUA" "HOP" "TBL" "CAL" "TOL" "PRI" "LEG" "JUA" "SEL" "REN" "SGI" "MAZ" "AK1" "AK2" "AK3" "AK4" "LAS" "JER" "TOF" "CRA" "RBY" "SHE" "MAL"

# Hellinger transformation of vector of abundances--> relative species composition rather than magnitude of abundance values... accounts for moderately skewed datasets
snp.hel <- decostand(snp.mat, method = "hellinger")
View(snp.hel)
## Environmental DATA

env1 <- read.table("env_predictors.txt", header = T)
rownames(env1) <- env1$Site
rownames(env1)
View(env1) #file looks good to me

# order rows to match snps that are correlated across environmental gradient
env1 <- env1[match(popnames, env1$Site), ] #need to capitalize "Site", be sure to check this on other scripts 
View(env1)
# separate lat/long and env vars
coords <- env1[,2:3]
env.vars <- env1[,4:ncol(env1)]

env.vars #okay good got data

# check variance inflation factor (VIF) of environmental variables
# VIF checks for for multicollinearity among variables
vif(env.vars)
#Variables        VIF
#1                  Mean_surface_salinity 489.123522
#2               Minimum_surface_salinity  41.865633
#3               Maximum_surface_salinity 285.882143
#4               Mean_surface_temperature  93.946589
#5            Minimum_surface_temperature  18.261385
#6            Maximum_surface_temperature  40.639257
#7  Mean_bottom_chlorophyll_concentration  26.620819
#8           Mean_bottom_current_velocity   5.438088
#9           Mean_bottom_dissolved_oxygen  11.829937
#10               Mean_bottom_temperature  79.323197 
#11            Maximum_bottom_temperature  39.718829
#12            Minimum_bottom_temperature  23.493223 
#13                  Mean_bottom_salinity   3.830173 
#14         Mean_surface_current_velocity   2.552170 

###make sure you are dropping VIF like mean surface temp####

# scale predictors
env.scale <- scale(env.vars, scale = T, center = T)

# run redundancy analysis (RDA) between hellinger transformed SNPs and all environmental variables (14 RDAs in total)
rda1 <- rda(snp.hel ~., data = as.data.frame(env.scale))
summary(rda1)
plot(rda1)

rda0 <- rda(snp.hel ~ 1, data = as.data.frame(env.scale))
plot(rda0)

summary(rda1) #summary output of all 6 RDAs
#Call: rda(formula = snp.hel ~ Mean_surface_salinity + Minimum_surface_salinity +      Maximum_surface_salinity + Mean_surface_temperature + Minimum_surface_temperature +      Maximum_surface_temperature + Mean_bottom_chlorophyll_concentration +      Mean_bottom_current_velocity + Mean_bottom_dissolved_oxygen +      Mean_bottom_temperature + Maximum_bottom_temperature + Minimum_bottom_temperature +      Mean_bottom_salinity + Mean_surface_current_velocity, data = as.data.frame(env.scale)) 

# considers partitioning of the variance, species (population here) scores and site scores (weighted sums of species scores)

RsquareAdj(rda1) #why take R squared--> correlation between the variables; percent of variance explained--> variation in y explained by variation in x
#$r.squared= 0.6493902 + $adj.r.squared- 0.1039971
anova(rda1) 
#999 permutations--> randomize the data to generate a null dist; ran model 1,000 times (randomly assign members to various groups)
#not assuming normality, just generating distribution
#is test statistic more significant in real data than the randomized dataset X999 runs

#readout:      Df   Variance    F     Pr(>F)  
#Model        14   0.064515  1.1907   0.015 *
#Residual     9    0.034832     

anova(rda1, by = "axis") 
#this considers all of the variation across each RDA (1-14) rather than one model with 14 Df

####################PART 2#

#### DETECT SNP OUTLIERS -->   "X84198_24= SNP"

# Function for detecting outliers (from Forester et al. 2018) #should we check this paper out?####
outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)                  
  # find +/- z sd from mean loading     
  x[x < lims[1] | x > lims[2]]   # locus names in these tails
}

# Outliers from the RDA with all env vars considered

## get the loadings
load.rda <- summary(rda1)$species[,1:3]
load.rda[,1]
load.rda

# for axes 1 and 2: 
cand1.3SD <- outliers(load.rda[,1], 3)
cand2.3SD <- outliers(load.rda[,2], 3)

ncand1.3SD <- length(cand1.3SD)
ncand2.3SD <- length(cand2.3SD)
ncand1.3SD #53  
ncand2.3SD #12

ncand <- ncand1.3SD+ncand2.3SD  
ncand #65

### 
cand1.3SD.df <- cbind.data.frame(rep(1, times = length(cand1.3SD)), names(cand1.3SD), unname(cand1.3SD)); colnames(cand1.3SD.df) <- c("axis", "snp", "loading")

cand2.3SD.df <- cbind.data.frame(rep(2, times = length(cand2.3SD)), names(cand2.3SD), unname(cand2.3SD)); colnames(cand2.3SD.df) <- c("axis", "snp", "loading")

cand <- rbind(cand1.3SD.df, cand2.3SD.df)
cand$snp <- as.character(cand$snp) #loading=how important an enviro axis is to SNP####

cand.mat <- matrix(nrow=(ncand), ncol=8)  # ncol = number of predictors/columns; 65 rows for cand SNPs
colnames(cand.mat) <- c("PCA1", "PCA2", "PCA3", "BO2_curvelmean_bdmean", "BO2_tempmean_bdmean", "BO2_tempmin_bdmean", "BO2_salinitymean_bdmean", "BO2_curvelmean_ss") #empty candidate matrix

#Joe's help####
length(cand.mat[i,]) #8
length(apply(env.scale,2,function(x) cor(x,snp.gen))) #14... should match... 
#also names of columns in env.scale file do not match cand.mat column names...

###ERROR:in cand.mat[i, ] <- apply(env.scale, 2, function(x) cor(x, snp.gen)) :  number of items to replace is not a multiple of replacement length####

apply(env.scale,2,function(x) cor(x,snp.gen))
#Mean_surface_salinity              Minimum_surface_salinity 
#0.18107736                            0.26385105 
#Maximum_surface_salinity              Mean_surface_temperature 
#0.12919691                           -0.17090813 
#Minimum_surface_temperature           Maximum_surface_temperature 
#-0.24348353                            0.15511505 
#Mean_bottom_chlorophyll_concentration          Mean_bottom_current_velocity 
#-0.28749631                            0.24724551 
#Mean_bottom_dissolved_oxygen               Mean_bottom_temperature 
#0.09651374                           -0.32276679 
#Maximum_bottom_temperature            Minimum_bottom_temperature 
#-0.14440751                           -0.15871449 
#Mean_bottom_salinity         Mean_surface_current_velocity 
#0.49810707                            0.10033776 

View(env.scale)
View(snp.mat)

#corr between enviro and allele freq @ specific SNP
#for (i in 1:length(cand$snp)) {
#  nam <- cand[i,2] #SNP name
#  snp.gen <- snp.mat[,nam]
#  cand.mat[i,] <- apply(env.scale,2,function(x) cor(x,snp.gen))
#} 

full.cand.df <-read.csv("cand.mat.csv")
View (full.cand.df)

#now that you have your new matrix, see correlations!
View(cand)

##SKIP?#
#full.cand.df <- cbind(cand, cand.mat)
#full.cand.df

cand$snp[duplicated(cand$snp)]  # check for duplicates
#no character repeats character(0)

# To determine which of the predictors each candidate SNP is most strongly correlated with:

for (i in 1:length(full.cand.df$snp)) {
  bar <- full.cand.df[i,]
  full.cand.df[i,12] <- names(which.max(abs(bar[4:11]))) # gives the variable
  full.cand.df[i,13] <- max(abs(bar[4:11]))   # gives the correlation
}
View(full.cand.df)

colnames(full.cand.df)[12] <- "predictor"
colnames(full.cand.df)[13] <- "correlation"

View(full.cand.df)

#how to plot the revised RDA?####


table(full.cand.df$predictor)  

#BO2_curvelmean_bdmean- 3 
#BO2_salinitymean_bdmean-3  
#BO2_tempmean_bdmean- 51
#BO2_tempmin_bdmean- 1
#PC2-1

table(full.cand.df$axis) 
# 1  2 
#51  8 

#creating the plot: ####
library(devtools)
install_github("vqv/ggbiplot")
sal.temp.pca <- prcomp(env1[,c(4:9)], center = TRUE,scale. = TRUE)
summary(sal.temp.pca)

#Importance of components:
#                         PC1    PC2    PC3     PC4     PC5     PC6
#Standard deviation     1.8945 1.1956 0.9132 0.30510 0.22632 0.05514
#Proportion of Variance 0.5982 0.2383 0.1390 0.01551 0.00854 0.00051
#Cumulative Proportion  0.5982 0.8364 0.9754 0.99096 0.99949 1.00000


library(ggbiplot)
ggbiplot(sal.temp.pca)


env2<-env1[,4:9]
bob<-prcomp(env2)
bob$x
bob$x[,1] #gets you values for PC1, etc.
bob<-prcomp(env2,scale=T,center=T)

bob$x[,1:3]
pc1to3 <-bob$x[,1:3]
View(pc1to3) #new file with PCA1, PCA2 and PCA3
d <- pc1to3
Site <- rownames(d)
rownames(d) <- NULL
data <- cbind(Site,d)
View(data) #put site column in
sub.env1 <- env1[,c(1,10:17)]
View(sub.env1) #enviro variables
total <- merge(data, sub.env1,by="Site") #merge enviro and PCAs
View(total) #yay dataset complete! :) 

