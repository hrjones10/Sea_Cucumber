## This script was used to compute polygenic scores across all candidate loci 

## 1. calculate allele frequencies of candidate loci using PLINK and organize data

allele_freq <- read.table("candidates.frq.strat", header = T, stringsAsFactors = F)
head(allele_freq)

temp_dat <- read.csv("meanbottomtemp_atsites.csv", header = T, stringsAsFactors = F)
head(temp_dat)

allele_freq$TEMP <- rep(temp_dat$TEMP, length(unique(allele_freq$SNP)))
allele_freq

allele_freq_temp <- select(allele_freq, c(SNP, MAF, TEMP))
head(allele_freq_temp)

## 2. perform correlations between allele frequencies at each site and environment (mean bottom temperature) per locus

correlations_by_snp <- allele_freq_temp %>%
  group_by(SNP) %>%
  summarize(COR=cor(TEMP, MAF))

correlations_by_snp


## 3. Convert vcf with 71 outlier SNPs to 012 format with vcftools then make .txt file and read in here 

snps <- read.table("candidates_012.txt", header = F, stringsAsFactors = F)
dim(snps)
snps <- snps[,-1]

# read in candidate and individual IDs

INDS <- read.table("individual_id.txt", header = T)
INDS

SNPS <- read.table("candidates_id.txt", header = T)
SNPS

colnames(snps) <- SNPS$ID
rownames(snps) <- INDS$IND


## 3. invert 0 and 2 for SNPs with a negative correlation 

head(correlations_by_snp)

snps[snps == -1] <- NA
cor_neg <- filter(correlations_by_snp, COR < 0)
cor_neg

snps_neg <- select(snps, cor_neg$SNP)
dim(snps_neg)  
for (i in 1:ncol(snps_neg)) {
  snps_neg[,i][which(snps_neg[,i] == 0)] <- 9 
  snps_neg[,i][which(snps_neg[,i] == 2)] <- 0
  snps_neg[,i][which(snps_neg[,i] == 9)] <- 2
}
snps_neg

cor_pos <- filter(correlations_by_snp, COR > 0)
snps_pos <- select(snps, cor_pos$SNP)
dim(snps_pos)

snps_recode <- cbind(snps_neg, snps_pos)
dim(snps_recode)
head(snps_recode)

## 4. calculate polygenic score (sum across rows) 

polygenic_score <- as.data.frame(rowSums(snps_recode, na.rm = T))
colnames(polygenic_score) <- "SCORE"
head(polygenic_score)
polygenic_score$IND <- rownames(polygenic_score)

polygenic_score$SITE <- substr(polygenic_score$IND, 1, 3)
head(polygenic_score)

site_temp <- temp_dat[, c(1,4)]
site_temp
colnames(site_temp) <- c("SITE", "TEMP")

polygenic_score <- merge(polygenic_score, site_temp, by = "SITE")
polygenic_score

## 5. Test relationship between polygenic score and temperature 

plot(polygenic_score$TEMP, polygenic_score$SCORE)

# linear:
mod1 <- lm(polygenic_score$SCORE ~ polygenic_score$TEMP)
summary(mod1)

# quadratic: 
mod2 <- lm(polygenic_score$SCORE ~ polygenic_score$TEMP + I(polygenic_score$TEMP^2))
summary(mod2)

AIC(mod1, mod2)
