## This script was used to compute polygenic scores across all candidate loci 

## 1. calculate allele frequencies of candidate loci using PLINK and organize data

allele_freq <- read.table("candidates.frq.strat", header = T, stringsAsFactors = F)
head(allele_freq)
#read out
#  CHR  SNP     CLST A1 A2 MAF MAC NCHROBS
#1   1 18276_19    1  C  T   0   0      62
#2   1 18276_19    2  C  T   0   0      62
#3   1 18276_19    3  C  T   0   0      56
#4   1 18276_19    4  C  T   0   0      62
#5   1 18276_19    5  C  T   0   0      40
#6   1 18276_19    6  C  T   0   0      56

temp_dat <- read.csv("meanbottomtemp_atsites.csv", header = T, stringsAsFactors = F)
head(temp_dat)
#minor allele frequency=0 across temp data
#     SNP     CLST MAF  TEMP   
#1 18502_23  OGD   0  9.22472
#2 18502_23  SGI   0 10.65276
#3 18502_23  LAS   0 11.29388
#4 18502_23  JER   0  9.95314
#5 18502_23  TOF   0 10.83483
#6 18502_23  CRA   0 10.90472

allele_freq$TEMP <- rep(temp_dat$TEMP, length(unique(allele_freq$SNP)))
allele_freq 

allele_freq_temp <- select(allele_freq, c(SNP, MAF, TEMP))
head(allele_freq_temp)
#   SNP      MAF  TEMP
#1 18276_19   0  9.22472
#2 18276_19   0 10.65276
#3 18276_19   0 11.29388
#4 18276_19   0  9.95314
#5 18276_19   0 10.83483
#6 18276_19   0 10.90472

## 2. perform correlations between allele frequencies at each site and environment (mean bottom temperature) per locus

correlations_by_snp <- allele_freq_temp %>%
  group_by(SNP) %>%
  summarize(COR=cor(TEMP, MAF))

correlations_by_snp
head(correlations_by_snp)
#list all SNPs
#    SNP      COR
#1 101_21   -0.323  (-)
#2 10406_30  0.541  (+)
#3 13062_15 -0.0887 (-)
#4 15144_18  0.548 (+)
#5 18276_19 -0.603 (-)
#6 18276_22 -0.622 (-)

## 3. Convert variant call format file (vcf) with 71 outlier SNPs to 012 format with vcftools then make .txt file and read in here 

snps <- read.table("candidates_012.txt", header = F, stringsAsFactors = F) ##dont have this file####
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
