# Evaluate association of CpGs with senescence

# libraries
library(dplyr)
library(ggplot2)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

# load files
load("1-functional-association/1-senescence/2-qc/out/pheno.Rdata")
load("1-functional-association/1-senescence/2-qc/out/beta.Rdata")

# sanity check
identical(pheno$basename, colnames(beta))

# Create DF for association with senescence
df <- data.frame(matrix(nrow = nrow(beta),
                        ncol = 3))
colnames(df) <- c("cg", "coef", "p")

# Construct model, accounting for set + cell type; loop over each CpG
pb <- txtProgressBar(min = 1, max = nrow(beta), style = 3)
for (i in 1:nrow(beta)){
  setTxtProgressBar(pb, i)
  fit <- lm(beta[i,] ~ pheno$type + as.factor(pheno$set) + as.factor(pheno$celltype))
  df$cg[i] <- rownames(beta)[i]
  df$coef[i] <- summary(fit)$coefficients[2, 1]
  df$p[i] <- summary(fit)$coefficients[2, 4]
}

# Check histograms
hist(df$p)
hist(df$coef)
# hist(p.adjust(df$p, method = "fdr"))

# Add annotation
anno <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19, lociNames = df$cg)
identical(df$cg, anno$Name)
df$chr <- anno$chr
save(df, file = "1-functional-association/1-senescence/3-output/coef.Rdata")