# Age correlation within epithelial fraction [cervical + buccal, inferred ic < 0.2]

# libraries
if(!require("ppcor")){
  install.packages("ppcor")
}
library(ppcor)

# Datasets: EGAD00010002079 and EGAD00010002080 (Controls only)
load("0-data/dat.Rdata")
dat <- dat |> 
  dplyr::filter(tissue == "epithelial")

# load beta matrices (after QC)
load("~/Dropbox/data/3c/beta_merged.Rdata")
beta_3c <- beta_merged
load("~/Dropbox/data/3c-buccal/beta_merged.Rdata")
beta_bu <- beta_merged
intersect <- intersect(rownames(beta_3c), rownames(beta_bu))
beta <- cbind(beta_3c[intersect,], beta_bu[intersect,])
beta <- beta[,dat$basename]
identical(dat$basename, colnames(beta))
rm(beta_merged, beta_3c, beta_bu);gc()

source("0-data/src/coerce_numeric.R")
beta <- coerce_numeric(beta)

# set up correlation test
corr.values <- data.frame(matrix(nrow = nrow(beta),
                                 ncol = 2))
rownames(corr.values) <- rownames(beta)
colnames(corr.values) <- c("corr", "p")

n <- nrow(beta)
pb <- txtProgressBar(min = 0, max = n, style = 3)

for (i in 1:n){
  tmp <- cor.test(dat$age, beta[i,], method = "spearman")
  corr.values$corr[i] <- tmp$estimate
  corr.values$p[i] <- tmp$p.value
  setTxtProgressBar(pb, i)
}
close(pb)
corr.values$p.adjust <- p.adjust(corr.values$p, n = nrow(corr.values),
                                 method = "fdr")
sum(corr.values$p.adjust<0.05 & abs(corr.values$corr)>0.4) # 3926 correlated.

# visualisation
hist(corr.values$p)
hist(corr.values$corr)
hist(corr.values$p.adjust)
save(corr.values, file = "2-age-association/2-output/age_epi.Rdata")
