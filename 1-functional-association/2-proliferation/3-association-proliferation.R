# Check asociation with proliferation

# libraries
library(ggplot2)
library(dplyr)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)

# data
load("1-functional-association/2-proliferation/1-output/pheno.Rdata")
directory <- "<your download directory/"
directory <- "~/Documents/Work/data.nosync/GEO/"
# load beta
load(paste0(directory, "proliferation/beta_merged.Rdata"))

# define sub-experiments: Mitomycin C and Serum withdrawal
p1 <- pheno |> 
  dplyr::filter(subexperiment %in% c("MMC") & culture_agent != 'NA')
table(p1$subculture)
hist(as.numeric(p1$days_in_culture))

p1 |> ggplot(aes(x = tissue,
                 y = as.numeric(days_in_culture),
                 colour = culture_agent,
                 shape = subculture))+
  geom_jitter() +
  facet_wrap(~coriell_id)

p1 |> ggplot(aes(x = tissue,
                 y = as.numeric(days_in_culture),
                 colour = culture_agent,
                 shape = batch_date))+
  geom_jitter() +
  facet_wrap(~coriell_id)

p2 <- pheno  |> 
  dplyr::filter(subexperiment %in% c("Serum withdrawal") & as.numeric(days_in_culture) > 30 & culture_pctfbs %in% c("0.50%",
                                                                                                                    "15%"))

pheno <- rbind(p1, p2) |> 
  dplyr::mutate(type = case_when(subexperiment == "MMC" & culture_agent == "MMC" ~ "non-proliferating",
                                 subexperiment == "Serum withdrawal" & culture_pctfbs == "0.50%" ~ "non-proliferating",
                                 TRUE ~ "Control"))
save(pheno, file = "1-functional-association/2-proliferation/3-output/pheno_filtered.Rdata")
beta <- beta_merged[,match(pheno$basename, colnames(beta_merged))]
identical(colnames(beta), pheno$basename)
rm(beta_merged);gc() # keep only small subset

# Delta-Betas
df <- data.frame(matrix(nrow = nrow(beta),
                        ncol = 3))
colnames(df) <- c("cg", "coef", "p")

# Construct model and loop; adjusting for subexperiment and tissue
pb <- txtProgressBar(min = 1, max = nrow(beta), style = 3)
for (i in 1:nrow(beta)){
  setTxtProgressBar(pb, i)
  fit <- lm(beta[i,] ~ pheno$type + as.factor(pheno$subexperiment) + as.factor(pheno$tissue))
  df$cg[i] <- rownames(beta)[i]
  df$coef[i] <- summary(fit)$coefficients[2, 1]
  df$p[i] <- summary(fit)$coefficients[2, 4]
}

# Add annotation
anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19, lociNames = df$cg)
identical(df$cg, anno$Name)
df$chr <- anno$chr
save(df, file = "1-functional-association/2-proliferation/3-output/coef.Rdata")

hist(df$p)
hist(p.adjust(df$p, method = "fdr"))
