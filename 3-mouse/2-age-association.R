# Age association, mouse
# Load mouse data
# Note we are using our own discovery set as GSE mouse methylation atlas uses multiple strains (129/Sv C57BL/6J    Mixed) and different ages that did not transfer well to our dataset. We use BALB/c mice that are from the same breed and housing as those used in our experiment.
load("~/Dropbox/data/mouse-project/beta_final.Rdata")
load("~/Dropbox/eca/mouse-project/2-uibk/1-preprocessing/1-output/pheno.Rdata")

intersect <- intersect(pheno$basename, colnames(beta_final))
pheno <- pheno[match(intersect, pheno$basename),]
beta_final <- beta_final[,intersect]

pheno <- droplevels(pheno[pheno$Experiment=='age' & !is.na(pheno$Experiment),])
table(pheno$Experiment)
beta <- beta_final[,pheno$basename]
beta <- na.omit(beta)

identical(pheno$basename, colnames(beta))
sum(is.na(beta))

# set age and tissue (for correction)
age <- as.numeric(pheno$age)
tissue <- droplevels(pheno$Tissue.type)

# Visualise age distribution
pheno |> 
  ggplot(aes(x = Tissue.type,
             y = age)) +
  geom_point() +
  coord_flip()

# compute correlations across all CpGs
cors <- data.frame(cor=rep(NA,nrow(beta)),
                   pval=rep(NA,nrow(beta)))
rownames(cors) <- rownames(beta)

n <- nrow(beta)
pb <- txtProgressBar(min = 0, max = n, style = 3)

for(i in 1:nrow(beta)){
  c <- lm(as.numeric(beta[i,]) ~ age + tissue)
  cors[i,1] <- summary(c)$coefficients[2,1]
  cors[i,2] <- summary(c)$coefficients[2,4]
  setTxtProgressBar(pb, i)
}

# save age correlations
save(cors, file = "3-mouse/2-output/age.Rdata")

# p value histogram
hist(cors$pval)
