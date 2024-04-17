# Check association of senescence score with p16 expression

# libraries
library(here)
library(AnnotationDbi)
library(ggplot2)

# First, download TCGA data on methylation and expression
# Folders: 1-functional-association/4-tcga-data/, scripts 1-3

directory <- "<your download directory>"

# Load data 
load("~/Dropbox/eca/sola/1-download/2-output/data.Rdata")
load("1-functional-association/4-tcga-data/2-output/data.Rdata")
normal <- data |>
  filter(type == "Control")

# Identify CDKN2A gene in ENSGDB
ensg <- AnnotationDbi::select(org.Hs.eg.db, keys = "CDKN2A",
                              columns = "ENSEMBL",
                              keytype = "SYMBOL")
normal$p16 <- NA

# get expression data
for (p in projects){
  if(file.exists(paste0(directory, "/TCGA/", p, "/expression"))){
    load(paste0(directory, "/TCGA/", p, "/expression/data_exp.Rdata"))
    expr <- assays(data)$tpm_unstrand
    ind <- match(substr(colnames(expr), 1, 16), normal$barcode2)
    normal[ind,]$p16 <- expr[grep(ensg$ENSEMBL, rownames(expr)),]
  }
}

normal <- normal[!is.na(normal$p16),] # 396 samples
hist(normal$p16)

# Compile beta for these samples
pb <- txtProgressBar(min = 1, max = length(unique(projects)), style = 3)

for (p in unique(normal$project)){
  setTxtProgressBar(pb, which(unique(normal$project)==p))
  load(paste0(directory, "/TCGA/", p, "/beta.Rdata"))
  
  tmp <- as.data.frame(beta[,normal[normal$project==p,]$barcode1])
  if(p == unique(normal$project)[1]){
    beta_merged <- tmp
  } else {
    intersect <- intersect(rownames(tmp), rownames(beta_merged))
    beta_merged <- beta_merged[intersect,]
    tmp <- tmp[intersect,]
    beta_merged <- cbind(beta_merged, tmp)
  }
  
  rm(tmp);gc()
}

setwd(here())
save(beta_merged, file = paste0(directory, "/TCGA/beta_expr.Rdata"))
save(normal, file = "1-functional-association/1-senescence/5-output/pheno_normal.Rdata")

# compute score
# load coefs
load("1-functional-association/1-senescence/3-output/coef.Rdata")
df$padj <- p.adjust(df$p, method = "fdr")
df <- df[df$padj<0.05 & !df$chr %in% c("chrX", "chrY"),]
df <- df[order(df$padj, decreasing = F),]

# double check matching of beta/pheno
intersect <- intersect(normal$barcode1, colnames(beta_merged))
beta <- beta_merged[,intersect]
pheno <- normal[match(intersect, normal$barcode1),]
identical(pheno$barcode1, colnames(beta))
rm(beta_merged, normal);gc()

# coerce to numeric
source("0-data/src/coerce_numeric.R")
beta <- coerce_numeric(beta)

# Score from df
hist(df$coef)
weights <- ifelse(df$coef>0, 1, -1)
names(weights) <- df$cg

intersect <- intersect(rownames(beta), names(weights))
b <- beta[intersect,]
w <- weights[intersect]

tmp <- apply(b, 2, function(i){
  return(sum(i*w)/length(w))
})

identical(names(tmp), pheno$barcode1)
pheno$idx <- tmp

t <- pheno |> 
  dplyr::group_by(tissue_simplified) |> 
  dplyr::count() |> 
  dplyr::filter(n > 3)

pheno |> 
  dplyr::filter(tissue_simplified %in% t$tissue_simplified) |> 
  ggplot(aes(x = idx,
             y = log2(p16+0.00001),
             colour = ic)) +
  geom_point() +
  geom_smooth(method = "lm",
              se = F) +
  ggpubr::stat_cor() +
  facet_wrap(~tissue_simplified) # positive association, although not significant for all tissues.

pheno |> 
  dplyr::filter(tissue_simplified %in% t$tissue_simplified) |> 
  ggplot(aes(x = idx,
             y = log2(p16+0.00001),
             colour = tissue_simplified)) +
  geom_point() +
  geom_smooth(method = "lm",
              se = F) +
  ggpubr::stat_cor() 

save(pheno, file = "1-functional-association/1-senescence/5-output/pheno.Rdata")