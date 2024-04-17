# Check association of senescence score with ki67 expression

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

# Identify MKi67 (Ki67) 
ensg <- AnnotationDbi::select(org.Hs.eg.db, keys = "MKI67",
                              columns = "ENSEMBL",
                              keytype = "SYMBOL")
normal$ki67 <- NA

# get expression data
for (p in projects){
  if(file.exists(paste0(directory, "/TCGA/", p, "/expression"))){
    load(paste0(directory, "/TCGA/", p, "/expression/data_exp.Rdata"))
    expr <- assays(data)$tpm_unstrand
    ind <- match(substr(colnames(expr), 1, 16), normal$barcode2)
    normal[ind,]$ki67 <- expr[grep(ensg$ENSEMBL, rownames(expr)),]
  }
}

normal <- normal[!is.na(normal$ki67),] # 396 samples
hist(normal$ki67)
save(normal, file = "1-functional-association/2-proliferation/5-output/pheno_normal.Rdata")

# Beta is same as for p16
load(paste0(directory, "/TCGA/beta_expr.Rdata"))
identical(colnames(beta_merged), normal$barcode1)
normal$barcode1

# compute score
# load coefs
load("1-functional-association/2-proliferation/3-output/coef.Rdata")
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

# Score from df (inverted)
hist(df$coef)
weights <- ifelse(df$coef>0, 1, -1)
names(weights) <- df$cg

intersect <- intersect(rownames(beta), names(weights))
b <- beta[intersect,]
w <- weights[intersect]

tmp <- apply(b, 2, function(i){
  return(abs(1-sum(i*w)/length(w)))
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
             y = log2(ki67+0.00001),
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

save(pheno, file = "1-functional-association/2-proliferation/5-output/pheno.Rdata")
     