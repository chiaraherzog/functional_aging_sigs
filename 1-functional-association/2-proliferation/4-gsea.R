# Exploratory plots and GSEA for Senescence-associated DMPs------------------

# libraries
library(ggplot2)
library(dplyr)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(methylGSA)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(stringr)

# load data
load("1-functional-association/2-proliferation/3-output/coef.Rdata")
df$padj <- p.adjust(df$p, method = "fdr")
df <- df[df$padj<0.05 & !df$chr %in% c("chrX", "chrY"),] # filter by pval and chr (no sex chromosomes)
df <- df[order(df$padj, decreasing = F),]
head(df)

# pheno
load("1-functional-association/2-proliferation/3-output/pheno_filtered.Rdata")
directory <- "<your download directory>"
directory <- "~/Documents/Work/data.nosync/GEO/"
# load beta
load(paste0(directory, "proliferation/beta_merged.Rdata"))
beta <- beta_merged[,pheno$basename]
identical(pheno$basename, colnames(beta))
rm(beta_merged);gc()

# exploratory plots:
# plot top 10
t <- beta[match(df$cg[1:10], rownames(beta)),]
tmp <- cbind(pheno, t(t))
tmp2 <- tmp |> 
  tidyr::pivot_longer(any_of(rownames(t)),
                      names_to = "cg",
                      values_to = "beta")

tmp2 |> 
  ggplot(aes(x = type,
             y = beta)) +
  geom_boxplot() + 
  geom_point(aes(colour = subexperiment)) +
  facet_wrap( ~ cg,
              scales = "free") # For most CpGs: with increasing methylation, decreasing proliferation (levels are lower in proliferating controls -> lower values indicate higher proliferation and should be considered for index.)

# plot a score
intersect <- intersect(rownames(beta), df$cg)
b <- beta[intersect,]
w <- df[match(intersect, df$cg),]
w <- ifelse(w$coef>0, 1, -1)
names(w) <- df$cg
tmp <- apply(b, 2, function(i){
  return(sum(i*w)/length(w))
})

identical(names(tmp), pheno$basename)
pheno$idx <- tmp

pheno |> 
  ggplot(aes(x = type,
             y = idx)) +
  geom_boxplot() + 
  geom_point(aes(colour = subexperiment)) # higher value indicates lower proliferation - index should be inverted

tmp <- apply(b, 2, function(i){
  return(abs(1-sum(i*w)/length(w)))
})
pheno$idx_inverted <- tmp

pheno |> 
  ggplot(aes(x = type,
             y = idx_inverted)) +
  geom_boxplot() + 
  geom_point(aes(colour = subexperiment)) # 'correct' directionality
save(pheno, file = "1-functional-association/2-proliferation/4-output/pheno_association.Rdata")


# Run GSEA for those that are significant
cpg.pval <- df$p
names(cpg.pval) <- df$cg

# GO with methylGSA
res1 = methylglm(cpg.pval = cpg.pval, minsize = 50, 
                 maxsize = 500, GS.type = "GO",array.type = "EPIC")
head(res1, 15)
hist(res1$padj)
x <- res1[res1$padj<0.05,]
x$Description # keratinization

# KEGG with methylGSA
res2 = methylglm(cpg.pval = cpg.pval, minsize = 50, 
                 maxsize = 500, GS.type = "KEGG",array.type = "450K")
x <- res2[res2$padj<0.05,]
x$Description # Olfactory transduction?
save(res1, res2, file = "1-functional-association/2-proliferation/4-output/gsa.Rdata")


# Manual: clusterProfiler
anno <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19, lociNames = df$cg)
cpgs <- anno |> 
  as.data.frame()

# 2. Patchway enrichment
genes <- unique(as.character(stringr::str_split(cpgs$UCSC_RefGene_Name, ";", simplify = T)))
genes <- genes[genes !=""]

allgenes <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19,
                          lociNames = rownames(beta))
allgenes <- unique(as.character(stringr::str_split(allgenes$UCSC_RefGene_Name, ";", simplify = T)))
allgenes <- allgenes[genes !=""]

enrichbp <- clusterProfiler::enrichGO(gene = genes,
                                      universe = allgenes,
                                      OrgDb = org.Hs.eg.db,
                                      keyType = 'SYMBOL',
                                      ont = "BP",
                                      pAdjustMethod = "BH")
# z <- as.data.frame(enrichbp)
# enrichplot::cnetplot(enrichbp, showCategory = 2)
# enrichplot::goplot(enrichbp, showCategory = 6)
# enrichplot::dotplot(enrichbp)

enrichmf <- clusterProfiler::enrichGO(gene = genes,
                                      universe = allgenes,
                                      OrgDb = org.Hs.eg.db,
                                      keyType = 'SYMBOL',
                                      ont = "MF",
                                      pAdjustMethod = "BH")
# z <- as.data.frame(enrichmf)
# enrichplot::cnetplot(enrichmf, showCategory = 20)
# enrichplot::goplot(enrichmf, showCategory = 20)

enrichcc <- clusterProfiler::enrichGO(gene = genes,
                                      universe = allgenes,
                                      OrgDb = org.Hs.eg.db,
                                      keyType = 'SYMBOL',
                                      ont = "CC",
                                      pAdjustMethod = "BH")
# z <- as.data.frame(enrichcc)

# Entrez
cpggenes <- bitr(genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
universegenes <- bitr(allgenes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

enrichpath <- enrichPathway(cpggenes$ENTREZID,
                            organism = "human",
                            universe = universegenes$ENTREZID,
                            pvalueCutoff = 0.05,
                            pAdjustMethod = "BH", 
                            readable = T,
                            minGSSize = 5)
x <- as.data.frame(enrichpath)

save(enrichbp, enrichmf, enrichcc,
     enrichpath, file = "1-functional-association/2-proliferation/4-output/clusterProfiler.Rdata")
rm(list=ls())

# Check association with ki67 expression? ----------------------

load("~/Dropbox/eca/sola/6-functional/1-proliferation/3-output/coef.Rdata")
df$padj <- p.adjust(df$p, method = "fdr")
df <- df[df$padj<0.05,]
df <- df[order(df$padj, decreasing = F),]

load("~/Documents/Work/data.nosync/tcga_sola/beta_ki67.Rdata")
load("2-expr-corr/1-output/pheno-corr.Rdata")
intersect <- intersect(normal$barcode1, colnames(beta_merged))
beta <- beta_merged[,intersect]
pheno <- normal[match(intersect, normal$barcode1),]
identical(pheno$barcode1, colnames(beta))
rm(beta_merged, normal);gc()

source("../../src/DNAm-tools/coerce_numeric.R")
beta <- coerce_numeric(beta)

# Score from df
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
pheno$idx_inverted <- tmp

t <- pheno |> 
  dplyr::group_by(tissue_simplified) |> 
  dplyr::count() |> 
  dplyr::filter(n > 3)

pheno |> 
  dplyr::filter(tissue_simplified %in% t$tissue_simplified) |> 
  ggplot(aes(x = idx_inverted,
             y = log2(ki67+0.00001),
             colour = ic)) +
  geom_point() +
  geom_smooth(method = "lm",
              se = F) +
  ggpubr::stat_cor() +
  facet_wrap(~tissue_simplified) # positive association, although not significant for all tissues.

save(pheno, file = "6-functional/1-proliferation/3-output/pheno_tcga.Rdata")
