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
load("1-functional-association/1-senescence/3-output/coef.Rdata")
df$padj <- p.adjust(df$p, method = "fdr")
df <- df[df$padj<0.05 & !df$chr %in% c("chrX", "chrY"),] # filter by pval and chr (no sex chromosomes)
df <- df[order(df$padj, decreasing = F),]
head(df)

# pheno
load("1-functional-association/1-senescence/2-qc/out/pheno.Rdata")
load("1-functional-association/1-senescence/2-qc/out/beta.Rdata")
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
  ggplot(aes(x = cg,
             y = beta,
             fill = type)) +
  geom_boxplot() + 
  # geom_point(aes(colour = celltype)) +
  theme_classic() +
  facet_wrap(~set,
             scales = "free") +
  theme(legend.position = "none") # senescent cells tend to lose methylation at top 10 sites

hist(df$coef)

df |> 
  ggplot(aes(x = coef)) +
  geom_histogram(bins = 100) +
  xlim(c(-0.4, 0.4)) +
  theme_minimal()

# Generate + plot a score
intersect <- intersect(rownames(beta), df$cg)
b <- beta[intersect,]
w <- df[match(intersect, df$cg),]
w <- ifelse(w$coef>0, 1, -1)
names(w) <- rownames(b)
tmp <- apply(b, 2, function(i){
  return(sum(i*w)/length(w))
})

identical(names(tmp), pheno$basename)
pheno$idx <- tmp

pheno |> 
  ggplot(aes(x = type,
             y = idx)) +
  geom_boxplot() + 
  geom_point(aes(colour = celltype)) +
  theme_bw()# higher value indicates senescence

pheno |> 
  ggplot(aes(x = type,
             y = idx)) +
  geom_boxplot() + 
  geom_point(aes(colour = celltype)) +
  facet_wrap(~set) +
  theme_bw()# higher value indicates senescence, consistent across sets although absolute values are different across sets

# Save pheno with idx
save(pheno, file = "1-functional-association/1-senescence/4-output/pheno_senescence_association.Rdata")

# Run GSEA for those that are significant
cpg.pval <- df$p
names(cpg.pval) <- df$cg

# GO with methylGSA
res1 = methylglm(cpg.pval = cpg.pval, minsize = 50, 
                 maxsize = 500, GS.type = "GO",array.type = "450K")
head(res1, 15)
hist(res1$padj)
x <- res1[res1$padj<0.05,]
x$Description # exoribonucleas, base excision repair, ...

# KEGG with methylGSA
res2 = methylglm(cpg.pval = cpg.pval, minsize = 50, 
                 maxsize = 500, GS.type = "KEGG",array.type = "450K")
x <- res2[res2$padj<0.05,]
x$Description # none remain significant
save(res1, res2, file = "1-functional-association/1-senescence/4-output/gsa.Rdata")

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
     enrichpath, file = "1-functional-association/1-senescence/4-output/clusterProfiler.Rdata")
rm(list=ls())