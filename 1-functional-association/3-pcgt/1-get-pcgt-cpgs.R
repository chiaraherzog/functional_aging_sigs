# Identification of CpGs in TSS200 of PCGT Genes that are fetally unmethylated (mean beta < 0.2)

# libraries
library(readxl)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)

# 1. Identify PCGT CpG sites -------------
# Working from Lee et al. 2006 paper based on promoter occupancy of genes [single]

# First, get EPIC annotation
data <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19, what = "everything")

# tidy
genes <- data$UCSC_RefGene_Name
genes <- stringr::str_split(genes, ";", simplify = T)[,1] # keep only first as illumina gives multiple
tss <- data$UCSC_RefGene_Group
tss <- stringr::str_split(tss, ";", simplify = T)[,1] # keep only first as illumina gives multiple
data$gene <- genes 
length(unique(genes)) # 26650 unique genes
data$tss <- tss

# Load paper data and filter single occupancy
dat <- read_xls("0-data/files/mmc10.xls", skip = 1, col_types = "text")
dat <- dat[(dat$`Suz12 Occupancy`==1 | dat$`Eed Occupancy` == 1 | dat$`H3K27me3 Occupancy` == 1),]

# Overlap with genes on EPIC array
pcgt_single <- intersect(data$gene, dat$`Gene Name`) # 1343 genes
cpgs_single <- data[data$gene %in% pcgt_single,] #57840 CpGs

# Filtering by TSS200
cpgs_single_tss200 <- cpgs_single[cpgs_single$tss == "TSS200",] # only 4625 left

# Check how many of these are in original epiTOC
load("0-data/files/13059_2016_1064_MOESM5_ESM.Rdata")
sum(rownames(cpgs_single) %in% epiTOCcpgs.v) # 316/385 in epiToc

save(pcgt_single, cpgs_single, cpgs_single_tss200, file = "1-functional-association/3-pcgt/1-output/single_occupancy.Rdata")
sum(rownames(cpgs_single_tss200) %in% epiTOCcpgs.v) # 316/385 in epiToc

# 2. Check methylation in fetal tissue -------------
# keep only those with methylation < 0.2 as per original epiTOC paper
# Fetal tissue:
#   - cord blood [GSE72867]
#   - stomach, heart, tongue, kidney, liver, brain, thymus, spleen, lung, adrenal gland [GSE31848].

library(GEOquery)
directory <- "<your download directory>"
directory <- "~/Documents/Work/data.nosync/GEO/"
setwd(directory)

# Here, we are using exprs () data, as raw idats only available for cord blood but not the remaining tissues, but we only want to extract CpGs which are methylated 0 anyway which are likely robust to various preprocessing pipelines
cb <- getGEO("GSE72867") # Cord blood
cb_beta <- na.omit(exprs(cb[[1]]))
cb <- pData(cb[[1]])
identical(rownames(cb), colnames(cb_beta))

other <- getGEO("GSE31848") # other samples
other_beta <- na.omit(exprs(other[[1]]))
# filtering out fetal tissue for 'other dataset'
other_pheno <- pData(other[[1]])
other_pheno <- other_pheno[other_pheno$`fetal vs adult tissue:ch1`=="Fetal" & !is.na(other_pheno$`fetal vs adult tissue:ch1`),]
ind <- match(rownames(other_pheno), colnames(other_beta))
fetal <- other_beta[,ind]
adult <- other_beta[,-ind]
identical(rownames(other_pheno), colnames(fetal))

# creating one "fetal" matrix (beta)
intersect <- intersect(rownames(cb_beta), rownames(fetal))
beta <- cbind(cb_beta[intersect,],
              fetal[intersect,])
# merge pheno
pheno <- plyr::rbind.fill(cb, other_pheno)

# keep adult for reference
adult <- adult[intersect,]

# getting rowmeans and selecting CpGs with mean methylation of 0.2
mean_fetal <- rowMeans(beta)
mean_adult <- rowMeans(adult)
hist(mean_fetal) # checking distribution
# hist(mean_adult)

y <- mean_fetal[mean_fetal < 0.2] # selecting CpGs with mean methylation below 0.2
length(y)
fetal_zero <- names(y)

# Get single occupancy CpGs (TSS200)
sum(fetal_zero %in% rownames(cpgs_single_tss200))
pcgt_cpgs_tss200_fetalzero <- cpgs_single_tss200[rownames(cpgs_single_tss200) %in% fetal_zero,]

setwd(here::here())
save(pcgt_cpgs_tss200_fetalzero, file = "1-functional-association/3-pcgt/1-output/pcgt_cpgs_tss200_fetalzero.Rdata")

# get genes within TSS200 which have fetal zero methylation
genes <- pcgt_cpgs_tss200_fetalzero$UCSC_RefGene_Name
for (i in 1:length(genes)){
  genes[i] <- strsplit(genes[i], ";")[[1]][1]
}
pcgt_single_fetalzero <- genes[(genes %in% pcgt_single)]
save(pcgt_single_fetalzero, file = "1-functional-association/3-pcgt/1-output/pcgt_single_fetalzero.Rdata")


