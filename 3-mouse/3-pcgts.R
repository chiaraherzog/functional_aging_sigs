# Define mouse PCGTs based on human genes
library(IlluminaMouseMethylationanno.12.v1.mm10)

# Load mouse array annotation
anno <- getAnnotation(IlluminaMouseMethylationanno.12.v1.mm10)
anno$mgene <- stringr::str_split(anno$GeneName_NCBI, ";", simplify = T)[,1] # keep only first as illumina gives multiple

# Load paper data and filter single occupancy
dat <- readxl::read_xls("0-data/files/mmc10.xls", skip = 1, col_types = "text")
dat <- dat[(dat$`Suz12 Occupancy`==1 | dat$`Eed Occupancy` == 1 | dat$`H3K27me3 Occupancy` == 1),]

# Get orthologues via Entrez
ortho <- dat$`EntrezGene ID` |> 
  babelgene::orthologs(species = 'mouse') |> 
  dplyr::mutate(`EntrezGene ID` = as.character(human_entrez)) |> 
  dplyr::rename(mgene = symbol)
dat <- dat |> 
  dplyr::left_join(dplyr::select(ortho, `EntrezGene ID`, mgene))

# Overlap with mouse array
pcgt_single <- intersect(anno$mgene, dat$mgene) # 1405 genes
cpgs_single <- anno[anno$mgene %in% pcgt_single & !is.na(anno$mgene),] #20114 CpGs

# Filtering by TSS200
tss <- stringr::str_split(cpgs_single$Feature_NCBI, ";", simplify = T)[,1]
cpgs_single_tss200 <- cpgs_single[tss == "tss_200" & !is.na(tss),] # 566 CpGs
sum(tss == "tss_200" & !is.na(tss))
nrow(cpgs_single_tss200)

save(pcgt_single, cpgs_single, cpgs_single_tss200, file = "3-mouse/3-output/pcgt.Rdata")

# Keep only fetal unmethylated CpGs

# Load and filter relevant samples
beta_final <- readRDS("~/Dropbox/eca/mouse-project/1-method-dev/data/1-mouse-methylation-atlas/beta_final.Rds")
load("~/Dropbox/eca/mouse-project/1-method-dev/1-preprocessing/1-mouse-methlation-atlas/1-output/pheno.Rdata")

fetal <- pheno |> 
  dplyr::filter(grepl("Fetal", tissue)) |> 
  droplevels()
beta <- beta_final[cpgs_single_tss200$Name,match(rownames(fetal), substr(colnames(beta_final), 1, 10))]
identical(rownames(fetal), substr(colnames(beta), 1, 10))
table(fetal$tissue)

adult <- pheno |> 
  dplyr::filter(!grepl("Fetal", tissue) & description %in% c("Tissue", "Tissue Biology")) |> 
  droplevels()
beta_adult <- beta_final[cpgs_single_tss200$Name,match(rownames(adult), substr(colnames(beta_final), 1, 10))]
identical(rownames(adult), substr(colnames(beta_adult), 1, 10))
table(adult$tissue)

# getting rowmeans and selecting CpGs with mean methylation of 0.2
mean_fetal <- rowMeans(beta)
mean_adult <- rowMeans(beta_adult)
hist(mean_fetal)
hist(mean_adult)

y <- mean_fetal[mean_fetal < 0.2] # selecting CpGs with mean methylation below 0.2
length(y) # 479 are unmethylated in fetal tissue
fetal_zero <- names(y)

# Get single occupancy CpGs (TSS200)
sum(fetal_zero %in% rownames(cpgs_single_tss200))
pcgt_cpgs_tss200_fetalzero <- cpgs_single_tss200[rownames(cpgs_single_tss200) %in% fetal_zero,]
save(pcgt_cpgs_tss200_fetalzero, file = "3-mouse/3-output/pcgt_cpgs_tss200_fetalzero.Rdata")
