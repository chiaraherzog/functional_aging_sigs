# Download additional cervical methylation to add to TCGA for cervical cancer - sparse dataset

# libraries
library(GEOquery)
library(dplyr)
library(ChAMP)
library(ggplot2)
library(uwot)

directory <- "<your download directory>"
directory <- "~/Documents/Work/data.nosync/GEO/"

# Pheno
gse <- getGEO("GSE211668")
pheno <- pData(gse[[1]])
save(pheno, file = "3-evaluation/2-additional-downloads/2-cervical-cancer/1-output/GSE211668.Rdata")

# Download signal intensities (no raw IDATs available); no exprs file
dir.create(paste0(directory, "/GSE211668/"))
setwd(paste0(directory, "/GSE211668/"))
getGEOSuppFiles("GSE211668", filter_regex = "GSE211668_UCL_Cervical_MatrixSignal.txt.gz", makeDirectory = F)
gunzip("GSE211668_UCL_Cervical_MatrixSignal.txt.gz")

# Process with GEORawFile
x <- readLines("GSE211668_UCL_Cervical_MatrixSignal.txt", n = 3)

# Reading in raw intensities from GEO file
# Using [.] around dots as otherwise grep within readGEORawFile will not recognise them as unique strings
Mset <- minfi::readGEORawFile(filename = "GSE211668_UCL_Cervical_MatrixSignal.txt",
                              sep = "\t",
                              Uname = "Unmethylated signal",
                              Mname = "Methylated signal",
                              array = "IlluminaHumanMethylation450k",
                              annotation = "ilmn12.hg19",
                              showProgress = TRUE)

# QC & Pipeline

# Define global thresholds
INTENSITY_THRESHOLD <- 9.5     # minimum median intensity required
DETECTION_P_THRESHOLD <- 0.01  # maximum detection p-value
FAILED_PROBE_THESHOLD <- 0.1   # maximum proportion of failed probes per sample

qc <- getQC(Mset)
plotQC(qc) # Most samples have homogenous intensity.

# Filter any samples with median (un)methylated intensity less than threshold
low_intensity_samples <- rownames(qc)[qc$mMed<INTENSITY_THRESHOLD | qc$uMed<INTENSITY_THRESHOLD] # one sample below threshold

# Read in detP
colnames <- strsplit(readLines("GSE211668_UCL_Cervical_MatrixSignal.txt", n = 1), "\t")[[1]]
select <- sort(grep("Detection Pval",colnames))
detP <- data.table::fread("GSE211668_UCL_Cervical_MatrixSignal.txt",
                          sep = "\t",
                          select = select)

# Filter samples with too many failed probes
failed_samples <- colnames(detP)[colSums(detP>DETECTION_P_THRESHOLD) > (nrow(detP) * FAILED_PROBE_THESHOLD)] # Sample with low intensity fails
rm(detP, colnames, select, qc);gc()

# Convert to RatioSet and then beta
RSet <- ratioConvert(Mset, what = "both", keepCN = TRUE)
rm(Mset);gc()

beta <- getBeta(RSet)
beta <- na.omit(beta)
rm(RSet);gc()

densityPlot(beta) # Plot densities # One odd sample; remove because otherwise ChAMP fails :(
z <- as.data.frame(t(beta[1:10000,])) |> 
  tibble::rownames_to_column("id") |> 
  tidyr::pivot_longer(cols = -c(id),
                      names_to = "cg",
                      values_to = "beta") 
library(ggplot2)
p <- z |> 
  ggplot(aes(x = beta,
             colour= id)) +
  geom_density()
plotly::ggplotly(p)  # find odd sample out. it is 376 AT which should have been removed before...

beta <- beta[,!colnames(beta) %in% c("376 A T")]
densityPlot(beta) # Plot densitiesnot looks ok.

# ChAMP Normalisation
norm <- champ.norm(beta, arraytype = "450k", method = "BMIQ", cores = 8)
beta <- norm
densityPlot(beta)

# rename columns
pheno$title <- gsub("^\\D+","",pheno$title)
pheno$title <- gsub("\\.", "", pheno$title)
ind <- match(pheno$title, colnames(beta)) # one sample missing -> removed sample.
pheno2 <- pheno[match(colnames(beta), pheno$title),]
identical(colnames(beta), pheno2$title)

# Rename colnames of beta
colnames(beta) <- rownames(pheno2)
save(beta, file = "beta.Rdata")
intersect <- intersect(colnames(beta), rownames(pheno))
pheno <- pheno[intersect,]
beta <- beta[,intersect]
identical(rownames(pheno), colnames(beta)) 

# Look at batch effects ~ UMAP ; prep pheno further--------
# load("beta.Rdata")
b <- beta[1:100000,]
umap <- uwot::umap(t(b))

ind <- match(rownames(umap), rownames(pheno2))
pheno2$umap1 <- umap[,1]
pheno2$umap2 <- umap[,2]

pheno |> 
  ggplot(aes(x = umap1, 
             y = umap2,
             colour = type)) +
  geom_point()

# compute indices
source("~/Dropbox/eca/sola/7-manuscript/manuscript/functional_aging_sigs/0-data/src/beta_params.R")
res <- beta_params(beta)
identical(rownames(res), rownames(pheno))
pheno <- cbind(pheno, res)

p <- pheno |> 
  dplyr::mutate(type = case_when(`disease state:ch1` == "Normal" ~ "Control",
                                 `disease state:ch1` == "Tumour" ~ "Cancer"),
                id = stringr::str_split(title, "\ ", simplify = T)[,1],
                accession = "GSE211668",
                celltype = "cervix",
                group = "cancer") |>
  dplyr::select(group, accession, celltype, type, id, any_of(colnames(as.data.frame(res))))
pheno <- p

save(pheno, file = "3-evaluation/2-additional-downloads/2-cervical-cancer/1-output/pheno_cervical.Rdata")
identical(colnames(beta), rownames(pheno))
setwd("~/Dropbox/eca/sola/7-manuscript/manuscript/functional_aging_sigs/")
source("0-data/src/beta_params.R")
tmp <- beta_params(beta)

pheno_min <- pheno |> dplyr::select(-any_of(colnames(tmp)))
pheno <- cbind(pheno_min, tmp)
save(pheno, file = "3-evaluation/2-additional-downloads/2-cervical-cancer/1-output/pheno_cervical.Rdata")
