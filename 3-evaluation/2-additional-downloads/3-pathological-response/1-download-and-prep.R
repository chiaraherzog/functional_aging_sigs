# Pathological response dataset

# libraries
library(GEOquery)
library(dplyr)
library(ChAMP)
library(ggplot2)
library(uwot)

directory <- "<your download directory>"
directory <- "~/Documents/Work/data.nosync/GEO/"

# Pheno
gse <- getGEO("GSE184159")
pheno <- pData(gse[[1]])
save(pheno, file = "3-evaluation/2-additional-downloads/3-pathological-response/pheno_GSE184159.Rdata")

# Beta
beta <- exprs(gse[[1]])


# Download signal intensities (no raw IDATs available); no exprs file
dir.create(paste0(directory, "/GSE184159/"))
setwd(paste0(directory, "/GSE184159/"))
getGEOSuppFiles("GSE184159", filter_regex = "GSE184159_raw_detp.csv.gz", makeDirectory = F)
gunzip("GSE184159_raw_detp.csv.gz")

# Process with GEORawFile
x <- readLines("GSE184159_raw_detp.csv", n = 3)
# Reading in raw intensities from GEO file
# Using [.] around dots as otherwise grep within readGEORawFile will not recognise them as unique strings
Mset <- minfi::readGEORawFile(filename = "GSE184159_raw_detp.csv",
                              sep = ",",
                              Uname = "^unmeth_",
                              Mname = "^meth_",
                              array = "IlluminaHumanMethylationEPIC",
                              annotation = "ilm10b4.hg19",
                              showProgress = TRUE)

# QC & Pipeline

# Define global thresholds
INTENSITY_THRESHOLD <- 9.5     # minimum median intensity required
DETECTION_P_THRESHOLD <- 0.01  # maximum detection p-value
FAILED_PROBE_THESHOLD <- 0.1   # maximum proportion of failed probes per sample

qc <- getQC(Mset)
plotQC(qc) 

# Filter any samples with median (un)methylated intensity less than threshold
low_intensity_samples <- rownames(qc)[qc$mMed<INTENSITY_THRESHOLD | qc$uMed<INTENSITY_THRESHOLD] # one sample below threshold

# Read in detP
colnames <- strsplit(readLines("GSE184159_raw_detp.csv", n = 1), ",")[[1]]
select <- sort(grep("^p_det",colnames))
detP <- data.table::fread("GSE184159_raw_detp.csv",
                          sep = ",",
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

densityPlot(beta) # densities do look a bit odd

# ChAMP Normalisation
norm <- champ.norm(beta, arraytype = "EPIC", method = "BMIQ", cores = 8)
densityPlot(norm) # looks improved
beta <- norm

colnames(beta)
pheno$title

# rename columns
pheno$title2 <- gsub("Tissue_|Normal_","", pheno$title)
ind <- match(pheno$title2, colnames(beta)) 
pheno2 <- pheno[match(colnames(beta), pheno$title2),]
identical(colnames(beta), pheno2$title2)

# Rename colnames of beta
colnames(beta) <- rownames(pheno2)
save(beta, file = "beta.Rdata")
intersect <- intersect(colnames(beta), rownames(pheno))
pheno <- pheno[intersect,]
beta <- beta[,intersect]
identical(rownames(pheno), colnames(beta)) 
save(pheno, file ="~/Dropbox/eca/sola/7-manuscript/manuscript/molecular-aging-sigs/3-evaluation/2-additional-downloads/3-pathological-response/1-output/pheno.Rdata")

# Look at batch effects ~ UMAP ; prep pheno further--------
# load("beta.Rdata")
b <- beta[1:30000,]
umap <- uwot::umap(t(b))

ind <- match(rownames(umap), rownames(pheno))
pheno2 <- pheno
pheno2$umap1 <- umap[,1]
pheno2$umap2 <- umap[,2]

pheno2 |> 
  ggplot(aes(x = umap1, 
             y = umap2,
             colour = `source_name_ch1`)) +
  geom_point()

# Compute scores and prep pheno (appending info from Supplementary Table 3 of original publication https://doi.org/10.1186/s13148-021-01210-6)

source("0-data/src/beta_params.R")
identical(colnames(beta), rownames(pheno))
res <- beta_params(beta)

# Load pheno info
dat <- readxl::read_xlsx("3-evaluation/2-additional-downloads/3-pathological-response/1-output/13148_2021_1210_MOESM3_ESM.xlsx", sheet = 2, skip = 2) |> 
  janitor::clean_names() |> 
  dplyr::filter(!is.na(age)) |> # baseline info 
  dplyr::mutate(id = stringr::str_split(sample_name, "_", simplify = T)[,1],
                prognosis = ifelse(prognosis == "Responder", "Responder", "Non-responder"),
                vital = ifelse(dead == "Yes", 1, NA)) |> 
  dplyr::select(id, age, relapse, time_to_relapse_days, prognosis, disease_free_survival_months,
                vital, days_on_study_without_death)


pheno <- pheno |> 
  dplyr::mutate(timepoint = `timepoint:ch1`,
                tissue = ifelse(source_name_ch1 == "Breast_Cancer", "Cancer tissue",
                                "Normal tissue"),
                response = case_when(`response:ch1` == "Partial" ~ "Partial response",
                                     `response:ch1` == "Non" ~ "No response",
                                     `response:ch1` == "Complete" ~ "Complete response"),
                id = stringr::str_split(title, "_", simplify = T)[,1]) |> 
  dplyr::select(tissue, timepoint, response, id) |> 
  dplyr::left_join(dat)

pheno <- cbind(pheno, res)
save(pheno, file = "3-evaluation/2-additional-downloads/3-pathological-response/1-output/pheno_clean.Rdata")
