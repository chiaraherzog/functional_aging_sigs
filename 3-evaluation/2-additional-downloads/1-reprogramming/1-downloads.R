# Download two reprogramming datasets

# libraries
library(GEOquery)

# source
source("0-data/src/beta_params.R")

directory <- "<your download directory>/"
directory <- "~/Documents/Work/data.nosync/GEO/"

# 1. GSE147436 ----------------------
# Neuronal data repair/reprogramming

# Pheno data
gse <- getGEO("GSE147436")
pheno <- pData(gse$GSE147436_series_matrix.txt.gz)
save(pheno, file = "3-evaluation/2-additional-downloads/1-reprogramming/1-output/GSE147436.Rdata")

# Beta
beta <- exprs(gse$GSE147436_series_matrix.txt.gz)
dir.create(paste0(directory, "GSE147436"))
save(beta, file = paste0(directory, "GSE147436/beta.Rdata"))
beta <- na.omit(beta)

dat <- beta_params(beta)
identical(rownames(dat), rownames(pheno))

# Compute
pheno <- cbind(pheno, dat)

# Tidy: keep only those without vincristine damage +- OSK
pheno <- pheno |> 
  dplyr::filter(`status/time point:ch1` == "no vincristine damage") |> 
  dplyr::mutate(accession = "GSE147436",
                celltype = "SH-SY5Y neurons",
                group = "OSK reprogramming",
                type = case_when(`osk:ch1` == "yes" ~ "OSK",
                              `osk:ch1` == "no" ~ "Control")) |> 
  dplyr::select(group, accession, celltype, type, any_of(colnames(as.data.frame(dat))))

save(pheno, file = "3-evaluation/2-additional-downloads/1-reprogramming/1-output/pheno_neuronal.Rdata")


# 2. GSE54848 ----------------------
# OSKM at different durations

# Pheno data
gse <- getGEO("GSE54848")
pheno <- pData(gse[[1]])
save(pheno, file = "3-evaluation/2-additional-downloads/1-reprogramming/1-output/GSE54848.Rdata")

# Beta
beta <- exprs(gse[[1]])
dir.create(paste0(directory, "GSE54848"))
save(beta, file = paste0(directory, "GSE54848/beta.Rdata"))
beta <- na.omit(beta)

dat <- beta_params(beta)
identical(rownames(dat), rownames(pheno))

# Compute
pheno <- cbind(pheno, dat)

pheno$title

# Tidy: 
pheno <- pheno |> 
  dplyr::mutate(type = case_when(`cell type:ch1` == "Human dermal fibroblast" ~ "Control",
                                 `cell type:ch1` %in% c("induced pluripotent stem cells")~ "iPSC",
                                 `cell type:ch1` %in% c("human Embryonic stem cells") ~ "ESC",
                                 grepl("Day", title) & source_name_ch1 == "Intermediate reprogrammed cells" ~ "intermediate reprogrammed cells"),
                t = case_when(type == "Control" ~ 0,
                              type == "intermediate reprogrammed cells" ~ readr::parse_number(stringr::str_split(title, "_", simplify = T)[,1]),
                              type %in% c("iPSC", "ESC") ~ 40),
                accession = "GSE54848",
                celltype = "fibroblast",
                group = "OSK reprogramming") |> # exemplary "t = 40" for iPSC
  dplyr::select(group, accession, celltype, type, t, any_of(colnames(as.data.frame(dat))))

save(pheno, file = "3-evaluation/2-additional-downloads/1-reprogramming/1-output/pheno_fibroblast.Rdata")
