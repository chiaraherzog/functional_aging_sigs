# Download senescence data

# libraries and settings
library(GEOquery)
library(here)

# set options
options(timeout = max(300, getOption("timeout")))
options(download.file.method.GEOquery = "wget")

# set data folder

directory <- "<your data folder>/"

# get source code 
source("0-data/src/downloadGEO.R")

# 1. Radiation/DNA-Damage-induced senescence -------------------------
# GSE112812: irradiated fibroblasts, strain 1; 450 K 
downloadGEO(accession = "GSE112812", directory)

# GSE112873: irradiated fibroblasts, strain 2; 450 K 
# "We did not observe epigenetic differences between stable and unstable irradiated clones. " - Flunkert et al. 2018
downloadGEO(accession = "GSE112873", directory)

# 2. Oncogene-induced senescence -------------------------
# GSE91069: transformed fibroblasts; 450k array 
downloadGEO(accession = "GSE91069", directory)

# 3. Replicative senescence -------------------------
# GSE82234: passaged cells; 450k array 
downloadGEO(accession = "GSE82234", directory)

# GSE131280: passaged cells; EPIC array 
downloadGEO(accession = "GSE131280", directory)

# GSE116375: passaged cells; 450k array
downloadGEO(accession = "GSE116375", directory)

rm(list = ls())

# 4. Merge pheno files  -------------------------

# Generate one single pheno file per array type that will be used for filtering files and QC
directory <- "<your data folder>/"


# GSE112812
accession = "GSE112812"
load(paste0(directory, accession, "/pheno.Rdata"))

GSE112812 <- pheno |> 
  dplyr::filter((`group:ch1` == "non-irradiated" & `stability:ch1` == "NA") |(`group:ch1` == "irradiated" & `stability:ch1` == "stable")) |> 
  dplyr::mutate(type = ifelse(`group:ch1` == "irradiated", "senescent", "control"),
                set = "GSE112812",
                celltype = "fetal fibroblast",
                induction = "irradiation",
                basename = gsub("_Grn.idat.gz", "", basename(supplementary_file))) |> 
  dplyr::select(basename, type, set, celltype, induction)

# GSE112873
accession = "GSE112873"
load(paste0(directory, accession, "/pheno.Rdata"))

GSE112873 <- pheno  |> 
  dplyr::filter((`group:ch1` == "non-irradiated" & `stability:ch1` == "NA") |(`group:ch1` == "irradiated" & `stability:ch1` == "stable")) |> 
  dplyr::mutate(type = ifelse(`group:ch1` == "irradiated", "senescent", "control"),
                set = "GSE112873",
                celltype = "fetal fibroblast",
                induction = "irradiation",
                basename = gsub("_Grn.idat.gz", "", basename(supplementary_file))) |> 
  dplyr::select(basename, type, set, celltype, induction)

# GSE91069
accession = "GSE91069"
load(paste0(directory, accession, "/pheno.Rdata"))

GSE91069 <- pheno |> 
  dplyr::filter(grepl("Early Passage|^Sene|^Oncogene", source_name_ch1)) |> 
  dplyr::mutate(type = ifelse(grepl("Early Passage", source_name_ch1), "control", "senescent"),
                set = "GSE91069",
                celltype = "fetal fibroblast",
                induction = "replicative-oncogene",
                basename = gsub("_Grn.idat.gz", "", basename(supplementary_file))) |> 
  dplyr::select(basename, type, set, celltype, induction)

# GSE82234
accession = "GSE82234"
load(paste0(directory, accession, "/pheno.Rdata"))

GSE82234 <- pheno |> 
  dplyr::mutate(type = ifelse(`passage:ch1` == "passage 4", "control", "senescent"),
                set = "GSE82234",
                celltype = "human umbilical vein endothelial cells",
                induction = "replicative",
                basename = gsub("_Grn.idat.gz", "", basename(supplementary_file))) |> 
  dplyr::select(basename, type, set, celltype, induction)

# GSE131280
accession = "GSE131280"
load(paste0(directory, accession, "/pheno.Rdata"))

GSE131280 <- pheno |> 
  dplyr::mutate(type = ifelse(as.numeric(`days_grown:ch1`) < 50, "control", 
                              ifelse(as.numeric(`days_grown:ch1`) > 110, "senescent", NA))) |> 
  dplyr::filter(!is.na(type) & `glucose_treatment:ch1` == "normal") |> 
  dplyr::mutate(set = "GSE131280",
                celltype = "adult fibroblast",
                induction = "replicative",
                basename = gsub("_Grn.idat.gz", "", basename(supplementary_file))) |> 
  dplyr::select(basename, type, set, celltype, induction)

# GSE116375
accession = "GSE116375"
load(paste0(directory, accession, "/pheno.Rdata"))

GSE116375 <- pheno |> 
  dplyr::filter(!grepl("OxBS", title)) |> 
  dplyr::mutate(type = ifelse(`passage:ch1` == "passage 4", "control", "senescent")) |> 
  dplyr::rename(celltype = `cell type:ch1`) |> 
  dplyr::mutate(set = "GSE116375",
                induction = "replicative",
                basename = gsub("_Grn.idat.gz", "", basename(supplementary_file))) |> 
  dplyr::select(basename, type, set, celltype, induction)

accession = "GSE116375"
load(paste0(directory, accession, "/pheno.Rdata"))

pheno <- rbind(GSE112812,
               GSE112873,
               GSE116375,
               GSE131280,
               GSE82234,
               GSE91069)

# Split by array
pheno_epic <- pheno |> 
  dplyr::filter(set %in% c("GSE131280"))

pheno_450k <- pheno |> 
  dplyr::filter(!set %in% c("GSE131280"))

save(pheno_epic, file = "1-functional-association/1-senescence/1-output/pheno_epic.Rdata")
save(pheno_450k, file = "1-functional-association/1-senescence/1-output/pheno_450k.Rdata")