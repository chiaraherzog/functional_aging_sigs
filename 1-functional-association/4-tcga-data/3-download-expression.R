# Download TCGA matched expression data for normal samples
# Solid normal tissue; methylation (already downloaded) is joined with RNA-seq data
# Downloading using TCGA biolinks

library(dplyr)
library(TCGAbiolinks)
library(SummarizedExperiment)
library(sesame)
library(sesameData)
library(here)

# already downloaded data:
load("~/Dropbox/eca/sola/1-download/2-output/data.Rdata")
load("1-functional-association/4-tcga-data/2-output/data.Rdata")

directory <- "<your download directory>"

normal <- data |>
  filter(type == "Control")
projects <- unique(normal$project)

# Find gene expression (legacy = FALSE; using preprocessed data)
for (p in projects){
  dir.create(paste0(directory, "/TCGA/", p, "/expression"))
  setwd(paste0(directory, "/TCGA/", p, "/expression"))
  
  query.exp <- GDCquery(project = p,
                        data.category = "Transcriptome Profiling",
                        data.type = "Gene Expression Quantification", 
                        workflow.type = "STAR - Counts",
                        sample.type = c("Solid Tissue Normal"),
                        barcode = normal[normal$project == p,]$barcode2)
  
  GDCdownload(query.exp, method = "client", files.per.chunk = 10)
  data <- GDCprepare(query.exp)
  save(data, file = "data_exp.Rdata")
  
  dat <- assays(data)[[1]]
  save(dat, file = "expression_matrix.Rds")
  
  file.remove("gdc_client_configuration.dtt")
  file.remove("gdc-client_v1.6.1_OSX_x64.zip")
  file.remove("gdc_manifest.txt")
  file.remove("gdc-client")
  unlink("GDCdata/", recursive = T)
  gc()
}

