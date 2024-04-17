# Download TCGA data for methylation + expression correlation
suppressPackageStartupMessages(library(TCGAbiolinks))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(SummarizedExperiment))
suppressPackageStartupMessages(library(here))

# TCGA projects -------------------------

# Get information on all projects
projects <- TCGAbiolinks::getGDCprojects()
projects <- projects[grepl("TCGA", projects$id),] # 33 projects

directory <- "<your download directory>"

# Get data for all projects
for (id in projects$id){
  
  cat("Starting ", id, " (", which(projects$id==id), "/", length(projects$id), ").\n", sep = "")
  
  # Query available files
  suppressMessages(query <- GDCquery(project = id,
                                     data.category = "DNA Methylation",
                                     platform = "Illumina Human Methylation 450",
                                     sample.type = c("Primary Tumor",
                                                     "Solid Tissue Normal",
                                                     "Blood Derived Normal",
                                                     "Primary Blood Derived Cancer - Peripheral Blood",
                                                     "Primary Blood Derived Cancer - Bone Marrow",
                                                     "Buccal Cell Normal",
                                                     "Bone Marrow Normal"),
                                     data.type = "Methylation Beta Value"))
  dat <- query$results[[1]]
  
  query_pheno <- GDCquery_clinic(project = id,
                                 type = "clinical")
  ind <- match(dat$cases.submitter_id, query_pheno$submitter_id)
  tcga_pheno <- query_pheno[ind,]
  save(tcga_pheno, dat, file = paste0(directory, "/TCGA/pheno_", substr(id, 6, nchar(id)), ".Rdata"))
  
  cat("Starting download (", as.character(Sys.time()), ").\n", sep = "")
  
  setwd(paste0(directory, "/TCGA/", id))
  GDCdownload(query, method = "client", files.per.chunk = 10) # RIP, this takes long!
  suppressMessages(data <- GDCprepare(query))
  beta <- assays(data)[[1]]
  beta <- na.omit(beta)
  cat("Done downloading, saving file.\n")
  save(beta, file = paste0(directory, "/TCGA/", id, "/beta.Rdata"))
  
  # delete files
  unlink("GDCdata/", recursive = T)
  suppressMessages(file.remove(c("gdc-client_v1.6.1_OSX_x64.zip",
                                   "gdc_manifest.txt",
                                   "gdc-client",
                                   "gdc_client_configuration.dtt")))
 setwd(here())
 # remove intrim files and clean
 rm(pheno, beta, data, tcga_pheno, dat, query_pheno, query);gc()
 
 
}



