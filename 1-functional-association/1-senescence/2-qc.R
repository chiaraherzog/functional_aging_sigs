# QC on downloaded senescence files + merging

if(!require("eutopsQC")){
  devtools::install_github("chiaraherzog/eutopsQC", force = T)
}

library(eutopsQC)
library(here)

directory <- "~/Documents/Work/data.nosync/GEO/"
directory <- "<your data folder>/"

# EPIC
eutopsQC::preprocessData(input = directory,
                         output = paste0(directory, "/senescence-epic/"),
                         array = "EPIC",
                         report = paste0(here(), "/1-functional-association/1-senescence/2-qc/epic/"),
                         pheno = paste0(here(), "/1-functional-association/1-senescence/1-output/pheno_epic.Rdata"),
                         find.files = T,
                         run.name = "Senescence-EPIC",
                         cores = 8)

# 450k
eutopsQC::preprocessData(input = directory,
                         output = paste0(directory, "/senescence-450k/"),
                         array = "450k",
                         report = paste0(here(), "/1-functional-association/1-senescence/2-qc/450k/"),
                         pheno = paste0(here(), "/1-functional-association/1-senescence/1-output/pheno_450k.Rdata"),
                         find.files = T,
                         run.name = "Senescence-450k",
                         cores = 8)

# Check outputs for good QC;

# Merge EPIC + 450k files
load(paste0(here(), "/1-functional-association/1-senescence/1-output/pheno_450k.Rdata"))
load(paste0(here(), "/1-functional-association/1-senescence/1-output/pheno_epic.Rdata"))
pheno <- rbind(pheno_epic, pheno_450k)

# Betas
load(paste0(directory, "/senescence-epic/beta_merged.Rdata"))
beta_epic <- beta_merged
load(paste0(directory, "/senescence-450k/beta_merged.Rdata"))
beta_450k <- beta_merged

intersect <- intersect(rownames(beta_epic), rownames(beta_450k))
beta <- cbind(beta_epic[intersect,],
              beta_450k[intersect,])
identical(colnames(beta), pheno$basename)
save(pheno, file = "1-functional-association/1-senescence/2-qc/out/pheno.Rdata")
save(beta, file = "1-functional-association/1-senescence/2-qc/out/beta.Rdata")
