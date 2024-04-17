if(!require("eutopsQC")){
  devtools::install_github("chiaraherzog/eutopsQC", force = T)
}

library(eutopsQC)
library(here)

directory <- "<your download directory>"

# Caveat - intermediate files are relatively large and should only be done if >26 GB are available in RAM.
# If not, recommend preprocessing in smaller batches (e.g. 1/3 of pheno file each) and merging afterwards

eutopsQC::preprocessData(input = directory,
                         output = paste0(directory, "/proliferation/"),
                         array = "EPIC",
                         report = paste0(here(), "/1-functional-association/2-proliferation/2-qc/"),
                         pheno = paste0(here(), "/1-functional-association/2-proliferation/1-output/pheno.Rdata"),
                         find.files = T,
                         run.name = "Proliferation",
                         cores = 8)
