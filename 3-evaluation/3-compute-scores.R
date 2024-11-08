# Compute scores

setwd("~/Dropbox/eca/sola/7-manuscript/manuscript/functional_aging_sigs/")
# libraries
library(dplyr)
source("0-data/src/beta_params.R")

# Compute scores in TCGA data ---------

# load pheno
load("~/Dropbox/eca/sola/1-download/2-output/data.Rdata")
# load("1-functional-association/4-tcga-data/2-output/data.Rdata")

# Keep those with at least 5 normal samples
n <- data |> 
  dplyr::group_by(project, type) |>
  dplyr::count() |> 
  tidyr::pivot_wider(id_cols = "project",
                     names_from = type, values_from = n) |> 
  dplyr::filter(Control >=5)

tmpdat <- data |> 
  dplyr::filter(project %in% n$project | project %in% c("TCGA-OV", "TCGA-CESC")) |> 
  dplyr::mutate(accession = project,
                celltype = tissue_simplified,
                group = "current cancer",
                type,
                basename = barcode1) |>
  dplyr::select(project, group, accession, celltype, type, basename, age)

# get beta
projects <- c(unique(n$project), "TCGA-OV", "TCGA-CESC")
# directory <- "<your download directory>/"
directory <- "~/Dropbox/data/tcga/"

for (p in projects){
  cat("Starting ", p, "...")
  load(paste0(directory, p, "/beta.Rdata"))
  
  tmp <- beta_params(beta)
  if(p == projects[1]){
    res <- as.data.frame(tmp)
  } else {
    res <- rbind(res, as.data.frame(tmp))
  }
  
  cat("done\n")
}

intersect <- intersect(tmpdat$basename, rownames(res))
data_tcga <- cbind(tmpdat[match(intersect, tmpdat$basename),],
                   res[intersect,])


# Compute scores in senescence datasets ---------
load("1-functional-association/1-senescence/2-qc/out/pheno.Rdata")
tmp <- pheno |> 
  tibble::rownames_to_column("id") |> 
  dplyr::mutate(type = case_when(type=="control" ~ "Control", 
                                 type!="control" & induction == "replicative" ~"Replicative senescence",
                                 type!="control" & induction == "irradiation" ~"Irradiation senescence",
                                 type!="control" & induction == "replicative-oncogene" & id %in% c("GSM2420518",
                                                                                                   "GSM2420519",
                                                                                                   "GSM2420520") ~ "Oncogene-induced senescence",
                                 type!="control" & induction == "replicative-oncogene" & !id %in% c("GSM2420518",
                                                                                                   "GSM2420519",
                                                                                                   "GSM2420520") ~ "Replicative senescence"),
                accession = set,
                celltype = celltype,
                group = "senescence") |>
  dplyr::select(group, accession, celltype, type, basename)
load("1-functional-association/1-senescence/2-qc/out/beta.Rdata")

beta <- beta[,tmp$basename]
identical(colnames(beta), tmp$basename)
res <- beta_params(beta)
data_sen <- cbind(tmp, as.data.frame(res))

# Compute scores in proliferation datasets ---------
load("1-functional-association/2-proliferation/3-output/pheno_filtered.Rdata")
directory <- "<your download directory/"
directory <- "~/Documents/Work/data.nosync/GEO/"
# load beta
load(paste0(directory, "proliferation/beta_merged.Rdata"))

tmp <- pheno |> 
  dplyr::filter(subexperiment != "MMC") |> 
  dplyr::mutate(type = case_when(type == "Control" ~ "Control",
                                 type != "Control" & subexperiment == "MMC" ~ "Reduced proliferation (mitomycin C)",
                                 type != "Control" & subexperiment == "Serum withdrawal" ~ "Reduced proliferation (serum withdrawal)"),
                t = days_in_culture, 
                pd = pd,
                accession = "GSE197512",
                celltype = tolower(tissue),
                group = "proliferation") |>
  dplyr::select(group, accession, celltype, type, basename, t, pd)

beta <- beta_merged[,tmp$basename]
identical(colnames(beta), tmp$basename)
res <- beta_params(beta)
data_prol <- cbind(tmp, as.data.frame(res))

# Add Reprogramming ---------
load("3-evaluation/2-additional-downloads/1-reprogramming/1-output/pheno_fibroblast.Rdata")
data_fib <- pheno |> 
  dplyr::filter(!is.na(type))

load("3-evaluation/2-additional-downloads/1-reprogramming/1-output/pheno_neuronal.Rdata")
data_neuron <- pheno

# Add Cervical ---------
load("3-evaluation/2-additional-downloads/2-cervical-cancer/1-output/pheno_cervical.Rdata")
data_cerv <- pheno
data_cerv$group <- "current cancer"

# Add pathological response --------
load("3-evaluation/2-additional-downloads/3-pathological-response/1-output/pheno_clean.Rdata")
data_prog <- pheno |>
  dplyr::filter(timepoint == "A") |> 
  dplyr::mutate(type = ifelse(prognosis == "Responder", "Control", "Non-responder"),
                group = "response",
                celltype = "breast",
                accession = "GSE184159") |> 
    dplyr::select(group, age, ic, accession, celltype, type, idx_AGE:Eosino)

# For all other datasets, values were equally computed as follows
# res <- as.data.frame(beta_params(beta))
load("0-data/data_assessment.Rdata")

# Merge all data
data <- plyr::rbind.fill(data_tcga, data_prol, data_sen,
                         data_fib, data_neuron,
                         data_cerv,
                         data_prog,
                         data)

save(data, file = "3-evaluation/3-output/data.Rdata")