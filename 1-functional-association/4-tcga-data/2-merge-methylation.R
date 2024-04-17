# Combine TCGA methylation data for all projects
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(EpiDISH))

# List TCGA projects -------------------------
projects <- TCGAbiolinks::getGDCprojects()
projects <- projects[grepl("TCGA", projects$id),] # 33 projects
pb <- txtProgressBar(min = 1, max = nrow(projects), style = 3)

directory <- "<your download directory>"

for (id in projects$id){
  setTxtProgressBar(pb, which(projects$id == id))
  
  if(!dir.exists(paste0(directory, "/TCGA/", id))){
    NULL
  } else {
    # Load Pheno
    load(paste0(directory, "/TCGA/pheno_", substr(id, 6, nchar(id)), ".Rdata"))
    
    # load Beta
    load(paste0(directory, "/TCGA/", id, "/beta.Rdata"))
    
    # Link submitter ids between beta/phenos
    dat <- dat[match(substr(colnames(beta), 1, 16), dat$sample.submitter_id),]
    tcga_pheno <- tcga_pheno[match(substr(dat$sample.submitter_id, 1, 12), tcga_pheno$submitter_id),]
    rownames(tcga_pheno) <- colnames(beta)
    tcga_pheno$barcode1 <- colnames(beta)
    tcga_pheno$barcode2 <- dat$sample.submitter_id
    tcga_pheno$sample_type <- dat$sample_type
    
    lookup <- c(age = "age_at_index",
                stage_ajcc = "ajcc_pathologic_stage",
                stage_figo = "figo_stage",
                stage_T = "ajcc_pathologic_t",
                stage_N = "ajcc_pathologic_n",
                stage_M = "ajcc_pathologic_m",
                tissue = "tissue_or_organ_of_origin",
                barcode3 = "submitter_id")
    
    cols <- c("project", "age", "type", "sample_type", "stage_ajcc", "stage_figo", "stage_T", "stage_N", "stage_M",
              "smoking_frequency", "tobacco_smoking_status", "years_smoked", "pack_years_smoked",
              "alcohol_intensity", "tissue",
              "barcode1", "barcode2", "barcode3") 
    
    pheno <- tcga_pheno |>
      mutate(type = case_when(sample_type %in% c("Solid Tissue Normal",
                                                 "Blood Derived Normal",
                                                 "Buccal Cell Normal",
                                                 "Bone Marrow Normal") ~ "Control",
                              sample_type %in% c("Primary Tumor",
                                                 "Primary Blood Derived Cancer - Peripheral Blood",
                                                 "Primary Blood Derived Cancer - Bone Marrow") ~ "Cancer")) |>
      rename(any_of(lookup)) |>
      dplyr::select(any_of(cols))
    
    # beta_param
    # tmp <- beta_parameters(beta)
    
    # identical(rownames(tmp), pheno$barcode1)
    pheno <- cbind(pheno, tmp)
    pheno <- pheno[!is.na(pheno$age),]
    
    if(which(projects$id == id) == 1){
      data <- pheno # if first one, pheno tmp is base
    } else {
      data <- plyr::rbind.fill(data, pheno) # append and compile all the data
    }
    
    rm(pheno, tcga_pheno, beta, dat, tmp);gc()
  }
}

save(data, file = "1-functional-association/4-tcga-data/2-output/data-interim.Rdata")

# simplified tissue and cancer descriptions (a bit unnecessarily complex but whatever...)
lookup_tissue_simplified <- c("adrenal gland" = "TCGA-ACC",
                              "bladder" = "TCGA-BLCA",
                              "breast" = "TCGA-BRCA",
                              "cervix" = "TCGA-CESC",
                              "bile duct" = "TCGA-CHOL",
                              "colon" = "TCGA-COAD",
                              "lymphoid tissues" = "TCGA-DLBC",
                              "oesophagus" = "TCGA-ESCA",
                              "brain" = "TCGA-GBM",
                              "oral region" = "TCGA-HNSC",
                              "kidney" = "TCGA-KICH",
                              "kidney" = "TCGA-KIRC",
                              "kidney" = "TCGA-KIRP",
                              "bone marrow" = "TCGA-LAML",
                              "brain" = "TCGA-LGG",
                              "liver" = "TCGA-LIHC",
                              "lung" = "TCGA-LUAD",
                              "lung" = "TCGA-LUSC",
                              "pleura" = "TCGA-MESO",
                              "ovary" = "TCGA-OV",
                              "pancreas" = "TCGA-PAAD",
                              "adrenergic tissue" = "TCGA-PCPG",
                              "prostate" = "TCGA-PRAD",
                              "rectum" = "TCGA-READ",
                              "connective tissue" = "TCGA-SARC",
                              "skin" = "TCGA-SKCM",
                              "stomach" = "TCGA-STAD",
                              "testis" = "TCGA-TGCT",
                              "thyroid" = "TCGA-THCA",
                              "thyroid" = "TCGA-THYM",
                              "endometrium" = "TCGA-UCEC",
                              "uterus" = "TCGA-UCS",
                              "uvea" = "TCGA-UVM")

lookup_cancer_type = c("adrenocortical carcinoma" = "TCGA-ACC",
                       "bladder urothelial carcinoma" = "TCGA-BLCA",
                       "breast invasive carcinoma" = "TCGA-BRCA",
                       "cervical squamous cell carcinoma and endocervical adenocarcinoma" = "TCGA-CESC",
                       "cholangiocarcinoma" = "TCGA-CHOL",
                       "lymphoid neoplasm" = "TCGA-DLBC",
                       "colon adenocarcinoma" = "TCGA-COAD",
                       "oesophageal carcinoma" = "TCGA-ESCA",
                       "glioblastoma multiforme" = "TCGA-GBM",
                       "head and neck squamous carcinoma" = "TCGA-HNSC",
                       "renal chromophobe cell carcinoma" = "TCGA-KICH",
                       "renal clear cell carcinoma" = "TCGA-KIRC",
                       "renal papillary cell carcinoma" = "TCGA-KIRP",
                       "acute myeloid leukemia" = "TCGA-LAML",
                       "lower grade glioma" = "TCGA-LGG",
                       "hepatocellular carcinoma" = "TCGA-LIHC",
                       "lung adenocarcinoma" = "TCGA-LUAD",
                       "lung squamous carcinoma" = "TCGA-LUSC",
                       "mesothelioma" = "TCGA-MESO",
                       "ovarian carcinoma" = "TCGA-OV",
                       "pancreatic adenocarcinoma" = "TCGA-PAAD",
                       "pheochromocytoma and paraganglioma" = "TCGA-PCPG",
                       "prostate adenocarcinoma" = "TCGA-PRAD",
                       "rectum adenocarcinoma" = "TCGA-READ",
                       "sarcoma" = "TCGA-SARC",
                       "skin cutaneous melanoma" = "TCGA-SKCM",
                       "stomach adenocarcinoma" = "TCGA-STAD",
                       "testicular germ cell tumour" = "TCGA-TGCT",
                       "thyroid carcinoma" = "TCGA-THCA",
                       "thymoma" = "TCGA-THYM",
                       "uterine corpus endometrial carcinoma" = "TCGA-UCEC",
                       "uterine carcinosarcoma" = "TCGA-UCS",
                       "uveal melanoma" = "TCGA-UVM")

lookup <- data.frame(project = lookup_cancer_type,
                     cancer_type = names(lookup_cancer_type),
                     tissue_simplified = names(lookup_tissue_simplified))

rownames(lookup) <- NULL
data <- data |>
  left_join(lookup)
save(data, file = "1-functional-association/4-tcga-data/2-output/data.Rdata")