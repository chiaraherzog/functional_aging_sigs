# Download and preprocess mouse methylation atlas data
library(GEOquery)

# download phenotype data from GEO
gds <- getGEO("GSE184410")
dat <- as.data.frame(gds$GSE184410_series_matrix.txt.gz)

# extract and clean relevant variables
pheno <- data.frame(matrix(NA,nrow=nrow(dat),ncol=12))
rownames(pheno) <- rownames(dat)

colnames(pheno) <- c('source',
                     'organism',
                     'site',
                     'sex',
                     'pathology',
                     'cell_line',
                     'age_days',
                     'description',
                     'tissue',
                     'sorted_cell_type',
                     'strain',
                     'basename')

pheno[,1] <- as.factor(dat$source_name_ch1)
pheno[,2] <- as.factor(dat$organism_ch1)
pheno[,3] <- as.factor(dat$site.ch1)
pheno[,4] <- as.factor(dat$Sex.ch1)
pheno[,5] <- as.factor(dat$pathology.ch1)
pheno[,6] <- as.factor(dat$cell.line.ch1)
pheno[,7] <- as.numeric(dat$age_days.ch1)
pheno[,8] <- as.factor(dat$description)

tissue <- character(nrow(dat))
for(i in 1:nrow(dat)){
  tissue[i] <- substr(dat$characteristics_ch1.5[i],
                      start=9,
                      stop=nchar(dat$characteristics_ch1.5[i]))
}
pheno[,9] <- as.factor(tissue)

sorted_cell_type <- character(nrow(dat))
for(i in 1:nrow(dat)){
  sorted_cell_type[i] <- substr(dat$characteristics_ch1.6[i],
                                start=19,
                                stop=nchar(dat$characteristics_ch1.6[i]))
}
pheno[,10] <- as.factor(sorted_cell_type)

strain <- character(nrow(dat))
for(i in 1:nrow(dat)){
  strain[i] <- substr(dat$characteristics_ch1.3[i],
                      start=9,
                      stop=nchar(dat$characteristics_ch1.3[i]))
}
pheno[,11] <- as.factor(strain)

pheno[,12] <- paste(rownames(dat),dat$description.1,sep='_')

summary(droplevels(pheno[pheno$description%in%c('Tissue Biology',
                                                'Sorted Blood Panel'),]))

save(pheno, file = "3-mouse/1-output/pheno.Rdata")


# Methylation beta file -------

library(minfi)
library(ChAMP)
library(IlluminaMouseMethylationmanifest)
library(IlluminaMouseMethylationanno.12.v1.mm10)
library(impute)

INTENSITY_THRESHOLD <- 9.5     # minimum median intensity required
DETECTION_P_THRESHOLD <- 0.01  # maximum detection p-value
FAILED_PROBE_THRESHOLD <- 0.1   # maximum proportion of failed probes per sample

# Set download directory
directory <- '<your download directory>'

# Load raw data
RGset <- read.metharray.exp(base = paste0(directory, '/GSE184410_RAW/'),
                            verbose = TRUE,
                            force = TRUE,
                            recursive = TRUE)

RGset@annotation <-  c(array = "IlluminaMouseMethylation", annotation = "12.v1.mm10")

detP <- minfi::detectionP(RGset, type = "m+u")
Mset <- preprocessRaw(RGset)
qc <- getQC(Mset)

# Background intensity correction
ssNOOB <- preprocessNoob(RGset, 
                         dyeCorr=TRUE, 
                         verbose=TRUE, 
                         dyeMethod='single')

# list low quality samples
low_intensity_samples <- rownames(qc)[qc$mMed<INTENSITY_THRESHOLD | qc$uMed<INTENSITY_THRESHOLD]
failed_samples <- colnames(detP)[colSums(detP>DETECTION_P_THRESHOLD) > (nrow(detP) * FAILED_PROBE_THRESHOLD)]

samples_to_remove <- unique(c(low_intensity_samples,
                              failed_samples))
rm(RGset);gc()

# extract beta values
beta_ssNOOB <- getBeta(ssNOOB)
rm(ssNOOB);gc()

# Remove failed probes (n=29)
beta_ssNOOB[detP > DETECTION_P_THRESHOLD] <- NA
r <- rowSums(is.na(beta_ssNOOB))
beta_ssNOOB <- beta_ssNOOB[(r < 0.8*ncol(beta_ssNOOB))>0, ]

# keep only samples needed
load("3-mouse/1-output/pheno.Rdata")

pheno <- droplevels(pheno[pheno$description%in%c('Tissue Biology',
                                                 'Sorted Blood Panel'),])
ind <- match(pheno$basename, colnames(beta_ssNOOB))
beta_ssNOOB_subset <- beta_ssNOOB[,ind]
rm(beta_ssNOOB);gc()
rm(detP);rm(Mset);gc()

# impute missing values
out <- capture.output(beta_ssNOOB_subset_imputed <- impute.knn(beta_ssNOOB_subset,k=10,rowmax=0.8)$data)
rm(beta_ssNOOB_subset);invisible(gc())

# probe type normalisation
beta_final <- IlluminaMouseMethylationmanifest::champ.normm(beta=beta_ssNOOB_subset_imputed,arraytype=array,cores=4)

# save output
saveRDS(beta_final, file= paste0(directory, '/beta_final.Rds'))
