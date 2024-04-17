# Proliferation blocking experiments
# GSE197512
dir <- "~/Documents/Work/data.nosync/GEO/"

# Download ------------------
# Issue with GEOquery/file (?) - manual download of txt and tar files for matrix file (pheno)

# set destination for download: matrix
dir <- "<your download directory>"
destfile <- paste0(dir, "/GSE197512/GSE197512_series_matrix.txt.gz")
url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE197nnn/GSE197512/matrix/GSE197512_series_matrix.txt.gz"
# Execute download
download.file(url, destfile)

# set destination for download: raw data
destfile <- paste0(dir, "/GSE197512/raw/GSE197512_RAW.tar")
url <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE197512&format=file" # raw data
# Execute download
download.file(url, destfile)

# unzip + tar
untar(paste0(dir, "/GSE197512/raw/GSE197512_RAW.tar"))

f <- list.files(paste0(dir, "/GSE197512/raw/"), full.names = T)
for (i in f){
  gunzip(i)
}

# Preprocess pheno ------------------
destfile <- paste0(dir, "/GSE197512/GSE197512_series_matrix.txt.gz")
gunzip(destfile)

# parse with GEOquery
fname <- gsub(".gz", "", destfile)
gse <- GEOquery::parseGEO(fname)
.na_strings = c('NA','null','NULL','Null')

# Read in the Beta file (code extracted from GEOquery parseGEO)
txt = readr:::read_lines(fname)
tbl_begin = grep('!\\w+_table_begin',txt,perl=TRUE)
if(length(tbl_begin>0)) {
  txt = txt[1:tbl_begin[1]]
  dat3 <- readr::read_tsv(fname, comment='!dataset_table_end', skip = tbl_begin[1],
                          guess_max = 10000, na = .na_strings)
}

tx <- txt[grep("^!Sample", txt)]
pheno <- data.frame(matrix(ncol = sum(grepl("^!Sample", txt)),
                           nrow = 372))

tmp <- strsplit(tx, "\t")
x <- do.call(cbind.data.frame, tmp)
colnames(x) <- gsub("!", "", x[1,])
x <- x[-1,]
colnames(x)[11:24] <- paste0(substr(colnames(x)[11:24], 1, nchar(colnames(x)[11:24])-1),
                             2:15)
colnames(x)[45] <- paste0(colnames(x[45]), 2)
colnames(x)
pheno <- x |> 
  dplyr::mutate_all(~ gsub("\"", "", .))

save(pheno, file = paste0(dir, "/GSE197512/pheno.Rdata"))
rm(dat3, tmp, x, fname, tbl_begin, tx, txt)
gc()

# prep pheno for qc
pheno_simplified <- pheno |> 
  dplyr::mutate(basename = gsub("_Grn.idat.gz", "", basename(Sample_supplementary_file))) |> 
  dplyr::mutate(across(contains("Sample_characteristics"), ~ stringr::str_trim(gsub(".*:","", .)))) |> 
  dplyr::rename(id = Sample_title,
                tissue = Sample_source_name_ch1,
                subexperiment = Sample_characteristics_ch1,
                coriell_id = Sample_characteristics_ch2,
                donor_age = Sample_characteristics_ch4,
                donor_sex = Sample_characteristics_ch5,
                culture_oxygen = Sample_characteristics_ch6,
                culture_agent = Sample_characteristics_ch7,
                culture_pctfbs = Sample_characteristics_ch8,
                expression_vector = Sample_characteristics_ch9,
                subculture = Sample_characteristics_ch10,
                pd = Sample_characteristics_ch11,
                days_in_culture = Sample_characteristics_ch12,
                batch_date = Sample_characteristics_ch13) |> 
  dplyr::select(basename, id,
                tissue, subexperiment, coriell_id,
                donor_age:batch_date)
pheno <- pheno_simplified    
save(pheno, file = "1-functional-association/2-proliferation/1-output/pheno.Rdata") 


