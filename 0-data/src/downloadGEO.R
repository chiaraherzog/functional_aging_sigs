downloadGEO <- function(accession,
                        directory,
                        getRaw = T,
                        regex = "_RAW.tar"){
  
  require(here)
  require(GEOquery)
  
  # make directory for the GEO data
  tmpdir <- paste0(directory, accession, "/")
  
  if(!dir.exists(tmpdir)){
    dir.create(tmpdir)
  }
  
  # download matrix file
  gse <- getGEO(accession)
  
  # extract and save pheno
  pheno <- pData(gse[[1]])
  save(pheno, file = paste0(tmpdir, "pheno.Rdata"))
  
  if(getRaw == T){
  # temporarily change wd so downloads will go in the right directory
  setwd(paste0(tmpdir))
    
  # Querying available GEO supplementary files and identify _RAW.tar
  files <- getGEOSuppFiles(accession, fetch_files = T, makeDirectory = T, filter_regex = regex)
  
  # Untar and unzip files
  untar(paste0(accession, "/", list.files(paste0(accession, "/"))))
  files <- list.files()
  files[grepl(".gz", files)]
      
  for (i in 1:length(files[grepl(".gz", files)])){
    gunzip(files[grepl(".gz", files)][i])
  } 
  
  setwd(here())
  }
}