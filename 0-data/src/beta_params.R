beta_params <- function(beta){
  
  require(EpiDISH)
  source("0-data/src/PROSENPCGT_age.R")
  source("0-data/src/hannum_clock.R")
  source("0-data/src/horvath_clock.R")
  source("0-data/src/PhenoAge.R")
  
  # Compute IC
  out.l <- epidish(beta.m = beta,
                   ref.m = centEpiFibIC.m,
                   method = "RPC")$estF
  
  ic <- out.l[,3]
  
  # Compute hEpidish
  # hepidish
  frac.m <- hepidish(beta.m = beta,
                     ref1.m = centEpiFibIC.m,
                     ref2.m = centBloodSub.m,
                     h.CT.idx = 3,
                     method = 'RPC',
                     maxit = 500)
  
  # Compute PROSENPCGT_age
  prosenpcgt <- PROSENPCGT_age(beta)
  hannum <- hannum_clock(beta)$index_primary
  horvath <- horvath_clock(beta)
  phenoage <- PhenoAge(beta)
  
  # Merge
  tmp <- cbind(as.data.frame(prosenpcgt),
               hannum,
               horvath,
               phenoage,
               ic,
               frac.m)
  
  return(tmp)
  
}