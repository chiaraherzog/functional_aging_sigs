# Age Signatures - PRO/SEN/PCGT_age

PROSENPCGT_age <- function(beta){
  
  # coef sheet
  coefs <- list.files("3-evaluation/1-output", full.names = T, pattern = "*.Rds")
  coefnames <- gsub(".Rds", "", basename(coefs))
  out <- vector(mode = "list", length = length(coefs))
  names(out) <- paste0("idx_", coefnames)
  
  for (i in 1:length(out)){
    # load coefs
    wx <- readRDS(coefs[i])
    
    # compute index
    intersect <- intersect(names(wx), rownames(beta))
    w <- wx[intersect]
    
    if(length(intersect) == 1) {
      tmp <- beta[intersect,]
      B <- as.data.frame(t(tmp))
      rownames(B) <- intersect
    } else {
      B <- beta[intersect,]
    }
    
    if(!grepl("idx_PRO", coefs[i])){
      tmp <- apply(B, 2, function(tmp){
        return(sum(w*tmp)/length(w))
      })
    } else {
    tmp <- apply(B, 2, function(i){
    return(abs(1-sum(i*w)/length(w)))
    })
    }
    
    out[[i]] <- tmp
    
  }
  
  return(out)
  
}
