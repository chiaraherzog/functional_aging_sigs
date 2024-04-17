hannum_clock <- function(beta){

  file_primary <- read.csv("~/Dropbox/src/data/hannum_clock_raw/mmc2_primary.csv")
  file_all <- read.csv("~/Dropbox/src/data/hannum_clock_raw/mmc2_all.csv")
  file_br <- read.csv("~/Dropbox/src/data/hannum_clock_raw/mmc3_breast.csv")
  
  # Primary index
  ind <- na.omit(match(file_primary$Marker, rownames(beta)))
  B <- beta[ind,]
  
  ind <- match(rownames(B), file_primary$Marker)
  w1 <- file_primary$Coefficient[ind]
  
  B1 <- B * w1
  index1 <- colSums(B1)
  
  # ModelAll index
  ind <- na.omit(match(file_all$Marker, rownames(beta)))
  B <- beta[ind,]
  
  ind <- match(rownames(B), file_all$Marker)
  w2 <- file_all$Coefficient[ind]
  
  B2 <- B * w2
  index2 <- colSums(B2)
  
  # calculating BR index 
  ind <- na.omit(match(file_br$Marker, rownames(beta)))
  B <- beta[ind,]
  
  ind <- match(rownames(B), file_br$Marker)
  w3 <- file_br$Coefficient[ind]
  
  B3 <- B * w3
  index3 <- colSums(B3) 
  
  return(list(index_primary=index1,
              w_primary=w1,
              index_all=index2,
              w_all=w2,
              index_br=index3,
              w_br=w3))
}
