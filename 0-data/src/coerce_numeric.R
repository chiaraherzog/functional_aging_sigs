coerce_numeric <- function(beta_merged){
    beta_merged_tmp <- matrix(NA, nrow=nrow(beta_merged), ncol=ncol(beta_merged))
    for(i in 1:ncol(beta_merged)){
      beta_merged_tmp[,i] <- as.numeric(beta_merged[,i])
    }
    colnames(beta_merged_tmp) <- colnames(beta_merged)
    rownames(beta_merged_tmp) <- rownames(beta_merged)
    beta_merged <- beta_merged_tmp
    rm(beta_merged_tmp);gc()
    return(beta_merged)
}
