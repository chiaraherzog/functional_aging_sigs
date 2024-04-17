# Adjust values by ic
# Indices can be freely defined

adjust_ic <- function(data,
                      indices = NULLtissue){
  
  # initiate columns
  for (i in indices){
    data[[paste0(i, "_adj")]] <- numeric(nrow(data))
    tmp <- data[data$type=="Control",]
    fit <- lm(tmp[[i]] ~ ic, data = tmp)
    data[,paste0(i, "_adj")] <- data[[i]] - as.numeric(predict(fit, newdata = data))
  }
  
  return(data)
}
