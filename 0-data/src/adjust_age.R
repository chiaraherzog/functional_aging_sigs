# Adjust values by age
# Indices can be freely defined

adjust_age <- function(data,
                      indices = NULLtissue){
  
  # initiate columns
  for (i in indices){
    data[[paste0(i, "_adj")]] <- numeric(nrow(data))
    tmp <- data[data$type=="Control",]
    fit <- lm(tmp[[i]] ~ age, data = tmp)
    data[,paste0(i, "_adj")] <- data[[i]] - as.numeric(predict(fit, newdata = data))
  }
  
  return(data)
}
