# Adjust values by age and ic

adjust_age_ic <- function(data,
                          indices = NULL){
  
  # initiate columns
  for (i in indices){
    data[[paste0(i, "_adj")]] <- numeric(nrow(data))
    tmp <- data[data$type=="Control",]
    fit <- lm(tmp[[i]] ~ age + ic, data = tmp)
    data[,paste0(i, "_adj")] <- data[[i]] - as.numeric(predict(fit, newdata = data))
  }
  
  return(data)
}