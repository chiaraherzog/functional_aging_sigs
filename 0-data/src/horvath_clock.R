horvath_clock <- function(beta){
  
  file <- read.csv(file="~/Dropbox/src/data/horvath_clock_raw/gb-2013-14-10-r115-S3.csv")
  intercept <- file$CoefficientTraining[1]
  
  file <- file[-1,]
  w <- file$CoefficientTraining
  
  ind <- na.omit(match(file$CpGmarker, rownames(beta)))
  
  beta <- beta[ind,]
  ind <- match(rownames(beta), file$CpGmarker)
  w <- w[ind]
  
  B1 <- beta * w
  index <- intercept + colSums(B1)
  
  anti.trafo = function(x,adult.age=20) {
    ifelse(x<0, (1+adult.age)*exp(x)-1, (1+adult.age)*x+adult.age)
    }
  index <- anti.trafo(index)
  
  return(index)
  
}
