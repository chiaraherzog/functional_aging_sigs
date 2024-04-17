# Computation of PhenoAge
# Data/paper by Levine et al., 2018 (Aging)

PhenoAge <- function(beta){
  
  #----------------#
  # load coefficients
  dat <- read.table("~/Dropbox/src/data/phenoage.csv", sep = ",",
                    header = T)
  intercept <- dat$Weight[1]
  w <- dat$Weight[-1]
  names(w) <- dat$CpG[-1]
  rm(dat)
  
  #------------------------------------------#
  # compute index
  intersect <- intersect(rownames(beta), names(w))
  B <- beta[intersect,]
  w <- w[intersect]
  
  B1 <- B*w
  PhenoAge <- intercept + apply(B1, MARGIN = 2, FUN = 'sum')
  
  return(PhenoAge)
}