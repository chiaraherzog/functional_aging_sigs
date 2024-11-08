#' Basic index function
#'
#' @param beta Beta-valued matrix with rownames labeling CpGs and columns labeling samples for which the CausAge, DamAge, and AdaptAge clocks should be computed
#' @return Dataframe with CausAge, DamAge, and AdaptAge clock values
#' @export


CausalClocks <- function(beta){
  
  load("0-data/src/coefAdaptAge.rda")
  load("0-data/src/coefCausAge.rda")
  load("0-data/src/coefDamAge.rda")
  
    
  # Compute clocks
  clocks <- lapply(c("DamAge", "AdaptAge", "CausAge"), function(c){
    
    # Get coefficients and intercept
    coef <- get(paste0('coef', c))
    intercept <- coef$estimate[1]
    coef <- coef[-1,]
    
    # find overlap with beta matrix and subset
    intersect <- intersect(rownames(beta), coef$term)
    coef <- coef[match(intersect, coef$term),]
    B <- beta[match(intersect, rownames(beta)),]
    cat("[", c, ": ", length(intersect), "/", nrow(coef), " (", round(length(intersect)/nrow(coef)*100, 1), "%) of CpGs in beta matrix]\n", sep = "")
    
    # Compute
    apply(B, 2, function(x){intercept + sum(x*coef$estimate)})
    
  })
  
  # Name and prepare df output
  names(clocks) <- c("DamAge", "AdaptAge", "CausAge")
  out <- as.data.frame(clocks)
  
  return(out)
  
}
