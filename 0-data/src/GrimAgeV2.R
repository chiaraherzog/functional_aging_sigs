#' Basic index function
#'
#' @param beta Beta-valued matrix with rownames labeling CpGs and columns labeling samples for which GrimAgeV2 are computed
#' @param age a numeric string with identical length to columns in beta, describing participant age.
#' @param sex a character string with identical length to rows in beta, describing participant sex. Default is 'f'
#' @return Dataframe with GrimAgeV2 values
#' @export


GrimAgeV2 <- function(beta, age, sex = 'f', sub_coefs = F){
  
  load("0-data/src/coefGrimAgeV2.rda")
  load("0-data/src/coefGrimAgeV2sub.rda")
  
  if(length(sex) == 1){
    sex <- rep(sex, ncol(beta))
  }
    
  if((!exists("age") | !exists("sex"))){
    stop("No age or sex supplied.")
  }
  
  if(length(age) != ncol(beta) | length(sex) !=  ncol(beta)){
    stop("Length of age or sex variable is not identical to the observations. Please check the data.")
  }
  
  # Append sex and age to beta
  beta <- rbind(beta,
                Age = age,
                Intercept = 1)
  
  # Subcoefficients
  sub <- lapply(unique(coefGrimAgeV2$Y.pred)[!unique(coefGrimAgeV2$Y.pred) %in% c("COX", 'transform')], function(s){
    
    # Get coefficients and intercept
    coef <- coefGrimAgeV2[coefGrimAgeV2$Y.pred==s,]
    
    intersect <- intersect(rownames(beta), coef$var)
    cat("[", s, ": ", length(intersect), "/", nrow(coef), " (", round(length(intersect)/nrow(coef)*100, 1), "%) of CpGs in beta matrix]\n", sep = "")
    coef <- coef[match(intersect, coef$var),]
    B <- beta[match(intersect, rownames(beta)),]
    
    apply(B, 2, function(x){sum(x*coef$beta)})
    
  })
  
  
  names(sub) <- unique(coefGrimAgeV2$Y.pred)[!unique(coefGrimAgeV2$Y.pred) %in% c("COX", 'transform')]
  sub <- as.data.frame(t(as.data.frame(sub)))
  
  # Append age + Sex again for cox model
  sub <- rbind(sub,
               Age = age,
               Female = ifelse(sex == 'f', 1, 0))
  
  cox <- apply(sub, 2, function(x){
    sum(x*coefGrimAgeV2[coefGrimAgeV2$Y.pred=='COX',]$beta)
  })
  
  # Transformation
  m_age <- coefGrimAgeV2[coefGrimAgeV2$Y.pred=='transform' & coefGrimAgeV2$var == 'm_age',]$beta
  sd_age <- coefGrimAgeV2[coefGrimAgeV2$Y.pred=='transform' & coefGrimAgeV2$var == 'sd_age',]$beta
  m_cox <- coefGrimAgeV2[coefGrimAgeV2$Y.pred=='transform' & coefGrimAgeV2$var == 'm_cox',]$beta
  sd_cox <- coefGrimAgeV2[coefGrimAgeV2$Y.pred=='transform' & coefGrimAgeV2$var == 'sd_cox',]$beta
  
  Y <- (cox - m_cox)/(sd_cox)
  DNAmGrimAgeV2 <- (Y * sd_age) + m_age
  AgeAccelGrimV2 <- DNAmGrimAgeV2-sub['Age',]
  
  if(sub_coefs == F){
    out <- data.frame(DNAmGrimAgeV2 = DNAmGrimAgeV2,
                      AgeAccelGrimV2 = as.numeric(AgeAccelGrimV2))
  } else {
    out <- data.frame(DNAmGrimAgeV2 = DNAmGrimAgeV2,
                      AgeAccelGrimV2 = as.numeric(AgeAccelGrimV2),
                      t(sub))
    
    out <- out |>
      dplyr::select(-c("Age", "Female"))
  }
  return(out)
  
}
