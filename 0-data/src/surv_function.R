surv_function <- function(data, adj = T){
  
  source("0-data/src/surv_heatmap.R")
  # generate dataframes
  
  # Overall survival
  os <- data.frame(matrix(nrow = length(unique(data$project))*length(unique(data$index)),ncol = 4))
  colnames(os) <- c("project", "index", "hr", "p")
  os$project <- rep(unique(data$project), each = length(unique(data$index)))
  os$index <- unique(data$index)
  
  # Progression
  pf <- os
  
  # Correlation with stage
  stage <- os
  colnames(stage) <- c("project", "index", "r", "p")

  # Correlation with lymphocyte
  lympho_cor <- stage
    
  for(p in unique(data$project)){
    
    for (i in unique(data$index)){
      
      # OS 
      if(adj == T){
        fit <- coxph(Surv(os_time, os) ~ value_scale + age + ic, data = data[data$index == i & data$project == p,])
      } else {
        fit <- coxph(Surv(os_time, os) ~ value_scale, data = data[data$index == i & data$project == p,])
      }
      os[os$project==p & os$index==i,]$hr <- exp(summary(fit)$coefficients[1,1])
      os[os$project==p & os$index==i,]$p <- summary(fit)$coefficients[1,5]
      
      # PFI 
      if(adj == T){
        fit <- coxph(Surv(pfi_time, pfi) ~ value_scale + age + ic, data = data[data$index == i & data$project == p,])
      } else {
        fit <- coxph(Surv(pfi_time, pfi) ~ value_scale, data = data[data$index == i & data$project == p,])
      }
      pf[pf$project==p & pf$index==i,]$hr <- exp(summary(fit)$coefficients[1,1])
      pf[pf$project==p & pf$index==i,]$p <- summary(fit)$coefficients[1,5]
      
      # Association with stage
      if(nrow(data[data$index == i & data$project == p & !is.na(data$stage),])>0){
        corrtest <- cor.test(data[data$index == i & data$project == p,]$stage,
                             data[data$index == i & data$project == p,]$value)
        stage[stage$project==p & stage$index==i,]$r <- corrtest$estimate
        stage[stage$project==p & stage$index==i,]$p <- corrtest$p.value
      }

      # Correlation with lympho: adjust for ic first
      if (adj == T){
        fit <- lm(value ~ ic, data = data[data$index == i & data$project == p,])
        data$value_adj <- NA
        data[data$index == i & data$project == p,]$value_adj <- data[data$index == i & data$project == p,]$value - as.numeric(predict(fit))
        corrtest <- cor.test(data[data$index == i & data$project == p,]$lymphocyte,
                             data[data$index == i & data$project == p,]$value_adj)
      } else {
        corrtest <- cor.test(data[data$index == i & data$project == p,]$lymphocyte,
                             data[data$index == i & data$project == p,]$value)
      }
      lympho_cor[lympho_cor$project==p & lympho_cor$index==i,]$r <- as.numeric(corrtest$estimate)
      lympho_cor[lympho_cor$project==p & lympho_cor$index==i,]$p <- corrtest$p.value
      
    }
  }
  
  
  # Plots
  out1 <- list(os = os,
                    pf = pf,
                    stage = stage,
                    lympho_cor = lympho_cor)
  
  out2 <- surv_heatmap(out1)
  
  return(list(data = out1, plots = out2))
  
  
}
