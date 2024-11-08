# ComputeAUC function
# Adjusting for ic, age, or both; 
# Paired samples (yes/no)
# Reference group is by default "Control", with type being the variable name.
# If multiple comparisons, "compgroups" can be specified (all will be compared to control reference.)

computeAUC <- function(data, 
                       refgroup = "Control",
                       refvar = "type",
                       celltype = "celltype",
                       group = "group",
                       adj.ic = T,
                       adj.age = T,
                       paired = F,
                       compgroups = NULL){
  
  # library
  library(pROC)
  library(rstatix)
  
  # indices
  indices <- c("idx_AGE", "idx_AGEepi", "idx_AGEimm",
               "idx_PCGTepi", "idx_PCGTimm", "idx_PCGTgen", "idx_PCGTpure", "idx_PCGTnonage",
               "idx_PROepi_pos", "idx_PROimm_pos", "idx_PROgen_pos", "idx_PROpure_pos",
               "idx_PROepi_neg", "idx_PROimm_neg", "idx_PROgen_neg", "idx_PROpure_neg",
               "idx_PROnonage",
               "idx_SENepi", "idx_SENimm", "idx_SENgen", "idx_SENpure", "idx_SENnonage",
               "hannum", "horvath", "phenoage",
               "DNAmGrimAgeV2", "CausAge", "DamAge", "AdaptAge")
  
  # helper functions
  source("0-data/src/adjust_age.R")
  source("0-data/src/adjust_ic.R")
  source("0-data/src/adjust_age_ic.R")
  
  # Format type
  if (refvar != "type"){
    data <- data |> 
      dplyr::mutate(type = get(UQ(refvar)))
  }
  
  # Rename controlgroup if not control
  if(refgroup != "Control"){
    data <- data |> 
      dplyr::mutate(type = ifelse(type == refgroup, "Control", as.character(type)))
  }
  
  # Relevel factors
  if(is.null(compgroups)){
    compgroups = as.character(unique(data$type))
    compgroups <- compgroups[compgroups != "Control"]
  } else {
    compgroups = compgroups
  }
  
  data <- data |> 
    dplyr::mutate(type = factor(type, levels = c("Control", compgroups)))
  
  # adjust ic and/or age with helper functions (or not, depending on model specification)
  if(adj.ic == T & adj.age == T){ 
    data <- data |>
      adjust_age_ic(indices = indices)
    
    # pivot long
    long <- data |> 
      tidyr::pivot_longer(all_of(paste0(indices, "_adj")),
                          names_to = "index",
                          values_to = "value") |> 
      droplevels() 
  } else if (adj.ic == T & adj.age == F){
    data <- data |>
      adjust_ic(indices = indices)
    
    # pivot long
    long <- data |> 
      tidyr::pivot_longer(all_of(paste0(indices, "_adj")),
                          names_to = "index",
                          values_to = "value") |> 
      droplevels() 
  } else if(adj.ic == F & adj.age == T){
    data <- data |>
      adjust_age(indices = indices)
    
    # pivot long
    long <- data |> 
      tidyr::pivot_longer(all_of(paste0(indices, "_adj")),
                          names_to = "index",
                          values_to = "value") |> 
      droplevels() 
  } else {
    long <- data |> 
      tidyr::pivot_longer(all_of(indices),
                          names_to = "index",
                          values_to = "value") |> 
      droplevels() 
  }
  
  # compute cohen's d
  if(paired == T){
    cohend <- long |> 
      dplyr::group_by(index) |>
      dplyr::arrange(id) |> 
      rstatix::cohens_d(value ~ type, ref.group = "Control", paired = T)
  } else {
    cohend <- long |> 
      dplyr::group_by(index) |>
      rstatix::cohens_d(value ~ type, ref.group = "Control", paired = F)
  }

  
  if(adj.ic == T | adj.age == T){
    indices = paste0(indices, "_adj")
  }
  
  for (c in unique(compgroups)){
    cat("Group ", c, "\n")
    # Create output DF
    tmpdf <- data.frame(matrix(nrow = length(indices),
                               ncol = 11))
    colnames(tmpdf) <- c("type", "index", "celltype", "group", "p", "d", "auc", "cilo", "cihi", "n_control", "n_case")
    
    # Compute P
    p <- long |> 
      dplyr::filter(type %in% c("Control", c)) |> 
      droplevels() |> 
      dplyr::group_by(index) |> 
      dplyr::do(tidy(wilcox.test(
        .$value[.$type == "Control"], 
        .$value[.$type == c], 
        paired = paired)))
    
    
    for (i in 1:length(indices)){
      
      tmpdf$type <- c
      tmpdf$celltype <- unique(data[[celltype]])
      tmpdf$group <- unique(data[[group]])
      tmpdf$index[i] <- paste0(indices[i])
      tmpdf$p[i] <- p[p$index == indices[i],]$p.value
      tmpdf$d[i] <- -as.numeric(cohend[cohend$group2 == c & cohend$index == indices[i],]$effsize)
      auc <- ci(auc(droplevels(data[data$type %in% c("Control", c),]$type), 
                    data[data$type %in% c("Control", c),][[indices[i]]],
                    direction = "<"))
      tmpdf$auc[i] <- auc[2]
      tmpdf$cilo[i] <- auc[1]
      tmpdf$cihi[i] <- auc[3]
      tmpdf$n_control[i] <- nrow(data[data$type == "Control",])
      tmpdf$n_case[i] <- nrow(data[data$type == c,])
    }
    
    if(c == unique(compgroups)[1]){
      out <- tmpdf
    } else {
      out <- rbind(out, tmpdf)
    }
    
  }
  
  return(out)
  
}
