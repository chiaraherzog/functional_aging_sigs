
# Heatmap

surv_heatmap <- function(data){
  
  library(ComplexHeatmap)
  
  scores <- c("PCGTpure", "PCGTgen", "PCGTepi", "PCGTimm",
              "SENpure", "SENgen", "AGE",
              "PROgen_neg", "PROpure_neg",
              "SENepi", "PROepi_neg",
              "SENimm", "PROimm_neg", "AGEimm",
              "SENnonage", "PROnonage", "PROepi_pos",
              "AGEepi", "PCGTnonage",
              "PROimm_pos", "PROgen_pos", "PROpure_pos",
              "hannum", "horvath", "phenoage")
  
  out <- vector(mode = "list", length = 4)
  names(out) <- names(data)
  
  for (i in 1:length(data)) {
    
    if(names(data)[i] %in% c("os", "pf")){
      tmp <- data[[i]] |> 
        tidyr::pivot_wider(id_cols = c(project),
                           names_from = 'index',
                           values_from = hr) |> 
        mutate(across(any_of(scores), ~ log(.x)))
      
      # create matrix for heatmap
      mat <- tmp[,scores]
      
      # set rownames
      rownames(mat) <- tmp$project
      
      # Create pval matrix
      tmp <- data[[i]] |> 
        tidyr::pivot_wider(id_cols = c(project),
                           names_from = 'index',
                           values_from = p) |> 
        dplyr::mutate(across(any_of(scores), ~ ifelse(.x < 0.05, "*", "")))
      matp <- tmp[,scores]
      rownames(matp) <- tmp$project
      
      out[[i]] <- Heatmap(as.matrix(mat),
              cluster_rows = F,
              cluster_columns = F,
              show_row_names = T,
              show_column_names = T,
              cell_fun = function(j, i, x, y, width, height, fill) {
                grid.text(matp[i, j], x, y, gp = gpar(fontsize = 10))
              },
              col = circlize::colorRamp2(breaks = seq(-2, 2, by = 0.5),
                                         colors = rev(brewer.pal(name = "RdBu", n = 9))),
              name = "log(HR)",
              row_names_gp = grid::gpar(fontsize = 8),
              column_names_gp = grid::gpar(fontsize = 9),
              show_row_dend = F,
              border = T)
      
    
    } else {
      tmp <- data[[i]] |> 
        tidyr::pivot_wider(id_cols = c(project),
                           names_from = 'index',
                           values_from = r) |> 
        tidyr::drop_na()
      
      # create matrix for heatmap
      mat <- tmp[,scores]
      
      # set rownames
      rownames(mat) <- tmp$project
      
      # Create pval matrix
      tmp <- data[[i]] |> 
        tidyr::pivot_wider(id_cols = c(project),
                           names_from = 'index',
                           values_from = p) |> 
        dplyr::mutate(across(any_of(scores), ~ ifelse(.x < 0.05, "*", ""))) |> 
        tidyr::drop_na()
      matp <- tmp[,scores]
      rownames(matp) <- tmp$project
      
      out[[i]] <- Heatmap(as.matrix(mat),
                          cluster_rows = F,
                          cluster_columns = F,
                          show_row_names = T,
                          show_column_names = T,
                          cell_fun = function(j, i, x, y, width, height, fill) {
                            grid.text(matp[i, j], x, y, gp = gpar(fontsize = 10))
                          },
                          col = circlize::colorRamp2(breaks = seq(-1, 1, by = 0.25),
                                                     colors = rev(brewer.pal(name = "RdBu", n = 9))),
                          name = "cor",
                          row_names_gp = grid::gpar(fontsize = 8),
                          column_names_gp = grid::gpar(fontsize = 9),
                          show_row_dend = F,
                          border = T)
    }
  }
  
  return(out)
}
