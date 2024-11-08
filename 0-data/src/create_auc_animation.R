library(ggplot2)
library(gganimate)
library(dplyr)
library(tidyr)
library(magick)

# Function to create AUC animation with customizable index and output filename
# and higher-quality output
create_auc_animation <- function(data_path, index_col, idx_label, output_filename, width = 800, height = 800, fps = 20,
                                 col = cols[12]) {
  # Load the data
  load(data_path)
  
  # Prepare the data
  dat <- data |> 
    dplyr::mutate(type = factor(type, levels = c("Control", "Cancer"))) |> 
    dplyr::filter(project == 'TCGA-BRCA')
  
  # Calculate the ROC AUC
  auc1 <- roc(dat$type, dat[[index_col]], direction = "<")
  pos <- data.frame(x = 'positive',
                    threshold = auc1$thresholds,
                    sens = auc1$sensitivities,
                    spec = auc1$specificities) |> 
    dplyr::filter(threshold != -Inf & threshold != Inf) |> 
    dplyr::arrange(dplyr::desc(threshold))
  
  # Prepare boxplot data by cross-joining thresholds
  box_dat <- dat |> dplyr::cross_join(dplyr::select(pos, threshold)) |> 
    dplyr::arrange(dplyr::desc(threshold))
  
  # Create the boxplot animation
  a <- box_dat |>
    ggplot(aes(x = type,
               y = .data[[index_col]])) +
    geom_half_boxplot(aes(fill = type),
                      alpha = 0.2,
                      outlier.shape = NA) +
    theme_bw() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.y = element_markdown(),
          plot.title = element_markdown()) +
    labs(x = '', y = paste0('<b style="color:', col, '">', idx_label, '</b>')) +
    scale_colour_manual(values = cols[c(9, 1)],
                        aesthetics = c("colour", "fill")) +
    geom_hline(aes(yintercept = threshold), linetype = "dashed", colour = col) +
    transition_reveal(along = threshold) +
    labs(title = 'Threshold: {frame_along}')
  
  # Create the ROC curve animation
  b <- pos |>
    ggplot(aes(x = 1 - spec,
               y = sens,
               colour = x)) +
    geom_path(size = 0.8) +
    geom_abline(slope = 1, intercept = 0,
                colour = 'grey60',
                linetype = 'dashed') +
    theme_bw() +
    theme(panel.grid = element_blank(),
          legend.position = 'none') + 
    labs(x = '1-Specificity',
         y = 'Sensitivity') +
    scale_colour_manual(values = col,
                        name = '') +
    transition_reveal(along = threshold)
  
  # Animate both plots with higher quality settings
  a_gif <- animate(a, width = width, height = height, fps = fps)
  b_gif <- animate(b, width = width, height = height, fps = fps)
  
  # Combine the animations using magick
  new_gif <- image_append(c(a_gif[1], b_gif[1]))
  for (i in 2:min(length(a_gif), length(b_gif))) {
    combined <- image_append(c(a_gif[i], b_gif[i]))
    new_gif <- c(new_gif, combined)
  }
  
  # Save the final animation
  anim_save(output_filename, animation = new_gif)
}

# Example usage with higher resolution and smoother fps
# create_auc_animation("3-evaluation/3-output/data.Rdata", "idx_PROpure_neg", "AUCanimation_highres.gif", width = 800, height = 800, fps = 30)
