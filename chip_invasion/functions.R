#functions, revised to use fastTps to build endothelial surface boundary
#2023-07-13

library(tidyverse)
library(here)
library(janitor)
library(readxl)
library(hablar)
library(ggbeeswarm)
library(plotly)
library(FactoMineR)
library(factoextra)
library(fields)
cwd <- here::here()

#function to read appropriate files into space
endothelial_epithelial_reader <- function(data_path) {
  require(readxl)
  require(dplyr)
  huvec_xl <- read_xlsx(file.path(data_path), sheet = c("Endothelial"), skip = 9) %>%
    mutate(Chip = as.double(str_extract(Compound, pattern = "(?<=chip)\\d*|(?<=Chip)\\s\\d*|(?<=Chip)\\d*")))
  
  huvec_xl <- huvec_xl %>% janitor::clean_names() %>%
    rename_with(~str_remove(.x, "_mm")) %>%
    rename_with(~str_remove(.x, "x546_upper_spot_bright_546_spot_bright_image_region_|x546_upper_spot_bright_")) %>%
    rename_with(~str_remove(.x, "\\d$"))
  
  gfp_xl <- read_xlsx(file.path(data_path), sheet = c("Epithelial"), skip = 9) %>%
    mutate(Chip = as.double(str_extract(Compound, pattern = "(?<=chip)\\d*|(?<=Chip)\\s\\d*|(?<=Chip)\\d*"))) 
  
  #cell_type <- "organoid"
  if("488 SPSB Nuclei Selected - 488 SPSB Nucleus Centroid Z [µm]" %in% colnames(gfp_xl)){
    cell_type <- "hct"
  } else {cell_type <- "organoid"}
  
  gfp_xl <- gfp_xl %>% janitor::clean_names() %>%
    rename_with(~str_remove(.x, "_mm|_mm\\d")) %>%
    rename_with(~str_remove(.x, "x488_spsb_nuclei_selected_488_spsb_nucleus_|x488_spsb_nuclei_selected_|x488_spot_bright_image_region_selected_488_spot_bright_image_region_|x488_spot_bright_image_region_selected_|x488_spot_bright_image_region_488_spot_bright_image_region_")) %>%
    rename_with(~str_remove_all(.x, "488_spot_bright_image_region_alexa_|nucleus_alexa_488_|x488_spot_bright_image_region_")) %>%
    rename_with(~str_remove(.x, "\\d$"))
  
  df_list <- list("huvec_df" = huvec_xl, "organoid_df" = gfp_xl, "cell_type" = cell_type)
  return(df_list)
}

#plot diagnostic histograms
diagnostic_plots_ee <- function(huvec_df, organoid_df) {
  huvec_plot <- ggplot(huvec_df, aes(x = centroid_z)) +
    geom_histogram() +
    facet_wrap(~ chip) +
    labs(title = "Endothelial cells at each z position")
  
  organoid_plot <- ggplot(organoid_df, aes(x = centroid_z))  +
    geom_histogram() +
    facet_wrap(~ chip) +
    labs(title = "Epithelial cells at each z position")
  
  print(huvec_plot)
  print(organoid_plot)
}

#define the z-boundary using endothelial cells
endo_boundzer <- function(huvec_df, organoid_df) {
  #browser()
  endo_chip_list <- split(huvec_df, f=huvec_df$chip)
  epi_chip_list <- split(organoid_df, f=organoid_df$chip)
  
  org_positions <- lapply(epi_chip_list, function(x) dplyr::select(x, position_x, position_y))
  
  endo_fits <- lapply(endo_chip_list, function(x) fastTps(x[, c("position_x", "position_y")], Y=x$centroid_z, aRange = 200)) 
  
  se_endo_fits <- lapply(endo_fits, function(x) predictSE(x))
  avg_se_fits <- lapply(se_endo_fits, function(x) mean(x)*qnorm(0.975)) #top half of the ci around the fitted surface
  ses_together <- bind_rows(avg_se_fits) %>%
    pivot_longer(cols = 1:6, names_to = "chip", values_to = "se") %>%
    mutate(chip = as.double(chip))
  
  pred_fits <- map2(.x=endo_fits, .y=org_positions, ~predict(.x, .y))
  org_fits <- map2(.x=epi_chip_list, .y=pred_fits, ~cbind(.x, .y))
  org_fits_clean <- lapply(org_fits, function(x) rename(x, "gfp_pred"=.y))
  
  ## split off here in some way
  org_fits_joined <- bind_rows(org_fits_clean) %>% 
    left_join(ses_together, by="chip") #%>%
  # group_by(chip) %>%
  # mutate(below_pred = centroid_z < (gfp_pred+se)) %>%
  # filter(below_pred == TRUE)
  
  #todo: make a summry table of zheights? idk how that would work here
  
  stuff_to_return <- list("huvec_list" = endo_chip_list, 
                          "org_list" = epi_chip_list, 
                          "fit_list" = endo_fits, 
                          "org_fits_df" = org_fits_joined)
  
  return(stuff_to_return)
}

#count cells in the bottom channel
epi_count <- function(org_fits_df) {
  
  filtered_gfp <- org_fits_df %>%
    group_by(chip) %>%
    mutate(boundary = gfp_pred+se,
           in_bottom_channel = centroid_z < boundary,
           in_top_channel = centroid_z > boundary) 
  
  # full_join(organoid_df, field_chip_huvec) %>% 
  # mutate(in_bottom_channel = centroid_z < mean_huvec_offset) %>% 
  # mutate(in_top_channel = centroid_z > mean_huvec_offset)
  
  gfp_in_bottom_channel <- filtered_gfp %>%
    filter(in_bottom_channel == TRUE) %>%
    convert(fct(field, chip))
  
  gfp_in_top_channel <- filtered_gfp %>%
    filter(in_top_channel == TRUE) %>%
    convert(fct(field, chip))
  
  gfp_per_chip_bottom <- gfp_in_bottom_channel %>%
    group_by(compound,chip) %>% 
    summarise(gfp_cells_in_bottom = n()) %>%
    arrange(chip)
  
  gfp_per_chip_top <- gfp_in_top_channel %>%
    group_by(compound, chip) %>% 
    summarise(gfp_cells_in_top = n()) %>%
    arrange(chip)
  
  gfp_per_field_bottom <- gfp_in_bottom_channel %>%
    group_by(chip, field, 
             .drop = FALSE
    ) %>% 
    summarise(gfp_cells_in_bottom = n()) %>% 
    ungroup() %>%
    arrange(chip)
  
  gfp_per_field_top <- gfp_in_top_channel %>%
    group_by(compound, chip, field) %>% 
    summarise(gfp_cells_in_top = n()) %>% 
    ungroup() %>%
    arrange(chip)
  
  gfp_count_list <- list("all_gfp_bottom" = gfp_in_bottom_channel, "all_gfp_top" = gfp_in_top_channel,
                         "gfp_per_chip_bottom" = gfp_per_chip_bottom, "gfp_per_chip_top" = gfp_per_chip_top,
                         "gfp_per_field_bottom" = gfp_per_field_bottom, "gfp_per_field_top" = gfp_per_field_top)
  
  return(gfp_count_list)
}

#make a table of cells in each field. group by top of the bottom channel and bottom of the bottom channel
epi_pf_count <- function(bottom_gfp_df, gfp_per_field_bottom){
  min_z_overall <- min(bottom_gfp_df$boundary)
  
  middle_filtered_gfp <- bottom_gfp_df %>%
    #group_by(compound, chip, field) %>%
    mutate(middle_bound = ((min_z_overall + mean(boundary))/2)) %>%
    mutate(in_top_half = centroid_z > middle_bound) %>%
    mutate(in_bottom_half = centroid_z < middle_bound)
  
  gfp_tophalf <- middle_filtered_gfp %>%
    filter(in_top_half == TRUE)
  
  gfp_bottomhalf <- middle_filtered_gfp %>%
    filter(in_bottom_half == TRUE)
  
  gfp_per_field_tophalf <- gfp_tophalf %>%
    group_by(chip, field, 
             # .drop = FALSE
    ) %>%
    summarise(gfp_cells_tophalf = n()) %>%
    ungroup() %>%
    arrange(chip)
  
  gfp_per_field_bottomhalf <- gfp_bottomhalf %>%
    group_by(chip, field, 
             #.drop = FALSE
    ) %>%
    summarise(gfp_cells_bottomhalf = n()) %>%
    ungroup() %>%
    arrange(chip)
  
  full_pf <- full_join(gfp_per_field_tophalf, gfp_per_field_bottomhalf) #%>% distinct(compound, chip, field)
  total_pf <- inner_join(gfp_per_field_bottom, full_pf)
  
  total_pf <- total_pf %>%
    mutate_all(~replace(., is.na(.), 0)) %>% 
    arrange(chip) %>%
    dplyr::rename("Chip" = chip, "Field" = field, "Total Epithelial Cells in the Bottom Channel" = gfp_cells_in_bottom, 
                  "Top of the Bottom" = gfp_cells_tophalf, "Bottom of the Bottom" = gfp_cells_bottomhalf)
  
  return(total_pf)
}

#plot the 3d endothelial surface
cute_surface_maker <- function(huvec_list, fit_list, chip) {
  
  x_seq <- seq(0, max(huvec_list[[chip]]$position_x), len = 100)
  y_seq <- seq(0, max(huvec_list[[chip]]$position_y), len = 100)  
  
  pos_grid <- expand.grid(x=x_seq, y=y_seq)
  
  fit_matrix <- matrix(predict(fit_list[[chip]], pos_grid), 100, 100)
  
  plot_ly(z=fit_matrix, colors = "plasma") %>% 
    add_surface()
  
  #return(cute_surface)
  
}

#create a 3d plot of the bottom channel
chip_3d_plot_ee <- function(bottom_gfp_df, chip_in) {
  chip_used <- bottom_gfp_df %>% filter(chip == chip_in)
  
  pos_x <- chip_used$position_x
  pos_y <- chip_used$position_y
  pos_z <- chip_used$centroid_z
  
  plot_ly(x = pos_x, y = pos_y, z = pos_z, type="scatter3d", mode="markers")
}

#boxplot of the number of cells in each field
epi_count_plot <- function(bottom_gfp_df) {
  
  plot_obj <- bottom_gfp_df %>%
    convert(fct(field, chip)) %>% 
    group_by(chip, field) %>%
    summarise(gfp_cells_in_bottom = n())
  
  ggplot(plot_obj, aes( x = chip, y = gfp_cells_in_bottom, label = field)) +
    geom_boxplot(outlier.colour = 'red') +
    geom_beeswarm() +
    geom_label(position = "jitter") +
    labs(title = "Epithelial cell counts in the bottom channel", 
         subtitle = "each point is a field")
}

#scatterplot of where each cell is on the x-y plane
position_chip_plots_xy <- function(organoid_df) {
  xyplot <- ggplot(organoid_df, aes(x = position_x, y = position_y)) + geom_point() + facet_wrap(~ chip) + labs(title = "XY Position Plot")
  
  return(xyplot)
}

#scatterplot of where each cell is on the x-z axis
position_chip_plots_xz <- function(organoid_df) {
  xzplot <- ggplot(organoid_df, aes(x=position_x, y=centroid_z)) + geom_point() + facet_wrap(~chip) + labs(title = "XZ Position Plot")
  
  return(xzplot)
}

#nicer-looking per-chip counts
per_chip_counter_ee <- function(gfp_in_bottom_channel) {
  
  gfp_per_chip <- gfp_in_bottom_channel %>%
    group_by(chip) %>% 
    summarise(gfp_cells_in_bottom = n()) %>%
    arrange(chip) %>%
    dplyr::rename("Chip" = chip, "Epithelial Cells in the Bottom Channel" = gfp_cells_in_bottom)
  
  return(gfp_per_chip)
}

#inital per-field counts
per_field_counter_ee <- function(gfp_in_bottom_channel) {
  
  gfp_per_field <- gfp_in_bottom_channel %>%
    group_by(chip, field, .drop = FALSE) %>% 
    summarise(gfp_cells_in_bottom = n()) %>% 
    ungroup() %>%
    #distinct(Chip, Field, .keep_all = TRUE) %>%
    arrange(chip)
  
  return(gfp_per_field)
}