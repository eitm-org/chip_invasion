#functions, revised to include analysis of HCTs
library(tidyverse)
library(here)
library(janitor)
library(readxl)
library(hablar)
library(ggbeeswarm)
library(plotly)
library(FactoMineR)
library(factoextra)
cwd <- here::here()

#function to read appropriate files into space
endothelial_epithelial_reader <- function(data_path) {
  require(readxl)
  require(dplyr)
  huvec_xl <- read_xlsx(file.path(data_path), sheet = c("Endothelial"), skip = 9) %>%
    mutate(Chip = as.double(str_extract(Compound, pattern = "(?<=chip)\\d*|(?<=Chip)\\s\\d*")))
  
  huvec_xl <- huvec_xl %>% janitor::clean_names() %>%
    rename_with(~str_remove(.x, "_mm")) %>%
    rename_with(~str_remove(.x, "x546_upper_spot_bright_546_spot_bright_image_region_|x546_upper_spot_bright_")) %>%
    rename_with(~str_remove(.x, "\\d$"))
  
  gfp_xl <- read_xlsx(file.path(data_path), sheet = c("Epithelial"), skip = 9) %>%
    mutate(Chip = as.double(str_extract(Compound, pattern = "(?<=chip)\\d*|(?<=Chip)\\s\\d*"))) 
  
  #cell_type <- "organoid"
  if("488 SPSB Nuclei Selected - 488 SPSB Nucleus Centroid Z [µm]" %in% colnames(gfp_xl)){
    cell_type <- "hct"
  } else {cell_type <- "organoid"}
  
  gfp_xl <- gfp_xl %>% janitor::clean_names() %>%
    rename_with(~str_remove(.x, "_mm|_mm\\d")) %>%
    rename_with(~str_remove(.x, "x488_spsb_nuclei_selected_488_spsb_nucleus_|x488_spsb_nuclei_selected_|x488_spot_bright_image_region_selected_488_spot_bright_image_region_|x488_spot_bright_image_region_selected_")) %>%
    rename_with(~str_remove(.x, "488_spot_bright_image_region_alexa_|nucleus_alexa_488_")) %>%
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

#define the z-boundary using endothelial cells. different for organoids vs hcts
endo_boundary <- function(huvec_df, cell_type) {
  
  mean_huvec_pop <- huvec_df %>% 
    # filter(`546 Upper Spot Bright - 546 Spot Bright Image Region Centroid Z [¬µm]` >= channel_boundary_number )  %>% 
    group_by(chip) %>% 
    summarise(mean_top_bottom = mean(centroid_z), 
              median_top_bottom = median(centroid_z),
              sd_top_bottom = sd(centroid_z)
    ) %>% 
    mutate(lower_cuttoff = mean_top_bottom - (2 * sd_top_bottom),
           upper_cuttoff = mean_top_bottom + (2 * sd_top_bottom))
  
  middle_chip_huvec <- huvec_df %>% 
    inner_join(mean_huvec_pop) %>%
    group_by(chip, field) %>%
    mutate(top_bottom_2sd_ll = centroid_z >= lower_cuttoff,
           top_bottom_2sd_ul = centroid_z <= upper_cuttoff)
  
  if(cell_type == "organoid"){
    middle_chip_huvec <- middle_chip_huvec %>%
      filter(
        # Field %in% c(4, 5, 6, 7, 8, 9, 10, 11, 12),
        top_bottom_2sd_ll == TRUE, 
        top_bottom_2sd_ul == TRUE)
  }
  
  field_chip_huvec <- middle_chip_huvec %>%
    dplyr::summarise(mean_huvec_z_per_field = mean(centroid_z), 
                     sd_huvec_z_per_field = sd(centroid_z),
                     huvec_count = n()) %>% 
    mutate(mean_huvec_offset = case_when(cell_type == "hct" ~ (mean_huvec_z_per_field + (0.7*(mean_huvec_z_per_field / sd_huvec_z_per_field))),
                                         cell_type == "organoid" ~ (mean_huvec_z_per_field - (0.25*sd_huvec_z_per_field)))) %>%
    mutate(mean_huvec_offset_old = (mean_huvec_z_per_field + (1.1*(mean_huvec_z_per_field / sd_huvec_z_per_field)))) %>%
    ungroup()
  
  field_chip_huvec_summary <- field_chip_huvec %>% 
    dplyr::select(chip, field, mean_huvec_z_per_field, mean_huvec_offset) %>%
    mutate(mean_huvec_z_per_field = round(mean_huvec_z_per_field, 2),
           mean_huvec_offset = round(mean_huvec_offset, 2)) %>%
    dplyr::rename("Chip" = chip, "Field" = field, "Mean Endothelial Z-height" = mean_huvec_z_per_field, "Endothelial Cutoff Z-height" = mean_huvec_offset)
  
  huvec_dfs <- list("field_chip_huvec" = field_chip_huvec, "field_chip_huvec_summary" = field_chip_huvec_summary)
  
  return(huvec_dfs)
}

#count cells in the bottom channel
epi_count <- function(organoid_df, field_chip_huvec) {
  
  filtered_gfp <- full_join(organoid_df, field_chip_huvec) %>% 
    mutate(in_bottom_channel = centroid_z < mean_huvec_offset) %>% 
    mutate(in_top_channel = centroid_z > mean_huvec_offset)
  
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
  min_z_overall <- min(bottom_gfp_df$centroid_z)
  
  middle_filtered_gfp <- bottom_gfp_df %>%
    #group_by(compound, chip, field) %>%
    mutate(middle_bound = (min_z_overall + mean_huvec_offset)/2) %>%
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

#filter out low-intensity cells (organoids only). collect intensity metrics
org_intensity_filtering <- function(gfp_in_bottom_channel) {
  
  gfp_high_intensity <- gfp_in_bottom_channel %>%
    group_by(chip) %>%
    mutate(median_int = median(intensity_488_mean)) %>% 
    ungroup() %>%
    mutate(above_avg_int = (intensity_488_mean > median_int)) %>%
    filter(above_avg_int == TRUE)
  
  gfp_intensity_metrics <- gfp_high_intensity %>%
    select(chip, median_int) %>%
    mutate(median_int = round(median_int, 2)) %>%
    dplyr::rename("Chip" = chip, "Median Intensity" = median_int) %>%
    distinct()
  
  intensity_dfs <- list("gfp_high_intensity" = gfp_high_intensity, "gfp_intensity_metrics" = gfp_intensity_metrics)
  
  return(intensity_dfs)
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

#select only columns needed for pca
ee_clean_for_pca <- function(intensity_filtered_df){
  clean_df <- intensity_filtered_df %>%
    select(-(5:8), -(11:17), -(19:27), -(37:39), -(44:46), -(48:52)) %>% #select relevant cols
    mutate(chip = as_factor(chip), field = as_factor(field))   #chip and field become factors

  return(clean_df)
}

#compute pca on intensity-filtered dataset
ee_invasion_pcar <- function(clean_df){
  nointensity.pca <- clean_df %>%
    mutate_if(is.numeric, scale) %>%
    PCA(quanti.sup = c(1:3, 5:7, 19, 21), quali.sup = c(4, 20), graph = FALSE)

  return(nointensity.pca)
}

#graph pca individuals plot to check for subgroups
invasion_pcar_graphs <- function(nointensity.pca, clean_df){
  indiv_plot <- fviz_pca_ind(nointensity.pca, label="none", habillage = as.factor(clean_df$chip))
  
  contrib_plot <- fviz_pca_var(nointensity.pca, col.var = "contrib") + 
    scale_color_gradient2(low="blue", mid="white", 
                          high="red", midpoint=5)+theme_bw()
  
  pca_out_plots <- list("individuals_plot" = indiv_plot, "variable_contribution_plot" = contrib_plot)
  
  return(pca_out_plots)
} 

#turn the pca into a density plot of Dim1 positions. optimize the density function to find a local minimum
#assign data as either cell or non-cell based on this Dim1 boundary
pcar_density_filter <- function(clean_df, nointensity.pca) {
  pca_density_df <- clean_df %>%
    mutate(dim1_coord = nointensity.pca$ind$coord[,1], 
           dim2_coord = nointensity.pca$ind$coord[,2]) 
  
  d <- density(pca_density_df$dim1_coord)
  global_max_density <- d$x[which.max(abs(diff(d$y)))]
  
  #set the boundary for cell v no-cell determination using the local min of the density plot
  local_min_density <- optimize(approxfun(d$x,d$y),interval=c(-6,global_max_density))$minimum
  
  #assign cell type based on the dim1 coord density minimum
  pca_density_df_cnc <- pca_density_df %>%
    mutate(cell_type = case_when(dim1_coord < local_min_density ~ "non_cell", 
                                 dim1_coord > local_min_density ~ "cell")) %>%
    filter(cell_type == "cell")
  
  density_output <- list("local_min_density" = local_min_density, "pca_filtered_df" = pca_density_df_cnc)
  
  return(density_output)
}