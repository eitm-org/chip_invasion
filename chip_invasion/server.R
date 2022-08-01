library(shiny)
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

source('functions.R')

shinyServer(function(input, output){
  options(shiny.maxRequestSize=100*1024^2)
  
  filedata <- reactive({
    infile <- input$gfp_huvec_upload
    if (is.null(infile)) {
      # User has not uploaded a file yet
      return(NULL)
    }
    endothelial_epithelial_reader(infile$datapath)
  })
  
  huvec_boundz <- reactive({
    celltype <- filedata() $cell_type
    huvecs <- filedata() $huvec_df
    if (is.null(huvecs)) {
      return(NULL)}
    endo_boundary(huvecs, celltype)
  })
  
  gfp_counts <- reactive({
    organoids <- filedata() $organoid_df
    huvec_border <- huvec_boundz() $field_chip_huvec
    if (is.null(organoids)) {
      return(NULL)}
    epi_count(organoids, huvec_border)
  })
  
  gfp_filterz <- reactive({
    req(input$analysis %in% c("Intensity filter", "PCA filter"))
    gfp_bottom <- gfp_counts() $all_gfp_bottom
    if (is.null(gfp_bottom)) {
      return(NULL)}
    org_intensity_filtering(gfp_bottom)
  })
  
  cleanerz <- reactive({
    req(input$analysis == "PCA filter")
    semi_filtered_gfps <- gfp_filterz() $gfp_high_intensity
    if (is.null(semi_filtered_gfps)) {
      return(NULL)}
    ee_clean_for_pca(semi_filtered_gfps)
  })
  
  pca_doer <- eventReactive(input$button2, {
    req(input$analysis == "PCA filter")
    cleaned_gfp_df <- cleanerz()
    if (is.null(cleaned_gfp_df)) {
      return(NULL)}
    ee_invasion_pcar(cleaned_gfp_df)
  })
  
  pca_density_filterz <- reactive({
    cleaned_gfp_df <- cleanerz()
    pca_output <- pca_doer()
    if (is.null(pca_output)) {
      return(NULL)}
    pcar_density_filter(cleaned_gfp_df, pca_output)
  }) 
  
  pca_plotz <- reactive({
    gfps_to_plot <- cleanerz()
    pca_output <- pca_doer()
    if (is.null(pca_output)) {
      return(NULL)}
    invasion_pcar_graphs(pca_output, gfps_to_plot)
  })
  
  chip_countz <- eventReactive(input$runner, {
    count_type <- input$analysis
    if(count_type == "Endothelial boundary only"){
      count_me <- gfp_counts() $all_gfp_bottom
    } else if(count_type == "Intensity filter"){
      count_me <- gfp_filterz() $gfp_high_intensity
    } else if(count_type == "PCA filter"){
      count_me <- pca_density_filterz() $pca_filtered_df
    }
    per_chip_counter_ee(count_me)
  })
  
  field_countz <- eventReactive(input$perfield_runner, {
    count_type <- input$analysis
    # pf_bottom <- gfp_counts() $gfp_per_field_bottom
    if(count_type == "Endothelial boundary only"){
      count_me <- gfp_counts() $all_gfp_bottom
    } else if(count_type == "Intensity filter"){
      count_me <- gfp_filterz() $gfp_high_intensity
    } else if(count_type == "PCA filter"){
      count_me <- pca_density_filterz() $pca_filtered_df
    }
    pfbot <- per_field_counter_ee(count_me)
    epi_pf_count(count_me, pfbot)
  })
  
  position_plotz <- eventReactive(input$posplotz, {
    count_type <- input$analysis
    if(count_type == "Endothelial boundary only"){
      count_me <- gfp_counts() $all_gfp_bottom
    } else if(count_type == "Intensity filter"){
      count_me <- gfp_filterz() $gfp_high_intensity
    } else if(count_type == "PCA filter"){
      count_me <- pca_density_filterz() $pca_filtered_df
    }
    position_chip_plots_xy(count_me)
  })
  
  position_plotzxz <- eventReactive(input$posplotzxz, {
    count_type <- input$analysis
    if(count_type == "Endothelial boundary only"){
      count_me <- gfp_counts() $all_gfp_bottom
    } else if(count_type == "Intensity filter"){
      count_me <- gfp_filterz() $gfp_high_intensity
    } else if(count_type == "PCA filter"){
      count_me <- pca_density_filterz() $pca_filtered_df
    }
    position_chip_plots_xz(count_me)
  })
  
  pf_count_plotz <- eventReactive(input$boxplotz, {
    count_type <- input$analysis
    if(count_type == "Endothelial boundary only"){
      count_me <- gfp_counts() $all_gfp_bottom
    } else if(count_type == "Intensity filter"){
      count_me <- gfp_filterz() $gfp_high_intensity
    } else if(count_type == "PCA filter"){
      count_me <- pca_density_filterz() $pca_filtered_df
    }
    epi_count_plot(count_me)
  })
  
  chip_3d_plotterz <- eventReactive(input$threedee, {
    chip_number <- input$chip_num
    count_type <- input$analysis
    if(count_type == "Endothelial boundary only"){
      count_me <- gfp_counts() $all_gfp_bottom
    } else if(count_type == "Intensity filter"){
      count_me <- gfp_filterz() $gfp_high_intensity
    } else if(count_type == "PCA filter"){
      count_me <- pca_density_filterz() $pca_filtered_df
    }
    chip_3d_plot_ee(count_me, chip_number)
  })
  
  observeEvent(req(input$analysis == "PCA filter"), {
    output$button2 <- renderUI({
      actionButton("button2", label = "PCA check", icon=icon("check-double"), class="btn btn-primary")
    })
  })
  
  output$huvec_z_border <- renderDataTable(huvec_boundz() $field_chip_huvec_summary) #x
  output$intensity_filt_valz <- renderDataTable(gfp_filterz() $gfp_intensity_metrics) #x
  output$pca_filt_val <- renderText({ 
    paste0("Objects less than ", pca_density_filterz() $local_min_density, " on the x-axis (Dim1) are non-cells and are removed.")
  })
  
  output$pca_plot <- renderPlot(pca_plotz() $individuals_plot)
  
  output$analysis_table_all <- renderDataTable(chip_countz()) #x
  output$field_table <- renderDataTable(field_countz()) #x
  
  output$position_plots <- renderPlot(position_plotz())
  output$position_plotsxz <- renderPlot(position_plotzxz())
  output$per_field_boxplots <- renderPlot(pf_count_plotz())
  output$plot_but_3d <- renderPlotly(chip_3d_plotterz())
})
