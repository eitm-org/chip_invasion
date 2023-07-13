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
library(fields)
library(waiter)
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
    waiter_show(html = tagList(spin_flower(), h4("Analyzing data...")), color = "#4E0B37")
    on.exit({waiter_hide()})
    huvecs <- filedata() $huvec_df
    orgs <- filedata() $organoid_df
    if (is.null(huvecs)) {
      return(NULL)}
    endo_boundzer(huvecs, orgs)
  })
  
  gfp_counts <- reactive({
    organoids <- huvec_boundz() $org_fits_df
    if (is.null(organoids)) {
      return(NULL)}
    epi_count(organoids)
  })
  
  chip_countz <- eventReactive(input$runner, {
    count_me <- gfp_counts() $all_gfp_bottom
    if (is.null(count_me)) {
      return(NULL)}
    per_chip_counter_ee(count_me)
  })
  
  field_countz <- eventReactive(input$perfield_runner, {
    count_me <- gfp_counts() $all_gfp_bottom
    if (is.null(count_me)) {
      return(NULL)}
    pfbot <- per_field_counter_ee(count_me)
    epi_pf_count(count_me, pfbot)
  })
  
  position_plotz <- eventReactive(input$posplotz, {
    count_me <- gfp_counts() $all_gfp_bottom
    if (is.null(count_me)) {
      return(NULL)}
    position_chip_plots_xy(count_me)
  })
  
  position_plotzxz <- eventReactive(input$posplotzxz, {
    count_me <- gfp_counts() $all_gfp_bottom
    if (is.null(count_me)) {
      return(NULL)}
    position_chip_plots_xz(count_me)
  })
  
  pf_count_plotz <- eventReactive(input$boxplotz, {
    count_me <- gfp_counts() $all_gfp_bottom
    if (is.null(count_me)) {
      return(NULL)}
    epi_count_plot(count_me)
  })
  
  chip_3d_plotterz <- eventReactive(input$threedee, {
    chip_number <- input$chip_num
    count_me <- gfp_counts() $all_gfp_bottom
    if (is.null(count_me)) {
      return(NULL)}
    chip_3d_plot_ee(count_me, chip_number)
  })
  
  cute_surfacez <- eventReactive(input$surfaceplot, {
    chippy <- input$chip_num
    huvz <- huvec_boundz() $huvec_list
    fitz <- huvec_boundz() $fit_list
    if (is.null(fitz)) {
      return(NULL)}
    cute_surface_maker(huvz, fitz, chippy)
  })
  
  #output$huvec_z_border <- renderDataTable(huvec_boundz() $field_chip_huvec_summary) #make a surface-compatible version?
  
  output$analysis_table_all <- renderDataTable(chip_countz()) #x
  output$field_table <- renderDataTable(field_countz()) #x
  
  output$position_plots <- renderPlot(position_plotz())
  output$position_plotsxz <- renderPlot(position_plotzxz())
  output$per_field_boxplots <- renderPlot(pf_count_plotz())
  
  output$plot_but_3d <- renderPlotly(chip_3d_plotterz())
  output$surf_plot <- renderPlotly(cute_surfacez())
})
