library(shiny)
library(bslib)
library(waiter)
source('functions.R')

shinyUI(fluidPage(
  theme = bs_theme(bootswatch = "lumen", bg = "#FDFDFD", fg = "#4E0B37", primary = "#FF81D3"),
  waiter::use_waiter(),
  
  fluidPage(
    titlePanel(
      "(CRC on a) Chip Invasion and Contour Analysis (ChICA)"
    ),
    tags$div(
      tags$hr(style="color:#8B1462;"),
      tags$b("Background"),
      # tags$br(),
      tags$p("The CRC on a chip model is a contained environment that can be used to study invasion of cells from the top channel to the bottom channel. 
              To study the amount of invasion from the top to the bottom channel, the border between them needs to be defined and anything below that can be counted as the bottom channel.
              Additional noise must be filtered out."),
      #tags$br(),
      tags$b("Problem"),
      tags$p("The chips don't always sit level on the mount and one z height can't be taken across the chip to use as a border to delineate top and bottom. 
              Cells that are on the bottom of the top channel can me mistakenly counted as being in the bottom channel. 
              Endothelial cells (such as HUVECS or HIMECS) can be used as a marker for the border since they line the bottom channel."),
      #tags$br(),
      tags$b("Objective"),
      tags$p("To model the endothelial cell layer as a topographical surface. This surface is then used to delineate invaded and invading objects from non-invaded objects.
              To count the number of invaded cells and provide descriptive and interactive data visualizations."),
      tags$b("Requirements"),
      tags$p("Data must be preprocessed using an image analysis tool. 
             Input files should contain the x, y, and z coordinates of the centroid of the endothelial objects and epithelial objects.
             Input files should cotain the Endothelial and Epithelial data as separate sheets in the same xlsx file.
             Example data below can be used as a template."),
      tags$a(href="test_data.xlsx", "Download Example", download=NA, target="_blank"),
     # downloadButton("downloadbutton", "Download Report")
    )
  ),
  
  sidebarLayout(
    sidebarPanel(
      style="color: #000; background-color: #FFE8F7; border-color: #8B1462",
      fileInput("gfp_huvec_upload", "Upload Endothelial and Epithelial data as separate sheets in the same xlsx file:"),
      tags$hr(style="color:#8B1462;"),
      tags$p("Cell Counts - Click 'Per-chip counts' to begin analysis. This can take up to several minutes."),
      actionButton("runner", "Per-chip counts", icon = icon("dna"), class="btn btn-primary"),
      actionButton("perfield_runner", "Per-field counts", icon = icon("flask"), class="btn btn-primary"),
      tags$hr(style="color:#8B1462;"),
      tags$p("Plots"),
      actionButton("boxplotz", "Per-field boxplot", icon = icon("vial"), class="btn btn-primary"),
      actionButton("posplotz", "XY position plot", icon = icon("microscope"), class="btn btn-primary"),
      actionButton("posplotzxz", "XZ position plot", icon = icon("disease"), class="btn btn-primary"),
      actionButton("threedee", "3D plot", icon = icon("cube"), class="btn btn-primary"), 
      actionButton("surfaceplot", "Surface plot", icon = icon("cube"), class="btn btn-primary"), 
      tags$br(),
      numericInput("chip_num", "Enter chip number to plot GFP objects or HUVEC surface in 3D:", 1)
      
    ), 
    
    mainPanel(
      tabsetPanel(
        tabPanel("Tables", fixedRow(column(4, dataTableOutput("analysis_table_all")),
                                    column(8, dataTableOutput("field_table"))
        )),
        tabPanel("Plots",  
                 fluidRow(plotlyOutput("surf_plot")),
                 fluidRow(plotlyOutput("plot_but_3d")),
                 fluidRow(plotOutput("per_field_boxplots")),
                 fluidRow(plotOutput("position_plots"),
                          plotOutput("position_plotsxz"))
        )
      )
    )
  )
)
)