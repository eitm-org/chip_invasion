library(shiny)
library(bslib)
source('functions.R')

shinyUI(fluidPage(
  theme = bs_theme(bootswatch = "lumen", bg = "#FDFDFD", fg = "#4E0B37", primary = "#FF81D3"),
  
  fluidPage(
    titlePanel(
      "CRC on Chip Invasion Image Analysis"
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
              Endothelial cells (such as HUVECS or HIMECS) can be used as a marker for the border since they line the bottom channel.
              When the inclusion criteria in the imaging software was broadened, low-intensity noise and dead cells could be included in the counts.
              Intensity and PCA distribution values can be used to filter out these populations."),
      #tags$br(),
      tags$b("Objective"),
      tags$p("To use the mean z height of the endothelial cells locally, on a field by field basis, to more accurately count the number of GFP expressing epithelial cells in the bottom channel.
             Intensity and cell v non-cell assignment based on PCA Dim1 density are used to refine the counts of epithelial cells.")
    )
  ),
  
  sidebarLayout(
    sidebarPanel(
      style="color: #000; background-color: #FFE8F7; border-color: #8B1462",
      fileInput("gfp_huvec_upload", "Upload Endothelial and Epithelial data as separate sheets in the same xlsx file:"),
      selectInput("analysis", label = "Select Analysis Type", choices = c("Endothelial boundary only", "Intensity filter", "PCA filter")),
      tags$hr(style="color:#8B1462;"),
      uiOutput("button2"), #button for PCA check, reacts to PCA filter from dropdown menu. details near bottom of server
      tags$hr(style="color:#8B1462;"),
      tags$p("Cell Counts"),
      actionButton("runner", "Per-chip counts", icon = icon("dna"), class="btn btn-primary"),
      actionButton("perfield_runner", "Per-field counts", icon = icon("flask"), class="btn btn-primary"),
      tags$hr(style="color:#8B1462;"),
      tags$p("Plots"),
      actionButton("boxplotz", "Per-field boxplot", icon = icon("vial"), class="btn btn-primary"),
      actionButton("posplotz", "XY position plot", icon = icon("microscope"), class="btn btn-primary"),
      actionButton("posplotzxz", "XZ position plot", icon = icon("disease"), class="btn btn-primary"),
      actionButton("threedee", "3D plot", icon = icon("cube"), class="btn btn-primary"), 
      tags$br(),
      numericInput("chip_num", "Enter chip number to plot in 3D:", 1)
      
    ), 
    
    mainPanel(
      tabsetPanel(
        tabPanel("Tables", fixedRow(column(4, dataTableOutput("analysis_table_all")),
                                    column(8, dataTableOutput("field_table"))
        )),
        tabPanel("Plots",  fluidRow(plotOutput("pca_plot")),
                 fluidRow(plotOutput("per_field_boxplots")),
                 fluidRow(plotOutput("position_plots"),
                          plotOutput("position_plotsxz")),
                 fluidRow(plotlyOutput("plot_but_3d"))
        ),
        tabPanel("Metrics", fixedRow(column(2, textOutput("pca_filt_val")),
                                     column(6, dataTableOutput("huvec_z_border")),
                                     column(4, dataTableOutput("intensity_filt_valz"))))
      )
    )
  )
)
)