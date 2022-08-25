context("server tests")

server <- function(input, output, session){
  count_me <- readxl::read_xlsx(paste0(cwd, "/tests/testthat/sheets/filtered_test_data.xlsx")) 
  
  output$position_plots <- renderPlot(position_chip_plots_xy(count_me))
}

shiny::testServer(server, {
  testthat::test_that("input fxn grabs all required components", {
    testthat::expect_lte(ncol(count_me), 50)
  }) 
    output$position_plots 
  }
)
