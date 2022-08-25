library(testthat)
library(shiny)
library(chip_invasion)

test_check("chip_invasion")

# test_dir(
#   "./testthat",
#   # Run in the app's environment containing all support methods.
#   env = shiny::loadSupport(),
#   # Display the regular progress output and throw an error if any test error is found
#   reporter = c("progress", "fail")
# )

#test_dir(paste0(cwd, "/tests/testthat/"))


