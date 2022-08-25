#test the functions of the chip invasion app
context("functions")

test_that("input fxn grabs all required components", {
 # testdata <- dataframe_i_will_make
  
  expect_length(endothelial_epithelial_reader(paste0(cwd, "/tests/testthat/sheets/organoid_test_data.xlsx")), 3)
  expect_length(endothelial_epithelial_reader(paste0(cwd, "/tests/testthat/sheets/hct_test_data.xlsx")), 3)
  #test that the cell_type assignment is correct -- hct vs org
}) #yeah? something like this

#expect output  - for pca fxn?
#is_a - make sure chip/field are factors (after which fxn tho?) -- after pca fxn; for counts make sure theyre numeric and char? or is it double and char
#also for pca??

# test_that("data is correct type for tables", {
#   #organoid_test_df <- epi_count()
#   bot_org_test <- readxl::read_xlsx(paste0(cwd, "/tests/testthat/sheets/filtered_test_data.xlsx"))
#   
#   expect_true(is.character(bot_org_test$chip)) #numeric? char? double? what are you
#   expect_true(is.character(bot_org_test$field))
# 
# })

# test_that("ggplots can print to ui", {
#   #ugh i need more test files -- after counts & stuff have been applied -- so that i can test plot fxns
#   expect_output(chip_3d_plot_ee(organoids_counted_df, 1))
# }) #rather that this ^, test that the plot generated can be printed? or that its able to print w/out error messages?
# #perhaps expect_true(is.ggplot(epi_count_plot(df))) <----- note, this wont work w the 3d plot (but i feel like thats ok)