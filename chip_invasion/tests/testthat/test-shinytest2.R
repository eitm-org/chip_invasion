library(shinytest2)
test_that("{shinytest2} recording: chip_invasion_test", {
  app <- AppDriver$new(name = "chip_invasion_test", height = 844, width = 1256)
  app$upload_file(gfp_huvec_upload = "test_endo.xlsx")
  app$expect_values()
  app$click("runner")
  app$set_inputs(waiter_shown = TRUE, allow_no_input_binding_ = TRUE, priority_ = "event")
  app$set_inputs(waiter_hidden = TRUE, allow_no_input_binding_ = TRUE, priority_ = "event")
  app$click("perfield_runner")
  app$click("boxplotz")
  app$click("posplotz")
  app$click("posplotzxz")
  app$click("threedee")
  app$set_inputs(`plotly_hover-A` = "[{\"curveNumber\":0,\"pointNumber\":0,\"x\":753.71,\"y\":2185.6,\"z\":130.452}]", 
      allow_no_input_binding_ = TRUE, priority_ = "event")
  app$set_inputs(`plotly_hover-A` = character(0), allow_no_input_binding_ = TRUE, 
      priority_ = "event")
  app$set_inputs(`plotly_hover-A` = "[{\"curveNumber\":0,\"pointNumber\":6,\"x\":897.16,\"y\":5228.99,\"z\":121.2}]", 
      allow_no_input_binding_ = TRUE, priority_ = "event")
  app$set_inputs(`plotly_hover-A` = "[{\"curveNumber\":0,\"pointNumber\":4,\"x\":101.64,\"y\":-3768.63,\"z\":-84.547}]", 
      allow_no_input_binding_ = TRUE, priority_ = "event")
  app$set_inputs(`plotly_hover-A` = character(0), allow_no_input_binding_ = TRUE, 
      priority_ = "event")
  app$click("surfaceplot")
  app$expect_values()
  app$set_inputs(chip_num = character(0))
  app$set_inputs(chip_num = 2)
  app$expect_values()
  app$set_inputs(`plotly_hover-A` = "[{\"curveNumber\":0,\"pointNumber\":[11,91],\"x\":11,\"y\":91,\"z\":135.6152361088585}]", 
      allow_no_input_binding_ = TRUE, priority_ = "event")
  app$set_inputs(`plotly_hover-A` = character(0), allow_no_input_binding_ = TRUE, 
      priority_ = "event")
  app$click("threedee")
  app$click("surfaceplot")
})
