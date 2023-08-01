library(shinytest2)

test_that("{shinytest2} recording: chip_invasion", {
  app <- AppDriver$new(variant = platform_variant(), name = "chip_invasion", height = 842, 
      width = 1330)
  app$upload_file(gfp_huvec_upload = "test_data.xlsx")
  app$click("runner")
  app$set_inputs(waiter_shown = TRUE, allow_no_input_binding_ = TRUE, priority_ = "event")
  app$set_inputs(waiter_hidden = TRUE, allow_no_input_binding_ = TRUE, priority_ = "event")
  app$expect_values()
  app$click("perfield_runner")
  app$click("boxplotz")
  app$click("posplotz")
  app$expect_screenshot()
  app$click("posplotzxz")
  app$click("threedee")
  app$click("surfaceplot")
  app$set_inputs(chip_num = character(0))
  app$set_inputs(chip_num = 2)
  app$set_inputs(`plotly_hover-A` = "[{\"curveNumber\":0,\"pointNumber\":[58,27],\"x\":58,\"y\":27,\"z\":130.34932462089313}]", 
      allow_no_input_binding_ = TRUE, priority_ = "event")
  app$set_inputs(`plotly_hover-A` = character(0), allow_no_input_binding_ = TRUE, 
      priority_ = "event")
  app$expect_values()
  app$click("threedee")
  app$click("surfaceplot")
})
