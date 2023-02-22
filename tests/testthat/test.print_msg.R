test_that("Just some check about print_msg", {
  set_verbosity(2)
  
  msg_info <- capture.output(print_msg("Hello world!", "INFO"))
  expect_equal(msg_info, "|-- INFO :  Hello world! ")

  msg_debug <- capture.output(print_msg("Hello world!", "DEBUG"))
  expect_equal(msg_debug, "|-- DEBUG :  Hello world! ")

  expect_warning(print_msg("Hello world!", "WARNING"),
               "|-- WARNING : Hello world!")
})
