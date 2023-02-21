test_that("Just some check about print_msg", {
  msg_info <- capture.output(print_msg("Hello world!", "INFO"))
  expect_equal(msg_info, "|-- INFO :  Hello world! ")

  msg_debug <- capture.output(print_msg("Hello world!", "DEBUG"))
  expect_equal(msg_debug, "|-- DEBUG :  Hello world! ")

  msg_warning <- suppressWarnings(print_msg("Hello world!", "WARNING"))
  expect_equal(msg_warning, "|-- WARNING : Hello world!")
})
