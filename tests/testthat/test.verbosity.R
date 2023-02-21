test_that("Checking set_verbosity() and get_verbosity()", {
  #==============================
  # Verbosity = 0
  set_verbosity(0)
  expect_equal(get_verbosity(), 0)
  
  msg_info <- capture.output(print_msg("Hello world!", "INFO"))
  expect_equal(msg_info, as.character())
  
  msg_debug <- capture.output(print_msg("Hello world!", "DEBUG"))
  expect_equal(msg_debug, as.character())
  
  msg_warning <- suppressWarnings(print_msg("Hello world!", "WARNING"))
  expect_equal(msg_warning, "|-- WARNING : Hello world!")
  
  
  #==============================
  # Verbosity = 1
  set_verbosity(1)
  expect_equal(get_verbosity(), 1)
  
  msg_info <- capture.output(print_msg("Hello world!", "INFO"))
  expect_equal(msg_info, "|-- INFO :  Hello world! ")
  
  msg_debug <- capture.output(print_msg("Hello world!", "DEBUG"))
  expect_equal(msg_debug, as.character())
  
  msg_warning <- suppressWarnings(print_msg("Hello world!", "WARNING"))
  expect_equal(msg_warning, "|-- WARNING : Hello world!")
  
  
  #==============================
  # Verbosity = 2
  set_verbosity(2)
  expect_equal(get_verbosity(), 2)
  
  msg_info <- capture.output(print_msg("Hello world!", "INFO"))
  expect_equal(msg_info, "|-- INFO :  Hello world! ")
  
  msg_debug <- capture.output(print_msg("Hello world!", "DEBUG"))
  expect_equal(msg_debug, "|-- DEBUG :  Hello world! ")
  
  msg_warning <- suppressWarnings(print_msg("Hello world!", "WARNING"))
  expect_equal(msg_warning, "|-- WARNING : Hello world!")
  
})
