test_that("Just some check about print_msg", {
  set_verbosity(2)
  
  expect_output(print_msg("Hello world!", "INFO"), 
                          "|-- INFO :  Hello world! ")

  expect_output(print_msg("Hello world!", "DEBUG"),
                          "|-- DEBUG :  Hello world! ")

  expect_warning(print_msg("Hello world!", "WARNING"), "|-- WARNING : Hello world!")
  expect_output(print_msg("Hello world!"), "|-- INFO :  Hello world! ")
})
