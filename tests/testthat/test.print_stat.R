test_that("Just some check about print_msg", {
  
  set_verbosity(2)
  set.seed(2)
  d <- rnorm(4)

  msg_info <- capture.output(print_stat("Hello world!",  d, 2,"INFO"))
  expect_equal(msg_info, "|-- INFO :  Hello world!: Min:-1.13 Q1:-0.96 Med:-0.36 Mean:-0.06 Q3:0.54 Max:1.59 ")

  msg_debug <- capture.output(print_stat("Hello world!",  d, 2,"DEBUG"))
  expect_equal(msg_debug, "|-- DEBUG :  Hello world!: Min:-1.13 Q1:-0.96 Med:-0.36 Mean:-0.06 Q3:0.54 Max:1.59 ")

  msg_warning <- suppressWarnings(print_msg("Hello world!", "WARNING"))
  expect_equal(msg_warning, NULL)
})
