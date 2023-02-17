test_that("Just some check about print_msg", {
  res <- print_msg("Hello world!", "INFO")
  expect_equal(res, NULL)
  res <- print_msg("Hello world!", "DEBUG")
  expect_equal(res, NULL)
  res <- print_msg("Hello world!", "WARNING")
  expect_equal(res, "|-- WARNING : Hello world!")  
})