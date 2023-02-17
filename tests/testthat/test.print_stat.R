test_that("Just some check about print_msg", {
  options(warn=-1)
  set.seed(123)
  data <- rnorm(3)
  res <- print_stat("Summary statistics for data", data, round_val = 2, msg_type = "WARNING")
  expect_equal(res, "|-- WARNING : Summary statistics for data: Min:-0.56 Q1:-0.4 Med:-0.23 Mean:0.26 Q3:0.66 Max:1.56")
  data <- rnorm(100)
  res <- print_stat("Summary statistics for data", data, round_val = 1, msg_type = "WARNING")
  
  data <- letters[1:10]
  res <- print_stat("Summary statistics for data", data, round_val = 1, msg_type = "WARNING")
  expect_equal(res, "|-- WARNING : Summary statistics for data: No Statistics")
  
})
