test_that("create_rand_str returns a string of length 10", {
  rand_str <- create_rand_str()
  expect_equal(nchar(rand_str), 10)
})

test_that("create_rand_str only contains letters and digits", {
  rand_str <- create_rand_str()
  expect_match(rand_str, "^[[:alnum:]]+$")
})

test_that("Check the result using seed", {
  set.seed(123)
  rand_str <- create_rand_str()
  expect_match(rand_str, "1eSN95tOk2")
})