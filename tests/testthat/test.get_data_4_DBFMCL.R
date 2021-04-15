test_that("Cheking get_data_4_DBFMCL returns a DBFMCLresult object", {
  m <- matrix(rnorm(80000), nc = 20)
  res <- get_data_4_DBFMCL(data=m)
  expect_equal(inherits(res, "list"), TRUE)
  expect_equal(inherits(res[[1]], "matrix"), TRUE)
  expect_equal(ncol(res[[1]]), 20)
  expect_equal(nrow(res[[1]]), 80000/20)
  expect_equal(is.null(res[[2]]), TRUE)
})