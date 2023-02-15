test_that("Cheking DBF is providing the right number of genes", {
  #Create matrix
  set.seed(123)
  m <- matrix(rnorm(80000), nc=20)
  m[1:100, 1:10] <- m[1:100, 1:10] + 3
  m[101:200, 11:20] <- m[101:200, 11:20] + 3
  m[301:400, 5:15] <- m[201:300, 5:15] - 3
  m[201:300, c(1,5,10,15,20)] <- m[201:300, c(1,5,10,15,20)] - 3

  res <- DBF(data = m,
             distance_method = "pearson",
             k = 75,
             row_sum=-Inf,
             highest=0.3,
             fdr = 1e-8)
  #Test results
  expect_equal(res@size, 438)
  #Remove output files
  file.remove("exprs.dbf_out.txt")
})

