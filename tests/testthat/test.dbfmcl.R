test_that("Cheking DBFMCL is providing the right number of genes", {
  set.seed(123)
  m <- matrix(rnorm(80000), nc=20)
  m[1:100,1:10] <- m[1:100,1:10] + 4
  m[101:200,11:20] <- m[101:200,11:20] + 3
  m[201:300,5:15] <- m[201:300,5:15] + -2
  res <- DBFMCL(data=m,
                distance_method="pearson",
                k=25)
  #expect_equal(length(res@size), 1)
  #expect_equal(res@size, 309)
})