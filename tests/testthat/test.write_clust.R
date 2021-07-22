test_that("Checking write_clust is working properly", {

  m <- matrix(rnorm(80000), nc=20)
  m[1:100,1:10] <- m[1:100,1:10] + 4
  m[101:200,11:20] <- m[101:200,11:20] + 3
  m[201:300,5:15] <- m[201:300,5:15] + -2
  res <- DBFMCL(data=m,
                distance_method="pearson",
                av_dot_prod_min = 0,
                inflation = 1.2,
                k=25,
                fdr = 10)
  write_clust(object = res,
              out_path = "/tmp",
              filename_out = "ALL.sign.txt")
  expect_equal(ncol(read.table("/tmp/ALL.sign.txt")), 21)
})