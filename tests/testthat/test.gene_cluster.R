test_that("Checking gene_cluster() is working...", { 
  
  # Set verbosity to 0
  set_verbosity(0)
  
  data("scigenex_test_I1.2")
  res <- scigenex_test_I1.2
  
  expect_equal(length(gene_cluster(res)), 359)
  expect_equal(paste0(table(gene_cluster(res)), collapse = " "), "123 88 81 67")
  expect_equal(as.vector(table(gene_cluster(res, cluster = 1))), 123)
  expect_equal(as.vector(table(gene_cluster(res, cluster = 2))), 88)
  expect_equal(as.vector(table(gene_cluster(res, cluster = 3))), 81)
  expect_equal(as.vector(table(gene_cluster(res, cluster = 4))), 67)
  expect_equal(as.vector(table(gene_cluster(res, cluster = c(1,4)))), c(123, 67))
  expect_equal(as.vector(table(gene_cluster(res, cluster = c(4,4)))), 67)
  expect_error(gene_cluster(res, cluster = c(-1,8)))
  expect_error(gene_cluster(res, cluster = c(0,9)))
  expect_error(gene_cluster(res, cluster = c(0:8)))
  expect_error(gene_cluster(res, cluster = c(1:10)))
  expect_equal(paste0(res@gene_clusters[[1]], collapse = " "), 
               paste0(names(gene_cluster(res, cluster = 1)), collapse = " "))
  expect_equal(paste0(unlist(res@gene_clusters[1:2]), collapse = " "), 
               paste0(names(gene_cluster(res, cluster = 1:2)), collapse = " "))
  expect_equal(paste0(unlist(res@gene_clusters[1:4]), collapse = " "), 
               paste0(names(gene_cluster(res, cluster = 1:4)), collapse = " "))
})
