test_that("Checking gene_cluster() is working...", { 
  
  # Set verbosity to 0
  set_verbosity(0)
  
  load_example_dataset("7871581/files/pbmc3k_medium_clusters")
  res <- pbmc3k_medium_clusters
  
  expect_equal(length(gene_cluster(res)), 291)
  expect_equal(paste0(table(gene_cluster(res)), collapse = " "), "51 49 45 24 20 18 14 14 14 12 7 7 6 5 5")
  expect_equal(as.vector(table(gene_cluster(res, cluster = 1))), 51)
  expect_equal(as.vector(table(gene_cluster(res, cluster = 2))), 49)
  expect_equal(as.vector(table(gene_cluster(res, cluster = 3))), 45)
  expect_equal(as.vector(table(gene_cluster(res, cluster = 4))), 24)
  expect_equal(as.vector(table(gene_cluster(res, cluster = c(1,4)))), c(51, 24))
  expect_equal(as.vector(table(gene_cluster(res, cluster = c(4,4)))), 24)
  expect_error(gene_cluster(res, cluster = c(-1,8)))
  expect_error(gene_cluster(res, cluster = c(0,9)))
  expect_error(gene_cluster(res, cluster = c(0:8)))
  expect_error(gene_cluster(res, cluster = c(1,40)))
  expect_equal(paste0(res@gene_clusters[[1]], collapse = " "), 
               paste0(names(gene_cluster(res, cluster = 1)), collapse = " "))
  expect_equal(paste0(unlist(res@gene_clusters[1:2]), collapse = " "), 
               paste0(names(gene_cluster(res, cluster = 1:2)), collapse = " "))
  expect_equal(paste0(unlist(res@gene_clusters[1:4]), collapse = " "), 
               paste0(names(gene_cluster(res, cluster = 1:4)), collapse = " "))
})
