test_that("Checking gene_cluster() is working...", { 
  
  set_verbosity(0)
  data("complex9Noisy")
  res <- find_gene_clusters(data=complex9Noisy[ ,1:2],
                            distance_method="euclidean",
                            inflation = 1.25,
                            k=20,
                            min_nb_supporting_cell=0,
                            highest = 0.95,
                            row_sum = -Inf,
                            fdr = 0.0001)
  
  expect_equal(length(gene_cluster(res)), 3286)
  expect_equal(paste0(table(gene_cluster(res)), collapse = " "), "1023 782 419 373 352 126 122 89")
  expect_equal(as.vector(table(gene_cluster(res, cluster = 1))), 1023)
  expect_equal(as.vector(table(gene_cluster(res, cluster = 2))), 782)
  expect_equal(as.vector(table(gene_cluster(res, cluster = 3))), 419)
  expect_equal(as.vector(table(gene_cluster(res, cluster = 8))), 89)
  expect_equal(as.vector(table(gene_cluster(res, cluster = 8))), 89)
  expect_equal(as.vector(table(gene_cluster(res, cluster = c(1,8)))), c(1023, 89))
  expect_equal(as.vector(table(gene_cluster(res, cluster = c(8,8)))), 89)
  expect_error(gene_cluster(res, cluster = c(-1,8)))
  expect_error(gene_cluster(res, cluster = c(0,9)))
  expect_error(gene_cluster(res, cluster = c(0:8)))
  expect_error(gene_cluster(res, cluster = c(1:10)))
  expect_equal(paste0(res@gene_clusters[[1]], collapse = " "), 
               paste0(names(gene_cluster(res, cluster = 1)), collapse = " "))
  expect_equal(paste0(unlist(res@gene_clusters[1:2]), collapse = " "), 
               paste0(names(gene_cluster(res, cluster = 1:2)), collapse = " "))
  expect_equal(paste0(unlist(res@gene_clusters[1:8]), collapse = " "), 
               paste0(names(gene_cluster(res, cluster = 1:8)), collapse = " "))
})