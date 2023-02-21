# Create matrix containing 3 signatures
m <- create_3_rnd_clust()
res_scigenex <- find_gene_clusters(
  data = m,
  name = "test",
  distance_method = "pearson",
  inflation = 2,
  k = 25,
  fdr = 10
)

test_that("Just some check for create_3_rnd_clust()", {
  expect_equal(ncol(res_scigenex), 20)
})

test_that("Just some check for create_3_rnd_clust()", {
  expect_equal(nrow(res_scigenex), 134)
})

test_that("Just some check for create_3_rnd_clust()", {
  expect_equal(dim(res_scigenex), c(134, 20))
})
