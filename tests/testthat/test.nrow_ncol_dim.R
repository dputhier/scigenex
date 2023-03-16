# Set verbosity to 0
set_verbosity(0)

#Create matrix containing 3 signatures
m <- create_4_rnd_clust()

## Select informative genes
res <- select_genes(data=m,
                    distance_method="kendall",
                    k=75,
                    row_sum=-Inf,
                    highest=0.3,
                    fdr = 1e-8)

## Cluster genes
res <- gene_clustering(object = res,
                       inflation = 1.2,
                       keep_nn = FALSE,
                       k = 5,
                       threads = 1)

test_that("Checking ncol()", {
  expect_equal(ncol(res), 20)
})

test_that("Checking nrow()", {
  expect_equal(nrow(res), 359)
})

test_that("Checking dim()", {
  expect_equal(dim(res), c(359, 20))
})
