# Set verbosity to 0
set_verbosity(0)

#Create matrix containing 3 signatures
m <- create_4_rnd_clust()

## Select informative genes
res <- select_genes(m,
                    distance = "pearson",
                    k = 80,
                    highest = 0.3,
                    fdr = 1e-8,
                    row_sum = -Inf)

gn <- gene_clustering(res, keep_nn = TRUE, inflation = 1.1)

test_that("Check 'which_clust' is working.", {
  expect_true(all(which_clust(gn, c('gene27', 'gene336', 'gene187')) == 4:2))
  expect_true(all(is.na(which_clust(gn, c('gene27', 'gene336', 'bla'))) == c(F,F,T)))
})