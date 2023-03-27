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

test_that("Check 'rename' is working.", {
  expect_true(nclust(rename(gn, 1:nclust(gn))) == nclust(gn))
  expect_true(all(names(rename(gn, nclust(gn):1)@gene_clusters) == 4:1))
  rn <- rename(gn, letters[1:nclust(gn)])
  expect_is(plot_profiles(rn, ident = c(rep(1, 10), rep(2, 10))), "gg")
  plot_heatmap(rn, cell_clusters = setNames(c(rep(1, 10), rep(2, 10)), col_names(rn)))
})

