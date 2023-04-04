# Set verbosity to 0
set_verbosity(0)

#Create matrix containing 3 signatures
m <- create_4_rnd_clust()

## Select informative genes

data("scigenex_test_I1.2")
gn <- scigenex_test_I1.2

test_that("Check 'rename' is working.", {
  expect_true(nclust(rename_clust(gn, 1:nclust(gn))) == nclust(gn))
  expect_true(all(names(rename_clust(gn, nclust(gn):1)@gene_clusters) == 4:1))
  rn <- rename_clust(gn, letters[1:nclust(gn)])
  expect_is(plot_profiles(rn, ident = c(rep(1, 10), rep(2, 10))), "gg")
  plot_heatmap(rn, cell_clusters = setNames(c(rep(1, 10), rep(2, 10)), col_names(rn)))
})

