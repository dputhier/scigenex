load_example_dataset("7871581/files/pbmc3k_medium_clusters")
set_verbosity(0)

test_that("Check clust_names() is working.", {
  a <- unique(unname(gene_cluster(pbmc3k_medium_clusters)))
  b <- clust_names(pbmc3k_medium_clusters)
  expect_true(all(as.character(a)==b))
})
