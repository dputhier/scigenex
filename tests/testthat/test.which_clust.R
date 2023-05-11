set_verbosity(0)
library(Seurat)
load_example_dataset("7871581/files/pbmc3k_medium_clusters")
res <- pbmc3k_medium_clusters


test_that("Check 'which_clust' is working.", {
  expect_true(all(which_clust(res, c('APOBEC3A', 'RPL11', 'PF4')) == c(3,1,2)))
  expect_true(all(is.na(which_clust(res, c('APOBEC3A', 'RPL11', 'bla'))) == c(F,F,T)))
})
