# Set verbosity to 0
set_verbosity(0)
library(Seurat)

load_example_dataset("7871581/files/pbmc3k_medium_clusters")
load_example_dataset("7871581/files/pbmc3k_medium")
gn <- pbmc3k_medium_clusters

test_that("Check 'rename' is working.", {
  expect_true(nclust(rename_clust(gn, 1:nclust(gn))) == nclust(gn))
  expect_true(all(names(rename_clust(gn, nclust(gn):1)@gene_clusters) == 15:1))
  rn <- rename_clust(gn, letters[1:nclust(gn)])
  expect_is(plot_profiles(rn, ident = Seurat::Idents(pbmc3k_medium)), "gg")
})

