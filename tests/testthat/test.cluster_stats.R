library(testthat)

# Set verbosity to 1 to display info messages only.
set_verbosity(0)

# Load a dataset
load_example_dataset('7871581/files/pbmc3k_medium_clusters')

# Compute some statistics about the clusters
df <- cluster_stats(pbmc3k_medium_clusters)  

testthat::test_that("Checking cluster_stats() #1", {
  testthat::expect_true(ncol(df) == 5)
  testthat::expect_true(nrow(df) == 15)
  testthat::expect_true(round(sum(df),0) == 3528)
})

  