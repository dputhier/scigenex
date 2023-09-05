library(testthat)

# Set verbosity to 0
set_verbosity(0)
load_example_dataset('7871581/files/pbmc3k_medium_clusters')
# DNA Binding: "GO:0003677"
pbmc3k_medium_clusters <- top_by_go(pbmc3k_medium_clusters, go_id = "GO:0003677")

testthat::test_that("Checking top_by_go() #1", {
  testthat::expect_true(length(unlist(pbmc3k_medium_clusters@top_genes)) > 15)
  testthat::expect_true("GATA2"  %in% unlist(pbmc3k_medium_clusters@top_genes))
})

testthat::test_that("Checking top_by_go() #2", {
  testthat::expect_true("MAFB"  %in% unlist(pbmc3k_medium_clusters@top_genes))
  testthat::expect_true("GFI1B"  %in% unlist(pbmc3k_medium_clusters@top_genes))  
})