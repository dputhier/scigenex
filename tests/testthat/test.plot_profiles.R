test_that("Check if enrich_go stops when species argument is invalid", {
  data(pbmc_small, package = "SeuratObject")
  # Compute the signatures using find_gene_clusters()
  clust_set <- find_gene_clusters(pbmc_small, k=50, no_dknn_filter=TRUE)
  p <- plot_profiles(clust_set, ident=Seurat::Idents(pbmc_small), size_text_y=5)
  expect_error(plot_profiles(clust_set, ident=Seurat::Idents(pbmc_small), color_cell_type = rainbow(2), size_text_y=5))
  expect_error(plot_profiles(clust_set, ident=Seurat::Idents(pbmc_small), color_cell_type = rainbow(4), size_text_y=5))
  expect_error(plot_profiles(clust_set, ident=1:5, color_cell_type = rainbow(5), size_text_y=5))
})
