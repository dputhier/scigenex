
testthat::test_that("Check if plot_profile is working properly", {
  # Set verbosity to 0
  set_verbosity(0)
  
  # Load seurat object
  data(pbmc_small, package = "SeuratObject")
  
  ## Select informative genes
  clust_set <- select_genes(data=pbmc_small,
                            distance_method="pearson",
                            k=10,
                            row_sum=-Inf,
                            noise_level=0.95,
                            fdr = 1e-6)
  
  ## Cluster genes
  clust_set <- gene_clustering(object = clust_set,
                               inflation = 1.2,
                               keep_nn = FALSE,
                               s = 5,
                               threads = 1)
  
  clust_set <- compute_centers(clust_set)
  p <- plot_profiles(clust_set, ident=Seurat::Idents(pbmc_small), size_text_y=5)
  testthat::expect_error(print(p), NA)
  
  p <- plot_profiles(clust_set[,colnames(clust_set@data)[1:80]], ident=Seurat::Idents(pbmc_small))
  testthat::expect_error(print(p), NA)
  
  p <- plot_profiles(clust_set[,colnames(clust_set@data)[1:10]], ident=Seurat::Idents(pbmc_small))
  testthat::expect_error(print(p), NA)
  
  p <- plot_profiles(clust_set[1,colnames(clust_set@data)[1:10]], ident=Seurat::Idents(pbmc_small))
  testthat::expect_error(print(p), NA)
  
  p <- plot_profiles(clust_set[2 ,colnames(clust_set@data)[1:10]], ident=Seurat::Idents(pbmc_small))
  testthat::expect_error(print(p), NA)
  
  p <- plot_profiles(clust_set[2 ,colnames(clust_set@data)[1:10]], ident=Seurat::Idents(pbmc_small))
  testthat::expect_error(print(p), NA)
  testthat::expect_error(plot_profiles(clust_set, ident=Seurat::Idents(pbmc_small), color_cell_type = rainbow(2), size_text_y=5))
  testthat::expect_error(plot_profiles(clust_set, ident=Seurat::Idents(pbmc_small), color_cell_type = rainbow(4), size_text_y=5))
  testthat::expect_error(plot_profiles(clust_set, ident=1:5, color_cell_type = rainbow(5), size_text_y=5))
  
})
