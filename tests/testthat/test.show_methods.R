test_that("Check clust_names() is working.", {
  a <- show_methods()
  expect_true(all(c("[", "clust_names", "clust_size", "cluster_stats", "col_names", 
                "dim", "enrich_go", "gene_cluster", "nclust", "rename_clust", 
                "row_names", "show", "top_genes", "viz_enrich", "which_clust")  %in% a))
})