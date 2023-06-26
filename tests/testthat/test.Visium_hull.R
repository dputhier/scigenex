library(testthat)
library(Seurat)
library(ggplot2)
set_verbosity(0)

test_that("Check getFlippedTissueCoordinates() is working.", {
  load_example_dataset("7870305/files/lymph_node_tiny_2")
  lymph_node_tiny <- getFlippedTissueCoordinates(lymph_node_tiny_2)
  df <- getFlippedTissueCoordinates(lymph_node_tiny_2, as_data_frame=TRUE)
  
  testthat::expect_true(ncol(df) == 2)
  testthat::expect_true(nrow(df) == 442)
})




test_that("Check display_hull() is working.", {
  load_example_dataset("7870305/files/lymph_node_tiny_2")
  load_example_dataset("7870305/files/lymph_node_tiny_clusters_2")
  lymph_node_tiny_2 <- Seurat::AddModuleScore(lymph_node_tiny_2, features = lymph_node_tiny_clusters_2@gene_clusters, nbin = 10)
  p <- Seurat::SpatialDimPlot(lymph_node_tiny_2, pt.size.factor = 4)
  hull <- display_hull(lymph_node_tiny_2, 
                       ident=ifelse(Seurat::Idents(lymph_node_tiny_2) %in% c(7, 8), 1, 0),
                       delta=1, size_x=3.4, size_y=3)
  hull <- display_hull(lymph_node_tiny_2, 
                       ident=ifelse(Seurat::  Idents(lymph_node_tiny_2) %in% c(7, 8), 1, 0),
                       delta=1, size_x=3.4, size_y=3, color="black")
  
  testthat::expect_true(inherits(hull, "LayerInstance"))
  testthat::expect_true(all(dim(hull$data) == c(162, 4)))
  testthat::expect_true(all(dim(hull$data) == c(162, 4)))
})
