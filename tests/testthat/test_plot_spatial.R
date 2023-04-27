set_verbosity(0)
load_example_dataset("7870305/files/lymph_node_tiny_2")


test_that("Check plot_spatial is working.", {
  
  expect_true(ggplot2::is.ggplot(plot_spatial(seurat_obj = lymph_node_tiny_2, 
                                              gene_name = "CCL21", 
                                              intensity_slot="data")))
  expect_true(ggplot2::is.ggplot(plot_spatial(seurat_obj = lymph_node_tiny_2, 
                                              metadata = "nCount_Spatial")))
  expect_true(ggplot2::is.ggplot(plot_spatial(seurat_obj = lymph_node_tiny_2, 
                                              metadata = "nCount_Spatial", 
                                              pt_size = 5)))
  expect_error(ggplot2::is.ggplot(plot_spatial(seurat_obj = lymph_node_tiny_2, 
                                               gene_name = "UNDEF")))
  expect_true(ggplot2::is.ggplot(plot_spatial(seurat_obj = lymph_node_tiny_2, 
                                              gene_name = "CCL21", pt_shape = 10))) 
  expect_true(ggplot2::is.ggplot(plot_spatial(seurat_obj = lymph_node_tiny_2, 
                                              gene_name = "CCL21", pt_shape = 10, 
                                              size_title=10, face_title = 'bold', 
                                              title = "toto")))
  expect_true(ggplot2::is.ggplot(plot_spatial(seurat_obj = lymph_node_tiny_2, 
                                              gene_name = "CCL21", pt_shape = 16, 
                                              pt_star = F,
                                              size_title=10, face_title = 'bold', 
                                              title = "toto")))
})

