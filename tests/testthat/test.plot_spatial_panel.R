set_verbosity(0)
load_example_dataset("7870305/files/lymph_node_tiny_2")
load_example_dataset("7870305/files/lymph_node_tiny_clusters_2")
lymph_node_tiny_2 <- AddModuleScore(lymph_node_tiny_2, features = lymph_node_tiny_clusters_2@gene_clusters, nbin = 15)

for(i in 1:nclust(lymph_node_tiny_clusters_2)){ # Normalizing module scores
 tmp <- lymph_node_tiny_2[[paste0("Cluster", i, sep="")]] 
 max_tmp <- max(tmp)
 min_tmp <- min(tmp)
 lymph_node_tiny_2[[paste0("Cluster", i, sep="")]]  <- (tmp[,1] - min(tmp))/(max_tmp - min_tmp)
}

test_that("Check plot_spatial is working.", {
  p <- plot_spatial_panel(lymph_node_tiny_2, metadata=paste0("Cluster", 1:4), ncol_layout=2,
                          guides='collect', pt_size=2.2, coord_flip=T)
  expect_true(ggplot2::is.ggplot(p))
  
  p <- plot_spatial_panel(lymph_node_tiny_2, gene=c('VPREB3', 'IGHG1', 'PRDX4', 
                                                    'LTB', 'CCL20', 'LYVE1', 
                                                    'IL7R', 'RGS9', 'MAPT'), 
                          ncol_layout=3,
                          pt_size=1.5, coord_flip=T, panel_names=LETTERS[1:9])
  
  expect_true(ggplot2::is.ggplot(p))
  
  p <- plot_spatial_panel(lymph_node_tiny_2, gene=c('VPREB3', 'IGHG1', 'PRDX4', 
                                                    'LTB', 'CCL20', 'LYVE1', 
                                                    'IL7R', 'RGS9', 'MAPT'), 
                          ncol_layout=3,
                          pt_size=1.5, coord_flip=T, panel_names=LETTERS[1:9], size_title = 10)
  
  expect_true(ggplot2::is.ggplot(p)) 
})

