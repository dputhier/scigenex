library(testthat)
library(Seurat)
library(SeuratObject)
library(scigenex)

set_verbosity(0)

 

test_that("Checking cluster_set_from_matrix() #1", {
  
  m <- create_3_rnd_clust()[1:300,] 
  rownames(m) <- paste0("gene", 1:300)
  m_sub <- m[1:200, ]
  
  marker_sub <- list(a=paste0("gene", 1:100), 
                     b=paste0("gene", 101:200))
  
  markers <- list(a=paste0("gene", 1:100), 
                  b=paste0("gene", 101:200),
                  c=paste0("gene", 201:300))
  
  marker_fake <- list(a=letters[1:10])
  
  x <- cluster_set_from_matrix(m, marker_sub)
  testthat::expect_true(all(dim(x) == c(200, 20))) 
  
  x <- cluster_set_from_matrix(m_sub, markers)
  testthat::expect_true(all(dim(x) == c(200, 20))) 
  
  x <- cluster_set_from_matrix(m_sub, marker_fake)
  testthat::expect_true(all(dim(x) == c(0, 0))) 

})

test_that("Check cluster_set_from_seurat is working.", {
    data("pbmc_small", package="SeuratObject")
    markers <- Seurat::FindAllMarkers(pbmc_small, 
                                       only.pos = TRUE)
    markers <- markers[markers$p_val_adj <= 0.01, ]
    cs <- cluster_set_from_seurat(pbmc_small, markers, p_val_adj = 0.01)
    expect_true(nrow(markers) == nrow(cs))
    
    expect_true(all(gsub("\\.[0-9]+$", "", cs@gene_clusters[[1]]) == gsub("~[0-9]+$", "", markers$gene[markers$cluster == "0"], )))
    expect_true(all(gsub("\\.[0-9]+$", "", cs@gene_clusters[[2]]) == gsub("~[0-9]+$", "", markers$gene[markers$cluster == "1"], )))
    expect_true(all(gsub("\\.[0-9]+$", "", cs@gene_clusters[[3]]) == gsub("~[0-9]+$", "", markers$gene[markers$cluster == "2"], )))
    expect_true(nrow(cs[1,]@data) == nrow(pbmc_small[["RNA"]]$data[markers$gene[markers$cluster == "0"],]))
    expect_true(nrow(cs["0",]@data) == nrow(pbmc_small[["RNA"]]$data[markers$gene[markers$cluster == "0"],]))
    expect_true(sum(cs[1,]@data) == sum(pbmc_small[["RNA"]]$data[markers$gene[markers$cluster == "0"],]))
    expect_true(nrow(cs[2,]@data) == nrow(pbmc_small[["RNA"]]$data[markers$gene[markers$cluster == "1"],]))
    expect_true(nrow(cs["1",]@data) == nrow(pbmc_small[["RNA"]]$data[markers$gene[markers$cluster == "1"],]))
    expect_true(sum(cs[2,]@data) == sum(pbmc_small[["RNA"]]$data[markers$gene[markers$cluster == "1"],]))
    expect_true(nrow(cs[3,]@data) == nrow(pbmc_small[["RNA"]]$data[markers$gene[markers$cluster == "2"],]))
    expect_true(sum(cs[3,]@data) == sum(pbmc_small[["RNA"]]$data[markers$gene[markers$cluster == "2"],]))
    expect_true(cs@gene_clusters_metadata$number == 3)
    expect_true(all(cs@gene_clusters_metadata$cluster_id == setNames(c("0","1","2"), c("0","1","2"))))
    expect_true(all(cs@gene_clusters_metadata$cluster_id == setNames(c("0","1","2"), c("0","1","2"))))
    p <- plot_heatmap(cs[2,])
    testthat::expect_error(print(p), NA)
    p <- plot_heatmap(cs, cell_clusters = Idents(pbmc_small))
    testthat::expect_error(print(p), NA)
    cs <- compute_centers(cs)
    p <- plot_profiles(cs, ident = Idents(pbmc_small))
    testthat::expect_error(print(p), NA) 
})

test_that("Check cluster_set_from_seurat is working with vector.", {
  
  data("pbmc_small", package="SeuratObject")
  markers <- Seurat::FindAllMarkers(pbmc_small, 
                                    only.pos = TRUE)
  markers <- markers[markers$p_val_adj <= 0.01, ]
  markers <- setNames(as.character(markers$cluster), markers$gene)
  cs <- cluster_set_from_seurat(pbmc_small, markers)
  expect_true(length(markers) == nrow(cs))
  expect_true(all(gsub("\\.[0-9]+$", "", cs@gene_clusters[[1]]) == names(gsub("~[0-9]+$", "", names(markers)[markers == "0"]))))
  expect_true(all(gsub("\\.[0-9]+$", "", cs@gene_clusters[[2]]) == names(gsub("~[0-9]+$", "", names(markers)[markers == "1"]))))
  expect_true(all(gsub("\\.[0-9]+$", "", cs@gene_clusters[[3]]) == names(gsub("~[0-9]+$", "", names(markers)[markers == "2"]))))
  
  expect_true(nrow(cs[1,]@data) == nrow(pbmc_small[["RNA"]]$data[names(markers)[markers == "0"],]))
  expect_true(nrow(cs["0",]@data) == nrow(pbmc_small[["RNA"]]$data[names(markers)[markers == "0"],]))
  expect_true(nrow(cs[2,]@data) == nrow(pbmc_small[["RNA"]]$data[names(markers)[markers == "1"],]))
  expect_true(nrow(cs["1",]@data) == nrow(pbmc_small[["RNA"]]$data[names(markers)[markers == "1"],]))  
  expect_true(nrow(cs[2,]@data) == nrow(pbmc_small[["RNA"]]$data[names(markers)[markers == "1"],]))  
  expect_true(cs@gene_clusters_metadata$number == 3)
  expect_true(all(cs@gene_clusters_metadata$cluster_id == setNames(c("0","1","2"), c("0","1","2"))))
  expect_true(all(cs@gene_clusters_metadata$cluster_id == setNames(c("0","1","2"), c("0","1","2"))))
  p <- plot_heatmap(cs[2,])
  testthat::expect_error(print(p), NA)
  p <- plot_heatmap(cs, cell_clusters = Idents(pbmc_small))
  testthat::expect_error(print(p), NA)
  cs <- compute_centers(cs)
  p <- plot_profiles(cs, ident = Idents(pbmc_small))
  testthat::expect_error(print(p), NA) 
})
