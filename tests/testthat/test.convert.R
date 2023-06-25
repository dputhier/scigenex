set_verbosity(0)

m <- create_3_rnd_clust()[1:300,] 
rownames(m) <- paste0("gene", 1:300)
m_sub <- m[1:200, ]

marker_sub <- list(a=paste0("gene", 1:100), 
                b=paste0("gene", 101:200))

markers <- list(a=paste0("gene", 1:100), 
                 b=paste0("gene", 101:200),
                 c=paste0("gene", 201:300))

marker_fake <- list(a=letters[1:10]) 

test_that("Checking cluster_set_from_matrix() #1", {
  x <- cluster_set_from_matrix(m, marker_sub)
  testthat::expect_true(all(dim(x) == c(200, 20))) 
  
  x <- cluster_set_from_matrix(m_sub, markers)
  testthat::expect_true(all(dim(x) == c(200, 20))) 
  
  x <- cluster_set_from_matrix(m_sub, marker_fake)
  testthat::expect_true(all(dim(x) == c(0, 0))) 

})
# Set verbosity to 0
set_verbosity(0)

library(Seurat)
library(SeuratObject)
data("pbmc_small", package="SeuratObject")
markers <- FindAllMarkers(pbmc_small)
cs <- cluster_set_from_seurat(pbmc_small, markers)

test_that("Check cluster_set_from_seurat is working.", {
 expect_true(all(cs@gene_clusters[[1]] == markers$gene[markers$cluster == "0"]))
 expect_true(all(cs@gene_clusters[[2]] == markers$gene[markers$cluster == "1"]))
 expect_true(all(cs@gene_clusters[[3]] == markers$gene[markers$cluster == "2"]))
 expect_true(nrow(cs[1,]@data) == nrow(pbmc_small@assays$RNA[markers$gene[markers$cluster == "0"],]))
 expect_true(nrow(cs["0",]@data) == nrow(pbmc_small@assays$RNA[markers$gene[markers$cluster == "0"],]))
 expect_true(sum(cs[1,]@data) == sum(pbmc_small@assays$RNA[markers$gene[markers$cluster == "0"],]))
 expect_true(nrow(cs[2,]@data) == nrow(pbmc_small@assays$RNA[markers$gene[markers$cluster == "1"],]))
 expect_true(nrow(cs["1",]@data) == nrow(pbmc_small@assays$RNA[markers$gene[markers$cluster == "1"],]))
 expect_true(sum(cs[2,]@data) == sum(pbmc_small@assays$RNA[markers$gene[markers$cluster == "1"],]))
 expect_true(nrow(cs[3,]@data) == nrow(pbmc_small@assays$RNA[markers$gene[markers$cluster == "2"],]))
 expect_true(sum(cs[3,]@data) == sum(pbmc_small@assays$RNA[markers$gene[markers$cluster == "2"],]))
 expect_true(cs@gene_clusters_metadata$number == 3)
 expect_true(all(cs@gene_clusters_metadata$cluster_id == setNames(c("0","1","2"), c("0","1","2"))))
 expect_true(all(cs@gene_clusters_metadata$cluster_id == setNames(c("0","1","2"), c("0","1","2"))))
 #' plot_heatmap(cs[1,])
 #' plot_heatmap(cs, cell_clusters = Idents(pbmc_small))
 #' plot_profiles(cs, ident = Idents(pbmc_small))
 p <- plot_heatmap(cs[1,])
 testthat::expect_error(print(p), NA)
 p <- plot_heatmap(cs, cell_clusters = Idents(pbmc_small))
 testthat::expect_error(print(p), NA)
 p <- plot_profiles(cs, ident = Idents(pbmc_small))
 testthat::expect_error(print(p), NA) 
})

