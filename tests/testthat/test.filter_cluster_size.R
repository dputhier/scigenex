set_verbosity(0)
obj <- new("ClusterSet")

n_smp <- 20
n_row <- 10

m1 <- matrix(1, ncol=n_smp, nrow=n_row)
m2 <- matrix(1, ncol=n_smp, nrow=n_row + 10)
m3 <- matrix(1, ncol=n_smp, nrow=n_row + 20)
m4 <- matrix(1, ncol=n_smp, nrow=n_row + 30)

rownames(m1) <- paste0("m1_", 1:n_row)
rownames(m2) <- paste0("m2_", 1:(n_row + 10))
rownames(m3) <- paste0("m3_", 1:(n_row + 20))
rownames(m4) <- paste0("m4_", 1:(n_row + 30))

obj@data <- as.matrix(rbind(m1, m2, m3, m4))
obj@gene_clusters <- list("1" = rownames(m1), "2"=rownames(m2), "3"=rownames(m3), "4"=rownames(m4))
obj@gene_clusters_metadata$cluster_id <- setNames(c("1", "2", "3", "4"), c("1", "2", "3", "4"))
obj@gene_clusters_metadata$size <- setNames(c(10, 20, 30, 40), c("1", "2", "3", "4"))
obj@dbf_output$center <- rbind(colMeans(m1), 
                               colMeans(m2),
                               colMeans(m3),
                               colMeans(m4))

test_that("Checking filter_cluster_size() #1", {
  x <- filter_cluster_size(obj, min_cluster_size = 9)
  testthat::expect_true(nclust(x) == 4) 
  x <- filter_cluster_size(obj, min_cluster_size = 10)
  testthat::expect_true(nclust(x) == 4) 
  x <- filter_cluster_size(obj, min_cluster_size = 20)
  testthat::expect_true(nclust(x) == 3) 
  x <- filter_cluster_size(obj, min_cluster_size = 21)
  testthat::expect_true(nclust(x) == 2) 
  testthat::expect_warning(filter_cluster_size(obj, min_cluster_size = 41)) 
})

