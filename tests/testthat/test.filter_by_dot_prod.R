
#' To do this, the function first converts the gene expression 
#' data for each cluster into a binary form (values greater than 
#' 1 are set to 1). Then it calculates the dot product for this 
#' binary matrix, which produces a gene-gene matrix showing the 
#' number of cells/spots where each pair of genes are expressed 
#' together. The function then calculates the median value of the 
#' maximum concordances across all genes, which can be used to 
#' determine whether a cluster should be filtered out or not.
#' 


obj <- new("ClusterSet")

n_smp <- 20
n_row <- 10

m1 <- m2 <- m3 <- m4 <- matrix(0, ncol=n_smp, nrow=n_row)

rownames(m1) <- paste0("m1_", 1:n_row)
rownames(m2) <- paste0("m2_", 1:n_row)
rownames(m3) <- paste0("m3_", 1:n_row)
rownames(m4) <- paste0("m4_", 1:n_row)

m1[, 1] <- 1

m2[, 1] <- 1
m2[, 2] <- 1

m3[, 8] <- 1
m3[, 7] <- 1
m3[, 6] <- 1

m4[, 3] <- 1
m4[, 4] <- 1
m4[, 9] <- 1
m4[, 10] <- 1

obj@data <- as.matrix(rbind(m1, m2, m3, m4))
obj@gene_clusters <- list("1" = rownames(m1), "2"=rownames(m2), "3"=rownames(m3), "4"=rownames(m4))
obj@gene_clusters_metadata$cluster_id <- setNames(c("1", "2", "3", "4"), c("1", "2", "3", "4"))
obj@gene_clusters_metadata$size <- setNames(c("1", "2", "3", "4"), rep(n_row, 4))
obj@dbf_output$center <- rbind(colMeans(m1), 
                               colMeans(m2),
                               colMeans(m3),
                               colMeans(m4))
rownames(obj@dbf_output$center) <- names(obj@gene_clusters)


res <- obj

test_that("Checking filter_by_dot_prod() when av_dot_prod_min argument is set to 2", {
    
  new_obj <- filter_by_dot_prod(res, av_dot_prod_min = 1)  
  expect_true(nclust(new_obj) == 3)

  new_obj <- filter_by_dot_prod(res, av_dot_prod_min = 2)  
  expect_true(nclust(new_obj) == 2) 
  
  new_obj <- filter_by_dot_prod(res, av_dot_prod_min = 3)  
  expect_true(nclust(new_obj) == 1)    
  
  new_obj <- filter_by_dot_prod(res, av_dot_prod_min = 4)  
  expect_true(nclust(new_obj) == 0)
})

# n_dbf_output$center -----------------------------------------------------

