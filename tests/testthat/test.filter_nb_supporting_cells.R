set_verbosity(0)

obj <- new("ClusterSet")

n_smp <- 20
n_row <- 10

m1 <- m2 <- m3 <- m4 <- matrix(0, ncol=n_smp, nrow=n_row)

rownames(m1) <- paste0("m1", 1:n_row)
rownames(m2) <- paste0("m2", 1:n_row)
rownames(m3) <- paste0("m3", 1:n_row)
rownames(m4) <- paste0("m4", 1:n_row)

m1[1:4, 1] <- 1

m2[1:4, 1] <- 1
m2[1:4, 2] <- 1

m3[1:3, 8] <- 1
m3[1:5, 7] <- 1
m3[1:6, 6] <- 1

m4[1:6, 3] <- 1
m4[1:6, 4] <- 1
m4[1:6, 9] <- 1
m4[1:6, 10] <- 1

obj@data <- as.matrix(rbind(m1, m2, m3, m4))
obj@gene_clusters <- list("1" = rownames(m1), "2"=rownames(m2), "3"=rownames(m3), "4"=rownames(m4))
obj@gene_clusters_metadata$cluster_id <- setNames(c("1", "2", "3", "4"), c("1", "2", "3", "4"))
obj@gene_clusters_metadata$size <- setNames(c("1", "2", "3", "4"), rep(n_row, 4))
obj@dbf_output$center <- rbind(colMeans(m1), 
                                  colMeans(m2),
                                  colMeans(m3),
                                  colMeans(m4))

test_that("Checking filter_nb_supporting_cells(), min_nb_supporting_cell=1, min_pct_gene_expressed = 40", {
            
  x <- filter_nb_supporting_cells(object = obj, 
                             min_nb_supporting_cell = 1, 
                             min_pct_gene_expressed = 40)
  testthat::expect_true(length(x@gene_clusters) == 4)
  testthat::expect_true(length(x@gene_clusters_metadata$cluster_id) == 4)
  testthat::expect_true(nrow(x@dbf_output$center) == 4)
  
  })

test_that("Checking filter_nb_supporting_cells(), min_nb_supporting_cell=2, min_pct_gene_expressed = 40", {
  
  x <- filter_nb_supporting_cells(object = obj, 
                                  min_nb_supporting_cell = 2, 
                                  min_pct_gene_expressed = 40)
  testthat::expect_true(length(x@gene_clusters) == 3)
  testthat::expect_true(length(x@gene_clusters_metadata$cluster_id) == 3)
  testthat::expect_true(nrow(x@dbf_output$center) == 3)
  
})

test_that("Checking filter_nb_supporting_cells(), min_nb_supporting_cell=3, min_pct_gene_expressed = 40", {
  
  x <- filter_nb_supporting_cells(object = obj, 
                                  min_nb_supporting_cell = 3, 
                                  min_pct_gene_expressed = 40)
  testthat::expect_true(length(x@gene_clusters) == 1)
  testthat::expect_true(length(x@gene_clusters_metadata$cluster_id) == 1)
  testthat::expect_true(nrow(x@dbf_output$center) == 1)
  
})

test_that("Checking filter_nb_supporting_cells(), min_nb_supporting_cell=1, min_pct_gene_expressed = 90", {
  
  x <- filter_nb_supporting_cells(object = obj, 
                                  min_nb_supporting_cell = 3, 
                                  min_pct_gene_expressed = 40)
  testthat::expect_true(length(x@gene_clusters) == 1)
  testthat::expect_true(length(x@gene_clusters_metadata$cluster_id) == 1)
  testthat::expect_true(nrow(x@dbf_output$center) == 1)
  
})

test_that("Checking filter_nb_supporting_cells(), min_nb_supporting_cell=2, min_pct_gene_expressed = 50", {
  
  x <- filter_nb_supporting_cells(object = obj, 
                                  min_nb_supporting_cell = 2, 
                                  min_pct_gene_expressed = 50)
  testthat::expect_true(length(x@gene_clusters) == 2)
  testthat::expect_true(length(x@gene_clusters_metadata$cluster_id) == 2)
  testthat::expect_true(nrow(x@dbf_output$center) == 2)
  
})

test_that("Checking filter_nb_supporting_cells(), min_nb_supporting_cell=10, min_pct_gene_expressed = 100", {
  
  x <- filter_nb_supporting_cells(object = obj, 
                                  min_nb_supporting_cell = 10, 
                                  min_pct_gene_expressed = 100)
  testthat::expect_true(length(x@gene_clusters) == 0)
  testthat::expect_true(length(x@gene_clusters_metadata$cluster_id) == 0)
  testthat::expect_true(nrow(x@dbf_output$center) == 0)
  
})
