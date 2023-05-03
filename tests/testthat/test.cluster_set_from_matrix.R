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
