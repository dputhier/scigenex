# Set verbosity to 0
set_verbosity(0)

#Create matrix containing 3 signatures
m <- create_4_rnd_clust()
# Create a signature with 3 samples expressing it
m <- rbind(m, matrix(rep(c(rep(5, 30*2),
                           rep(0, 30*17),
                           rep(5, 30*1)), 20), nrow = 30, ncol = 20))
m <- rbind(m, matrix(rep(c(rep(0, 30*19),
                           rep(5, 30*1)), 20), nrow = 30, ncol = 20))

## Select informative genes
res <- select_genes(data=m,
                    distance_method="kendall",
                    k=75,
                    row_sum=-Inf,
                    dist_threads = 6,
                    highest=0.3,
                    fdr = 1e-8)

# Cluster informative genes
res <- gene_clustering(object = res,
                       inflation = 1.2,
                       keep_nn = FALSE,
                       threads = 6,
                       k = 5)



test_that("Checking filter_by_dot_prod() when av_dot_prod_min argument is set to 2", {
            
            nb_clusters_before_filtering <- length(res@gene_clusters)
            
            res_filtered <- filter_by_dot_prod(object = res,
                                               av_dot_prod_min = 2)
            
            # Check non-filtered clusters
            expect_equal(length(res_filtered@gene_clusters), 3)
            expect_equal(res_filtered@gene_clusters_metadata$number, 3)
            expect_equal(res_filtered@gene_clusters_metadata$cluster_id, c(1, 2, 3))
            
            # Check number of filtered out clusters
            expect_equal(nb_clusters_before_filtering - length(res_filtered@gene_clusters), 3)
            
            # Check data slot
            expect_equal(nrow(res_filtered@data), 203)
            expect_equal(rownames(res_filtered@data), unlist(res_filtered@gene_clusters, use.names = FALSE))
            
            # Check centers slot
            expect_equal(round(sum(res_filtered@dbf_output$center), 4), 67.6282)
            expect_equal(round(mean(res_filtered@dbf_output$center), 4), 1.1271)
            expect_equal(round(median(res_filtered@dbf_output$center), 4), 0)
            expect_equal(round(sd(res_filtered@dbf_output$center), 4), 1.6292)
          })






test_that("Checking filter_by_dot_prod() when av_dot_prod_min argument is set to 4", {
  
  nb_clusters_before_filtering <- length(res@gene_clusters)
  
  res_filtered <- filter_by_dot_prod(object = res,
                                     av_dot_prod_min = 4)
  
  # Check non-filtered clusters
  expect_equal(length(res_filtered@gene_clusters), 2)
  expect_equal(res_filtered@gene_clusters_metadata$number, 2)
  expect_equal(res_filtered@gene_clusters_metadata$cluster_id, c(1, 2))
  
  # Check number of filtered out clusters
  expect_equal(nb_clusters_before_filtering - length(res_filtered@gene_clusters), 4)
  
  # Check data slot
  expect_equal(nrow(res_filtered@data), 173)
  expect_equal(rownames(res_filtered@data), unlist(res_filtered@gene_clusters, use.names = FALSE))
  
  # Check centers slot
  expect_equal(round(sum(res_filtered@dbf_output$center), 4), 52.6282)
  expect_equal(round(mean(res_filtered@dbf_output$center), 4), 1.3157)
  expect_equal(round(median(res_filtered@dbf_output$center), 4), 1.0745)
  expect_equal(round(sd(res_filtered@dbf_output$center), 4), 1.5072)
})


test_that("Checking filter_by_dot_prod() stops when there is no conserved clusters", {
  
  nb_clusters_before_filtering <- length(res@gene_clusters)
  
  expect_error(filter_by_dot_prod(object = res, av_dot_prod_min = 40000))
})
