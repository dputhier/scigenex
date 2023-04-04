# Set verbosity to 0
set_verbosity(0)

#Create matrix containing 3 signatures
m <- create_4_rnd_clust()

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
                       inflation = 4,
                       keep_nn = FALSE,
                       k = 5)

test_that("Checking filter_cluster_size() when min_cluster_size\
          argument is set to 5", {
            
            nb_clusters_before_filtering <- length(res@gene_clusters)

            # Args set to min_cluster_size = 5
            res_size5 <- filter_cluster_size(object = res,
                                             min_cluster_size = 5)
            
            # Check non-filtered clusters
            expect_equal(length(res_size5@gene_clusters), 1)
            expect_equal(res_size5@gene_clusters_metadata$number, 1)
            expect_equal(res_size5@gene_clusters_metadata$cluster_id, 1)
            
            # Check number of filtered out clusters
            expect_equal(nb_clusters_before_filtering - length(res_size5@gene_clusters), 243)
            
            # Check data slot
            expect_equal(nrow(res_size5@data), 7)
            expect_equal(rownames(res_size5@data), unlist(res_size5@gene_clusters, use.names = FALSE))
            
            # Check centers slot
            expect_equal(round(sum(res_size5@dbf_output$center), 4), -15.0361)
            expect_equal(round(mean(res_size5@dbf_output$center), 4), -0.7518)
            expect_equal(round(median(res_size5@dbf_output$center), 4), -0.3019)
            expect_equal(round(sd(res_size5@dbf_output$center), 4), 1.5521)
          })


test_that("Checking filter_cluster_size() when min_cluster_size\
          argument is set to 1", {
            
            nb_clusters_before_filtering <- length(res@gene_clusters)
            # Args set to min_cluster_size = 1
            res_size1 <- filter_cluster_size(object = res,
                                             min_cluster_size = 1)
            
            # Check non-filtered clusters
            expect_equal(length(res_size1@gene_clusters), 74)
            expect_equal(res_size1@gene_clusters_metadata$number, 74)
            expect_equal(res_size1@gene_clusters_metadata$cluster_id, seq(1, 74))
            
            # Check number of filtered out clusters
            expect_equal(nb_clusters_before_filtering - length(res_size1@gene_clusters), 170)
            
            # Check data slot
            expect_equal(nrow(res_size1@data), 189)
            expect_equal(rownames(res_size1@data), unlist(res_size1@gene_clusters, use.names = FALSE))
            
            # Check centers slot
            expect_equal(round(sum(res_size1@dbf_output$center), 4), 522.9372)
            expect_equal(round(mean(res_size1@dbf_output$center), 4), 0.3533)
            expect_equal(round(median(res_size1@dbf_output$center), 4), 0.1936)
            expect_equal(round(sd(res_size1@dbf_output$center), 4), 1.8418)
          })




test_that("Checking if all clusters are conserved when filter_cluster_size() \
          with min_cluster_size set to 0", {
            
            nb_clusters_before_filtering <- length(res@gene_clusters)
            # Args set to min_cluster_size = 0
            res_size0 <- filter_cluster_size(object = res,
                                             min_cluster_size = 0)
            
            # Check non-filtered clusters
            expect_equal(length(res_size0@gene_clusters), 244)
            expect_equal(res_size0@gene_clusters_metadata$number, 244)
            expect_equal(res_size0@gene_clusters_metadata$cluster_id, seq(1, 244, 1))
            
            # Check number of filtered out clusters
            expect_equal(nb_clusters_before_filtering - length(res_size0@gene_clusters), 0)
            
            # Check data slot
            expect_equal(nrow(res_size0@data), 359)
            expect_equal(rownames(res_size0@data), unlist(res_size0@gene_clusters, use.names = FALSE))
            
            # Check centers slot
            expect_equal(round(sum(res_size0@dbf_output$center), 3), 1846.763)
            expect_equal(round(mean(res_size0@dbf_output$center), 4), 0.3784)
            expect_equal(round(median(res_size0@dbf_output$center), 4), 0.1386)
            expect_equal(round(sd(res_size0@dbf_output$center), 4), 2.1379)
            
            # Check if clusters are the same in the inital ClusterSet object
            expect_equal(rownames(res_size0@data), rownames(res@data))
            expect_equal(unlist(res_size0@gene_clusters, use.names = FALSE),
                         unlist(res@gene_clusters, use.names = FALSE))
          })





test_that("Checking if filter_cluster_size() stops when all clusters are filtered out", {
            
            nb_clusters_before_filtering <- length(res@gene_clusters)
            # Args set to min_cluster_size = 50
            expect_error(filter_cluster_size(object = res,
                                             min_cluster_size = 50))
          })

