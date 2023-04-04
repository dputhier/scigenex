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



test_that("Checking filter_nb_supporting_cells() when min_nb_supporting_cell\
          argument is set to 10 and min_pct_gene_expressed to 80", {
            
            nb_clusters_before_filtering <- length(res@gene_clusters)
            
            res_filtered <- filter_nb_supporting_cells(object = res,
                                                       min_nb_supporting_cell = 10,
                                                       min_pct_gene_expressed = 80)
            
            # Check non-filtered clusters
            expect_equal(length(res_filtered@gene_clusters), 4)
            expect_equal(res_filtered@gene_clusters_metadata$number, 4)
            expect_equal(res_filtered@gene_clusters_metadata$cluster_id, c(1, 2, 3, 4))
            
            # Check number of filtered out clusters
            expect_equal(nb_clusters_before_filtering - length(res_filtered@gene_clusters), 2)
            
            # Check data slot
            expect_equal(nrow(res_filtered@data), 306)
            expect_equal(rownames(res_filtered@data), unlist(res_filtered@gene_clusters, use.names = FALSE))
            
            # Check centers slot
            expect_equal(round(sum(res_filtered@dbf_output$center), 4), 19.1442)
            expect_equal(round(mean(res_filtered@dbf_output$center), 4), 0.2393)
            expect_equal(round(median(res_filtered@dbf_output$center), 4), 0.0201)
            expect_equal(round(sd(res_filtered@dbf_output$center), 4), 1.7346)
          })




test_that("Checking filter_nb_supporting_cells() when min_nb_supporting_cell\
          argument is set to 2 and min_pct_gene_expressed to 80", {
            
            nb_clusters_before_filtering <- length(res@gene_clusters)
            # Args set to min_cluster_size = 5
            res_filtered <- filter_nb_supporting_cells(object = res,
                                                       min_nb_supporting_cell = 2,
                                                       min_pct_gene_expressed = 80)
            
            # Check non-filtered clusters
            expect_equal(length(res_filtered@gene_clusters), 5)
            expect_equal(res_filtered@gene_clusters_metadata$number, 5)
            expect_equal(res_filtered@gene_clusters_metadata$cluster_id, c(1, 2, 3, 4, 5))
            
            # Check number of filtered out clusters
            expect_equal(nb_clusters_before_filtering - length(res_filtered@gene_clusters), 1)
            
            # Check data slot
            expect_equal(nrow(res_filtered@data), 336)
            expect_equal(rownames(res_filtered@data), unlist(res_filtered@gene_clusters, use.names = FALSE))
            
            # Check centers slot
            expect_equal(round(sum(res_filtered@dbf_output$center), 4), 34.1442)
            expect_equal(round(mean(res_filtered@dbf_output$center), 4), 0.3414)
            expect_equal(round(median(res_filtered@dbf_output$center), 4), 0)
            expect_equal(round(sd(res_filtered@dbf_output$center), 4), 1.757)
          })

test_that("Checking filter_nb_supporting_cells() stops when min_nb_supporting_cell\
          argument is set to 21 and min_pct_gene_expressed to 100", {
            
            nb_clusters_before_filtering <- length(res@gene_clusters)
            # Args set to min_cluster_size = 5
            expect_error(filter_nb_supporting_cells(object = res,
                                                    min_nb_supporting_cell = 21,
                                                    min_pct_gene_expressed = 100))
            
          })
