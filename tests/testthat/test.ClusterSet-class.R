testthat::test_that("Test the indexing operator of the ClusterSet", {
  
  clusters <- 1:5
  cluster_size <- c(12,10,8,6,4)
  m <- matrix(rnorm(20000), nc=500, nr=40)
  colnames(m) <- paste0("S_", 1:500)
  gn <- paste0("G_", 1:40)
  cls <- unlist(mapply(rep, clusters, cluster_size))
  gene_cluster <- split(gn, cls) 
  rownames(m) <- paste0("G_", 1:40)
  cls <- unlist(mapply(rep, clusters, cluster_size))
  gene_clusters_metadata <- list()
  gene_clusters_metadata$cluster_id <- as.numeric(names(gene_cluster))
  gene_clusters_metadata$number <- length(clusters)
  gene_clusters_metadata$size <- cluster_size
  names(gene_clusters_metadata$size) <- clusters
  gene_cluster_annotations <- list()
  top_genes <- lapply(gene_cluster, "[", 1:3)
  dbf_output <- list()
  dbf_output$dknn <- runif(1000)
  dbf_output$simulated_dknn <- runif(1000)
  dbf_output$critical_distance <- 0.5
  dbf_output$fdr <- 0.05
  dbf_output$center <- do.call(rbind, 
                               lapply(split.data.frame(as.data.frame(m), 
                                                       as.factor(cls)), 
                                      apply, 2, mean))
  
  colnames(dbf_output$center) <- colnames(m)
  rownames(dbf_output$center) <- names(gene_cluster)
  
  res <- new("ClusterSet",
              data = m,
              gene_clusters = split(gn, cls),
              top_genes = top_genes,
              gene_clusters_metadata = gene_clusters_metadata,
              gene_cluster_annotations = gene_cluster_annotations,
              cells_metadata = data.frame(cells_barcode=colnames(m), 
                                          row.names = colnames(m)),
              dbf_output = dbf_output,
              parameters = list(filename="oP0H7t6Sc0", distance_method="pearson", 
                                k=50, inflation=2, highest=0.3, fdr=0.05, 
                                min_nb_supporting_cell=0, min_pct_gene_expressed=40, 
                                min_cluster_size=10, row_sum=-Inf, seed=123)
            )
  testthat::expect_equal(dim(res[1:3,4:5]), c(30, 2))
  testthat::expect_equal(dim(res[c(2,4),1:5]), c(16, 5))
  testthat::expect_equal(length(res[c(2,4),1:5]@gene_clusters), 2)
  testthat::expect_equal(length(res[c(2,4),1:5]@gene_clusters), 2)
  testthat::expect_equal(length(res[c(1,2,4),]@gene_clusters), 3)
  testthat::expect_equal(dim(res[logical(), logical()]), c(0,0))
  testthat::expect_equal(dim(res[c("1", "2"), ]), c(22, 500))
  testthat::expect_equal(res[c("1", "2"), ]@gene_clusters_metadata$number, 2)
  testthat::expect_equal(res[2:4, ]@gene_clusters_metadata$number, 3)
  testthat::expect_equal(res[c("1", "2", "3"), ]@gene_clusters_metadata$number, 3)
  testthat::expect_equal(res[logical(), ]@gene_clusters_metadata$number, 0)
  testthat::expect_equal(dim(res[, ]), c(40, 500))
  testthat::expect_equal(dim(res[, 10:20]), c(40, 11))
  testthat::expect_equal(dim(res[1:3, 10:20]@cells_metadata), c(11, 1))
  testthat::expect_equal(dim(res[1:3, 10:20]@cells_metadata), c(11, 1))
  testthat::expect_equal(dim(res[1, 1]@cells_metadata), c(1, 1))
  testthat::expect_equal(dim(res[1, 1]), c(12, 1))
  testthat::expect_equal(length(res[1:2, 1]@top_genes), 2)
  testthat::expect_equal(length(res[1:2, c()]@top_genes), 2)
  
  testthat::expect_identical(res[c(1, 3), ], res[c(1, 3), seq_len(ncol(res))])
  testthat::expect_identical(res[, c("S_408", "S_409")], res[, c(408, 409)])
  testthat::expect_identical(res[c(1, 3), c("S_3", "S_4")], res[c(1, 3), 3:4])
  testthat::expect_identical(res[-c(1, 2), ], res[c(3:nclust(res)), ])
  testthat::expect_identical(res[as.character(2:5), c("S_3", "S_4")], res[2:5, 3:4])
  testthat::expect_identical(res[as.character(2:5), ], res[2:5, ])
  testthat::expect_identical(res[,paste0("S_",2:5)], res[,2:5])
  testthat::expect_identical(res[c("2", "4"), ]@gene_clusters, res[c(2, 4), ]@gene_clusters)
  testthat::expect_identical(res[c("1", "3"), ]@gene_clusters_metadata$number, res[c(1, 3), ]@gene_clusters_metadata$number)
  testthat::expect_identical(res[3, ]@gene_clusters_metadata, res["3", ]@gene_clusters_metadata)
  testthat::expect_error(res["a", "b"])
  })

