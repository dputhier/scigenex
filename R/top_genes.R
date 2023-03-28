#################################################################
##    Define top_genes function for ClusterSet object
#################################################################

#' @title Most co-expressed genes from each gene cluster
#' @description
#' Extract the most highly co-expressed genes of each gene cluster.
#' @param object A \code{ClusterSet} object.
#' @param top A value for the number of genes to select from each cluster.
#' @param cluster A vector of gene cluster identity.
#'
#' @return A \code{ClusterSet} object.
#' @export top_genes
#'
#' @examples
#' # Set verbosity to 1 to only display info messages.
#' set_verbosity(1)
#' 
#' # Create a matrix with 4 signatures
#' m <- create_4_rnd_clust()
#' 
#' # Select informative genes
#' res <- select_genes(m,
#'                     distance = "kendall",
#'                     k = 75,
#'                     highest = 0.3,
#'                     fdr = 1e-8,
#'                     row_sum = -Inf)
#'                     
#' # Cluster informative features
#' res <- gene_clustering(res, 
#'                        inflation = 1.2,
#'                        keep_nn = FALSE,
#'                        k = 5)
#'                        
#' # Identify top 5 genes of each gene cluster
#' res <- top_genes(res, top = 5)
#' 
#' # Plot heatmap of gene clusters
#' plot_heatmap(res, use_top_genes = TRUE)
#'
top_genes <- function(object,
                      top = 20,
                      cluster = "all") {
  
  ## Check format object arg
  check_format_cluster_set(object)
  
  if (unique(cluster == "all")) {
    cluster <- object@gene_clusters_metadata$cluster_id
  }
  
  # Display a warning message if there is less than n top genes in a gene cluster
  loop <- 0
  for (size in object@gene_clusters_metadata$size[cluster]) {
    loop <- loop + 1

    if (top > size) {
      print_msg(paste0("Number of top genes is greater than the number of genes in cluster ", 
                       loop, 
                       ". All genes will be used and ordered by similarity rank."),
                msg_type = "WARNING")
    }
  }
  
  # Initialization for the for loop
  clusters <- object@gene_clusters
  genes_top <- matrix(ncol = top)
  
  # Extract top co-expressed genes for each gene cluster
  for (i in cluster) {
    # Extract gene names in cluster i
    genes <- clusters[[i]]
    
    # Compute distances between genes in cluster i
    # Use the same distance used by SciGeneX
    dist_method <- object@parameters$distance_method
    
    if (dist_method == "pearson") {
      dist <- cor(t(object@data[genes, ]), method = "pearson")
      dist <- 1 - dist
    }
    
    if (dist_method == "cosine") {
      dist <- as.matrix(qlcMatrix::cosSparse(t(object@data[genes, ])))
      dist <- 1 - dist
    }
    
    if (dist_method == "euclidean") {
      dist <- as.matrix(dist(object@data[genes, ],
                             method = "euclidean",
                             upper = TRUE,
                             diag = TRUE
      ))
      diag(dist) <- NA
    }
    
    if (dist_method == "kendall") {
      dist <- as.matrix(amap::Dist(object@data[genes, ],
                                   method = "kendall"
      ))
    }
    
    
    # Compute mean correlation for each gene in cluster i
    dist_means <- colMeans(dist, na.rm = TRUE)
    dist_means <- dist_means[order(dist_means, decreasing = FALSE)]
    
    # Extract top genes with the highest correlation mean
    genes_top <- rbind(genes_top, names(dist_means[1:top]))
  }
  
  # Prepare top gene matrix
  if (length(cluster) > 1) {
    genes_top <- genes_top[2:(length(cluster) + 1), ]
  } else {
    genes_top <- as.matrix(t(genes_top[2:(length(cluster) + 1), ]))
  }
  
  colnames(genes_top) <- paste0("gene_top_", 1:top)
  rownames(genes_top) <- paste0("cluster_", cluster)
  
  # Put top gene matrix in object@top_genes
  object@top_genes <- split(x = unname(genes_top), f = cluster)
  
  for (clust in cluster) {
    object@top_genes[[clust]] <- object@top_genes[[as.character(clust)]][!is.na(object@top_genes[[as.character(clust)]])]
  }
  
  # Print
  print_msg(
    msg = "Results are stored in top_genes slot of ClusterSet object.",
    msg_type = "INFO"
  )
  
  return(object)
}
