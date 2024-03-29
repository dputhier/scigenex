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
#' @keywords internal
setGeneric("top_genes", 
           function(object,
                    top = 20,
                    cluster = "all")
             standardGeneric("top_genes")
)

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
#' # Set verbosity to 1 to display info messages only.
#' set_verbosity(1)
#' 
#' # Load a dataset
#' load_example_dataset('7871581/files/pbmc3k_medium_clusters')
#' 
#' # Store top genes in the object
#' pbmc3k_medium_clusters <- top_genes(pbmc3k_medium_clusters)
#'
setMethod("top_genes", 
          signature("ClusterSet"), 
              function(object,
                      top = 20,
                      cluster = "all") {
  
  ## Check format object arg
  check_format_cluster_set(object)
  
  if (unique(cluster == "all")) {
    cluster <- object@gene_clusters_metadata$cluster_id
  }
  
  # Display an info message if there is less than n top genes in a gene cluster
  loop <- 0
  for (size in object@gene_clusters_metadata$size[cluster]) {
    loop <- loop + 1

    if (top > size) {
      print_msg(paste0("Number of top genes is greater than the number of genes in cluster ", 
                       loop, 
                       ". All genes will be used."),
                msg_type = "INFO")
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
    
    if (dist_method %in% c("kendall", "pearson", "spearman")) {
      dist <- cor(t(object@data[genes, ]), method = dist_method)
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
})
