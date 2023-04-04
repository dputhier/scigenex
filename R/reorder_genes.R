#################################################################
##    Define the reorder_genes function
#################################################################

#' @title
#' reorder_genes
#' @description
#' Reorder genes in gene cluster.
#' @param object A ClusterSet object.
#' @param order_by A character string related to the reordering type. 
#' Default is "correlation" which reorder genes based on their correlation 
#' within their respective gene cluster.
#'
#' @return ClusterSet-class object
#' @export
#'
#' @examples
#' # Set verbosity to 1 to display info messages only.
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
#' # Reorder genes using hierarchical clustering
#' res <- reorder_genes(res, order_by = "hclust")
#' 
#' # Plot heatmap of gene clusters
#' plot_heatmap(res)
#' 
#' @rdname reorder_genes

reorder_genes <- function(object,
                          order_by = c("gene_names", "hclust", "correlation")) {
  
  ## Check format object arg
  check_format_cluster_set(object)
  
  if (!order_by %in% c("gene_names", "hclust", "correlation")) {
    stop("Provided order_by parameter does not exist.")
  }
  
  
  reordered_genes <- c()
  for (gene_cluster in object@gene_clusters_metadata$cluster_id) {
    
    # Reorder by gene names
    if (order_by == "gene_names") {
      gene_cluster_names <- object@gene_clusters[[gene_cluster]]
      gene_cluster_names <- sort(gene_cluster_names)
      
      # Reorder genes in gene_clusters slot
      object@gene_clusters[[gene_cluster]] <- gene_cluster_names
      
      # Reorder genes in data slot
      object@data <- object@data[unlist(object@gene_clusters, use.names = FALSE),]
    }
    
    # Reorder using hierarchical clustering (amap::hcluster)
    if (order_by == "hclust") {
      set.seed(123)
      gene_cluster_names <- object@gene_clusters[[gene_cluster]]
      expression_matrix <- object@data[gene_cluster_names,]
      
      
      hclust_res <- hcluster(expression_matrix,
                             method = "pearson",
                             diag = FALSE,
                             upper = FALSE,
                             link = "average",
                             nbproc = 2,
                             doubleprecision = TRUE)
      
      # Reorder genes in gene_clusters slot
      object@gene_clusters[[gene_cluster]] <- hclust_res$labels[hclust_res$order]
      
      # Reorder genes in data slot
      object@data <- object@data[unlist(object@gene_clusters, use.names = FALSE),]
    }
    
    # Reorder using top_genes
    if (order_by == "correlation") {
      gene_cluster_names <- object@gene_clusters[[gene_cluster]]
      gene_cluster_names <- top_genes(object,
                                      cluster = gene_cluster,
                                      top = length(gene_cluster_names))@top_genes[[gene_cluster]]
      
      # Reorder genes in gene_clusters slot
      object@gene_clusters[[gene_cluster]] <- gene_cluster_names
      
      # Reorder genes in data slot
      object@data <- object@data[unlist(object@gene_clusters, use.names = FALSE),]
    }
  }
  
  
  return(object)
}
