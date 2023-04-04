#################################################################
##    Define get_genes function for ClusterSet object
#################################################################

#' @title
#' Extract genes list from each gene clusters.
#' @description
#' Extract genes list from each gene clusters. Can extract genes from all clusters (default) or from specified clusters.
#' @param object A \code{ClusterSet} object.
#' @param cluster A vector of gene cluster identity (all by default).
#' @param top A logical to indicate whether to provide top genes.
#'
#' @return A character vector.
#' @export get_genes
#'
#' @examples
#' # Set verbosity to 1 to display info messages only.
#' set_verbosity(1)
#' 
#' # Create a matrix with 4 signatures
#' m <- create_4_rnd_clust()
#' 
#' # Select informative genes
#' clust_set <- select_genes(m,
#'                           distance = "kendall",
#'                           k = 75,
#'                           highest = 0.3,
#'                          fdr = 1e-8,
#'                          row_sum = -Inf)
#'                     
#' # Cluster informative features
#' clust_set <- gene_clustering(clust_set, 
#'                             inflation = 1.2,
#'                             keep_nn = FALSE,
#'                             k = 5)
#' 
#' # Get all selected genes
#' get_genes(clust_set)
#' 
#' # Get genes from gene cluster 1
#' get_genes(clust_set, cluster = 1)
#' 
#' # Get top 5 genes from cluster 1
#' clust_set <- top_genes(clust_set, top = 5)
#' get_genes(clust_set, cluster = 1, top = TRUE)
#' 
get_genes <- function(object,
                      cluster = "all",
                      top = FALSE) {
  
  ## Check format object arg
  check_format_cluster_set(object)
  
  if(unique(cluster == "all")) {
    cluster <- object@gene_clusters_metadata$cluster_id
  }

  if(top){
    if (length(object@top_genes) == 0) {
      stop(paste0("The slot top_genes of the input ClusterSet object is empty. Be sure to run top_genes() before."))
    }
    gene_names <- unlist(object@top_genes[cluster], use.names = FALSE)
  }else{
    gene_names <- unlist(object@gene_clusters[cluster], use.names = FALSE)
  }
  
  return(gene_names)
}
