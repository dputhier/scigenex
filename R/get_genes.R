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
#'
#' @examples
#' load_example_dataset('7871581/files/pbmc3k_medium_clusters')
#' # Get all selected genes
#' x <- get_genes(pbmc3k_medium_clusters)
#' 
#' # Get genes from gene cluster 1
#' x<- get_genes(pbmc3k_medium_clusters, cluster = 1)
#' 
#' # Get top 5 genes from cluster 1
#' pbmc3k_medium_clusters <- top_genes(pbmc3k_medium_clusters, top = 5)
#' x <- get_genes(pbmc3k_medium_clusters, cluster = 1, top = TRUE)
#' @export
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
