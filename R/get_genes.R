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
#' 
#' m <- create_3_rnd_clust()
#' 
#' res <- find_gene_clusters(data=m,
#'                              distance_method="pearson",
#'                              inflation = 2,
#'                              k=75,
#'                              row_sum=-Inf,
#'                              highest=0.3,
#'                              min_nb_supporting_cell = 0,
#'                              fdr = 1e-8)
#'               
#'genes <- get_genes(object = res, cluster="1")
#'genes
#' 
#' \dontrun{
#' ## with an artificial dataset
#'
#' m <- matrix(rnorm(80000), nc = 20)
#' res <- get_data_for_scigenex(data=m)
#' }
#' 

get_genes <- function(object,
                      cluster = "all",
                      top = FALSE) {
  
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
