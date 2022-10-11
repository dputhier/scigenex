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
#' @return 
#' @export get_genes
#'
#' @examples
#' 
#' m <- matrix(rnorm(80000), nc=20)
#' m[1:100,1:10] <- m[1:100,1:10] + 4
#' m[101:200,11:20] <- m[101:200,11:20] + 3
#' m[201:300,5:15] <- m[201:300,5:15] + -2
#' 
#' res <- find_gene_clusters(data=m,
#'                           distance_method="pearson",
#'                           av_dot_prod_min = 0,
#'                           inflation = 1.2,
#'                           k=25,
#'                           fdr = 10)
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
    cluster <- c(1:length(object@size))
  }

  if(top){
    if (nrow(object@top_genes) == 1 &
        ncol(object@top_genes) == 1 &
        is.na(object@top_genes[1,1])) {
      stop(paste0("The slot top_genes of the input ClusterSet object is empty. Be sure to run top_genes() before."))
    }
    gene_names <- as.vector(t(object@top_genes[paste0("cluster_", cluster),]))
  }else{
    gene_names <- names(object@gene_patterns[object@gene_patterns %in% cluster])
  }
  
  return(gene_names)
}
