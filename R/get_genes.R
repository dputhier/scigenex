#################################################################
##    Define get_genes function for ClusterSet object
#################################################################

#' @title
#' Extract genes list from each gene clusters.
#' @description
#' Extract genes list from each gene clusters. Can extract genes from all clusters (default) or from specified clusters.
#' @param object A \code{ClusterSet} object.
#' @param cluster A vector of gene cluster identity (all by default).
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
#' res <- DBFMCL(data=m,
#'               distance_method="pearson",
#'               av_dot_prod_min = 0,
#'               inflation = 1.2,
#'               k=25,
#'               fdr = 10)
#'               
#'genes <- get_genes(object = res, cluster="1")
#'genes
#' 
#' \dontrun{
#' ## with an artificial dataset
#'
#' m <- matrix(rnorm(80000), nc = 20)
#' res <- get_data_4_DBFMCL(data=m)
#' }
#' 

get_genes <- function(object,
                      cluster = "all") {
  
  if(unique(cluster == "all")) {
    cluster <- c(1:length(object@size))
  }
  
  gene_names <- names(object@gene_patterns[object@gene_patterns %in% cluster])
  
  return(gene_names)
}
