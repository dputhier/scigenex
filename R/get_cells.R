#################################################################
##    Define get_cells function for ClusterSet object
#################################################################

#' @title
#' Extract cell identities from each cell clusters.
#' @description
#' Extract cell identities from each cell clusters. Can extract cells from all clusters (default) or from specified clusters.
#' @param object A \code{ClusterSet} object.
#' @param cluster A vector of cell cluster identity (all by default).
#' @param core A logical to indicate whether to provide core cells.
#'
#' @return 
#' @export get_cells
#'
#' @examples
#' \dontrun{
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
#' res <- cell_clust(res, min_cluster_size = 3)
#' cells <- get_cells(object = res, cluster="1")
#' cells
#' }
#' 

get_cells <- function(object,
                      cluster = "all",
                      core = FALSE) {
  
  if(unique(cluster == "all")) {
    cluster <- order(unique(object@cell_clusters$labels))
  }
  
  if(!core){
    cell_names <- names(object@cell_clusters$labels[object@cell_clusters$labels %in% cluster])
  } else {
    cell_names <- names(object@cell_clusters$cores[object@cell_clusters$cores %in% cluster])
  }
  
  
  return(cell_names)
}
