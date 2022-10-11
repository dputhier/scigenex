#################################################################
##    REDEFINE dim/ncol/nrow METHOD FOR CLASS OBJECT : ClusterSet
#################################################################

#' @title
#' ncol
#' @description
#' The number of columns of a ClusterSet object.
#' @param x A ClusterSet object.
#' @return The number of columns.
#' @export
#'
#' @examples
#'
#' \dontrun{
#'   m <- matrix(rnorm(80000), nc=20)
#'    m[1:100,1:10] <- m[1:100,1:10] + 4
#'    m[101:200,11:20] <- m[101:200,11:20] + 3
#'    m[201:300,5:15] <- m[201:300,5:15] + -2
#'    res <- find_gene_clusters(data=m,
#'                              distance_method="pearson",
#'                              av_dot_prod_min = 0,
#'                              inflation = 1.2,
#'                              k=25,
#'                              fdr = 10)
#'   ncol(res)
#' }
#'
setMethod(
  "ncol", signature("ClusterSet"),
  
  function(x) {
    
    return(ncol(x@data))
    
  }
)

#' @title
#' nrow
#' @description
#' The number of rows of a ClusterSet object.
#' @param x A ClusterSet object.
#' @return The number of rows.
#' @export
#'
#' @examples
#'
#' \dontrun{
#'   m <- matrix(rnorm(80000), nc=20)
#'    m[1:100,1:10] <- m[1:100,1:10] + 4
#'    m[101:200,11:20] <- m[101:200,11:20] + 3
#'    m[201:300,5:15] <- m[201:300,5:15] + -2
#'    res <- find_gene_clusters(data=m,
#'                              distance_method="pearson",
#'                              av_dot_prod_min = 0,
#'                              inflation = 1.2,
#'                              k=25,
#'                              fdr = 10)
#'   nrow(res)
#' }
#'
setMethod(
  "nrow", signature("ClusterSet"),
  
  function(x) {
    
    return(nrow(x@data))
    
  }
)


#' @title
#' dim
#' @description
#' The number of rows/columns of a ClusterSet object.
#' @param x A ClusterSet object.
#' @return The number of rows/columns.
#' @export
#'
#' @examples
#'
#' \dontrun{
#'   m <- matrix(rnorm(80000), nc=20)
#'    m[1:100,1:10] <- m[1:100,1:10] + 4
#'    m[101:200,11:20] <- m[101:200,11:20] + 3
#'    m[201:300,5:15] <- m[201:300,5:15] + -2
#'   res <- find_gene_clusters(data=m,
#'                             distance_method="pearson",
#'                             av_dot_prod_min = 0,
#'                             inflation = 1.2,
#'                             k=25,
#'                             fdr = 10)
#'   dim(res)
#' }
#'
setMethod(
  "dim", signature("ClusterSet"),
  
  function(x) {
    
    return(dim(x@data))
    
  }
)