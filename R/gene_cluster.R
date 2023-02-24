#################################################################
##    Define get_cells function for ClusterSet object
#################################################################

#' @title
#' Extract gene clusters from a cluster set object
#' @description
#' Extract gene clusters from a cluster set object.
#' @param object A \code{ClusterSet} object.
#' @param cluster A vector of gene cluster number (e.g. c(1,3))
#'
#' @return A named character vector.
#'
#' @examples
#'
#' set_verbosity(0)
#' data("complex9Noisy")
#' res <- find_gene_clusters(data=complex9Noisy[ ,1:2],
#'                           distance_method="euclidean",
#'                           inflation = 1.25,
#'                           k=20,
#'                           min_nb_supporting_cell=0,
#'                           highest = 0.95,
#'                           row_sum = -Inf,
#'                           fdr = 0.0001)
#'
#' plot(res@data, col=gene_cluster(res))
#'
#' @export gene_cluster
gene_cluster <- function(object,
                         cluster = 0,
                         core = FALSE) {
  
  if(!is.null(object@gene_clusters)){
    nb_clust <- length(object@gene_clusters)
  }else{
    print_msg("There is no cluster in this object.",
              msg_type = 'STOP')
  }
  
  if(!is.numeric(cluster))
    print_msg("Please provide a numeric value.",
              msg_type = 'STOP')
  
  cluster <- unique(cluster)
  
  if (!all(cluster-floor(cluster)==0) | any(cluster < 0 | any(cluster > nb_clust)))
    print_msg("Please provide a zero (all clusters) or positive integer in the required range.",
              msg_type = 'STOP')
  
  
  if(length(cluster) == 1){
    if (cluster == 0)
      cluster <- 1:length(object@gene_clusters)
  }

  if(length(cluster) > 1){
    if(length(cluster[cluster == 0]))
      print_msg("Zero is out of range.",
                msg_type = 'STOP')
  }
  
  
  if (nb_clust) {
    cluster_as_int <- unlist(mapply(rep,
                                    cluster,
                                    lapply(object@gene_clusters[cluster], length),
                                    SIMPLIFY = TRUE))
    cluster_as_int <- as.vector(as.matrix(cluster_as_int))
    names(cluster_as_int) <-
      unlist(object@gene_clusters[cluster])
    return(cluster_as_int)
    
  } else{
    return(NULL)
  }
}
