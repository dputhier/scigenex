#' @title Compute statistics about the clusters
#' @description Compute statistics about the clusters
#' @param object a ClusterSet object.
#' @importFrom stats var
#' @details
#' For each cluster this function computes:
#' - The sum of (normalized) counts divided by the number of genes (cluster size). 
#' This gives a good estimates of whether the genes are weakly or highly expressed (*i.e* housekeeping)
#' - The variance. This provides an estimate of the dispertion of the values inside a cluster.
#' - The standard deviation. This provides an estimate of the dispertion of the values inside a cluster.
#' - The coefficient of variation. This provides an estimate of the dispertion of the values inside a cluster (normalized by the mean expression).
#' @examples
#' # load a dataset
#' load_example_dataset('7871581/files/pbmc3k_medium_clusters')
#' df <- cluster_stats(pbmc3k_medium_clusters)
#' @export cluster_stats
#' @keywords internal
setGeneric("cluster_stats", 
           function(object)
             standardGeneric("cluster_stats")
)

#' @title Compute statistics about the clusters
#' @description Compute statistics about the clusters
#' @param object a ClusterSet object.
#' @importFrom stats var
#' @details
#' For each cluster this function computes:
#' - The sum of (normalized) counts divided by the number of genes (cluster size). 
#' This gives a good estimates of whether the genes are weakly or highly expressed (*i.e* housekeeping)
#' - The variance. This provides an estimate of the dispertion of the values inside a cluster.
#' - The standard deviation. This provides an estimate of the dispertion of the values inside a cluster.
#' - The coefficient of variation. This provides an estimate of the dispertion of the values inside a cluster (normalized by the mean expression).
#' @examples
#' # load a dataset
#' load_example_dataset('7871581/files/pbmc3k_medium_clusters')
#' df <- cluster_stats(pbmc3k_medium_clusters)
#' @export cluster_stats
setMethod(
  "cluster_stats", signature("ClusterSet"),
  function(object) {
    #object <- scigenex::check_format_cluster_set(object)
    
    df <- data.frame(size=clust_size(object), 
                     row.names = names(clust_size(object)))
    
    print_msg('Splitting expression data', msg_type="DEBUG")
    gene_clust <- as.factor(gene_cluster(object))
    df_split <- split(as.data.frame(object@data), gene_clust)
    
    print_msg('Computing sum', msg_type="DEBUG")
    tmp <- unlist(lapply(df_split, sum)) 
    tmp_2 <-  clust_size(object)


    df$sum_by_row <- tmp[row.names(df)] / tmp_2[row.names(df)]
    
    print_msg('Computing sd.', msg_type="DEBUG")
    df_split <- split(object@data, gene_clust)
    
    print_msg('Computing var.', msg_type="DEBUG")
    tmp <- unlist(lapply(df_split, stats::var)) 
    df$var <- tmp[row.names(df)]
    

    tmp <- unlist(lapply(df_split, stats::sd))
    df$sd <- tmp[row.names(df)]
 
    print_msg('Computing cv.', msg_type="DEBUG")
    tmp <- unlist(lapply(df_split, stats::sd)) / unlist(lapply(df_split, mean))
    df$cv <- tmp[row.names(df)]
    
    return(df)
    
  }
)


#' @title Plot the result of cluster_stats()
#' @description Plot the result of cluster_stats()
#' @param x A \code{ClusterSet} object.
#' @param highlight a vector with two modalities indicating which clusters to highlight
#' @return A \code{ClusterSet} object.
#' @importFrom ggplot2 ggplot geom_col facet_grid coord_flip theme_bw
#' @importFrom reshape2 melt
#' @export plot_cluster_stats
#'
#' @examples
#' # Set verbosity to 1 to display info messages only.
#' set_verbosity(1)
#' 
#' # Load a dataset
#' load_example_dataset('7871581/files/pbmc3k_medium_clusters')
#' 
#'                        
#' # Compute some statistics about the clusters
#' df <- cluster_stats(pbmc3k_medium_clusters)  
#' plot_cluster_stats(df)
#' plot_cluster_stats(df, highlight=df$size > 8) 
#' 
#' # Select the cluster of interest
#' pbmc3k_medium_clusters_filtered <- pbmc3k_medium_clusters[df$size > 8, ] 
#' nclust(pbmc3k_medium_clusters_filtered)          
plot_cluster_stats <- function(x, highlight=NULL){
  
  if(!is.null(highlight)){
    if(length(highlight) != nrow(x))
      print_msg("The number of clusters to highlight should have same length as nrow(x).",
                msg_type = "STOP")
    
    highlight <- as.factor(highlight)
    
    if(length(levels(highlight)) != 2)
      print_msg("The 'highlight' vector should have two modalities.",
                msg_type = "STOP")
  }
  
  x <- reshape2::melt(as.matrix(x))
  colnames(x) <- c("cluster", "stat", "value")
  
  if(!is.null(highlight))
    x$highlight <- highlight 
  
  x$cluster <- factor(x$cluster, 
                      levels=sort(unique(x$cluster)),
                      ordered=T)
  x$stat <- factor(x$stat, 
                   levels=sort(unique(x$stat)),
                   ordered=T)
  if(!is.null(highlight)){
    ggplot2::ggplot(data=x, 
                    mapping=ggplot2::aes(x=cluster,
                                         y=value,
                                         fill=highlight)) +
      ggplot2::geom_col() + 
      ggplot2::facet_grid(~stat, scales = "free") +
      ggplot2::coord_flip() +
      ggplot2::theme_bw()
  }else{
    ggplot2::ggplot(data=x, 
                    mapping=ggplot2::aes(x=cluster,
                                         y=value)) +
      ggplot2::geom_col() + 
      ggplot2::facet_grid(~stat, scales = "free") +
      ggplot2::coord_flip()+
      ggplot2::theme_bw()
  }
}
