#' Filter cluster from a ClusterSet object using dot product. 
#' 
#' This function filters clusters of gene expression data based 
#' on their dot products. It aims to remove clusters that have a 
#' lot of zeros and are supported by only a few cells or spots. 
#' To do this, the function first converts the gene expression 
#' data for each cluster into a binary form (values greater than 
#' 1 are set to 1). Then it calculates the dot product for this 
#' binary matrix, which produces a gene-gene matrix showing the 
#' number of cells/spots where each pair of genes are expressed 
#' together. The function then calculates the median value of the 
#' maximum concordances across all genes, which can be used to 
#' determine whether a cluster should be filtered out or not.
#' 
#' @param data A ClusterSet object.
#' @param av_dot_prod_min Any cluster with average dot product below this value is discarded. This allow to delete
#' clusters in which correlation is influenced/supported by very few samples (typically 1).
#' @export filter_by_dot_prod
filter_by_dot_prod <- function(data = NULL,
                               av_dot_prod_min = 2) {
  
  if (is.null(data) | !inherits(data, "ClusterSet"))
    print_msg("Please provide a ClusterSet objet.", 
              msg_type = "STOP")
  
  if (!is.numeric(av_dot_prod_min) | av_dot_prod_min < 0)
    print_msg("The av_dot_prod_min argument should be a positive numeric value.", 
              msg_type = "STOP")

  selected_cluster <- names(res@gene_clusters)
  all_dot_prod <- vector()

  for (i in 1:length(data@gene_clusters)) {
    
    print_msg(paste0("Computing dot product for cluster: ", i), 
              msg_type = "DEBUG")
    cur_clust <- data@data[data@gene_clusters[[i]],]
    cur_clust[cur_clust > 0] <- 1
    cur_dot_prod <- cur_clust %*% t(cur_clust)
    diag(cur_dot_prod) <- NA
    cur_dot_prod_median_of_max <-median(apply(cur_dot_prod, 
                                              1, 
                                              max, 
                                              na.rm = T))
    all_dot_prod[i] <- cur_dot_prod_median_of_max
    
    # Dot product filtering
    if (cur_dot_prod_median_of_max <= av_dot_prod_min) {
      selected_cluster[i] <- NA
    }
  }

  selected_cluster <- selected_cluster[!is.na(selected_cluster)]
  print_msg(paste0("Number of selected clusters: ", length(selected_cluster)), 
            msg_type = "DEBUG")
  n_data <- data[selected_cluster, ]
  n_data@gene_clusters_metadata$dot_product <- all_dot_prod
  names(n_data@gene_clusters_metadata$dot_product) <- names(n_data@gene_clusters)
  return(data[selected_cluster, ])

}
