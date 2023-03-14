############################## filters all fct ##############################
#' Filter out gene clusters
#'
#' This function filters out gene clusters from a ClusterSet object
#' 
#' @param object A ClusterSet object.
#' @param min_cluster_size An integer indicating the minimum size of clusters to keep.
#' @param min_nb_supporting_cell An integer indicating the minimum number of cell supporting a cluster.
#' A cell supports a cluster if it expresses at least min_pct_gene_expressed \% of the genes from the cluster.
#' @param min_pct_gene_expressed See min_nb_supporting_cell argument.
#' 
#' @return A ClusterSet object where clusters that did not pass the filter have been removed.
#' 
#' @export filter_gene_clusters

filter_gene_clusters <- function(object = NULL,
                                 min_cluster_size = 5,
                                 min_nb_supporting_cell = 3,
                                 min_pct_gene_expressed = 50) {
  
  object <- filter_cluster_size(object = object,
                                min_cluster_size = min_cluster_size)
  
  object <- filter_nb_supporting_cells(object = object,
                                       min_nb_supporting_cell = min_nb_supporting_cell,
                                       min_pct_gene_expressed = min_pct_gene_expressed)
  
  return(object)
}





############################## filter_cluster_size ##############################
#' Filter out gene clusters based on their size
#'
#' This function filters out gene clusters from a ClusterSet object based on their size.
#'
#' @param object A ClusterSet object.
#' @param min_cluster_size An integer indicating the minimum size of clusters to keep.
#'
#' @return A ClusterSet object where clusters that did not pass the filter have been removed.
#'
#' @export filter_cluster_size

filter_cluster_size <- function(object = NULL,
                                min_cluster_size = 5) {
  
  # Store the initial number of clusters (used to compute the number of cluster filtered out)
  nb_clusters_before_filtering <- names(object@gene_clusters)
  
  cluster_to_keep <- c()
  
  for (i in seq_along(object@gene_clusters)) {
    gene_in_clust <- object@gene_clusters[[i]]
    
    if (length(gene_in_clust) > min_cluster_size) {
      cluster_to_keep <- c(cluster_to_keep, i)
    }
  }
  
  # Update ClusterSet object
  ## Remove clusters that did not pass the filter
  object@gene_clusters <- object@gene_clusters[cluster_to_keep]
  names(object@gene_clusters) <- seq(1, length(object@gene_clusters))
  
  ## Update gene_cluster_metadata slots
  object@gene_clusters_metadata <- list("cluster_id" = as.numeric(names(object@gene_clusters)),
                                        "number" = max(as.numeric(names(object@gene_clusters))),
                                        "size" = unlist(lapply(object@gene_clusters, length)))

  ## Remove filtered out genes in data slot
  object@data <- object@data[unlist(object@gene_clusters, use.names = FALSE),]
  
  ## Compute centers
  nb_clusters = length(names(object@gene_clusters))
  centers <- matrix(ncol = ncol(object@data),
                    nrow = nb_clusters)
  colnames(centers) <- colnames(object@data)
  rownames(centers) <- names(object@gene_clusters)
  
  ## calcul of the mean profils
  for (i in 1:nb_clusters) {
    if(is(object@data[object@gene_clusters[[i]], ])[2] == "vector") {
      centers[i, ] <- apply(t(as.matrix(object@data[object@gene_clusters[[i]], ])),
                            2, mean,
                            na.rm = TRUE)
    } else {
      centers[i, ] <- apply(object@data[object@gene_clusters[[i]], ],
                            2, mean,
                            na.rm = TRUE)
    }
  }
  
  object@dbf_output$center <- centers
  
  # Print number of cluster filtered out
  nb_cluster_out <- length(nb_clusters_before_filtering) - length(cluster_to_keep)
  print_msg(
    paste(
      nb_cluster_out,
      "clusters with less than", min_cluster_size, "genes are filtered out."
    ),
    msg_type = "INFO"
  )
  
  return(object)
}




############################## filter_nb_supporting_cells ##############################
#' Filter out gene clusters based on the number of cells expressing at least a provided percentage of genes
#' 
#' This function filters out the gene clusters based on the number of cells expressing a certain percentage of genes present in this gene cluster.
#' 
#' @param min_nb_supporting_cell An integer indicating the minimum number of cell supporting a cluster.
#' A cell supports a cluster if it expresses at least min_pct_gene_expressed \% of the genes from the cluster.
#' @param min_pct_gene_expressed See min_nb_supporting_cell argument.
#' 
#' @return A ClusterSet object where clusters that did not pass the filter have been removed.
#'
#' @export filter_nb_supporting_cells

filter_nb_supporting_cells <- function(object = NULL,
                                       min_nb_supporting_cell = 3,
                                       min_pct_gene_expressed = 50) {
  
  # Store the initial number of clusters (used to compute the number of cluster filtered out)
  nb_clusters_before_filtering <- names(object@gene_clusters)
  
  cluster_to_keep <- c()
  
  for (i in seq_along(object@gene_clusters)) {
    cur_clust <- object@data[object@gene_clusters[[i]], ]
    cur_clust[cur_clust > 0] <- 1
    min_pct_gene_expressed_cur <- apply(cur_clust, 2, sum) / nrow(cur_clust) * 100
    min_nb_supporting_cell_cur <- length(min_pct_gene_expressed_cur[min_pct_gene_expressed_cur >= min_pct_gene_expressed])
    
    if (min_nb_supporting_cell_cur >= min_nb_supporting_cell) {
      cluster_to_keep <- c(cluster_to_keep, i)
    }
  }
  
  # Update ClusterSet object
  ## Remove clusters that did not pass the filter
  object@gene_clusters <- object@gene_clusters[cluster_to_keep]
  names(object@gene_clusters) <- seq(1, length(object@gene_clusters))
  
  ## Update gene_cluster_metadata slots
  object@gene_clusters_metadata <- list("cluster_id" = as.numeric(names(object@gene_clusters)),
                                        "number" = max(as.numeric(names(object@gene_clusters))),
                                        "size" = unlist(lapply(object@gene_clusters, length)))
  
  ## Remove filtered out genes in data slot
  object@data <- object@data[unlist(object@gene_clusters, use.names = FALSE),]
  
  ## Compute centers
  nb_clusters = length(names(object@gene_clusters))
  centers <- matrix(ncol = ncol(object@data),
                    nrow = nb_clusters)
  colnames(centers) <- colnames(object@data)
  rownames(centers) <- names(object@gene_clusters)
  
  ## calcul of the mean profils
  for (i in 1:nb_clusters) {
    if(is(object@data[object@gene_clusters[[i]], ])[2] == "vector") {
      centers[i, ] <- apply(t(as.matrix(object@data[object@gene_clusters[[i]], ])),
                            2, mean,
                            na.rm = TRUE)
    } else {
      centers[i, ] <- apply(object@data[object@gene_clusters[[i]], ],
                            2, mean,
                            na.rm = TRUE)
    }
  }
  
  object@dbf_output$center <- centers
  
  # Print number of cluster filtered out
  nb_cluster_out <- length(nb_clusters_before_filtering) - length(cluster_to_keep)
  print_msg(
    paste0(
      nb_cluster_out,
      " clusters with less than ", min_nb_supporting_cell, " cells expressing at least ", min_pct_gene_expressed, "% of genes are filtered out."
    ),
    msg_type = "INFO"
  )
  
  return(object)
}




############################## filter_av_dot_prod ##############################
# filter_av_dot_prod <- function(object = NULL) {
#   
# }


############################## filter_manual ##############################
# filter_manual <- function(object = NULL) {
#   
# }


