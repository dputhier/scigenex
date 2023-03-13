################################################################################
filter_gene_clusters <- function(object = NULL,
                                 min_cluster_size = 5) {
  
}








################################################################################
filter_cluster_size <- function(object = NULL,
                                min_cluster_size = 5) {
  
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
  
  ## Update gene_cluster_metadata slots
  object@gene_clusters_metadata <- list("cluster_id" = as.numeric(names(object@gene_clusters)),
                                        "number" = max(as.numeric(names(object@gene_clusters))),
                                        "size" = unlist(lapply(object@gene_clusters, length)))
  
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
  nb_cluster_out <- length(obj@gene_clusters) - length(cluster_to_keep)
  print_msg(
    paste(
      nb_cluster_out,
      "clusters with less than", min_cluster_size, "genes are filtered out."
    ),
    msg_type = "INFO"
  )
  
  return(object)
}




################################################################################
filter_nb_supporting_cells <- function(object = NULL,
                                       min_nb_supporting_cell = 3,
                                       min_pct_gene_expressed = 50) {
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
  
  ## Update gene_cluster_metadata slots
  object@gene_clusters_metadata <- list("cluster_id" = as.numeric(names(object@gene_clusters)),
                                        "number" = max(as.numeric(names(object@gene_clusters))),
                                        "size" = unlist(lapply(object@gene_clusters, length)))
  
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
  nb_cluster_out <- length(obj@gene_clusters) - length(cluster_to_keep)
  print_msg(
    paste0(
      nb_cluster_out,
      " clusters with less than ", min_nb_supporting_cell, " cells expressing at least ", min_pct_gene_expressed, "% of genes are filtered out."
    ),
    msg_type = "INFO"
  )
  
  return(object)
}





################################################################################
filter_av_dot_prod <- function(object = NULL) {
  
}

filter_manual <- function(object = NULL) {
  
}


