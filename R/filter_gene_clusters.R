############################## filter_cluster_size ##############################
#' @title Filter out gene clusters containing few genes
#' @description This function filters out gene clusters from a ClusterSet object 
#' based on their size.
#' @param object A ClusterSet object.
#' @param min_cluster_size An integer indicating the minimum size for a clusters to be kept.
#'
#' @return A ClusterSet object where clusters not passing the filter have been removed.
#' 
#' @examples 
#' # Set verbosity to 1 to display info messages only.
#' set_verbosity(1)
#' 
#' # Create a matrix with 4 signatures
#' m <- create_4_rnd_clust()
#' 
#' # Select informative genes
#' clust_set <- select_genes(m,
#'                          distance = "pearson",
#'                          k = 75,
#'                          highest = 0.3,
#'                          fdr = 1e-8,
#'                          row_sum = -Inf)
#'                     
#' # Cluster informative features
#' clust_set <- gene_clustering(clust_set, 
#'                             inflation = 1.2,
#'                             keep_nn = FALSE,
#'                             k = 5)
#' clust_size(clust_set)
#' 
#' # Remove the cluster with less than 100 genes
#' clust_set <- filter_cluster_size(clust_set,
#'                                  min_cluster_size = 100)
#' clust_size(clust_set)
#' 
#' @export filter_cluster_size

filter_cluster_size <- function(object = NULL,
                                min_cluster_size = 5) {
  
  ## Check format object arg
  check_format_cluster_set(object)
  
  # Store the initial number of clusters (used to compute the number of cluster filtered out)
  nb_clusters_before_filtering <- names(object@gene_clusters)
  
  cluster_to_keep <- object@gene_clusters_metadata$size > min_cluster_size
  
  print_msg(
    paste(
      sum(!cluster_to_keep),
      " clusters with less than", min_cluster_size, "genes will be filtered out."
    ),
    msg_type = "INFO"
  )
  
  object <- object[cluster_to_keep, ]
  object <- rename_clust(object)

  return(object)
}
  
############################## filter_nb_supporting_cells ##############################
#' @title Filter out cluster supported by few cells.
#' 
#' @description
#'  This function filters out the gene clusters based on the number of cells expressing a 
#' certain percentage of genes present in this gene cluster.
#' 
#' @param object A ClusterSet object.
#' @param min_nb_supporting_cell An integer indicating the minimum number of cell supporting a cluster.
#' A cell supports a cluster if it expresses at least min_pct_gene_expressed \% of the genes from the cluster.
#' @param min_pct_gene_expressed See min_nb_supporting_cell argument.
#' 
#' @examples
#' 
#' 
#' @return A ClusterSet object where clusters that did not pass the filter have been removed.
#'
#' @export filter_nb_supporting_cells

filter_nb_supporting_cells <- function(object = NULL,
                                       min_nb_supporting_cell = 3,
                                       min_pct_gene_expressed = 50) {
  
  if(is.null(object))
    print_msg('Please provide a valid object.')

  ## Check format object arg
  check_format_cluster_set(object)
  
  # Store the initial number of clusters (used to compute the number of cluster filtered out)
  nb_clusters_before_filtering <- names(object@gene_clusters)
  
  cluster_to_keep <- vector()
  
  for (i in seq_along(object@gene_clusters)) {
    
    cur_clust <- object@data[object@gene_clusters[[i]], , drop=FALSE]
    
    # Compute the percentage of genes expressed per cell
    cur_clust[cur_clust < 1 ] <- 0
    cur_clust[cur_clust > 1 ] <- 1
    
    pct_gene_expressed_cur <- apply(cur_clust, 2, sum) / nrow(cur_clust) * 100

    min_nb_supporting_cell_cur <- length(pct_gene_expressed_cur[pct_gene_expressed_cur >= min_pct_gene_expressed])
    
    if (min_nb_supporting_cell_cur >= min_nb_supporting_cell) {
      cluster_to_keep <- append(cluster_to_keep, i)
    }
  }
  
  # Print number of cluster filtered out
  nb_cluster_out <- length(nb_clusters_before_filtering) - length(cluster_to_keep)
  print_msg(
    paste0(
      nb_cluster_out,
      " clusters with less than ", min_nb_supporting_cell, " cells expressing at least ", min_pct_gene_expressed, "% of genes were filtered out."
    ),
    msg_type = "INFO"
  )
  
  # Return NULL if all clusters were filtered out
  if(length(cluster_to_keep)==0){
    return(NULL)
  }else{
    object <- object[cluster_to_keep, ]
    object <- rename_clust(object)
    return(object)
  }
}




############################## filter_by_dot_prod ##############################
#' @title Filter cluster from a ClusterSet object using dot product. 
#' 
#' @description
#'  This function filters clusters of gene expression data based 
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
#' @param object A ClusterSet object.
#' @param av_dot_prod_min Any cluster with average dot product below this value is discarded. This allow to delete
#' clusters in which correlation is influenced/supported by very few samples (typically 1).
#' 
#' @examples
#' # Load a Seurat object
#' data(pbmc_small, package = "SeuratObject")
#' 
#' # Compute the signatures using find_gene_clusters()
#' clust_set <- select_genes(pbmc_small,
#'                           distance = "kendall",
#'                           k = 75,
#'                           highest = 0.3,
#'                           fdr = 1e-8,
#'                           row_sum = -Inf,
#'                           no_dknn_filter = TRUE)
#'                           
#' # Cluster informative features
#' clust_set <- gene_clustering(clust_set, 
#'                             inflation = 1.6,
#'                             keep_nn = FALSE,
#'                             k = 5)
#'                             
#' # Plot heatmap of gene clusters
#' plot_heatmap(clust_set)
#' 
#' # Remove the cluster 5 containing less than 4 cells expressing 100% of the genes in cluster 4
#' clust_set <- filter_by_dot_prod(clust_set,
#'                                 av_dot_prod_min = 5)
#' plot_heatmap(clust_set)
#' 
#' @export filter_by_dot_prod
filter_by_dot_prod <- function(object = NULL,
                               av_dot_prod_min = 2) {
  
  ## Check format object arg
  check_format_cluster_set(object)
  
  if (is.null(object) | !inherits(object, "ClusterSet"))
    print_msg("Please provide a ClusterSet objet.", 
              msg_type = "STOP")
  
  if (!is.numeric(av_dot_prod_min) | av_dot_prod_min < 0)
    print_msg("The av_dot_prod_min argument should be a positive numeric value.", 
              msg_type = "STOP")
  
  selected_cluster <- names(object@gene_clusters)
  all_dot_prod <- vector()
  
  for (i in 1:length(object@gene_clusters)) {
    
    print_msg(paste0("Computing dot product for cluster: ", i), 
              msg_type = "DEBUG")
    cur_clust <- object@data[object@gene_clusters[[i]],]
    cur_clust[cur_clust > 1] <- 1
    cur_clust[cur_clust < 1] <- 0
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
  
  # Extract number of filtered out clusters
  nb_cluster_out <- length(selected_cluster[is.na(selected_cluster)])
  
  # Extract selected clusters
  selected_cluster <- selected_cluster[!is.na(selected_cluster)]
  
  # Stop if all the clusters are filtered out
  if(length(selected_cluster) == 0){
    print_msg("No clusters conserved.", msg_type = "STOP")
  }
  
  print_msg(paste0(nb_cluster_out, " clusters are filtered out based on their dot product."),
            msg_type = "INFO")
  print_msg(paste0("Number of selected clusters: ", length(selected_cluster)), 
            msg_type = "DEBUG")
  
  
  # Update ClusterSet object
  ## Add dot product values
  object@gene_clusters_metadata$all_dot_prod <- all_dot_prod
  names(object@gene_clusters_metadata$all_dot_prod) <- object@gene_clusters_metadata$cluster_id
  
  ## Remove clusters that did not pass the filter
  object@gene_clusters <- object@gene_clusters[selected_cluster]
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
  
  return(object)
}

############################## filter_cluster_sd ##############################
#' @title Filter out gene clusters with small standard deviation.
#' @description Filter out gene clusters with small standard deviation.
#' @param object A ClusterSet object.
#' @param min_sd An integer indicating the minimum standard deviation for a clusters to be kept.
#'
#' @return A ClusterSet object where clusters not passing the filter have been removed.
#' 
#' @examples 
#' # Set verbosity to 1 to display info messages only.
#' set_verbosity(1)
filter_cluster_sd <- function(object = NULL,
                              min_sd = 0.1) {
  
  ## Check format object arg
  check_format_cluster_set(object)
  
  # Store the initial number of clusters (used to compute the number of cluster filtered out)
  nb_clusters_before_filtering <- names(object@gene_clusters)
  
  gene_clust <- as.factor(gene_cluster(object))
  df_split <- split(object@data, gene_clust)
  sd_total <- unlist(lapply(df_split, sd))

  cluster_to_keep <- sd_total < min_sd
  
  print_msg(
    paste(
      sum(!cluster_to_keep),
      " clusters with std dev lower than", min_sd, " will be filtered out."
    ),
    msg_type = "INFO"
  )
  
  object <- object[cluster_to_keep, ]
  object <- rename_clust(object)
  
  return(object)
}


############################## filter_manual ##############################
# filter_manual <- function(object = NULL) {
#   
# }


#' ############################## filters all fct ##############################
#' #' Filter out gene clusters
#' #'
#' #' This function filters out gene clusters from a ClusterSet object
#' #' 
#' #' @param object A ClusterSet object.
#' #' @param min_cluster_size An integer indicating the minimum size of clusters to keep.
#' #' @param min_nb_supporting_cell An integer indicating the minimum number of cell supporting a cluster.
#' #' A cell supports a cluster if it expresses at least min_pct_gene_expressed \% of the genes from the cluster.
#' #' @param min_pct_gene_expressed See min_nb_supporting_cell argument.
#' #' 
#' #' @return A ClusterSet object where clusters that did not pass the filter have been removed.
#' #' 
#' #' @export filter_gene_clusters
#' 
#' filter_gene_clusters <- function(object = NULL,
#'                                  min_cluster_size = 5,
#'                                  min_nb_supporting_cell = 3,
#'                                  min_pct_gene_expressed = 50) {
#'   
#'   object <- filter_cluster_size(object = object,
#'                                 min_cluster_size = min_cluster_size)
#'   
#'   object <- filter_nb_supporting_cells(object = object,
#'                                        min_nb_supporting_cell = min_nb_supporting_cell,
#'                                        min_pct_gene_expressed = min_pct_gene_expressed)
#'   
#'   return(object)
#' }