############################## filter_cluster_size ##############################
#' @title Filter out gene clusters containing few genes
#' @description This function filters out gene clusters from a ClusterSet object 
#' based on their size.
#' @param object A ClusterSet object.
#' @param min_cluster_size An integer indicating the minimum size for a clusters to be kept.
#' @param rename Should cluster be renamed.
#'
#' @return A ClusterSet object where clusters not passing the filter have been removed.
#' 
#' @examples 
#' # Set verbosity to 1 to display info messages only.
#' set_verbosity(1)
#' 
#' # Load a dataset
#' load_example_dataset('7871581/files/pbmc3k_medium_clusters')
#' clust_size(pbmc3k_medium_clusters)
#' # Remove the cluster with less than 20 genes
#' clust_set <- filter_cluster_size(pbmc3k_medium_clusters,
#'                                  min_cluster_size = 20)
#' clust_size(clust_set)
#' 
#' @export
filter_cluster_size <- function(object = NULL,
                                min_cluster_size = 5,
                                rename=FALSE) {
  
  ## Check format object arg
  check_format_cluster_set(object)
  
  # Store the initial number of clusters (used to compute the number of cluster filtered out)
  clust_names_before_filtering <- names(object@gene_clusters)
  
  cluster_to_keep <- lapply(object@gene_clusters, length) >= min_cluster_size
  
  print_msg(
    paste(
      sum(!cluster_to_keep),
      " clusters with less than", min_cluster_size, "were filtered out."
    ),
    msg_type = "INFO"
  )
  
  object <- object[cluster_to_keep, ]
  
  if(rename)
    object <- rename_clust(object)

  print_msg(paste0("Number of clusters left ", nclust(object)), msg_type = "INFO")
  
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
#' @param rename Should cluster be renamed.
#' 
#' @examples
#' load_example_dataset("7871581/files/pbmc3k_medium_clusters")
#' new_cs <- filter_nb_supporting_cells(pbmc3k_medium_clusters, 3, 50)
#' 
#' @return A ClusterSet object where clusters that did not pass the filter have been removed.
#'
#' @export
filter_nb_supporting_cells <- function(object = NULL,
                                       min_nb_supporting_cell = 3,
                                       min_pct_gene_expressed = 50,
                                       rename=TRUE) {
  
  if(is.null(object))
    print_msg('Please provide a valid object.')

  ## Check format object arg
  check_format_cluster_set(object)
  
  # Store the initial number of clusters (used to compute the number of cluster filtered out)
  clust_names_before_filtering <- names(object@gene_clusters)
  
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
  nb_cluster_out <- length(clust_names_before_filtering) - length(cluster_to_keep)
  print_msg(
    paste0(
      nb_cluster_out,
      " clusters with less than ", min_nb_supporting_cell, " cells expressing at least ", min_pct_gene_expressed, "% of genes were filtered out."
    ),
    msg_type = "INFO"
  )
  

  if(length(cluster_to_keep)==0){
    object <- object[-c(1:nclust(object)),]
  }else{
    object <- object[cluster_to_keep, ]
    if(rename)
      object <- rename_clust(object)
  }
  
  print_msg(paste0("Number of clusters left ", nclust(object)), msg_type = "INFO")
  
  return(object)
}


############################## filter_by_dot_prod ##############################

# median_of_max_dot_prod()
# Expect a binary matrix as input. Calculates the dot product for 
# this matrix, which produces a gene-gene matrix showing the 
# number of cells/spots where each pair of genes are expressed 
# together. The function then calculates the median value of the 
# maximum concordances across all genes, which can be used to 
# determine whether a cluster should be filtered out or not.
median_of_max_dot_prod <- function(cur_clust){
  cur_dot_prod <- cur_clust %*% t(cur_clust)
  diag(cur_dot_prod) <- NA
  cur_dot_prod_median_of_max <- median(apply(cur_dot_prod, 
                                            1, 
                                            max, 
                                            na.rm = T))
  return(cur_dot_prod_median_of_max)
}

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
#' @param rename Should cluster be renamed.
#' clusters in which correlation is influenced/supported by very few samples (typically 1).
#' 
#' @examples
#' load_example_dataset("7871581/files/pbmc3k_medium_clusters")
#' pbmc3k_medium_clusters <- top_genes(pbmc3k_medium_clusters)
#' nclust(pbmc3k_medium_clusters)
#' obj <- filter_by_dot_prod(pbmc3k_medium_clusters, av_dot_prod_min=5)
#' nclust(obj) 
#' @export
filter_by_dot_prod <- function(object = NULL, 
                               av_dot_prod_min = 2,
                               rename=TRUE) {

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
    all_dot_prod[i] <- median_of_max_dot_prod(cur_clust)
    
    # Dot product filtering
    if (all_dot_prod[i] <= av_dot_prod_min) {
      selected_cluster[i] <- NA
    }
  }
  
  # Extract number of filtered out clusters
  nb_cluster_out <- length(selected_cluster[is.na(selected_cluster)])
  
  # Extract selected clusters
  selected_cluster <- selected_cluster[!is.na(selected_cluster)]
  
  print_msg(paste0("Number of selected clusters: ", length(selected_cluster)), 
            msg_type = "DEBUG")
  
  
  # Update ClusterSet object
  ## Add dot product values
  object@gene_clusters_metadata$all_dot_prod <- all_dot_prod[selected_cluster]
  names(object@gene_clusters_metadata$all_dot_prod) <- object@gene_clusters_metadata$cluster_id[selected_cluster]

   
  if(length(selected_cluster)==0){
    return(object[-c(1:nclust(object)),])
  }else{
    print_msg(paste0("Selected clusters : ", selected_cluster), msg_type = "DEBUG")
    object <- object[selected_cluster, ]
    
    print_msg("Renaming.", msg_type = "DEBUG")
    if(rename)
      object <- rename_clust(object)
    return(object)
  }

}

############################## filter_cluster_sd ##############################
#' @title Filter out gene clusters with small standard deviation.
#' @description Filter out gene clusters with small standard deviation.
#' @param object A ClusterSet object.
#' @param min_sd An integer indicating the minimum standard deviation for a clusters to be kept.
#' @param rename Should cluster be renamed.
#'
#' @return A ClusterSet object where clusters not passing the filter have been removed.
#' @examples 
#' load_example_dataset("7871581/files/pbmc3k_medium_clusters")
#' df <- cluster_stats(pbmc3k_medium_clusters)  
#' plot_cluster_stats(df)
#' plot_cluster_stats(df, highlight=df$sd > 0.3) 
#' pbmc3k_medium_clusters_sub <- filter_cluster_sd(pbmc3k_medium_clusters, min_sd=0.3, rename=FALSE)
#' plot_cluster_stats(cluster_stats(pbmc3k_medium_clusters_sub))
#' @export
filter_cluster_sd <- function(object = NULL,
                              min_sd = 0.2,
                              rename=TRUE) {
  
  ## Check format object arg
  check_format_cluster_set(object)
  
  # Store the initial number of clusters (used to compute the number of cluster filtered out)
  clust_names_before_filtering <- names(object@gene_clusters)
  
  gene_clust <- object@gene_clusters
  sd_clust <- vector()

  for(pos in 1:length(gene_clust)){
    sd_clust[pos] <- stats::sd(object@data[gene_clust[[pos]], ]) 
  }

  cluster_to_keep <- sd_clust >= min_sd
  
  print_msg(
    paste(
      sum(!cluster_to_keep),
      " clusters with std dev lower than", min_sd, " will be filtered out."
    ),
    msg_type = "INFO"
  )
  
  if(length(cluster_to_keep)==0){
    object <- object[-c(1:nclust(object)),]
  }else{
    object <- object[cluster_to_keep, ]
    if(rename)
      object <- rename_clust(object)
  }
  
  print_msg(paste0("Number of clusters left ", nclust(object)), msg_type = "INFO")
  
  return(object)
}

