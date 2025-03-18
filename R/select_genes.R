################################################################################
#' Selects informative genes based on k-nearest neighbour analysis.
#'
#' This function selects genes based on k-nearest neighbour analysis.
#' The function takes a seurat object or gene expression matrix as input 
#' and compute distance to k-nearest neighbour for each gene/feature.
#' A threshold is set based on permutation analysis and FDR computation. 
#'
#' @param data A matrix, data.frame or Seurat object.
#' @param distance_method a character string indicating the method for computing distances (one of "pearson", "cosine", 
#' "euclidean", spearman or "kendall").
#' @param noise_level This parameter controls the fraction of genes with high dknn (ie. noise) whose neighborhood (i.e associated distances) 
#' will be used to compute simulated DKNN values. A value of 0 means to use all the genes. A value close to 1 means  to use only gene 
#' with high dknn (i.e close to noise).
#' @param k An integer specifying the size of the neighborhood.
#' @param row_sum A feature/gene whose row sum is below this threshold will be discarded. Use -Inf to keep all genes. 
#' @param fdr A numeric value indicating the false discovery rate threshold (range: 0 to 100).
#' @param which_slot a character string indicating which slot to use from the input scRNA-seq object (one of "data", "sct" or "counts") to compute distances. 
#' @param no_dknn_filter a logical indicating whether to skip the k-nearest-neighbors (KNN) filter. If FALSE, all genes are kept for the next steps.
#' @param no_anti_cor If TRUE, correlation below 0 are set to zero ("pearson", "cosine", "spearman" "kendall"). This may increase the 
#' relative weight of positive correlation (as true anti-correlation may be rare).
#' @param final_slot This is the slot to export in the final object. In case distance is set to binary, the function will force the usage of 'counts' for distance
#' computation. However, the final object may contain normalized data (i.e log transformed) extract from the 'data' slot.
#' @param seed An integer specifying the random seed to use.
#'
#' @return a ClusterSet class object
#' 
#' @author Julie Bavais, Sebastien Nin, Lionel Spinelli and Denis Puthier
#' 
#' @references
#' - Lopez F.,Textoris J., Bergon A., Didier G., Remy E., Granjeaud
#' S., Imbert J. , Nguyen C. and Puthier D. TranscriptomeBrowser: a powerful
#' and flexible toolbox to explore productively the transcriptional landscape
#' of the Gene Expression Omnibus database. PLoSONE, 2008;3(12):e4001.
#'
#' @examples
#' 
#' # Restrict vebosity to info messages only.
#' set_verbosity(1)
#' 
#' # Load a dataset
#' load_example_dataset("7871581/files/pbmc3k_medium")
#' 
#' # Select informative genes
#' res <- select_genes(pbmc3k_medium,
#'                     distance = "pearson",
#'                     row_sum=5)
#' 
#' # Result is a ClusterSet object
#' is(res)
#' slotNames(res)
#' 
#' # The selected genes
#' nrow(res)
#' head(row_names(res))
#' 
#' @export select_genes

select_genes <- function(data = NULL,
                         distance_method = c("pearson",
                                             "cosine",
                                             "euclidean",
                                             "spearman",
                                             "kendall",
                                             "binary"),
                         noise_level = 0.00005,
                         k = 80,
                         row_sum = 1,
                         fdr = 0.005,
                         which_slot = c("data", "sct", "counts"),
                         no_dknn_filter = FALSE,
                         no_anti_cor=FALSE,
                         final_slot=c("data", "sct", "counts"),
                         seed = 123) {


  # Set parameters ----------------------------------------------------------

  ## set a seed for reproductibility
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  final_slot <- match.arg(final_slot)
  distance_method <- match.arg(distance_method)
  which_slot <- match.arg(which_slot)
  
  print_msg(paste0("Slot for computing distances was set to: ", which_slot), msg_type = "INFO")
  print_msg(paste0("Final expression slot was set to: ", final_slot), msg_type = "INFO")
  
  if(distance_method == "binary" & which_slot != 'counts'){
    print_msg("Method binary was selected forcing the use of 'counts', slot", msg_type = "INFO")
    which_slot <- 'counts'
  }
  
  data <- get_data_for_scigenex(data = data, which_slot=which_slot)
  

  # Check if noise_level is between 0 and 1
  if (noise_level < 0 | noise_level > 1) {
    print_msg("noise_level argument should be >= 0 and <= 1.",
              msg_type = "STOP")              
  }
  

  # Computing the distance matrix -------------------------------------------

  # Remove genes with small counts
  
  if(!inherits(data, "dgCMatrix")){
    genes_to_keep <- rowSums(data) > row_sum
  }else{
    genes_to_keep <- Matrix::rowSums(data) > row_sum
  }
  
  mtx_sel <- data[genes_to_keep, ]
  
  # Compute gene-gene correlation/distance matrix
  
  print_msg(
    paste0(
      "Computing distances using selected method: ",
      distance_method
    ),
    msg_type = "INFO"
  )
  print_msg(paste0("Number of selected rows/genes (row_sum): ", nrow(mtx_sel)), 
            msg_type = "DEBUG")
  
  # Compute correlations and corresponding 
  # distance matrix. Note that for pearson
  # and cosine, 0 < distance < 2.
  if (distance_method == "pearson") {
    dist_matrix <- qlcMatrix::corSparse(t(mtx_sel))
    if(no_anti_cor)
      dist_matrix[dist_matrix < 0] <- 0
    dist_matrix <- 1 - dist_matrix
  }else  if (distance_method == "kendall") {
    dist_matrix <- as.matrix(cor(t(mtx_sel), method = "kendall"))
    if(no_anti_cor)
      dist_matrix[dist_matrix < 0] <- 0
    dist_matrix <- 1 - dist_matrix  
  }else  if (distance_method == "spearman") {
      dist_matrix <- as.matrix(cor(t(mtx_sel), method = "spearman"))
      if(no_anti_cor)
        dist_matrix[dist_matrix < 0] <- 0
      dist_matrix <- 1 - dist_matrix                      
  } else if (distance_method == "cosine") {
    dist_matrix <- as.matrix(qlcMatrix::cosSparse(t(mtx_sel)))
    if(no_anti_cor)
      dist_matrix[dist_matrix < 0] <- 0
    dist_matrix <- 1 - dist_matrix
  } else if (distance_method == "euclidean") {
    dist_matrix <- as.matrix(dist(mtx_sel))
  }else if (distance_method == "binary") {
    dist_matrix <- ((1 - mtx_sel) %*% t(mtx_sel) + mtx_sel %*% t(1 - mtx_sel)) / (mtx_sel %*% t(mtx_sel) + (1 - mtx_sel) %*% t(mtx_sel) + (mtx_sel) %*% t(1 - mtx_sel))
  }
  
  print_stat("Distance matrix stats", 
             data = dist_matrix, msg_type = "DEBUG")
  
  # Compute min and max distances.
  # They will be used later to compute weights.
  min_dist <- min(dist_matrix)
  max_dist <- max(dist_matrix)
  
  # Set the rownames / colnames of the distance matrix
  rownames(dist_matrix) <- rownames(mtx_sel)
  colnames(dist_matrix) <- rownames(mtx_sel)
  
  # The distance from a gene to itself is 'hidden'
  diag(dist_matrix) <- NA
  
  

  # Compute distance to KNN -------------------------------------------------
  
  print_msg(paste0("Computing distances to KNN."), msg_type = "INFO")
  
  # Create a dataframe to store the DKNN values.
  # Gene_id appear both as rownames and column
  # for coding convenience
  df_dknn <- data.frame(
    dknn_values = rep(NA, nrow(dist_matrix)),
    row.names = rownames(dist_matrix),
    gene_id = rownames(dist_matrix)
  )
  
  # This list will be used to store the
  # dknn values
  l_knn <- list()
  
  # Extract the DKNN for each gene
  for (pos in seq_len(nrow(dist_matrix))) {
    gene <- rownames(df_dknn)[pos]
    gene_dist <- dist_matrix[gene, ]
    # Reorder the distance values in increasing order.
    # The distance from a gene to itself (previously set to NA)
    # is placed at the end of the vector.
    gene_dist <- sort(gene_dist, na.last = T)[1:k]
    
    # Add the neigbhors to the list
    l_knn[[gene]] <- gene_dist
    
    # Select the kth distance values. 
    # This value corresponds to the DKNN of the gene(i)
    df_dknn[gene, "dknn_values"] <- gene_dist[k]
  }
  
  print_stat("Observed DKNN stats", 
             data = df_dknn$dknn_values, msg_type = "DEBUG")
  

  # Compute simulated distance to KNN ---------------------------------------

  sim_dknn <- vector()
  critical_distance <- vector()
  
  if(!no_dknn_filter){
    
    #################### DKNN simulation
    print_msg(paste0("Computing simulated distances to KNN."), msg_type = "INFO")
    
    nb_selected_genes <- nrow(mtx_sel)
    
    if (noise_level == 0) {
      tresh_dknn <- min(df_dknn$dknn_values)
    } else if (noise_level == 1) {
      tresh_dknn <- max(df_dknn$dknn_values)
    } else {
      tresh_dknn <- stats::quantile(df_dknn$dknn_values, noise_level)
    }
    
    gene_with_low_dknn <- df_dknn$gene_id[df_dknn$dknn_values > tresh_dknn]
    dist_values_sub <- dist_matrix[gene_with_low_dknn, ]
    dist_values_sub <- dist_values_sub[!is.na(dist_values_sub)]
    
    
    for (sim_nb in 1:nb_selected_genes) {
      # Randomly sample distances for one simulated gene
      dist_sim <- sample(dist_values_sub,
                         size = nb_selected_genes,
                         replace = FALSE)
      
      # Extract the k nearest neighbors of these simulated gene
      dist_sim <- sort(dist_sim)
      sim_dknn[sim_nb] <- dist_sim[k]
    }
    
    print_stat("Simulated DKNN stats", 
               data = sim_dknn, msg_type = "DEBUG")
    
    # The simulated DKNN values follow a normal distribution.
    # Compute the parameters of this distribution
    mean_sim <- mean(sim_dknn)
    sd_sim <- sd(sim_dknn)
    

    # Compute the DKNN threshold (or critical distance) -----------------------

    # Order genes by DKNN values
    print_msg(paste0("Computing threshold of distances to KNN (DKNN threshold)."),
              msg_type = "INFO")
    
    # Order the dknn values from low to high
    df_dknn <- df_dknn[order(df_dknn$dknn_values), ]
    df_dknn[, "FDR"] <- NA
    
    # Compute the FDR
    for (i in 1:nb_selected_genes) {
      gene <- df_dknn$gene_id[i]
      df_dknn[gene, "FDR"] <-
        stats::pnorm(df_dknn[gene, "dknn_values"],
                     mean = mean_sim,
                     sd = sd_sim,
                     lower.tail = T) / (i / nb_selected_genes) * 100
    }
    
    df_dknn$FDR[df_dknn$FDR > 100] <- 100
    
    # Select genes with a distance value under critical distance ------------

    print_msg(paste0("Selecting informative genes."), msg_type = "INFO")
    fdr_tresh_pos <- which(df_dknn$FDR > fdr)[1]
    selected_genes <- df_dknn[1:fdr_tresh_pos,]$gene_id
    critical_distance <-df_dknn$dknn_values[fdr_tresh_pos]
    
  }else{
    selected_genes <- rownames(dist_matrix)
  }
  


  # Create the ClusterSet object --------------------------------------------

  print_msg(paste0("Creating the ClusterSet object."), msg_type = "INFO")

  if(final_slot != which_slot){
    data <- get_data_for_scigenex(data = data, which_slot=final_slot)
  }  
  
  obj <- new("ClusterSet")

  if (length(selected_genes) > 0) {
    obj@data <- as.matrix(data[selected_genes, ])
    obj@gene_clusters <- list("1" = rownames(obj@data))
    obj@gene_clusters_metadata <- list("cluster_id" = as.numeric(names(obj@gene_clusters)),
                                       "number" = max(names(obj@gene_clusters)),
                                       "size" = nrow(obj@data))
    
    obs_dknn <- as.vector(df_dknn[, "dknn_values"])
    names(obs_dknn) <- df_dknn[, "gene_id"]
    if(length(selected_genes) > 1){
      center = matrix(apply(obj@data[obj@gene_clusters$`1`, ],
                            2,
                            mean,
                            na.rm = TRUE),
                      nrow = 1)
    }else{
      center = matrix(obj@data[obj@gene_clusters$`1`, ],
                      nrow = 1)
    }
    
    obj@dbf_output <- list("dknn" = obs_dknn,
                           "simulated_dknn" = sim_dknn,
                           "critical_distance" = critical_distance,
                           "fdr" = df_dknn$FDR,
                           "center" = center,
                           "all_gene_expression_matrix" = data,
                           "all_neighbor_distances" = l_knn[selected_genes])
    # obj@normal_model_mean <- mean_sim
    # obj@normal_model_sd <- mean_sd
  }else{
    
    obs_dknn <- as.vector(df_dknn[, "dknn_values"])
    names(obs_dknn) <- df_dknn[, "gene_id"]
    obj@dbf_output <- list("dknn" = obs_dknn,
                           "simulated_dknn" = sim_dknn,
                           "critical_distance" = critical_distance,
                           "fdr" = df_dknn$FDR,
                           "center" = NULL,
                           "all_gene_expression_matrix" = data,
                           "all_neighbor_distances" = NULL)
  }
  
  obj@parameters <- list("distance_method" = distance_method,
                         "k" = k,
                         "noise_level" = noise_level,
                         "fdr" = fdr,
                         "row_sum" = row_sum,
                         "no_dknn_filter" = no_dknn_filter,
                         "seed" = seed)
  
  return(obj)
}
