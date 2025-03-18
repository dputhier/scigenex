#################################################################
##    Define top_genes function for ClusterSet object
#################################################################
#' @title Most co-expressed genes from each gene cluster
#' @description
#' Extract the most highly co-expressed genes of each gene cluster.
#' @param object A \code{ClusterSet} object.
#' @param top A value for the number of genes to select from each cluster.
#' @param cluster A vector of gene cluster identity.
#' @param fast Use qlcMatrix::corSparse for pearson computation. Default to cor(). Default to FALSE for retro-compatibility.
#' @param distance_method Overright the object distance_method (slot parameters$distance_method). Useful when object was created using cluster_set_from_matrix() for instance. One of c("kendall", "spearman", "cosine", "euclidean", "pearson") or NULL.
#' @return A \code{ClusterSet} object.
#' @export top_genes
#' @keywords internal
setGeneric("top_genes", 
           function(object,
                    top = 20,
                    cluster = "all",
                    fast=FALSE,
                    distance_method=NULL)
             standardGeneric("top_genes")
)

#' @title Most co-expressed genes from each gene cluster
#' @description
#' Extract the most highly co-expressed genes of each gene cluster.
#' @param object A \code{ClusterSet} object.
#' @param top A value for the number of genes to select from each cluster.
#' @param cluster A vector of gene cluster identity.
#' @param fast Use qlcMatrix::corSparse for pearson computation. Default to cor(). Default to FALSE for retro-compatibility.
#' @param distance_method Overright the object distance_method (slot parameters$distance_method). Useful when object was created using cluster_set_from_matrix() for instance. One of c("kendall", "spearman", "cosine", "euclidean", "pearson") or NULL.   
#' @return A \code{ClusterSet} object.
#' @export top_genes
#'
#' @examples
#' # Set verbosity to 1 to display info messages only.
#' set_verbosity(1)
#' 
#' # Load a dataset
#' load_example_dataset('7871581/files/pbmc3k_medium_clusters')
#' 
#' # Store top genes in the object
#' pbmc3k_medium_clusters <- top_genes(pbmc3k_medium_clusters)
#'
setMethod("top_genes", 
          signature("ClusterSet"), 
              function(object,
                      top = 20,
                      cluster = "all",
                      fast=FALSE,
                      distance_method=NULL) {
  
  ## Check format object arg
  check_format_cluster_set(object)
  
  if (unique(cluster == "all")) {
    cluster <- object@gene_clusters_metadata$cluster_id
  }
  
  # Display an info message if there is less than n top genes in a gene cluster
  loop <- 0
  for (size in object@gene_clusters_metadata$size[cluster]) {
    loop <- loop + 1

    if (top > size) {
      print_msg(paste0("Number of top genes is greater than the number of genes in cluster ", 
                       loop, 
                       ". All genes will be used."),
                msg_type = "INFO")
    }
  }
  
  
  clusters <- object@gene_clusters
  genes_top <- matrix(ncol = top)
  
  print_msg("Extracting top co-expressed genes for each gene cluster", 
            msg_type="DEBUG")
  
  # Compute distances between genes in cluster i
  # Use the same distance used by SciGeneX
  if(!is.null(distance_method)){
    if(!distance_method %in% c("kendall", "spearman", "cosine", "euclidean", "pearson")){
      print('dist_method should be one of "kendall", "spearman", "cosine", "euclidean", "pearson".', 
            msg_type="STOP")
    }

    dist_method <-distance_method

  }else{
    dist_method <- object@parameters$distance_method
  }
  
  for (i in cluster) {
    
    print_msg(paste0("Processing cluster ", i), msg_type="DEBUG")
    
    # Extract gene names in cluster i
    genes <- clusters[[i]]
    
    if (dist_method == "pearson"){
      if(!fast){
        dist <- cor(t(object@data[genes, ]), method = dist_method)
      }else{
        print_msg("Using fast computation of pearson correlations.", msg_type="DEBUG")
        dist <- as.matrix(qlcMatrix::corSparse(t(object@data[genes, ])))
        colnames(dist) <- genes
        rownames(dist) <- genes
      }
      dist <- 1 - dist
    }
    
    if (dist_method == "binary_assym"){
      
      mtx_sel <- object@data[genes, ,drop=FALSE]
      dist <- ((1 - mtx_sel) %*% t(mtx_sel) + mtx_sel %*% t(1 - mtx_sel)) / (mtx_sel %*% t(mtx_sel) + (1 - mtx_sel) %*% t(mtx_sel) + (mtx_sel) %*% t(1 - mtx_sel))
      colnames(dist) <- genes
      rownames(dist) <- genes
      }
    
    if(dist_method %in% c("kendall", "spearman")) {
      dist <- cor(t(object@data[genes, , drop=FALSE]), method = dist_method)
      dist <- 1 - dist
    }
    
    if (dist_method == "cosine") {
      dist <- as.matrix(qlcMatrix::cosSparse(t(object@data[genes, , drop=FALSE])))
      dist <- 1 - dist
      colnames(dist) <- genes
      rownames(dist) <- genes
    }
    
    if (dist_method == "euclidean") {
      dist <- as.matrix(dist(object@data[genes, , drop=FALSE],
                             method = "euclidean",
                             upper = TRUE,
                             diag = TRUE
      ))
      diag(dist) <- NA
    }
    
    if (dist_method == "unknown") {
      print_msg("Distance type is unknown", msg_type = "STOP")
    }

    # Compute mean correlation for each gene in cluster i
    dist_means <- colMeans(dist, na.rm = TRUE)
    dist_means <- dist_means[order(dist_means, decreasing = FALSE)]
    
    # Extract top genes with the highest correlation mean
    genes_top <- rbind(genes_top, names(dist_means[1:top]))
  }
  
  # Prepare top gene matrix
  if (length(cluster) > 1) {
    genes_top <- genes_top[2:(length(cluster) + 1), ]
  } else {
    genes_top <- as.matrix(t(genes_top[2:(length(cluster) + 1), , drop=FALSE]))
  }
  
  colnames(genes_top) <- paste0("gene_top_", 1:top)
  rownames(genes_top) <- paste0("cluster_", cluster)
  
  # Put top gene matrix in object@top_genes
  object@top_genes <- split(x = unname(genes_top), f = cluster)
  
  for (clust in cluster) {
    object@top_genes[[clust]] <- object@top_genes[[as.character(clust)]][!is.na(object@top_genes[[as.character(clust)]])]
  }
  
  # Print
  print_msg(
    msg = "Results are stored in top_genes slot of ClusterSet object.",
    msg_type = "INFO"
  )
  
  return(object)
})
