#################################################################
##    Define top_genes function for ClusterSet object
#################################################################
#' @title Most co-expressed genes from each gene cluster
#' @description
#' Extract the most highly co-expressed genes of each gene cluster.
#' @param object A \code{ClusterSet} object.
#' @param top A value for the number of genes to select from each cluster.
#' @param fast Use qlcMatrix::corSparse for pearson computation. Default to cor(). Default to FALSE for retro-compatibility.
#' @param distance_method Overright the object distance_method (slot parameters$distance_method). Useful when object was created using cluster_set_from_matrix() for instance. One of c("kendall", "spearman", "cosine", "euclidean", "pearson") or NULL.
#' @return A \code{ClusterSet} object.
#' # Set verbosity to 1 to display info messages only.
#' set_verbosity(1)
#' 
#' # Load a dataset
#' load_example_dataset('7871581/files/pbmc3k_medium_clusters')
#' 
#' # Store top genes in the object
#' pbmc3k_medium_clusters <- top_genes(pbmc3k_medium_clusters)
#' @export
#' @noRd
setGeneric("top_genes", 
           function(object,
                    top = 20,
                    fast=FALSE,
                    distance_method=NULL)
             standardGeneric("top_genes")
)

#' @title Most co-expressed genes from each gene cluster
#' @description
#' Extract the most highly co-expressed genes of each gene cluster.
#' @param object A \code{ClusterSet} object.
#' @param top A value for the number of genes to select from each cluster.
#' @param fast This argument is no more functional. Fast computation is used by default now.
#' @param distance_method Overright the object distance_method (slot parameters$distance_method). Useful when object was created using cluster_set_from_matrix() for instance. 
#' @return A \code{ClusterSet} object.
#' @export top_genes
#' @importFrom Matrix t
#' @importFrom qlcMatrix corSparse
#' @importFrom qlcMatrix cosSparse
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
                      fast=TRUE,
                      distance_method=c(NULL, "pearson", "kendall", "spearman", "cosine", "euclidean")) {
    
    distance_method <- match.arg(distance_method)
  
    ## Check format object arg
    check_format_cluster_set(object)
    
    # Display an info message if there is less than 
    # n top genes in a gene cluster
    if(any(object@gene_clusters_metadata$size < top)){
      print_msg(paste0("One or several clusters contain less than ", 
                       top, " genes. Retrieving all genes"),
                msg_type = "INFO")
    }              
    
    # Prepare some vars
    clusters <- object@gene_clusters
    genes_top <- list()
    
    # Extract top co-expressed genes for each gene cluster
    for (cl_name in names(clusters)) {
      # Extract gene names in cluster i
      genes <- clusters[[cl_name]]
      
      # Compute distances between genes in cluster i
      # Use the same distance used by SciGeneX
      if(is.null(distance_method)){
        dist_method <- object@parameters$distance_method
      }else{
        dist_method <- distance_method
      }
        
      
      if (dist_method == "unknown") {
        print_msg("Distance type is unknown", msg_type = "INFO")
        print_msg("Using Pearson", msg_type = "INFO")
        dist_method<- "pearson"
      }
      
      if (dist_method == "pearson"){
        print_msg("Using fast computation of pearson correlations.", msg_type="DEBUG")
        dist <- as.matrix(qlcMatrix::corSparse(Matrix::t(object@data[genes, ])))
        colnames(dist) <- genes
        rownames(dist) <- genes
        dist <- 1 - dist
      }
      
      if(dist_method %in% c("kendall", "spearman")) {
        dist <- cor(t(as.matrix(object@data[genes, ])), method = dist_method)
        dist <- 1 - dist
      }
      
      if (dist_method == "cosine") {
        dist <- as.matrix(qlcMatrix::cosSparse(Matrix::t(object@data[genes, ])))
        dist <- 1 - dist
        colnames(dist) <- genes
        rownames(dist) <- genes
      }
      
      if (dist_method == "euclidean") {
        dist <- as.matrix(stats::dist(as.matrix(object@data[genes, ]),
                               method = "euclidean",
                               upper = TRUE,
                               diag = TRUE
        ))
        diag(dist) <- NA
        
      }
      
      
      dist_means <- colMeans(dist, na.rm = TRUE)
      dist_means <- dist_means[order(dist_means, decreasing = FALSE)]  


      
      # Extract top genes with the highest correlation mean
      max_g <- ifelse(top > length(dist_means), length(dist_means), top) 
      genes_top[[cl_name]] <- names(dist_means[1:max_g])
    }
    
    
    object@top_genes <- genes_top
    
    # Print
    print_msg(
      msg = "Results are stored in 'top_genes' slot of the object.",
      msg_type = "INFO"
    )
    
    return(object)
  })

