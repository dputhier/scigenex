#################################################################
##    Define top_genes function for ClusterSet object
#################################################################
#' @title Most co-expressed genes from each gene cluster
#' @description
#' Extract the most highly co-expressed genes of each gene cluster.
#' @param object A \code{ClusterSet} object.
#' @param top A value for the number of genes to select from each cluster.
#'
#' @return A \code{ClusterSet} object.
#' @export top_genes
#' @keywords internal
setGeneric("top_genes", 
           function(object,
                    top = 20)
             standardGeneric("top_genes")
)

#' @title Most co-expressed genes from each gene cluster
#' @description
#' Extract the most highly co-expressed genes of each gene cluster.
#' @param object A \code{ClusterSet} object.
#' @param top A value for the number of genes to select from each cluster.
#'
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
                   top = 20) {
            
            ## Check format object arg
            check_format_cluster_set(object)
            
            # Display an info message if there is less than 
            # n top genes in a gene cluster
            if(any(object@gene_clusters_metadata$size < top)){
              print_msg(paste0("One our several cluster contains less than ", 
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
              dist_method <- object@parameters$distance_method
              
              if (dist_method %in% c("kendall", "pearson", "spearman")) {
                dist <- cor(t(object@data[genes, ]), method = dist_method)
                dist <- 1 - dist
              }
              
              if (dist_method == "cosine") {
                dist <- as.matrix(qlcMatrix::cosSparse(t(object@data[genes, ])))
                dist <- 1 - dist
              }
              
              if (dist_method == "euclidean") {
                dist <- as.matrix(dist(object@data[genes, ],
                                       method = "euclidean",
                                       upper = TRUE,
                                       diag = TRUE
                ))
                diag(dist) <- NA
              }
              
              
              # Compute mean correlation for each gene in cluster i
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

