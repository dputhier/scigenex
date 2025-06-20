#' @title Compute Gene Module Quality Scores
#' @description This function computes various quality metrics (connectivity, dunn index and silhouette coefficient) for gene modules in a ClusterSet object.
#' @param object a ClusterSet object.
#' @return A data frame containing the following cluster quality scores:
#' \itemize{
#'   \item \code{connectivity}: Connectivity score.
#'   \item \code{dunn_index}: Dunn index.
#'   \item \code{silhouette_coefficient}: Mean silhouette coefficient.
#' }
#' @examples
#' # Load your data object
#' load_example_dataset('7871581/files/pbmc3k_medium_clusters')
#' # Compute quality scores of gene modules
#' quality_scores <- module_quality_scores(pbmc3k_medium_clusters)
#' quality_scores
#' @export
#' @noRd
setGeneric("module_quality_scores", 
           function(object)
             standardGeneric("module_quality_scores")
)


#' @title Compute Gene Module Quality Scores
#' @description This function computes various quality metrics (connectivity, dunn index and silhouette coefficient) for gene modules in a ClusterSet object.
#' @param object a ClusterSet object.
#' @return A data frame containing the following cluster quality scores:
#' \itemize{
#'   \item \code{connectivity}: Connectivity score.
#'   \item \code{dunn_index}: Dunn index.
#'   \item \code{silhouette_coefficient}: Mean silhouette coefficient.
#' }
#' @examples
#' # Load your data object
#' load_example_dataset('7871581/files/pbmc3k_medium_clusters')
#' # Compute quality scores of gene modules
#' quality_scores <- module_quality_scores(pbmc3k_medium_clusters)
#' quality_scores
#' @export
setMethod(
  "module_quality_scores", signature("ClusterSet"),
  function(object){
    #Compute correlation matrix
    cor <- corSparse(t(object@data))
    
    # Transform correlation matrix to distance matrix
    diag(cor) <- NA
    dist <- (2 - (cor + 1)) / 2
    rownames(dist) <- rownames(object@data)
    colnames(dist) <- rownames(object@data)
    dist <- as.dist(dist)
    
    # Extract gene clusters
    gene_cl <- gene_cluster(object)
    
    ### Connectivity
    connectivity_score <- clValid::connectivity(dist, gene_cl)
    
    ### Silhouette score
    # Compute silhouette similarity score
    sil <- cluster::silhouette(gene_cl, dist = dist)
    # Extract mean silhouette similarity score
    sil_mean <- mean(sil[,"sil_width"])
    
    ### Dunn index
    du <- clValid::dunn(clusters = gene_cl, distance = dist)
    
    
    ### Store scores in a dataframe
    scores <- data.frame("connectivity" = connectivity_score,
                         "dunn_index" = du,
                         "silhouette_coefficient" = sil_mean)
    
    return(scores)
  }
)
