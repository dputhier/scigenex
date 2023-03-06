#################################################################
##    DEFINITION OF A SPECIFIC CLASS OBJECT : ClusterSet
#################################################################

#' @title
#' ClusterSet-class
#' @description
#' This class is a representation of a partitioning algorithm and is intented to store gene clusters.
#' @slot data A matrix containing the filtered and partitioned data.
#' @slot gene_clusters A list contains the partitioned genes of the dataset. Each element of the list corresponds to a cluster, and contains the indices of the genes assigned to that cluster.
#' @slot top_genes A list contains the top genes from the gene clusters. Each element of the list corresponds to a cluster, and contains the indices of the genes assigned to that cluster ranked by their correlation value within their cluster.
#' @slot gene_clusters_metadata A list contains metadata related to the gene clusters such as the number of gene clusters, their ID, and the number of genes contained in each of them.
#' @slot gene_cluster_annotations A list contains the result obtained from the GO enrichment analysis of gene clusters.
#' @slot cells_metadata A list containing metadata related to the cell clusters such as the clustering results the number of cell clusters, their order, colors associated to each cluster,... 
#' @slot dbf_output A list containing the intermediates outputs of the DBF function : dknn, simulated distances, critical distance and fdr values. 
#' @slot parameters A list containing the parameter used. Each element of the list correspond to a parameter.
#' 
#' @return A ClusterSet object.
#' @export
#'
#' @examples
#' 
#' \dontrun{
#' ## An artificial example
#' ## with continuous values
#' m <- create_3_rnd_clust()
#' res <- find_gene_clusters(data=m,
#'                              distance_method="pearson",
#'                              inflation = 2,
#'                              k=75,
#'                              row_sum=-Inf,
#'                              highest=0.3,
#'                              min_nb_supporting_cell = 0,
#'                              fdr = 1e-8)
#' is(res)
#' res
#' plot_heatmap(res, row_labels = F, line_size_horizontal = 2)
#' 
#' plot_profiles(res)
#' ## The subset operator may be used to select clusters/rows
#' ## and cells/columns of interest
#' res[1,]
#' res[1:2,]
#' }          
setClass("ClusterSet",
         representation = list(
           data = "matrix",
           gene_clusters = "list",
           top_genes = "list",
           gene_clusters_metadata = "list",
           gene_cluster_annotations = "list",
           cells_metadata = "data.frame",
           dbf_output = "list",
           parameters = "list"
         ),
         prototype = list(
           data = matrix(nr = 0, nc = 0),
           gene_clusters = list(),
           top_genes = list(),
           gene_clusters_metadata = list(),
           gene_cluster_annotations = list(),
           cells_metadata = data.frame(),
           dbf_output = list(),
           parameters = list()
         )
)

#################################################################
##    REDEFINE SHOW METHOD FOR CLASS OBJECT : ClusterSet
#################################################################

setMethod(
  "show", signature("ClusterSet"),
  
  function(object) {
    
    cat("\t\tAn object of class ClusterSet\n")
    cat("\t\tName:", slot(object, "parameters")$name, "\n")
    cat("\t\tMemory used: ", object.size(object), "\n")
    cat("\t\tNumber of cells: ", ncol(slot(object, "data")), "\n")
    cat(
      "\t\tNumber of informative genes: ",
      nrow(slot(object, "data")), "\n"
    )
    cat("\t\tNumber of gene clusters: ", slot(object, "gene_clusters_metadata")$number, "\n")
    cat("\t\tThis object contains the following informations:\n")
    
    for (i in slotNames(object)) {
      cat("\t\t\t - ", i, "\n")
    }
    if (length(slot(object, "parameters")) > 0) {
      for (i in 1:length(slot(object, "parameters"))) {
        cat(
          "\t\t\t\t * ", names(slot(object, "parameters"))[[i]],
          " = ", slot(object, "parameters")[[i]], "\n"
        )
      }
    }
  }
)


################################################################################
##      NCOL/NROW/DIM METHOD FOR CLASS OBJECT : ClusterSet
################################################################################

ncol.ClusterSet <- function (x) {
  ncol(x@data)
}


nrow.ClusterSet <- function (x) {
  nrow(x@data)
}


if (!isGeneric("dim"))
  setGeneric("dim", 
             function(x)
               standardGeneric("dim")
  )


setMethod("dim",signature(x="ClusterSet"),
          function(x) {
            dim(x@data)
          }
) 


colnames.ClusterSet <- function (x, do.NULL = TRUE, prefix = "col") {
            colnames(x@data)
}

rownames.ClusterSet <- function (x, do.NULL = TRUE, prefix = "row") {
  names(x@data)
}



################################################################################
##      Method for function"[". Subsetting 
##      ClusterSet object
################################################################################

setMethod("[", signature(x = "ClusterSet"),
          function (x, i, j, ..., drop = FALSE) {
            
            if(is.null(names(x@gene_clusters_metadata$cluster_id)))
               names(x@gene_clusters_metadata$cluster_id) <-  names(x@gene_clusters)
            
            if(is.null(names(x@gene_clusters_metadata$size)))
               names(x@gene_clusters_metadata$size) <-  names(x@gene_clusters)
               
            if (missing(j)) {
              if (missing(i)) {
                n_data <- x@data
                n_gene_clusters <- x@gene_clusters
                n_top_genes <- x@top_genes
                n_gene_clusters_metadata <- x@gene_clusters_metadata
                n_gene_cluster_annotations <- x@gene_cluster_annotations
                n_cells_metadata <- x@cells_metadata
                n_dbf_output <- x@dbf_output
              }else {
                n_data <- x@data[unlist(x@gene_clusters[i]), , drop=FALSE]
                n_gene_clusters <- x@gene_clusters[i]
                
                if(length(x@top_genes)){
                  n_top_genes <- x@top_genes[i]
                }else{
                  n_top_genes <- x@top_genes
                }
                
                n_gene_clusters_metadata <- x@gene_clusters_metadata
                n_gene_clusters_metadata$cluster_id <- 
                  x@gene_clusters_metadata$cluster_id[i]
                n_gene_clusters_metadata$number <- length(n_gene_clusters)
                n_gene_clusters_metadata$size <- x@gene_clusters_metadata$size[i]
                
                if(length(x@gene_cluster_annotations) > 0){
                  n_gene_cluster_annotations <- x@gene_cluster_annotations[i]
                }else{
                  n_gene_cluster_annotations <- x@gene_cluster_annotations
                }
                  
                n_cells_metadata <- x@cells_metadata
                n_dbf_output <- x@dbf_output
                n_dbf_output$center <- n_dbf_output$center[i, , drop=FALSE]
              }
            } else {
              if (missing(i)) {
                n_data <- x@data[,j]
                n_gene_clusters <- x@gene_clusters
                n_top_genes <- x@top_genes
                n_gene_clusters_metadata <- x@gene_clusters_metadata
                n_gene_cluster_annotations <- x@gene_cluster_annotations
                n_cells_metadata <- x@cells_metadata[j, , drop=FALSE]
                n_dbf_output <- x@dbf_output
                n_dbf_output$center <- n_dbf_output$center[ , j, drop=FALSE]
              }else {
                n_data <- x@data[unlist(x@gene_clusters[i]), j, drop=FALSE]
                n_gene_clusters <- x@gene_clusters[i]
                
                if(length(x@top_genes)){
                  n_top_genes <- x@top_genes[i]
                }else{
                  n_top_genes <- x@top_genes
                }
                
                n_gene_clusters_metadata <- x@gene_clusters_metadata
                n_gene_clusters_metadata$cluster_id <- x@gene_clusters_metadata$cluster_id[i]
                n_gene_clusters_metadata$number <- length(i)
                n_gene_clusters_metadata$size <- x@gene_clusters_metadata$size[i]
                if(length(x@gene_cluster_annotations) > 0){
                  n_gene_cluster_annotations <- x@gene_cluster_annotations[i]
                }else{
                  n_gene_cluster_annotations <- x@gene_cluster_annotations
                }
                n_cells_metadata <- x@cells_metadata[j, , drop=FALSE]
                n_dbf_output <- x@dbf_output
                n_dbf_output$center <- n_dbf_output$center[i, j, drop=FALSE]
              }
            }
            
            
            new(
              "ClusterSet",
              data = n_data,
              gene_clusters = n_gene_clusters,
              top_genes = n_top_genes,
              gene_clusters_metadata = n_gene_clusters_metadata,
              gene_cluster_annotations = n_gene_cluster_annotations,
              cells_metadata = n_cells_metadata,
              dbf_output = n_dbf_output,
              parameters = x@parameters
            )
          })


################################################################################
##      Method nclust for a 
##      ClusterSet object
################################################################################


setGeneric("nclust", 
           function(x)
             standardGeneric("nclust")
)

#' @export
setMethod(
  "nclust", signature("ClusterSet"),
  function(x) {
    x@gene_clusters_metadata$number
  }
)
