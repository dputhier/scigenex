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
#' # Set verbosity to 1 to display info messages only.
#' set_verbosity(1)
#' 
#' # Create a matrix with 4 signatures
#' m <- create_3_rnd_clust()
#' 
#' # Select informative genes
#' res <- select_genes(m,
#'                     distance = "kendall",
#'                     k = 75,
#'                     highest = 0.3,
#'                     fdr = 1e-8,
#'                     row_sum = -Inf)
#'                     
#' # Cluster informative features
#' res <- gene_clustering(res, 
#'                        inflation = 1.2,
#'                        keep_nn = FALSE,
#'                        k = 5)
#' is(res)
#' 
#' # Plot heatmap of gene clusters
#' plot_heatmap(res, row_labels = FALSE, line_size_horizontal = 2)
#' plot_heatmap(res[1,], row_labels = FALSE, line_size_horizontal = 2)
#' plot_heatmap(res[1:2, ], row_labels = FALSE, line_size_horizontal = 2)
#' plot_heatmap(res[1:2, 1:15], row_labels = FALSE, line_size_horizontal = 2)
#' 
#' # plot the profiles
#' idents <- c(rep(1,10), rep(2,10))
#' names(idents) <- colnames(res@data)
#' plot_profiles(res,
#'               ident = idents,
#'              color_cell_type = c("1" = "#4E79A7", "2" = "#A0CBE8"))
#'              
#' # Some methods of the ClusterSet object
#' x <- ncol(res)
#' x <- nrow(res)
#' x <- dim(res)
#' x <- col_names(res)
#' x <- row_names(res)
#' x <- get_genes(res)
#' x<-  clust_size(res)
#' x <- c("gene400", "gene141") %in% res
#' x <- which_clust(res, genes = c("gene400", "gene141"))
#' res <- top_genes(res, top=5)
#' res <- res[2:3, ]
#' clust_names(res)
#' scigenex_test_I1.2 <- rename_clust(res)
#' clust_names(res)
#' res <- res[,col_names(res)[1:10]]
#' show(res)    
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

#' @title
#' The show method of a ClusterSet.
#' @description
#' The show method of a ClusterSet.
#' @param object A ClusterSet object.
#' @export show
#' @keywords internal
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

#' @title
#' ncol.ClusterSet
#' @description
#' The number of column of a ClusterSet object.
#' @param x The ClusterSet object
#' @keywords internal
ncol.ClusterSet <- function (x) {
  ncol(x@data)
}

#' @title
#' nrow.ClusterSet 
#' @description
#' The number of rows of a ClusterSet object.
#' @param x The ClusterSet object
#' @keywords internal
nrow.ClusterSet <- function (x) {
  nrow(x@data)
}

#' @title Names of gene clusters stored in the ClusterSet object
#' @description
#' The names of the gene clusters stored in the ClusterSet object.
#' @param x The ClusterSet object.
#' @export clust_names
#' @keywords internal
setGeneric("clust_names", 
           function(x)
             standardGeneric("clust_names")
)

#' @title Names of gene clusters stored in the ClusterSet object
#' @description
#' The names of the gene clusters stored in the ClusterSet object.
#' @param x The ClusterSet object
#' @export clust_names
setMethod("clust_names", "ClusterSet",
          function(x){
            clust <- unlist(mapply(rep, 
                                   names(x@gene_clusters), 
                                   lapply(x@gene_clusters, length)))
            names(clust) <- unlist(mapply(rep, 
                                          names(x@gene_clusters), 
                                          lapply(x@gene_clusters, length)))

            return(clust)
          }
            
)

#' @title Dimension of a ClusterSet object.
#' dim
#' @description
#' The dimension of a ClusterSet object.
#' @param x The ClusterSet object
#' @keywords internal
setMethod("dim",signature(x="ClusterSet"),
          function(x) {
            dim(x@data)
          }
) 


#' @title Column names of an object
#' @description
#' The column names of a ClusterSet object.
#' @param x The ClusterSet object
#' @export col_names
#' @keywords internal
setGeneric(name="col_names", def=function(x) standardGeneric("col_names"))

#' @title Column names of a ClusterSet object.
#' @description
#' The column names of a ClusterSet object.
#' @param x The ClusterSet object
#' @export col_names
setMethod(f="col_names", 
          signature="ClusterSet", 
          definition=function(x) colnames(x@data))


#' @title Row names of a ClusterSet object.
#' row_names
#' @description 
#' The row names of a ClusterSet object.
#' @param x The ClusterSet object
#' @export row_names
#' @keywords internal
setGeneric("row_names", 
             function(x)
               standardGeneric("row_names"))

#' @title Row names of a ClusterSet object.
#' row_names
#' @description 
#' The row names of a ClusterSet object.
#' @param x The ClusterSet object
#' @export row_names
setMethod("row_names", "ClusterSet",
          function(x)
            rownames(x@data)
)


################################################################################
##      Method for function"[". Subsetting 
##      ClusterSet object
################################################################################

#' @title Subsetting operator of a ClusterSet object
#' Extract
#' @description
#' The subsetting operator of a ClusterSet object.
#' The i axis correspond to clusters and j to column/cells
#' @param i indices specifying rows to extract. Indices are numeric or character vectors or empty (missing) or NULL.
#' @param j indices specifying column to extract. Indices are numeric or character vectors or empty (missing) or NULL.
#' @param ... See ?'['. Not functionnal here.
#' @param drop For matrices and arrays. If TRUE the result is coerced to the lowest possible dimension. Not functionnal here.
#' @keywords internal
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
##      Method nclust/clust_size for a 
##      ClusterSet object
################################################################################
#' @title Number of clusters in a ClusterSet object.
#' nclust
#' @description
#' The number of clusters in a ClusterSet object.
#' @param x The ClusterSet object
#' @export nclust
#' @examples 
#' # load a dataset
#' data("scigenex_test_I1.2")
#' n_clust <- nclust(scigenex_test_I1.2)
#' @keywords internal
setGeneric("nclust", 
           function(x)
             standardGeneric("nclust")
)


#' @title Number of clusters in a ClusterSet object.
#' nclust
#' @description
#' The number of clusters in a ClusterSet object.
#' @param x The ClusterSet object
#' @export nclust
#' @examples 
#' # load a dataset
#' data("scigenex_test_I1.2")
#' n_clust <- nclust(scigenex_test_I1.2)
setMethod(
  "nclust", signature("ClusterSet"),
  function(x) {
    x@gene_clusters_metadata$number
  }
)

#' @title Sizes of the clusters stored in a ClusterSet object
#' clust_size
#' @description
#' The sizes of the clusters stored in a ClusterSet object.
#' @param x A ClusterSet object.
#' @export clust_size
#' @examples
#' # load a dataset
#' data("scigenex_test_I1.2")
#' x <- clust_size(scigenex_test_I1.2)
#' @keywords internal
setGeneric("clust_size", 
           function(x)
             standardGeneric("clust_size")
)

#' @title Sizes of the clusters stored in a ClusterSet object
#' clust_size
#' @description
#' The sizes of the clusters stored in a ClusterSet object.
#' @param x A ClusterSet object.
#' @export
#' @examples
#' # load a dataset
#' data("scigenex_test_I1.2")
#' x <- clust_size(scigenex_test_I1.2)
setMethod(
  "clust_size", signature("ClusterSet"),
  function(x) {
    x@gene_clusters_metadata$size
  }
)

#################################################################
##    Define gene_cluster function for ClusterSet object
#################################################################

#' @title The gene clusters stored in a ClusterSet.
#' gene_cluster
#' @description
#' Returns a named vector (gene as names) and cluster
#' as value.
#' @param object a ClusterSet object.
#' @param cluster The cluster of interest. 0 means all cluster. Otherwise a non-null integer value.
#' @examples
#' # load a dataset
#' data("scigenex_test_I1.2")
#' g_clust <- gene_cluster(scigenex_test_I1.2)
#' @export gene_cluster
#' @keywords internal
setGeneric("gene_cluster", 
           function(object,
                    cluster = 0)
             standardGeneric("gene_cluster")
)

#' @title The gene clusters stored in a ClusterSet.
#' gene_cluster
#' @description
#' Returns a named vector (gene as names) and cluster
#' as value.
#' @param object a ClusterSet object.
#' @param cluster The cluster of interest. 0 means all cluster. Otherwise a non-null integer value.
#' @examples
#' # load a dataset
#' data("scigenex_test_I1.2")
#' g_clust <- gene_cluster(scigenex_test_I1.2)
#' @export gene_cluster
setMethod(
  "gene_cluster", signature("ClusterSet"),
      function(object,
      cluster = 0) {
  
  if(!is.null(object@gene_clusters)){
    nb_clust <- length(object@gene_clusters)
  }else{
    print_msg("There is no cluster in this object.",
              msg_type = 'STOP')
  }
  
  if(!is.numeric(cluster))
    print_msg("Please provide a numeric value.",
              msg_type = 'STOP')
  
  cluster <- unique(cluster)
  
  if (!all(cluster-floor(cluster)==0) | any(cluster < 0 | any(cluster > nb_clust)))
    print_msg("Please provide a zero (all clusters) or positive integer in the required range.",
              msg_type = 'STOP')
  
  if(length(cluster) == 1){
    if (cluster == 0)
      cluster <- 1:length(object@gene_clusters)
  }
  
  if(length(cluster) > 1){
    if(length(cluster[cluster == 0]))
      print_msg("Zero is out of range.",
                msg_type = 'STOP')
  }
  
  
  if (nb_clust) {
    cluster_as_int <- unlist(mapply(rep,
                                    cluster,
                                    lapply(object@gene_clusters[cluster], length),
                                    SIMPLIFY = TRUE))
    cluster_as_int <- as.vector(as.matrix(cluster_as_int))
    names(cluster_as_int) <-
      unlist(object@gene_clusters[cluster])
    return(cluster_as_int)
    
  } else{
    return(NULL)
  }
})


################################################################################
##      Method for function matching genes in a ClusterSet.
################################################################################

#' @title Match operator of a ClusterSet object
#' @description The match operator of a ClusterSet object
#' @param x The gene to be searched;
#' @param table The ClusterSet object.
#' @keywords internal
setMethod("%in%", signature(x = "character", table = "ClusterSet"), function(x, table) {
  x %in% names(gene_cluster(table))
})


#' @title Which clusters contain a set of genes.
#' @description Which clusters contain a set of genes.
#' @param object a ClusterSet object.
#' @param genes The genes to be searched in the ClusterSet.
#' @export which_clust
#' @examples
#' # load a dataset
#' data("scigenex_test_I1.2")
#' hit <- which_clust(scigenex_test_I1.2, genes = c("gene400", "gene141"))
#' @keywords internal
setGeneric("which_clust", 
           function(object,
                    genes = NULL)
             standardGeneric("which_clust")
)

#' @title Which clusters contain a set of genes.
#' @description Which clusters contain a set of genes.
#' @param object a ClusterSet object.
#' @param genes The genes to be searched in the ClusterSet.
#' @export which_clust
#' @examples
#' # load a dataset
#' data("scigenex_test_I1.2")
#' hit <- which_clust(scigenex_test_I1.2, genes = c("gene400", "gene141"))
setMethod("which_clust", 
          signature("ClusterSet"), 
          function(object, genes) {
            check_format_cluster_set(object)
            gc <- gene_cluster(object)
            tmp <- gc[which(names(gc) %in% genes)]
            tmp[genes]
})

#' @title Rename the gene clusters of a ClusterSet
#' @description Rename the gene clusters of a ClusterSet.
#' @param object a ClusterSet object.
#' @param new_names The new names for the clusters.
#' @export rename_clust
#' @examples
#' # load a dataset
#' data("scigenex_test_I1.2")
#' new_obj <- rename_clust(scigenex_test_I1.2, new_names = letters[1:nclust(scigenex_test_I1.2)])
#' @keywords internal
setGeneric("rename_clust", 
           function(object, new_names=NULL)
             standardGeneric("rename_clust")
)

#' @title Rename the gene clusters of a ClusterSet
#' @description Rename the gene clusters of a ClusterSet.
#' @param object a ClusterSet object.
#' @param new_names The new names for the clusters.
#' @export rename_clust
#' @examples
#' # load a dataset
#' data("scigenex_test_I1.2")
#' new_obj <- rename_clust(scigenex_test_I1.2, new_names = letters[1:nclust(scigenex_test_I1.2)])
setMethod("rename_clust", 
          signature("ClusterSet"), 
          function(object, 
                   new_names=NULL) {
            
            check_format_cluster_set(object)
            
            if(is.null(new_names))
              new_names <- 1:nclust(object)
            
            if(length(new_names) != nclust(object))
              print_msg("The number of labels should be the same a the number of clusters.")
            
            names(object@gene_clusters) <- new_names
            
            if(length(object@top_genes) > 0)
              names(object@top_genes) <- new_names
            
            object@gene_clusters_metadata$cluster_id <- new_names
            
            names(object@gene_clusters_metadata$size) <- new_names
            
            if(length(object@gene_cluster_annotations) > 0)
              names(object@gene_cluster_annotations) <- new_names
            
            rownames(object@dbf_output$center) <- new_names
            
            return(object)
            
})

