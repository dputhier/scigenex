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
#' @importFrom SeuratObject SparseEmptyMatrix
#' @export
#'
#' @examples
#' library(Seurat)
#' load_example_dataset("7871581/files/pbmc3k_medium")
#'
#' # Select informative genes
#' res <- select_genes(pbmc3k_medium)
#'
#' # Cluster informative features
#' res <- gene_clustering(res, inflation=1.6)
#' is(res)
#'
#' # Plot heatmap of gene clusters
#' plot_heatmap(res, row_labels = FALSE, line_size_horizontal = 2)
#' plot_heatmap(res[1,], row_labels = FALSE, line_size_horizontal = 2)
#' plot_heatmap(res[1:2, ], row_labels = FALSE, line_size_horizontal = 2)
#' plot_heatmap(res[1:2, 1:15], row_labels = FALSE, line_size_horizontal = 2)
#'
#' # plot the profiles
#' idents <- Seurat::Idents(pbmc3k_medium)
#' plot_profiles(res,
#'               ident = idents)
#'
#' # Some methods of the ClusterSet object
#' x <- ncol(res)
#' x <- nrow(res)
#' x <- dim(res)
#' x <- col_names(res)
#' x <- row_names(res)
#' x <- get_genes(res)
#' x<-  clust_size(res)
#' x <- c("IL32", "CCL5") %in% res
#' x <- which_clust(res, genes = c("IL32", "CCL5"))
#' res <- top_genes(res, top=5)
#' res <- res[2:3, ]
#' res <- rename_clust(res)
#' clust_names(res)
#' res <- res[, col_names(res)[1:10]]
#' show(res)
#' show_methods(res)
setClass(
  "ClusterSet",
  representation = list(
    data = "dgCMatrix",
    gene_clusters = "list",
    top_genes = "list",
    gene_clusters_metadata = "list",
    gene_cluster_annotations = "list",
    cells_metadata = "data.frame",
    dbf_output = "list",
    parameters = "list"
  ),
  prototype = list(
    data = SeuratObject::SparseEmptyMatrix(ncol=0, nrow=0),
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
#' @export
#' @noRd
setMethod("show", signature("ClusterSet"),
          function(object) {
            cat("\t\tAn object of class ClusterSet\n")
            cat("\t\tName:", slot(object, "parameters")$name, "\n")
            cat("\t\tMemory used: ", object.size(object), "\n")
            cat("\t\tNumber of cells: ", ncol(slot(object, "data")), "\n")
            cat("\t\tNumber of informative genes: ",
                nrow(slot(object, "data")), "\n")
            cat(
              "\t\tNumber of gene clusters: ",
              slot(object, "gene_clusters_metadata")$number,
              "\n"
            )
            cat("\t\tThis object contains the following informations:\n")
            
            for (i in slotNames(object)) {
              cat("\t\t\t - ", i, "\n")
            }
            if (length(slot(object, "parameters")) > 0) {
              for (i in 1:length(slot(object, "parameters"))) {
                cat("\t\t\t\t * ",
                    names(slot(object, "parameters"))[[i]],
                    " = ",
                    slot(object, "parameters")[[i]],
                    "\n")
              }
            }
          })

################################################################################
##      NCOL/NROW/DIM METHOD FOR CLASS OBJECT : ClusterSet
################################################################################

#' @title
#' ncol.ClusterSet
#' @description
#' The number of column of a ClusterSet object.
#' @param x The ClusterSet object
#' @noRd
ncol.ClusterSet <- function (x) {
  ncol(x@data)
}

#' @title
#' nrow.ClusterSet
#' @description
#' The number of rows of a ClusterSet object.
#' @param x The ClusterSet object
#' @noRd
nrow.ClusterSet <- function (x) {
  nrow(x@data)
}

#' @title Names of gene clusters stored in the ClusterSet object
#' @description
#' The names of the gene clusters stored in the ClusterSet object.
#' @param x The ClusterSet object.
#' @export
#' @noRd
setGeneric("clust_names",
           function(x)
             standardGeneric("clust_names"))

#' @title The names of the gene clusters stored in the ClusterSet object
#' @description
#' The names of the gene clusters stored in the ClusterSet object.
#' @param x The ClusterSet object
#' @export
#' @examples
#' load_example_dataset('7871581/files/pbmc3k_medium_clusters')
#' clust_names(pbmc3k_medium_clusters)
#' 
setMethod("clust_names",
          "ClusterSet",
          function(x) {
            return(names(x@gene_clusters))
          })

#' @title Dimension of a ClusterSet object.
#' dim
#' @description
#' The dimension of a ClusterSet object.
#' @param x The ClusterSet object
#' @noRd
setMethod("dim", signature(x = "ClusterSet"),
          function(x) {
            dim(x@data)
          })


#' @title Column names of an object
#' @description
#' The column names of a ClusterSet object.
#' @param x The ClusterSet object
#' @export
#' @noRd
setGeneric(
  name = "col_names",
  def = function(x)
    standardGeneric("col_names")
)

#' @title Column names of a ClusterSet object.
#' @description
#' The column names of a ClusterSet object.
#' @param x The ClusterSet object
#' @examples
#' load_example_dataset('7871581/files/pbmc3k_medium_clusters')
#' col_names(pbmc3k_medium_clusters)
#' @export
setMethod(
  f = "col_names",
  signature = "ClusterSet",
  definition = function(x)
    colnames(x@data)
)


#' @title Row names of a ClusterSet object.
#' row_names
#' @description
#' The row names of a ClusterSet object.
#' @param x The ClusterSet object
#' @export
#' @noRd
setGeneric("row_names",
           function(x)
             standardGeneric("row_names"))

#' @title Row names of a ClusterSet object.
#' row_names
#' @description
#' The row names of a ClusterSet object.
#' @param x The ClusterSet object
#' @examples
#' load_example_dataset('7871581/files/pbmc3k_medium_clusters')
#' row_names(pbmc3k_medium_clusters)
#' @export
setMethod("row_names", "ClusterSet",
          function(x)
            rownames(x@data))


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
#' @noRd
#' @importFrom Matrix Matrix
setMethod("[", signature(x = "ClusterSet"),
          function (x, i, j, ..., drop = FALSE) {
            if (is.null(names(x@gene_clusters_metadata$cluster_id)))
              names(x@gene_clusters_metadata$cluster_id) <-
                names(x@gene_clusters)
            
            
            if(is.null(names(x@gene_clusters_metadata$cluster_id)))
               names(x@gene_clusters_metadata$cluster_id) <-  names(x@gene_clusters)
            
            
            if (missing(j)) {
              if (missing(i)) {
                n_data <- x@data
                n_gene_clusters <- x@gene_clusters
                n_top_genes <- x@top_genes
                n_gene_clusters_metadata <- x@gene_clusters_metadata
                n_gene_cluster_annotations <-
                  x@gene_cluster_annotations
                n_cells_metadata <- x@cells_metadata
                n_dbf_output <- x@dbf_output
              } else {
                n_data <- x@data[unlist(x@gene_clusters[i]), , drop = FALSE]
                n_gene_clusters <- x@gene_clusters[i]
                
                if (length(x@top_genes)) {
                  n_top_genes <- x@top_genes[i]
                } else{
                  n_top_genes <- x@top_genes
                }
                
                n_gene_clusters_metadata <- x@gene_clusters_metadata
                n_gene_clusters_metadata$cluster_id <-
                  x@gene_clusters_metadata$cluster_id[i]
                n_gene_clusters_metadata$number <-
                  length(n_gene_clusters)
                n_gene_clusters_metadata$size <-
                  x@gene_clusters_metadata$size[i]
                
                if (length(x@gene_cluster_annotations) > 0) {
                  n_gene_cluster_annotations <- x@gene_cluster_annotations[i]
                } else{
                  n_gene_cluster_annotations <- x@gene_cluster_annotations
                }
                
                n_cells_metadata <- x@cells_metadata
                n_dbf_output <- x@dbf_output
                
                if(!is.null(n_dbf_output$center))
                  n_dbf_output$center <-
                    n_dbf_output$center[i, , drop = FALSE]
              }
            } else {
              if (missing(i)) {
                n_data <- x@data[, j, drop = FALSE]
                n_gene_clusters <- x@gene_clusters
                n_top_genes <- x@top_genes
                n_gene_clusters_metadata <- x@gene_clusters_metadata
                n_gene_cluster_annotations <-
                  x@gene_cluster_annotations
                n_cells_metadata <-
                  x@cells_metadata[j, , drop = FALSE]
                n_dbf_output <- x@dbf_output
                
                if(!is.null(n_dbf_output$center))
                  n_dbf_output$center <-
                    n_dbf_output$center[, j, drop = FALSE]
                
              } else {
                n_data <- x@data[unlist(x@gene_clusters[i]), j, drop = FALSE]
                n_gene_clusters <- x@gene_clusters[i]
                
                if (length(x@top_genes)) {
                  n_top_genes <- x@top_genes[i]
                } else{
                  n_top_genes <- x@top_genes
                }
                
                n_gene_clusters_metadata <- x@gene_clusters_metadata
                n_gene_clusters_metadata$cluster_id <-
                  x@gene_clusters_metadata$cluster_id[i]
                n_gene_clusters_metadata$number <- length(i)
                n_gene_clusters_metadata$size <-
                  x@gene_clusters_metadata$size[i]
                if (length(x@gene_cluster_annotations) > 0) {
                  n_gene_cluster_annotations <- x@gene_cluster_annotations[i]
                } else{
                  n_gene_cluster_annotations <- x@gene_cluster_annotations
                }
                n_cells_metadata <-
                  x@cells_metadata[j, , drop = FALSE]
                n_dbf_output <- x@dbf_output
                
                if(!is.null(n_dbf_output$center))
                  n_dbf_output$center <-
                    n_dbf_output$center[i, j, drop = FALSE]
              }
            }
            
            cname <- clust_names(x)
            csize <- clust_size(x)

            new(
              "ClusterSet",
              data = Matrix::Matrix(n_data, sparse=TRUE),
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
#' @export
#' @examples
#' load_example_dataset('7871581/files/pbmc3k_medium_clusters')
#' n_clust <- nclust(pbmc3k_medium_clusters)
#' @noRd
setGeneric("nclust",
           function(x)
             standardGeneric("nclust"))


#' @title Number of clusters in a ClusterSet object.
#' nclust
#' @description
#' The number of clusters in a ClusterSet object.
#' @param x The ClusterSet object
#' @export
#' @examples
#' # load a dataset
#' load_example_dataset('7871581/files/pbmc3k_medium_clusters')
#' n_clust <- nclust(pbmc3k_medium_clusters)
setMethod("nclust", signature("ClusterSet"),
          function(x) {
            length(x@gene_clusters)
          })

#' @title Sizes of the clusters stored in a ClusterSet object
#' clust_size
#' @description
#' The sizes of the clusters stored in a ClusterSet object.
#' @param x A ClusterSet object.
#' @export
#' @examples
#' # load a dataset
#' load_example_dataset('7871581/files/pbmc3k_medium_clusters')
#' clust_size(pbmc3k_medium_clusters)
#' @noRd
setGeneric("clust_size",
           function(x)
             standardGeneric("clust_size"))

#' @title Sizes of the clusters stored in a ClusterSet object
#' clust_size
#' @description
#' The sizes of the clusters stored in a ClusterSet object.
#' @param x A ClusterSet object.
#' @export
#' @examples
#' # load a dataset
#' load_example_dataset('7871581/files/pbmc3k_medium_clusters')
#' clust_size(pbmc3k_medium_clusters)
setMethod("clust_size", signature("ClusterSet"),
          function(x) {
            x@gene_clusters_metadata$size
          })

#' @title Extract centers from a ClusterSet object.
#' @description
#' Extract centers from a ClusterSet object.
#' @param x A ClusterSet object.
#' @export
#' @examples
#' # load a dataset
#' load_example_dataset('7871581/files/pbmc3k_medium_clusters')
#' centers(pbmc3k_medium_clusters)
#' @noRd
setGeneric("centers",
           function(x)
             standardGeneric("centers"))

#' @title Extract centers from a ClusterSet object.
#' @description
#' Extract centers from a ClusterSet object.
#' @param x A ClusterSet object.
#' @export
#' @examples
#' # load a dataset
#' load_example_dataset('7871581/files/pbmc3k_medium_clusters')
#' centers(pbmc3k_medium_clusters)
setMethod("centers", signature("ClusterSet"),
          function(x) {
            x@dbf_output$center
          })
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
#' @param as_string Return cluster names as strings.
#' @examples
#' # load a dataset
#' load_example_dataset('7871581/files/pbmc3k_medium_clusters')
#' g_clust <- gene_cluster(pbmc3k_medium_clusters)
#' @export
#' @noRd
setGeneric("gene_cluster",
           function(object,
                    cluster = 0,
                    as_string = FALSE)
             standardGeneric("gene_cluster"))

#' @title The gene clusters stored in a ClusterSet.
#' gene_cluster
#' @description
#' Returns a named vector (gene as names) and cluster
#' as value.
#' @param object a ClusterSet object.
#' @param cluster The cluster of interest. 0 means all cluster. Otherwise a non-null integer value.
#' @param as_string Return cluster names as strings.
#' @examples
#' # load a dataset
#' load_example_dataset('7871581/files/pbmc3k_medium_clusters')
#' g_clust <- gene_cluster(pbmc3k_medium_clusters)
#' @export
setMethod("gene_cluster", signature("ClusterSet"),
          function(object,
                   cluster = 0,
                   as_string = FALSE) {
            if (!is.null(object@gene_clusters)) {
              nb_clust <- length(object@gene_clusters)
            } else{
              print_msg("There is no cluster in this object.",
                        msg_type = 'STOP')
            }
            
            if (!is.numeric(cluster))
              print_msg("Please provide a numeric value.",
                        msg_type = 'STOP')
            
            cluster <- unique(cluster)
            
            if (!all(cluster - floor(cluster) == 0) |
                any(cluster < 0 | any(cluster > nb_clust)))
              print_msg(
                "Please provide a zero (all clusters) or positive integer in the required range.",
                msg_type = 'STOP'
              )
            
            if (length(cluster) == 1) {
              if (cluster == 0) {
                cluster <- 1:length(object@gene_clusters)
              }
            }
            
            if (length(cluster) > 1) {
              if (length(cluster[cluster == 0]))
                print_msg("Zero is out of range.",
                          msg_type = 'STOP')
            }
            
            
            if (nb_clust) {
              if (!as_string) {
                cluster_as_int <- unlist(mapply(
                  rep,
                  cluster,
                  lapply(object@gene_clusters[cluster], length),
                  SIMPLIFY = TRUE
                ))
                cluster_as_int <-
                  as.vector(as.matrix(cluster_as_int))
                names(cluster_as_int) <-
                  unlist(object@gene_clusters[cluster])
                return(cluster_as_int)
              } else{
                cluster_as_str <- unlist(mapply(
                  rep,
                  clust_names(object)[cluster],
                  lapply(object@gene_clusters[cluster], length),
                  SIMPLIFY = TRUE
                ))
                cluster_as_str <-
                  as.vector(as.matrix(cluster_as_str))
                names(cluster_as_str) <-
                  unlist(object@gene_clusters[cluster])
                return(cluster_as_str)
              }
              
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
#' @noRd
#' @export
setMethod("%in%", signature(x = "character", table = "ClusterSet"), function(x, table) {
  x %in% names(gene_cluster(table))
})


#' @title Which clusters contain a set of genes.
#' @description Which clusters contain a set of genes.
#' @param object a ClusterSet object.
#' @param genes The genes to be searched in the ClusterSet.
#' @export
#' @examples
#' # load a dataset
#' load_example_dataset('7871581/files/pbmc3k_medium_clusters')
#' hit <- which_clust(pbmc3k_medium_clusters, genes = c("TJP2", "GLA", "UNKNOWN"))
#' @noRd
setGeneric("which_clust",
           function(object,
                    genes = NULL)
             standardGeneric("which_clust"))

#' @title Which clusters contain a set of genes.
#' @description Which clusters contain a set of genes.
#' @param object a ClusterSet object.
#' @param genes The genes to be searched in the ClusterSet.
#' @examples
#' # load a dataset
#' load_example_dataset('7871581/files/pbmc3k_medium_clusters')
#' hit <- which_clust(pbmc3k_medium_clusters, genes = c("TJP2", "GLA", "UNKNOWN"))
#' @export
setMethod("which_clust",
          signature("ClusterSet"),
          function(object, genes) {
            check_format_cluster_set(object)
            gc <- gene_cluster(object)
            tmp <- gc[which(names(gc) %in% genes)]
            tmp[genes]
          })

################################################################################
##      Method for searching genes using REGEXP
################################################################################

#' @title Search gene module using a regular expression.
#' @description Search gene module using a regular expression.
#' @param object a ClusterSet object.
#' @param reg_exp The regular expression indicating the genes to be found.
#' @param as_list Whether to return the result as a list.
#' @param val if FALSE, a vector containing the (integer) indices of the matches determined
#' by grep is returned, and if TRUE, a vector containing the matching elements themselves
#' is returned.
#' @examples
#' # load a dataset
#' load_example_dataset('7871581/files/pbmc3k_medium_clusters')
#' hit <- grep_clust(pbmc3k_medium_clusters, reg_exp="^T.*[0-9]$")
#' @noRd
#' @export
setGeneric("grep_clust", 
           function(object,
                    reg_exp = NULL,
                    as_list=FALSE,
                    val=FALSE)
             standardGeneric("grep_clust")
)

#' @title Search gene module using a regular expression.
#' @description Search gene module using a regular expression.
#' @param object a ClusterSet object.
#' @param reg_exp The regular expression indicating the genes to be found.
#' @param as_list Whether to return the result as a list.
#' @param val if FALSE, a vector containing the (integer) indices of the matches determined
#' by grep is returned, and if TRUE, a vector containing the matching elements themselves
#' is returned.
#' @examples
#' # load a dataset
#' load_example_dataset('7871581/files/pbmc3k_medium_clusters')
#' hit <- grep_clust(pbmc3k_medium_clusters, reg_exp="^T.*[0-9]$")
#' @export
setMethod("grep_clust", 
          signature("ClusterSet"), 
          function(object=NULL, 
                   reg_exp=NULL,
                   as_list=FALSE,
                   val=TRUE) {
            check_format_cluster_set(object)
            grep_term <- function(x, y, val=TRUE){ grep(y, x, val=val, perl = TRUE)}
            hits <- lapply(object@gene_clusters, grep_term, reg_exp, val=val)
            if(as_list){
              return(hits)
            }else{
              hits <- stack(hits)
              hits <- setNames(hits$values, hits$ind)
              return(hits)
            }
            
            
          })


################################################################################
##      Method for renaming clusters from a clusterSet
################################################################################


#' @title Rename the gene clusters of a ClusterSet
#' @description Rename the gene clusters of a ClusterSet.
#' @param object a ClusterSet object.
#' @param new_names The new names for the clusters.
#' @export
#' @examples
#' # load a dataset
#' load_example_dataset('7871581/files/pbmc3k_medium_clusters')
#' new_obj <- rename_clust(pbmc3k_medium_clusters, new_names = letters[1:nclust(pbmc3k_medium_clusters)])
#' @noRd
setGeneric("rename_clust",
           function(object, new_names = NULL)
             standardGeneric("rename_clust"))

#' @title Rename the gene clusters of a ClusterSet
#' @description Rename the gene clusters of a ClusterSet.
#' @param object a ClusterSet object.
#' @param new_names The new names for the clusters.
#' @export
#' @examples
#' # load a dataset
#' load_example_dataset('7871581/files/pbmc3k_medium_clusters')
#' new_obj <- rename_clust(pbmc3k_medium_clusters, new_names = letters[1:nclust(pbmc3k_medium_clusters)])
setMethod("rename_clust",
          signature("ClusterSet"),
          function(object,
                   new_names = NULL) {
            check_format_cluster_set(object)
            
            if (is.null(new_names)) {
              if (length(object@gene_clusters)) {
                new_names <- 1:nclust(object)
              } else{
                new_names <- NULL
              }
            }
            
            
            if (length(new_names) != nclust(object))
              print_msg("The number of labels should be the same a the number of clusters.")
            
            
            names(object@gene_clusters) <- new_names
            
            if (length(object@top_genes) > 0)
              names(object@top_genes) <- new_names
            
            object@gene_clusters_metadata$cluster_id <- new_names
            
            names(object@gene_clusters_metadata$size) <- new_names
            
            if (length(object@gene_cluster_annotations) > 0)
              names(object@gene_cluster_annotations) <- new_names
            
            if(!is.null(object@dbf_output$center))
              rownames(object@dbf_output$center) <- new_names
            
            return(object)
            
          })


################################################################################
##      Method for writing gene list into an excel sheet.
################################################################################

#' @title Write Cluster-Set gene lists into an excel sheet.
#' @description  Write gene lists from a Cluster-Set object into an excel sheet.
#' @param object The ClusterSet object.
#' @param file_path The file path.
#' @noRd
#' @export
setGeneric("cluster_set_to_xls",
           function(object,
                    file_path = NULL)
             standardGeneric("cluster_set_to_xls"))

#' @title Write Cluster-Set gene lists into an excel sheet.
#' @description  Write gene lists from a Cluster-Set object into an excel sheet.
#' @param object The ClusterSet object.
#' @param file_path The file path.
#' @importFrom WriteXLS WriteXLS
#' @examples 
#' #' Load an example dataset
#' load_example_dataset('7871581/files/pbmc3k_medium_clusters')
#'
#' # Write gene lists to file
#' tp_dir <- tempdir()
#' dir.create(tp_dir, showWarnings = FALSE)
#' cluster_set_to_xls(pbmc3k_medium_clusters, file.path(tp_dir, "test.xls"))
#' @export
setMethod("cluster_set_to_xls",
          signature("ClusterSet"),
          function(object,
                   file_path = NULL) {
            check_format_cluster_set(object)
            object <- reorder_genes(object)
            dir_n <- dirname(file_path)
            
            if (!dir.exists(dir_n))
              print_msg("Directory does not exist. Exiting.", msg_type = "STOP")
            
            
            if (file.exists(file_path))
              print_msg("File  already exist. Exiting.", msg_type = "STOP")
            
            
            gnc <- gene_cluster(object)
            df_list <- list(x=data.frame(All_modules = unname(gnc),
                                         "official_gene_symbol" = names(gnc)))
            
            tmp <- lapply(object@gene_clusters, as.data.frame)
            
            for(i in 1:length(tmp)){
              colnames(tmp[[i]]) <- paste0("Module ", i)
              
            }
            
            df_list <- append(df_list, tmp)                           

            WriteXLS::WriteXLS(
              x=df_list,
              ExcelFileName = file_path,
              SheetNames = c("All_modules", paste0("Module ", 1:length(object@gene_clusters)))
            )
            
})



################################################################################
##      Method for reordering clusters from a clusterSet
################################################################################

#' @title Reorder the clusters from a ClusterSet
#' @description Reorder the clusters from a ClusterSet based on their names
#' @param object a ClusterSet object.
#' @param new_order The names from the clusterSet in an alternative order.
#' @export
#' @examples
#' # load a dataset
#' load_example_dataset('7871581/files/pbmc3k_medium_clusters')
#' clust_size(pbmc3k_medium_clusters)
#' new_obj <- reorder_clust(pbmc3k_medium_clusters, new_order = 15:1)
#' clust_size(pbmc3k_medium_clusters)
#' @noRd
setGeneric("reorder_clust",
           function(object, new_order = NULL)
             standardGeneric("reorder_clust"))

#' @title Reorder the clusters from a ClusterSet
#' @description Reorder the clusters from a ClusterSet based on their names
#' @param object a ClusterSet object.
#' @param new_order The names from the clusterSet in an alternative order.
#' @export
#' @examples
#' # load a dataset
#' load_example_dataset('7871581/files/pbmc3k_medium_clusters')
#' clust_size(pbmc3k_medium_clusters)
#' new_obj <- reorder_clust(pbmc3k_medium_clusters, new_order = 15:1)
#' clust_size(new_obj)
setMethod("reorder_clust",
          signature("ClusterSet"),
          function(object,
                   new_order = NULL) {
            check_format_cluster_set(object)
            
            if (is.null(new_order)) {
              print_msg('Please provide new_order argument.')
            }
            
            if (length(new_order) != nclust(object))
              print_msg("The number of labels should be the same a the number of clusters.")
            
            if (!all(sort(new_order) == sort(clust_names(object))))
              print_msg("The labels should be the same in an alternative order.")
            
            new_pos <- match(new_order, clust_names(object))
            object@gene_clusters <- object@gene_clusters[new_pos]
            
            if (length(object@top_genes) > 0)
              object@top_genes <- object@top_genes[new_pos]
            
            object@gene_clusters_metadata$cluster_id <-
              object@gene_clusters_metadata$cluster_id[new_pos]
            
            object@gene_clusters_metadata$size <-
              object@gene_clusters_metadata$size[new_pos]
            
            if (length(object@gene_cluster_annotations) > 0)
              object@gene_cluster_annotations <-
              object@gene_cluster_annotations[new_pos]
            
            if(!is.null(object@dbf_output$center))
              object@dbf_output$center <-
                object@dbf_output$center[new_pos,]
            
            return(object)
            
          })


################################################################################
##      Method for selecting a subset of column/cell for each identity
################################################################################

#' @title Given ncell, a target number, select ncell from each class of cell/column.
#' @description Given ncell, a target number, select ncell from each class of cell/column.
#' @param object a ClusterSet object.
#' @param ident A named vector. Names are cell/column names, values are classes/identity. 
#' Typically the result of the Seurat::Ident() function.
#' @param nbcell The number of cell to select.
#' @param seed A seed for subsampling.
#' @export
#' @examples
#' # load a dataset
#' @noRd
setGeneric("subsample_by_ident", 
           function(object, 
                    ident=NULL, 
                    nbcell=TRUE,
                    seed=123)
             standardGeneric("subsample_by_ident")
)

#' @title Given ncell, a target number, select ncell from each class of cell/column.
#' @description Given ncell, a target number, select ncell from each class of cell/column.
#' @param object a ClusterSet object.
#' @param ident A named vector. Names are cell/column names, values are classes/identity. 
#' Typically the result of the Seurat::Ident() function.
#' @param nbcell The number of cell to select per population (cell identity).
#' @param seed A seed for subsampling.
#' @export
#' @examples
#' # Set verbosity to 1 to display info messages only.
#' set_verbosity(1)
#' 
#' # Load a dataset
#' load_example_dataset('7871581/files/pbmc3k_medium_clusters')
#' idents <- sample(1:10, size=ncol(pbmc3k_medium_clusters), rep=TRUE)
#' names(idents) <- col_names(pbmc3k_medium_clusters)
#' sub <- subsample_by_ident(pbmc3k_medium_clusters, 
#'          ident=idents,
#'          nbcell=10)
setMethod("subsample_by_ident", 
          signature("ClusterSet"), 
          function(object, 
                   ident=NULL, 
                   nbcell=1e6,
                   seed=123) {
            
            
            check_format_cluster_set(object)
            
            print_msg(paste0("Number of cell in the object: ", ncol(object)), msg_type = "DEBUG")
            
            if(is.null(ident)){
              print_msg('Please set the ident argument.')
            }else{
              
              name_idents <- names(ident)
              
              if(is.null(name_idents)){
                print_msg("The 'ident' argument needs a named vector.")
              }
              
              if(length(which(names(ident) == "")) != 0){
                print_msg("The 'ident' argument needs a named vector or a named list of named vector.")
              }
              
              
            }
            
            cell_ident <- split(names(ident),  ident)
            
            subsample <- function(x, y, seed){
                            if(length(x) < y){
                              print_msg("Not enough cells for sampling, returning max", msg_type = "DEBUG")
                                return(x)
                            }else{
                              set.seed(seed)
                              print_msg("Sampling requested number of cells", msg_type = "DEBUG")
                              return(sample(x, size = y, replace = FALSE))
                            }
                        }
            
            cell_ident <- unlist(lapply(cell_ident, subsample, nbcell, seed))
            
            print_msg(paste0("Number of cell left: ", length(cell_ident)), msg_type = "DEBUG")
            
            object <- object[, cell_ident]

            return(object)
            
          })

################################################################################
##      Method for printing gene clusters
################################################################################

#' @title Write the cluster to files.
#' @description  Write the cluster to files.
#' @param object a ClusterSet object.
#' @param sep The separator
#' @param file_prefix A file prefix.
#' @param file_suffix A file suffix.
#' @param path A directory to store the files.
#' @param single_file Logical. Whether to write all clusters in a single file (one cluster / line). Need to change the default separator (e.g to ","). The file_prefix is used as file name.
#' @param write_cname Whether to add the cluster name. The cluster name is written as a prefix of each line in file(s) and followed by two pipes ("||"). 
#' @export
#' @examples
#' # load a dataset
#' load_example_dataset('7871581/files/pbmc3k_medium_clusters')
#' write_clust(pbmc3k_medium_clusters[1:3,], path="/tmp")
#' @noRd
setGeneric("write_clust", 
           function(object,
                    sep = "\n",
                    file_prefix="scigenex_clust",
                    file_suffix=".csv",
                    path=NULL,
                    single_file=FALSE,
                    write_cname=FALSE)
             standardGeneric("write_clust")
)

#' @title Write the cluster to files.
#' @description  Write the cluster to files.
#' @param object a ClusterSet object.
#' @param sep The separator
#' @param file_prefix A file prefix.
#' @param file_suffix A file suffix.
#' @param path A directory to store the files.
#' @param single_file Logical. Whether to write all clusters in a single file (one cluster / line). Need to change the default separator (e.g to ","). The file_prefix is used as file name.
#' @param write_cname Whether to add the cluster name. The cluster name is written as a prefix of each line in file(s) and followed by two pipes ("||"). 
#' @export
#' @examples
#' # load a dataset
#' load_example_dataset('7871581/files/pbmc3k_medium_clusters')
#' write_clust(pbmc3k_medium_clusters[1:3,], path="/tmp")
setMethod("write_clust", 
          signature("ClusterSet"), 
          function(object,
                   sep = ",",
                   file_prefix="scigenex_clust",
                   file_suffix=".csv",
                   path=NULL,
                   single_file=TRUE,
                   write_cname=TRUE) {
            
            if(is.null(path)){
              path <- getwd()
            }else{
              if(!dir.exists(path)){
                print_msg(paste0("Creating a path for output: ", 
                                 path), 
                          msg_type = "INFO")
                dir.create(path, showWarnings = FALSE, recursive = TRUE) 
                
              }
            }
            
            check_format_cluster_set(object)
            
            if(!single_file){
              
              cat_fun <- function(x, sep=NULL, file=NULL, write_cname=FALSE, clust_names){
                if(!write_cname){
                  cat(paste0(sort(x), collapse = sep), file=file)
                }else{
                  cat(paste0(clust_names, "||", paste0(sort(x), collapse = sep), sep=""), file=file)
                }
              }  
              
              for(i in 1:length(object@gene_clusters)){
                file_out <- paste0(file_prefix, "_", i, file_suffix)
                cat_fun(object@gene_clusters[[i]], 
                        sep=sep, 
                        file=file.path(path, file_out),
                        write_cname=write_cname,
                        clust_names=names(object@gene_clusters)[i])
              }
            }else{
              for(i in 1:length(object@gene_clusters)){
                if(!write_cname){
                  cat(paste0(object@gene_clusters[[i]], collapse = sep),
                      file=file.path(path, paste0(file_prefix, file_suffix)), 
                      append = TRUE,
                      sep="\n")
                }else{
                  clust_names <- names(object@gene_clusters)[i]
                  cat(paste0(clust_names, "||", paste0(object@gene_clusters[[i]], collapse = sep)),
                      file=file.path(path, paste0(file_prefix, file_suffix)), 
                      append = TRUE,
                      sep="\n")
                }
                
                
              }
            }
            
          })



#################################################################
##    Define top_by_intersect function for a ClusterSet object
#################################################################
#' @title Select top_genes based on intersection with a list.
#' @description
#' The clusterSet object contains a top_genes slot that can be used to display 
#' genes in heatmaps (see \code{plot_heatmap} function). Here the function select 
#' top_genes based on intersection with a list.
#' @param object A \code{ClusterSet} object.
#' @param set A list to compare clusters to.
#' @param as_list Return a list of clusters not a ClusterSet object.
#' @return A \code{ClusterSet} object or a list (see as_list).
#' @export
#' @noRd
#' @examples
#' # Set verbosity to 1 to display info messages only.
#' set_verbosity(1)
#' 
#' # Load a dataset
#' load_example_dataset('7871581/files/pbmc3k_medium_clusters')
#' 
#' set <- c('MS4A1', 'ISG20', 'CD3D', 'SEC14L5', 'RPL11', 'RPL32')
#' pbmc3k_medium_clusters <- top_by_intersect(pbmc3k_medium_clusters, set=set)
#' pbmc3k_medium_clusters@top_genes
setGeneric("top_by_intersect", 
           function(object,
                    set=NULL,
                    as_list=FALSE)
             standardGeneric("top_by_intersect")
)

#################################################################
##    Define top_by_intersect function for a ClusterSet object
#################################################################
#' @title Select top_genes based on intersection with a list.
#' @description
#' The clusterSet object contains a top_genes slot that can be used to display 
#' genes in heatmaps (see \code{plot_heatmap} function). Here the function select 
#' top_genes based on intersection with a list.
#' @param object A \code{ClusterSet} object.
#' @param set A list to compare clusters to.
#' @param as_list Return a list of clusters not a ClusterSet object.
#' @return A \code{ClusterSet} object or a list (see as_list).
#' @export
#' @examples
#' # Set verbosity to 1 to display info messages only.
#' set_verbosity(1)
#' 
#' # Load a dataset
#' load_example_dataset('7871581/files/pbmc3k_medium_clusters')
#' 
#' set <- c('MS4A1', 'ISG20', 'CD3D', 'SEC14L5', 'RPL11', 'RPL32')
#' pbmc3k_medium_clusters <- top_by_intersect(pbmc3k_medium_clusters, set=set)
#' pbmc3k_medium_clusters@top_genes
setMethod("top_by_intersect", 
          signature("ClusterSet"), 
          function(object,
                   set=NULL,
                   as_list=FALSE) {
            
            if(is.null(set))
              print_msg("Please provide a set ('set' argument)", 
                        msg_type = "STOP")
            
            top_gn <- lapply(object@gene_clusters, intersect, set)
            
            if(as_list){
              return(top_gn)
            }else{
              object@top_genes <- top_gn
              return(object)
            }
            
          })


#################################################################
##    Define top_by_grep function for a ClusterSet object
#################################################################
#' @title Select top_genes based on a regular expression search
#' @description
#' The clusterSet object contains a top_genes slot that can be used to display 
#' genes in heatmaps (see \code{plot_heatmap} function). Here the function select 
#' top_genes based on a regular expression search
#' @param object A \code{ClusterSet} object.
#' @param regexp A regular expression
#' @param as_list Return a list of clusters not a ClusterSet object.
#' @return A \code{ClusterSet} object or a list (see as_list).
#' @export
#' @noRd
#' @examples
#' # Set verbosity to 1 to display info messages only.
#' set_verbosity(1)
#' 
#' # Load a dataset
#' load_example_dataset('7871581/files/pbmc3k_medium_clusters')
#' 
#' pbmc3k_medium_clusters <- top_by_grep(pbmc3k_medium_clusters, regexp="^CD")
#' pbmc3k_medium_clusters@top_genes
setGeneric("top_by_grep", 
           function(object,
                    regexp=NULL,
                    as_list=FALSE)
             standardGeneric("top_by_grep")
)

#################################################################
##    Define top_by_grep function for a ClusterSet object
#################################################################
#' @title Select top_genes based on a regular expression search
#' @description
#' The clusterSet object contains a top_genes slot that can be used to display 
#' genes in heatmaps (see \code{plot_heatmap} function). Here the function select 
#' top_genes based on a regular expression search
#' @param object A \code{ClusterSet} object.
#' @param regexp A regular expression
#' @param as_list Return a list of clusters not a ClusterSet object.
#' @return A \code{ClusterSet} object or a list (see as_list).
#' @export
#' @examples
#' # Set verbosity to 1 to display info messages only.
#' set_verbosity(1)
#' 
#' # Load a dataset
#' load_example_dataset('7871581/files/pbmc3k_medium_clusters')
#' 
#' pbmc3k_medium_clusters <- top_by_grep(pbmc3k_medium_clusters, regexp="^CD")
#' pbmc3k_medium_clusters@top_genes
setMethod("top_by_grep", 
          signature("ClusterSet"), 
          function(object,
                   regexp=NULL,
                   as_list=FALSE) {
            
            if(is.null(regexp))
              print_msg("Please provide a regexp ('regexp' argument)", 
                        msg_type = "STOP")
            
            fun_grep <- function(x, regexp, val=TRUE, perl=TRUE) grep(regexp, x, value=val, perl=perl)
            top_gn <- lapply(object@gene_clusters, fun_grep, regexp, val=TRUE, perl=TRUE)
            
            if(as_list){
              return(top_gn)
            }else{
              object@top_genes <- top_gn
              return(object)
            }
            
          })


#################################################################
##    Define compute_centers() function for a ClusterSet object
#################################################################
#' @title Compute centers of gene modules (e.g mean profiles)
#' @description
#' Compute centers of gene modules (e.g mean profiles)
#' @param object A \code{ClusterSet} object.
#' @export
#' @noRd
#' @examples
#' # Set verbosity to 1 to display info messages only.
#' set_verbosity(1)
#' 
#' # Load a dataset
#' load_example_dataset('7871581/files/pbmc3k_medium_clusters')
#' pbmc3k_medium_clusters <- compute_centers(pbmc3k_medium_clusters)
#' pbmc3k_medium_clusters@dbf_output$center
setGeneric("compute_centers", 
           function(object)
             standardGeneric("compute_centers")
)

#' @title Compute centers of gene modules (e.g mean profiles)
#' @description
#' Compute centers of gene modules (e.g mean profiles)
#' @param object A \code{ClusterSet} object.
#' @export
#' @examples
#' # Set verbosity to 1 to display info messages only.
#' set_verbosity(1)
#' 
#' # Load a dataset
#' load_example_dataset('7871581/files/pbmc3k_medium_clusters')
#' pbmc3k_medium_clusters <- compute_centers(pbmc3k_medium_clusters)
#' pbmc3k_medium_clusters@dbf_output$center
#' @importFrom Matrix Matrix
#' @importFrom Matrix colMeans
setMethod("compute_centers", 
          signature("ClusterSet"), 
          function(object) {
            
            print_msg("Computing centers.", msg_type = "INFO")
            
            nb_clusters <- nclust(object)
            
            print_msg("Preparing matrix.", msg_type = "DEBUG")
            
            centers <- matrix(NA, 
                              ncol = ncol(object),
                              nrow = nb_clusters)
            
            print_msg("Renaming columns / rows.", msg_type = "DEBUG")
            
            colnames(centers) <- colnames(object@data)
            rownames(centers) <- names(object@gene_clusters)
            
            print_msg("Looping over clusters...", msg_type = "DEBUG")
            
            for (i in 1:nb_clusters) {
              print_msg(paste0("Computing cluster ", i, " center."), msg_type = "DEBUG")
              tmp <- Matrix::colMeans(object@data[object@gene_clusters[[i]], , drop=FALSE],
                                      na.rm = TRUE)
              centers[i, ] <- setNames(tmp, NULL)
              
            }

            centers <- Matrix::Matrix(centers)
            object@dbf_output$center <- centers
            
            return(object)
})
