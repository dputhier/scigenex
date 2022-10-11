#################################################################
##    DEFINITION OF A SPECIFIC CLASS OBJECT : ClusterSet
#################################################################

library(ggplot2)
library(reshape2)
library(dplyr)
library(igraph)
library(iheatmapr)

#' @title
#' ClusterSet
#' @description
#' This class is a representation of a partitioning algorithm and is intented to store gene clusters.
#' @slot name character. The original input file name (if applicable).
#' @slot data matrix. The matrix containing the filtered/partitionned data.
#' @slot distances vector. The observed distance values with the knn.
#' @slot simulated_distances vector. The simulated distance values with the knn.
#' @slot critical_distance vector. The critical threshold distance to select informative genes.
#' @slot gene_patterns vector. Mapping of row/genes to gene_patterns.
#' @slot size vector. The size of each cluster.
#' @slot dot_product vector. The median dot product of each gene clusters.
#' @slot top_genes matrix The highly co-expressed genes of each gene clusters.
#' @slot center matrix. The centers of each clusters.
#' @slot parameters list. The parameter used.
#' @slot algorithm vector. The algorithm used to produce the clusters.
#' @slot cell_clusters list The cell clusters.
#' @slot cell_types vector. The cell types.
#' @slot cell_colors vector. The cell types to color mapping.
#' @slot cell_order vector. How cell should be ordered.
#' @slot cluster_annotations list. Functional annotation of clusters.
#' @return A ClusterSet object.
#' @export
#'
#' @examples
#' 
#' \dontrun{
#'   m <- matrix(rnorm(80000), nc=20)
#'    m[1:100,1:10] <- m[1:100,1:10] + 4
#'    m[101:200,11:20] <- m[101:200,11:20] + 3
#'    m[201:300,5:15] <- m[201:300,5:15] + -2
#'    res <- find_gene_clusters(data=m,
#'                              distance_method="pearson",
#'                              av_dot_prod_min = 0,
#'                              inflation = 1.2,
#'                              k=25,
#'                              fdr = 10)
#' plot_clust(res, ceil = 10, floor = -10)
#' plot_clust(res, type="tile", ceil = 10, floor = -10)
#' write_clust(res, filename_out = "ALL.sign.txt")
#'   is(res)
#' }
#'               
setClass("ClusterSet",
         representation = list(
           name = "character",
           data = "matrix",
           distances = "vector",
           simulated_distances = "vector",
           critical_distance = "vector",
           gene_patterns = "vector",
           size = "vector",
           dot_product = "vector",
           top_genes = "matrix",
           center = "matrix",
           parameters = "list",
           algorithm = "vector",
           cell_clusters = "list",
           cell_types = "vector",
           cell_colors = "vector",
           cell_order = "vector",
           cluster_annotations = "list"
         ),
         prototype = list(
           name = character(),
           data = matrix(nr = 0, nc = 0),
           distances = vector(),
           simulated_distances = vector(),
           critical_distance = vector(),
           gene_patterns = numeric(),
           size = numeric(),
           dot_product = numeric(),
           top_genes = matrix(),
           center = matrix(nc = 0, nr = 0),
           parameters = list(),
           algorithm = character(),
           cell_clusters = list(),
           cell_types = vector(),
           cell_colors = vector(),
           cell_order = vector(),
           cluster_annotations = list()
         )
)

#################################################################
##    REDEFINE SHOW METHOD FOR CLASS OBJECT : ClusterSet
#################################################################

setMethod(
  "show", signature("ClusterSet"),
  
  function(object) {
    
    cat("\t\tAn object of class ClusterSet\n")
    cat("\t\tName:", slot(object, "name"), "\n")
    cat("\t\tMemory used: ", object.size(object), "\n")
    cat("\t\tNumber of cells: ", ncol(slot(object, "data")), "\n")
    cat(
      "\t\tNumber of informative genes: ",
      nrow(slot(object, "data")), "\n"
    )
    cat("\t\tNumber of clusters: ", length(slot(object, "size")), "\n")
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