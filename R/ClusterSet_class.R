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
#' @slot data A matrix containing the filtered and partitioned data.
#' @slot gene_clusters A list contains the partitioned genes of the dataset. Each element of the list corresponds to a cluster, and contains the indices of the genes assigned to that cluster.
#' @slot top_genes A list contains the top genes from the gene clusters. Each element of the list corresponds to a cluster, and contains the indices of the genes assigned to that cluster ranked by their correlation value within their cluster.
#' @slot gene_clusters_metadata A list contains metadata related to the gene clusters such as the number of gene clusters, their ID, and the number of genes contained in each of them.
#' @slot gene_cluster_annotations A list contains the result obtained from the GO enrichment analysis of gene clusters.
#' @slot cell_metadata A list containing metadata related to the cell clusters such as the clustering results the number of cell clusters, their order, colors associated to each cluster,... 
#' @slot dbf_output A list containing the intermediates outputs of the DBF function : dknn, simulated distances, critical distance and fdr values. 
#' @slot parameters A list containing the parameter used. Each element of the list correspond to a parameter.
#' 
#' @return A ClusterSet object.
#' @export
#'
#' @examples
#' 
#' \dontrun{
#'  m <- matrix(rnorm(80000), nc=20)
#'  m[1:100,1:10] <- m[1:100,1:10] + 4
#'  m[101:200,11:20] <- m[101:200,11:20] + 3
#'  m[201:300,5:15] <- m[201:300,5:15] + -2
#'  res <- find_gene_clusters(data=m,
#'                            distance_method="pearson",
#'                            inflation = 1.2,
#'                            k=25,
#'                            fdr = 10)
#'  is(res)
#'  res
#'               
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
if (!isGeneric("ncol"))
  setGeneric("ncol", 
             function(x)
               standardGeneric("ncol")
  )

setMethod("ncol",signature(x="ClusterSet"),
          function(x) {
            ncol(x@data)
          }
) 

if (!isGeneric("nrow"))
  setGeneric("nrow", 
             function(x)
               standardGeneric("nrow")
  )

setMethod("nrow",signature(x="ClusterSet"),
          function(x) {
            nrow(x@data)
          }
)

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
