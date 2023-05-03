# -------------------------------------------------------------------------
# Define the reorder_genes function ---------------------------------------
# -------------------------------------------------------------------------
#' @title Reorder genes in a ClusterSet object.
#' @description
#' Reorder genes in a ClusterSet object according to various
#' ordering rules.
#' @param object A ClusterSet object.
#' @param order_by A character string specifying the ordering rule. 
#' Default is "gene_names" which reorder genes in alphanumerical 
#' order. If 'hclust', hierarchical clustering is applied. If 
#' 'correlation' the gene as ordered based on there similarity 
#' with the cluster center.
#' @param nb_proc The number of processor to be used if order_by=='hclust'.
#' @param method The distance measure to be used order_by=='hclust'.
#' @param link The agglomeration method to be used if order_by=='hclust'.
#'
#' @return ClusterSet-class object
#' @export reorder_genes
#' @keywords internal
setGeneric("reorder_genes", 
           function(object=NULL, 
                    order_by = NULL,
                    method=NULL,
                    link=NULL,
                    nb_proc=2)
             standardGeneric("reorder_genes")
)

#' @title Reorder genes in a ClusterSet object.
#' @description
#' Reorder genes in a ClusterSet object according to various
#' ordering rules.
#' @param object A ClusterSet object.
#' @param order_by A character string specifying the ordering rule. 
#' Default is "gene_names" which reorder genes in alphanumerical 
#' order. If 'hclust', hierarchical clustering is applied. If 
#' 'correlation' the gene as ordered based on there similarity 
#' with the cluster center.
#' @param nb_proc The number of processor to be used if order_by=='hclust'.
#' @param method The distance measure to be used order_by=='hclust'.
#' @param link The agglomeration method to be used if order_by=='hclust'.
#'
#' @return ClusterSet-class object
#' @importFrom amap hcluster
#' @export reorder_genes
#'
#' @examples
#' # Set verbosity to 1 to display info messages only.
#' set_verbosity(1)
#' # load a dataset
#' load_example_dataset('7871581/files/pbmc3k_medium_clusters')
#'                        
#' # Reorder genes using hierarchical clustering
#' pbmc3k_medium_clusters <- reorder_genes(pbmc3k_medium_clusters, order_by = "hclust")
#' 
#' # Plot heatmap of gene clusters
#' plot_heatmap(pbmc3k_medium_clusters)
#' 
#' # Reorder and plot
#' plot_heatmap(reorder_genes(pbmc3k_medium_clusters[1,], 
#'                            order_by = "gene_names"), 
#'              label_size = 5)
#' 
setMethod("reorder_genes", 
          signature("ClusterSet"), 
          function(object, 
                   order_by = c("gene_names", "hclust", "correlation"),
                   method=c("pearson", "maximum", "manhattan", "canberra", "binary", "euclidean", 
                            "abspearson", "correlation", "abscorrelation", "spearman", "kendall"),
                   link=c("average", "single", "complete", "ward", "mcquitty", "median", "centroid",
                          "centroid2"),
                   nb_proc=2) {
  
  order_by <- match.arg(order_by)
  method <- match.arg(method)
  link <- match.arg(link)
  
  ## Check format object arg
  check_format_cluster_set(object)

  reordered_genes <- vector()
  for (gene_cluster in object@gene_clusters_metadata$cluster_id) {
    
    # Reorder by gene names
    if (order_by == "gene_names") {
      gene_cluster_names <- object@gene_clusters[[gene_cluster]]
      gene_cluster_names <- sort(gene_cluster_names)
      
      # Reorder genes in gene_clusters slot
      object@gene_clusters[[gene_cluster]] <- gene_cluster_names
      
      # Reorder genes in data slot
      object@data <- object@data[unlist(object@gene_clusters, use.names = FALSE),]
    }
    
    # Reorder using hierarchical clustering (amap::hcluster)
    if (order_by == "hclust") {
      set.seed(123)
      gene_cluster_names <- object@gene_clusters[[gene_cluster]]
      expression_matrix <- object@data[gene_cluster_names,]
      
      hclust_res <- amap::hcluster(expression_matrix,
                             method = method,
                             diag = FALSE,
                             upper = FALSE,
                             link = link,
                             nbproc = nb_proc,
                             doubleprecision = TRUE)
      
      # Reorder genes in gene_clusters slot
      object@gene_clusters[[gene_cluster]] <- hclust_res$labels[hclust_res$order]
      
      # Reorder genes in data slot
      object@data <- object@data[unlist(object@gene_clusters, use.names = FALSE),]
    }
    
    # Reorder using top_genes
    if (order_by == "correlation") {
      gene_cluster_names <- object@gene_clusters[[gene_cluster]]
      gene_cluster_names <- top_genes(object,
                                      cluster = gene_cluster,
                                      top = length(gene_cluster_names))@top_genes[[gene_cluster]]
      
      # Reorder genes in gene_clusters slot
      object@gene_clusters[[gene_cluster]] <- gene_cluster_names
      
      # Reorder genes in data slot
      object@data <- object@data[unlist(object@gene_clusters, use.names = FALSE),]
    }
  }
  
  
  return(object)
})
