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
#' @param method The distance measure to be used order_by=='hclust'.
#' @param link The agglomeration method to be used if order_by=='hclust'.
#' @param decreasing Whether the sorting should be decreasing (FALSE by default)
#' @param nb_proc The number of processor to be used if order_by=='hclust'.
#'
#' @return ClusterSet-class object
#' @export reorder_genes
#' @keywords internal
setGeneric("reorder_genes", 
           function(object=NULL, 
                    order_by = NULL,
                    method=NULL,
                    link=NULL,
                    decreasing=FALSE,
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
#' @param method The distance measure to be used order_by=='hclust'.
#' @param link The agglomeration method to be used if order_by=='hclust'.
#' @param decreasing Whether the sorting should be decreasing (FALSE by default)
#' @param nb_proc The number of processor to be used if order_by=='hclust'.
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
                   method=c("pearson", "maximum", "manhattan", 
                            "canberra", "binary", "euclidean", 
                            "abspearson", "correlation", 
                            "abscorrelation", "spearman", "kendall"),
                   link=c("average", "single", "complete", "ward", 
                          "mcquitty", "median", "centroid",
                          "centroid2"),
                   decreasing=FALSE,
                   nb_proc=2) {
  
  order_by <- match.arg(order_by)
  method <- match.arg(method)
  link <- match.arg(link)
  
  ## Check format object arg
  check_format_cluster_set(object)

  if(order_by == "gene_names"){
    for (cur_clust in names(object@gene_clusters)) {

    gene_cluster_names <- object@gene_clusters[[cur_clust]]
    gene_cluster_names <- sort(gene_cluster_names, decreasing=decreasing)
    
    object@gene_clusters[[cur_clust]] <- gene_cluster_names
    }
    
  }else if(order_by == "hclust"){
    for (cur_clust in names(object@gene_clusters)) {
      
    gene_cluster_names <- object@gene_clusters[[cur_clust]]
    expression_matrix <- object@data[gene_cluster_names,]
    
    set.seed(123)
    hclust_res <- amap::hcluster(expression_matrix,
                                 method = method,
                                 diag = FALSE,
                                 upper = FALSE,
                                 link = link,
                                 nbproc = nb_proc,
                                 doubleprecision = TRUE)
    
    # Reorder genes in gene_clusters slot
    object@gene_clusters[[cur_clust]] <- hclust_res$labels[hclust_res$order]
    
    }
  }else if (order_by == "correlation") {
    clust_ordered <- top_genes(object, top=max(clust_size(object)))@top_genes
    object@gene_clusters <- clust_ordered
  }
  
  # Reorder genes in data slot
  object@data <- object@data[unlist(object@gene_clusters, use.names = FALSE),]
  
  return(object)
})
