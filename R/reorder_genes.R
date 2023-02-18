#################################################################
##    Define the reorder_genes function
#################################################################

#' @title
#' reorder_genes
#' @description
#' Reorder genes in gene cluster.
#' @param object A ClusterSet object.
#' @param order_by A character string related to the reordering type. 
#' Default is "correlation" which reorder genes based on their correlation 
#' within their respective gene cluster.
#'
#' @return ClusterSet-class object
#' @export
#'
#' @examples
#' \dontrun{
#' m <- create_3_rnd_clust()
#' 
#' res <- find_gene_clusters(data=m,
#'                              distance_method="pearson",
#'                              inflation = 2,
#'                              k=75,
#'                              row_sum=-Inf,
#'                              highest=0.3,
#'                              min_nb_supporting_cell = 0,
#'                              fdr = 1e-8)
#'                              
#' res <- reorder_genes(res, "correlation")
#' 
#' plot_heatmap(res,
#'              use_core_cells = FALSE,
#'              use_top_genes = FALSE)
#' 
#' }
#' 
#' @rdname reorder_genes

reorder_genes <- function(object,
                          order_by = c("gene_names", "hclust", "correlation")) {
  
  
  
  
  if (!order_by %in% c("gene_names", "hclust", "correlation")) {
    stop("Provided order_by parameter does not exist.")
  }
  
  
  reordered_genes <- c()
  for (gene_cluster in object@gene_clusters_metadata$cluster_id) {
    
    # Reorder by gene names
    if (order_by == "gene_names") {
      gene_cluster_names <- object@gene_clusters[[gene_cluster]]
      gene_cluster_names <- sort(gene_cluster_names)
      object@gene_clusters[[gene_cluster]] <- gene_cluster_names
    }
    
    # Reorder using hierarchical clustering (amap::hcluster)
    if (order_by == "hclust") {
      gene_cluster_names <- object@gene_clusters[[gene_cluster]]
      expression_matrix <- object@data[gene_cluster_names,]
      
      
      hclust_res <- hcluster(expression_matrix,
                             method = "pearson",
                             diag = FALSE,
                             upper = FALSE,
                             link = "average",
                             nbproc = 2,
                             doubleprecision = TRUE)
      
      object@gene_clusters[[gene_cluster]] <- hclust_res$labels[hclust_res$order]
    }
    
    # Reorder using top_genes
    if (order_by == "correlation") {
      gene_cluster_names <- object@gene_clusters[[gene_cluster]]
      gene_cluster_names <- top_genes(object,
                                      cluster = gene_cluster,
                                      top = length(gene_cluster_names))@top_genes[[gene_cluster]]
      
      object@gene_clusters[[gene_cluster]] <- gene_cluster_names
    }
  }
  
  
  return(object)
}

