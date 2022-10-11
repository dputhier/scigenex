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
#' m <- matrix(rnorm(80000), nc=20)
#' m[1:100,1:10] <- m[1:100,1:10] + 4
#' m[101:200,11:20] <- m[101:200,11:20] + 3
#' m[201:300,5:15] <- m[201:300,5:15] + -2
#' res <- DBFMCL(data=m,
#'               distance_method="pearson",
#'               av_dot_prod_min = 0,
#'               inflation = 2,
#'               k=25,
#'              fdr = 10)
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
  for (gene_cluster in unique(object@gene_patterns)) {
    
    # Reorder by gene names
    if (order_by == "gene_names") {
      gene_cluster_names <- names(object@gene_patterns[object@gene_patterns == gene_cluster])
      gene_cluster_names <- sort(gene_cluster_names)
      reordered_genes <- c(reordered_genes, gene_cluster_names)
    }
    
    # Reorder using hierarchical clustering (amap::hcluster)
    if (order_by == "hclust") {
      gene_cluster_names <- names(object@gene_patterns[object@gene_patterns == gene_cluster])
      expression_matrix <- object@data[gene_cluster_names,]
      
      
      hclust_res <- hcluster(expression_matrix,
                             method = "pearson",
                             diag = FALSE,
                             upper = FALSE,
                             link = "average",
                             nbproc = 2,
                             doubleprecision = TRUE)
      
      reordered_genes <- c(reordered_genes, hclust_res$labels[hclust_res$order])
    }
    
    # Reorder using top_genes
    if (order_by == "correlation") {
      gene_cluster_names <- names(object@gene_patterns[object@gene_patterns == gene_cluster])
      gene_cluster_names <- c(top_genes(object,
                                        cluster = gene_cluster,
                                        top = length(gene_cluster_names))@top_genes)
      
      reordered_genes <- c(reordered_genes, gene_cluster_names)
    }
  }
  
  object@gene_patterns <- object@gene_patterns[reordered_genes]
  
  return(object)
}

