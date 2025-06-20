#' @title Transform a Seurat objects / FindAllMarkers result into a ClusterSet.
#' @description Transform a Seurat objects into a ClusterSet.
#' @param object A Seurat object.
#' @param markers A Seurat::FindAllMarkers() result or a named vector (clusters with gene_names as named).
#' @param layer One of 'data', 'counts' or 'scale.data'. The slot to extract from the seurat object to perform clustering analysis.
#' SCT is the recommended method from Seurat package when working with spatial transcriptomics data.
#' @param assay The type of assay (e.g. "RNA", "Spatial", "Sketch", "SCT"...).
#' @importFrom SeuratObject LayerData
#' @importFrom Matrix Matrix
#' @examples
#' ## From a scRNA-seq/Seurat object
#' library(SeuratObject)
#' library(Seurat)
#' data("pbmc_small", package="SeuratObject")
#' cs <- cluster_set_from_seurat(pbmc_small, Seurat::FindAllMarkers(pbmc_small))
#' plot_heatmap(cs)
#' plot_heatmap(cs)
#' plot_heatmap(cs[1,])
#' plot_heatmap(cs, cell_clusters = Seurat::Idents(pbmc_small))
#' plot_heatmap(cs[1,Idents(pbmc_small) == "0"], 
#'              cell_clusters = Seurat::Idents(pbmc_small), label_size = 6)
#' plot_profiles(cs, ident = Seurat::Idents(pbmc_small))
#' @export
cluster_set_from_seurat <- function(object=NULL, 
                                    markers=NULL,
                                    layer=c('data', 'counts', 'scale.data'),
                                    assay='RNA'){
  
  layer <- match.arg(layer)
  
  object <- SeuratObject::LayerData(object, assay=assay, layer=layer)
  
  if(inherits(markers, "data.frame")){
    
    gn <- markers$gene
    clusters <- markers$cluster
    names(clusters) <- gn
    object <- object[markers$gene, ]
    
  }else if(is.vector(markers)){
    if(is.null(names(markers)))
      print_msg("The 'markers' argument should be a named vector.",
                msg_type="STOP")
    clusters <- markers
  }else{
    print_msg("The 'markers' argument should be a data.frame or named vector.",
              msg_type="STOP")
  }
  
  gn <- split(names(clusters), clusters)

  obj_out <- new(Class = "ClusterSet")
  obj_out@gene_clusters <- gn
  obj_out@data <- Matrix::Matrix(object, sparse = TRUE)
  obj_out@gene_clusters_metadata <- list("cluster_id" = setNames(as.character(unique(clusters)), 
                                                                 as.character(unique(clusters))),
                                          "number" = length(table(clusters)),
                                          "size" = table(clusters))
  obj_out@parameters    <- list("distance_method" = "pearson",
                                "k" = vector(),
                                "noise_level" = vector(),
                                "fdr" = vector(),
                                "row_sum" = vector(),
                                "no_dknn_filter" = vector(),
                                "seed" = vector())
  
  obj_out@dbf_output <- list("dknn" = vector(),
                         "simulated_dknn" = vector(),
                         "critical_distance" = vector(),
                         "fdr" = vector(),
                         "center" = NULL,
                         "all_gene_expression_matrix" = vector(),
                         "all_neighbor_distances" = vector())
  
  return(obj_out)
}



#' @title Transform any matrix and list into a ClusterSet object.
#' @description 
#' Transform any matrix (e.g expression matrix) and list (e.g markers obtained through kmeans 
#' or any partitioning algorithm)  into a ClusterSet object.
#' @param object A matrix or data.frame.
#' @param markers A list of vector containing the gene sets
#' @examples
#' m <- create_3_rnd_clust()[1:300,] 
#' rownames(m) <- paste0("gene", 1:300)
#' markers <- list(a=paste0("gene", 1:100), 
#'                 b=paste0("gene", 101:200),
#'                 c=paste0("gene", 201:300))
#' cs <- cluster_set_from_matrix(m, markers)
#' plot_heatmap(cs, interactive = FALSE)
#' 
#' @export
cluster_set_from_matrix <- function(object=NULL, 
                                    markers=NULL){
  
  if(!(inherits(object, "data.frame") |  inherits(object, "matrix"))){
    print_msg("The 'object' argument should be a data.frame or matrix.",
              msg_type="STOP")
  }
  
 
  if(!inherits(markers, "list")){
    print_msg("The 'marker' argument should be a list.",
              msg_type="STOP")
  }
  
  print_msg("Checking markers.", msg_type = "INFO")
  marker_found <- unlist(markers)[unlist(markers) %in% rownames(object)]
  
  if(length(marker_found) == 0){
    print_msg("No marker were found in the matrix.")
    return(new(Class = "ClusterSet"))
  }
  
  print_msg("Subsetting object.", msg_type = "INFO")
  object <- object[marker_found, ]
  
  print_msg("Subsetting object.", msg_type = "INFO")
  select_markers <- function(x, y){ x[x %in% y] }
  markers <- lapply(markers, "select_markers", marker_found)
  clusters <- unlist(mapply(FUN = "rep", 1:length(markers), 
                     unlist(lapply(markers, length))))

  print_msg("Computing centers.", msg_type = "INFO")

  print_msg("Creating ClusterSet object.", msg_type = "INFO")
  obj_out <- new(Class = "ClusterSet")
  obj_out@gene_clusters <- markers
  obj_out@data <- Matrix::Matrix(object, sparse = TRUE)
  obj_out@gene_clusters_metadata <- list("cluster_id" = setNames(as.character(unique(clusters)), 
                                                                 as.character(unique(clusters))),
                                         "number" = length(table(clusters)),
                                         "size" = table(clusters))

  obj_out@parameters    <- list("distance_method" = "unknown",
                                "k" = vector(),
                                "noise_level" = vector(),
                                "fdr" = vector(),
                                "row_sum" = vector(),
                                "no_dknn_filter" = vector(),
                                "seed" = vector())
  
  obj_out@dbf_output <- list("dknn" = vector(),
                             "simulated_dknn" = vector(),
                             "critical_distance" = vector(),
                             "fdr" = vector(),
                             "center" = NULL,
                             "all_gene_expression_matrix" = vector(),
                             "all_neighbor_distances" = vector())
  
  return(obj_out)
}

