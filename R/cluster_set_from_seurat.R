#' @title Transform a Seurat objects into a ClusterSet.
#' @description Transform a Seurat objects into a ClusterSet.
#' @param object A Seurat object.
#' @param markers A Seurat::FindAllMarkers() result or a named vector (clusters with gene_names as named).
#' @param which_slot One of 'data', 'counts' or 'sct'. The slot to extract from the seurat object to perform clustering analysis.
#' SCT is the recommended method from Seurat package when working with spatial transcriptomics data.
#' @importFrom SeuratObject GetAssayData
#' @examples
#' ## From a scRNA-seq/Seurat object
#' data("pbmc_small")
#' cs <- cluster_set_from_seurat(pbmc_small, FindAllMarkers(pbmc_small))
#' plot_heatmap(cs)
#' plot_heatmap(cs)
#' plot_heatmap(cs[1,])
#' plot_heatmap(cs, cell_clusters = Idents(pbmc_small))
#' plot_heatmap(cs[1,Idents(pbmc_small) == "0"], 
#' cell_clusters = Idents(pbmc_small), label_size = 6)
#' plot_profiles(cs, ident = Idents(pbmc_small))
#' @export cluster_set_from_seurat
cluster_set_from_seurat <- function(object=NULL, 
                                    markers=NULL,
                                    which_slot=c('data', 'counts', 'sct')){
  
  which_slot <- match.arg(which_slot)
  
  if (which_slot %in% c("data", "counts")) {
    object <- SeuratObject::GetAssayData(object, slot = which_slot)
  } else if (which_slot == "sct") {
    if ("SCT" %in% names(object@assays)) {
      object <- object@assays$SCT@data
    } else{
      print_msg("This object has no 'SCT' slot. Use SCTransform() before.",
                msg_type = 'STOP')
    }
  }
  
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
  centers <- split(data.frame(object[names(clusters), ]), clusters)
  centers <- lapply(centers, apply, 2, mean)
  name_center <- names(centers)
  centers <- do.call("rbind", centers)
  rownames(centers) <- name_center
  colnames(centers) <- colnames(object)

  
  obj_out <- new(Class = "ClusterSet")
  obj_out@gene_clusters <- gn
  obj_out@data <- as.matrix(object)
  obj_out@gene_clusters_metadata <- list("cluster_id" = setNames(as.character(unique(clusters)), 
                                                                 as.character(unique(clusters))),
                                          "number" = length(table(clusters)),
                                          "size" = table(clusters))
  obj_out@parameters    <- list("distance_method" = "FindAllMarkers",
                                "k" = vector(),
                                "highest" = vector(),
                                "fdr" = vector(),
                                "row_sum" = vector(),
                                "no_dknn_filter" = vector(),
                                "seed" = vector())
  
  obj_out@dbf_output <- list("dknn" = vector(),
                         "simulated_dknn" = vector(),
                         "critical_distance" = vector(),
                         "fdr" = vector(),
                         "center" = centers,
                         "all_gene_expression_matrix" = vector(),
                         "all_neighbor_distances" = vector())
  
  return(obj_out)
}