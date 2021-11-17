#################################################################
##    Define the load_seurat function for class ClusterSet
#################################################################

#' @title
#' load_seurat
#' @description
#' Load a seurat object into a ClusterSet object. At the moment the objective is mainly
#' to store cell identity (i.e cell types/groups to barcode mapping) and
#' cell type to color mapping.
#' @param object A ClusterSet object.
#' @param seurat_obj A seurat object to extract cell to group mapping.
#' @param dimplot_obj To display cell clusters on the profile diagram with colors extracted from Dimplot output (provide also a Seurat object).
#
#' @return A ClusterSet object.
#' @export
#'
#' @examples
#' # see online examples
setGeneric("load_seurat",
           function(object,
                    seurat_obj=NULL,
                    dimplot_obj=NULL) {
             standardGeneric("load_seurat")
           })

#' @rdname load_seurat
setMethod("load_seurat",
          signature(object = "ClusterSet"),
          function(object,
                   seurat_obj=NULL,
                   dimplot_obj=NULL) {
            if (!inherits(seurat_obj, "Seurat")) {
              stop("Please provide a Seurat and patchwork object.")
            }
            if (!inherits(dimplot_obj, "patchwork")) {
              stop("Please provide a Seurat and patchwork object.")
            }
            if(ncol(object) != ncol(seurat_obj)){
              stop("The number of cells is not the same in ClusterSet object and Seurat object.")
            }
            
            g <- ggplot_build(dimplot_obj)
            tmp_mat <- distinct(as.data.frame(cbind(g$data[[1]]$colour,
                                                    as.character(g$plot$data$ident))))
            cell_col_tmp <- tmp_mat[,1]
            cell_grp_tmp <- tmp_mat[,2]
            object@cell_colors <- setNames(as.character(cell_col_tmp), cell_grp_tmp)
            
            tmp_mat <- distinct(as.data.frame(cbind(rownames(g$plot$data), as.character(g$plot$data$ident))))
            object@cell_types <- setNames(as.character(tmp_mat[,2]), tmp_mat[,1])
            
            object@cell_order <- rownames(g$plot$data)[order(g$plot$data$ident)]
            
            return(object)
          }
)