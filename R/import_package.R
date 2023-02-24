# @importClassesFrom methods character list matrix vector numeric
#'  
#' @importMethodsFrom methods show
#' 
#' @importFrom methods new slot slotNames setClass signature setMethod setGeneric
#' @importFrom grDevices colorRampPalette
#' @importFrom stats sd dist hclust heatmap median rnorm setNames as.dist hclust cor pnorm quantile
#' @importFrom utils object.size read.table write.table read.csv
#' @importFrom dplyr summarise group_by distinct
#' @importFrom ggplot2 ggplot geom_line aes geom_tile geom_vline geom_text
#' @importFrom ggplot2 theme theme_bw element_text element_blank element_line facet_grid guide_legend stat_density xlab ylab
#' @importFrom ggplot2 scale_fill_gradientn scale_color_manual scale_linetype_manual scale_size_manual 
#' @importFrom ggplot2 ggplot_build ggtitle geom_histogram scale_fill_manual
#' @importFrom plotly layout
#' @importFrom reshape2 melt 
#' @importFrom magrittr %>%
#' @importFrom testthat expect_equal
#' @importFrom igraph graph_from_data_frame as_adj
#' @importFrom iheatmapr main_heatmap modify_layout add_row_labels add_col_labels add_row_title add_col_title add_col_annotation add_col_dendro
#' @importFrom enrichplot dotplot
#' @importFrom AnnotationDbi select
#' @importFrom clusterProfiler enrichGO
#' @importFrom qlcMatrix corSparse cosSparse
#' @importFrom dynamicTreeCut cutreeHybrid
#' @import org.Hs.eg.db
#' @import org.Mm.eg.db
#' @importFrom amap hcluster Dist
#' @import SeuratObject
#' @importFrom Seurat DimPlot
#' @importFrom graphics barplot
#' @importFrom pheatmap pheatmap
#' @importFrom Matrix rowSums
#' @importFrom SparseM t
#' 
NULL