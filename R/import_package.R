# @importClassesFrom methods character list matrix vector numeric
#' 
#' @useDynLib scigenex
#' 
#' @importMethodsFrom methods show
#' 
#' @importFrom methods new slot slotNames setClass signature setMethod setGeneric
#' @importFrom grDevices colorRampPalette
#' @importFrom stats sd dist hclust heatmap median rnorm setNames as.dist hclust
#' @importFrom utils object.size read.table write.table read.csv
#' @importFrom dplyr summarise group_by distinct
#' @importFrom ggplot2 ggplot geom_line aes geom_tile geom_vline geom_text
#' @importFrom ggplot2 theme theme_bw element_text element_blank element_line facet_grid guide_legend stat_density xlab ylab
#' @importFrom ggplot2 scale_fill_gradientn scale_color_manual scale_linetype_manual scale_size_manual 
#' @importFrom reshape2 melt 
#' @importFrom magrittr %>%
#' @importFrom testthat expect_equal
#' @importFrom igraph graph_from_data_frame as_adj
#' @importFrom MCL mcl
#' @importFrom iheatmapr main_heatmap modify_layout add_row_labels add_col_labels add_row_title add_col_title
#' @importFrom gprofiler2 gost gostplot
#' 
#' 
NULL