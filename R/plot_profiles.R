#' Plot mean expression profiles of each cluster from a ClusterSet object.
#'
#' This function generates a barplot showing the expression profiles of
#' cell type-specific genes across different cell types.
#'
#' @param data A ClusterSet object.
#' @param ident A named vector containing the cell type identities for each cell.
#' Typically the result from the Idents() function on a Seurat object (see Seurat library).
#' @param nb_column The number of columns in the facet grid of the plot. If not provided,
#' it is automatically computed as the square root of the number of cell types.
#' @param color_cell_type A named vector of colors (with cell type as names) used to 
#' distinguish between different cell types in the plot. If not provided, the default 
#' hue color palette is used.
#' @param size_text_y The font size of the y-axis tick labels.
#' @param size_label The font size of the cluster labels.
#' @param legend_name A name for the legend.
#' @return A ggplot object showing the expression profiles of cell type-specific genes.
#'
#' @examples
#' # Load a Seurat object
#' data(pbmc_small, package = "SeuratObject")
#' library(Seurat)
#' # Compute the signatures using find_gene_clusters()
#' clust_set <- find_gene_clusters(pbmc_small, k=50, no_dknn_filter=TRUE)
#' plot_profiles(clust_set, ident=Seurat::Idents(pbmc_small))
#' pal <- c("#4E79A7", "#A0CBE8", "#F28E2B")
#' names(pal) <- levels(Idents(pbmc_small))
#' plot_profiles(clust_set, ident=Seurat::Idents(pbmc_small), color_cell_type = pal)
#' plot_profiles(clust_set[2:4,], ident=Seurat::Idents(pbmc_small), color_cell_type = pal)
#' @importFrom ggplot2 geom_col facet_wrap theme_minimal geom_text scale_color_manual
#' @importFrom scales hue_pal
#' @importFrom reshape2 melt
#' @importFrom ggplot2 .data
#' @export
plot_profiles <- function(data = NULL,
                          ident = NULL,
                          nb_column = NULL,
                          color_cell_type = NULL,
                          size_text_y = 5,
                          size_label = 2,
                          legend_name="Cell\ntype") {
  
  if (is.null(data) | !inherits(data, "ClusterSet"))
    print_msg("Please provide a ClusterSet objet.", msg_type = "STOP")
  
  ident <- ident[colnames(data@data)]
  
  centers <- data@dbf_output$center
  
  if (is.null(ident))
    print_msg("Please provide cell identification.", msg_type = "STOP")
  
  ident <- sort(ident)
  ident <- factor(ident, levels=levels(ident), ordered = T)
  nb_cell_type <- length(unique(ident))
  
  
  
  if (is.null(nb_column))
    nb_column <- round(sqrt(nrow(data@dbf_output$center)), 0)
  
  if (is.null(color_cell_type)){
    color_cell_type <- scales::hue_pal()(nb_cell_type)
  }else{
    if(nb_cell_type != length(color_cell_type))
      print_msg("The number of colors should be the same as the number of cell types.", 
                msg_type = "STOP")
    
    if(is.null(names(color_cell_type)))
      print_msg("The color_cell_type argument should be a named vector.", 
                msg_type = "STOP")
    
    if(!all(names(color_cell_type) %in% unique(ident)))
      print_msg("The color_cell_type argument contains unknown cell type.", 
                msg_type = "STOP")
  }
    
  print_msg(paste0("Number of cells types: ", nb_cell_type),
            msg_type = "INFO")
  
  nb_cells <- ncol(centers)
  
  print_msg(paste0("Number of cells: ", nb_cells),
            msg_type = "INFO")
  
  centers <- centers[, names(ident), drop=FALSE]
  
  m <- reshape2::melt(as.matrix(centers))
  
  colnames(m) <- c("Cluster", "Cell", "Intensity")
  m$Cluster <- factor(
    paste0("Cluster: ", m$Cluster),
    levels = paste0("Cluster: ", unique(m$Cluster)),
    ordered = T
  )
  
  m$Ident <- ident[m$Cell]
  
  print_msg(paste0("Centers dimension: ", paste0(dim(centers), collapse = " ")),
            msg_type = "DEBUG")
  
  y_text <- apply(centers, 1, max)
  
  y_text <- y_text + 0.1 * y_text
  df_text <- data.frame(x = colnames(centers)[round(nb_cells / 3, 0)], 
               y = y_text)
  df_text$Cluster <- factor(paste0("Cluster: ", 
                                   rownames(centers)),
                            ordered = T)
  
  ggplot2::ggplot(data= m,
                  ggplot2::aes(
                    x = .data[["Cell"]],
                    y = .data[["Intensity"]],
                    group = 1,
                    fill = .data[["Ident"]]
                  )) + 
    ggplot2::geom_col()  +
    ggplot2::facet_wrap(~Cluster, scales = "free_y",
               ncol = nb_column) +
    ggplot2::theme_minimal() +
    ggplot2::geom_text(
      data = df_text,
      mapping = ggplot2::aes(x = .data[["x"]], 
                             y =.data[["y"]], 
                             label = .data[["Cluster"]]),
      size = size_label,
      inherit.aes = F
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      plot.background = ggplot2::element_rect(fill = "white"),
      axis.text.y = ggplot2::element_text(size = size_text_y),
      strip.text = ggplot2::element_blank()
    ) +
    ggplot2::scale_fill_manual(values=color_cell_type, name=legend_name)
  
}
