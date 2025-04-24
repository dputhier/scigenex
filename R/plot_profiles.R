#' @title Plot mean expression profiles of each cluster from a ClusterSet object.
#'
#' @description This function generates a barplot showing the expression profiles of
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
#' @param to_lin Transform in linear scale (i.e. 2^x).
#' @return A ggplot object showing the expression profiles of cell type-specific genes.
#'
#' @examples
#' # Load a Seurat object
#' load_example_dataset("7870305/files/lymph_node_tiny_2")
#' # Load a ClusterSet object
#' load_example_dataset("7870305/files/lymph_node_tiny_clusters_2")
#'                             
#' plot_profiles(lymph_node_tiny_clusters_2, ident=Seurat::Idents(lymph_node_tiny_2))
#' lv <- levels(Seurat::Idents(lymph_node_tiny_2))
#' pal <- discrete_palette(n=length(lv))
#' names(pal) <- lv
#' plot_profiles(lymph_node_tiny_clusters_2, ident=Seurat::Idents(lymph_node_tiny_2), color_cell_type = pal)
#' plot_profiles(lymph_node_tiny_clusters_2[2:4,], ident=Seurat::Idents(lymph_node_tiny_2), color_cell_type = pal)
#' 
#' @importFrom ggplot2 geom_col facet_wrap theme_minimal geom_text scale_color_manual
#' @importFrom scales hue_pal
#' @importFrom reshape2 melt
#' @importFrom ggplot2 .data
#' @export plot_profiles
plot_profiles <- function(data = NULL,
                          ident = NULL,
                          nb_column = NULL,
                          color_cell_type = NULL,
                          size_text_y = 5,
                          size_label = 2,
                          legend_name="Cell\ntype",
                          to_lin=FALSE) {
  
  if (is.null(data) | !inherits(data, "ClusterSet"))
    print_msg("Please provide a ClusterSet objet.", msg_type = "STOP")

  centers <- data@dbf_output$center
  
  if(is.null(centers))
    print_msg("Please run compute_centers() on ClusterSet.", msg_type = "STOP")
  
  if (is.null(ident))
    print_msg("Please provide cell identification.", msg_type = "STOP")
  
  if (length(ident) < ncol(data@data))
    print_msg("The length of the 'ident' argument should be equal or greater to ncol(data).", msg_type = "STOP")
  
  if (is.null(names(ident)))
    names(ident) <- colnames(data@data)
    
  ident <- ident[colnames(data@data)]
  
  ident <- sort(ident)
  ident <- factor(ident, levels=levels(as.factor(ident)), ordered = T)
  nb_cell_type <- length(unique(ident))
  
  
  if (is.null(nb_column))
    nb_column <- round(sqrt(nrow(centers)), 0)
  
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
  
  if(to_lin)
    m$Intensity <- 2^m$Intensity 
  
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
      inherit.aes = FALSE
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

#' @title Plot a summarized view of expression profiles.
#'
#' @description This function generates line plots displaying mean expression profiles accross
#' cell types and modules.
#'
#' @param data A ClusterSet object.
#' @param ident A named vector containing the cell type identities for each cell.
#' Typically the result from the Idents() function on a Seurat object (see Seurat library).
#' @param color_cell_type A named vector of colors (with cell type as names) used to 
#' distinguish between different clusters.
#' @param size_text_y The font size of the y-axis tick labels.
#' @param size_label The font size of the cluster labels.
#' @param legend_name A name for the legend.
#' @param center Whether to center (substract mean) each row.
#' @return A ggplot object showing the expression profiles of cell type-specific genes.
#'
#' @examples
#' # Load a Seurat object
#' load_example_dataset("7871581/files/pbmc3k_medium")
#' # Load a ClusterSet object
#' load_example_dataset("7871581/files/pbmc3k_medium_clusters")
#'                             
#' plot_multi_profiles(pbmc3k_medium_clusters, ident=Seurat::Idents(pbmc3k_medium))
#' pal <- discrete_palette(nclust(pbmc3k_medium_clusters))
#' names(pal) <- names(pbmc3k_medium_clusters@gene_clusters)
#' plot_multi_profiles(pbmc3k_medium_clusters, 
#'         ident=Seurat::Idents(pbmc3k_medium), color_cluster = pal)
#' plot_multi_profiles(pbmc3k_medium_clusters[2:4,], 
#'       ident=Seurat::Idents(pbmc3k_medium), 
#'       color_cluster = pal[2:4],
#'       center=TRUE)
#' 
#' @importFrom ggplot2 geom_col facet_wrap theme_minimal geom_text scale_color_manual
#' @importFrom scales hue_pal
#' @importFrom reshape2 melt
#' @importFrom ggplot2 .data
#' @export plot_multi_profiles
plot_multi_profiles <- function(data = NULL,
                                ident = NULL,
                                color_cluster = NULL,
                                size_text_y = 5,
                                size_label = 2,
                                legend_name="Cell\ntype",
                                nb_column=NULL,
                                center=FALSE) {
  

  if (is.null(data) | !inherits(data, "ClusterSet"))
    print_msg("Please provide a ClusterSet objet.", msg_type = "STOP")

  centers <- data@dbf_output$center
  
  if(is.null(centers))
    print_msg("Please run compute_centers() on ClusterSet.", msg_type = "STOP")
  
  if (is.null(nb_column))
    nb_column <- round(sqrt(nrow(data@dbf_output$center)), 0)
  
  if(center){
    centers <- sweep(centers, 1, STATS=rowMeans(centers), FUN="-")
  }

  if (is.null(ident))
    print_msg("Please provide cell identification.", msg_type = "STOP")
  
  if (length(ident) < ncol(data@data))
    print_msg("The length of the 'ident' argument should be equal or greater to ncol(data).", msg_type = "STOP")
  
  if (is.null(names(ident)))
    names(ident) <- colnames(data@data)
  
  ident <- ident[colnames(data@data)]
  
  ident <- sort(ident)
  ident <- factor(ident, levels=levels(as.factor(ident)), ordered = T)
  
  if (is.null(color_cluster)){
    color_cluster <- scales::hue_pal()(nrow(centers))
  }else{
    if(nrow(centers) != length(color_cluster)){
      print_msg(paste0("Number of clusters: ", nrow(centers)), msg_type = "DEBUG")
      print_msg(paste0("Number of colors: ", length(color_cluster)), msg_type = "DEBUG")
      print_msg("The number of colors should be the same as the number of gene clusters.", 
                msg_type = "STOP")
    }

    if(is.null(names(color_cluster)))
      print_msg("The color_cluster argument should be a named vector.", 
                msg_type = "STOP")
    
    if(!all(names(color_cluster) %in% unique(rownames(centers))))
      print_msg("The color_cluster argument contains unknown clusters.", 
                msg_type = "STOP")
  }
  
  nb_cells <- ncol(centers)
  
  print_msg(paste0("Number of cells: ", nb_cells),
            msg_type = "INFO")
  
  centers <- centers[, names(ident), drop=FALSE]
  
  centers <- centers[, names(ident), drop=FALSE]
  
  centers_summarized_by_cell_type <- matrix(NA,
                                            nr=nrow(centers),
                                            ncol=length(levels(ident)))
  for(rown in 1:nrow(centers)){
    centers_summarized_by_cell_type[rown, ] <- tapply(centers[rown,], ident, mean)
  }
  
  colnames(centers_summarized_by_cell_type) <- as.character(levels(ident))
  rownames(centers_summarized_by_cell_type) <- rownames(centers)
  
  print_msg(paste0("Centers dimension: ", paste0(dim(centers), collapse = " ")),
            msg_type = "DEBUG")
  
  m <- reshape2::melt(centers_summarized_by_cell_type)

  colnames(m) <- c("Cluster", "Cell_type", "Intensity")
  m$Cell_type <- factor(m$Cell_type, levels=levels(ident), ordered = TRUE)
  m$Cluster <-   factor(m$Cluster, levels=unique(m$Cluster), ordered = TRUE)
  
  ggplot2::ggplot(data= m,
                  ggplot2::aes(
                    x = .data[["Cell_type"]],
                    y = .data[["Intensity"]],
                    group = Cluster,
                    color = .data[["Cluster"]]
                  )) + 
    ggplot2::geom_line() +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      plot.background = ggplot2::element_rect(fill = "white"),
      axis.text.y = ggplot2::element_text(size = size_text_y),
      strip.text = ggplot2::element_blank()
    ) +
    ggplot2::scale_color_manual(values=color_cluster, name=legend_name)
  
}
