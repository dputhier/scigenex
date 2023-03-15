#################################################################
##    Define the plot_heatmap
#################################################################

#' @title
#' plot_heatmap
#' @description
#' Plot the observed and simulated distance with the Kth nearest neighbors.
#' @param object A ClusterSet object.
#' @param center A logical to indicate whether to center row.. 
#' @param ceil A value for ceiling (NULL for no ceiling). Ceiling is performed centering.
#' @param floor A value for flooring (NULL for no flooring). Flooring is performed after centering.
#' @param cell_clusters A vector of cell clusters with cell barcodes as names.
#' @param show_dendro A logical to indicate whether to show column dendrogram.
#' @param use_top_genes A logical to indicate whether to use highly similar genes in the slot top_genes of ClusterSet.
#' @param interactive A logical to indicate if the heatmap should be interactive.
#' @param name A title for the heatmap.
#' @param xlab A title for the x axis.
#' @param ylab A title for the y axis.
#' @param colorbar_name A title for the colorbar.
#' @param show_legend A logical to indicate whether to show colorbar.
#' @param colors A vector of colors.
#' @param colors_cell_clusters A named vector of colors for cell identity annotations.
#' @param row_labels A logical to indicate whether to show row labels.
#' @param col_labels A logical to indicate whether to show col labels.
#' @param label_size A value for label font size.
#' @param line_size_vertical An integer for the size of horizontal white line which separate gene clusters.
#' @param line_size_horizontal An integer for the size of vertical white line  which separate cell clusters.
#' @param link The aggloremative criterion for hierarchical clustering. One of "average", "complete" or "single". 
#' Default to average.
#' @return Iheatmap-class object.
#' @export
#'
#' @examples
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
#' plot_heatmap(object = res)
#' plot_heatmap(object = res[1, ])
#' 
#' @rdname plot_heatmap
plot_heatmap <- function(object,
                         center = TRUE,
                         ceil = 1,
                         floor = -1,
                         cell_clusters = NULL,
                         show_dendro = TRUE,
                         use_top_genes = FALSE,
                         interactive = TRUE,
                         name = NULL,
                         xlab = NULL,
                         ylab = NULL,
                         colorbar_name = "Exp. level",
                         show_legend = TRUE,
                         colors = colors_for_gradient("Ju1"),
                         colors_cell_clusters = NULL,
                         row_labels = TRUE,
                         col_labels = FALSE,
                         label_size = 10,
                         line_size_vertical = 3,
                         line_size_horizontal = 3,
                         link=c("average", "complete", "single")) {
  
  link <- match.arg(link) 
  
  check_format_cluster_set(object)
  
  if(show_dendro & is.null(cell_clusters))
    print_msg("cell_clusters is not NULL. Setting show_dendro to FALSE",
              msg_type = "INFO")
  
  if(!is.null(cell_clusters) & is.null(colors_cell_clusters))
    colors_cell_clusters <- discrete_palette(n=length(unique(cell_clusters)), "ggplot")
  
  if(!is.null(cell_clusters)){
    if(is.null(names(cell_clusters)))
      print_msg("The cell_clusters argument should be a named vector", msg_type="STOP")
    print_msg("Extracting cell identity.", msg_type="INFO")
    cell_clusters <- cell_clusters[names(cell_clusters) %in% colnames(object@data)] 
    if(inherits(cell_clusters, "factor"))
       cell_clusters <- droplevels(cell_clusters)
  }
  
  # Ensure there are enough colors
  if(length(as.character(colors_cell_clusters)) < length(table(as.character(cell_clusters))))
    print_msg("Not enough colors for cell_clusters (see colors_cell_clusters).", 
              msg_type="STOP")
  
  # Ensure there is no class equal to zero
  if(0 %in% cell_clusters){
    names_cell_clusters <- names(cell_clusters)
    cell_clusters <- as.numeric(cell_clusters)
    if(0 %in% cell_clusters){
      cell_clusters <- cell_clusters + 1
      print_msg("Found 0 in cell_clusters.", msg_type = "DEBUG")
    }
    names(cell_clusters) <- names_cell_clusters
    cell_clusters <- as.factor(cell_clusters)
  }

  print_msg(paste0("Color palette for cells: ", paste0(colors_cell_clusters, collapse=", ")),
            msg_type = "DEBUG")
  
  m <- as.matrix(object@data)
  
  # Centering
  if(center){
    print_msg("Centering matrix.", msg_type="INFO")
    m <- t(scale(t(m), center = TRUE, scale = FALSE))
  }
  
  # Ceiling and flooring
  if(!is.null(ceil)){
    print_msg("Ceiling matrix.", msg_type="INFO")
    m[m > ceil] <- ceil
  }
  
  if(!is.null(floor)){
    print_msg("Flooring matrix.", msg_type="INFO")
    m[m < floor] <- floor
  }  
  
  gene_to_clust <- gene_cluster(object)
  
  if(is.null(cell_clusters)){
    print_msg("Ordering cells/columns using hierarchical clustering.",
              msg_type = "INFO")
    
    dist_cells <- cor(m, method = "pearson")
    dist_cells <- as.dist((1-dist_cells)/2)
    hclust_cells <- hclust(dist_cells, method = link)
    hclust_cells_order <- colnames(m)[hclust_cells$order]
  }
  
  if(use_top_genes) {
    print_msg("Only top genes will be used.", 
              msg_type = "INFO")
    if (length(object@top_genes) == 0) 
      print_msg("If use_top_genes is TRUE, run top_genes() before.", 
                msg_type = "STOP")
    gene_top <- unlist(object@top_genes, use.names = FALSE)
    m <- m[gene_top, ]
    gene_to_clust <- gene_to_clust[gene_top]
  }
  
  # Add blank row to separate gene clusters in the heatmap

  if(length(table(as.character(gene_to_clust))) > 1){
    
    blank_row <- matrix(NA, nrow = line_size_horizontal, ncol = ncol(m))
    
    print_msg(paste0("line_size_horizontal: ", line_size_horizontal), msg_type = "DEBUG")
    print_msg(paste0("Dim[1] blank matrix: ", nrow(blank_row)), msg_type = "DEBUG")
    print_msg(paste0("Dim[2] blank matrix: ", ncol(blank_row)), msg_type = "DEBUG")
    
    colnames(blank_row) <- colnames(m)
    m_split <- split(as.data.frame(m), gene_to_clust)
    
    nb_NA_row <- 1
    
    for(i in 1:(length(m_split)-1)){
      print_msg(paste0("Adding blank lines for cluster ", i), 
                msg_type = "DEBUG")
      rownames(blank_row) <- paste("NA.", 
                                   nb_NA_row:(nb_NA_row + line_size_horizontal - 1),
                                   sep="")
      nb_NA_row <- nb_NA_row + line_size_horizontal
      m_split[[i]] <- as.matrix(rbind(m_split[[i]], blank_row))
    }
    
    m_blank <-  do.call(rbind, m_split)
    m <- as.matrix(m_blank)
  }

  if(!is.null(cell_clusters) & length(table(as.character(cell_clusters))) > 1){
    
    # Add blank col to separate cell clusters in the heatmap
    blank_col <- matrix(NA, ncol = nrow(m), nrow = line_size_vertical)
    
    print_msg(paste0("line_size_vertical: ", line_size_vertical), msg_type = "DEBUG")
    print_msg(paste0("Dim[1] blank matrix: ", nrow(blank_col)), msg_type = "DEBUG")
    print_msg(paste0("Dim[2] blank matrix: ", ncol(blank_col)), msg_type = "DEBUG")
    
    colnames(blank_col) <- rownames(m)
    m_split <- split(as.data.frame(t(m)), cell_clusters)
    
    cell_clusters <- sort(cell_clusters)
    m <- m[ , names(cell_clusters)]
    nb_NA_row <- 1

    for(i in 1:(length(m_split)-1)){
      print_msg(paste0("Adding blank lines for cell cluster ", i), 
                msg_type = "DEBUG")
      rownames(blank_col) <- paste("NA.", 
                                   nb_NA_row:(nb_NA_row + line_size_vertical - 1),
                                   sep="")
      nb_NA_row <- nb_NA_row + line_size_vertical
      m_split[[i]] <- as.matrix(rbind(m_split[[i]], blank_col))
    }
      
      m_blank <-  do.call(rbind, m_split)
      # Now correct also the names that have been 
      # changed by do.call()...
      row.names(m_blank) <- gsub("^[0-9]+\\.", "", rownames(m_blank), perl=T)
      m <- as.matrix(m_blank)
      # These rownames are in fact colnames (it will be transposed).
      # later on we will need the colnames to be used as rownames 
      # of the annotation_col dataframe. So we need colnames
      # to be unique. We wil add 1, 2 (...) space char... (" ").
      pos_NA <- grep("^NA\\.[0-9]+$", rownames(m), perl=T)
      row.names(m)[pos_NA] <- unlist(lapply(mapply(rep, ' ', 1:length(pos_NA)), paste0, collapse=""))
      m <- t(m)
  }
  
  # Now correct also the rownames that have been 
  # changed by do.call()
  if(length(table(as.character(gene_to_clust))) > 1){
    row.names(m) <- gsub("^[0-9]+\\.", "", rownames(m), perl=T)
    row.names(m)[grep("^NA\\.[0-9]+$", rownames(m), perl=T)] <- " "
  }

  if(is.null(cell_clusters)) {
    # This condition, if(interactive), 
    # may seem weird. In fact interactive 
    # and non interactive do not use the 
    # same plotting function. They have
    # different requirement regarding 
    # column sorting
    if(interactive)
    m <- m[, hclust_cells_order]
  }
  
  # Preparing a data.frame containing cell annotations
  if(!is.null(cell_clusters)){
    column <- as.factor(as.numeric(as.character(cell_clusters[colnames(m)])))
    cell_clusters_anno <- data.frame("Indent."=column)
    rownames(cell_clusters_anno) <- colnames(m)
  }

  ####### Heatmap #######
  # Main heatmap
  print_msg("Plotting heatmap.", msg_type="INFO")
  
  if(interactive){
    
    #Flip rows
    m <- m[order(nrow(m):1),]
    
    print_msg("Plot is interactive...")
    
    htmp <- iheatmapr::main_heatmap(data = m,
                                   name = colorbar_name,
                                   show_colorbar = show_legend,
                                   colors = colors)
    
    htmp <- htmp %>% iheatmapr::modify_layout(list(margin = list(t=20, 
                                                                 r=10, 
                                                                 b=20, 
                                                                 l=20)))
    
    print_msg("Adding labels.", msg_type="DEBUG")
    
    if(row_labels){
      htmp <- htmp %>% iheatmapr::add_row_labels(font = list(size = label_size))}
    if(col_labels){
      htmp <- htmp %>% iheatmapr::add_col_labels(font = list(size = label_size))}
    
    if(!is.null(cell_clusters)){
      
      # Here, there is a bug with add_col_annotation()
      # if a single color is passed then it calls brewer.pal()
      # with the color as the palette name. To fix, we will
      # pass it 3 times the same color wich fix the bug...
      if(length(colors_cell_clusters) == 1)
        colors_cell_clusters <- rep(colors_cell_clusters, 3)

      htmp <- htmp %>% iheatmapr::add_col_annotation(cell_clusters_anno, 
                                                     colors = list("Indent." = colors_cell_clusters))
    }
    
    if(show_dendro & is.null(cell_clusters)) {
      htmp <- htmp %>% iheatmapr::add_col_dendro(hclust_cells, reorder = FALSE)
    }

    print_msg("Adding Titles.", msg_type="DEBUG")
    
    if(!is.null(ylab)){
      htmp <- htmp %>% add_row_title(ylab, side="right", font = list(size = 12))}
    if(!is.null(xlab)){
      htmp <- htmp %>% add_col_title(xlab, side="top", font = list(size = 12))}
    if(!is.null(name)){
      htmp <- htmp %>% add_col_title(name, side="top", font = list(size = 24))}
  }else{
    # Reorder rows to get the same order as the interactive heatmap
    #m <- m[order(nrow(m):1),]
    
    cluster_cols <- F
    annotation_col <- NA
    
    if(is.null(cell_clusters)){
      if(show_dendro){
        cluster_cols <- hclust_cells
      }
    }

    if(!is.null(cell_clusters)){
      annotation_col <- cell_clusters_anno
    }

    htmp <- pheatmap(mat = m, 
                     annotation_legend = show_legend, legend = show_legend,
                     color = colors,
                     cluster_rows = F, 
                     cluster_cols = cluster_cols, 
                     fontsize_row= label_size,
                     show_rownames = row_labels, 
                     show_colnames = col_labels, 
                     annotation_col = annotation_col,
                     border_color = NA, 
                     scale = "none", 
                     annotation_colors=list("Indent."=setNames(colors_cell_clusters, 
                                                                    unique(as.character(sort(cell_clusters))))),
                     na_col = "white")
  }

  return(htmp)
}

#################################################################
##    Define the plot_dist function
#################################################################

#' Plot Distribution of distances
#' 
#' This function creates a histogram of the simulated and observed distances 
#' of a ClusterSet object.
#' 
#' @param object A ClusterSet object.
#' @param bins The number of bins to use in the histogram (default is 150).
#' @param alpha The level of transparency for the bars in the histogram
#'  (default is 0.5).
#' @param colors A named vector of colors to use for the histograms, 
#' with "Simulated" and "Observed" as the named elements 
#' (default is c("Simulated" = "#FB8500", "Observed" = "#36949D")).
#' @param xlim Limits for the x axis of the histogram (e.g. c(0.8, 1)).
#' @param vline_color Color of the vertical line indicating 
#' the critical distance with KNN.
#' @param text_size Font size for the label of the critical distance.
#' @param text_hjust Horizontal justification for the label 
#' of the critical distance.
#' @param text_vjust Vertical justification for the label 
#' of the critical distance.
#' 
#' @return A ggplot object.
#' 
#' @examples
#' 
#' set_verbosity(0)
#' m <- create_4_rnd_clust()
#' 
#' res <- find_gene_clusters(data=m,
#'                           distance_method="pearson",
#'                           inflation = 2,
#'                           k=75,
#'                           row_sum=-Inf,
#'                           highest=0.3,
#'                           min_nb_supporting_cell = 0,
#'                           fdr = 1e-8)
#' p <- plot_dist(res)
#' 
#' @export plot_dist
#' 
plot_dist <- function(object,
                      bins=150,
                      alpha=0.5,
                      colors=c("Simulated"="#FB8500", "Observed"="#36949D"),
                      xlim=NULL,
                      vline_color = "#4f6d7a",
                      text_size = 4,
                      text_hjust = -0.8,
                      text_vjust = -0.5) {
  
  check_format_cluster_set(object)
  
  DKNN = c(object@dbf_output$simulated_dknn,
           object@dbf_output$dknn)
  Type = c(rep("Simulated",
               length(object@dbf_output$simulated_dknn)),
           rep("Observed",
               length(object@dbf_output$dknn)))
  df <- data.frame(DKNN, Type)
  p <-  ggplot(df, aes(x=DKNN, fill=Type)) +
    geom_histogram(bins=bins, 
                   position="identity", 
                   alpha=alpha, 
                   color="white") + 
    theme_bw() + 
    scale_fill_manual(values=colors) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.border = element_blank(),
          legend.title = element_text("Distance with KNN"),
          axis.title = element_text(size = 15),
          axis.text = element_text(size = 10)) +
    xlab(label = "Distance with KNN") +
    ylab(label = "Count")  +
    geom_vline(aes(xintercept = object@dbf_output$critical_distance),
               color = vline_color,
               linetype = "dotdash") +
    geom_text(mapping = aes(x = object@dbf_output$critical_distance,
                            y = 0,
                            label = "Critical distance with KNN",
                            hjust = text_hjust,
                            vjust = text_vjust,
                            angle = 90),
              color = vline_color,
              size = text_size)
  
  if(!is.null(xlim))
    p <- p + xlim(xlim)
  
  return(p)
}
