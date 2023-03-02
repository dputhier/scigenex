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
#' @param gene_clusters A cluster id to plot. Default is NULL for plotting all cluster.
#' @param use_top_genes A logical to indicate whether to use highly similar genes in the slot top_genes of ClusterSet.
#' @param interactive A logical to indicate if the heatmap should be interactive.
#' @param name A title for the heatmap.
#' @param xlab A title for the x axis.
#' @param ylab A title for the y axis.
#' @param colorbar_name A title for the colorbar.
#' @param show_legend A logical to indicate whether to show colorbar.
#' @param colors A vector of colors.
#' @param colors_cell_clusters A vector of colors for column annotations.
#' @param row_labels A logical to indicate whether to show row labels.
#' @param col_labels A logical to indicate whether to show col labels.
#' @param label_size A value for label font size.
#' @param line_size_vertical An integer for the size of horizontal white line which separate gene clusters.
#' @param line_size_horizontal An integer for the size of vertical white line  which separate cell clusters.
#'
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
#' plot_heatmap(object = res, cluster = "1")
#' 

#' @rdname plot_heatmap

plot_heatmap <- function(object,
                         center = TRUE,
                         ceil = 1,
                         floor = -1,
                         cell_clusters = NULL,
                         show_dendro = TRUE,
                         gene_clusters = NULL,
                         use_top_genes = FALSE,
                         interactive = TRUE,
                         name = NULL,
                         xlab = NULL,
                         ylab = NULL,
                         colorbar_name = "Expression level",
                         show_legend = TRUE,
                         colors = c("#A9D6E5", "#2166AC", "#000000", "#B2182B", "#FFCA3A"),
                         colors_cell_clusters = c("#9F1717", "#AE5B11", "#C48D00", "#517416", "#115C8A", "#584178", "#9D1C70",
                                                  "#E96767", "#EC9342", "#FFCA3A", "#8AC926", "#4DADE8", "#9579B9", "#E25CB4", 
                                                  "#DB2020", "#DA7316", "#F0AE00", "#6D9D1E", "#1882C0", "#71529A", "#D02494",
                                                  "#EF9292", "#F2B57D", "#FFDA77", "#B6E36A", "#7BC4EE", "#AD98C9", "#EA8AC9"),
                         row_labels = TRUE,
                         col_labels = FALSE,
                         label_size = 9,
                         line_size_vertical = 15,
                         line_size_horizontal = 15) {
  
  m <- object@data
  
  # # Config
  if (is.null(gene_clusters)) {
    gene_clusters <- object@gene_clusters_metadata$cluster_id
  }
  
  # Cell order
  if (!is.null(cell_clusters)){
    print_msg("Ordering cells.", msg_type="INFO")
    m <- m[,names(sort(cell_clusters))]
  } else {
    print_msg("Ordering cells/columns using hierarchical clustering.",
              msg_type = "INFO")
    dist_cells <- cor(m, method = "pearson")
    dist_cells <- as.dist((1-dist_cells)/2)
    hclust_cells <- hclust(dist_cells, method = "complete")
    if(interactive){
      # Reorder cells based on hierarchical clustering
      m <- m[,hclust_cells$order]  
    }
  }
  
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
  
  
  # Reduce matrix to gene clusters in gene_clusters parameter
  if(!is.null(gene_clusters)){
    gene_cl_int <- unlist(object@gene_clusters[gene_clusters], use.names = FALSE)
    m <- m[gene_cl_int,]
  }
  
  # Reduce m rows to only keep genes from top_genes
  if(use_top_genes) {
    if (length(object@top_genes) == 0) {
      print_msg("The slot top_genes of the input ClusterSet object is empty. Be sure to run top_genes() before.", msg_type = "STOP")
    }
    
    if(is.null(gene_clusters)){
      genes_top <- unlist(object@top_genes, use.names = FALSE)
      m <- m[genes_top,]
    } else {
      genes_top <- unlist(object@top_genes[gene_clusters], use.names = FALSE)
      m <- m[genes_top,]
    }
  }
  
  
  # Add blank row to separate gene clusters in heatmap
  if(!length(gene_clusters) == 0){
    ## Create blank row
    blank_row <- matrix(nrow = line_size_horizontal, ncol = ncol(m))
    
    ## Insert blank row in matrix
    m_blank <- matrix(ncol = ncol(m))
    m <- rbind(m, matrix(nrow = 1, ncol = ncol(m)))
    
    for (i in gene_clusters) {
      if(!use_top_genes) {
        size <- object@gene_clusters_metadata$size[gene_clusters]
        row_start <- sum(size[names(size) %in% 0:(i-1)])+1
        row_end <- sum(size[names(size) %in% 1:i])
        
        m_blank_loop <- rbind(m[row_start:row_end,], blank_row)
      } else {
        nb_top_genes <- length(object@top_genes[[i]])
        
        m_top_i <- m[0:nb_top_genes,]
        m <- m[(nb_top_genes+1):nrow(m),]
        
        m_blank_loop <- rbind(m_top_i, blank_row)
        #row_start <- i*ncol(object@top_genes) - ncol(object@top_genes) + 1
        #row_end <- i*nb_top_genes
      }
      
      
      #rownames(test)[(nrow(test)-line_size+1):nrow(test)] <- paste(rep(" ", 2), collapse = '')
      m_blank <- rbind(m_blank, m_blank_loop)
    }
    m <- m_blank
  }
  
  
  
  # Add blank col to separate cell clusters in heatmap
  if(!is.null(cell_clusters)){
    ## Create blank row
    blank_col <- matrix(nrow = nrow(m), ncol = line_size_vertical)
    
    ## Insert blank row in matrix
    m_blank <- matrix(nrow = nrow(m))
    m <- cbind(m, matrix(nrow = nrow(m), ncol = 1))
    
    cell_names_blank <- c(NA)
    
    cell_names <- names(sort(cell_clusters))
    
    if (0 %in% cell_clusters){
      cell_clusters_tmp <- as.numeric(as.character(cell_clusters)) + 1 #Add one to each element in the vector
    } else {
      cell_clusters_tmp <- cell_clusters
    }
    
    for (i in sort(unique(cell_clusters_tmp))) {
      col_start <- sum(table(cell_clusters_tmp)[0:(i-1)])+1
      col_end <- sum(table(cell_clusters_tmp)[1:i])
      
      m_blank_loop <- cbind(m[,col_start:col_end], blank_col)
      cell_names_loop <- c(cell_names[col_start:col_end], rep(NA, line_size_vertical))
      
      #rownames(test)[(nrow(test)-line_size+1):nrow(test)] <- paste(rep(" ", 2), collapse = '')
      m_blank <- cbind(m_blank, m_blank_loop)
      cell_names_blank <- c(cell_names_blank, cell_names_loop)
    }
    m <- m_blank
  }
  
  #Flip rows
  m <- m[order(nrow(m):1),]
  
  
  ####### Heatmap #######
  # Main heatmap
  print_msg("Plotting heatmap.", msg_type="INFO")
  
  if(interactive){
    htmp <- main_heatmap(data = m,
                         name = colorbar_name,
                         show_colorbar = show_legend,
                         colors = colors,
                         row_order = seq_len(nrow(m)))
    
    htmp <- htmp %>% modify_layout(list(margin = list(t=20, r=10, b=20, l=10)))
    
    print_msg("Adding labels.", msg_type="DEBUG")
    
    if(row_labels){
      htmp <- htmp %>% add_row_labels(font = list(size = label_size))}
    if(col_labels){
      htmp <- htmp %>% add_col_labels(font = list(size = label_size))}
    
    if(!is.null(cell_clusters)){
      cell_clusters_anno <- cell_clusters[match(cell_names_blank, names(cell_clusters))]
      cell_clusters_anno <- as.factor(cell_clusters_anno)
      cell_clusters_anno <- as.data.frame(cell_clusters_anno)
      colnames(cell_clusters_anno) <- "cell_clusters"
      htmp <- htmp %>% add_col_annotation(cell_clusters_anno, colors = list("cell_clusters" = colors_cell_clusters))
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
    m <- m[order(nrow(m):1),]
    
    if(show_dendro){
      htmp <- pheatmap(mat = m, 
                       annotation_legend = show_legend, legend = show_legend,
                       color = colorRampPalette(colors)(50),
                       cluster_rows = F, cluster_cols = hclust_cells, 
                       show_rownames = row_labels, show_colnames = col_labels, 
                       border_color = NA, scale = "none", na_col = "white")
    }else{
      htmp <- pheatmap(mat = m[, hclust_cells$order], 
                       annotation_legend = show_legend, legend = show_legend,
                       color = colorRampPalette(colors)(50),
                       cluster_rows = F, cluster_cols = F, 
                       show_rownames = row_labels, show_colnames = col_labels, 
                       border_color = NA, scale = "none", na_col = "white")
    }
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
#' plot_dist(res)
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
  
  if (!inherits(object, "ClusterSet")) {
    stop("Please provide ClusterSet object.")
  }
  
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
               linetype = "dotdash",
               size = 0.5) +
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
