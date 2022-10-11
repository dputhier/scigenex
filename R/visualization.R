#################################################################
##    Define the viz_dist function
#################################################################

#' @title
#' viz_dist
#' @description
#' Plot the observed and simulated distance with the Kth nearest neighbors.
#' @param object A ClusterSet object.
#' @param line_type A vector of character defining the line type for the observed distance line and the simulated distance line.
#' @param line_color A vector of numeric defining the color for the observed distance line and the simulated distance line.
#' @param line_size A vector of numeric defining the size for the observed distance line and the simulated distance line.
#' @param vline_type Type of vertical line.
#' @param vline_color Color of vertical line.
#' @param vline_size Size of vertical line
#' @param text_size Size of vertical line text
#' @param text_hjust Horizontal position adjustment of vertical line text.
#' @param text_vjust Vertical position adjustment of vertical line text
#' input file containing distances and cutting threshold is stored. 
#' Default to current working directory.
#'
#' @return A ggplot diagram.
#' @export
#'
#' @examples
#' # see online examples

#' @rdname viz_dist


viz_dist <-  function(object,
                       line_type = c("solid", "longdash"),
                       line_color = c("#006D77", "#83C5BE"),
                       line_size = c(1, 0.8),
                       vline_type = "dotdash",
                       vline_color = "#E29578",
                       vline_size = 0.5,
                       text_size = 4,
                       text_hjust = -0.8,
                       text_vjust = -0.5) {
  
  # Extract observed and simulated distances
  dist_obs <- data.frame(distance_value = object@distances, type = "Observed")
  dist_sim <- data.frame(distance_value = object@simulated_distances, type = "Simulated")
  dist_p <- rbind(dist_obs, dist_sim)
  
  # plot density of distance values for the observed and simulated conditions
  p <- ggplot(data = dist_p, aes(x = distance_value, color = type)) +
    stat_density(aes(linetype = type, size = type), geom = "line", position = "identity") +
    scale_linetype_manual(breaks = c("Observed", "Simulated"), values = line_type) +
    scale_size_manual(values = line_size) +
    scale_color_manual(values = line_color) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.border = element_blank(),
          legend.title = element_text("Distance with KNN"),
          axis.title = element_text(size = 15),
          axis.text = element_text(size = 10)) +
    xlab(label = "Distance with KNN") +
    ylab(label = "Density") +
    geom_vline(aes(xintercept = object@critical_distance),
               color = vline_color,
               linetype = vline_type,
               size = vline_size) +
    geom_text(mapping = aes(x = object@critical_distance,
                            y = 0,
                            label = "Critical distance with KNN",
                            hjust = text_hjust,
                            vjust = text_vjust,
                            angle = 90),
              color = vline_color,
              size = text_size)
  
  return(p)
  
}




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
#' @param cell_order A vector of cell names already ordered.
#' @param cell_ordering_method The clustering method to be used. This must be "hclust", ADD OTHER CLUSTERING METHOD.
#' @param show_dendro A logical to indicate whether to show column dendrogram.
#' @param gene_cluster A cluster id to plot. Default is NULL for plotting all cluster.
#' @param use_top_genes A logical to indicate whether to use highly similar genes in the slot top_genes of ClusterSet.
#' @param use_core_cells A logical to indicate whether to use core cells obtained by find_cell_clusters function.
#' @param name A title for the heatmap.
#' @param xlab A title for the x axis.
#' @param ylab A title for the y axis.
#' @param colorbar_name A title for the colorbar.
#' @param show_legend A logical to indicate whether to show colorbar.
#' @param colors A vector of colors.
#' @param row_labels A logical to indicate whether to show row labels.
#' @param col_labels A logical to indicate whether to show col labels.
#' @param label_size A value for label font size.
#' @param line_size_vertical An integer for the size of horizontal white line which separate gene clusters.
#' @param line_size_horizontal An integer for the size of vertical white line  which separate cell clusters.
#'
#' @return Iheatmap-class object
#' @export
#'
#' @examples
#' m <- matrix(rnorm(40000), nc=20)
#' m[1:100,1:10] <- m[1:100,1:10] + 4
#' m[101:200,11:20] <- m[101:200,11:20] + 3
#' m[201:300,5:15] <- m[201:300,5:15] + -2
#' 
#' res <- find_gene_clusters(data=m,
#'                           distance_method="pearson",
#'                           av_dot_prod_min = 0,
#'                           inflation = 2,
#'                           k = 25,
#'                           fdr = 10)
#' 
#' plot_heatmap(object = res)
#' plot_heatmap(object = res, cluster = "1")
#' 

#' @rdname plot_heatmap

plot_heatmap <- function(object,
                         center = TRUE,
                         ceil = 1,
                         floor = -1,
                         cell_order = NULL,
                         #cell_ordering_method = "hclust",
                         show_dendro = FALSE,
                         gene_cluster = NULL,
                         use_top_genes = FALSE,
                         use_core_cells = FALSE,
                         name = NULL,
                         xlab = NULL,
                         ylab = NULL,
                         colorbar_name = "Expression level",
                         show_legend = TRUE,
                         colors = c("#A9D6E5", "#2166AC", "#000000", "#B2182B", "#FFCA3A"),
                         row_labels = TRUE,
                         col_labels = FALSE,
                         label_size = 9,
                         line_size_vertical = 15,
                         line_size_horizontal = 15) {
  
  m <- object@data[names(object@gene_clusters),]
  
  # Config
  if(use_core_cells | use_top_genes) {
    show_dendro = FALSE
  }
  
  # Centering
  if(center){
    print_msg("Centering matrix.", msg_type="DEBUG")
    m <- t(scale(t(m), center = TRUE, scale = FALSE))
  }
  
  
  # Cell order
  if (!is.null(cell_order)){
    print_msg("Ordering cells.", msg_type="DEBUG")
    m <- m[,cell_order]
  } else {
    if (length(object@cell_clusters) == 0) {
      object@cell_clusters$labels <- object@cell_clusters$labels[colnames(m)]
      object@cell_clusters$cores <- object@cell_clusters$cores[colnames(m)]
      
    } else {
      
      if(!show_dendro) {
        m <- m[,names(sort(object@cell_clusters$labels))]
        object@cell_clusters$labels <- object@cell_clusters$labels[colnames(m)]
        object@cell_clusters$cores <- object@cell_clusters$cores[colnames(m)]
      } else {
        print_msg("Ordering cells based on hierarchical clustering.", msg_type="DEBUG")
        m <- m[,object@cell_clusters$hclust_res$order]
      }
    }
  }
  
  
  # Ceiling and flooring
  if(!is.null(ceil)){
    print_msg("Ceiling matrix.", msg_type="DEBUG")
    m[m > ceil] <- ceil
  }
  
  if(!is.null(floor)){
    print_msg("Flooring matrix.", msg_type="DEBUG")
    m[m < floor] <- floor
  }  
  
  
  # Reduce matrix to one cluster
  if(!is.null(gene_cluster)){
    gene_cl_int <- names(which(object@gene_clusters == gene_cluster))
    m <- m[gene_cl_int,]
  }
  
  # Reduce m rows to only keep genes from top_genes
  if(use_top_genes) {
    if (nrow(object@top_genes) == 1 &
        ncol(object@top_genes) == 1 &
        is.na(object@top_genes[1,1])) {
      stop(paste0("The slot top_genes of the input ClusterSet object is empty. Be sure to run top_genes() before."))
    }
    
    if(is.null(gene_cluster)){
      genes_top <- unlist(as.data.frame(t(object@top_genes)), use.names = FALSE)
      genes_top <- genes_top[!is.na(genes_top)]
      m <- m[genes_top,]
    } else {
      genes_top <- unlist(as.data.frame(t(object@top_genes[gene_cluster,])), use.names = FALSE)
      genes_top <- genes_top[!is.na(genes_top)]
      m <- m[genes_top,]
    }
  }
  
  # Reduce m cols to only keep core cells from find_cell_clusters function
  if(use_core_cells) {
    if(length(object@cell_clusters) == 0){
      stop(paste0("The slot cell_clusters of the input ClusterSet object is empty. Be sure to run find_cell_clusters() before."))
    } else {
      cell_names <- names(which(sort(object@cell_clusters$cores) !=0))
      m <- m[,cell_names]
    }
  } else {
    if(is.null(cell_order)){
      cell_names <- object@cell_clusters$labels
    } else {
      cell_names <- sort(object@cell_clusters$labels)
    }
  }
  
  
  # Add blank row to separate gene clusters in heatmap
  if(is.null(gene_cluster)){
    ## Create blank row
    blank_row <- matrix(nrow = line_size_horizontal, ncol = ncol(m))
    
    ## Insert blank row in matrix
    m_blank <- matrix(ncol = ncol(m))
    m <- rbind(m, matrix(nrow = 1, ncol = ncol(m)))
    
    for (i in 1:length(object@size)) {
      if(!use_top_genes) {
        row_start <- sum(object@size[0:(i-1)])+1
        row_end <- sum(object@size[1:i])
        
        m_blank_loop <- rbind(m[row_start:row_end,], blank_row)
      } else {
        nb_top_genes <- length(object@top_genes[i,!is.na(object@top_genes[i,])])
        
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
  if(!length(object@cell_clusters) == 0 & is.null(cell_order) & !show_dendro){
    ## Create blank row
    blank_col <- matrix(nrow = nrow(m), ncol = line_size_vertical)
    
    ## Insert blank row in matrix
    m_blank <- matrix(nrow = nrow(m))
    m <- cbind(m, matrix(nrow = nrow(m), ncol = 1))
    
    cell_names_blank <- c(NA)
    #cell_names_blank <- cell_names
    
    for (i in sort(unique(object@cell_clusters$labels))) {
      if(!use_core_cells) {
        col_start <- sum(table(object@cell_clusters$labels)[0:(i-1)])+1
        col_end <- sum(table(object@cell_clusters$labels)[1:i])
        
        m_blank_loop <- cbind(m[,col_start:col_end], blank_col)
        cell_names_loop <- c(cell_names[col_start:col_end], rep(NA, line_size_vertical))
      } else {
        nb_core_cells <- length(which(object@cell_clusters$cores == i))
        
        m_core_i <- m[,0:nb_core_cells]
        m <- m[,(nb_core_cells+1):ncol(m)]
        
        m_blank_loop <- cbind(m_core_i, blank_col)
        
        cell_names_loop <- c(object@cell_clusters$cores[which(object@cell_clusters$cores == i)], rep(NA, line_size_vertical) )
      }
      
      
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
  print_msg("Plotting heatmap.", msg_type="DEBUG")
  
  htmp <- main_heatmap(data = m,
                       name = colorbar_name,
                       show_colorbar = show_legend,
                       colors = colors,
                       row_order = seq_len(nrow(m)))
  
  htmp <- htmp %>% modify_layout(list(margin = list(t=20, r=10, b=20, l=10)))
  
  # Add labels
  if(row_labels){
    htmp <- htmp %>% add_row_labels(font = list(size = label_size))}
  if(col_labels){
    htmp <- htmp %>% add_col_labels(font = list(size = label_size))}
  
  
  # Add titles
  if(!is.null(ylab)){
    htmp <- htmp %>% add_row_title(ylab, side="right", font = list(size = 12))}
  if(!is.null(xlab)){
    htmp <- htmp %>% add_col_title(xlab, side="top", font = list(size = 12))}
  if(!is.null(name)){
    htmp <- htmp %>% add_col_title(name, side="top", font = list(size = 24))}
  
  # Show cell clusters
  if(!length(object@cell_clusters) == 0 & is.null(cell_order) & !show_dendro) {
    htmp <- htmp %>% add_col_annotation( data.frame("Clusters" = as.factor(cell_names_blank)),
                                         colors = list("Clusters"= c("#9F1717", "#AE5B11", "#C48D00", "#517416", "#115C8A", "#584178", "#9D1C70",
                                                                     "#E96767", "#EC9342", "#FFCA3A", "#8AC926", "#4DADE8", "#9579B9", "#E25CB4", 
                                                                     "#DB2020", "#DA7316", "#F0AE00", "#6D9D1E", "#1882C0", "#71529A", "#D02494",
                                                                     "#EF9292", "#F2B57D", "#FFDA77", "#B6E36A", "#7BC4EE", "#AD98C9", "#EA8AC9")))
  } else {
    if (!length(object@cell_clusters) == 0 & !show_dendro) {
      htmp <- htmp %>% add_col_annotation( data.frame("Clusters" = as.factor(object@cell_clusters$labels[names(cell_names)])),
                                           colors = list("Clusters"= c("#9F1717", "#AE5B11", "#C48D00", "#517416", "#115C8A", "#584178", "#9D1C70",
                                                                       "#E96767", "#EC9342", "#FFCA3A", "#8AC926", "#4DADE8", "#9579B9", "#E25CB4", 
                                                                       "#DB2020", "#DA7316", "#F0AE00", "#6D9D1E", "#1882C0", "#71529A", "#D02494",
                                                                       "#EF9292", "#F2B57D", "#FFDA77", "#B6E36A", "#7BC4EE", "#AD98C9", "#EA8AC9")))
    } else {
      if(!length(object@cell_clusters) == 0) {
        htmp <- htmp %>% add_col_annotation( data.frame("Clusters" = as.factor(object@cell_clusters$labels[object@cell_clusters$hclust_res$order])),
                                             colors = list("Clusters"= c("#9F1717", "#AE5B11", "#C48D00", "#517416", "#115C8A", "#584178", "#9D1C70",
                                                                         "#E96767", "#EC9342", "#FFCA3A", "#8AC926", "#4DADE8", "#9579B9", "#E25CB4", 
                                                                         "#DB2020", "#DA7316", "#F0AE00", "#6D9D1E", "#1882C0", "#71529A", "#D02494",
                                                                         "#EF9292", "#F2B57D", "#FFDA77", "#B6E36A", "#7BC4EE", "#AD98C9", "#EA8AC9")))
      }
    }
  }
  
  # Show dendrogram from hclust
  if(show_dendro & !(use_core_cells) & is.null(cell_order) & !length(object@cell_clusters) == 0) {
    htmp <- htmp %>% add_col_dendro(object@cell_clusters$hclust_res, reorder = FALSE)
  }
  
  
  return(htmp)
}
