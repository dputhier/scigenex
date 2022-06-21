#################################################################
##    Define the plot_profile for class ClusterSet
#################################################################

#' @title
#' plot_clust
#' @description
#' Plot the results (heatmap or profiles) contained in a ClusterSet object.
#' @param object A ClusterSet object.
#' @param type The type of diagram ("line or "tile").
#' @param to_log2 Whether data should be transform in logarithm base 2 (+ 1 as a pseudocount).
#' @param average_only Only display the average gene expression profile (if type="line").
#' @param average_line_color The color of the average profile.
#' @param standardizing Whether rows should be divided by standard deviation.
#' @param ceil A value for ceiling (NULL for no ceiling). Ceiling is performed after log transformation, centering and standardization.
#' @param floor A value for flooring (NULL for no flooring). Flooring is performed after log transformation, centering and standardization.
#' @param centering Whether rows should be centered. 
#'
#' @return A ggplot diagram.
#' @export
#'
#' @examples
#' # see online examples
setGeneric("plot_clust",
           
           function(object,
                    type = c("line", "tile"),
                    to_log2 = FALSE,
                    average_only=FALSE,
                    average_line_color="skyblue4",
                    standardizing = FALSE,
                    ceil=1,
                    floor=-1,
                    centering = TRUE) {
             
             standardGeneric("plot_clust")
           })

#' @rdname plot_clust
setMethod(
  "plot_clust",
  signature(object = "ClusterSet"),
  function(object,
           type = c("line", "tile"),
           to_log2 = FALSE,
           average_only=FALSE,
           average_line_color="skyblue4",
           standardizing = FALSE,
           ceil=1,
           floor=-1,
           centering = TRUE) {
    
    # The type of diagram
    type <- type[1]
    
    print_msg("getting matrix", msg_type="DEBUG")
    
    nb <- length(object@size)
    m <- object@data
    if(length(object@cell_order) != 0){
      m <- m[, object@cell_order]
    }
    
    if (to_log2) {
      m <- log2(m + 1)
    }
    
    ## median-centering of row
    if (centering) {
      print_msg("Median-centering rows.", msg_type="DEBUG")
      mean_row <- apply(m, 1, mean)
      m <- sweep(m, MARGIN = 1, STATS = mean_row, FUN = "-")
    }
    
    ## Standardizing row
    if (standardizing) {
      print_msg("Standardizing rows.", msg_type="DEBUG")
      sd_row <- apply(m, 1, sd)
      m <- sweep(m, MARGIN = 1, STATS = sd_row, FUN = "/")
    }
    
    ## Ceiling / flooring
    if(!is.null(ceil)){
      print_msg("Ceiling matrix.", msg_type="DEBUG")
      m[m > ceil] <- ceil
    }
    
    if(!is.null(floor)){
      print_msg("Flooring matrix.", msg_type="DEBUG")
      m[m < floor] <- floor
    }
    
    ## melting
    print_msg("Melting matrix.", msg_type="DEBUG")
    m_melt <- as.data.frame(m)
    m_melt$cluster <- object@gene_patterns
    m_melt$gene <- row.names(object@data)
    
    m_melt <- melt(m_melt,
                   id.vars = c("cluster", "gene"),
                   variable.name = "samples"
    )
    
    print_msg("Storing cell types.", msg_type="DEBUG")
    ## Storing cell types:
    
    if(!is.null(object@cell_types)){
      m_melt$cell_types <- as.character(object@cell_types[as.character(m_melt$samples)])
    }else{
      print_msg("Warning: cell type is undefined.", msg_type="WARNING")
      m_melt$cell_types <- "unknown_cell_type"
    }
    
    ## plotting
    # Note that samples, value, gene, cluster
    # may appear as undefined variable to "R check" command.
    # A workaround is to define them as NULL first...
    samples <- value <- gene <- cluster <- cluster_mean <- NULL
    if (type == "line") {
      print_msg("Preparing diagram.", msg_type="DEBUG")
      p <- ggplot(data = m_melt, aes(
        x = samples,
        y = value
      ))
      
      ## displaying cell types:
      if(length(object@cell_types) > 0 & length(object@cell_colors) > 0){
        print_msg("Adding cell populations to the diagram.", msg_type="DEBUG")
        cell_types <- NULL # Avoid "no visible binding for global variable" inn R check.
        p <- p + geom_vline(aes(xintercept= samples, color=cell_types))
        p <- p + scale_color_manual(values=object@cell_colors,  guide = guide_legend(override.aes = list(size = 5)))
      }
      
      if(! average_only){
        print_msg("Adding gene profile.", msg_type="DEBUG")
        p <- p + geom_line(color = "azure3", aes(group = gene), size=0.1)
      }
      
      print_msg("Adding average profile.", msg_type="DEBUG")
      
      p <- p + geom_line(
        data = m_melt %>%
          group_by(cluster, samples) %>%
          summarise(cluster_mean = mean(value)),
        aes(
          x = samples,
          y = cluster_mean,
          group = cluster
        ),
        color = average_line_color,
        size=0.2
      )
      
      print_msg("Faceting.", msg_type="DEBUG")
      p <- p + facet_grid(cluster ~ .)
      
      print_msg("Theming.", msg_type="DEBUG")
      p <- p + theme_bw()
      p <- p + theme(
        strip.text.y = element_text(angle = 0),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
      )
      
      
    } else if (type == "tile") {
      print_msg("Preparing diagram (tile).", msg_type="DEBUG")
      p <- ggplot(
        data = m_melt,
        aes(
          x = samples,
          y = gene,
          fill = value
        )
      )
      
      print_msg("Preparing color palette.", msg_type="DEBUG")
      col <- unlist(strsplit("#67001f,#b2182b,#d6604d,#f4a582,#fddbc7,#f7f7f7,#d1e5f0,#92c5de,#4393c3,#2166ac,#053061", ","))
      color.ramp <- colorRampPalette(col)(10)
      p <- p + geom_tile()
      p <- p + theme_bw()
      p <- p + scale_fill_gradientn(
        colours = color.ramp,
        name = "Signal"
      )
      
      print_msg("Theming.", msg_type="DEBUG")
      p <- p + theme(
        strip.text.y = element_text(angle = 0),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank()
      )
    }
    
    print_msg("Adding facets.", msg_type="DEBUG")
    p <- p + facet_grid(cluster ~ ., scales = "free_y")
    
    return(p)
  }
)




#################################################################
##    Define the viz_dist
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
#' @param cluster A cluster id to plot. Default is NULL for plotting all cluster.
#' @param use_top_genes A logical to indicate whether to use highly similar genes in the slot top_genes of ClusterSet.
#' @param use_core_cells A logical to indicate whether to use core cells obtained by cell_clust function.
#' @param name A title for the heatmap.
#' @param xlab A title for the x axis.
#' @param ylab A title for the y axis.
#' @param colorbar_name A title for the colorbar.
#' @param show_legend A logical to indicate whether to show colorbar.
#' @param colors A vector of colors.
#' @param row_labels A logical to indicate whether to show row labels.
#' @param col_labels A logical to indicate whether to show col labels.
#' @param label_size A value for label font size.
#' @param line_size An integer for the horizontal white line size.
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
#' res <- DBFMCL(data=m,
#' distance_method="pearson",
#' av_dot_prod_min = 0,
#' inflation = 2,
#' k=25,
#' fdr = 10)
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
                         show_dendro = TRUE,
                         cluster = NULL,
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
                         line_size = 15) {
  
  m <- object@data
  
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
      
      show_dendro = FALSE
      
    } else {
#      if (cell_ordering_method == "hclust"){
        print_msg("Ordering cells based on hierarchical clustering.", msg_type="DEBUG")
        # m_dist <- as.dist(1 - cor(m, method = 'pearson'))
        # m_clust <- hclust(m_dist, method = 'average')
        m <- m[,object@cell_clusters$hclust_res$order]
        
        show_dendro = TRUE
        
        if(length(object@cell_clusters$labels) != 0) {
          object@cell_clusters$labels <- object@cell_clusters$labels[colnames(m)]
          object@cell_clusters$cores <- object@cell_clusters$cores[colnames(m)]
        }
        
#      }
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
  if(!is.null(cluster)){
    gene_cl_int <- names(which(object@gene_patterns == cluster))
    m <- m[gene_cl_int,]
  }
  
  # Reduce m rows to only keep genes from top_genes
  if(use_top_genes) {
    if (nrow(object@top_genes) == 1 &
        ncol(object@top_genes) == 1 &
        is.na(object@top_genes[1,1])) {
      stop(paste0("The slot top_genes of the input ClusterSet object is empty. Be sure to run top_genes() before."))
    }
    
    if(is.null(cluster)){
      genes_top <- unlist(as.data.frame(t(object@top_genes)), use.names = FALSE)
      genes_top <- genes_top[!is.na(genes_top)]
      m <- m[genes_top,]
    } else {
      genes_top <- unlist(as.data.frame(t(object@top_genes[cluster,])), use.names = FALSE)
      genes_top <- genes_top[!is.na(genes_top)]
      m <- m[genes_top,]
    }
  }
  
  # Reduce m cols to only keep core cells from cell_clust
  if(use_core_cells) {
    if(length(object@cell_clusters) == 0){
      stop(paste0("The slot cell_clusters of the input ClusterSet object is empty. Be sure to run cell_clust() before."))
    } else {
      cell_names <- names(which(sort(object@cell_clusters$cores) !=0))
      m <- m[,cell_names]
    }
  } else {
    if(is.null(cell_order)){
      cell_names <- names(object@cell_clusters$labels)
    } else {
      cell_names <- names(sort(object@cell_clusters$labels))
    }
  }
  
  
  # Add blank row to separate feature clusters in heatmap
  if(is.null(cluster)){
    ## Create blank row
    blank_row <- matrix(nrow = line_size, ncol = ncol(m))
    
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
  if(!length(object@cell_clusters) == 0) {
    htmp <- htmp %>% add_col_annotation( data.frame("Clusters" = as.factor(object@cell_clusters$labels[cell_names])), colors = list("Clusters"= c("#DB2020", "#DA7316", "#F0AE00", "#6D9D1E", "#1882C0", "#71529A", "#D02494",
                                                                                                                                                  "#9F1717", "#AE5B11", "#C48D00", "#517416", "#115C8A", "#584178", "#9D1C70")))
  }
  
  # Show dendrogram from hclust
  if(show_dendro & !(use_core_cells) & is.null(cell_order)) {
    htmp <- htmp %>% add_col_dendro(object@cell_clusters$hclust_res, reorder = FALSE)
  }
  
  
  return(htmp)
}
