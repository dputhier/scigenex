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
    m_melt$cluster <- object@cluster
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
##    Define the plot_distance
#################################################################

#' @title
#' plot_dist
#' @description
#' Plot the observed and simulated distance with the Kth nearest neighbors.
#' @param object A ClusterSet object.
#' @param path A character string representing the data directory where
#' input file containing distances and cutting threshold is stored. 
#' Default to current working directory.
#'
#' @return A ggplot diagram.
#' @export
#'
#' @examples
#' # see online examples

#' @rdname plot_dist


plot_dist <-  function(object,
                       path = ".") {
  #Extract name of the file containing distance and cutting threshold.
  file_name <- object@opt_name
  
  #Path to the file containg distance and cutting threshold.
  if (path == ".") path <- getwd()
  file_path <- file.path(path, file_name)
  file_path <- gsub(pattern = "//", replacement = "/", x = file_path)
  
  #Read file containing distance and cutting threshold.
  opt_data <- readLines(file_path)
  
  # Extract cutting threshold value
  dknn <- opt_data[(which(opt_data == ">>thresholds")+1):length(opt_data)]
  dknn <- as.numeric(dknn[1])
  
  # Extract distances values 
  dist <- opt_data[(which(opt_data == ">>dists")+1):(which(opt_data == ">>thresholds")-1)]
  dist <- strsplit(dist, "\t")
  dist <- do.call(rbind, dist)
  # Convert it in dataframe
  dist <- as.data.frame(dist)
  
  # Modify column names 
  dist_name <- c("Observed")
  for (i in 1:(ncol(dist)-1)) {
    dist_name_temp <- paste0("simulation_", i, "_distance")
    dist_name <- c(dist_name, dist_name_temp)
  }
  names(dist) <- dist_name
  
  # Convert column from factor to numeric
  dist <- sapply(dist[1:4], function(x) as.numeric(as.character(x)))
  dist <- as.data.frame(dist)
  
  # Prepare dataframe for ggplot
  dist_p <- melt(dist, variable.name = "type", value.name = "distance_value", id.vars = NULL)
  
  dist_p[,"type"] <- as.character(dist_p[,"type"])
  dist_p[grep(dist_p[,"type"], pattern = "sim*"), "type"] <- "Simulated"
  
  # plot density of distance values for the observed and simulated conditions
  p <- ggplot(data = dist_p, aes(x = distance_value, color = type)) +
    stat_density(aes(linetype = type, size = type), geom = "line", position = "identity") +
    scale_linetype_manual(breaks = c("Observed", "Simulated"), values = c("solid", "longdash")) +
    scale_size_manual(values = c(1, 0.8)) +
    scale_color_manual(values = c("#006D77", "#83C5BE")) +
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
    geom_vline(aes(xintercept = dknn),
               color = "#E29578",
               linetype = "dotdash",
               size = 0.5) +
    geom_text(mapping = aes(x = dknn,
                            y = 0,
                            label = "Critical distance with KNN",
                            hjust = -0.5,
                            vjust = -1,
                            angle = 90),
              color = "#E29758")
  
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
#' @param use_top_genes A logical to indicate whether to use highly similar genes in the slot top_genes of ClusterSet.
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
#' m <- matrix(rnorm(6000), nc=20)
#' m[1:100,1:10] <- m[1:100,1:10] + 4
#' m[101:200,11:20] <- m[101:200,11:20] + 3
#' m[201:300,5:15] <- m[201:300,5:15] + -2
#' 
#' res <- DBFMCL(data=m,
#' distance_method="pearson",
#' av_dot_prod_min = 0,
#' inflation = 1.2,
#' k=25,
#' fdr = 10)
#' 
#' plot_heatmap(object = res)
#' 

#' @rdname plot_heatmap

plot_heatmap <- function(object,
                         center = TRUE,
                         ceil = 1,
                         floor = -1,
                         cell_order = NULL,
                         cell_ordering_method = "hclust",
                         use_top_genes = FALSE,
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
    if (cell_ordering_method == "hclust"){
      print_msg("Ordering cells based on hierarchical clustering.", msg_type="DEBUG")
      m_dist <- as.dist(1 - cor(m, method = 'pearson'))
      m_clust <- hclust(m_dist, method = 'average')
      m <- m[,m_clust$order]
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
  
  
  # Add blank row to separate feature clusters in heatmap
  ## Create blank row
  blank_row <- data.frame()
  blank_row[1:line_size, 1:ncol(object@data)] <- NA
  row.names(blank_row) <- make.names(rep("w_line", line_size), unique=TRUE)
  colnames(blank_row) <- colnames(object@data)
  
  ## Insert blank row in matrix
  m_blank <- data.frame()
  for (i in 1:length(object@size)) {
    row_start <- sum(object@size[0:(i-1)])+1
    row_end <- sum(object@size[1:i])
    m_blank_loop <- rbind(m[row_start:row_end,], blank_row)
    #rownames(test)[(nrow(test)-line_size+1):nrow(test)] <- paste(rep(" ", 2), collapse = '')
    m_blank <- rbind(m_blank, m_blank_loop)
  }
  m <- as.matrix(m_blank)
  
  # Reduce m rows to only keep genes from top_genes
  if(use_top_genes) {
    if (nrow(object@top_genes) == 1 &
        ncol(object@top_genes) == 1 &
        is.na(object@top_genes[1,1])) {
      stop(paste0("The slot top_genes of the input ClusterSet object is empty. Be sure to run top_genes() before."))
    }
    genes_top <- unlist(as.data.frame(t(object@top_genes)), use.names = FALSE)
    genes_top <- genes_top[!is.na(genes_top)]
    m <- m[genes_top,]
  }
  
  
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
  
  
  
  return(htmp)
}