################################################################
##        Main script of the R PACKAGE : scigenex
##
## Authors : J. BAVAIS, BERGON A, 
##  with the collaboration of LOPEZ F., TEXTORIS J. and PUTHIER D.
##
##
## R CMD SHLIB dbf.c -o dbf
#################################################################
##    UTILS (a set of useful function)
#################################################################
# 0 : No message
# 1 : info
# 2 : DEBUG

VERBOSITY_SCIGENEX = 3

print_msg <- function(msg, msg_type="INFO"){
  if(msg_type == "INFO")
    if(VERBOSITY_SCIGENEX > 0)
      cat(paste("|-- ", msg, "\n"))
  if(msg_type == "DEBUG")
    if(VERBOSITY_SCIGENEX > 1)
      cat(paste("|-- ", msg, "\n"))
  if(msg_type == "WARNING")
      cat(paste("|-- ", msg, "\n"))
}
#################################################################
##    DEFINITION OF A SPECIFIC CLASS OBJECT : ClusterSet
#################################################################

library(ggplot2)
library(reshape2)
library(dplyr)
library(igraph)
library(iheatmapr)

#' @title
#' ClusterSet
#' @description
#' This class is a representation of a partitioning algorithm and is intented to store gene clusters.
#' @slot name character. The original input file name (if applicable).
#' @slot opt_name character. The name is file containing distance and cutting threshold.
#' @slot data matrix. The matrix containing the filtered/partitionned data.
#' @slot cluster vector. Mapping of row/genes to clusters.
#' @slot size vector. The size of each cluster.
#' @slot top_genes matrix The highly co-expressed genes of each gene clusters.
#' @slot center matrix. The centers of each clusters.
#' @slot parameters list. The parameter used.
#' @slot algorithm vector. The algorithm used to produce the clusters.
#' @slot cell_types vector. The cell types.
#' @slot cell_colors vector. The cell types to color mapping.
#' @slot cell_order vector. How cell should be ordered.
#' @slot cluster_annotations list. Functional annotation of clusters.
#' @return A ClusterSet object.
#' @export
#'
#' @examples
#' 
#' \dontrun{
#'   m <- matrix(rnorm(80000), nc=20)
#'    m[1:100,1:10] <- m[1:100,1:10] + 4
#'    m[101:200,11:20] <- m[101:200,11:20] + 3
#'    m[201:300,5:15] <- m[201:300,5:15] + -2
#'    res <- DBFMCL(data=m,
#'                  distance_method="pearson",
#'                  av_dot_prod_min = 0,
#'                  inflation = 1.2,
#'                  k=25,
#'                  fdr = 10)
#' plot_clust(res, ceil = 10, floor = -10)
#' plot_clust(res, type="tile", ceil = 10, floor = -10)
#' write_clust(res, filename_out = "ALL.sign.txt")
#'   is(res)
#' }
#'               
setClass("ClusterSet",
  representation = list(
    name = "character",
    opt_name = "character",
    data = "matrix",
    cluster = "vector",
    size = "vector",
    top_genes = "matrix",
    center = "matrix",
    parameters = "list",
    algorithm = "vector",
    cell_types = "vector",
    cell_colors = "vector",
    cell_order = "vector",
    cluster_annotations = "list"
  ),
  prototype = list(
    name = character(),
    data = matrix(nr = 0, nc = 0),
    cluster = numeric(),
    size = numeric(),
    top_genes = matrix(),
    center = matrix(nc = 0, nr = 0),
    parameters = list(),
    algorithm = character(),
    cell_types = vector(),
    cell_colors = vector(),
    cell_order = vector(),
    cluster_annotations = list()
  )
)

#################################################################
##    REDEFINE SHOW METHOD FOR CLASS OBJECT : ClusterSet
#################################################################

setMethod(
  "show", signature("ClusterSet"),

  function(object) {
    
    cat("\t\tAn object of class ClusterSet\n")
    cat("\t\tName:", slot(object, "name"), "\n")
    cat("\t\tMemory used: ", object.size(object), "\n")
    cat("\t\tNumber of cells: ", ncol(slot(object, "data")), "\n")
    cat(
      "\t\tNumber of informative genes: ",
      nrow(slot(object, "data")), "\n"
    )
    cat("\t\tNumber of clusters: ", length(slot(object, "size")), "\n")
    cat("\t\tThis object contains the following informations:\n")

    for (i in slotNames(object)) {
      cat("\t\t\t - ", i, "\n")
    }
    if (length(slot(object, "parameters")) > 0) {
      for (i in 1:length(slot(object, "parameters"))) {
        cat(
          "\t\t\t\t * ", names(slot(object, "parameters"))[[i]],
          " = ", slot(object, "parameters")[[i]], "\n"
        )
      }
    }
  }
)


#################################################################
##    REDEFINE dim/ncol/nrow METHOD FOR CLASS OBJECT : ClusterSet
#################################################################

#' @title
#' ncol
#' @description
#' The number of columns of a ClusterSet object.
#' @param x A ClusterSet object.
#' @return The number of columns.
#' @export
#'
#' @examples
#'
#' \dontrun{
#'   m <- matrix(rnorm(80000), nc=20)
#'    m[1:100,1:10] <- m[1:100,1:10] + 4
#'    m[101:200,11:20] <- m[101:200,11:20] + 3
#'    m[201:300,5:15] <- m[201:300,5:15] + -2
#'    res <- DBFMCL(data=m,
#'                 distance_method="pearson",
#'                  av_dot_prod_min = 0,
#'                  inflation = 1.2,
#'                  k=25,
#'                  fdr = 10)
#'   ncol(res)
#' }
#'
setMethod(
  "ncol", signature("ClusterSet"),

  function(x) {

      return(ncol(x@data))

      }
)

#' @title
#' nrow
#' @description
#' The number of rows of a ClusterSet object.
#' @param x A ClusterSet object.
#' @return The number of rows.
#' @export
#'
#' @examples
#'
#' \dontrun{
#'   m <- matrix(rnorm(80000), nc=20)
#'    m[1:100,1:10] <- m[1:100,1:10] + 4
#'    m[101:200,11:20] <- m[101:200,11:20] + 3
#'    m[201:300,5:15] <- m[201:300,5:15] + -2
#'    res <- DBFMCL(data=m,
#'                  distance_method="pearson",
#'                  av_dot_prod_min = 0,
#'                  inflation = 1.2,
#'                  k=25,
#'                  fdr = 10)
#'   nrow(res)
#' }
#'
setMethod(
  "nrow", signature("ClusterSet"),

  function(x) {

      return(nrow(x@data))

      }
)


#' @title
#' dim
#' @description
#' The number of rows/columns of a ClusterSet object.
#' @param x A ClusterSet object.
#' @return The number of rows/columns.
#' @export
#'
#' @examples
#'
#' \dontrun{
#'   m <- matrix(rnorm(80000), nc=20)
#'    m[1:100,1:10] <- m[1:100,1:10] + 4
#'    m[101:200,11:20] <- m[101:200,11:20] + 3
#'    m[201:300,5:15] <- m[201:300,5:15] + -2
#'   res <- DBFMCL(data=m,
#'                  distance_method="pearson",
#'                  av_dot_prod_min = 0,
#'                  inflation = 1.2,
#'                  k=25,
#'                 fdr = 10)
#'   dim(res)
#' }
#'
setMethod(
  "dim", signature("ClusterSet"),

  function(x) {

      return(dim(x@data))

      }
)


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




#################################################################
##    Define the write function for class ClusterSet
#################################################################


#' @title
#' write_clust
#' @description
#' Write a ClusterSet into a flat file.
#' @param object ClusterSet. 
#' @param filename_out The outfile name.
#' @param out_path The path to the file.
#' @param nb_na_row Number of separating rows (containing NAs).
#' @return Write a file.
#' @export
#'
#' @examples
#' # see online help.
setGeneric("write_clust",
           
           function(object, 
                    filename_out = NULL,
                    out_path = ".",
                    nb_na_row=3) {
                standardGeneric("write_clust")
})


#' @rdname write_clust
setMethod(
  "write_clust",
  signature(object = "ClusterSet"),
  function(object,
           filename_out = NULL,
           out_path = ".",
           nb_na_row=5) {
           
    if (out_path == ".") out_path <- getwd()

    if (is.null(filename_out)) {
      filename_out <- "exprs.dataMods.txt"
    }

    data <- object@data
    nb <- 0
    dataT <- c("clusters", colnames(data))

    ## processing data
    for (i in 1:length(object@size)) {

      print_msg(paste("Cluster ", i, " --> ", object@size[i], " probes"), msg_type="INFO")

      subData <- data[object@cluster == i, ]
      subData <- cbind(rownames(subData), subData)
      if (nb_na_row > 0){
        intLine <- matrix(rep(NA, (ncol(data) + 1)*nb_na_row), nrow = nb_na_row)
        dataT <- rbind(dataT, subData, intLine)
      }
      nb <- nb + 1

    }


    ## exporting results
    print_msg("Exporting results", msg_type="DEBUG")
    write.table(dataT, file.path(out_path, filename_out),
      col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE
    )

      print_msg(paste("--> Creating file : ",
                      file.path(out_path, filename_out)),
                msg_type="DEBUG")
  }
)

#################################################################
##    DBF-MCL
#################################################################

#' @title
#' The "Density Based Filtering and Markov CLustering" algorithm (DBF-MCL).
#' @description
#' DBFMCL is a tree-steps adaptative algorithm that \emph{(i)} find elements
#' located in dense areas (DBF), \emph{(ii)} uses selected items to construct a
#' graph, \emph{(iii)} performs graph partitioning using the Markov CLustering
#' Algorithm (MCL).
#'
#' This function requires installation of the mcl program
#' (\url{http://www.micans.org/mcl}). See "Warnings" section for more
#' informations.
#'
#' When analyzing a noisy dataset, one is interested in isolating dense regions
#' as they are populated with genes/elements that display weak distances to
#' their nearest neighbors (i.e. strong profile similarities). To isolate these
#' regions DBF-MCL computes, for each gene/element, the distance with its kth
#' nearest neighbor (DKNN).In order to define a critical DKNN value that will
#' depend on the dataset and below which a gene/element will be considered as
#' falling in a dense area, DBF-MCL computes simulated DKNN values by using an
#' empirical randomization procedure. Given a dataset containing n genes and p
#' samples, a simulated DKNN value is obtained by sampling n distance values
#' from the gene-gene distance matrix D and by extracting the kth-smallest
#' value. This procedure is repeated n times to obtain a set of simulated DKNN
#' values S. Computed distributions of simulated DKNN are used to compute a FDR
#' value for each observed DKNN value. The critical value of DKNN is the one
#' for which a user-defined FDR value (typically 10\%) is observed. Genes with
#' DKNN value below this threshold are selected and used to construct a graph.
#' In this graph, edges are constructed between two genes (nodes) if one of
#' them belongs to the k-nearest neighbors of the other. Edges are weighted
#' based on the respective coefficient of correlation (\emph{i.e.}, similarity)
#' and the graph obtained is partitioned using the Markov CLustering Algorithm
#' (MCL).
#'
#' @param data a \code{matrix}, \code{data.frame} or \code{Seurat} object.
#' @param filename a character string representing the file name.
#' @param name a prefix for the names of the intermediary files created by DBF
#' and MCL.
#' @param path a character string representing the data directory where
#' intermediary files are to be stored. Default to current working directory.
#' @param output_path a character string representing the data directory where
#' output files will be stored. Default to current working directory.
#' @param optional_output if TRUE then DBF generate optional output file in the
#' specified output_path directory. This file contains observed and simulated 
#' distances, cutting threshold, number of kept genes and FDR value.
#' @param mcl_threads An integer to determine number of threads for mcl algorithm.
#' @param distance_method a method to compute the distance to the k-th nearest
#' neighbor. One of "pearson" (Pearson's correlation coefficient-based
#' distance), "spearman" (Spearman's rho-based distance), "euclidean".
#' @param av_dot_prod_min Any cluster with average dot product below this value is discarded. This allow to delete
#' clusters in which correlation is influenced/supported by very few samples (typically 1).
#' @param min_cluster_size Minimum number of element inside a cluster. MCL tend to create lots of clusters with
#' very few (e.g 2) objects.
#' @param silent if set to TRUE, the progression of distance matrix calculation
#' is not displayed.
#' @param k the neighborhood size.
#' @param random the number of simulated distributions S to compute. By default
#' \code{random = 3}.
#' @param memory_used size of the memory used to store part of the distance
#' matrix. The subsequent sub-matrix is used to computed simulated distances to
#' the k-th nearest neighbor (see detail section).
#' @param fdr an integer value corresponding to the false discovery rate
#' (range: 0 to 100).
#' @param inflation the main control of MCL. Inflation affects cluster
#' granularity. It is usually chosen somewhere in the range \code{[1.2-5.0]}.
#' \code{inflation = 5.0} will tend to result in fine-grained clusterings, and
#' whereas \code{inflation = 1.2} will tend to result in very coarse grained
#' clusterings. By default, \code{inflation = 2.0}. Default setting gives very
#' good results for microarray data when k is set between 70 and 250.
#' @param set.seed specify seeds for random number generator.
#' @return a ClusterSets class object.
#' @section Warnings: With the current implementation, this function only works
#' only on UNIX-like plateforms.
#'
#' MCL should be installed. One can used the following command lines in a
#' terminal:
#'
#' \code{# Download the latest version of mcl (the script has been tested
#' successfully with the 06-058 version).}
#' \code{wget http://micans.org/mcl/src/mcl-latest.tar.gz}
#' \code{# Uncompress and install mcl}
#' \code{tar xvfz mcl-latest.tar.gz}
#' \code{cd mcl-xx-xxx}
#' \code{./configure}
#' \code{make}
#' \code{sudo make install}
#' \code{# You should get mcl in your path}
#' \code{mcl -h}
#' @author Bergon A., Bavais J., Textoris J., Granjeaud S., Lopez F and Puthier
#' D.
#' @references
#' - Van Dongen S. (2000) A cluster algorithm for graphs. National
#' Research Institute for Mathematics and Computer Science in the 1386-3681.
#' - Lopez F.,Textoris J., Bergon A., Didier G., Remy E., Granjeaud
#' S., Imbert J. , Nguyen C. and Puthier D. TranscriptomeBrowser: a powerful
#' and flexible toolbox to explore productively the transcriptional landscape
#' of the Gene Expression Omnibus database. PLoSONE, 2008;3(12):e4001.
#' @keywords clustering, gene expression, classification, MCL.
#' @examples
#'
#' \dontrun{
#'   m <- matrix(rnorm(80000), nc=20)
#'   m[1:100,1:10] <- m[1:100,1:10] + 4
#'    m[101:200,11:20] <- m[101:200,11:20] + 3
#'    m[201:300,5:15] <- m[201:300,5:15] + -2
#'    res <- DBFMCL(data=m,
#'                  distance_method="pearson",
#'                  av_dot_prod_min = 0,
#'                  inflation = 1.2,
#'                  k=25,
#'                  fdr = 10)
#' plot_clust(res, ceil = 10, floor = -10)
#' plot_clust(res, type="tile", ceil = 10, floor = -10)
#' write_clust(res, filename_out = "ALL.sign.txt")
#' }
#'
#' @export DBFMCL
DBFMCL <- function(data = NULL, 
                   filename = NULL, 
                   path = ".",
                   output_path = ".",
                   optional_output = TRUE,
                   mcl_threads=1,
                   name = NULL,
                   distance_method = c("pearson", "spearman", "euclidean"),
                   av_dot_prod_min=2,
                   min_cluster_size=10,
                   silent = FALSE,
                   k = 50,
                   random = 3, 
                   memory_used = 1024,
                   fdr = 10,
                   inflation = 8,
                   set.seed = 123) {
  
  ## testing the system
  if (.Platform$OS.type == "windows") {
    stop("\t--> A unix-like OS is required to launch mcl and cluster programs.")
  }
  
  ## getting parameters
  data_source <- get_data_4_DBFMCL(data = data, filename = filename, path = path)
  data_matrix <- data_source$data
  
  # A simple function to create a random string
  create_rand_str <- function() {
    v = c(sample(LETTERS, 3, replace = TRUE),
          sample(0:9, 4, replace = TRUE),
          sample(letters, 3, replace = TRUE))
    return(paste0(sample(v),collapse = ""))
  }
  
  if (is.null(name)) name <- data_source$name
  if (is.null(name)) name <- create_rand_str()
  
  # Put the current working directory in output_path or path
  if(path == ".") {
    path <- getwd()
  }
  if(output_path == ".") {
    output_path <- getwd()
  }
  # Check if output directory exists. If not stop the command.
  if(!file.exists(output_path)){
    stop("Output directory provided does not exist.")
  }
  
  
  
  distance_method <- match.arg(distance_method)
  txt <- paste("\n\tInflation: ", inflation, sep = "")
  
  ## writting all parameters
  
  cat(
    "The following parameters will be used :",
    "\n\tWorking directory: ", path,
    "\n\tOuput directory: ", output_path,
    "\n\tName: ", name,
    "\n\tDistance method: ", distance_method,
    "\n\tMinimum average dot product for clusters: ", av_dot_prod_min,
    "\n\tMinimum cluster size: ", min_cluster_size,
    "\n\tNumber of neighbors: ", k,
    "\n\tNumber of randomizations: ", random,
    "\n\tFDR: ", fdr, "%", txt,
    "\n\tVisualize standard outputs from both mcl and cluster",
    "commands: ", silent,
    "\n\tMemory used : ", memory_used, "\n\n"
  )
  
  
  ## DBF algorithm, returns a ClusterSet object
  obj <- DBF(data_matrix,
             output_path = output_path,
             name,
             optional_output = optional_output,
             distance_method = distance_method,
             silent = silent,
             k = k,
             random = random,
             memory_used = memory_used,
             fdr = fdr,
             set.seed = set.seed
  )
  
  dbf_out_file <- paste0(output_path, "/", name, ".dbf_out.txt")
  dbf_out_file <- gsub(pattern = "//", replacement = "/", x = dbf_out_file)
  mcl_out_file <- paste0(output_path, "/", name, ".mcl_out.txt")
  mcl_out_file <- gsub(pattern = "//", replacement = "/", x = mcl_out_file)
  
  print_msg("DBF completed. Starting MCL step.", msg_type="DEBUG")
  
  if (length(readLines(dbf_out_file)) > 0) {
    
    ## Launching mcl (command line)
    mcl_system_cmd(name, inflation = inflation, input_path = output_path, silent = silent, threads = mcl_threads)
    
    
    print_msg(paste0("Reading MCL output: ", mcl_out_file), msg_type="DEBUG")
    
    ## getting mcl results into the ClusterSet object
    mcl_cluster <- readLines(mcl_out_file)
    gene_list <- NULL
    clusters <- NULL
    size <- NULL
    nb <- 0
    nb_cluster_deleted <- 0
    
    for (i in 1:length(mcl_cluster)) {
      h <- unlist(strsplit(mcl_cluster[i], "\t"))
      cur_clust <- data_matrix[h,]
      cur_clust[cur_clust > 0 ] <- 1
      cur_dot_prod <- cur_clust %*% t(cur_clust)
      
      if(mean(cur_dot_prod) > av_dot_prod_min & length(h) > min_cluster_size){
        
        nb <- nb + 1
        gene_list <- c(gene_list, h)
        clusters <- c(clusters, rep(nb, length(h)))
        if (is.null(size)) {
          size <- length(h)
        }
        else {
          size <- c(size, length(h))
        }
      }else{
        nb_cluster_deleted <- nb_cluster_deleted + 1
      }
    }
    print_msg(paste(nb, " clusters conserved after MCL partitioning."),
              msg_type="INFO")
    print_msg(paste(nb_cluster_deleted,
                    " clusters filtered out from MCL partitioning (size and mean dot product)."),
              msg_type="INFO")
    
    
    ## build ClusterSet object
    if (nb > 0) {
      obj@name <- name
      obj@data <- as.matrix(data_matrix[gene_list, ])
      names(clusters) <- rownames(obj@data)
      obj@cluster <- clusters
      obj@size <- size
      
      centers <- matrix(ncol = ncol(data_matrix), nrow = nb)
      ## calcul of the mean profils
      for (i in 1:nb) {
        centers[i, ] <- apply(obj@data[obj@cluster == i, ],
                              2, mean,
                              na.rm = TRUE
        )
      }
      obj@center <- centers
      
      ## add DBFMCL parameters used to build this object
      obj@parameters <- list(
        distance_method = distance_method,
        k = k,
        random = random,
        fdr = fdr,
        set.seed = set.seed,
        inflation = inflation
      )
    }
  }
  else {
    stop("\t--> There is no conserved gene.\n\n")
  }
  return(obj)
  
}

###############################################################
##    COMPUTE DBF algorithm
###############################################################

#' @title
#' DBF
#' @description
#' This function is an internal function used by \code{\link{DBFMCL}} to detect
#' informative elements (\emph{i.e.}, those that belong to dense regions). User
#' should not use this function. Instead they can use the \code{\link{DBFMCL}}
#' function with \code{clustering} argument set to \code{FALSE}.
#'
#' See \code{\link{DBFMCL}}
#'
#' @param data a matrix or data.frame
#' @param output_path a character string representing the data directory where
#' output files will be stored. Default to current working directory.
#' @param name a prefix for the file name
#' @param optional_output if TRUE then DBF generate optional output file in the
#' specified output_path directory. This file contains observed and simulated 
#' distances, cutting threshold, number of kept genes and FDR value.
#' @param distance_method a method to compute the distance to the k-th nearest
#' neighbor. One of "pearson" (Pearson's correlation coefficient-based
#' distance), "spearman" (Spearman's rho-based distance) or "euclidean".
#' @param silent if set to TRUE (default), the progression of distance matrix
#' calculation is not displayed.
#' @param k the neighborhood size.
#' @param random the number of simulated distributions S to compute. By default
#' \code{random = FALSE}.
#' @param fdr a value for the false discovery rate.
#' @param memory_used size of the memory used to store part of the distance
#' matrix. The subsequent sub-matrix is used to computed simulated distances to
#' the k-th nearest neighbor (see detail section).
#' @param set.seed specify seeds for random number generator.
#' @section Warnings: Works only on UNIX-alikes platforms.
#' @author Bergon A., Bavais J., Textoris J., Granjeaud S., Lopez F and Puthier
#' D.
#' @seealso \code{\link{DBFMCL}}
#' @references Lopez F.,Textoris J., Bergon A., Didier G., Remy E., Granjeaud
#' S., Imbert J. , Nguyen C. and Puthier D. TranscriptomeBrowser: a powerful
#' and flexible toolbox to explore productively the transcriptional landscape
#' of the Gene Expression Omnibus database. PLoSONE, 2008;3(12):e4001.
#' @keywords manip
#' @export DBF
DBF <- function(data,
                output_path = ".",
                name = NULL,
                distance_method = c("spearman", "pearson", "euclidean"),
                silent = FALSE,
                k = 100,
                fdr = 10,
                set.seed = 123) {

  ## testing the system
  if (.Platform$OS.type != "windows") {
    if (!is.null(data)) {
      ## getting data and parameters
      if (output_path == ".") {output_path <- getwd()}
      if (is.null(name)) name <- "exprs"
      data <- get_data_4_DBFMCL(data = data)$data
      row <- rownames(data)
      col <- colnames(data)
      distance_method <- match.arg(distance_method)

      ## transforming data into double
      data <- apply(data, 2, as.double)

      if (silent) {
        print_msg(paste0("Computing distances to the kth-nearest neighbors ",
                         "and associated FDR values... \n"),
                  msg_type = "INFO")

      }
      
      # Directory and name of the principal output
      outfile <- paste(output_path, "/", name, ".dbf_out.txt", sep = "")
      outfile <- gsub(pattern = "//", replacement = "/", x = outfile)
      
      
      
      #################### Correlation and distance matrices
      # Remove genes with 0 values for all cells
      select_for_correlation <- data_matrix[-c(which(rowSums(data_matrix) == 0)),]
      
      # Compute gene-gene pearson correlation matrix
      if(distance_method == "pearson") {
        cor_matrix <- corSparse(t(select_for_correlation))
        rownames(cor_matrix) <- rownames(select_for_correlation)
        colnames(cor_matrix) <- rownames(select_for_correlation)
      } else {
        print_msg(msg_type = "WARNING",
                  msg = "Distance used is not provided.")
      }
      
      # Remove the diagonal and the upper triangular matrix (values is replaced by NA)
      cor_matrix[lower.tri(cor_matrix, diag=TRUE)] <- NA
      
      # Transform correlation matrix into distance matrix (values between 0 and 2)
      dist_matrix <- 2 - (cor_matrix + 1)
      
      
      #Create a dataframe with a column that contains all the gene ID
      df_dknn <- data.frame("gene_id" = rownames(dist_matrix))
      l_knn <- list()
      
      #################### DKNN for each genes
      # Extract the DKNN for each gene
      for (gene in df_dknn[,"gene_id"]){
        #Create a vector with all the correlation values for one gene(i)
        gene_dist <- c(subset(dist_matrix[gene,], !is.na(dist_matrix[gene,])),
                       subset(dist_matrix[,gene], !is.na(dist_matrix[,gene]))) 
        
        #Reorder the pearson correlation values (increasing order)
        row_dknn <- order(gene_dist, decreasing=F)[1:k]
        gene_dknn <- gene_dist[row_dknn]
        
        #Store the results in a list
        l_knn[[gene]] <- gene_dknn
        
        #Select the kth pearson correlation values. This value corresponds to the DKNN of the gene(i)
        df_dknn[which(df_dknn[,"gene_id"] == gene) ,"DKNN_cor"] <- gene_dknn[k]
      }
      
      
      #################### DKNN simulation
      # Extract the distance values from the distance matrix
      dist_values <- dist_matrix[!is.na(dist_matrix)]
      nb_of_gene_sim <- nrow(select_for_correlation) * 1.5 #REMPLACER PAR NB DE GENE x2
      sim_dknn <- vector()
      
      # Generate simulated distances
      for (sim_nb in 1:nb_of_gene_sim) {
        
        # Randomly sample distances for one simulated gene
        nb_gene <- nrow(select_for_correlation)
        dist_sim <- sample(dist_values, size=nb_gene, replace=FALSE)
        
        # Extract the k nearest neighbors of these simulated gene
        dist_sim <- dist_sim[order(dist_sim)]
        sim_dknn[sim_nb] <- dist_sim[k]
        
      }
      
      
      #################### Determine the DKNN threshold (or critical distance)
      # Order genes by DKNN values
      df_dknn <- df_dknn[order(df_dknn[,"DKNN_cor"]),]
      df_dknn[, c("nb_dknn_sim", "nb_dknn_obs", "ratio_sim_obs")] <- 0
      
      
      #Recuperer les valeurs sim_dknn < obs_dknn(i)
      for(i in 2:nb_gene) {
        
        #Compute the number of simulated DKNN values under DKNN value at rank i
        #If there is no simulated DKNN values under DKNN value at rank i, put 0 in nb_dknn_sim
        if (min(sim_dknn) < df_dknn[i,"DKNN_cor"]) {
          nb_dknn_sim <- sim_dknn[which(sim_dknn < df_dknn[i,"DKNN_cor"])]
          nb_dknn_sim <- length(nb_dknn_sim)
          df_dknn[i,"nb_dknn_sim"] <- nb_dknn_sim
          
        }else {
          nb_dknn_sim <- 0 
          df_dknn[i,"nb_dknn_sim"] <- nb_dknn_sim 
          
        }
        #Compute the number of observed DKNN values under DKNN value at rank i
        nb_dknn_obs <- length(df_dknn[which(df_dknn[,"DKNN_cor"] < df_dknn[i,"DKNN_cor"]), "DKNN_cor"]) #Nombre de valeurs de dknn observees inferieures a la valeur dknn_obs(i)
        df_dknn[i,"nb_dknn_obs"] <- nb_dknn_obs
        
        #Compute the ratio between number of simulated DKNN values and the number of observed DKNN values under DKNN value at rank i
        df_dknn[i,"ratio_sim_obs"] <- nb_dknn_sim/(nb_dknn_obs + nb_dknn_sim)
      }
      
      #################### Select genes with a distance value under critical distance
      selected_genes <- df_dknn[which(df_dknn[,"ratio_sim_obs"] < fdr*0.01),]
      
      
      #################### Outputs
      # Create the ClusterSet object
      obj <- new("ClusterSet")
      obj@algorithm <- "DBFMCL"
      
      if (length(selected_genes[,"gene_id"]) > 0) {
        obj@opt_name <- "extra_output-dbfAll.txt"
        obj@data <- as.matrix(data[selected_genes[,"gene_id"],])
        obj@cluster <- rep(1, nrow(obj@data))
        obj@size <- nrow(obj@data)
        obj@center <- matrix(
          apply(obj@data[obj@cluster == 1, ],
                2,
                mean,
                na.rm = TRUE
          ),
          nrow = 1
        )
      }
      
      
      
      
      
      
      return(obj)
    } else {
      stop("\t--> Please provide a matrix...\n\n")
    }
  }
  else {
    stop("\t--> A unix-like OS is required to launch mcl and cluster programs.")
  }
}

##############################################################
##    MCL
##############################################################


#' @title
#' Invokes the command line version of Markov CLustering algorithm (MCL).
#' @description
#'  This function invokes the mcl system command. MCL is a clustering algorithm
#' for graphs that was developped by Stijn van Dongen (see references for
#' further informations).
#' @param name a character string corresponding to the file name.
#' @param inflation the main control of MCL. Inflation affects cluster
#' granularity. It is usually chosen somewhere in the range \code{[1.2-5.0]}.
#' \code{inflation = 5.0} will tend to result in fine-grained clusterings, and
#' whereas \code{inflation = 1.2} will tend to result in very coarse grained
#' clusterings. By default, \code{inflation = 2.0}. Default setting gives very
#' good results for microarray data when k is set around 100.
#' @param input_path a character string representing the directory path of 
#' the input file used by mcl. Default is the current working directory.
#' @param silent if set to TRUE, the progression of the MCL partitionning is
#' not displayed.
#' @param threads The number of threads to use.
#' @return Returns a file with the ".mcl\_out.txt" extension.
#' @section warning: Works only on UNIX-like plateforms. MCL should be
#' installed. The following command lines can be used for installation.
#'
#' \code{# Download the latest version of mcl (RTools4TB has been tested
#' successfully with the 06-058 version).}
#' \code{wget http://micans.org/mcl/src/mcl-latest.tar.gz}
#' \code{# Uncompress and install mcl}
#' \code{tar xvfz mcl-latest.tar.gz}
#' \code{cd mcl-xx-xxx}
#' \code{./configure}
#' \code{make}
#' \code{sudo make install}
#' \code{# You should get mcl in your path}
#' \code{mcl -h}
#' @author Bergon A., Lopez F., Textoris J., Granjeaud S. and Puthier D.
#' @references Stijn van Dongen. A cluster algorithm for graphs.  Technical
#' Report INS-R0010, National Research Institute for Mathematics and Computer
#' Science in the Netherlands, Amsterdam, May 2000.
#' \url{http://www.cwi.nl/ftp/CWIreports/INS/INS-R0010.ps.Z}
#' @keywords manip
#' @export mcl_system_cmd
mcl_system_cmd <- function(name, inflation = 2.0, input_path = ".", silent = FALSE, threads=1) {
  ## testing the system
  if (.Platform$OS.type != "windows") {

    ## Testing mcl installation
    if (system("mcl --version | grep 'Stijn van Dongen'", intern = TRUE) > 0) {
      if (!silent) {
        cat("Running mcl (graph partitioning)... \n")
        verb <- ""
      }
      else {
        verb <- "-V all "
      }
      if (inflation != 2) {
        i <- paste("-I ", as.character(round(inflation, 1)), sep = "")
      } else {
        i <- "-I 2.0"
      }
      threads <- paste("-te", threads, sep = " ")
      ## launching mcl program
      cmd <- paste0("mcl ",
                   input_path, "/", name, ".dbf_out.txt ",
                   i,
                   " --abc -o ",
                   input_path, "/", name, ".mcl_out.txt ",
                   verb,
                   threads)
      cmd <- gsub(pattern = "//", replacement = "/", x = cmd)
      system(cmd)

      if (!silent) {
        print_msg("Done", msg_type="INFO")
        print_msg(paste0("creating file : ",
                        file.path(getwd(), paste(name, ".mcl_out.txt", sep = ""))),
                  msg_type="INFO")
      }
    } else {
      stop(
        "\t--> Please install mcl on your computer...\n",
        "\t--> You can download it from : 'http://www.micans.org/mcl/'\n\n"
      )
    }
  }
  else {
    stop("--> A unix-like OS is required to launch mcl and cluster programs.")
  }
}


#################################################################
##    Getting data
#################################################################

#' @title
#' Fetch an expression matrix from a file, dataframe or Seurat object.
##' @description
#' This function fetchs an expression matrix from a file, dataframe or Seurat object.
#' @param data A \code{matrix}, \code{data.frame} or \code{Seurat} object.
#' @param filename A character string representing the file name.
#' @param path A character string representing the data directory where
#' intermediary files are to be stored. Default to current working directory.
#'
#' @return A list containing a matrix and the filename (if filename argument is used).
#' @export get_data_4_DBFMCL
#'
#' @examples
#' 
#' \dontrun{
#' ## with an artificial dataset
#'
#' m <- matrix(rnorm(80000), nc = 20)
#' res <- get_data_4_DBFMCL(data=m)
#' }
#' 
get_data_4_DBFMCL <- function(data = NULL, filename = NULL, path = ".") {

  ## getting matrix (probesID vs SamplesID)
  if (!is.null(data)) {
    if (inherits(data, "Seurat")) {
      data <- as.matrix(data@assays$RNA@data)
    }
    else if (is.data.frame(data)) {
      data <- as.matrix(data)
    }
    if (!is.matrix(data)) {
      stop(
        "\t--> Please provide a Seurat Object, a data.frame",
        " or a matrix.\n"
      )
    }
    name <- NULL
  }
  else {
    if (!is.null(filename)) {
      data <- as.matrix(read.table(file.path(path, filename),
        sep = "\t", header = TRUE, row.names = 1, quote = ""
      ))
      name <- unlist(strsplit(filename, "\\."))[1]
    }
    else {
      stop(
        "\t--> Please provide an ExpressionSet, a data.frame, ",
        "a matrix or a tabular file\n"
      )
    }
  }
  ## adding dimnames if not provided
  if (is.null(rownames(data))) {
    print_msg("Row names not provided. Adding.", msg_type = "DEBUG")
    rownames(data) <- paste("gene", 1:nrow(data), sep = "")
  }
  if (is.null(colnames(data))) {
    print_msg("Colum names not provided. Adding.", msg_type = "DEBUG")
    colnames(data) <- paste("sample", 1:ncol(data), sep = "")
  }

  return(list(data = data, name = name))
}

#################################################################
##    Define top_genes function for ClusterSet object
#################################################################

#' @title
#' Best co-expressed genes from each gene cluster
#' @description
#' Extract the highly co-expressed genes of each gene cluster.
#' @param object A \code{ClusterSet} object.
#' @param top A value for the number of the most similar genes in gene clusters.
#' @param cluster A vector of gene cluster identity.
#'
#' @return 
#' @export top_genes
#'
#' @examples
#' 
#' m <- matrix(rnorm(80000), nc=20)
#' m[1:100,1:10] <- m[1:100,1:10] + 4
#' m[101:200,11:20] <- m[101:200,11:20] + 3
#' m[201:300,5:15] <- m[201:300,5:15] + -2
#' 
#' res <- DBFMCL(data=m,
#'               distance_method="pearson",
#'               av_dot_prod_min = 0,
#'               inflation = 1.2,
#'               k=25,
#'               fdr = 10)
#'               
#'res <- top_genes(object = res, top=500, cluster="all")
#'res@top_genes
#' 
#' \dontrun{
#' ## with an artificial dataset
#'
#' m <- matrix(rnorm(80000), nc = 20)
#' res <- get_data_4_DBFMCL(data=m)
#' }
#' 

top_genes <- function(object,
                      top = 20,
                      cluster = "all") {
  
   if(unique(cluster == "all")) {
    cluster <- c(1:length(object@size))
   }
  
  # Display a warning message if there is less than n top genes in a gene cluster
  loop <- 0
  for (size in object@size[cluster]){
    loop <- 1 + loop
    if(top > size) {
      warning(paste0("Number of top genes is greater than the number of genes in cluster ", loop, ". All genes will be use and order by similarity rank."))
    }
  }
  
  # Initialization for the for loop
  clusters <- object@cluster
  l_cor_means <- list()
  genes_top <- matrix(ncol = top)
  
  # Extract top co-expressed genes for each gene cluster
  for (i in cluster) {
    #Extract gene names in cluster i
    genes <- names(clusters[which(clusters == i)])
    
    #Compute distances between genes in cluster i
    #Use the same distance used by DBFMCL
    dist_method <- object@parameters$distance_method
    
    if (dist_method == "pearson") {
      cor_genes <- cor(t(object@data[genes,]), method = "pearson")
      dist <- 1 - cor_genes
    }
    
    if(dist_method == "spearman") {
      cor_genes <- cor(t(object@data[genes,]), method = "spearman")
      dist <- 1 - cor_genes
    }
    
    if(dist_method == "euclidean") {
      dist <- as.matrix(dist(object@data[genes,], method = "euclidean", upper = TRUE, diag = TRUE))
      diag(dist) <- NA
    }
    
    
    #Compute mean correlation for each gene in cluster i
    dist_means <- colMeans(dist, na.rm = TRUE)
    dist_means <- dist_means[order(dist_means, decreasing = FALSE)]
    
    #Extract top genes with the highest correlation mean
    genes_top <- rbind(genes_top, names(dist_means[1:top]))
  }
  
  # Prepare genes_top matrix
  if (length(cluster) > 1) {
    genes_top <- genes_top[2:(length(cluster)+1),]
  } else {
    genes_top <- as.matrix(t(genes_top[2:(length(cluster)+1),]))
  }
  colnames(genes_top) <- make.names(rep("gene.top", top+1), unique=TRUE)[2:(top+1)]
  rownames(genes_top) <- make.names(rep("cluster", length(cluster)+1), unique=TRUE)[2:(length(cluster)+1)]
  
  # Put genes_top matrix in object@top_genes
  object@top_genes <- genes_top
  
  # Print
  print_msg(msg = "Results are stored in top_genes slot of ClusterSet object.",
            msg_type = "INFO")
  
  return(object)
}

#################################################################
##    Running gprofiler2 to perform enrichment analysis
#################################################################

#' @title
#' Perform gene enrichment analysis.
##' @description
#' Perform enrichment analysis on all MCL clusters indepentently and store the results in the cluster_annotations slot of the ClusterSet object.
#' @param object A \code{ClusterSet} object.
#' @param specie Specie name, as a concatenation of the first letter of the name and the family name, e.g human - hsapien
#'
#' @return A \code{ClusterSet} object
#' @export enrich_analysis
#'
#' @examples
#' 
#' \dontrun{
#' ## Assuming myobject is a ClusterSet object with at least 1 cluster.
#'
#' gores <- enrich_analysis(myobject)
#' }

setGeneric("enrich_analysis",
    function(object,
            specie="hsapiens") {
      standardGeneric("enrich_analysis")
})

#' @rdname enrich_analysis
setMethod("enrich_analysis",
    signature(object = "ClusterSet"),
    function(object,
            specie="hsapiens") {

      for(cluster in unique(object@cluster)){
        print(paste0("Enrichment analysis for cluster ", cluster))
        cluster_name = paste0("Cluster_", cluster)
        query = rownames(object@data[object@cluster == cluster,])
        gostres <- gost(query, organism = "hsapiens", ordered_query = FALSE, significant = TRUE, exclude_iea = T)
        object@cluster_annotations[[cluster]] = list(result = gostres$result, meta = gostres$meta)
      }
      return(object)
  }
)

#################################################################
##    Manhattan-like-plot for enrichment analysis on ClusterSet object
#################################################################

#' @title
#' Manhattan-like-plot of enrichment analysis results
##' @description
#' Retrieve enrichment analysis results from a ClusterSet object and draw a Manhattan-like-plot.
#' @param object A \code{ClusterSet} object.
#' @param clusters  A vector of cluster id to plot.
#' @param verbose Whether or not to print progression in the console.
#'
#' @return A \code{ClusterSet} object
#' @export enrich_viz
#'
#' @examples
#' 
#' \dontrun{
#' ## Assuming myobject is a ClusterSet object with at least 1 cluster.
#'
#' enrich_viz(myobject)
#' }

setGeneric("enrich_viz",
           function(object,
                    clusters = "all",
                    verbose = TRUE) {
             standardGeneric("enrich_viz")
           })

#' @rdname enrich_viz
setMethod("enrich_viz",
          signature(object = "ClusterSet"),
          function(object,
                   clusters = "all",
                   verbose = TRUE) {
            
            if (length(clusters) == 1){
              if (clusters == "all"){
                clusters <- unique(object@cluster)
              }
            }
            
            
            for (cur_cluster in clusters) {
              # Check if the current cluster id provided exists
              if (!(cur_cluster %in% unique(object@cluster))) {
                stop(paste0("Cluster ", cur_cluster, " doesn't exist."))
              }
              
              # Check if there is a result provided by enrich_analysis function for the current cluster
              if(is.null(object@cluster_annotations[[cur_cluster]]$result)){
                print_msg(msg_type = "WARNING",
                          msg = paste0("No functional enrichment analysis results for cluster ", cur_cluster, ".")) #Continue through the next cluster without plotting
              } else {
                
                # Print an informative message to announce plot of the results of the current cluster
                if (verbose) {
                  print_msg(msg_type = "INFO",
                            msg = paste0("Plot enrichment analysis results for cluster ", cur_cluster))
                }
                
                # Create a plotly result plot
                plot_enrich <- gostplot(object@cluster_annotations[[cur_cluster]],
                                        interactive = TRUE)
                plot_enrich <- plot_enrich %>% plotly::layout(title = paste0("Cluster ", cur_cluster),
                                                              xaxis = list(title = 'Database'))
                
                # Store the plot in the cluster_annotation slot
                object@cluster_annotations[[cur_cluster]]$plot <- plot_enrich
                print(plot_enrich)
              }
            }
            
            # Print a message to inform which slot contains the results
            if (verbose){
              print_msg(msg_type = "INFO",
                        msg = paste("Plots are stored in object@cluster_annotations[[<cluster>]]$plot"))
            }
            
            return(object)
          }
)

#########################################################
##      END PACKAGE scigenex
#########################################################
