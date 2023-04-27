#################################################################
##    Compare two lists of gene sets
#################################################################

# An internal function to compute p-value
# of hypergeometric test
hyper_fun <- function(a, b, N) {

  a <- unique(a[a %in% N])
  b <- unique(b[b %in% N])
  N <- length(unique(N))
  x <- intersect(a,b)
  q <- ifelse(length(x)-1 < 0, 0, length(x)-1)
  m <- length(a)
  k <- length(b)
  n <- N - m 

  print_msg(paste0("Calling phyper(q=", 
                   q, 
                   ", m=", 
                   m, 
                   ", n=",
                   n,
                   ", k=",
                   k,
                   ", lower.tail = FALSE)"), msg_type = "DEBUG")
  
  pval <- stats::phyper(q=q, 
                 m=m, 
                 n=n, 
                 k=k, 
                 lower.tail = FALSE)
  return(pval)
}

vhyper <- Vectorize(hyper_fun)


# An internal function to compute jaccard
# intersectio,, union...
check_inter <- function(a, b, 
                        stat=c("union", 
                               "intersection", 
                               "jaccard",
                               "size_set_1",
                               "size_set_2",
                               "diff_set_1",
                               "diff_set_2"),
                        background=NULL) {
  
  a <- unique(a[a %in% background])
  b <- unique(b[b %in% background])
  
  intersection <- length(intersect(a, b))
  union <- length(unique(c(a, b)))
  
  # from: https://tinyurl.com/h6dsuc26
  if(stat=="intersection"){
    return(intersection)
  }else if(stat=="union"){
    return(union)
  }else if(stat=="jaccard"){
    return (intersection/union)
  }else if(stat=="size_set_1"){
    return(length(a))
  }else if(stat=="size_set_2"){
    return(length(b))
  }else if(stat=="diff_set_1"){
    return(length(a) - intersection)
  }else if(stat=="diff_set_2"){
    return(length(b) - intersection)
  }
}

vcheck_inter <- Vectorize(check_inter)

#' @title Compare two lists of gene sets
#' @description  
#' This function compares two lists of gene sets using either the Jaccard index or the hypergeometric test.
#'
#' @param set_1 A list containing gene sets to be compared.
#' @param set_2 A list containing gene sets to be compared.
#' @param stat The statistics to be computed between gene sets. It can be either "jaccard", "hypergeom", "intersection"
#' "size_set_1", "size_set_2", "diff_set_1" (specific to set_1), "diff_set_2" (specific to set_2). The background
#' is taken into account. Note that hypergeometric tests check for enrichment.
#' @param background The background (universe) to consider. Default to the non-redundant list of elements 
#' merged from set_1 and set2. You may provide a vector with all genes of the genome for instance. 
#' @return A matrix of comparison results where each row corresponds to a gene set in set_1, 
#' and each column corresponds to a gene set in set_2.
#' @details The Jaccard index is a measure of similarity between two sets defined as the size of the 
#' intersection divided by the size of the union of the sets. The hypergeometric test is used to 
#' determine whether the overlap between two sets is more significant than expected by chance. The 
#' 'intersection' method, simply computes the size of the intersection between a and b. The "union", 
#' "size_set_1", "size_set_2", "diff_set_1" and "diff_set_2" compute the union of the two sets, the size 
#' of gene sets from set_1, the size of gene sets from set_2, the gene that are specific to set_1,  the 
#' gene that are specific to set_2, respectively. 
#' @export compare_genesets
#' @examples
#' set.seed(123)
#' set_1 <- list(letters[1:10], letters[11:20])
#' x <- sample(letters[1:20])
#' set_2 <- list(x[1:5], x[6:20])
#' comp <- compare_genesets(set_1, set_2, stat = "jaccard")
#' comp <- compare_genesets(set_1, set_2, stat = "hypergeom")
#' 
compare_genesets <- function(set_1=NULL, 
                             set_2=NULL,
                             stat=c("jaccard",
                                    "hypergeom",
                                    "intersection",
                                    "union",
                                    "size_set_1",
                                    "size_set_2",
                                    "diff_set_1",
                                    "diff_set_2"),
                             background=NULL){

  if(is.null(set_1) | is.null(set_2))
    print_msg('Please provide two lists as input.', msg_type = "STOP")
  
  if(!is.list(set_1) | !is.list(set_2))
    print_msg('Please provide two lists as input.', msg_type = "STOP")
  
  stat <- match.arg(stat)
  
  if(is.null(background)){
    print_msg("Computing background from the union of set_1 and set_2.",
              msg_type = "INFO")
    background <- unique(c(unlist(set_1), unlist(set_2)))
  }else{
    background <- unique(background)
  }

  if(stat != "hypergeom"){
    res <- outer(X=set_1, 
                 Y=set_2, 
                 FUN = "vcheck_inter", 
                 stat=stat, background=list(background))

  }else if(stat == "hypergeom"){

    res <- outer(X=set_1, Y=set_2, FUN = "vhyper", N=list(background))
    
  }

  if(is.null(names(set_1))){
    rownames(res) <- paste0("Set_1_", 1:length(set_1))
  }else{
    rownames(res) <- names(set_1)
  }
  
  if(is.null(names(set_2))){
    colnames(res) <- paste0("Set_2_", 1:length(set_2))
  }else{
    colnames(res) <- names(set_2)
  }
  
  return(res)
  
} 

#################################################################
##    Plot the result of compare_genesets
#################################################################
#' @title Plot the result of compare_genesets
#'
#' @description
#' This function plots the results of compare_genesets.
#'
#' @param set_1 A list containing gene sets to be compared.
#' @param set_2 A list containing gene sets to be compared.
#' @param stat The statistics to be computed between gene sets. It can be either "jaccard", "hypergeom", "intersection"
#' "size_set_1", "size_set_2", "diff_set_1" (specific to set_1), "diff_set_2" (specific to set_2). The background
#' is taken into account. Note that hypergeometric tests check for enrichment.
#' @param transform The transformation to be applied to the values. It can be either "NA", "log10", or "log2", "-log10", "-log2.
#' @param colors The color palette to be used in the plot.
#' @param layout The type of diagram. Either "raster" (a scatter plot showing the statistics of interest) or "square".
#' The "square" layout shows the hypergeometric pvalue (color) and the Jaccard result (size of the square). 
#' @param background The background (universe) to consider. Default to the non-redundant list of elements 
#' merged from set_1 and set2. You may provide a vector with all genes of the genome for instance. 
#' @param coord_equal make sure that the an equal length on both axis represents the same change in units
#' @details see compare_genesets. 
#' @importFrom ggplot2 geom_hline scale_x_continuous expansion scale_y_continuous scale_size_area coord_equal
#' @return A ggplot object representing the comparison results.
#'
#' @export plot_cmp_genesets
#'
#' @examples
#' set.seed(123)
#' set_1 <- list(letters[1:10], letters[11:20], letters[21:30])
#' x <- sample(letters[1:30])
#' set_2 <- list(x[1:5], x[6:20], letters[21:30])
#' res <- compare_genesets(set_1, set_2, stat = "jaccard")
#' plot_cmp_genesets(set_1, set_2, stat = "jaccard")
#' plot_cmp_genesets(set_1, set_2, stat = "hypergeom", transform = "log10")
#' plot_cmp_genesets(set_1, set_2, layout="square", transform = "-log10")
#'
plot_cmp_genesets <- function(set_1=NULL, 
                              set_2=NULL,
                              stat=c("jaccard",
                                     "hypergeom",
                                     "intersection",
                                     "union",
                                     "size_set_1",
                                     "size_set_2",
                                     "diff_set_1",
                                     "diff_set_2"),
                              transform=c("None", "log10", "log2", "-log10", "-log2"),
                              colors=colors_for_gradient("Ju1"),
                              layout=c("raster", "square"),
                              background=NULL,
                              coord_equal=TRUE){

  if(is.null(set_1) | is.null(set_2))
    print_msg('Please provide two lists as input.', msg_type = "STOP")
  
  if(!is.list(set_1) | !is.list(set_2))
    print_msg('Please provide two lists as input.', msg_type = "STOP")
  
  stat <- match.arg(stat)
  layout <- match.arg(layout)
  transform <- match.arg(transform)
  
  print_msg(paste0("Using transformation: ", transform),
            msg_type="INFO")
  
  if(layout == "raster"){
    res <- compare_genesets(set_1=set_1, 
                          set_2=set_2, 
                          stat=stat)
  }else{
    res <- compare_genesets(set_1=set_1, 
                            set_2=set_2, 
                            stat="hypergeom")
  }
  
  res_melt <- reshape2::melt(as.matrix(res))
  colnames(res_melt) <- c("Set_1", "Set_2", "stat")
  
  if(transform != "None"){
    res_melt$stat <- do.call(gsub("-","", transform), list(res_melt$stat))
    if(transform %in% c("-log10", "-log2")){
      res_melt$stat <- -res_melt$stat
    }
  }
  
  if(coord_equal){
    coord_equal <- ggplot2::coord_equal()
  }else{
    coord_equal <- NULL
  }

  if(layout == "raster"){
    ggplot(res_melt, mapping=aes(x=Set_1, 
                                 y=Set_2, 
                                 fill=stat)) + 
      geom_tile(color="white", linewidth=2) +
      theme_bw() +
      theme(axis.text.x = element_text(angle=45, vjust = 0.5)) +
      scale_fill_gradientn(colours = colors) +
      coord_equal
  }else if(layout == "square"){
    
    if(is.na(transform)){
      legend_label <- "p-value"
    }else{
      legend_label <- paste0(transform, "(p-value)")
    }
    inter <-compare_genesets(set_1=set_1, 
                             set_2=set_2, 
                             stat="jaccard")

    inter_melt <- reshape2::melt(as.matrix(inter))    
    res_melt$jaccard <- inter_melt$value
    res_melt$jaccard_2 <- inter_melt$value
    
    res_melt$label_set_1 <- as.character(res_melt$Set_1)
    res_melt$label_set_2 <- as.character(res_melt$Set_2)
    res_melt$Set_1 <- as.numeric(as.factor(res_melt$Set_1))
    res_melt$Set_2 <- as.numeric(as.factor(res_melt$Set_2))
    
    ggplot(res_melt, mapping=aes(x=Set_1, 
                                 y=Set_2, 
                                 fill=stat,
                                 width=jaccard,
                                 height=jaccard)) + 
      geom_vline(xintercept = seq(0.5, length(set_1), by=1), 
                 col="black",
                 linewidth=0.3) + 
      geom_hline(yintercept = seq(0.5, length(set_2), by=1), 
                 col="black",
                 linewidth=0.3) + 
      geom_tile() +
      theme_bw() +
      theme(panel.grid = element_blank(),
            axis.text.x = element_text(angle=45, vjust = 0.9, hjust=1)) +
      scale_x_continuous(expand=expansion(mult = c(0, 0), 
                                             add = c(0, 0.5)), 
                         breaks = res_melt$Set_1, 
                         labels = res_melt$label_set_1) +
      scale_y_continuous(expand=expansion(mult = c(0, 0), 
                                          add = c(0, 0.5)), 
                         breaks = res_melt$Set_2, 
                         labels = res_melt$label_set_2) +
      scale_size_area("Jaccard", max_size = 1, guide = "legend") +
      scale_fill_gradientn(legend_label, colours = colors, ) +
      coord_equal 
  }

}
