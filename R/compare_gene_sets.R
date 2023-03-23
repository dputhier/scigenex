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
  print(N)
  print_msg(paste0("Calling phyper(q=", 
                   q, 
                   ", m=", 
                   m, 
                   ", n=",
                   n,
                   ", k=",
                   k,
                   ", lower.tail = FALSE)"), msg_type = "DEBUG")
  
  pval <- phyper(q=q, 
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

#' Compare two lists of gene sets
#'
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
#' @export
#' @examples
#' set.seed(123)
#' set_1 <- list(letters[1:10], letters[11:20])
#' x <- sample(letters[1:20])
#' set_2 <- list(x[1:5], x[6:20])
#' compare_genesets(set_1, set_2, stat = "jaccard")
#' compare_genesets(set_1, set_2, stat = "hypergeom")
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
    colnames(res) <- paste0("Set_2_", 1:length(set_1))
  }else{
    colnames(res) <- names(set_2)
  }
  
  return(res)
  
} 

#################################################################
##    Plot the result of compare_genesets
#################################################################
# Plot the result of compare_genesets

plot_cmp_geneset <- function(set_1=NULL, 
                             set_2=NULL,
                             method=c("Jaccard",
                                      "hypergeom"),
                             transform=c(NA,"log10", "log2"),
                             colors=colors_for_gradient("Ju1")){
  if(is.null(set_1) | is.null(set_2))
    print_msg('Please provide two lists as input.', msg_type = "STOP")
  
  if(!is.list(set_1) | !is.list(set_2))
    print_msg('Please provide two lists as input.', msg_type = "STOP")
  
  method <- match.arg(method)
  
  
  transform <- match.arg(transform)
  
  print(paste0("Using transformation: ", transform))
  
  res <- compare_genesets(set_1=set_1, 
                          set_2=set_2, 
                          method=method)
  
  res_melt <- reshape2::melt(as.matrix(res))
  colnames(res_melt) <- c("Set_1", "Set_2", "value")
  
  if(!is.na(transform))
    res_melt$value <- do.call(transform, list(res_melt$value))
  
  ggplot(res_melt, mapping=aes(x=Set_1, 
                               y=Set_2, 
                               fill=value)) + 
    geom_tile(color="white", size=2) +
    theme_bw() +
    theme(axis.text.x = element_text(angle=45, vjust = 0.5)) +
    scale_fill_gradientn(colours = colors)
}
