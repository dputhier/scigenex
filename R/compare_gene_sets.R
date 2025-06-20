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
#' @examples
#' set.seed(123)
#' set_1 <- list(letters[1:10], letters[11:20])
#' x <- sample(letters[1:20])
#' set_2 <- list(x[1:5], x[6:20])
#' comp <- compare_genesets(set_1, set_2, stat = "jaccard")
#' comp <- compare_genesets(set_1, set_2, stat = "hypergeom")
#' @export
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
  
  for(i in set_1){
    if(length(i)==0)
      print_msg("Empty sets not allowed.", msg_type = "STOP")
  }
  
  for(i in set_2){
    if(length(i)==0)
      print_msg("Empty sets not allowed.", msg_type = "STOP")
  }
  
  if(is.null(background)){
    print_msg("Computing background from the union of set_1 and set_2.",
              msg_type = "INFO")
    background <- unique(c(unlist(set_1), unlist(set_2)))
  }else{
    background <- unique(background)

  }

  print_msg(paste0("Background size : ", length(background)),
            msg_type = "INFO")
  
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
#' @export
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
  
  if(stat == "hypergeom" | layout == "square"){
    print_msg("Ceiling pvalue to 1e-320 (R limit).")
    res[res <= 1e-320] <- 1e-320
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
      geom_tile(color="white", linewidth=0) +
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
    
    # Some 'magick' (or not...) for axes 

    v1 <- length(unique(res_melt$label_set_1))
    v2 <- length(unique(res_melt$label_set_2))
    add_axe <- ifelse( v1 == v2,
                      0,
                      0.5)
    
    # Create the diagram
    ggplot(res_melt, mapping=aes(x=Set_1, 
                                 y=Set_2, 
                                 fill=stat,
                                 width=jaccard,
                                 height=jaccard)) + 
      geom_tile() +
      theme_bw() +
      scale_size_area("Jaccard", max_size = 1, guide = "legend") +
      scale_fill_gradientn(legend_label, colours = colors) +
      geom_vline(xintercept = seq(0.5, length(set_1), by=1), 
                 col="black",
                 linewidth=0.3) + 
      geom_hline(yintercept = seq(0.5, length(set_2), by=1), 
                 col="black",
                 linewidth=0.3) + 
      theme(panel.grid = element_blank(),
            axis.text.x = element_text(angle=45, vjust = 0.9, hjust=1)) +
      scale_x_continuous(expand=expansion(mult = c(0, 0), 
                                             add = c(0, add_axe)), 
                         breaks = res_melt$Set_1, 
                         labels = res_melt$label_set_1) +
      scale_y_continuous(expand=expansion(mult = c(0, 0), 
                                          add = c(0, add_axe)), 
                         breaks = res_melt$Set_2, 
                         labels = res_melt$label_set_2) +
      coord_equal 
  }

}



#################################################################
##    Map a user list to a ClusterSet and
##    Display various metrics related to overlap 
##    with internaly stored clusters   
#################################################################
#' @title Map a user list to a ClusterSet and display various metrics  related to overlap
#' with internaly stored clusters 
#' @description
#' Map a user list to a ClusterSet and display various metrics  related to overlap
#' with internaly stored clusters 
#' @param object A ClusterSet object
#' @param user_list A set of elements to be compared to clusters stored in the object.
#' @param name_user_list A name for the user list.
#' @param background A background(e.g. all the gene of the genome).
#' @param colors A set of colors.
#' @param module The analyses to be performed.
#' @param as.list Whether a list of plot should be returned (default to a patchwork).
#' @return A patchwork representing the various metrics (or a list of diagram). 
#' @importFrom ggplot2 theme_minimal theme element_blank geom_col aes coord_flip 
#' @importFrom ggplot2 scale_fill_manual scale_y_continuous geom_bar ggtitle coord_polar
#' @importFrom scales percent
#'
#' @examples
#' set.seed(124)
#' load_example_dataset("7871581/files/pbmc3k_medium_clusters")
#' user_list <- sample(row_names(pbmc3k_medium_clusters), 100)
#' cmp_to_a_list(pbmc3k_medium_clusters, user_list, background=user_list)
#' @export
cmp_to_a_list <- function(object=NULL, 
                          user_list=NULL,
                          name_user_list="user list",
                          background=NULL,
                          colors=c("#ED7931", "#3C78AE"),
                          module=c("clust_size", 
                                   "overlap_size",
                                   "percent_overlap",
                                   "percent_covered",
                                   "jacard_stat",
                                   "hypergeom"),
                         as.list=FALSE){
      
      check_format_cluster_set(object)
      
      if(is.null(user_list))
        print_msg('Please provide a lists as input.', msg_type = "STOP")
      
      if(is.null(background))
        print_msg("Please provide background genes (e.g all the known genes).",
                  msg_type = "STOP")
    
      
      print_msg(paste0("Background size : ", length(background)),
                msg_type = "INFO")
    
      plot_list <- list()
      
      theming <- theme_minimal() +
                 theme(panel.grid.minor = element_blank()) 
      
      if("clust_size" %in% module){

        print_msg("Computing cluster sizes.", msg_type = "INFO")
    
        tmp <- clust_size(object)
        df <- data.frame(gene_set=c(name_user_list, names(tmp)), 
                         size=c(length(user_list), tmp),
                         source=c(name_user_list, 
                                  rep("ClusterSet", nclust(object))))
        
        df$gene_set <- factor(df$gene_set, levels=rev(c(name_user_list, 
                                                        clust_names(object))),
                              ordered = TRUE)

        gene_set <- size <- source <- NULL
        p1 <- ggplot(df, aes(x=gene_set, 
                             y=size, 
                             fill = source)) +
          geom_col(color="white") +
          scale_fill_manual(values=colors) + 
          coord_flip() +
          ylab("Number of\ngenes") +
          xlab("Gene set") +
          theming
       
        plot_list[["clust_size"]] <-  p1
      }
      
      if("overlap_size" %in% module){
        
        print_msg("Computing 'overlap_size'.")
        
        print_msg("Computing overlapping sizes.", msg_type = "INFO")
        
        intersection <- compare_genesets(object@gene_clusters, list(user_list), stat = "intersection")[,1]
        
        df <- data.frame(gene_set=names(intersection), 
                         intersection=intersection)
        
        df$gene_set <- factor(df$gene_set, levels=rev(clust_names(object)),
                              ordered = TRUE)
    
        p2 <- ggplot(df, aes(x=gene_set, 
                             y=intersection)) +
          geom_col(fill=colors[1], color="white") +
          scale_fill_manual(values=colors) + 
          coord_flip() +
          ylab(paste0("Number of genes\nintersecting with\n", name_user_list)) +
          xlab("") +
          theming
        
        plot_list[["overlap_size"]] <-  p2
      }
      
      if("percent_overlap" %in% module){
        
        print_msg("Computing 'percent_overlap'.")
        
        tmp <- gene_cluster(object)
        
        df <- data.frame(gene_set=tmp, 
                         gene=names(tmp),
                         overlap=ifelse(names(tmp) %in% user_list, 
                                        name_user_list, 
                                        paste0("Not in ", name_user_list)))
        
        df$gene_set <- factor(df$gene_set, 
                              levels=rev(clust_names(object)),
                              ordered = TRUE)
        
        this_col <- colors[1:2]
        
        names(this_col) <- c(name_user_list, paste0("Not in ", name_user_list))
        
        p3 <- ggplot(df, aes(x=gene_set, 
                             fill=overlap)) +
          geom_bar(aes(y = (..count..)/sum(..count..)), color="white", position="fill") + 
          scale_y_continuous(labels=scales::percent) +
          scale_fill_manual(values=this_col) + 
          coord_flip() +
          ylab(paste0("Fraction of Gene Set covered \nby ",  name_user_list)) +
          xlab("") +
          theming
        
        plot_list[["percent_overlap"]] <-  p3
      }
      
      if("percent_covered" %in% module){
        
        names(this_col) <- c("In gene set", "Not in gene set" )
        
        print_msg("Computing 'percent_covered'.")
        
        intersection <- compare_genesets(object@gene_clusters, list(user_list), 
                                         stat = "intersection", 
                                         background = background)[,1]

        df_1 <- data.frame(gene_set=names(intersection), 
                         intersection=intersection/length(user_list) * 100)
        
        df_1$source <- "In gene set" 
        
        df_2 <- df_1
        
        df_2$intersection <- 100 - df_1$intersection
        
        df_2$source <-  "Not in gene set" 
        
        
        df_3 <- rbind(df_2, df_1)
    
        df_3$source <- factor(df_3$source, 
                              levels=rev(c("In gene set", "Not in gene set" )), 
                              ordered = TRUE)
        
        df_3$gene_set <- factor(df_3$gene_set, 
                                levels=rev(clust_names(object)),
                                ordered = TRUE)
        
        p4 <- ggplot(df_3, aes(x=gene_set, 
                               y=intersection,
                               fill=source)) +
          geom_col(color="white", position = "stack") + 
          coord_flip() +
          scale_fill_manual(values=this_col) +
          ylab(paste0("Fraction of\n",  name_user_list, " covered")) +
          xlab("") +
          theming
        
        plot_list[["percent_covered"]] <-  p4
      }
      
      if("jacard_stat" %in% module){
        
        print_msg("Computing 'jacard_stat'.")
        
        print_msg("Computing overlapping sizes.", msg_type = "INFO")
        
        jaccard <- compare_genesets(object@gene_clusters, list(user_list), stat = "jaccard", background=background)[,1]
        
        df <- data.frame(gene_set=names(jaccard), 
                         jaccard=jaccard)
        
        df$gene_set <- factor(df$gene_set, 
                              levels=rev(clust_names(object)),
                              ordered = TRUE)

        p5 <- ggplot(df, aes(x=gene_set, 
                             y=jaccard)) +
          geom_col(fill=colors[1], color="white") +
          scale_fill_manual(values=colors) + 
          ggtitle("Jaccard value") +
          xlab("Clusters") +
          ylab("") +
          coord_polar() +
          theming 
        
        
        plot_list[["jacard_stat"]] <-  p5
      }
      
      if("hypergeom" %in% module){
        
        print_msg("Computing 'hypergeometric'.")
        
        print_msg("Computing overlapping sizes.", msg_type = "INFO")
        
        hypergeometric <- compare_genesets(object@gene_clusters, list(user_list), stat = "hypergeom", background=background)[,1]
        hypergeometric <- -log10(hypergeometric)
        
        df <- data.frame(gene_set=names(hypergeometric), 
                         jaccard=hypergeometric)
        
        df$gene_set <- factor(df$gene_set, 
                              levels=rev(clust_names(object)),
                              ordered = TRUE)
        
        p6 <- ggplot(df, aes(x=gene_set, 
                             y=hypergeometric)) +
          geom_col(fill=colors[1], color="white") +
          scale_fill_manual(values=colors) + 
          ggtitle("Hypergeometric test\n-log10(p-value)") +
          xlab("Clusters") +
          ylab("") +
          coord_polar() + 
          theming 
        
        plot_list[["hypergeom"]] <-  p6
      }
      
      if(as.list){
        return(plot_list)
      }else{
        return((p1 | p2 | p3 ) / (p4 | p5 | p6 ))
      }

}

  
    
    
    
