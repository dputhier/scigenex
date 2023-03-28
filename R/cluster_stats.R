
setGeneric("cluster_stats", 
           function(object)
             standardGeneric("cluster_stats")
)

setMethod(
  "cluster_stats", signature("ClusterSet"),
  function(object) {
    #object <- scigenex::check_format_cluster_set(object)
    
    df <- data.frame(size=clust_size(object), 
                     row.names = names(clust_size(object)))
    
    gene_clust <- as.factor(gene_cluster(object))
    df_split <- split(as.data.frame(object@data), gene_clust)
    
    tmp <- unlist(lapply(df_split, sum))
    df$sum_count <- tmp
    
    tmp <- unlist(lapply(df_split, var))
    df$var_total <- tmp
    
    tmp <- unlist(lapply(df_split, sd))
    df$sd_total <- tmp
    
    tmp <- unlist(lapply(df_split, sd))
    df$sd_total <- tmp
    
    tmp <- unlist(lapply(df_split, sd))
    df$sd_total <- tmp
    
    all_dot_prod <- vector()
    
    for (i in 1:length(object@gene_clusters)) {
      
      print_msg(paste0("Computing dot product for cluster: ", i), 
                msg_type = "DEBUG")
      cur_clust <- object@data[object@gene_clusters[[i]],]
      cur_clust[cur_clust >= 1] <- 1
      cur_clust[cur_clust < 1] <- 0
      cur_dot_prod <- cur_clust %*% t(cur_clust)
      diag(cur_dot_prod) <- NA
      cur_dot_prod_median_of_max <-median(apply(cur_dot_prod, 
                                                1, 
                                                max, 
                                                na.rm = T))
      all_dot_prod[i] <- cur_dot_prod_median_of_max
      
    }
    
    df$dot_prod <- all_dot_prod
  
    return(df)
    
}
)

plot_cluster_stats <- function(x, highlight=NULL){

  if(!is.null(highlight)){
    if(length(highlight) != nrow(x))
      print_msg("The number of clusters to highlight should have same length as nrow(x).",
                msg_type = "STOP")
    
    highlight <- as.factor(highlight)
    
    if(length(levels(highlight)) != 2)
      print_msg("The 'highlight' vector should have two modalities.",
                msg_type = "STOP")
  }

  x <- reshape2::melt(as.matrix(x))
  colnames(x) <- c("cluster", "stat", "value")
  
  if(!is.null(highlight))
    x$highlight <- highlight 
  
  x$cluster <- factor(x$cluster, 
                         levels=sort(unique(x$cluster)),
                         ordered=T)
  x$stat <- factor(x$stat, 
                      levels=sort(unique(x$stat)),
                      ordered=T)
  if(!is.null(highlight)){
    ggplot2::ggplot(data=x, 
                    mapping=ggplot2::aes(x=cluster,
                                         y=value,
                                         fill=highlight)) +
            ggplot2::geom_col() + 
            ggplot2::facet_grid(~stat, scales = "free") +
            coord_flip() +
            theme_bw()
  }else{
    ggplot2::ggplot(data=x, 
                    mapping=ggplot2::aes(x=cluster,
                                         y=value)) +
            ggplot2::geom_col() + 
            ggplot2::facet_grid(~stat, scales = "free") +
            coord_flip()+
            theme_bw()
  }
}
