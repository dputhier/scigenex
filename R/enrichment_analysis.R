#################################################################
##    Running gprofiler2 to perform enrichment analysis
#################################################################

#' @title
#' Enrichment analysis using GO database of gene clusters from a ClusterSet object 
##' @description
#' This function performs enrichment analysis of all gene clusters from a ClusterSet object
#'  and store the results in the cluster_annotations slot.
#' @param object an object of class \code{ClusterSet}.
#' @param species a character string indicating the species name,
#'  as a concatenation of the first letter of the name (uppercase) and the family name (lowercase),
#'   e.g human -> Hsapiens, mouse -> Mmusculus
#' @param ontology a character string indicating the ontology to use for the enrichment analysis. One of "BP", "MF", and "CC" ontology, or "ALL".
#'
#' @return A \code{ClusterSet} object
#' @method enrich_go ClusterSet
#' @export enrich_go
#' @keywords internal
setGeneric("enrich_go",
           function(object,
                    species="Hsapiens",
                    ontology="all") {
             standardGeneric("enrich_go")
           })

#' @title
#' Enrichment analysis using GO database of gene clusters from a ClusterSet object 
##' @description
#' This function performs enrichment analysis of all gene clusters from a ClusterSet object
#'  and store the results in the cluster_annotations slot.
#' @param object an object of class \code{ClusterSet}.
#' @param species a character string indicating the species name,
#'  as a concatenation of the first letter of the name (uppercase) and the family name (lowercase),
#'   e.g human -> Hsapiens, mouse -> Mmusculus
#' @param ontology a character string indicating the ontology to use for the enrichment analysis. One of "BP", "MF", and "CC" ontology, or "ALL".
#'
#' @return A \code{ClusterSet} object
#' @method enrich_go ClusterSet
#' @export enrich_go
#'
#' @examples
#' load_example_dataset("7871581/files/pbmc3k_medium_clusters")
#' load_example_dataset("8028126/files/pbmc3k_medium_clusters_enr")
#' # pbmc3k_medium_clusters_enr <- enrich_go(pbmc3k_medium_clusters)
#' pbmc3k_medium_clusters_enr
#' is(pbmc3k_medium_clusters_enr@gene_cluster_annotations)
setMethod("enrich_go",
          signature(object = "ClusterSet"),
          function(object,
                   species="Hsapiens",
                   ontology="all") {
            
            ## Check format object arg
            check_format_cluster_set(object)
            
            if (!species %in% c("Hsapiens", "Mmusculus")){
              print_msg("This species name is not supported.", 
                        msg_type = "STOP")}
            
            if (species == "Hsapiens") {
              org_db <- org.Hs.eg.db::org.Hs.eg.db
              go_name <- "org.Hs.eg.db"
              print_msg(msg_type = "INFO",
                        msg = "Species used : Homo sapiens")
            }
            
            if (species == "Mmusculus") {
              org_db <- org.Mm.eg.db::org.Mm.eg.db
              go_name <- "org.Mm.eg.db"
              print_msg(msg_type = "INFO",
                        msg = "Species used : Mus musculus")
            }
            
            
            for(cluster in unique(names(object@gene_clusters))){
              print_msg(paste0("Enrichment analysis for cluster ", cluster),
                        msg_type = "INFO")
              cluster_name = paste0("Cluster_", cluster)
              query = object@gene_clusters[[cluster]]
              
              # Convert gene id to EntrezId format
              suppressMessages(query_entrezid <- AnnotationDbi::select(org_db, 
                                                                       keys = query,
                                                                       columns = c("ENTREZID", "SYMBOL"),
                                                                       keytype = "SYMBOL"))
              
              print_msg(msg_type = "DEBUG",
                        msg = paste0("Cluster ",
                                     cluster,
                                     "\n",
                                     "\t|--Number of gene names converted to EntrezId : ", 
                                     length(which(!is.na(query_entrezid[, "ENTREZID"]))),
                                     "\n", 
                                     "\t|--Number of gene names not converted to EntrezId : ", 
                                     length(which(is.na(query_entrezid[, "ENTREZID"]))),
                                     "\n\t|-- Calling enrichGO."))
              
              
              # Enrichment analysis using GO database
              enrich_go_res <- clusterProfiler::enrichGO(query_entrezid[,"ENTREZID"],
                                                         OrgDb = go_name,
                                                         ont = ontology,
                                                         readable = TRUE)
              
              # Store results in the ClusterSet object
              object@gene_cluster_annotations[[cluster]] <- enrich_go_res
            }
            return(object)
          }
)


############################################################################
##    Visualization of enrichment analyses results from a ClusterSet object
############################################################################
#' @title
#' Visualization of enrichment analyses results from a ClusterSet object
##' @description
#' Display enrichment analyses results for a ClusterSet object.
#' @param object A \code{ClusterSet} object.
#' @param clusters  A vector of cluster id to plot.
#' @param type A vector providing the type of plot.
#' @param nb_terms An integer indicating the number of terms in the plot.
#' @param terms_size An interger indicating the wrap length of terms
#' @param font_size The font size for the x axis.
#'
#' @return A \code{ClusterSet} object
#' @method viz_enrich ClusterSet
#' @export viz_enrich
#' @keywords internal
setGeneric("viz_enrich",
           function(object,
                    clusters = "all",
                    type = "dotplot",
                    nb_terms = 20,
                    terms_size = 50,
                    font_size=4) {
             standardGeneric("viz_enrich")
           })


#' @title
#' Visualization of enrichment analyses results from a ClusterSet object
##' @description
#' Display enrichment analyses results for a ClusterSet object.
#' @param object A \code{ClusterSet} object.
#' @param clusters  A vector of cluster id to plot.
#' @param type A vector providing the type of plot.
#' @param nb_terms An integer indicating the number of terms in the plot.
#' @param terms_size An interger indicating the wrap length of terms
#' @param font_size The font size for the x axis.
#'
#' @return A \code{ClusterSet} object
#' @method viz_enrich ClusterSet
#' @export viz_enrich
#'
#' @examples
#' load_example_dataset("8028226/files/pbmc3k_medium_clusters_enr_sub")
#' viz_enrich(pbmc3k_medium_clusters_enr_sub, cluster = 1, nb_terms = 5)
setMethod("viz_enrich",
          signature(object = "ClusterSet"),
          function(object,
                   clusters = "all",
                   type = c("dotplot", "barplot"),
                   nb_terms = 20,
                   terms_size = 50,
                   font_size=4) {
            
            ## Check format object arg
            check_format_cluster_set(object)
            
            if (length(clusters) == 1){
              if (clusters == "all"){
                clusters <- object@gene_clusters_metadata$cluster_id
              }
            }
            
            list_plot <- list()
            list_dotplot <- list()
            list_barplot <- list()
            
            for (cur_cluster in clusters) {
              # Check if the current cluster id provided exists
              if (!(cur_cluster %in% object@gene_clusters_metadata$cluster_id)) {
                stop(paste0("Cluster ", cur_cluster, " doesn't exist."))
              }
              
              # Check if there is a result provided by enrich_go function for the current cluster
              if(nrow(object@gene_cluster_annotations[[cur_cluster]]@result) == 0){
                print_msg(msg_type = "WARNING",
                          msg = paste0("No functional enrichment analysis results for cluster ", cur_cluster, ".")) #Continue through the next cluster without plotting
              } else {
                
                # Print an informative message to announce plot of the results of the current cluster
                print_msg(msg_type = "INFO",
                          msg = paste0("Plot enrichment analysis results for cluster ", cur_cluster))
                
                
                
                # Create a ggplot - dotplot
                if ("dotplot" %in% type){
                  if(object@gene_cluster_annotations[[cur_cluster]]@ontology == "GOALL"){
                    dot_plot <- enrichplot::dotplot(object@gene_cluster_annotations[[cur_cluster]],
                                                    split="ONTOLOGY",
                                                    showCategory=nb_terms,
                                                    label_format = terms_size)
                    dot_plot <- dot_plot + facet_grid(ONTOLOGY~., scales="free")
                    dot_plot <- dot_plot + ggtitle(paste0("Gene cluster ", cur_cluster))
                    dot_plot <- dot_plot + theme(axis.text = element_text(size=font_size))
                  } else {
                    dot_plot <- enrichplot::dotplot(object@gene_cluster_annotations[[cur_cluster]],
                                                    showCategory=nb_terms,
                                                    label_format = terms_size)
                    dot_plot <- dot_plot + ggtitle(paste0("Gene cluster ", cur_cluster))
                    dot_plot <- dot_plot + theme(axis.text = element_text(size=font_size))
                  }
                  list_dotplot <- append(list_dotplot, list(dot_plot))
                }
                
                # Create a ggplot - barplot
                if ("barplot" %in% type){
                  if(object@gene_cluster_annotations[[cur_cluster]]@ontology == "GOALL"){
                    bar_plot <- graphics::barplot(object@gene_cluster_annotations[[cur_cluster]],
                                                  split="ONTOLOGY",
                                                  showCategory=nb_terms,
                                                  label_format = terms_size)
                    bar_plot <- bar_plot + facet_grid(ONTOLOGY~., scales="free")
                    bar_plot <- bar_plot + ggtitle(paste0("Gene cluster ", cur_cluster)) 
                  } else {
                    bar_plot <- graphics::barplot(object@gene_cluster_annotations[[cur_cluster]],
                                                  showCategory=nb_terms,
                                                  label_format = terms_size)
                    bar_plot <- bar_plot + ggtitle(paste0("Gene cluster ", cur_cluster))
                  }
                  list_barplot <- append(list_barplot, list(bar_plot))
                }
                list_plot <- append(list_barplot, list_dotplot)
              }
            }
            
            return(invisible(lapply(list_plot, print)))
          }
)


############################################################################
##    plot_clust_enrichments()
############################################################################

#' @title Plot Cluster Enrichments
#'
#' @title Display cluster enrichments for gene GO terms.
#'
#' @param object An object of class "ClusterSet".
#' @param stat_shown A character string specifying the statistic to be shown. Must be one of "qvalue", "p.adjust", or "pvalue". Defaults to "qvalue".
#' @param nb_go_term An integer specifying the number of top GO terms to select per cluster. Defaults to 6. Note that some terms may be shared between clusters.
#' @param gradient_palette A vector of colors specifying the gradient palette for the color scale. Defaults to colors_for_gradient("Magma").
#' @param label_fun A function to customize the labels of the GO terms. Defaults to NULL.
#' @param term_order A vector specifying the desired order of the GO terms. Defaults to NULL.
#' @param return_mat Logical. Do not return a ggplot object but the prepared matrix.
#' @param floor A value for flooring the statistic.
#' 
#' @return A ggplot object displaying the cluster enrichments for the GO terms.
#' 
#' @seealso
#' \code{\link{enrich_go}}
#' \code{\link{colors_for_gradient}}
#'  
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 scale_color_gradientn
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 element_rect
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 guide_colourbar
#' @importFrom ggplot2 facet_grid
#' @importFrom reshape2 melt
#' @importFrom grid unit
#' @export plot_clust_enrichments
#' @examples
#' load_example_dataset("8028226/files/pbmc3k_medium_clusters_enr_sub")
#' plot_clust_enrichments(pbmc3k_medium_clusters_enr_sub)
#' @keywords internal
setGeneric("plot_clust_enrichments",
           function(object,
                    stat_shown=c("qvalue", "p.adjust", 
                                 "pvalue"),
                    nb_go_term=3,
                    gradient_palette=colors_for_gradient("Ju1"),
                    label_fun=NULL,
                    term_order=NULL,
                    return_mat=FALSE,
                    floor=NULL) {
             standardGeneric("plot_clust_enrichments")
           })


#' @title Plot Cluster Enrichments
#'
#' @title Display cluster enrichments for gene GO terms.
#'
#' @param object An object of class "ClusterSet".
#' @param stat_shown A character string specifying the statistic to be shown. Must be one of "qvalue", "p.adjust", or "pvalue". Defaults to "qvalue".
#' @param nb_go_term An integer specifying the number of top GO terms to select per cluster. Defaults to 6. Note that some terms may be shared between clusters.
#' @param gradient_palette A vector of colors specifying the gradient palette for the color scale. Defaults to colors_for_gradient("Magma").
#' @param label_fun A function to customize the labels of the GO terms. Defaults to NULL.
#' @param term_order A vector specifying the desired order of the GO terms. Defaults to NULL.
#' @param return_mat Logical. Do not return a ggplot object but the prepared matrix.
#' @param floor A value for flooring the statistic.
#' 
#' @return A ggplot object displaying the cluster enrichments for the GO terms.
#' 
#' @seealso
#' \code{\link{enrich_go}}
#' \code{\link{colors_for_gradient}}
#'  
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 scale_color_gradientn
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 element_rect
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 guide_colourbar
#' @importFrom ggplot2 facet_grid
#' @importFrom reshape2 melt
#' @importFrom grid unit
#' @export plot_clust_enrichments
#' @examples
#' load_example_dataset("8028226/files/pbmc3k_medium_clusters_enr_sub")
#' plot_clust_enrichments(pbmc3k_medium_clusters_enr_sub)
setMethod("plot_clust_enrichments",
          signature(object = "ClusterSet"),
          function(object,
                   stat_shown=c("qvalue", "p.adjust", 
                                "pvalue"),
                   nb_go_term=6,
                   gradient_palette=colors_for_gradient("Magma"),
                   label_fun=NULL,
                   term_order=NULL,
                   return_mat=FALSE,
                   floor=NULL) {
            
            stat_shown <- match.arg(stat_shown)
            print_msg(paste0("Using ", stat_shown, " as ordering statistic"), msg_type = "DEBUG")
            
            print_msg("Checking object format", msg_type = "DEBUG")
            check_format_cluster_set(object)
            
           
            if(is.null(object@gene_cluster_annotations) | length(object@gene_cluster_annotations)==0){
                print_msg("Please run enrich_go() first.", msg_type = "STOP")
            }
            
            print_msg("Extracting and formating data", msg_type = "DEBUG")
            print_msg(paste0("Argument nb_go_term = ", nb_go_term),  msg_type = "DEBUG")
            
            result_go <- lapply(object@gene_cluster_annotations, "slot", "result")
            
            order_mat <- function(x, column){x[order(x[,column]), ]}
            result_go <- lapply(result_go, order_mat, stat_shown)
            add_rank <- function(x, nb_go_term){x$rank <- 1:nrow(x); 
                                                x$rank[(nb_go_term + 1):nrow(x)] <- NA; return(x)}

            result_go <- lapply(result_go, add_rank, nb_go_term)
            
            for(i in 1:length(result_go)){
              result_go[[i]]$cluster <- names(result_go)[i]
            }
            
            print_msg("Merging results", msg_type = "DEBUG")
            
            m <- do.call('rbind', result_go)
          
            go_retained <- unique(m$Description[!is.na(m$rank)])
            m <- m[m$Description %in% go_retained, ]

            if(!is.null(term_order)){
              term_order_selected <- intersect(term_order, m$Description)
              m <- m[m$Description %in% term_order_selected, ]
            }
            
            m$stat <- -log10(m[, stat_shown])
            m <- m[order(m$cluster, m$rank), ]
            
            m$is_top <- m$stat
            m$is_top[is.na(m$rank)] <- "No"
            m$is_top[!is.na(m$rank)] <- "Yes"
            m$is_top <- as.factor(m$is_top)
            
            if(return_mat){
              return(m)
            }
            
            if(!is.null(label_fun)){
              m$Description <- sapply(m$Description, label_fun)
            }
            
            if(!is.null(term_order)){
              desc_levels <- sapply(term_order_selected, label_fun)
            }else{
              desc_levels <- rev(unique(m$Description))
            }

            m$Description <- factor(m$Description, levels=desc_levels, ordered = TRUE)
            m$cluster <- factor(m$cluster, 
                                levels=as.character(sort(as.numeric(unique(m$cluster)))), 
                                ordered=TRUE)
            
            cluster <- go_term <- stat <- size <- rank <- Counts <- NULL
            
            print_msg("Melting...", msg_type = "DEBUG")
            
            m_melt <- reshape2::melt(m, id.vars = c("Description", "Count", "stat", "cluster", "rank", "is_top"))
            colnames(m_melt) <- c("go_term", "Counts", "stat", "cluster", "rank", "is_top", "ID", "value")
            
            m_gg <- m_melt[order(m_melt[,"cluster"]), ]
            
            print_msg("Running ggplot...", msg_type = "DEBUG")
            
            if(!is.null(floor))
              m_gg$stat[m_gg$stat > floor] <- floor
            
            ggplot(data=m_gg, mapping=aes(x=cluster, y = go_term, color=stat, size=Counts)) +
              ggplot2::theme_bw() + 
              geom_point(shape=16) +
              ggplot2::scale_color_gradientn(colours = gradient_palette,
                                             guide = guide_colourbar(title=paste0("-log10(", stat_shown, ")"), barwidth = 0.75, barheight = 5)) + 
              ggplot2::theme(
                axis.ticks.y = ggplot2::element_blank(),
                axis.ticks.x = ggplot2::element_blank(),
                panel.spacing = grid::unit(0.0, "lines"),
                panel.border = ggplot2::element_blank(),
                strip.background.x = ggplot2::element_rect(fill = "#444444", colour = "white"),
                strip.text.x = ggplot2::element_text(colour = "white", angle = 0),
                axis.text.x = ggplot2::element_blank()
              ) + 
              ggplot2::xlab('Clusters') +
              ggplot2::ylab('Go terms')  +
              ggplot2::facet_grid(~cluster, scale="free_x") 
          
})

############################################################################
##    plot_markers_to_clusters
############################################################################

#' @title Map cell markers to clusters (hypergeometric/jaccard)
#'
#' @description
#' Use jaccard index and hypergeometric to map a list of cell markers to gene clusters.
#'
#' @param object A ClusterSet object.
#' @param markers A list of cell markers.
#' @param gradient_palette A vector of colors for the gradient palette.
#' @param background The background for hypergeometric test. Typically the list of genes annotated in the considered gene ontology  (e.g. BP).
#' @param cell_type_order The order of cell types. Can be "None" or "hclust". Not supported at the moment.
#' @param only_specific Logical value indicating whether to show only specific markers.
#' @param jaccard_scale The scale factor for Jaccard coefficient.
#'
#' @return A plot showing the relationship between markers and clusters.
#'
#' @export
#' @import ggplot2
#' @importFrom reshape2 melt
#' @importFrom grid unit
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 scale_color_gradientn
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 element_line
#' @importFrom ggplot2 element_rect
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 facet_grid
#' @importFrom ggplot2 guide_colourbar
#' @importFrom ggplot2 geom_point
#'
#' @examples
#' load_example_dataset("8028126/files/pbmc3k_medium_clusters_enr")
#' library(clustermole)
#' m <- clustermole::clustermole_markers(species = "hs")
#' markers <- split(m$gene[m$organ == "Blood" & m$species == "Human"], 
#'                  m$celltype[m$organ == "Blood" & m$species == "Human"])
#' plot_markers_to_clusters(pbmc3k_medium_clusters_enr, markers=markers)
#' @keywords internal
setGeneric("plot_markers_to_clusters",
           function(object, 
                    markers=NULL,
                    gradient_palette=colors_for_gradient("Magma"),
                    background=NULL,
                    cell_type_order=c("None", "hclust"),
                    only_specific=FALSE,
                    jaccard_scale=10) {
             standardGeneric("plot_markers_to_clusters")
           })

#' @title Map cell markers to clusters (hypergeometric/jaccard)
#'
#' @description
#' Use jaccard index and hypergeometric to map a list of cell markers to gene clusters.
#'
#' @param object A ClusterSet object.
#' @param markers A list of cell markers.
#' @param gradient_palette A vector of colors for the gradient palette.
#' @param background The background for hypergeometric test. Typically the list of genes annotated in the considered gene ontology  (e.g. BP).
#' @param cell_type_order The order of cell types. Can be "None" or "hclust". Not supported at the moment.
#' @param only_specific Logical value indicating whether to show only specific markers.
#' @param jaccard_scale The scale factor for Jaccard coefficient.
#'
#' @return A plot showing the relationship between markers and clusters.
#'
#' @export
#' @import ggplot2
#' @importFrom reshape2 melt
#' @importFrom grid unit
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 scale_color_gradientn
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 element_line
#' @importFrom ggplot2 element_rect
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 facet_grid
#' @importFrom ggplot2 guide_colourbar
#' @importFrom ggplot2 geom_point
#' @importFrom stats as.formula
#'
#' @examples
#' load_example_dataset("8028126/files/pbmc3k_medium_clusters_enr")
#' library(clustermole)
#' m <- clustermole_markers(species = "hs")
#' markers <- split(m$gene[m$organ == "Blood" & m$species == "Human"], 
#'                  m$celltype[m$organ == "Blood" & m$species == "Human"])
#' plot_markers_to_clusters(pbmc3k_medium_clusters_enr, markers=markers)
setMethod("plot_markers_to_clusters",
          signature(object = "ClusterSet"),
          function(object, 
                   markers=NULL,
                   gradient_palette=c("white", colors_for_gradient("De2")),
                   background=NULL,
                   cell_type_order=c("None", "hclust"),
                   only_specific=FALSE,
                   jaccard_scale=10) {

            cell_type_order <- match.arg(cell_type_order)
            
            if(is.null(background)){
              print_msg("WARNING: Running without a background but you should probably provide one.", msg_type = "INFO")
            }
            
            print_msg("Checking object format", msg_type = "DEBUG")
            check_format_cluster_set(object)
            
            if(length(object@gene_clusters) == 0){
              print_msg("No cluster found in the object.", msg_type = "STOP")
            }
            
            if(!is.list(markers))
              print_msg("The markers argument should be a list.", msg_type = "STOP")
            
            if(is.null(markers) | length(markers)==0){
              print_msg("No marker found.", msg_type = "STOP")
            }
            
            if(only_specific){
              for(i in 1:length(markers)){
                cell_marker_curr <- markers[[i]] 
                all_other_marker <- unique(unlist(cell_marker_curr[(1:length(markers))[-i]]))
                cell_marker_curr <- cell_marker_curr[!cell_marker_curr %in% all_other_marker]
                markers[[i]] <- cell_marker_curr
              }
              
              nb_cell_markers <- lapply(markers, length)
              markers <- markers[nb_cell_markers > 0]
            }
            
            comparison_jaccard <- compare_genesets(object@gene_clusters, 
                                                markers, 
                                                background, 
                                                stat = "jaccard")
            comparison_jaccard_m <- reshape2::melt(comparison_jaccard)
            
            comparison_hyper <- compare_genesets(object@gene_clusters, 
                                                   markers, 
                                                   background, 
                                                   stat = "hypergeom")
            
            comparison_hyper_m <- reshape2::melt(comparison_hyper)
            
            comparisons_m <- comparison_jaccard_m
            colnames(comparisons_m) <- c("cluster", "marker", "jaccard")

            comparisons_m$pvalue <- comparison_hyper_m$value
            comparisons_m$pvalue[comparisons_m$pvalue < 1e-320] <-  1e-320

            comparisons_m$pvalue <- -log10(comparisons_m$pvalue)
            
            comparisons_m$cluster <- factor(comparisons_m$cluster)
            comparisons_m$marker <- factor(comparisons_m$marker)
            
            m_gg <- comparisons_m[order(comparisons_m[,"cluster"], 
                                        comparisons_m[,"marker"], 
                                        comparisons_m[,"pvalue"]), ]
            
            
            cluster <- marker <- pvalue <- jaccard <- NULL
            print_msg("Running ggplot...", msg_type = "DEBUG")
            
            ggplot(data=m_gg, mapping=aes(x=cluster, y = marker, color=pvalue, size=jaccard*jaccard_scale)) +
              ggplot2::theme_bw() + 
              geom_point(shape=15) +
              ggplot2::scale_color_gradientn(colours = gradient_palette,
                                             guide = guide_colourbar(title=paste0("-log10(pvalue)"), barwidth = 0.75, barheight = 5)) + 
              ggplot2::theme(
                axis.ticks.x = ggplot2::element_blank(),
                axis.ticks.y = ggplot2::element_line(color="gray"),
                panel.spacing = grid::unit(0.0, "lines"),
                panel.border = element_rect(size=0.1, color="gray"), 
                panel.grid = element_blank(),
                strip.background.x = ggplot2::element_rect(fill = "#444444", colour = "white"),
                strip.text.x = ggplot2::element_text(colour = "white", angle = 0),
                axis.text.x = ggplot2::element_blank()
              ) + 
              ggplot2::xlab('Clusters') +
              ggplot2::ylab('Markers')  +
              ggplot2::facet_grid(~cluster, scale="free_x") 
})
