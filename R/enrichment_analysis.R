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
#' # Load a Seurat object
#' data(pbmc_small, package = "SeuratObject")
#' 
#' # Compute the signatures using find_gene_clusters()
#' clust_set <- select_genes(pbmc_small,
#'                           distance = "kendall",
#'                           k = 75,
#'                           highest = 0.3,
#'                           fdr = 1e-8,
#'                           row_sum = -Inf,
#'                           no_dknn_filter = TRUE)
#'                           
#' # Cluster informative features
#' clust_set <- gene_clustering(clust_set, 
#'                             inflation = 1.6,
#'                             keep_nn = FALSE,
#'                             k = 5)
#' 
#' # Do enrichment analysis using GO ontology
#' clust_set <- enrich_go(clust_set[1:2,])
#' 
#' # Results of enrichment analysis are 'gene_cluster_annotations' slot
#' print(clust_set@gene_cluster_annotations)
#' 
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
              org_db <- org.Hs.eg.db
              go_name <- "org.Hs.eg.db"
              print_msg(msg_type = "INFO",
                        msg = "Species used : Homo sapiens")
            }
            
            if (species == "Mmusculus") {
              org_db <- org.Mm.eg.db
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
#' # Load a Seurat object
#' data(pbmc_small, package = "SeuratObject")
#' 
#' # Compute the signatures using find_gene_clusters()
#' clust_set <- select_genes(pbmc_small,
#'                           distance = "kendall",
#'                           k = 75,
#'                           highest = 0.3,
#'                           fdr = 1e-8,
#'                           row_sum = -Inf,
#'                           no_dknn_filter = TRUE)
#'                           
#' # Cluster informative features
#' clust_set <- gene_clustering(clust_set, 
#'                             inflation = 1.6,
#'                             keep_nn = FALSE,
#'                             k = 5)
#' 
#' # Do enrichment analysis using GO ontology
#' # (here also subset the two first clusters)
#' clust_set <- enrich_go(clust_set[1:2,])
#' viz_enrich(clust_set, cluster = 1, nb_terms = 5)
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
