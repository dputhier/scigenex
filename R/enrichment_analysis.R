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
#' @export enrich_go
#'
#' @examples
#' 
#' \dontrun{
#' ## Assuming myobject is a ClusterSet object with at least 1 cluster.
#'
#' go_res <- enrich_go(myobject)
#' }

setGeneric("enrich_go",
           function(object,
                    species="Hsapiens",
                    ontology="all") {
             standardGeneric("enrich_go")
           })

#' @rdname enrich_go
setMethod("enrich_go",
          signature(object = "ClusterSet"),
          function(object,
                   species="Hsapiens",
                   ontology="all") {
            
            
            if (!species %in% c("Hsapiens", "Mmusculus")){
              stop("Species name provided doesn't exists.")}
            
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
            
            
            for(cluster in unique(object@gene_clusters)){
              print(paste0("Enrichment analysis for cluster ", cluster))
              cluster_name = paste0("Cluster_", cluster)
              query = rownames(object@data[object@gene_clusters == cluster,])
              
              # Convert gene id to EntrezId format
              suppressMessages(query_entrezid <- AnnotationDbi::select(org_db, 
                                                                       keys = query,
                                                                       columns = c("ENTREZID", "SYMBOL"),
                                                                       keytype = "SYMBOL"))
              
              print_msg(msg_type = "DEBUG",
                        msg = paste0("Cluster ",
                                     cluster,
                                     "\n",
                                     "Number of gene names converted to EntrezId : ", 
                                     length(which(!is.na(query_entrezid[, "ENTREZID"]))),
                                     "\n", 
                                     "Number of gene names not converted to EntrezId : ", 
                                     length(which(is.na(query_entrezid[, "ENTREZID"])))))
              
              # Enrichment analysis using GO database
              enrich_go_res <- clusterProfiler::enrichGO(query_entrezid[,"ENTREZID"],
                                        OrgDb = go_name,
                                        ont = ontology,
                                        readable = TRUE)
              # Store results in the ClusterSet object
              object@cluster_annotations[[cluster]] <- enrich_go_res
            }
            return(object)
          }
)

#################################################################
##    Visualization of enrichment analyses results from a ClusterSet object
#################################################################

#' @title
#' Visualization of enrichment analyses results from a ClusterSet object
##' @description
#' Display enrichment analyses results for a ClusterSet object.
#' @param object A \code{ClusterSet} object.
#' @param clusters  A vector of cluster id to plot.
#' @param type A vector providing the type of plot.
#' @param nb_terms An integer indicating the number of terms in the plot.
#' @param terms_size An interger indicating the wrap length of terms
#'
#' @return A \code{ClusterSet} object
#' @export viz_enrich
#'
#' @examples
#' 
#' \dontrun{
#' ## Assuming myobject is a ClusterSet object with at least 1 cluster.
#'
#' viz_enrich(myobject)
#' }

setGeneric("viz_enrich",
           function(object,
                    clusters = "all",
                    type = "dotplot",
                    nb_terms = 20,
                    terms_size = 50) {
             standardGeneric("viz_enrich")
           })

#' @rdname viz_enrich
setMethod("viz_enrich",
          signature(object = "ClusterSet"),
          function(object,
                   clusters = "all",
                   type = c("dotplot", "barplot"),
                   nb_terms = 20,
                   terms_size = 50) {
            
            if (length(clusters) == 1){
              if (clusters == "all"){
                clusters <- unique(object@gene_clusters)
              }
            }
            
            
            for (cur_cluster in clusters) {
              # Check if the current cluster id provided exists
              if (!(cur_cluster %in% unique(object@gene_clusters))) {
                stop(paste0("Cluster ", cur_cluster, " doesn't exist."))
              }
              
              # Check if there is a result provided by enrich_go function for the current cluster
              if(nrow(object@cluster_annotations[[as.numeric(cur_cluster)]]@result) == 0){
                print_msg(msg_type = "WARNING",
                          msg = paste0("No functional enrichment analysis results for cluster ", cur_cluster, ".")) #Continue through the next cluster without plotting
              } else {
                
                # Print an informative message to announce plot of the results of the current cluster
                print_msg(msg_type = "INFO",
                          msg = paste0("Plot enrichment analysis results for cluster ", cur_cluster))
                
                
                # Create a ggplot - dotplot
                if ("dotplot" %in% type){
                  if(object@cluster_annotations[[as.numeric(cur_cluster)]]@ontology == "GOALL"){
                    dot_plot <- enrichplot::dotplot(object@cluster_annotations[[as.numeric(cur_cluster)]],
                                                    split="ONTOLOGY",
                                                    showCategory=nb_terms,
                                                    label_format = terms_size)
                    dot_plot <- dot_plot + facet_grid(ONTOLOGY~., scales="free")
                    dot_plot <- dot_plot + ggtitle(paste0("Gene cluster ", cur_cluster))
                  } else {
                    dot_plot <- enrichplot::dotplot(object@cluster_annotations[[cur_cluster]],
                                                    showCategory=nb_terms,
                                                    label_format = terms_size)
                    dot_plot <- dot_plot + ggtitle(paste0("Gene cluster ", cur_cluster))
                  }
                  object@cluster_annotations[[paste0("plot_cl",cur_cluster)]]$dotplot <- dot_plot
                }
                
                # Create a ggplot - barplot
                if ("barplot" %in% type){
                  if(object@cluster_annotations[[as.numeric(cur_cluster)]]@ontology == "GOALL"){
                    bar_plot <- graphics::barplot(object@cluster_annotations[[as.numeric(cur_cluster)]],
                                                  split="ONTOLOGY",
                                                  showCategory=nb_terms,
                                                  label_format = terms_size)
                    bar_plot <- bar_plot + facet_grid(ONTOLOGY~., scales="free")
                    bar_plot <- bar_plot + ggtitle(paste0("Gene cluster ", cur_cluster))
                  } else {
                    bar_plot <- graphics::barplot(object@cluster_annotations[[cur_cluster]],
                                                  showCategory=nb_terms,
                                                  label_format = terms_size)
                    bar_plot <- bar_plot + ggtitle(paste0("Gene cluster ", cur_cluster))
                  }
                  object@cluster_annotations[[paste0("plot_cl",cur_cluster)]]$barplot <- bar_plot
                }
              }
            }
            
            # Print a message to inform which slot contains the results
            print_msg(msg_type = "INFO",
                      msg = paste("Plots are stored in object@cluster_annotations[[plot_cl<cluster>]]$<type of plot> \n",
                                  "For example : object@cluster_annotations[[plot_cl1]]$dotplot"))
            return(object)
          }
)