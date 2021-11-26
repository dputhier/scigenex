#################################################################
##    Running gprofiler2 to perform enrichment analysis
#################################################################

#' @title
#' Perform gene enrichment analysis.
##' @description
#' Perform enrichment analysis on all MCL clusters indepentently and store the results in the cluster_annotations slot of the ClusterSet object.
#' @param object A \code{ClusterSet} object.
#' @param specie Specie name, as a concatenation of the first letter of the name (uppercas) and the family name, e.g human - Hsapien
#' @param ontology One of "BP", "MF", and "CC" subontologies, or "ALL" for all three.
#' @param verbose Whether or not to print progression in the console.
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
                    specie="Hsapiens",
                    ontology="ALL",
                    verbose = TRUE) {
             standardGeneric("enrich_go")
           })

#' @rdname enrich_go
setMethod("enrich_go",
          signature(object = "ClusterSet"),
          function(object,
                   specie="Hsapiens",
                   ontology="ALL",
                   verbose = TRUE) {
            
            
            if (!specie %in% c("Hsapiens", "Mmusculus")){
              stop("Specie name provided doesn't exists.")}
            
            if (specie == "Hsapiens") {
              hs <- org.Hs.eg.db
              go_name <- "org.Hs.eg.db"
              if(verbose) {print_msg(msg_type = "INFO",
                                     msg = "Specie used : Homo sapiens")}
            }
            
            if (specie == "Mmusculus") {
              mm <- org.Mm.eg.db
              go_name <- "org.Mm.eg.db"
              if(verbose) {print_msg(msg_type = "INFO",
                                     msg = "Specie used : Mus musculus")}
            }
            
            
            for(cluster in unique(object@cluster)){
              if (verbose) {print(paste0("Enrichment analysis for cluster ", cluster))}
              cluster_name = paste0("Cluster_", cluster)
              query = rownames(object@data[object@cluster == cluster,])
              
              # Convert gene id to EntrezId format
              suppressMessages(query_entrezid <- select(hs, 
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
                enrich_go_res <- enrichGO(de_entrezid[,"ENTREZID"],
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
                    type = "dotplot",
                    verbose = TRUE) {
             standardGeneric("enrich_viz")
           })

#' @rdname enrich_viz
setMethod("enrich_viz",
          signature(object = "ClusterSet"),
          function(object,
                   clusters = "all",
                   type = "dotplot",
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
              
              # Check if there is a result provided by enrich_go function for the current cluster
              if(is.null(object@cluster_annotations[[cur_cluster]]@result)){
                print_msg(msg_type = "WARNING",
                          msg = paste0("No functional enrichment analysis results for cluster ", cur_cluster, ".")) #Continue through the next cluster without plotting
              } else {
                
                # Print an informative message to announce plot of the results of the current cluster
                if (verbose) {
                  print_msg(msg_type = "INFO",
                            msg = paste0("Plot enrichment analysis results for cluster ", cur_cluster))
                }
                
                
                # Create a ggplot - dotplot
                if (type = "dotplot"){
                  dotplot(object@cluster_annotations[[cur_cluster]], showCategory=30)
                }
                
                # Create a ggplot - barplot
                if (type = "dotplot"){
                  barplot(object@cluster_annotations[[cur_cluster]], showCategory=20)
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