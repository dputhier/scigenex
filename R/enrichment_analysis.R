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
#' go_res <- enrich_analysis(myobject)
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