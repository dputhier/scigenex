#################################################################
##    Define top_by_go function for a ClusterSet object
#################################################################
#' @title Select top_genes based on go terms
#' @description
#' The clusterSet object contains a top_genes slot that can be used to display 
#' genes in heatmaps (see \code{plot_heatmap} function). Here the function select 
#' top_genes based on GO terms.
#' @param object A \code{ClusterSet} object.
#' @param go_id A GO ID.
#' @param species The target species according to biomaRt::listDatasets(mart = ensembl)$dataset 
#' (e.g hsapiens, mmusculus, rnorvegicus...).
#' @param gene_id The type of gene identifier in the clusterSet object (default hgnc_symbol). 
#' use biomaRt::listFilters() to get the vailable list.
#' @param host Host to connect to. Defaults to www.ensembl.org.
#' @return A \code{ClusterSet} object.
#' @export top_by_go
#' @keywords internal
#' @examples
#' # Set verbosity to 1 to display info messages only.
#' set_verbosity(1)
#' 
#' # Load a dataset
#' load_example_dataset('7871581/files/pbmc3k_medium_clusters')
#' 
#' # DNA Binding: "GO:0003677"
#' pbmc3k_medium_clusters <- top_by_go(pbmc3k_medium_clusters, go_id = "GO:0003677")
#' pbmc3k_medium_clusters@top_genes
#' 
#' # Cell surface receptor signaling pathway: "GO:0007166"
#' pbmc3k_medium_clusters <- top_by_go(pbmc3k_medium_clusters, go_id =c("GO:0007166"))
#'                                                                     
#' # Signaling receptor binding: GO:0005102 
#' pbmc3k_medium_clusters <- top_by_go(pbmc3k_medium_clusters, go_id="GO:0005102")
#' 
#' # ECM genes
#' # GO:0031012 - extracellular matrix; 
#' # GO:0005578 - proteinacious extracellular matrix; 
#' # GO:0005201 - extracellular matrix structural constituent; 
#' # GO:1990430 - extracellular matrix protein binding; and 
#' # GO:0035426 - extracellular matrix cell signalling).
#' pbmc3k_medium_clusters <- top_by_go(pbmc3k_medium_clusters, go_id=c("GO:0031012", "GO:0005578",
#'                                                                    "GO:0005201", "GO:1990430",
#'                                                                    "GO:0035426"))
setGeneric("top_by_go", 
           function(object,
                    go_id = "GO:0003677",
                    species = "hsapiens",
                    gene_id = "hgnc_symbol",
                    host="https://www.ensembl.org")
             standardGeneric("top_by_go")
)

#' @title Select top_genes based on go terms
#' @description
#' The clusterSet object contains a top_genes slot that can be used to display 
#' genes in heatmaps (see \code{plot_heatmap} function). Here the function select 
#' top_genes based on GO terms.
#' @param object A \code{ClusterSet} object.
#' @param go_id A GO ID.
#' @param species The target species according to biomaRt::listDatasets(mart = ensembl)$dataset 
#' (e.g hsapiens, mmusculus, rnorvegicus...).
#' @param gene_id The type of gene identifier in the clusterSet object (default hgnc_symbol). 
#' use biomaRt::listFilters() to get the vailable list.
#' @param host Host to connect to. Defaults to www.ensembl.org.
#' @return A \code{ClusterSet} object.
#' @export top_by_go
#'
#' @examples
#' # Set verbosity to 1 to display info messages only.
#' set_verbosity(1)
#' 
#' # Load a dataset
#' load_example_dataset('7871581/files/pbmc3k_medium_clusters')
#' 
#' # DNA Binding: "GO:0003677"
#' pbmc3k_medium_clusters <- top_by_go(pbmc3k_medium_clusters, go_id = "GO:0003677")
#' pbmc3k_medium_clusters@top_genes
#' 
#' # Cell surface receptor signaling pathway: "GO:0007166"
#' pbmc3k_medium_clusters <- top_by_go(pbmc3k_medium_clusters, go_id =c("GO:0007166"))
#'                                                                     
#' # Signaling receptor binding: GO:0005102 
#' pbmc3k_medium_clusters <- top_by_go(pbmc3k_medium_clusters, go_id="GO:0005102")
#' 
#' # ECM genes
#' # GO:0031012 - extracellular matrix; 
#' # GO:0005578 - proteinacious extracellular matrix; 
#' # GO:0005201 - extracellular matrix structural constituent; 
#' # GO:1990430 - extracellular matrix protein binding; and 
#' # GO:0035426 - extracellular matrix cell signalling).
#' pbmc3k_medium_clusters <- top_by_go(pbmc3k_medium_clusters, go_id=c("GO:0031012", "GO:0005578",
#'                                                                    "GO:0005201", "GO:1990430",
#'                                                                    "GO:0035426"))
setMethod("top_by_go", 
          signature("ClusterSet"), 
          function(object,
                    go_id = "GO:0003677",
                    species = "hsapiens",
                    gene_id = "hgnc_symbol",
                    host="https://www.ensembl.org") {
            
            ensembl <- biomaRt::useEnsembl(biomart = "ensembl")
            
            if(!species  %in% gsub("_gene_ensembl", "", listDatasets(mart = ensembl)$dataset))
              print_msg("Unknow species in biomaRt", msg_type = "STOP")
            
            print_msg("Connecting to Biomart...", msg_type = "INFO")
            
            mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", 
                            dataset = paste0(species, "_gene_ensembl"),
                            host=host) 
            
            top_gn <- list()
            
            for(i in 1:length(object@gene_clusters)){
              
              print_msg(paste0("Processing cluster ", i), msg_type = "INFO")
              
              go_list <- biomaRt::getBM(attributes=c("go_id",
                                            "hgnc_symbol"),
                               filters = "hgnc_symbol",
                               values = object@gene_clusters[[i]],
                               mart=mart)
              
              to_keep <- go_list$hgnc_symbol[go_list$go_id %in% go_id]
              top_gn[[i]] <- sort(to_keep)
            }
            
            top_gn <- lapply(top_gn, unique)
            top_gn <- lapply(top_gn, sort)
            
            object@top_genes <- top_gn
            
            return(object)
})
