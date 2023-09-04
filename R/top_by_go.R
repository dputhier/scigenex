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
#'\dontrun{
#' # Set verbosity to 1 to display info messages only.
#' set_verbosity(1)
#' 
#' # Load a dataset
#' load_example_dataset('7871581/files/pbmc3k_medium_clusters')
#' 
#' # DNA Binding: "GO:0003677"
#' pbmc3k_medium_clusters <- top_by_go(pbmc3k_medium_clusters, go_id = "GO:0003677")
#' pbmc3k_medium_clusters@top_genes
#' }
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
#' @importFrom biomaRt listDatasets useMart getBM
setMethod("top_by_go", 
          signature("ClusterSet"), 
          function(object,
                    go_id = "GO:0003677",
                    species = "hsapiens",
                    gene_id = "hgnc_symbol",
                    host="https://www.ensembl.org") {
            
            ensembl <- biomaRt::useEnsembl(biomart = "ensembl")
            
            if(!species  %in% gsub("_gene_ensembl", "", biomaRt::listDatasets(mart = ensembl)$dataset))
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
