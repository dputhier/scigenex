#################################################################
##    Getting data
#################################################################

#' @title
#' Fetch an expression matrix from a file, dataframe or Seurat object.
##' @description
#' This function fetchs an expression matrix from a file, dataframe or Seurat object.
#' @param data A \code{matrix}, \code{data.frame} or \code{Seurat} object.
#' @param layer One of "SCT", "data" or count. The slot to extract from the seurat object to perform clustering analysis.
#' @param assay The assay to use in the Seurat object. If NULL, the function will try to guess.
#' SCT is the recommended method from Seurat package when working with spatial transcriptomics data.
#'
#' @return An expression matrix.
#'
#' @examples
#' ## with an artificial dataset
#' m <- create_3_rnd_clust()
#'
#' res <- get_data_for_scigenex(data=m)
#' 
#' @export
get_data_for_scigenex <- function(data = NULL,
                                  layer = c("data", "sct", "counts"),
                                  assay=NULL) {
  layer <- match.arg(layer)
  
  # Stop the function if no data were provided
  if (is.null(data)) {
    print_msg(paste0(
                    "Please provide a Seurat Object, a data.frame",
                    ", a matrix or dgCMatrix.\n",
                  msg_type = "STOP"
    ))
  }
  
  if (inherits(data, "Seurat")) {
    
    all_assays <- names(data@assays)
    
    if(is.null(assay)){
      assay <- all_assays[1]
    }else{
      if(!assay %in% all_assays)
        print_msg("Assay was not found in Seurat object.", msg_type="STOP")
    }
    
    print_msg(paste0("Extracting '", layer, "' layer from Seurat object"),
              msg_type = "DEBUG")
    
    if(!layer %in% SeuratObject::Layers(data[[assay]]))
      print_msg("Layer not found in Seurat object assay.", msg_type = "STOP")
    
    rn  <- rownames(data[[assay]])
    cn <- colnames(data[[assay]])
    data <- SeuratObject::LayerData(data, assay=assay, layer=layer)
    colnames(data) <- cn
    rownames(data) <- rn

  } else if (is.data.frame(data)) {
    print_msg("Converting dataframe to matrix", msg_type = "DEBUG")
    data <- as.matrix(data)
  } else if (is.matrix(data)) {
    data <- data
  } else if (inherits(data, "dgCMatrix")) {
    data <- data
  } else {
    print_msg(
      paste0(
        "\t--> Please provide a Seurat Object, a data.frame",
        " or a matrix."
      ),
      msg_type = "STOP"
    )
  }
  
  ## adding dimnames if not provided
  if (is.null(rownames(data))) {
    print_msg("Row names not provided. Adding.", msg_type = "DEBUG")
    rownames(data) <- paste("gene", 1:nrow(data), sep = "")
  }
  
  if (is.null(colnames(data))) {
    print_msg("Column names not provided. Adding.", msg_type = "DEBUG")
    colnames(data) <- paste("sample", 1:ncol(data), sep = "")
  }
  
  return(data)
}
