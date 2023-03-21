#################################################################
##    Getting data
#################################################################

#' @title
#' Fetch an expression matrix from a file, dataframe or Seurat object.
##' @description
#' This function fetchs an expression matrix from a file, dataframe or Seurat object.
#' @param data A \code{matrix}, \code{data.frame} or \code{Seurat} object.
#' @param which_slot One of "SCT", "data" or count. The slot to extract from the seurat object to perform clustering analysis.
#' SCT is the recommended method from Seurat package when working with spatial transcriptomics data.
#'
#' @return An expression matrix.
#' @export get_data_for_scigenex
#'
#' @examples
#' ## with an artificial dataset
#' m <- create_3_rnd_clust()
#'
#' res <- get_data_for_scigenex(data=m)
#' 
#'
get_data_for_scigenex <- function(data = NULL,
                                  which_slot = c("data", "sct", "counts")) {
  which_slot <- match.arg(which_slot)
  
  # Stop the function if data not provided
  if (is.null(data)) {
    print_msg(
      paste0(
        "Please provide a Seurat Object, a data.frame",
        ", a matrix or dgCMatrix.\n"
      ),
      msg_type = "STOP"
    )
  }
  
  if (inherits(data, "Seurat")) {
    print_msg(paste0("Extracting '", which_slot, "' slot from Seurat object"),
              msg_type = "DEBUG")
    
    if (which_slot %in% c("data", "counts")) {
      data <- SeuratObject::GetAssayData(data, slot = which_slot)
    } else if (which_slot == "sct") {
      if ("SCT" %in% names(data@assays)) {
        data <- data@assays$SCT@data
      } else{
        print_msg("This object has no 'SCT' slot. Use SCTransform() before.",
                  msg_type = 'STOP')
      }
    }
    
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
