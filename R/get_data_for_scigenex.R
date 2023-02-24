#################################################################
##    Getting data
#################################################################

#' @title
#' Fetch an expression matrix from a file, dataframe or Seurat object.
##' @description
#' This function fetchs an expression matrix from a file, dataframe or Seurat object.
#' @param data A \code{matrix}, \code{data.frame} or \code{Seurat} object.
#'
#' @return An expression matrix.
#' @export get_data_for_scigenex
#'
#' @examples
#' 
#' \dontrun{
#' ## with an artificial dataset
#'
#' m <- create_3_rnd_clust()
#' 
#' res <- get_data_for_scigenex(data=m)
#' }
#' 
get_data_for_scigenex <- function(data = NULL) {
  
  ## getting matrix (probesID vs SamplesID)
  
  # Stop the function if data not provided
  if (is.null(data)) {
    print_msg(paste0("Please provide a Seurat Object, a data.frame",
                      ", a matrix or dgCMatrix.\n"), 
              msg_type = "STOP")
  }
  
  if (inherits(data, "Seurat")) { 
    # Extract normalized count matrix in Seurat object
    data <- SeuratObject::GetAssayData(data, slot = 'data')
  } else if (is.data.frame(data)) {
    # Convert dataframe to a matrix
    data <- as.matrix(data)
  } else if (is.matrix(data)) {
    data <- data
  } else if (inherits(data, "dgCMatrix")) {
    data <- data
  } else {
        print_msg(paste0("\t--> Please provide a Seurat Object, a data.frame",
                         " or a matrix."), 
                  msg_type = "STOP")
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
