#################################################################
##    Getting data
#################################################################

#' @title
#' Fetch an expression matrix from a file, dataframe or Seurat object.
##' @description
#' This function fetchs an expression matrix from a file, dataframe or Seurat object.
#' @param data A \code{matrix}, \code{data.frame} or \code{Seurat} object.
#' @param filename A character string representing the file name.
#' @param path A character string representing the data directory where
#' intermediary files are to be stored. Default to current working directory.
#'
#' @return A list containing a matrix and the filename (if filename argument is used).
#' @export get_data_4_DBFMCL
#'
#' @examples
#' 
#' \dontrun{
#' ## with an artificial dataset
#'
#' m <- matrix(rnorm(80000), nc = 20)
#' res <- get_data_4_DBFMCL(data=m)
#' }
#' 
get_data_4_DBFMCL <- function(data = NULL, filename = NULL, path = ".") {
  
  ## getting matrix (probesID vs SamplesID)
  if (!is.null(data)) {
    if (inherits(data, "Seurat")) {
      data <- as.matrix(GetAssayData(data, slot = 'data'))
    }
    else if (is.data.frame(data)) {
      data <- as.matrix(data)
    }
    if (!is.matrix(data)) {
      stop(
        "\t--> Please provide a Seurat Object, a data.frame",
        " or a matrix.\n"
      )
    }
    name <- NULL
  }
  else {
    if (!is.null(filename)) {
      data <- as.matrix(read.table(file.path(path, filename),
                                   sep = "\t", header = TRUE, row.names = 1, quote = ""
      ))
      name <- unlist(strsplit(filename, "\\."))[1]
    }
    else {
      stop(
        "\t--> Please provide an ExpressionSet, a data.frame, ",
        "a matrix or a tabular file\n"
      )
    }
  }
  ## adding dimnames if not provided
  if (is.null(rownames(data))) {
    print_msg("Row names not provided. Adding.", msg_type = "DEBUG")
    rownames(data) <- paste("gene", 1:nrow(data), sep = "")
  }
  if (is.null(colnames(data))) {
    print_msg("Colum names not provided. Adding.", msg_type = "DEBUG")
    colnames(data) <- paste("sample", 1:ncol(data), sep = "")
  }
  
  return(list(data = data, name = name))
}
