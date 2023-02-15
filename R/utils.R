#################################################################
##    set_verbosity
#################################################################
#' Set the verbosity level for the SciGeneX package
#'
#' This function sets the verbosity level for the SciGeneX package, 
#' which controls the amount of information that is printed to the console by 
#' the \code{\link{print_msg}} function. The verbosity level can be set to
#'  any non-negative integer, with higher values indicating more detailed output. 
#'  By default, the verbosity level is set to 1.
#'
#' @param verbosity_value A non-negative integer indicating the verbosity level to be set.
#'
#' @return NULL
#'
#' @export
#'
#' @examples
#' # Set verbosity level to 2
#' set_verbosity(2)
#'
#' # Set verbosity level to 0
#' set_verbosity(0)

# 0 : No message
# 1 : Display only INFO type message
# 2 : Display both INFO and DEBUG type message

set_verbosity <- function(verbosity_value) {
  if(!is.null(verbosity_value) & verbosity_value >=0 & is.numeric(verbosity_value)){
    options(scigenex_verbosity = verbosity_value)
  }
}

#################################################################
##    get_verbosity()
#################################################################
#' Get the current verbosity level.
#'
#' This function get the verbosity level of the SciGeneX package which 
#' controls the amount of information that is printed to the console by 
#' the \code{\link{print_msg}} function. 
#'
#'
#' @return A vector
#'
#' @export
#'
#' @examples
#' get_verbosity()
#'

get_verbosity <- function() {
  if(is.null(unlist(options()["scigenex_verbosity"]))) {
    options(scigenex_verbosity = 1)
  }
  return(options()$scigenex_verbosity)
}

#################################################################
##    print_msg
#################################################################
#' Print a message based on the level of verbosity
#'
#' @param msg The message to be printed
#' @param msg_type The type of message, one of "INFO", "DEBUG", or "WARNING"
#'
#' @return None
#'
#' @export
#' @examples
#' print_msg("Hello world!", "INFO")
#' print_msg("Debugging message", "DEBUG")
#' print_msg("Warning message", "WARNING")
print_msg <- function(msg, msg_type = "INFO") {
  if(is.null(unlist(options()["scigenex_verbosity"]))) {
    options(scigenex_verbosity = 1)
  }
  if (msg_type == "INFO")
    if (unlist(options()["scigenex_verbosity"]) > 0)
      cat(paste("|-- INFO : ", msg, "\n"))
  if (msg_type == "DEBUG")
    if (unlist(options()["scigenex_verbosity"]) > 1)
      cat(paste("|-- DEBUG : ", msg, "\n"))
  if (msg_type == "WARNING")
    warning("|-- WARNING : ", msg, call. = FALSE)
}