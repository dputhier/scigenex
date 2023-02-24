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
  if (!is.null(verbosity_value) &
      verbosity_value >= 0 & is.numeric(verbosity_value)) {
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
  if (is.null(unlist(options()["scigenex_verbosity"]))) {
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
#' set_verbosity(1)
#' print_msg("Hello world!", "INFO")
#' set_verbosity(3)
#' print_msg("Debugging message", "DEBUG")
#' set_verbosity(0)
#' print_msg("Hello world!", "INFO")
#' print_msg("Debugging message", "DEBUG")
#' options(warn=0)
#' print_msg("Warning message", "WARNING")
#' options(warn=-1)
#' print_msg("A warning message not displayed", "WARNING")
#' options(warn=opt_warn)
print_msg <-
  function(msg,
           msg_type = c("INFO", "DEBUG", "WARNING", "STOP")) {
    if (is.null(unlist(options()["scigenex_verbosity"]))) {
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
    if (msg_type == "STOP")
      stop(paste0("|-- STOP : ", msg), call. = FALSE)
  }

#################################################################
##    print_stat
#################################################################
#' Mostly a debugging function that will print some summary
#' statistics about a numeric vector, matrix or dataframe.
#'
#' @param msg The message to users.
#' @param data  The vector (numeric) for which the stats are to be produced
#' (a vector, )
#' @param msg_type The type of message, one of "INFO", "DEBUG", or "WARNING"
#' @param round_val Round the values in its first argument to the specified number
#'  of decimal. Set argument to -1 for no rounding
#' @return None
#'
#' @export
#' @examples
#' print_stat("My data", 1:10, msg_type="INFO")
#' set_verbosity(3)
#' print_stat("My data", matrix(rnorm(10), nc=2), msg_type="DEBUG")
#' set_verbosity(0)
#' print_stat("My data", matrix(rnorm(10), nc=2), msg_type="DEBUG")
#' (opt_warn <- options()$warn)
#' print_stat("My data", matrix(rnorm(10), nc=2), msg_type="WARNING")
#' options(warn=-1)
#' print_stat("My data", matrix(rnorm(10), nc=2), msg_type="WARNING")
#' options(warn=opt_warn)
print_stat <-
  function(msg,
           data,
           round_val = 2,
           msg_type = c("INFO", "DEBUG", "WARNING")) {
    
    msg_type <- match.arg(arg = msg_type, c("DEBUG", "WARNING", "INFO"))
    
    if (inherits(data, "data.frame")) {
      data <- as.matrix(data)
    }
    
    data <- as.vector(data)
    
    if (!is.numeric(data)) {
      
      print_msg("Can't print stats from non numeric object", msg_type = "WARNING")
      stats="No Statistics"
    }else{
      stats <- summary(data)
      names(stats) <- c("Min", "Q1", "Med", "Mean", "Q3", "Max")
      
      if (round_val > 0 & is.numeric(round_val)) {
        stats <- round(stats, round_val)
      }
      
      stats <- paste(names(stats), stats, sep = ":", collapse = " ")
    }
    
    print_msg(paste0(msg, ": ", stats), msg_type = msg_type)
    
  }


#################################################################
##    A simple function to create a random string
#################################################################
#' Generate a random string of letters and numbers
#'
#' This function generates a random string of 10 characters, consisting of 3 uppercase letters,
#' 4 digits, and 3 lowercase letters.
#'
#' @return A character string of length 10, consisting of random letters and numbers
#' 
#' @export
#' 
#' @examples
#' create_rand_str()
create_rand_str <- function() {
  v <- c(
    sample(LETTERS, 3, replace = TRUE),
    sample(0:9, 4, replace = TRUE),
    sample(letters, 3, replace = TRUE)
  )
  return(paste0(sample(v), collapse = ""))
}

#################################################################
##    A simple function that return a dataset (1)
#################################################################
#' Generate an example dataset with three clusters of profiles 
#'
#' Generate an artificial matrix with with random noise and 3 
#' clusters of 'expression' profiles (as row). 
#'
#' @return a matrix.
#' 
#' @export
#' 
#' @examples
#' m <- create_3_rnd_clust()
#' 
create_3_rnd_clust <- function(){
  set.seed(123)
  m <- matrix(rnorm(80000), nc=20)
  m[1:100, 1:10] <- m[1:100, 1:10] + 4
  m[101:200, 11:20] <- m[101:200, 11:20] + 3
  m[201:300, 5:15] <- m[201:300, 5:15] -2
  return(m)
}

#################################################################
##    A simple function that return a dataset (2)
#################################################################
#' Generate an example dataset with four clusters of profiles 
#'
#' Generate an artificial matrix with with random noise and 4 
#' clusters of 'expression' profiles (as row). 
#'
#' @return a matrix.
#' 
#' @export
#' 
#' @examples
#' m <- create_4_rnd_clust()
#' 
#' 
create_4_rnd_clust <- function(){
  set.seed(123)
  m <- matrix(rnorm(80000), nc=20)
  m[1:100, 1:10] <- m[1:100, 1:10] + 2
  m[101:200, 11:20] <- m[101:200, 11:20] + 4
  m[301:400, 5:15] <- m[201:300, 5:15] - 2
  m[201:300, c(1,5,10,15,20)] <- m[201:300, c(1,5,10,15,20)] - 3
  return(m)
}




