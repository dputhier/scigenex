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
#' opt_warn <- options()$warn
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
    
    msg_type <- match.arg(msg_type)
    
    if (is.null(unlist(options()["scigenex_verbosity"]))) {
      options(scigenex_verbosity = 1)
    }
    
    if (msg_type == "DEBUG"){
      if (unname(unlist(options()["scigenex_verbosity"]) > 1))
        cat(paste("|-- DEBUG : ", msg, "\n"))
      
    }else if (msg_type == "WARNING"){
      warning("|-- WARNING : ", msg, call. = FALSE)
      
    }else if (msg_type == "STOP"){
      stop(paste0("|-- STOP : ", msg), call. = FALSE)
      
    }else{
        if (unname(unlist(options()["scigenex_verbosity"]) > 0))
          cat(paste("|-- INFO : ", msg, "\n"))
    }
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
#' opt_warn <- options()$warn
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
  n <- 80000
  m <- 20
  m <- matrix(rnorm(n), ncol=m)
  print_msg(paste0("Creating a matrix of size: ", 
                   paste0(dim(m), collapse = " x "), 
                   collapse = ""))
  
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
  m <- matrix(rnorm(80000), ncol=20)
  m[1:100, 1:10] <- m[1:100, 1:10] + 2
  m[101:200, 11:20] <- m[101:200, 11:20] + 4
  m[301:400, 5:15] <- m[201:300, 5:15] - 2
  m[201:300, c(1,5,10,15,20)] <- m[201:300, c(1,5,10,15,20)] - 3
  return(m)
}


#################################################################
##    Returns a palette for gradients
#################################################################
#' Generate a vector of colors for a gradient
#'
#' This function generates a vector of colors for a gradient, given 
#' a specified palette name.
#'
#' @param palette A character vector specifying the palette to use. One of: "Je1", 
#' "Seurat_like", "Ju1", "De1",  "De2", "De3", "De4", "De5", "De6", "De7", "De8", "De9.
#' @return A character vector of color codes.
#' @export colors_for_gradient
#' @examples
#' colors_for_gradient()
#' colors_for_gradient(palette = "Ju1")
#' 
colors_for_gradient <- function(palette=c("Je1", "Seurat_Like", "Ju1", "De1", 
                                          "De2", "De3", "De4", "De5",
                                          "De6", "De7", "De8", "De9")){
  palette <- match.arg(palette)
  
  if(palette == "Seurat_Like"){
    return(c("#5D50A3", "#9FD7A4", "#FBFDBA", "#FEB163", "#A80B44"))
  }else if(palette == "Ju1"){
    return(c("#A9D6E5", "#2166AC", "#000000", "#B2182B", "#FFCA3A"))
  }else if(palette == "De1"){
    return(c("#d73027","#fc8d59","#fee090","#e0f3f8","#91bfdb","#253494"))
  }else if(palette == "De2"){
    return(c("#FFF7FB","#ECE2F0","#D0D1E6","#A6BDDB","#67A9CF","#3690C0","#02818A","#016450"))
  }else if(palette == "De3"){
    return(c("#1A1835","#15464E","#2B6F39","#757B33","#C17A70","#D490C6","#C3C1F2","#CFEBEF"))
  }else if(palette == "De4"){
    return(c("#0000FF","#00FFFF","#80FF80","#FFFF00","#FF0000"))
  }else if(palette == "De5"){
    return(c("#0000AA","#0000FF","#00FFFF","#80FF80","#FFFF00","#FF0000","#AA0000"))
  }else if(palette == "De6"){
    return(c("#4575b4","#74add1","#abd9e9","#e0f3f8","#fee090","#fdae61","#f46d43","#d73027"))
  }else if(palette == "De7"){
    return(c("#67001f","#b2182b","#d6604d","#f4a582","#fddbc7","#f7f7f7","#d1e5f0","#92c5de","#4393c3","#2166ac","#053061"))
  }else if(palette == "De8"){
    return(c("#2b83ba","#abdda4","#fdae61","#d7191c"))
  }else if(palette == "De9"){
    return(c("#0000BF","#0000FF","#0080FF","#00FFFF","#40FFBF","#80FF80","#BFFF40","#FFFF00","#FF8000","#FF0000","#BF0000"))
  }else if(palette == "Je1"){
    return(c("#27408B", "#3A5FCD", "#3288BD", "#66C2A5","#ABDDA4", "#E6F598","#FEE08B", "#FDAE61","#F46D43","#D53E4F","#8B2323"))
  }
}

#################################################################
##    Returns a discrete color palette
#################################################################
#' Generate a discrete color palette
#'
#' This function generates a vector of colors for a discrete variable, 
#' given a specified palette name.
#'
#' @param n An integer specifying the number of colors to generate. 
#' @param palette A character vector specifying the palette to use. 
#' @return A character vector of color codes.
#' @export discrete_palette
#' @examples
#' discrete_palette()
#' discrete_palette(n=20, palette = "ggplot")
#' 
discrete_palette <- function(n=10, palette=c("Ju1", "ggplot")){
  
  palette <- match.arg(palette)
  
  if(palette == "Ju1"){
    palette <- colorRampPalette(c(  "#9F1717", "#AE5B11", "#C48D00", "#517416", 
                                    "#115C8A", "#584178", "#9D1C70", "#E96767", 
                                    "#EC9342", "#FFCA3A", "#8AC926", "#4DADE8", 
                                    "#9579B9", "#E25CB4", "#DB2020", "#DA7316", 
                                    "#F0AE00", "#6D9D1E", "#1882C0", "#71529A", 
                                    "#D02494", "#EF9292", "#F2B57D", "#FFDA77", 
                                    "#B6E36A", "#7BC4EE", "#AD98C9", "#EA8AC9"))(n)
  }else if(palette == "ggplot"){
    gg_color_hue <- function(n) {
      hues = seq(15, 375, length = n + 1)
      grDevices::hcl(h = hues, l = 65, c = 100)[1:n]
    }
    palette <- gg_color_hue(n)
  }
  
  names(palette) <- 1:n
  
  return(palette)
  
}


#################################################################
##    Check the format of a Clusterset object
#################################################################
#' Check the format of a Clusterset object
#'
#' Check the format of a Clusterset object (ncol, nrow, inherits...)
#'
#' @param object the clusterSet object to be tested 
#' @return None.
#' @examples 
#' \dontrun{
#'  set_verbosity(0)
#'  data(pbmc_small, package = "SeuratObject")
#'  clust_set <- select_genes(data=pbmc_small,
#'                            distance_method="pearson",
#'                            k=10,
#'                            row_sum=-Inf,
#'                            highest=0.95,
#'                            fdr = 1e-6)
#'  clust_set <- gene_clustering(object = clust_set,
#'                               inflation = 1.2,
#'                               keep_nn = FALSE,
#'                               k = 5,
#'                               threads = 1)
#'  check_format_cluster_set(clust_set)
#' }
#' @export check_format_cluster_set
check_format_cluster_set <- function(object) {

  if(!inherits(object, "ClusterSet"))
    print_msg("Please provide a ClusterSet as input.", 
              msg_type = "STOP")
  
  if(nrow(object) == 0)
    print_msg("The ClusterSet object does not contain any gene.", 
              msg_type = "STOP")
    
  if(ncol(object) == 0)
    print_msg("The ClusterSet object does not contain any cell.", 
              msg_type = "STOP")
}


#################################################################
##    Install MCL
#################################################################
#' Install MCL
#'
#' This function installs the MCL (Markov Cluster) program, which is required for
#' some of the functions in the \code{\link{scigenex}} package. MCL is a cluster
#' algorithm that uses stochastic flow simulation to cluster graphs.
#'
#' @param force logical indicating whether to force installation even if MCL is
#'   already installed. Default is \code{FALSE}.
#'
#' @examples
#' @importFrom utils download.file
#' # Install MCL
#' install_mcl()
#' @export install_mcl
install_mcl <- function(force=FALSE){
  if (.Platform$OS.type == "windows") {
    print_msg("A unix-like OS is required to launch the MCL program.",
              msg_type = "ERROR")
  }else{
    if(nchar(Sys.which("mcl")) == 0 | force ){
      
      if(is.null(unlist(options()["scigenex_mcl_path"])) | force){
        dir_path <- file.path(path.expand('~'), ".scigenex")
        print_msg(paste0("Creating a path for mcl installation: ", 
                         dir_path), 
                  msg_type = "INFO")
        dir.create(dir_path, showWarnings = FALSE)
        setwd(dir_path)
        utils::download.file("http://micans.org/mcl/src/mcl-latest.tar.gz",
                      destfile="mcl-latest.tar.gz")
        system("tar xvfz mcl-latest.tar.gz")
        system("rm -f mcl-latest.tar.gz")
        mcl_version <- dir()
        setwd(mcl_version)
        print_msg("Installing MCL.", 
                  msg_type = "INFO")
        system("./configure")
        system("make")
        mcl_install_path <- file.path(getwd(), "src/shmcl/mcl")
        print_msg("MCL program installed.", 
                  msg_type = "INFO")
        options(scigenex_mcl_path = mcl_install_path)
      }
    }
  }
}
