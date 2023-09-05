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
#' options(warn=-1)
#' print_msg("A warning message not displayed", "WARNING")
#' options(warn=opt_warn)
#' @keywords internal
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
#' @title Debugging statistics (vector, matrix or dataframe)
#' @description
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
#' @keywords internal
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
#' @title Generate a random string of letters and numbers
#' @description
#' This function generates a random string of 10 characters, consisting of 3 uppercase letters,
#' 4 digits, and 3 lowercase letters.
#'
#' @returns A character string of length 10, consisting of random letters and numbers
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
#' @title Generate an example dataset with three clusters of profiles 
#'
#' @description 
#' Generate an artificial matrix with with random noise and 3 
#' clusters of 'expression' profiles (as row). 
#'
#' @returns a matrix.
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
#' @title Generate an example dataset with four clusters of profiles 
#' @description
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
#' @title Generate a vector of colors for a gradient
#'
#' @description 
#' This function generates a vector of colors for a gradient, given 
#' a specified palette name.
#'
#' @param palette A character vector specifying the palettes to use. One of: "Je1", 
#' "Seurat_like", "Ju1", "De1",  "De2", "De3", "De4", "De5", "De6", "De7", "De8", "De9", 
#' "Magma", "viridis", "magma2", "plasma"
#' @return A character vector of color codes.
#' @export colors_for_gradient
#' @examples
#' colors_for_gradient()
#' colors_for_gradient(palette = "Ju1")
#' 
colors_for_gradient <- function(palette=c("Je1", "Seurat_Like", "Ju1", "De1", 
                                          "De2", "De3", "De4", "De5",
                                          "De6", "De7", "De8", "De9", "Magma", 
                                          "viridis", "magma2", "plasma")){
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
  }else if(palette == "Magma"){
    return(c("#ffdb00", "#ffa904", "#ee7b06", "#a12424", "#400b0b"))
  }  else if(palette == "viridis"){
    return(c("#fde725", "#5ec962", "#21918c", "#3b528b", "#440154"))
  } else if(palette == "magma2"){
    return(c("#fcfdbf", "#fc8961", "#b73779", "#51127c", "#000004"))
  } else if(palette == "plasma"){
    return(c("#f0f921", "#f89540", "#cc4778", "#7e03a8", "#0d0887"))
  }
  
}

# -------------------------------------------------------------------------
# Returns a discrete color palette ----------------------------------------
# -------------------------------------------------------------------------
#' @title Generate a discrete color palette
#'
#' @description
#' This function generates a vector of colors for a discrete variable, 
#' given a specified palette name.
#'
#' @param n An integer specifying the number of colors to generate. 
#' @param palette A character vector specifying the palette to use. 
#' @return A character vector of color codes.
#' @export discrete_palette
#' @importFrom grDevices hcl
#' @examples
#' discrete_palette()
#' discrete_palette(n=20, palette = "ggplot")
#' 
discrete_palette <- function(n=10, palette=c("Ju1", "De1", "ggplot")){
  
  palette <- match.arg(palette)
  
  if(palette == "Ju1"){
    palette <- colorRampPalette(c(  "#9F1717", "#AE5B11", "#C48D00", "#517416", 
                                    "#115C8A", "#584178", "#9D1C70", "#E96767", 
                                    "#EC9342", "#FFCA3A", "#8AC926", "#4DADE8", 
                                    "#9579B9", "#E25CB4", "#DB2020", "#DA7316", 
                                    "#F0AE00", "#6D9D1E", "#1882C0", "#71529A", 
                                    "#D02494", "#EF9292", "#F2B57D", "#FFDA77", 
                                    "#B6E36A", "#7BC4EE", "#AD98C9", "#EA8AC9"))(n)
  }else if(palette == "De1"){
    
    palette <- colorRampPalette(c("grey","#EA8AC9", "#9F1717"))(n)
    
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

# -------------------------------------------------------------------------
# Check the format of a Clusterset object ---------------------------------
# ------------------------------------------------------------------------
#' @title Check the format of a Clusterset object
#' @description
#' Check the format of a Clusterset object (ncol, nrow, inherits...)
#'
#' @param object the clusterSet object to be tested 
#' @returns None.
#' @examples 
#' # load a dataset
#' load_example_dataset('7871581/files/pbmc3k_medium_clusters')
#' check_format_cluster_set(pbmc3k_medium_clusters)
#' @export check_format_cluster_set
#' @keywords internal
check_format_cluster_set <- function(object) {
  
  print_msg("Checking ClusterSet format.", msg_type = "DEBUG")
  
  if(is.null(object)){
    print_msg("The ClusterSet is a NULL object.", 
              msg_type = "STOP")
  }
  
  if(!inherits(object, "ClusterSet"))
    print_msg("Please provide a ClusterSet as input.", 
              msg_type = "STOP")
  
  if(nrow(object) == 0)
    print_msg("The ClusterSet object does not contain any gene.", 
              msg_type = "WARNING")
  
  if(ncol(object) == 0)
    print_msg("The ClusterSet object does not contain any cell.", 
              msg_type = "WARNING")
}

# -------------------------------------------------------------------------
# Install MCL -------------------------------------------------------------
# -------------------------------------------------------------------------
#' @title Install MCL
#' @description
#' This function installs the MCL (Markov Cluster) program, which is required for
#' the main functions of the scigenex package. MCL is a cluster
#' algorithm that uses stochastic flow simulation to cluster graphs.
#'
#' @importFrom utils download.file
#' @examples
#' # Install MCL
#' install_mcl()
#' @export install_mcl
install_mcl <- function(){
  
  
  if (.Platform$OS.type == "windows") {
    print_msg("A unix-like OS is required to launch the MCL program.",
              msg_type = "STOP")
  }else{
        old_path <- getwd()
        dir_path <- file.path(path.expand('~'), ".scigenex")
        print_msg(paste0("Creating a path for mcl installation: ", 
                         dir_path), 
                  msg_type = "INFO")
        dir.create(dir_path, showWarnings = FALSE, recursive = TRUE)
        setwd(dir_path)
        utils::download.file("http://micans.org/mcl/src/mcl-latest.tar.gz",
                             destfile="mcl-latest.tar.gz")
        system("tar xvfz mcl-latest.tar.gz")
        system("rm -f mcl-latest.tar.gz")
        mcl_version <- dir(pattern = "mcl*")
        setwd(mcl_version)
        print_msg("Installing MCL.", msg_type = "INFO")
        system("sh configure; make")
        
        print_msg("MCL program installed.", 
                  msg_type = "INFO")
      
        setwd(old_path)
  }
  
}

# -------------------------------------------------------------------------
# Show available methods --------------------------------------------------
# -------------------------------------------------------------------------
#' @title List the method for the ClusterSet object.
#' @description
#' List the method for the ClusterSet object.
#' @param class The name of the class. Default to "ClusterSet".
#' @param where The package name. Default to "scigenex".
#' @returns A vector with the available methods.
#' @importFrom methods showMethods
#' @importFrom utils capture.output
#' @examples
#' show_methods()
#' @export show_methods
show_methods <- function(class="ClusterSet",
                         where="scigenex"){
  class_method <- utils::capture.output(methods::showMethods(class="ClusterSet", 
                                                      where = paste0("package:", where)))
  
  class_method <- class_method[grep("Function", class_method, perl=T)]
  class_method <- gsub("Function: ([^ ]+) \\(.*", "\\1", class_method)
  
  return(class_method)
}

# -------------------------------------------------------------------------
# load_example_dataset()   ------------------------------------------------
# -------------------------------------------------------------------------
#' @title Load/download a Seurat or ClusterSet example dataset.
#' @description
#' Load/download a Seurat or ClusterSet example dataset.
#' @param dataset The name of the dataset.
#' @param timeout Set the timout for download (options(timeout=timeout))
#' @param force Reload the dataset even if already in globalEnv. 
#' @returns Load the dataset.
#' @examples
#' # An example Seurat/Visium dataset
#' load_example_dataset("7870305/files/lymph_node_tiny_2")
#' lymph_node_tiny_2
#' # An example clusterSet dataset
#' load_example_dataset("7870305/files/lymph_node_tiny_clusters_2")
#' lymph_node_tiny_clusters_2
#' @export load_example_dataset
load_example_dataset <- function(dataset=c("7871581/files/pbmc3k_medium",
                                           "7871581/files/pbmc3k_medium_clusters",
                                           "8028126/files/pbmc3k_medium_clusters_enr",
                                           "8028226/files/pbmc3k_medium_clusters_enr_sub",
                                           "7870305/files/lymph_node_tiny_2", 
                                           "7870305/files/lymph_node_tiny_clusters_2",
                                           "7869307/files/lymph_node_tiny",
                                         "7869307/files/lymph_node_tiny_clusters"),
                                 timeout=NULL,
                                 force=FALSE){
  
  dataset <- match.arg(dataset)
  
  if(!is.null(timeout))
    options(timeout=timeout)
  
  file_data <- gsub(".*\\/", "", dataset)
  dir_path <- file.path(path.expand('~'), ".scigenex", "datasets")
  
  if(!dir.exists(dir_path)){
    print_msg(paste0("Creating a path for dataset installation: ", 
                     dir_path), 
              msg_type = "INFO")
    dir.create(dir_path, showWarnings = FALSE, recursive = TRUE) 

  }
  
  old_path <- getwd()
  setwd(dir_path)
  
  if(!file.exists(paste0(file_data, ".rda"))){
    download.file(url=paste0("https://zenodo.org/record/", 
                             dataset, ".rda"), 
                  destfile = paste0(file_data, ".rda"))
  }

  dataset_short <- gsub(".*/", "", dataset)
  
  if(!dataset_short %in% ls(envir = globalenv()) | force){
    load(file.path(dir_path, paste0(file_data, ".rda")), envir = .GlobalEnv)
    print_msg(paste0("Dataset ", dataset, " has been loaded."),
              msg_type = "INFO")
    
  }else{
    print_msg(paste0("Dataset ", dataset, " was already loaded."),
              msg_type = "INFO")
  }
  
  setwd(old_path)

}

