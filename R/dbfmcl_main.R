################################################################
##        Main script of the R PACKAGE : DBFMCL
##
## Authors : BERGON A, J. BAVAIS
##  with the collaboration of LOPEZ F., TEXTORIS J. and PUTHIER D.
##
##
## R CMD SHLIB dbf.c -o dbf
##
#################################################################
##    DEFINITION OF A SPECIFIC CLASS OBJECT : DBFMCLresult
#################################################################
 
library(Biobase)
library(ggplot2)
library(reshape2)
library(dplyr)

#' Title
#' This class represents the results of the \code{\link{DBFMCL}} algorithm.
#' @slot name character. The original input file name (if applicable).
#' @slot data matrix. The matrix containing the filtered/partitionned data.
#' @slot cluster vector. Mapping of row/genes to clusters.
#' @slot size vector. The size of each cluster.
#' @slot center matrix. The centers of each clusters.
#' @slot parameters list. The parameter used.
#'
#' @return A DBFMCLresult object.
#' @export
#'
#' @examples
#' 
#' \dontrun{
#'   m <- matrix(rnorm(80000), nc=20)
#'   m[1:100,1:10] <- m[1:100,1:10] + 4
#'   m[101:200,11:20] <- m[101:200,11:20] + 3
#'   m[201:300,5:15] <- m[201:300,5:15] + -2
#'   res <- DBFMCL(data=m,
#'                 distance_method="pearson",
#'                 k=25)
#' is(res) 
#' }
#'               
setClass("DBFMCLresult",
  representation = list(
    name = "character",
    data = "matrix",
    cluster = "vector",
    size = "vector",
    center = "matrix",
    parameters = "list"
  ),
  prototype = list(
    name = character(),
    data = matrix(nr = 0, nc = 0),
    cluster = numeric(),
    size = numeric(),
    center = matrix(nc = 0, nr = 0),
    parameters = list()
  )
)

#################################################################
##    REDEFINE SHOW METHOD FOR CLASS OBJECT : DBFMCLresult
#################################################################

setMethod(
  "show", signature("DBFMCLresult"),

  function(object) {
    
    cat("\t\tAn object of class DBFMCLresult\n")
    cat("\t\tName:", slot(object, "name"), "\n")
    cat("\t\tMemory used: ", object.size(object), "\n")
    cat("\t\tNumber of samples: ", ncol(slot(object, "data")), "\n")
    cat(
      "\t\tNumber of informative genes: ",
      nrow(slot(object, "data")), "\n"
    )
    cat("\t\tNumber of clusters: ", length(slot(object, "size")), "\n")
    cat("\t\tThis object contains the following informations:\n")

    for (i in slotNames(object)) {
      cat("\t\t\t - ", i, "\n")
    }
    if (length(slot(object, "parameters")) > 0) {
      for (i in 1:length(slot(object, "parameters"))) {
        cat(
          "\t\t\t\t * ", names(slot(object, "parameters"))[[i]],
          " = ", slot(object, "parameters")[[i]], "\n"
        )
      }
    }
  }
)


#################################################################
##    REDEFINE dim/ncol/nrow METHOD FOR CLASS OBJECT : DBFMCLresult
#################################################################


setMethod(
  "ncol", signature("DBFMCLresult"),

  function(x) {

      return(ncol(x@data))

      }
)

setMethod(
  "nrow", signature("DBFMCLresult"),

  function(x) {

      return(nrow(x@data))

      }
)
setMethod(
  "dim", signature("DBFMCLresult"),

  function(x) {

      return(dim(x@data))

      }
)

#################################################################
##    Define the nclust function for class DBFMCLresult
#################################################################

#' nclust
#' returns the number of clusters contained in a DBFMCLresult object.
#' @param object A DBFMCLresult object.
#'
#' @return the number of clusters.
#' @export
#'
#' @examples
setGeneric("nclust",
           function(object) {
              standardGeneric("nclust")
})

#' @rdname nclust
setMethod(
  "nclust",
  signature(object = "DBFMCLresult"),
  function(object) {
    return(sort(unique(object@cluster)))
  }
)
#################################################################
##    Define the plot_profile for class DBFMCLresult
#################################################################

#' plot_dbf
#' Plot the results (heatmap or profiles) contained in a DBFMCLresult object.
#' @param object A DBFMCLresult object.
#' @param type The type of diagram.
#' @param to_log2 Whether data should be transform in logarithm base 2 (+ 1 as a pseudocount).
#' @param colo  r_palette A color palette.
#' @param standardizing Whether rows should be divided by standard deviation.
#' @param centering Whether rows should be centered. 
#'
#' @return A ggplot diagram.
#' @export
#'
#' @examples
setGeneric("plot_dbf", 
           
           function(object,
                    type = c("line", "tile"),
                    to_log2 = FALSE,
                    color_palette = NULL,
                    standardizing = FALSE,
                    centering = TRUE) {
             
              standardGeneric("plot_dbf")
})

#' @rdname plot_dbf
setMethod(
  "plot_dbf",
  signature(object = "DBFMCLresult"),
  function(object,
           type = c("line", "tile"),
           to_log2 = FALSE,
           color_palette = "#0000BF,#0000FF,#0080FF,#00FFFF,#40FFBF,#80FF80,#BFFF40,#FFFF00,#FF8000,#FF0000,#BF0000",
           standardizing = FALSE,
           centering = TRUE) {

    # The type of diagram
    type <- type[1]


    # color_palette_list = color_palette.split(",")
    # if len(color_palette_list) < 2:
    #    message("Need more than 2 colors for heatmap color palette.",
    #            type="ERROR")

    ## getting matrix
    nb <- length(object@size)
    m <- object@data

    if (to_log2) {
      m <- log2(m + 1)
    }

    ## median-centering of row
    if (centering) {
      mean_row <- apply(m, 1, mean)
      m <- sweep(m, MARGIN = 1, STATS = mean_row, FUN = "-")
    }

    ## median-centering of row
    if (standardizing) {
      sd_row <- apply(m, 1, sd)
      m <- sweep(m, MARGIN = 1, STATS = sd_row, FUN = "/")
    }

    ##  hclust on samples
    dis <- dist(t(m))
    h <- hclust(dis, method = "av")
    m <- m[, h$order]
    myOut <- ""


    ## melting
    m_melt <- as.data.frame(m)
    m_melt$cluster <- object@cluster
    m_melt$gene <- row.names(object@data)

    m_melt <- melt(m_melt,
      id.vars = c("cluster", "gene"),
      variable.name = "samples"
    )


    ## ploting
    col <- unlist(strsplit("#67001f,#b2182b,#d6604d,#f4a582,#fddbc7,#f7f7f7,#d1e5f0,#92c5de,#4393c3,#2166ac,#053061", ","))
    color.ramp <- colorRampPalette(col)(10)
    # Note that samples, value, gene, cluster
    # may appear as undefined variable to R check.
    # A workaround is to define them as NULL first...
    samples <- value <- gene <- cluster <- cluster_mean <- NULL
    if (type == "line") {
      p <- ggplot(data = m_melt, aes(
        x = samples,
        y = value
      ))
      p <- p + geom_line(color = "azure3", aes(group = gene), size=0.1)
      p <- p + theme_bw()
      p <- p + facet_grid(cluster ~ .)
      p <- p + geom_line(
        data = m_melt %>%
          group_by(cluster, samples) %>%
          summarise(cluster_mean = mean(value)),
        aes(
          x = samples,
          y = cluster_mean,
          group = cluster
        ),
        color = "skyblue4",
	size=0.8
      )

      p <- p + theme(
        strip.text.y = element_text(angle = 0),
        axis.text.x = element_blank()
      )
    } else if (type == "tile") {
      p <- ggplot(
        data = m_melt,
        aes(
          x = samples,
          y = gene,
          fill = value
        )
      )

      p <- p + geom_tile()
      p <- p + theme_bw()
      p <- p + scale_fill_gradientn(
        colours = color.ramp,
        name = "Signal"
      )
      p <- p + theme(
        strip.text.y = element_text(angle = 0),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank()
      )
    }

    p <- p + facet_grid(cluster ~ ., scales = "free_y")

    return(p)
  }
)

#################################################################
##    Define the write function for class DBFMCLresult
#################################################################

#' write_dbf
#' Write a DBFMCLresult into a flat file.
#' @param object DBFMCLresult. 
#' @param filename_out The outfile name.
#' @param path The path to the file.
#' @param nb_na_row Number of separating rows (containing NAs).
#' @param verbose Whether verbosity should be activated.
#' @return Write a file.
#' @export
#'
#' @examples
setGeneric("write_dbf", 
           
           function(object, 
                    filename_out = NULL,
                    path = ".",
                    nb_na_row=3,
                    verbose = TRUE) {
                standardGeneric("write_dbf")
})


#' @rdname write_dbf
setMethod(
  "write_dbf",
  signature(object = "DBFMCLresult"),
  function(object,
           filename_out = NULL,
           path = ".",
           nb_na_row=5,
           verbose = TRUE) {
           
    if (path == ".") path <- getwd()

    if (is.null(filename_out)) {
      filename_out <- "exprs.dataMods.txt"
    }

    data <- object@data
    nb <- 0
    dataT <- c("clusters", colnames(data))

    ## processing data
    for (i in 1:length(object@size)) {
      if (verbose) {
        cat("\n\tCluster ", i, " --> ", object@size[i], " probes")
      }

      subData <- data[object@cluster == i, ]
      subData <- cbind(rownames(subData), subData)
      if (nb_na_row > 0){
        intLine <- matrix(rep(NA, (ncol(data) + 1)*nb_na_row), nrow = nb_na_row)
        dataT <- rbind(dataT, subData, intLine)
      }
      nb <- nb + 1

    }

    ## exporting results
    write.table(dataT, file.path(path, filename_out),
      col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE
    )

    if (verbose) {
      cat("\n\t--> Creating file : ",
        file.path(path, filename_out), "\n\n")
    }
  }
)

#################################################################
##    DBF-MCL
#################################################################

#' The "Density Based Filtering and Markov CLustering" algorithm (DBF-MCL).
#'
#' DBF-MCL is a tree-steps adaptative algorithm that \emph{(i)} find elements
#' located in dense areas (DBF), \emph{(ii)} uses selected items to construct a
#' graph, \emph{(iii)} performs graph partitioning using the Markov CLustering
#' Algorithm (MCL).
#'
#' This function requires installation of the mcl program
#' (\url{http://www.micans.org/mcl}). See "Warnings" section for more
#' informations.
#'
#' When analyzing a noisy dataset, one is interested in isolating dense regions
#' as they are populated with genes/elements that display weak distances to
#' their nearest neighbors (i.e. strong profile similarities). To isolate these
#' regions DBF-MCL computes, for each gene/element, the distance with its kth
#' nearest neighbor (DKNN).In order to define a critical DKNN value that will
#' depend on the dataset and below which a gene/element will be considered as
#' falling in a dense area, DBF-MCL computes simulated DKNN values by using an
#' empirical randomization procedure. Given a dataset containing n genes and p
#' samples, a simulated DKNN value is obtained by sampling n distance values
#' from the gene-gene distance matrix D and by extracting the kth-smallest
#' value. This procedure is repeated n times to obtain a set of simulated DKNN
#' values S. Computed distributions of simulated DKNN are used to compute a FDR
#' value for each observed DKNN value. The critical value of DKNN is the one
#' for which a user-defined FDR value (typically 10\%) is observed. Genes with
#' DKNN value below this threshold are selected and used to construct a graph.
#' In this graph, edges are constructed between two genes (nodes) if one of
#' them belongs to the k-nearest neighbors of the other. Edges are weighted
#' based on the respective coefficient of correlation (\emph{i.e.}, similarity)
#' and the graph obtained is partitioned using the Markov CLustering Algorithm
#' (MCL).
#'
#' @param data a \code{matrix}, \code{data.frame} or \code{Seurat} object.
#' @param filename a character string representing the file name.
#' @param name a prefix for the names of the intermediary files created by DBF
#' and MCL.
#' @param path a character string representing the data directory where
#' intermediary files are to be stored. Default to current working directory.
#' @param distance_method a method to compute the distance to the k-th nearest
#' neighbor. One of "pearson" (Pearson's correlation coefficient-based
#' distance), "spearman" (Spearman's rho-based distance), "euclidean".
#' @param av_dot_prod_min Any cluster with average dot product below this value is discarded. This allow to delete
#' clusters in which correlation is influenced/supported by very few samples (typically 1).
#' @param min_cluster_size Minimum number of element inside a cluster. MCL tend to create lots of clusters with
#' very few (e.g 2) objects.
#' @param silent if set to TRUE, the progression of distance matrix calculation
#' is not displayed.
#' @param verbose if set to TRUE the function runs verbosely.
#' @param k the neighborhood size.
#' @param random the number of simulated distributions S to compute. By default
#' \code{random = 3}.
#' @param memory_used size of the memory used to store part of the distance
#' matrix. The subsequent sub-matrix is used to computed simulated distances to
#' the k-th nearest neighbor (see detail section).
#' @param fdr an integer value corresponding to the false discovery rate
#' (range: 0 to 100).
#' @param inflation the main control of MCL. Inflation affects cluster
#' granularity. It is usually chosen somewhere in the range \code{[1.2-5.0]}.
#' \code{inflation = 5.0} will tend to result in fine-grained clusterings, and
#' whereas \code{inflation = 1.2} will tend to result in very coarse grained
#' clusterings. By default, \code{inflation = 2.0}. Default setting gives very
#' good results for microarray data when k is set between 70 and 250.
#' @param set.seed specify seeds for random number generator.
#' @return a DBFMCLresults class object.
#' @section Warnings: With the current implementation, this function only works
#' only on UNIX-like plateforms.
#'
#' MCL should be installed. One can used the following command lines in a
#' terminal:
#'
#' \code{# Download the latest version of mcl (the script has been tested
#' successfully with the 06-058 version).}
#' \code{wget http://micans.org/mcl/src/mcl-latest.tar.gz}
#' \code{# Uncompress and install mcl}
#' \code{tar xvfz mcl-latest.tar.gz}
#' \code{cd mcl-xx-xxx}
#' \code{./configure}
#' \code{make}
#' \code{sudo make install}
#' \code{# You should get mcl in your path}
#' \code{mcl -h}
#' @author Bergon A., Bavais J., Textoris J., Granjeaud S., Lopez F and Puthier
#' D.
#' @references Van Dongen S. (2000) A cluster algorithm for graphs. National
#' Research Institute for Mathematics and Computer Science in the 1386-3681.
#' @keywords manip
#' @examples
#'
#' \dontrun{
#' ## with an artificial dataset
#'
#' m <- matrix(rnorm(80000), nc = 20)
#' m[1:100, 1:10] <- m[1:100, 1:10] + 4
#' m[101:200, 11:20] <- m[101:200, 11:20] + 3
#' m[201:300, 5:15] <- m[201:300, 5:15] + -2
#'   res <- DBFMCL(data=m,
#'                 distance_method="pearson",
#'                 k=25)
#' plot_profile(res)
#' write_clusters(res, filename_out = "ALL.sign.txt")
#' }
#'
#' @export DBFMCL
DBFMCL <- function(data = NULL, 
                   filename = NULL, 
                   path = ".", 
                   name = NULL,
                   distance_method = c("pearson", "spearman", "euclidean"),
                   av_dot_prod_min=2,
                   min_cluster_size=10,
                   silent = FALSE, 
                   verbose = TRUE, 
                   k = 150,
                   random = 3, 
                   memory_used = 1024,
                   fdr = 10, 
                   inflation = 2.0,
                   set.seed = 123) {

  ## testing the system
  if (.Platform$OS.type == "windows") {
    stop("\t--> A unix-like OS is required to launch mcl and cluster programs.")
  }

  ## getting parameters
  data_source <- get_data_4_DBFMCL(data = data, filename = filename, path = path)
  m <- data_source$data

  # A simple function to create a random string
  create_rand_str <- function() {
      v = c(sample(LETTERS, 3, replace = TRUE),
            sample(0:9, 4, replace = TRUE),
            sample(letters, 3, replace = TRUE))
      return(paste0(sample(v),collapse = ""))
  }

  if (is.null(name)) name <- data_source$name
  if (is.null(name)) name <- create_rand_str()

  distance_method <- match.arg(distance_method)
  txt <- paste("\n\tInflation: ", inflation, sep = "")

  ## writting all parameters
  if (verbose) {
    cat(
      "The following parameters will be used :",
      "\n\tWorking directory: ", getwd(),
      "\n\tName: ", name,
      "\n\tDistance method: ", distance_method,
      "\n\tMinimum average dot product for clusters: ", av_dot_prod_min,
      "\n\tMinimum cluster size: ", min_cluster_size,
      "\n\tNumber of neighbors: ", k,
      "\n\tNumber of randomizations: ", random,
      "\n\tFDR: ", fdr, "%", txt,
      "\n\tVisualize standard outputs from both mcl and cluster",
      "commands: ", silent,
      "\n\tMemory used : ", memory_used, "\n\n"
    )
  }

  ## DBF algorithm, returns a DBFMCLresult object
  obj <- DBF(m,
             name,
             distance_method = distance_method,
             silent = silent,
             k = k,
             random = random,
             memory_used = memory_used,
             fdr = fdr,
             set.seed = set.seed
  )

  if (length(readLines(paste(name, ".dbf_out.txt", sep = ""))) > 0) {

      ## Launching mcl
      if (is.null(inflation)) inflation <- 2.0
      MCL(name, inflation = inflation, silent = silent)

      ## getting mcl results into the DBFMCLresult object
      mcl_cluster <- readLines(paste(name, ".mcl_out.txt", sep = ""))
      gene_list <- NULL
      clusters <- NULL
      size <- NULL
      nb <- 0
      nb_cluster_deleted <- 0

      for (i in 1:length(mcl_cluster)) {
        h <- unlist(strsplit(mcl_cluster[i], "\t"))
        cur_clust <- m[h,]
        cur_clust[cur_clust > 0 ] <- 1
        cur_dot_prod <- cur_clust %*% t(cur_clust)

        if(mean(cur_dot_prod) > av_dot_prod_min & length(h) > min_cluster_size){

          nb <- nb + 1
          gene_list <- c(gene_list, h)
          clusters <- c(clusters, rep(nb, length(h)))
          if (is.null(size)) {
            size <- length(h)
          }
          else {
            size <- c(size, length(h))
          }
        }else{
          nb_cluster_deleted <- nb_cluster_deleted + 1
        }
      }
      if (verbose) {
        cat(
          "\t--> ", nb, " clusters conserved after MCL partitioning\n\n"
        )
        cat(
          "\t--> ", nb_cluster_deleted, " clusters deleted after MCL partitioning\n\n"
        )
      }

      ## build DBFMCLresult object
      if (nb > 0) {
        obj@name <- name
        obj@data <- as.matrix(m[gene_list, ])
        obj@cluster <- clusters
        obj@size <- size

        centers <- matrix(ncol = ncol(m), nrow = nb)
        ## calcul of the mean profils
        for (i in 1:nb) {
          centers[i, ] <- apply(obj@data[obj@cluster == i, ],
            2, mean,
            na.rm = TRUE
          )
        }
        obj@center <- centers

        ## add DBFMCL parameters used to build this object
        obj@parameters <- list(
          distance_method = distance_method,
          k = k,
          random = random,
          fdr = fdr,
          set.seed = set.seed,
          inflation = inflation
        )
      }
  }
  else {
    stop("\t--> There is no conserved gene.\n\n")
  }
  return(obj)

}


###############################################################
##    COMPUTE DBF algorithm
###############################################################

#' Density-Based Filtering.
#'
#' This function is an internal function used by \code{\link{DBFMCL}} to detect
#' informative elements (\emph{i.e.}, those that belong to dense regions). User
#' should not use this function. Instead they can use the \code{\link{DBFMCL}}
#' function with \code{clustering} argument set to \code{FALSE}.
#'
#' See \code{\link{DBFMCL}}
#'
#' @param data a matrix or data.frame
#' @param name a prefix for the file name
#' @param distance_method a method to compute the distance to the k-th nearest
#' neighbor. One of "pearson" (Pearson's correlation coefficient-based
#' distance), "spearman" (Spearman's rho-based distance) or "euclidean".
#' @param silent if set to TRUE (default), the progression of distance matrix
#' calculation is not displayed.
#' @param k the neighborhood size.
#' @param random the number of simulated distributions S to compute. By default
#' \code{random = FALSE}.
#' @param fdr a value for the false discovery rate.
#' @param memory_used size of the memory used to store part of the distance
#' matrix. The subsequent sub-matrix is used to computed simulated distances to
#' the k-th nearest neighbor (see detail section).
#' @param set.seed specify seeds for random number generator.
#' @section Warnings: Works only on UNIX-alikes platforms.
#' @author Bergon A., Bavais J., Textoris J., Granjeaud S., Lopez F and Puthier
#' D.
#' @seealso \code{\link{DBFMCL}}
#' @references Lopez F.,Textoris J., Bergon A., Didier G., Remy E., Granjeaud
#' S., Imbert J. , Nguyen C. and Puthier D. TranscriptomeBrowser: a powerful
#' and flexible toolbox to explore productively the transcriptional landscape
#' of the Gene Expression Omnibus database. PLoSONE, 2008;3(12):e4001.
#' @keywords manip
#' @export DBF
DBF <- function(data, name = NULL,
                distance_method = c("spearman", "pearson", "euclidean"),
                silent = FALSE,
                k = 100,
                random = 3,
                fdr = 10,
                memory_used = 1024,
                set.seed = 123) {

  ## testing the system
  if (.Platform$OS.type != "windows") {
    if (!is.null(data)) {
      ## getting data and parameters
      if (is.null(name)) name <- "exprs"
      data <- get_data_4_DBFMCL(data = data)$data
      row <- rownames(data)
      col <- colnames(data)
      distance_method <- match.arg(distance_method)

      ## transforming data into double
      data <- apply(data, 2, as.double)

      if (silent) {
        cat(
          "\t--> Computing distances to the kth-nearest neighbors",
          " and associated FDR values... \n"
        )
      }
      outfile <- paste(name, ".dbf_out.txt", sep = "")

      ## launching DBF
      a <- .C("DBF",
        data,
        as.integer(nrow(data)),
        as.integer(ncol(data)),
        row,
        col,
        distance_method,
        as.integer(k),
        as.integer(random),
        as.integer(!silent),
        as.integer(memory_used),
        as.integer(fdr),
        as.integer(!silent),
        m2 = vector(length = nrow(data), mode = "character"),
        outfile,
        as.integer(set.seed)
      )

      ## creation of the DBFMCLresult object
      obj <- new("DBFMCLresult")

      informative <- a$m2[a$m2 != ""]
      if (length(informative) > 0) {
        obj@data <- as.matrix(data[a[[4]] %in% informative, ])
        obj@cluster <- rep(1, nrow(obj@data))
        obj@size <- nrow(obj@data)
        obj@center <- matrix(
          apply(obj@data[obj@cluster == 1, ],
            2,
            mean,
            na.rm = TRUE
          ),
          nrow = 1
        )
      }
      return(obj)
    } else {
      stop("\t--> Please provide a matrix...\n\n")
    }
  }
  else {
    stop("\t--> A unix-like OS is required to launch mcl and cluster programs.")
  }
}

##############################################################
##    MCL
##############################################################

#' Invokes the Markov CLustering algorithm (MCL).
#'
#' This function invokes the mcl system command. MCL is a clustering algorithm
#' for graphs that was developped by Stijn van Dongen (see references for
#' further informations).
#'
#'
#' @param name a character string corresponding to the file name.
#' @param inflation the main control of MCL. Inflation affects cluster
#' granularity. It is usually chosen somewhere in the range \code{[1.2-5.0]}.
#' \code{inflation = 5.0} will tend to result in fine-grained clusterings, and
#' whereas \code{inflation = 1.2} will tend to result in very coarse grained
#' clusterings. By default, \code{inflation = 2.0}. Default setting gives very
#' good results for microarray data when k is set around 100.
#' @param silent if set to TRUE, the progression of the MCL partitionning is
#' not displayed.
#' @return Returns a file with the ".mcl\_out.txt" extension.
#' @section warning: Works only on UNIX-like plateforms. MCL should be
#' installed. The following command lines can be used for installation.
#'
#' \code{# Download the latest version of mcl (RTools4TB has been tested
#' successfully with the 06-058 version).}
#' \code{wget http://micans.org/mcl/src/mcl-latest.tar.gz}
#' \code{# Uncompress and install mcl}
#' \code{tar xvfz mcl-latest.tar.gz}
#' \code{cd mcl-xx-xxx}
#' \code{./configure}
#' \code{make}
#' \code{sudo make install}
#' \code{# You should get mcl in your path}
#' \code{mcl -h}
#' @author Bergon A., Lopez F., Textoris J., Granjeaud S. and Puthier D.
#' @references Stijn van Dongen. A cluster algorithm for graphs.  Technical
#' Report INS-R0010, National Research Institute for Mathematics and Computer
#' Science in the Netherlands, Amsterdam, May 2000.
#' \url{http://www.cwi.nl/ftp/CWIreports/INS/INS-R0010.ps.Z}
#' @keywords manip
#' @export MCL
MCL <- function(name, inflation = 2.0, silent = FALSE) {
  ## testing the system
  if (.Platform$OS.type != "windows") {

    ## Testing mcl installation
    if (system("mcl --version | grep 'Stijn van Dongen'", intern = TRUE) > 0) {
      if (!silent) {
        cat("Running mcl (graph partitioning)... ")
        verb <- ""
      }
      else {
        verb <- "-V all"
      }
      if (inflation != 2) {
        i <- paste("-I ", as.character(round(inflation, 1)), sep = "")
      } else {
        i <- "-I 2.0"
      }
      ## launching mcl program
      system(paste("mcl ", name, ".dbf_out.txt ", i, " --abc -o ",
        name, ".mcl_out.txt ", verb,
        sep = ""
      ))

      if (!silent) {
        cat("\t-->  Done.\n\n")
        cat(
          "\t--> creating file : ",
          file.path(getwd(), paste(name, ".mcl_out.txt", sep = "")),
          "\n\n"
        )
      }
    } else {
      stop(
        "\t--> Please install mcl on your computer...\n",
        "\t--> You can download it from : 'http://www.micans.org/mcl/'\n\n"
      )
    }
  }
  else {
    stop("A unix-like OS is required to launch mcl and cluster programs.")
  }
}


#################################################################
##    Getting data
#################################################################

#' Title
#' Fetch an expression matrix from a file, dataframe or Seurat object.
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
      data <- as.matrix(data@assays$RNA@data)
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
    rownames(data) <- paste("gene", 1:nrow(data), sep = "")
  }
  if (is.null(colnames(data))) {
    colnames(data) <- paste("sample", 1:ncol(data), sep = "")
  }

  return(list(data = data, name = name))
}
#########################################################
##      END PACKAGE DBFMCL
#########################################################
