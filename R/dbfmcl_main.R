################################################################
##        Main script of the R PACKAGE : scigenex
##
## Authors : J. BAVAIS, BERGON A, 
##  with the collaboration of LOPEZ F., TEXTORIS J. and PUTHIER D.
##
##
#################################################################
##    UTILS (a set of useful function)
#################################################################
# 0 : No message
# 1 : info
# 2 : DEBUG

VERBOSITY_SCIGENEX = 3

print_msg <- function(msg, msg_type="INFO"){
  if(msg_type == "INFO")
    if(VERBOSITY_SCIGENEX > 0)
      cat(paste("|-- ", msg, "\n"))
  if(msg_type == "DEBUG")
    if(VERBOSITY_SCIGENEX > 1)
      cat(paste("|-- ", msg, "\n"))
  if(msg_type == "WARNING")
      cat(paste("|-- ", msg, "\n"))
}



#################################################################
##    DBF-MCL
#################################################################

#' @title
#' The "Density Based Filtering and Markov CLustering" algorithm (DBF-MCL).
#' @description
#' DBFMCL is a tree-steps adaptative algorithm that \emph{(i)} find elements
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
#' @param output_path a character string representing the data directory where
#' output files will be stored. Default to current working directory.
#' @param mcl_threads An integer to determine number of threads for mcl algorithm.
#' @param distance_method a method to compute the distance to the k-th nearest neighbor.
#' @param av_dot_prod_min Any cluster with average dot product below this value is discarded. This allow to delete
#' clusters in which correlation is influenced/supported by very few samples (typically 1).
#' @param min_cluster_size Minimum number of element inside a cluster. MCL tend to create lots of clusters with
#' very few (e.g 2) objects.
#' @param silent if set to TRUE, the progression of distance matrix calculation
#' is not displayed.
#' @param k the neighborhood size.
#' @param fdr an integer value corresponding to the false discovery rate
#' (range: 0 to 100).
#' @param inflation the main control of MCL. Inflation affects cluster
#' granularity. It is usually chosen somewhere in the range \code{[1.2-5.0]}.
#' \code{inflation = 5.0} will tend to result in fine-grained clusterings, and
#' whereas \code{inflation = 1.2} will tend to result in very coarse grained
#' clusterings. By default, \code{inflation = 2.0}. Default setting gives very
#' good results for microarray data when k is set between 70 and 250.
#' @param seed specify seeds for random number generator.
#' @return a ClusterSets class object.
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
#' @references
#' - Van Dongen S. (2000) A cluster algorithm for graphs. National
#' Research Institute for Mathematics and Computer Science in the 1386-3681.
#' - Lopez F.,Textoris J., Bergon A., Didier G., Remy E., Granjeaud
#' S., Imbert J. , Nguyen C. and Puthier D. TranscriptomeBrowser: a powerful
#' and flexible toolbox to explore productively the transcriptional landscape
#' of the Gene Expression Omnibus database. PLoSONE, 2008;3(12):e4001.
#' @keywords clustering, gene expression, classification, MCL.
#' @examples
#'
#' \dontrun{
#'   m <- matrix(rnorm(80000), nc=20)
#'   m[1:100,1:10] <- m[1:100,1:10] + 4
#'    m[101:200,11:20] <- m[101:200,11:20] + 3
#'    m[201:300,5:15] <- m[201:300,5:15] + -2
#'    res <- DBFMCL(data=m,
#'                  distance_method="pearson",
#'                  av_dot_prod_min = 0,
#'                  inflation = 1.2,
#'                  k=25,
#'                  fdr = 10)
#' plot_clust(res, ceil = 10, floor = -10)
#' plot_clust(res, type="tile", ceil = 10, floor = -10)
#' write_clust(res, filename_out = "ALL.sign.txt")
#' }
#'
#' @export DBFMCL
DBFMCL <- function(data = NULL, 
                   filename = NULL, 
                   path = ".",
                   output_path = ".",
                   mcl_threads=1,
                   name = NULL,
                   distance_method = "pearson",
                   av_dot_prod_min=2,
                   min_cluster_size=10,
                   silent = FALSE,
                   k = 50,
                   fdr = 10,
                   inflation = 8,
                   seed = 123) {
  
  ## testing the system
  if (.Platform$OS.type == "windows") {
    stop("\t--> A unix-like OS is required to launch mcl and cluster programs.")
  }
  
  ## getting parameters
  data_source <- get_data_4_DBFMCL(data = data, filename = filename, path = path)
  data_matrix <- data_source$data
  
  # A simple function to create a random string
  create_rand_str <- function() {
    v = c(sample(LETTERS, 3, replace = TRUE),
          sample(0:9, 4, replace = TRUE),
          sample(letters, 3, replace = TRUE))
    return(paste0(sample(v),collapse = ""))
  }
  
  if (is.null(name)) name <- data_source$name
  if (is.null(name)) name <- create_rand_str()
  
  # Put the current working directory in output_path or path
  if(path == ".") {
    path <- getwd()
  }
  if(output_path == ".") {
    output_path <- getwd()
  }
  # Check if output directory exists. If not stop the command.
  if(!file.exists(output_path)){
    stop("Output directory provided does not exist.")
  }
  
  
  
  distance_method <- match.arg(distance_method)
  txt <- paste("\n\tInflation: ", inflation, sep = "")
  
  ## writting all parameters
  
  cat(
    "The following parameters will be used :",
    "\n\tWorking directory: ", path,
    "\n\tOuput directory: ", output_path,
    "\n\tName: ", name,
    "\n\tDistance method: ", distance_method,
    "\n\tMinimum average dot product for clusters: ", av_dot_prod_min,
    "\n\tMinimum cluster size: ", min_cluster_size,
    "\n\tNumber of neighbors: ", k,
    "\n\tFDR: ", fdr, "%", txt,
    "\n\tVisualize standard outputs from both mcl and cluster",
    "commands: ", silent,"\n\n"
  )
  
  
  ## DBF algorithm, returns a ClusterSet object
  obj <- DBF(data_matrix,
             output_path = output_path,
             name,
             distance_method = distance_method,
             silent = silent,
             k = k,
             fdr = fdr,
             seed = seed
  )
  
  dbf_out_file <- paste0(output_path, "/", name, ".dbf_out.txt")
  dbf_out_file <- gsub(pattern = "//", replacement = "/", x = dbf_out_file)
  mcl_out_file <- paste0(output_path, "/", name, ".mcl_out.txt")
  mcl_out_file <- gsub(pattern = "//", replacement = "/", x = mcl_out_file)
  
  print_msg("DBF completed. Starting MCL step.", msg_type="DEBUG")
  
  if (length(readLines(dbf_out_file)) > 0) {
    
    ## Launching mcl (command line)
    mcl_system_cmd(name, inflation = inflation, input_path = output_path, silent = silent, threads = mcl_threads)
    
    
    print_msg(paste0("Reading and filtering MCL output: ", mcl_out_file), msg_type="DEBUG")
    
    ## getting mcl results into the ClusterSet object
    mcl_cluster <- readLines(mcl_out_file)
    gene_list <- NULL
    clusters <- NULL
    size <- NULL
    nb <- 0
    nb_cluster_deleted <- 0
    median_cur_dot_prod <- c()
    
    for (i in 1:length(mcl_cluster)) {
      h <- unlist(strsplit(mcl_cluster[i], "\t"))
      cur_clust <- data_matrix[h,]
      cur_clust[cur_clust > 0 ] <- 1
      cur_dot_prod <- cur_clust %*% t(cur_clust)
      diag(cur_dot_prod) <- NA
      cur_dot_prod_median_of_max <-median(apply(cur_dot_prod, 1, max, na.rm=T))
      
      if(cur_dot_prod_median_of_max >= av_dot_prod_min & length(h) > min_cluster_size){
        # Extract median value of dot product for gene signature i
        median_cur_dot_prod[i] <- median(cur_dot_prod)
        
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
    print_msg(paste(nb, " clusters conserved after MCL partitioning."),
              msg_type="INFO")
    print_msg(paste(nb_cluster_deleted,
                    " clusters filtered out from MCL partitioning (size and mean dot product)."),
              msg_type="INFO")
    
    
    ## build ClusterSet object
    if (nb > 0) {
      obj@name <- name
      obj@data <- as.matrix(data_matrix[gene_list, ])
      names(clusters) <- rownames(obj@data)
      obj@gene_patterns <- clusters
      obj@size <- size
      
      centers <- matrix(ncol = ncol(data_matrix), nrow = nb)
      ## calcul of the mean profils
      for (i in 1:nb) {
        centers[i, ] <- apply(obj@data[obj@gene_patterns == i, ],
                              2, mean,
                              na.rm = TRUE
        )
      }
      obj@center <- centers
      
      # Add median values of dot product for each gene cluster
      names(median_cur_dot_prod) <- paste0("cluster_", seq(1:length(size)))
      obj@dot_prodcut <- median_cur_dot_prod
      
      ## add DBFMCL parameters used to build this object
      obj@parameters <- list(
        distance_method = distance_method,
        k = k,
        fdr = fdr,
        seed = seed,
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

#' @title
#' DBF
#' @description
#' This function is an internal function used by \code{\link{DBFMCL}} to detect
#' informative elements (\emph{i.e.}, those that belong to dense regions). User
#' should not use this function. Instead they can use the \code{\link{DBFMCL}}
#' function with \code{clustering} argument set to \code{FALSE}.
#'
#' See \code{\link{DBFMCL}}
#'
#' @param data a matrix or data.frame
#' @param output_path a character string representing the data directory where
#' output files will be stored. Default to current working directory.
#' @param name a prefix for the file name
#' @param distance_method a method to compute the distance to the k-th nearest
#' neighbor.
#' @param silent if set to TRUE (default), the progression of distance matrix
#' calculation is not displayed.
#' @param k the neighborhood size.
#' @param fdr a value for the false discovery rate.
#' @param seed specify seeds for random number generator.
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
DBF <- function(data,
                output_path = ".",
                name = NULL,
                distance_method = "pearson",
                silent = FALSE,
                k = 100,
                fdr = 10,
                seed = 123) {

  ## testing the system
  if (.Platform$OS.type != "windows") {
    if (!is.null(data)) {
      ## set a seed for reproductibility
      if (!is.null(seed)) {
        set.seed(seed)
      }
      
      ## getting data and parameters
      if (output_path == ".") {output_path <- getwd()}
      if (is.null(name)) name <- "exprs"
      data <- get_data_4_DBFMCL(data = data)$data
      row <- rownames(data)
      col <- colnames(data)
      distance_method <- match.arg(distance_method)

      if (silent) {
        print_msg(paste0("Computing distances to the kth-nearest neighbors ",
                         "and associated FDR values... \n"),
                  msg_type = "INFO")

      }
      
      # Directory and name of the principal output
      outfile <- paste(output_path, "/", name, ".dbf_out.txt", sep = "")
      outfile <- gsub(pattern = "//", replacement = "/", x = outfile)
      
      
      
      #################### Correlation and distance matrices
      # Remove genes with 0 values for all cells
      genes_to_rm <- c(which(rowSums(data) == 0))
      if (length(genes_to_rm) > 0) {
        select_for_correlation <- data[-c(which(rowSums(data) == 0)),]
      } else {
        select_for_correlation <- data
      }
      
      # Compute gene-gene pearson correlation matrix
      if(distance_method == "pearson") {
        cor_matrix <- corSparse(t(select_for_correlation))
        rownames(cor_matrix) <- rownames(select_for_correlation)
        colnames(cor_matrix) <- rownames(select_for_correlation)
      } else {
        print_msg(msg_type = "WARNING",
                  msg = "Distance used is not provided.")
      }
      
      # Store correlation matrix in input_mcl. It will be used to create the input file for mcl.
      input_mcl <- cor_matrix
      
      # Remove the diagonal and the upper triangular matrix (values is replaced by NA)
      cor_matrix[lower.tri(cor_matrix, diag=TRUE)] <- NA
      
      # Transform correlation matrix into distance matrix (values between 0 and 2)
      dist_matrix <- 2 - (cor_matrix + 1)
      
      
      #Create a dataframe with a column that contains all the gene ID
      df_dknn <- data.frame("gene_id" = rownames(dist_matrix))
      l_knn <- list()
      
      #################### DKNN for each genes
      # Extract the DKNN for each gene
      for (gene in df_dknn[,"gene_id"]){
        #Create a vector with all the correlation values for one gene(i)
        gene_dist <- c(subset(dist_matrix[gene,], !is.na(dist_matrix[gene,])),
                       subset(dist_matrix[,gene], !is.na(dist_matrix[,gene]))) 
        
        #Reorder the pearson correlation values (increasing order)
        row_dknn <- order(gene_dist, decreasing=F)[1:k]
        gene_dknn <- gene_dist[row_dknn]
        
        #Store the results in a list
        l_knn[[gene]] <- gene_dknn
        
        #Select the kth pearson correlation values. This value corresponds to the DKNN of the gene(i)
        df_dknn[which(df_dknn[,"gene_id"] == gene) ,"dknn_values"] <- gene_dknn[k]
      }
      
      
      #################### DKNN simulation
      # Extract the distance values from the distance matrix
      dist_values <- dist_matrix[!is.na(dist_matrix)]
      nb_of_gene_sim <- nrow(select_for_correlation) * 1 #REMPLACER PAR NB DE GENE x2 POUR RAJOUTER UNE SIMULATION
      sim_dknn <- vector()
      
      # Generate simulated distances
      for (sim_nb in 1:nb_of_gene_sim) {
        
        # Randomly sample distances for one simulated gene
        nb_gene <- nrow(select_for_correlation)
        dist_sim <- sample(dist_values, size=nb_gene, replace=FALSE)
        
        # Extract the k nearest neighbors of these simulated gene
        dist_sim <- dist_sim[order(dist_sim)]
        sim_dknn[sim_nb] <- dist_sim[k]
        
      }
      
      
      #################### Determine the DKNN threshold (or critical distance)
      # Order genes by DKNN values
      df_dknn <- df_dknn[order(df_dknn[,"dknn_values"]),]
      df_dknn[, c("nb_dknn_sim", "nb_dknn_obs", "ratio_sim_obs")] <- 0
      
      
      #Recuperer les valeurs sim_dknn < obs_dknn(i)
      for(i in 2:nb_gene) {
        
        #Compute the number of simulated DKNN values under DKNN value at rank i
        #If there is no simulated DKNN values under DKNN value at rank i, put 0 in nb_dknn_sim
        if (min(sim_dknn) < df_dknn[i,"dknn_values"]) {
          nb_dknn_sim <- sim_dknn[which(sim_dknn < df_dknn[i,"dknn_values"])]
          nb_dknn_sim <- length(nb_dknn_sim)
          df_dknn[i,"nb_dknn_sim"] <- nb_dknn_sim
          
        }else {
          nb_dknn_sim <- 0 
          df_dknn[i,"nb_dknn_sim"] <- nb_dknn_sim 
          
        }
        #Compute the number of observed DKNN values under DKNN value at rank i
        nb_dknn_obs <- length(df_dknn[which(df_dknn[,"dknn_values"] < df_dknn[i,"dknn_values"]), "dknn_values"]) #Nombre de valeurs de dknn observees inferieures a la valeur dknn_obs(i)
        df_dknn[i,"nb_dknn_obs"] <- nb_dknn_obs
        
        #Compute the ratio between number of simulated DKNN values and the number of observed DKNN values under DKNN value at rank i
        df_dknn[i,"ratio_sim_obs"] <- nb_dknn_sim/nb_dknn_obs
      }
      
      #################### Select genes with a distance value under critical distance
      critical_distance <- df_dknn[which(df_dknn[,"ratio_sim_obs"] > fdr*0.01)[1], "dknn_values"]
      selected_genes <- df_dknn[df_dknn[,"dknn_values"] < critical_distance,]
      # Remove genes not selected on the previously created list including observed distance values
      l_knn_selected <- l_knn[as.character(selected_genes[,"gene_id"])]
      
      
      
      ####################  Create the input file for mcl algorithm
      # Distance matrix of selected genes
      input_mcl <- (input_mcl+1)/2
      input_mcl <- input_mcl[as.character(selected_genes[,"gene_id"]), as.character(selected_genes[,"gene_id"])]
      
      # Remove distances values not conserved
      for (gene in as.character(selected_genes[,"gene_id"])) {
        gene_comb_loop <- names(l_knn_selected[[gene]])
        gene_comb_loop <- gene_comb_loop[which(gene_comb_loop %in% selected_genes[,"gene_id"])]
        
        input_mcl[-which(rownames(input_mcl) %in% gene_comb_loop), gene] <- 0
      }
      
      # Melt the matrix to adapt it the the MCL format
      input_mcl <- melt(input_mcl)
      
      # Remove NA values and zero values
      input_mcl <- input_mcl[-which(is.na(input_mcl[,"value"]) | input_mcl[,"value"] == 0),]
      
      write.table(input_mcl, 
                  file=file.path( output_path, paste0(name, ".dbf_out.txt")), 
                  row.names=FALSE, 
                  col.names=FALSE, 
                  sep='\t', 
                  quote=FALSE)
      
      
      
      #################### Create the ClusterSet object
      obj <- new("ClusterSet")
      obj@algorithm <- "DBFMCL"
      
      if (length(selected_genes[,"gene_id"]) > 0) {
        obj@data <- as.matrix(data[selected_genes[,"gene_id"],])
        obj@gene_patterns <- rep(1, nrow(obj@data))
        obj@size <- nrow(obj@data)
        obj@center <- matrix(
          apply(obj@data[obj@gene_patterns == 1, ],
                2,
                mean,
                na.rm = TRUE
          ),
          nrow = 1
        )
        obs_dknn <- as.vector(df_dknn[,"dknn_values"])
        names(obs_dknn) <- df_dknn[,"gene_id"]
        obj@distances <- obs_dknn
        obj@simulated_distances <- sim_dknn
        obj@critical_distance <- critical_distance
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


#' @title
#' Invokes the command line version of Markov CLustering algorithm (MCL).
#' @description
#'  This function invokes the mcl system command. MCL is a clustering algorithm
#' for graphs that was developped by Stijn van Dongen (see references for
#' further informations).
#' @param name a character string corresponding to the file name.
#' @param inflation the main control of MCL. Inflation affects cluster
#' granularity. It is usually chosen somewhere in the range \code{[1.2-5.0]}.
#' \code{inflation = 5.0} will tend to result in fine-grained clusterings, and
#' whereas \code{inflation = 1.2} will tend to result in very coarse grained
#' clusterings. By default, \code{inflation = 2.0}. Default setting gives very
#' good results for microarray data when k is set around 100.
#' @param input_path a character string representing the directory path of 
#' the input file used by mcl. Default is the current working directory.
#' @param silent if set to TRUE, the progression of the MCL partitionning is
#' not displayed.
#' @param threads The number of threads to use.
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
#' @export mcl_system_cmd
mcl_system_cmd <- function(name, inflation = 2.0, input_path = ".", silent = FALSE, threads=1) {
  ## testing the system
  if (.Platform$OS.type != "windows") {

    ## Testing mcl installation
    if (system("mcl --version | grep 'Stijn van Dongen'", intern = TRUE) > 0) {
      if (!silent) {
        cat("Running mcl (graph partitioning)... \n")
        verb <- ""
      }
      else {
        verb <- "-V all "
      }
      if (inflation != 2) {
        i <- paste("-I ", as.character(round(inflation, 1)), sep = "")
      } else {
        i <- "-I 2.0"
      }
      threads <- paste("-te", threads, sep = " ")
      ## launching mcl program
      cmd <- paste0("mcl ",
                   input_path, "/", name, ".dbf_out.txt ",
                   i,
                   " --abc -o ",
                   input_path, "/", name, ".mcl_out.txt ",
                   verb,
                   threads)
      cmd <- gsub(pattern = "//", replacement = "/", x = cmd)
      system(cmd)

      if (!silent) {
        print_msg("Done", msg_type="INFO")
        print_msg(paste0("creating file : ",
                        file.path(getwd(), paste(name, ".mcl_out.txt", sep = ""))),
                  msg_type="INFO")
      }
    } else {
      stop(
        "\t--> Please install mcl on your computer...\n",
        "\t--> You can download it from : 'http://www.micans.org/mcl/'\n\n"
      )
    }
  }
  else {
    stop("--> A unix-like OS is required to launch mcl and cluster programs.")
  }
}

#########################################################
##      END PACKAGE scigenex
#########################################################
