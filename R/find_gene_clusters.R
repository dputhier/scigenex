################################################################
##        Main script of the R PACKAGE : scigenex
##
## Authors : J. BAVAIS, BERGON A,
##  with the collaboration of LOPEZ F., TEXTORIS J. and PUTHIER D.
##
##
################################################################


#################################################################
##    find_gene_clusters function
#################################################################

#' @title
#' The central function of Scigenex package used to select informative genes
#' and create partitions of co-expressed genes.
#' @description
#' The find_gene_clusters function used an implementation of the  DBFMCL algorithm.
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
#' falling in a dense area, DBF-MCL computes simulated DKNN values based on 
#' distances corresponding to genes with the genes highest dknn values. Computed 
#' distributions of simulated DKNN are used to compute a FDR value for each observed 
#' DKNN value. Genes with corresponding FDR value <= user defined FDR are selected 
#' and used to construct a graph. In this graph, edges are constructed between 
#' two genes (nodes) if one of them belongs to the k-nearest neighbors of the other. 
#' Edges are weighted based on the respective coefficient of correlation (\emph{i.e.}, 
#' similarity) and the graph obtained is partitioned using the Markov CLustering 
#' Algorithm (MCL).
#'
#' @param data a \code{matrix}, \code{data.frame} or \code{Seurat} object.
#' @param name a prefix for the names of the intermediary files created by DBF
#' and MCL.
#' @param output_path a character string representing the data directory where
#' output files will be stored. Default to current working directory.
#' @param mcl_threads An integer to determine number of threads for mcl algorithm.
#' @param distance_method A method to compute the distance to the k-th nearest neighbor 
#' (one of "pearson", "cosine" or "euclidean").
#' @param min_nb_supporting_cell Use to filter out the clusters. The minimum number 
#' of cell supporting a cluster. A cell supports a cluster if it expresses at least 
#' min_pct_gene_expressed \% of the genes from the cluster.
#' @param min_pct_gene_expressed See min_nb_supporting_cell argument.
#' @param min_cluster_size Minimum number of element inside a cluster. 
#' MCL tend to create lots of clusters with very few (e.g 2) objects.
#' @param highest During the process, genes will be ordered by their distance to their 
#' k nearest neighbors (dknn). This parameter controls the fraction of genes with high dknn (ie. noise)
#' whose neighborhood (i.e associated distances) will be used to compute simulated values. 0 would mean to 
#' use all the genes.  
#' @param k the neighborhood size.
#' @param row_sum Select only lines whose row sum is greater than row_sum (may be 
#' -Inf if no filter needed).
#' @param fdr an integer value corresponding to the false discovery rate (range: 0 to 100).
#' @param inflation the main control of MCL. Inflation affects cluster
#' granularity. It is usually chosen somewhere in the range \code{[1.2-3]}.
#' \code{inflation = 5.0} will tend to result in fine-grained clusterings, and
#' whereas \code{inflation = 1.2} will tend to result in very coarse grained
#' clusterings. By default, \code{inflation = 2.0}. Default setting gives very
#' good results for microarray data when k is set between 70 and 250.
#' @param no_dknn_filter Do not perform distance to k nearest neighbor (DBF) step. This may produce
#' a very large graph.
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
#' @author Bergon A., Bavais J., Textoris J., Lopez F and Puthier D.
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
#' set_verbosity(2)
#' m <- create_4_rnd_clust()
#' 
#' res <- find_gene_clusters(data=m,
#'                              distance_method="pearson",
#'                              inflation = 2,
#'                              k=75,
#'                              row_sum=-Inf,
#'                              highest=0.3,
#'                              min_nb_supporting_cell = 0,
#'                              fdr = 1e-8)
#' plot_heatmap(res, row_labels = F)
#'
#' @export find_gene_clusters
find_gene_clusters <- function(data = NULL,
                               output_path = tempdir(),
                               mcl_threads = 1,
                               name = NULL,
                               distance_method = c("pearson",
                                                   "cosine",
                                                   "euclidean"),
                               min_nb_supporting_cell = 2,
                               min_pct_gene_expressed = 40,
                               min_cluster_size = 10,
                               highest = 0.25,
                               k = 50,
                               row_sum = 0,
                               fdr = 0.05,
                               inflation = 2,
                               no_dknn_filter=FALSE,
                               seed = 123) {
  ## testing the system
  if (.Platform$OS.type == "windows") {
    stop("\t--> A unix-like OS is required to launch mcl and cluster programs.")
  }
  
  ## getting parameters
  data_matrix <- get_data_for_scigenex(data = data)
  
  if (highest < 0 || highest > 1) {
    stop("highest argument should be >= 0 and <= 1.")
  }
  
  if (min_pct_gene_expressed < 0 || min_pct_gene_expressed > 100) {
    stop("min_pct_gene_expressed argument should be >= 0 and <= 100.")
  }
  
  if (is.null(name))
    name <- create_rand_str()
  
  # Put the current working directory in output_path or path
  if (output_path == ".") {
    output_path <- getwd()
  }
  # Check if output directory exists. If not stop the command.
  if (!file.exists(output_path)) {
    stop("Output directory provided does not exist.")
  }
  
  distance_method <- match.arg(distance_method)
  
  print_msg(paste0(
    "The following parameters will be used :",
    "\n\tOuput directory: ", output_path,
    "\n\tName: ", name,
    "\n\tDistance method: ", distance_method,
    "\n\tMinimum number of cell supporting each cluster: ", min_nb_supporting_cell,
    "\n\tThe min fraction of genes expressed by supporting cells : ", min_pct_gene_expressed,
    "\n\tRow sum threshold: ", row_sum,
    "\n\tMinimum cluster size: ", min_cluster_size,
    "\n\tNumber of neighbors: ", k,
    "\n\tFDR: ", fdr, "%",
    "\n\tInflation: ", inflation,
    "\n\tVerbosity Level: ", get_verbosity(),
    "\n\n"),
    msg_type = "INFO"
  )
  

    ## DBF algorithm, returns a ClusterSet object
    obj <- DBF(
      data_matrix,
      output_path = output_path,
      name = name,
      distance_method = distance_method,
      highest = highest,
      k = k,
      row_sum=row_sum,
      fdr = fdr,
      no_dknn_filter=no_dknn_filter,
      seed = seed
    )

  dbf_out_file <- paste0(output_path, "/", name, ".dbf_out.txt")
  dbf_out_file <-
    gsub(pattern = "//",
         replacement = "/",
         x = dbf_out_file)
  mcl_out_file <- paste0(output_path, "/", name, ".mcl_out.txt")
  mcl_out_file <-
    gsub(pattern = "//",
         replacement = "/",
         x = mcl_out_file)
  
  print_msg("DBF completed. Starting MCL step.", msg_type = "INFO")
  
  if (length(readLines(dbf_out_file)) > 0) {
    ## Launching mcl (command line)
    mcl_system_cmd(
      name,
      inflation = inflation,
      input_path = output_path,
      threads = mcl_threads
    )
    
    
    print_msg(paste0("Reading and filtering MCL output: ", mcl_out_file),
              msg_type = "INFO")
    
    ## getting mcl results into the ClusterSet object
    mcl_cluster <- readLines(mcl_out_file)
    gene_list <- NULL
    clusters <- NULL
    size <- NULL
    nb <- 0
    nb_cluster_deleted_size <- 0
    nb_cluster_deleted_pct_gene <- 0
    
    for (i in seq_along(mcl_cluster)) {
      gene_in_clust <- unlist(strsplit(mcl_cluster[i], "\t"))
      
      # Cluster size filtering
      if (length(gene_in_clust) > min_cluster_size) {
        # Filter on the number of cells
        # expressing at least a certain purcentage of genes.
        if (min_nb_supporting_cell > 0) {
          cur_clust <- data_matrix[gene_in_clust, ]
          cur_clust[cur_clust > 0] <- 1
          min_pct_gene_expressed_cur <-
            apply(cur_clust, 2, sum) / nrow(cur_clust) * 100
          min_nb_supporting_cell_cur <-
            length(min_pct_gene_expressed_cur[min_pct_gene_expressed_cur >= min_pct_gene_expressed])
          
          if (min_nb_supporting_cell_cur >= min_nb_supporting_cell) {
            nb <- nb + 1
            gene_list <- c(gene_list, gene_in_clust)
            clusters <- c(clusters, rep(nb, length(gene_in_clust)))
            size <- c(size, length(gene_in_clust))
          } else {
            nb_cluster_deleted_pct_gene <- nb_cluster_deleted_pct_gene + 1
          }
        } else {
          nb <- nb + 1
          gene_list <- c(gene_list, gene_in_clust)
          clusters <- c(clusters, rep(nb, length(gene_in_clust)))
          
          if (is.null(size)) {
            size <- length(gene_in_clust)
          } else {
            size <- c(size, length(gene_in_clust))
          }
        }
      } else {
        nb_cluster_deleted_size <- nb_cluster_deleted_size + 1
      }
    }
    
    print_msg(
      paste(
        nb_cluster_deleted_pct_gene,
        " clusters filtered out after MCL partitionning (min_nb_supporting_cell)."
      ),
      msg_type = "INFO"
    )
    print_msg(
      paste(
        nb_cluster_deleted_size,
        " clusters filtered out after MCL partitionning (cluster size)."
      ),
      msg_type = "INFO"
    )
    print_msg(paste(length(clusters),
                    " genes selected after analysis."),
              msg_type = "INFO")
    print_msg(paste(nb, " clusters conserved after MCL partitioning."),
              msg_type = "INFO")
    
    ## build ClusterSet object
    if (nb > 0) {
      obj@data <- as.matrix(data_matrix[gene_list, ])
      clusters_df <- data.frame("gene_clusters" = clusters,
                                "gene_id" = rownames(obj@data))
      obj@gene_clusters <- split(clusters_df[,"gene_id"], f=clusters_df$gene_clusters)
      obj@gene_clusters_metadata <- list("cluster_id" = as.numeric(names(obj@gene_clusters)),
                                         "number" = max(names(obj@gene_clusters)),
                                         "size" = unlist(lapply(obj@gene_clusters, length)))
      
      centers <- matrix(ncol = ncol(data_matrix), nrow = nb)
      ## calcul of the mean profils
      for (i in 1:nb) {
        centers[i, ] <- apply(obj@data[obj@gene_clusters[[i]], ],
                              2, mean,
                              na.rm = TRUE)
      }
      obj@dbf_output$center <- centers
      
      obj@cells_metadata <- data.frame("cells_barcode" = colnames(obj@data))
      
      ## add DBFMCL parameters used to build this object
      obj@parameters <- list(
        filename = name,
        distance_method = distance_method,
        k = k,
        inflation = inflation,
        highest = highest,
        fdr = fdr,
        min_nb_supporting_cell = min_nb_supporting_cell,
        min_pct_gene_expressed, min_pct_gene_expressed,
        min_cluster_size = min_cluster_size,
        row_sum = row_sum,
        seed = seed
      )
    }
  } else {
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
#' This function is an internal function used by \code{\link{find_gene_clusters}} to detect
#' informative elements (\emph{i.e.}, those that belong to dense regions). User
#' should not use this function. Instead they can use the \code{\link{find_gene_clusters}}
#' function with \code{clustering} argument set to \code{FALSE}.
#'
#' See \code{\link{find_gene_clusters}}
#'
#' @param data a \code{matrix}, \code{data.frame} or \code{Seurat} object.
#' @param output_path a character string representing the data directory where
#' output files will be stored. Default to current working directory.
#' @param name a prefix for the names of the intermediary files created by DBF
#' and MCL.
#' @param distance_method A method to compute the distance to the k-th nearest neighbor 
#' (one of "pearson", "cosine" or "euclidean").
#' @param silent if set to TRUE (default), the progression of distance matrix
#' calculation is not displayed.
#' @param highest During the process, genes will be ordered by their distance to their 
#' k nearest neighbors (dknn). This parameter controls the fraction of genes with high dknn (ie. noise)
#' whose neighborhood (i.e associated distances) will be used to compute simulated values. 0 would mean to 
#' use all the genes.  
#' @param k the neighborhood size.
#' @param row_sum Select only lines whose row sum is greater than row_sum (may be 
#' -Inf if no filter needed).
#' @param fdr a value for the false discovery rate.
#' @param no_dknn_filter Do not perform distance to k nearest neighbor (DBF) step.
#' @param seed specify seeds for random number generator.
#' @section Warnings: Works only on UNIX-alikes platforms.
#' @author Bergon A., Bavais J., Textoris J., Lopez F and Puthier D.
#' @seealso \code{\link{find_gene_clusters}}
#' @references Lopez F.,Textoris J., Bergon A., Didier G., Remy E., Granjeaud
#' S., Imbert J. , Nguyen C. and Puthier D. TranscriptomeBrowser: a powerful
#' and flexible toolbox to explore productively the transcriptional landscape
#' of the Gene Expression Omnibus database. PLoSONE, 2008;3(12):e4001.
#' @keywords manip
#' @export DBF
DBF <- function(data,
                output_path = ".",
                name = NULL,
                distance_method = c("pearson", "cosine", "euclidean"),
                silent = FALSE,
                highest = 0.5,
                k = 100,
                row_sum = 0,
                fdr = 0.05,
                no_dknn_filter=FALSE,
                seed = 123) {
  ## testing the system
  if (.Platform$OS.type == "windows") {
    stop("\t--> A unix-like OS is required to launch mcl program.")
  }
  
  if (is.null(data) ||
      !inherits(data, c("matrix", "data.frame", "Seurat"))) {
    stop("\t--> Please provide a matrix...\n\n")
  }
  
  ## set a seed for reproductibility
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  if (highest < 0 || highest > 1) {
    stop("highest argument should be >= 0 and <= 1.")
  }
  
  ## getting data and parameters
  if (output_path == ".")
    output_path <- getwd()
  if (is.null(name))
    name <- "exprs"
  
  data <- get_data_for_scigenex(data = data)
  
  distance_method <- match.arg(distance_method)
  
  if (silent) {
    print_msg(
      paste0(
        "Computing distances to the kth-nearest neighbors ",
        "and associated FDR values... \n"
      ),
      msg_type = "INFO"
    )
  }
  
  # Directory and name of the principal output
  outfile <- paste(output_path, "/", name, ".dbf_out.txt", sep = "")
  outfile <- gsub(pattern = "//",
                  replacement = "/",
                  x = outfile)
  
  #################### Correlation and distance matrices
  # Remove genes with 0 values for all cells
  
  genes_to_keep <- rowSums(data) > row_sum
  select_for_correlation <- data[genes_to_keep, ]
  
  # Compute gene-gene correlation/distance matrix
  
  print_msg(
    paste0(
      "Computing distances using selected method: ",
      distance_method
    ),
    msg_type = "INFO"
  )
  
  # Compute correlations and corresponding 
  # distance matrix. Note that for pearson
  # and cosine, 0 < distance < 2.
  if (distance_method == "pearson") {
    dist_matrix <- corSparse(t(select_for_correlation))
    dist_matrix <- 1 - dist_matrix
  } else if (distance_method == "cosine") {
    dist_matrix <- as.matrix(corSparse(t(select_for_correlation)))
    dist_matrix <- 1 - dist_matrix
  } else if (distance_method == "euclidean") {
    dist_matrix <- as.matrix(dist(select_for_correlation))
  }
  
  print_stat("Distance matrix stats", 
             data = dist_matrix, msg_type = "DEBUG")
  
  # Compute min and max distances.
  # They will be used later to compute weights.
  min_dist <- min(dist_matrix)
  max_dist <- max(dist_matrix)
  
  # Set the rownames / colnames of the distance matrix
  rownames(dist_matrix) <- rownames(select_for_correlation)
  colnames(dist_matrix) <- rownames(select_for_correlation)
  
  # The distance from a gene to itself is 'hidden'
  diag(dist_matrix) <- NA
  

  # Create a dataframe to store the DKNN values.
  # Gene_id appear both as rownames and column
  # for coding convenience
  df_dknn <- data.frame(
    dknn_values = rep(NA, nrow(dist_matrix)),
    row.names = rownames(dist_matrix),
    gene_id = rownames(dist_matrix)
  )
  
  # This list will be used to store the
  # dknn values
  l_knn <- list()
  
  #################### DKNN for each genes
  # Extract the DKNN for each gene
  print_msg(paste0("Computing distances to KNN."), msg_type = "INFO")
  
  for (pos in seq_len(nrow(dist_matrix))) {
    gene <- rownames(df_dknn)[pos]
    gene_dist <- dist_matrix[gene, ]
    # Reorder the distance values in increasing order.
    # The distance from a gene to itself (previously set to NA)
    # is placed at the end of the vector.
    gene_dist <- sort(gene_dist, na.last = T)[1:k]
    
    # Add the neigbhors to the list
    l_knn[[gene]] <- gene_dist
    
    # Select the kth pearson correlation values. 
    # This value corresponds to the DKNN of the gene(i)
    df_dknn[gene, "dknn_values"] <- gene_dist[k]
  }
    
  print_stat("Observed DKNN stats", 
             data = df_dknn$dknn_values, msg_type = "DEBUG")

  sim_dknn <- vector()
  critical_distance <- vector()
  
  if(!no_dknn_filter){
    
      #################### DKNN simulation
      print_msg(paste0("Computing simulated distances to KNN."), msg_type = "INFO")
      
      nb_selected_genes <- nrow(select_for_correlation)
      
      if (highest == 0) {
        tresh_dknn <- min(df_dknn$dknn_values)
      } else if (highest == 1) {
        tresh_dknn <- max(df_dknn$dknn_values)
      } else {
        tresh_dknn <- stats::quantile(df_dknn$dknn_values, highest)
      }
      
      gene_with_low_dknn <- df_dknn$gene_id[df_dknn$dknn_values > tresh_dknn]
      dist_values_sub <- dist_matrix[gene_with_low_dknn, ]
      dist_values_sub <- dist_values_sub[!is.na(dist_values_sub)]
      
      
      for (sim_nb in 1:nb_selected_genes) {
        # Randomly sample distances for one simulated gene
        dist_sim <- sample(dist_values_sub,
                           size = nb_selected_genes,
                           replace = FALSE)
        
        # Extract the k nearest neighbors of these simulated gene
        dist_sim <- sort(dist_sim)
        sim_dknn[sim_nb] <- dist_sim[k]
      }
    
      print_stat("Simulated DKNN stats", 
                  data = sim_dknn, msg_type = "DEBUG")
    
      # The simulated DKNN values follow a normal distribution.
      # Compute the parameters of this distribution
      mean_sim <- mean(sim_dknn)
      sd_sim <- sd(sim_dknn)
    
      #################### Determine the DKNN threshold (or critical distance)
      # Order genes by DKNN values
      print_msg(paste0("Computing distances to KNN (DKNN) threshold."),
                msg_type = "INFO")
      
      # Order the dknn values from low to high
      df_dknn <- df_dknn[order(df_dknn$dknn_values), ]
      df_dknn[, "FDR"] <- NA
    
      # Compute the FDR
      for (i in 1:nb_selected_genes) {
        gene <- df_dknn$gene_id[i]
        df_dknn[gene, "FDR"] <-
          stats::pnorm(df_dknn[gene, "dknn_values"],
                       mean = mean_sim,
                       sd = sd_sim,
                       lower.tail = T) / (i / nb_selected_genes) * 100
      }
    
      df_dknn$FDR[df_dknn$FDR > 100] <- 100
    
      #################### Select genes with a distance value under critical distance
      print_msg(paste0("Selecting informative genes."), msg_type = "INFO")
      fdr_tresh_pos <- which(df_dknn$FDR > fdr)[1]
      selected_genes <- df_dknn[1:fdr_tresh_pos,]$gene_id
      critical_distance <-df_dknn$dknn_values[fdr_tresh_pos]

  }else{
    selected_genes <- rownames(dist_matrix)
  }
  
  # Remove genes not selected on the previously created list including observed distance values
  
  l_knn_selected <- l_knn[selected_genes]
  
  ####################  Create the input file for mcl algorithm
  print_msg(paste0("Creating the input file for MCL algorithm."), msg_type = "INFO")
  
  mcl_out_as_list_of_df <- list()
  
  for (g in names(l_knn_selected)) {
    mcl_out_as_list_of_df[[g]] <- data.frame(
      src = g,
      dest = names(l_knn_selected[[g]]),
      weight = l_knn_selected[[g]]
    )
  }
  
  mcl_out_as_df <- do.call(rbind, mcl_out_as_list_of_df)
  
  #############  Convert distances into weights
  # scale dist between 0..1
  mcl_out_as_df$weight <- 
    (mcl_out_as_df$weight - min_dist) / (max_dist - min_dist)
  # Convert scaled dist to weight
  mcl_out_as_df$weight <- abs(mcl_out_as_df$weight - 1)

  print_stat("Graph weights (after convertion)", 
             data = mcl_out_as_df$weight, 
             msg_type = "DEBUG")
  
  # A and B are added if A is in the
  # neighborhood of B and B in the neighborhood
  # of A
  mcl_out_as_df <-
    mcl_out_as_df[duplicated(t(apply(mcl_out_as_df[, c("src", "dest")], 1, sort))), ]
  
  # Ensure an edge (A->B and B->A)
  # is not defined twice
  mcl_out_as_df <-
    mcl_out_as_df[!duplicated(t(apply(mcl_out_as_df[, c("src", "dest")], 1, sort))), ]
  
  print_msg(paste0("Writing table."), msg_type = "INFO")
  data.table::fwrite(
    mcl_out_as_df,
    file = file.path(output_path,
                     paste0(name, ".dbf_out.txt")),
    sep = "\t",
    eol = "\n",
    row.names = F,
    col.names = FALSE
  )
  
  #################### Create the ClusterSet object
  obj <- new("ClusterSet")
  obj@parameters <- list(algorithm = "DBFMCL",
                         name = name,
                         distance_method = distance_method,
                         k = k,
                         highest = highest,
                         fdr = fdr,
                         row_sum = row_sum,
                         seed = seed
                         )
  
  if (length(selected_genes) > 0) {
    obj@data <- as.matrix(data[selected_genes, ])
    obj@gene_clusters <- list("1" = rownames(obj@data))
    obj@gene_clusters_metadata <- list("cluster_id" = as.numeric(names(obj@gene_clusters)),
                                       "number" = max(names(obj@gene_clusters)),
                                       "size" = nrow(obj@data))
    
    obs_dknn <- as.vector(df_dknn[, "dknn_values"])
    names(obs_dknn) <- df_dknn[, "gene_id"]
    center = matrix(apply(obj@data[obj@gene_clusters$`1`, ],
                          2,
                          mean,
                          na.rm = TRUE),
                    nrow = 1)
    obj@dbf_output <- list("dknn" = obs_dknn,
                           "simulated_dknn" = sim_dknn,
                           "critical_distance" = critical_distance,
                           "fdr" = df_dknn$FDR,
                           "center" = center)
    # obj@normal_model_mean <- mean_sim
    # obj@normal_model_sd <- mean_sd
  }
  
  return(obj)
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
mcl_system_cmd <- function(name,
                           inflation = 2.0,
                           input_path = ".",
                           silent = FALSE,
                           threads = 1) {
  ## testing the system
  if (.Platform$OS.type == "windows") {
    stop("--> A unix-like OS is required to launch the MCL program.")
  }else {
    print_msg("Running mcl_system_cmd() under a unix-like system.", msg_type = "DEBUG")
  }
  
  ## Testing mcl installation
  if (system("mcl --version | grep 'Stijn van Dongen'", intern = TRUE) > 0) {
    print_msg("Found MCL program in the path...", msg_type = "DEBUG")
  } else {
    stop(
      "\t--> Please install mcl on your computer...\n",
      "\t--> You can download it from : 'http://www.micans.org/mcl/'\n\n"
    )
  }
  
  
  print_msg("Running MCL (graph partitioning)...", msg_type = "DEBUG")
  
  if (get_verbosity() > 0) {
    verb <- ""
  }
  else {
    verb <- "-V all "
  }
  
  i <- paste("-I ", as.character(round(inflation, 1)), sep = "")
  
  threads <- paste("-te", threads, sep = " ")
  
  ## launching mcl program
  cmd <- paste0("mcl ",
                input_path,
                "/",
                name,
                ".dbf_out.txt ",
                i,
                " --abc -o ",
                input_path,
                "/",
                name,
                ".mcl_out.txt ",
                verb,
                threads)
  
  cmd <- gsub(pattern = "//",
              replacement = "/",
              x = cmd)
  
  print_msg(paste0("Running command: ", cmd), msg_type = "DEBUG")
  
  system(cmd)
  
  print_msg("MCL step is finished.", msg_type = "DEBUG")
  print_msg(paste0("creating file : ",
                   file.path(
                     getwd(), paste(name, ".mcl_out.txt", sep = "")
                   )),
            msg_type = "DEBUG")
  
}

#########################################################
##      END PACKAGE scigenex
#########################################################
