################################################################################
#' Gene clustering using the Markov Cluster Algorithm (MCL) method.
#'
#' This function performs gene clustering using the MCL algorithm. 
#' The method starts by creating a graph with genes as nodes and edges connecting each gene to its nearest neighbors.
#' Then the method use the MCL algorithm to detect clusters of co-expressed genes (method argument).
#'
#' @param object A ClusterSet object.
#' @param s If method="closest_neighborhood", s is an integer value indicating the size of the neighbourhood used for graph construction. Default is 5.
#' @param inflation A numeric value indicating the MCL inflation parameter. Default is 2.
#' @param method Which method to use to build the graph. If "closest_neighborhood", creates an edge between 
#' two selected genes a and b if b is part of the kg closest nearest neighbors of a (with kg < k). If "reciprocal_neighborhood"),
#' inspect the neighborhood of size k of all selected genes and put an edge between two genes a and b if they are reciprocally 
#' in the neighborhood of the other
#' @param threads An integer value indicating the number of threads to use for MCL.
#' @param output_path a character indicating the path where the output files will be stored.
#' @param name a character string giving the name for the output files. If NULL, a random name is generated.
#' @param keep_nn Deprecated. Use 'method' instead.
#' @return A ClusterSet object
#' @references
#' - Van Dongen S. (2000) A cluster algorithm for graphs. National
#' Research Institute for Mathematics and Computer Science in the 1386-3681.
#' @examples
#' # Restrict vebosity to info messages only.
#' library(Seurat)
#' set_verbosity(1)
#' 
#' # Load a dataset
#' load_example_dataset("7871581/files/pbmc3k_medium")
#' 
#' # Select informative genes
#' res <- select_genes(pbmc3k_medium,
#'                     distance = "pearson",
#'                     row_sum=5)
#'                     
#' # Cluster informative features
#' 
#' ## Method 1 - Construct a graph with a 
#' ## novel neighborhood size
#' res <- gene_clustering(res, method="closest_neighborhood",
#'                        inflation = 1.5, threads = 4)
#'                        
#' # Display the heatmap of gene clusters
#' res <- top_genes(res)
#' plot_heatmap(res)
#' plot_heatmap(res, cell_clusters = Seurat::Idents(pbmc3k_medium))
#' 
#' 
#' ## Method 2 - Conserve the same neighborhood 
#' ## size
#' res <- gene_clustering(res, 
#'                        inflation = 2.2,
#'                        method="reciprocal_neighborhood")
#'                        
#' # Display the heatmap of gene clusters
#' res <- top_genes(res)
#' plot_heatmap(res)
#' plot_heatmap(res, cell_clusters = Seurat::Idents(pbmc3k_medium))
#' 
#' @export gene_clustering

gene_clustering <- function(object = NULL,
                            s = 5,
                            inflation = 2,
                            method=c("closest_neighborhood", "reciprocal_neighborhood"),
                            threads = 1,
                            output_path = tempdir(),
                            name = NULL,
                            keep_nn = FALSE) {
  
  
  method <- match.arg(method)
  
  if(keep_nn){
    print_msg("The use of keep_nn=TRUE is deprecated. Use 'method' argument set to 'reciprocal_neighborhood'.")
    method <- "reciprocal_neighborhood"
  }
    
  # Construct graph for mcl and save it in a new file
  if (method == "reciprocal_neighborhood") {
    object <- keep_dbf_graph(object = object,
                             output_path = output_path,
                             name = name)
  } else {
    object <- construct_new_graph(object = object,
                                  k = s,
                                  output_path = output_path,
                                  name = name)
  }
  
  # Run MCL algorithm
  object <- mcl_system_cmd(object = object,
                           inflation = inflation,
                           threads = threads)
  
  print_msg(paste0("Adding results to a ClusterSet object."), msg_type = "INFO")
  
  # Update ClusterSet object
  ## Read mcl out file
  mcl_out_file <- file.path(object@parameters$output_path,
                            paste0(object@parameters$name, ".mcl_out.txt"))
  
  ## Extract gene clusters
  print_msg(paste0("Read MCL output file."), msg_type = "DEBUG")
  mcl_cluster <- readLines(mcl_out_file)
  mcl_cluster <- strsplit(mcl_cluster, "\t")
  names(mcl_cluster) <- seq(1, length(mcl_cluster))
  
  ## Add gene clusters to ClusterSet object
  object@gene_clusters <- mcl_cluster
  
  ## Update gene_cluster_metadata slots
  object@gene_clusters_metadata <- list("cluster_id" = as.numeric(names(object@gene_clusters)),
                                        "number" = max(as.numeric(names(object@gene_clusters))),
                                        "size" = unlist(lapply(object@gene_clusters, length)))
  
  ## Update data slot
  object@data <- object@data[unlist(object@gene_clusters, use.names = FALSE), ]
  
  ## Compute centers
  print_msg(paste0("Compute centers."), msg_type = "DEBUG")
  nb_clusters = length(names(object@gene_clusters))
  centers <- matrix(ncol = ncol(object@data),
                    nrow = nb_clusters)
  colnames(centers) <- colnames(object@data)
  rownames(centers) <- names(object@gene_clusters)
  
  ## calcul of the average profiles
  for (i in 1:nb_clusters) {

      centers[i, ] <- apply(object@data[object@gene_clusters[[i]], , drop=FALSE],
                            2, mean,
                            na.rm = TRUE)
    
  }
  
  object@dbf_output$center <- centers
  rownames(object@dbf_output$center) <- names(object@gene_clusters)
  
  object@cells_metadata <- data.frame("cells_barcode" = colnames(object@data),
                                      row.names = colnames(object@data))
  
  return(object)
}



################################################################################
#' Construct a new graph for a ClusterSet object
#'
#' This function constructs a new graph for a given ClusterSet object based on selected genes. 
#' The graph is constructed using the Density K-Nearest Neighbor (DKNN) method.
#'
#' @param object A ClusterSet object.
#' @param k The number of nearest neighbors to consider for each gene.
#' @param output_path a character indicating the path where the output files will be stored.
#' @param name a character string giving the name for the output files. If NULL, a random name is generated.
#'
#' @return A ClusterSet object.
#'
#' @examples
#' # Set verbosity to 1 to display info messages only.
#' set_verbosity(1)
#' 
#' # Create a matrix with 4 signatures
#' m <- create_4_rnd_clust()
#' 
#' # Select informative features
#' res <- select_genes(m[1:600,],
#'                     distance = "pearson",
#'                     k = 75,
#'                     noise_level = 0.3,
#'                     fdr = 1e-8,
#'                     row_sum = -Inf)
#'                     
#' # Construct a new graph based on genes selected with select_genes()
#' res <- construct_new_graph(object = res, k = 5)
#' @keywords internal
#' @export construct_new_graph
construct_new_graph <- function(object = NULL,
                                k = 5,
                                output_path = tempdir(),
                                name = NULL) {
  ## Check format object arg
  check_format_cluster_set(object)
  
  ## Generate output path
  # Get working directory if output_path is "."
  if (output_path == ".") {
    output_path <- getwd()
  }
  # Check if output directory exists. If not stop the command.
  if (!file.exists(output_path)) {
    stop("Output directory provided does not exist.")
  }
  
  # Create a random string if name is not provided
  if (is.null(name)){
    name <- create_rand_str()
  }
  
  # Directory and name of the principal output
  path_input_mcl <- paste(output_path, "/", name, ".input_mcl.txt", sep = "")
  path_input_mcl <- gsub(pattern = "//",
                         replacement = "/",
                         x = path_input_mcl)
  
  
  
  
  # Extract selected genes
  selected_genes <- object@gene_clusters$`1`
  # Extract normalized gene expression matrix
  data_selected_genes <- object@data
  
  
  
  # ======================
  #### Compute distances ####
  # ======================
  # Compute distances between selected genes
  print_msg(paste0("Compute distances between selected genes."), msg_type = "INFO")
  dist_matrix_selected_genes <- qlcMatrix::corSparse(t(data_selected_genes))
  dist_matrix_selected_genes <- 1 - dist_matrix_selected_genes
  
  # Set the rownames / colnames of the distance matrix
  rownames(dist_matrix_selected_genes) <- rownames(data_selected_genes)
  colnames(dist_matrix_selected_genes) <- rownames(data_selected_genes)
  
  # The distance from a gene to itself is 'hidden'
  diag(dist_matrix_selected_genes) <- NA
  
  
  # ======================
  #### Compute distances to KNN ####
  # ======================
  print_msg(paste0("Compute distances to KNN between selected genes."), msg_type = "INFO")
  # Create a dataframe to store the DKNN values.
  # Gene_id appear both as rownames and column
  # for coding convenience
  df_dknn_selected_genes <- data.frame(
    dknn_values = rep(NA, nrow(dist_matrix_selected_genes)),
    row.names = rownames(dist_matrix_selected_genes),
    gene_id = rownames(dist_matrix_selected_genes)
  )
  
  # This list will be used to store the
  # dknn values
  l_knn_selected_genes <- list()
  
  #################### DKNN for each genes
  # Extract the DKNN for each gene
  print_msg(paste0("Computing distances to KNN."), msg_type = "INFO")
  
  for (pos in seq_len(nrow(dist_matrix_selected_genes))) {
    gene <- rownames(df_dknn_selected_genes)[pos]
    gene_dist <- dist_matrix_selected_genes[gene, ]
    # Reorder the distance values in increasing order.
    # The distance from a gene to itself (previously set to NA)
    # is placed at the end of the vector.
    gene_dist <- sort(gene_dist, na.last = T)[1:k]
    
    # Add the neigbhors to the list
    l_knn_selected_genes[[gene]] <- gene_dist
    
    # Select the kth pearson correlation values. 
    # This value corresponds to the DKNN of the gene(i)
    df_dknn_selected_genes[gene, "dknn_values"] <- gene_dist[k]
  }
  
  
  # ======================
  #### Create the input file for mcl algorithm ####
  # ======================
  print_msg(paste0("Creating the input file for MCL algorithm."), msg_type = "INFO")
  
  mcl_out_as_list_of_df <- list()
  
  for (g in names(l_knn_selected_genes)) {
    mcl_out_as_list_of_df[[g]] <- data.frame(
      src = g,
      dest = names(l_knn_selected_genes[[g]]),
      weight = l_knn_selected_genes[[g]]
    )
  }
  
  mcl_out_as_df <- do.call(rbind, mcl_out_as_list_of_df)
  
  # # A and B are added if A is in the
  # # neighborhood of B and B in the neighborhood
  # # of A
  # mcl_out_as_df <-
  #   mcl_out_as_df[duplicated(t(apply(mcl_out_as_df[, c("src", "dest")], 1, sort))), ]
  # 
  # length(unique(c(mcl_out_as_df[,"src"], mcl_out_as_df[,"dest"])))
  
  # Ensure an edge (A->B and B->A)
  # is not defined twice
  mcl_out_as_df <-
    mcl_out_as_df[!duplicated(t(apply(mcl_out_as_df[, c("src", "dest")], 1, sort))), ]
  
  
  ############# Convert distances into weights
  # scale dist between 0..1
  min_dist <- min(mcl_out_as_df$weight)
  max_dist <- max(mcl_out_as_df$weight)
  mcl_out_as_df$weight <- 
    (mcl_out_as_df$weight - min_dist) / (max_dist - min_dist)
  # Convert scaled dist to weight
  mcl_out_as_df$weight <- abs(mcl_out_as_df$weight - 1)
  
  print_stat("Graph weights (after convertion)", 
             data = mcl_out_as_df$weight, 
             msg_type = "DEBUG")
  
  
  ############# Write input files for mcl
  print_msg(paste0("Writing the input file."), msg_type = "INFO")
  data.table::fwrite(
    mcl_out_as_df,
    file = path_input_mcl,
    sep = "\t",
    eol = "\n",
    row.names = FALSE,
    col.names = FALSE
  )
  
  # Add k used to construct graph to ClusterSet object
  object@parameters <- append(object@parameters,
                              list("keep_nn" = FALSE,
                                   "k_mcl_graph" = k,
                                   "output_path" = output_path,
                                   "name" = name))
  
  print_msg(paste0("The input file saved in '", path_input_mcl, "'."), msg_type = "INFO")
  
  return(object)
}






################################################################################
#' Construct a graph based on selected genes in a ClusterSet object and their neighbors
#'
#' This function creates a graph based on the gene selected by select_genes() and their neighbors.
#' It generates an input file required for the MCL (Markov Cluster) algorithm. 
#'  
#'
#' @param object A ClusterSet object.
#' @param output_path a character indicating the path where the output files will be stored.
#' @param name a character string giving the name for the output files. If NULL, a random name is generated.
#'
#' @return The output of this function is a tab-separated file, which serves as input for the MCL algorithm.
#'
#' @examples
#' # Set verbosity to 1 to display info messages only.
#' set_verbosity(1)
#' 
#' # Create a matrix with 4 signatures
#' m <- create_4_rnd_clust()
#' 
#' # Select informative features
#' res <- select_genes(m[1:600,],
#'                     distance = "pearson",
#'                     k = 75,
#'                     noise_level = 0.3,
#'                     fdr = 1e-8,
#'                     row_sum = -Inf)
#'                     
#' # Construct a graph based on genes selected with select_genes() and their neighbors
#' keep_dbf_graph(object = res)
#'
#' @keywords internal
#' @export keep_dbf_graph
keep_dbf_graph <- function(object = NULL,
                           output_path = tempdir(),
                           name = NULL) {
  
  ## Generate output path
  # Get working directory if output_path is "."
  if (output_path == ".") {
    output_path <- getwd()
  }
  # Check if output directory exists. If not stop the command.
  if (!file.exists(output_path)) {
    stop("Output directory provided does not exist.")
  }
  
  # Create a random string if name is not provided
  if (is.null(name)){
    name <- create_rand_str()
  }
  
  # Directory and name of the principal output
  path_input_mcl <- paste(output_path, "/", name, ".input_mcl.txt", sep = "")
  path_input_mcl <- gsub(pattern = "//",
                         replacement = "/",
                         x = path_input_mcl)
  
  
  # Extract list of distances for each gene
  l_knn_selected <- object@dbf_output$all_neighbor_distances
  
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
  min_dist <- min(mcl_out_as_df$weight)
  max_dist <- max(mcl_out_as_df$weight)
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
  
  # Add parameters to ClusterSet object
  object@parameters <- append(object@parameters,
                              list("keep_nn" = TRUE,
                                   "output_path" = output_path,
                                   "name" = name))
  
  ############# Write input files for mcl
  print_msg(paste0("Writing the input file."), msg_type = "INFO")
  data.table::fwrite(
    mcl_out_as_df,
    file = path_input_mcl,
    sep = "\t",
    eol = "\n",
    row.names = FALSE,
    col.names = FALSE
  )
  print_msg(paste0("The input file saved in '", path_input_mcl, "'."), msg_type = "INFO")
  
  return(object)
}







################################################################################
#' Call MCL program for graph partitioning (internal).
#'
#' @details This function call the MCL program for graph partitioning. 
#' MCL is a graph clustering algorithm that detects clusters of nodes in a graph 
#' based on the flow of information through the edges. 
#' It is a fast and efficient algorithm that can be used for clustering large-scale graphs.
#' 
#' @param object a ClusterSet object.
#' @param inflation numeric. The inflation parameter used to control the granularity of the clustering.
#' @param threads integer. The number of threads to use.
#'
#' @section Warnings: With the current implementation, this function only works
#' only on UNIX-like plateforms.
#'
#' MCL should be installed. Run the \code{install_mcl} command to install it.
#'
#' @references
#' - Van Dongen S. (2000) A cluster algorithm for graphs. National
#' Research Institute for Mathematics and Computer Science in the 1386-3681.
#'
#' @return a ClusterSet object with the updated parameters.
#' 
#' @export mcl_system_cmd
mcl_system_cmd <- function(object = NULL,
                           inflation = inflation,
                           threads = 1) {
  ## testing the system
  if (.Platform$OS.type == "windows") {
    stop("--> A unix-like OS is required to launch the MCL program.")
  }else {
    print_msg("Running mcl_system_cmd() under a unix-like system.", msg_type = "DEBUG")
  }
  
  ## Testing mcl installation
  if (Sys.which("mcl") != "") {
    print_msg("Found MCL program in the path...", msg_type = "DEBUG")
    mcl_dir <- ""
  } else {
    mcl_dir <- Sys.glob(file.path(path.expand('~'), 
                                  ".scigenex", "mcl*",
                                  "src", "shmcl"))
    if(length(mcl_dir) == 0){
      print_msg("MCL was not found in the PATH nor in  ~/.scigenex. Installing in ~/.scigenex")  
      install_mcl()
    }else{
      print_msg("MCL was found in the ~/.scigenex directory.")
    }
    
  }

  name <- object@parameters$name
  input_path <- object@parameters$output_path
  
  if (get_verbosity() > 0) {
    verb <- ""
  } else {
    verb <- "-V all "
  }
  
  i <- paste("-I ", as.character(round(inflation, 1)), sep = "")
  
  threads <- paste("-te", threads, sep = " ")
  
  ## launching mcl program
  if(mcl_dir == "" | mcl_dir == " "){
    mcl_path <- "mcl"
  }else{
    mcl_path <- file.path(mcl_dir, "mcl")
  }
  
  cmd <- paste0(mcl_path, 
                " ",
                input_path,
                "/",
                name,
                ".input_mcl.txt ",
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

  print_msg(paste0("Running mcl command: ", cmd), msg_type = "DEBUG")
  
  system(cmd)
  
  print_msg("MCL step is finished.", msg_type = "DEBUG")
  print_msg(paste0("creating file : ",
                   file.path(
                     getwd(), paste(name, ".mcl_out.txt", sep = "")
                   )),
            msg_type = "DEBUG")
  
  # Add inflation parameter to ClusterSet object
  object@parameters <- append(object@parameters,
                              list("inflation" = inflation))
  
  return(object)
}
