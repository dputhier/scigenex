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
#' @param algorithm The algorihm to use for clustering: MCL", "infomap", "louvain" or "walktrap".
#' @param threads An integer value indicating the number of threads to use for MCL.
#' @param output_path a character indicating the path where the output files will be stored.
#' @param name a character string giving the name for the output files. If NULL, a random name is generated.
#' @param keep_nn Deprecated. Use 'method' instead.
#' @param louv_resolution Resolution of Louvain algorithm if chosen.
#' @param walktrap_step The length of the random walks to perform for walktrap algorithm.
#' @param infomap_nb nb.trials parameter for igraph::cluster_infomap().
#' @param infomap_modularity modularity parameter for igraph::cluster_infomap().
#' @param mcl_dir A path to MCL (if issue with automated installation).
#' @param delete_output_file Whether to delete output files (.graph_input.txt, .graph_out.txt).
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
                            algorithm=c("MCL", "infomap", "louvain", "walktrap"),
                            threads = 1,
                            output_path = NULL,
                            name = NULL,
                            keep_nn = FALSE,
                            louv_resolution=5,
                            walktrap_step=4,
                            infomap_nb=10,
                            infomap_modularity=TRUE,
                            mcl_dir=NULL,
                            delete_output_file=FALSE) {
  
  print_msg("Retrieving args.", msg_type = "DEBUG")
  
  method <- match.arg(method)
  algorithm <- match.arg(algorithm)
  
  print_msg("Preparing output path", msg_type = "DEBUG")
  
  if(is.null(output_path)){
    output_path <- tempdir()
    dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
  }else{
    if (!dir.exists(output_path)) {
      dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
    }
  }
  
  if(keep_nn){
    print_msg("The use of keep_nn=TRUE is deprecated. Use 'method' argument set to 'reciprocal_neighborhood'.")
    method <- "reciprocal_neighborhood"
  }
  
  print_msg("Preparing graph as a file", msg_type = "DEBUG")
  
  if (method == "reciprocal_neighborhood") {
    object <- do_reciprocal_neighbor_graph(object = object,
                             output_path = output_path,
                             name = name)
  } else {
    object <- do_closest_neighbor_graph(object = object,
                                        k = s,
                                        output_path = output_path,
                                        name = name)
  }
  
  #---------------- Run graph partitionning
  if(algorithm == "MCL"){

    print_msg("MCL algorithm has been selected.")
    object <- mcl_system_cmd(object = object,
                             inflation = inflation,
                             threads = threads,
                             mcl_dir=mcl_dir)
    
    clust_out_file <- file.path(object@parameters$output_path,
                              paste0(object@parameters$name, ".graph_out.txt"))
    
    
  }else if(algorithm == "louvain"){
  
    print_msg("Louvain algorithm has been selected.")
    object <- call_louvain_clusterset(object = object,
                                      resolution = louv_resolution)
    
    clust_out_file <- file.path(object@parameters$output_path,
                                paste0(object@parameters$name, ".graph_out.txt"))
    
  }else if(algorithm == "walktrap"){
    
    print_msg("The walktrap algorithm has been selected.")

    object <- call_walktrap_clusterset(object, step=walktrap_step)
    
    return(object)
    
    clust_out_file <- file.path(object@parameters$output_path,
                                paste0(object@parameters$name, ".graph_out.txt"))
  
  }else if(algorithm == "infomap"){
    
    print_msg("The infomap algorithm has been selected.")
    
    object <- call_infomap_clusterset(object, 
                                      nb.trials=infomap_nb,
                                      modularity=infomap_modularity)
    return(object)
    
    clust_out_file <- file.path(object@parameters$output_path,
                                paste0(object@parameters$name, ".graph_out.txt"))
  }

  print_msg(paste0("Reading graph clustering output file."), msg_type = "DEBUG")
  algo_cluster <- readLines(clust_out_file)
  algo_cluster <- strsplit(algo_cluster, "\t")
  names(algo_cluster) <- seq(1, length(algo_cluster))

  print_msg("Adding clusters to a ClusterSet object.", msg_type = "DEBUG")
  object@gene_clusters <- algo_cluster
  
  print_msg("Update gene_cluster_metadata slots.", msg_type = "DEBUG")
  object@gene_clusters_metadata <- list("cluster_id" = as.numeric(names(object@gene_clusters)),
                                        "number" = max(as.numeric(names(object@gene_clusters))),
                                        "size" = unlist(lapply(object@gene_clusters, length)))
  
  print_msg("Updating @data slot.", msg_type = "DEBUG")
  object@data <- object@data[unlist(object@gene_clusters, use.names = FALSE), ]
  
  print_msg("Computing centers.", msg_type = "DEBUG")
  nb_clusters = length(names(object@gene_clusters))
  centers <- matrix(ncol = ncol(object@data),
                    nrow = nb_clusters)
  colnames(centers) <- colnames(object@data)
  rownames(centers) <- names(object@gene_clusters)
  
  for (i in 1:nb_clusters) {

      centers[i, ] <- apply(object@data[object@gene_clusters[[i]], , drop=FALSE],
                            2, mean,
                            na.rm = TRUE)
    
  }
  
  object@dbf_output$center <- centers
  rownames(object@dbf_output$center) <- names(object@gene_clusters)
  
  object@cells_metadata <- data.frame("cells_barcode" = colnames(object@data),
                                      row.names = colnames(object@data))
  
  print_msg(paste0("Number of clusters found: ", nclust(object)))
  
  if(delete_output_file){
    unlink(Sys.glob(paste0(file.path(object@parameters$output_path, 
                                     object@parameters$name),
                           "*")), 
           force = TRUE)
    object@parameters$output_path <- ""
    object@parameters$name <- ""
  }
    
  return(object)
}


################################################################################
#' Construct a new graph for a ClusterSet object (internal)
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
#' res <- do_closest_neighbor_graph(object = res, k = 5)
#' @keywords internal
#' @export do_closest_neighbor_graph
do_closest_neighbor_graph <- function(object = NULL,
                                      k = 5,
                                      output_path = NULL,
                                      name = NULL) {
  
  print_msg("Checking object format.", msg_type = "DEBUG")
  
  check_format_cluster_set(object)
  
  print_msg("Preparing output file.", msg_type = "DEBUG")
  
  if (is.null(name)){
    name <- create_rand_str()
  }
  
  if(is.null(output_path)){
    output_path <- tempdir()
    dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
  }
  
  if (output_path == ".") {
    output_path <- getwd()
  }

  if (!dir.exists(output_path)) {
    dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
  }
  
  path_input_graph <- file.path(output_path, 
                                paste0(name, ".graph_input.txt"))

  print_msg("Retrieving selected features.", msg_type = "DEBUG")
  
  selected_genes <- object@gene_clusters$`1`
  data_selected_genes <- object@data
  
  print_msg("Computing distances between selected genes.", msg_type = "INFO")
  dist_matrix_selected_genes <- qlcMatrix::corSparse(t(data_selected_genes))
  dist_matrix_selected_genes <- 1 - dist_matrix_selected_genes
  rownames(dist_matrix_selected_genes) <- rownames(data_selected_genes)
  colnames(dist_matrix_selected_genes) <- rownames(data_selected_genes)
  
  # The distance from a gene to itself is 'hidden'
  diag(dist_matrix_selected_genes) <- NA
  
  print_msg("Compute distances to KNN between selected genes.", msg_type = "INFO")
  df_dknn_selected_genes <- data.frame(
    dknn_values = rep(NA, nrow(dist_matrix_selected_genes)),
    row.names = rownames(dist_matrix_selected_genes),
    gene_id = rownames(dist_matrix_selected_genes)
  )
  
  # A list to store the dknn values
  l_knn_selected_genes <- list()
  
  print_msg("Computing distances to KNN.", msg_type = "INFO")
  
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
  
  print_msg("Creating the input file for graph partitioning.", msg_type = "INFO")
  
  mcl_out_as_list_of_df <- list()
  
  for (g in names(l_knn_selected_genes)) {
    mcl_out_as_list_of_df[[g]] <- data.frame(
      src = g,
      dest = names(l_knn_selected_genes[[g]]),
      weight = l_knn_selected_genes[[g]]
    )
  }
  
  mcl_out_as_df <- do.call(rbind, mcl_out_as_list_of_df)
  
  print_msg("Deleting reciprocal edges.", msg_type = "INFO")
  
  mcl_out_as_df <-
    mcl_out_as_df[!duplicated(t(apply(mcl_out_as_df[, c("src", "dest")], 1, sort))), ]
  
  print_msg("Converting distances into weights.", msg_type = "INFO")
  
  min_dist <- min(mcl_out_as_df$weight)
  max_dist <- max(mcl_out_as_df$weight)
  mcl_out_as_df$weight <- 
    (mcl_out_as_df$weight - min_dist) / (max_dist - min_dist)
  
  mcl_out_as_df$weight <- abs(mcl_out_as_df$weight - 1)
  
  print_stat("Graph weights (after convertion)", 
             data = mcl_out_as_df$weight, 
             msg_type = "DEBUG")
  

  print_msg("Writing the input file.", msg_type = "INFO")
  
  data.table::fwrite(
    mcl_out_as_df,
    file = path_input_graph,
    sep = "\t",
    eol = "\n",
    row.names = FALSE,
    col.names = FALSE
  )
  
  print_msg("Storing analysis parameters in ClusterSet object.", msg_type = "DEBUG")
  
  object@parameters <- append(object@parameters,
                              list("keep_nn" = FALSE,
                                   "k_graph" = k,
                                   "output_path" = output_path,
                                   "name" = name))
 
   
  print_msg(paste0("Input file saved in :", path_input_graph), msg_type = "INFO")

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
#' do_reciprocal_neighbor_graph(object = res)
#'
#' @keywords internal
#' @export do_reciprocal_neighbor_graph
do_reciprocal_neighbor_graph <- function(object = NULL,
                                         output_path = NULL,
                                         name = NULL) {
  
  print_msg("Checking object format.", msg_type = "DEBUG")
  
  check_format_cluster_set(object)
  
  print_msg("Preparing output file.", msg_type = "DEBUG")
  
  if(is.null(output_path)){
    output_path <- tempdir()
    dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
  }

  if (output_path == ".") {
    output_path <- getwd()
  }

  if (!dir.exists(output_path)) {
    dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
  }
  
  if (is.null(name)){
    name <- create_rand_str()
  }
  
  path_input_graph <- file.path(output_path, 
                                paste0(name, ".graph_input.txt"))
  
  print_msg("Extracting distances to neighbors.", msg_type = "DEBUG")
  l_knn_selected <- object@dbf_output$all_neighbor_distances
  
  print_msg("Creating the input file for MCL algorithm.", msg_type = "INFO")
  
  mcl_out_as_list_of_df <- list()
  
  for (g in names(l_knn_selected)) {
    mcl_out_as_list_of_df[[g]] <- data.frame(
      src = g,
      dest = names(l_knn_selected[[g]]),
      weight = l_knn_selected[[g]]
    )
  }
  
  mcl_out_as_df <- do.call(rbind, mcl_out_as_list_of_df)
  
  print_msg("Convert distances into weights.", msg_type = "INFO")
  
  min_dist <- min(mcl_out_as_df$weight)
  max_dist <- max(mcl_out_as_df$weight)
  mcl_out_as_df$weight <- 
    (mcl_out_as_df$weight - min_dist) / (max_dist - min_dist)
  mcl_out_as_df$weight <- abs(mcl_out_as_df$weight - 1)
  
  print_stat("Graph weights (after convertion)", 
             data = mcl_out_as_df$weight, 
             msg_type = "DEBUG")
  
  print_msg("Selecting only reciprocal neighborhood.", msg_type = "INFO")
  mcl_out_as_df <-
    mcl_out_as_df[duplicated(t(apply(mcl_out_as_df[, c("src", "dest")], 1, sort))), ]
  
  print_msg("Deleting reciprocal edges.", msg_type = "DEBUG")
  mcl_out_as_df <-
    mcl_out_as_df[!duplicated(t(apply(mcl_out_as_df[, c("src", "dest")], 1, sort))), ]
  
  print_msg("Adding parameters to the ClusterSet object.", msg_type = "DEBUG")
  object@parameters <- append(object@parameters,
                              list("keep_nn" = TRUE,
                                   "output_path" = output_path,
                                   "name" = name))
  
  print_msg(paste0("Writing the input file."), msg_type = "INFO")
  
  data.table::fwrite(
    mcl_out_as_df,
    file = path_input_graph,
    sep = "\t",
    eol = "\n",
    row.names = FALSE,
    col.names = FALSE
  )
  print_msg(paste0("Input file saved in :", path_input_graph), msg_type = "INFO")
  
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
#' @param mcl_dir The directory were mcl is installed (guessed otherwise).
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
                           threads = 1,
                           mcl_dir=NULL) {
  ## testing the system
  if (.Platform$OS.type == "windows") {
    stop("--> A unix-like OS is required to launch the MCL program.")
  }else {
    print_msg("Running mcl_system_cmd() under a unix-like system.", msg_type = "DEBUG")
  }
  
  ## Testing mcl installation
  if(is.null(mcl_dir)){
    
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
        mcl_dir <- Sys.glob(file.path(path.expand('~'), 
                                      ".scigenex", "mcl*",
                                      "src", "shmcl"))
      }else{
        print_msg("MCL was found in the .scigenex home directory.")
      }
      
    }
    
  }
  
  
  print(paste0("MCL path : ", mcl_dir))

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
  print_msg(paste0("mcl_dir:", mcl_dir))
  
  if(mcl_dir == ""){
    mcl_dir <- "mcl"
  }else if(mcl_dir == " "){
    mcl_dir <- "mcl"
  }else{
    mcl_dir <- file.path(mcl_dir, "mcl")
  }
  
  cmd <- paste0(mcl_dir, 
                " ",
                input_path,
                "/",
                name,
                ".graph_input.txt ",
                i,
                " --abc -o ",
                input_path,
                "/",
                name,
                ".graph_out.txt ",
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
                     getwd(), paste(name, ".graph_out.txt", sep = "")
                   )),
            msg_type = "DEBUG")
  
  # Add inflation parameter to ClusterSet object
  object@parameters <- append(object@parameters,
                              list("inflation" = inflation,
                                   "algorithm"="mcl"))
  
  return(object)
}

#' Perform Community Detection Using the Louvain Algorithm (internal function)
#'
#' Applies the Louvain algorithm for community detection on a graph derived 
#' from the input data in a `ClusterSet` object. Outputs the resulting clusters 
#' to a specified file and updates the parameters in the `ClusterSet` object.
#'
#' @param object A `ClusterSet` object. 
#' 
#' @param resolution The level of resolution.
#'
#' @return The modified `ClusterSet` object with updated parameters, including 
#'   the Walktrap algorithm results
#'
#' @details This function reads a graph from the input file defined in the `ClusterSet` 
#'   object, processes the graph using the Louvain community detection algorithm 
#'   from the `igraph` package, and return the clusters enclose in the `ClusterSet`.
#'
#' @examples
#' # Restrict vebosity to info messages only.
#' set_verbosity(1)
#' # Load a dataset
#' load_example_dataset("7871581/files/pbmc3k_medium")
#' # Select informative genes
#' res <- select_genes(pbmc3k_medium,
#'                     distance = "pearson",
#'                     row_sum=5)
#' 
#' # Cluster informative features
#' res <- gene_clustering(res, method="closest_neighborhood",
#'                        algorithm="louvain")
#' 
#' @importFrom data.table fread
#' @importFrom igraph graph_from_data_frame cluster_louvain
#' @export
call_louvain_clusterset <- function(object, 
                                    resolution = 1){

  name <- object@parameters$name
  input_dir <- object@parameters$output_path
  input_path <- file.path(input_dir, 
                           paste0(name, ".graph_input.txt"))
  
  g <- data.table::fread(input_path, 
                         header = FALSE)
  g <- as.data.frame(g)
  colnames(g) <- c("src", "dest", "weight")
  
  g <- g[!duplicated(paste0(g$src, g$dest)), ]
  
  g <- igraph::graph_from_data_frame(g, directed = FALSE)
   
  print_msg("Calling Louvain algorithm.")
  clust <- igraph::cluster_louvain(g, resolution=resolution)
  
  clust <- split(clust$names, clust$membership)
  clust <- lapply(clust, paste0, collapse="\t")
  clust <- unlist(clust)
  
  out_path <- file.path(input_dir, 
                        paste0(name, ".graph_out.txt"))
  
  print_msg(paste0("Output writen to ", out_path), msg_type = "DEBUG")
  
  cat(clust, file=out_path, sep = "\n")
  
  # Add inflation parameter to ClusterSet object
  params <- object@parameters
  params$resolution <- resolution
  params$algorithm <- "louvain"
  object@parameters <- params 
  
 return(object)
}

#' Perform Community Detection Using the Walktrap Algorithm (internal function)
#'
#' Applies the Walktrap algorithm for community detection on a graph derived 
#' from the input data in a `ClusterSet` object. Outputs the resulting clusters 
#' to a specified file and updates the parameters in the `ClusterSet` object.
#'
#' @param object A `ClusterSet` object.
#' @param steps The length of the random walks to perform.
#' @param merges See igraph::cluster_walktrap().
#' @param  modularity See igraph::cluster_walktrap().
#' @param  membership See igraph::cluster_walktrap().
#' @return The modified `ClusterSet` object with updated parameters, including 
#'   the Walktrap algorithm results
#'
#' @details This function reads a graph from the input file defined in the `ClusterSet` 
#'   object, processes the graph using the Walktrap community detection algorithm 
#'   from the `igraph` package, and return the clusters enclose in the `ClusterSet`.
#'
#' @examples
#' # Restrict vebosity to info messages only.
#' set_verbosity(1)
#' # Load a dataset
#' load_example_dataset("7871581/files/pbmc3k_medium")
#' # Select informative genes
#' res <- select_genes(pbmc3k_medium,
#'                     distance = "pearson",
#'                     row_sum=5)
#' 
#' # Cluster informative features
#' res <- gene_clustering(res, method="closest_neighborhood",
#'                        algorithm="walktrap")
#' 
#' @importFrom data.table fread
#' @importFrom igraph graph_from_data_frame cluster_walktrap
#' @export
call_walktrap_clusterset <- function(object, 
                                     steps = 4,
                                     merges = FALSE,
                                     modularity = FALSE,
                                     membership = TRUE){
  
  name <- object@parameters$name
  input_dir <- object@parameters$output_path
  input_path <- file.path(input_dir, 
                          paste0(name, ".graph_input.txt"))
  
  g <- data.table::fread(input_path, 
                         header = FALSE)
  g <- as.data.frame(g)
  colnames(g) <- c("src", "dest", "weight")
  
  g <- g[!duplicated(paste0(g$src, g$dest)), ]
  
  g <- igraph::graph_from_data_frame(g, directed = FALSE)
  
  print_msg("Calling walktrap algorithm.")

  clust <- igraph::cluster_walktrap(g, 
                                    steps=steps,
                                    merges=merges,
                                    modularity=modularity,
                                    membership=membership)

  clust <- split(clust$names, clust$membership)
  clust <- lapply(clust, paste0, collapse="\t")
  clust <- unlist(clust)
  
  out_path <- file.path(input_dir, 
                        paste0(name, ".graph_out.txt"))
  
  print_msg(paste0("Output writen to ", out_path), msg_type = "DEBUG")
  
  cat(clust, file=out_path, sep = "\n")
  
  # Add inflation parameter to ClusterSet object
  params <- object@parameters
  params$steps <- steps
  params$merges <- merges
  params$modularity <- modularity
  params$membership <- membership
  params$algorithm <- "walktrap"
  object@parameters <- params 
  
  return(object)
}



#' Perform Community Detection Using the infomap Algorithm (internal function)
#'
#' Applies the infomap algorithm for community detection on a graph derived 
#' from the input data in a `ClusterSet` object. Outputs the resulting clusters 
#' to a specified file and updates the parameters in the `ClusterSet` object.
#'
#' @param object A `ClusterSet` object.
#' @param modularity See igraph::cluster_infomap().
#' @param nb.trials See igraph::cluster_infomap().
#'
#' @return The modified `ClusterSet` object with updated parameters, including 
#'   the infomap algorithm results
#'
#' @details This function reads a graph from the input file defined in the `ClusterSet` 
#'   object, processes the graph using the infomap community detection algorithm 
#'   from the `igraph` package, and return the clusters enclose in the `ClusterSet`.
#'
#' @examples
#' # Restrict vebosity to info messages only.
#' set_verbosity(1)
#' # Load a dataset
#' load_example_dataset("7871581/files/pbmc3k_medium")
#' # Select informative genes
#' res <- select_genes(pbmc3k_medium,
#'                     distance = "pearson",
#'                     row_sum=5)
#' 
#' # Cluster informative features
#' res <- gene_clustering(res, method="closest_neighborhood",
#'                        algorithm="infomap")
#' 
#' @importFrom data.table fread
#' @importFrom igraph graph_from_data_frame cluster_walktrap
#' @export
call_infomap_clusterset <- function(object, 
                                    modularity = FALSE,
                                    nb.trials = 10){
  
  name <- object@parameters$name
  input_dir <- object@parameters$output_path
  input_path <- file.path(input_dir, 
                          paste0(name, ".graph_input.txt"))
  
  g <- data.table::fread(input_path, 
                         header = FALSE)
  g <- as.data.frame(g)
  colnames(g) <- c("src", "dest", "weight")
  
  g <- g[!duplicated(paste0(g$src, g$dest)), ]
  
  g <- igraph::graph_from_data_frame(g, directed = FALSE)
  
  print_msg("Calling infomap algorithm.")
  
  clust <- igraph::cluster_infomap(g, 
                                   modularity=modularity,
                                   nb.trials=nb.trials)
  

  clust <- split(clust$names, clust$membership)
  clust <- lapply(clust, paste0, collapse="\t")
  clust <- unlist(clust)
  
  out_path <- file.path(input_dir, 
                        paste0(name, ".graph_out.txt"))
  
  print_msg(paste0("Output writen to ", out_path), msg_type = "DEBUG")
  
  cat(clust, file=out_path, sep = "\n")
  
  # Add inflation parameter to ClusterSet object
  params <- object@parameters
  params$modularity <- modularity
  params$nb.trials <- nb.trials
  params$algorithm <- "infomap"
  object@parameters <- params 
  
  return(object)
}
