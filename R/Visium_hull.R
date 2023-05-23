###########################################################
#
# The getFlippedTissueCoordinates() function
#
###########################################################

#' @name getFlippedTissueCoordinates
#' @title Seurat object internally store spot coordinates (see Seurat::GetTissueCoordinates()). 
#' However, at least in the case of Visium, data are flipped and rotated before SpatialDimPlot. 
#' This function  return the rotated/flipped tissue Coordinates from a Seurat object.
#' @param seurat_obj a seurat object with tissue coordinates.
#' @param as_data_frame return x/y coords as data.frame. Default to SeuratObject.
#' @return a seurat object with slots/metadata $x_coord and $y_coord (or a dataframe if as_data_frame is TRUE).
#' @seealso display_visium_hull
#' @examples
#' load_example_dataset("7870305/files/lymph_node_tiny_2")
#' lymph_node_tiny <- getFlippedTissueCoordinates(lymph_node_tiny_2)
#' df <- getFlippedTissueCoordinates(lymph_node_tiny_2, as_data_frame=TRUE)
#' @export getFlippedTissueCoordinates
getFlippedTissueCoordinates <- function(seurat_obj, 
                                        as_data_frame=FALSE){
  # Get coord and flip
  coord_spot <- GetTissueCoordinates(seurat_obj)[,2:1] # rotation
  colnames(coord_spot) <- c("x", "y")
  min_coord_y <- min(coord_spot$y)
  max_coord_y <- max(coord_spot$y)
  coord_spot$y <- -coord_spot$y + 2*min_coord_y + max_coord_y-min_coord_y
  
  # prepare return
  coord_spot_x <- coord_spot$x
  names(coord_spot_x) <- colnames(seurat_obj)
  coord_spot_y <- coord_spot$x
  names(coord_spot_y) <- colnames(seurat_obj)
  if(!as_data_frame){
    seurat_obj$x_coord <- coord_spot$x 
    seurat_obj$y_coord <-coord_spot$x   
    return(seurat_obj)
  }else{
    m <- data.frame("x"=coord_spot$x, 
                    "y"=coord_spot$y,
                    "row.names" = colnames(seurat_obj))
  }
  
}

###########################################################
#
# The center_and_rotate() function
#
###########################################################
# Centers and rotates a dataset around a specified center point and angle
# (internal function)
#
# This function centers and rotates a given dataset around a specified center point
# and angle. The resulting dataset is returned as a data frame.
#
# param data A data frame to be centered and rotated
# param center_x The x-coordinate of the center point
# param center_y The y-coordinate of the center point
# param angle The angle (in degrees) to rotate the dataset
# param x The name of the column representing the x-coordinate in the data frame
# param y The name of the column representing the y-coordinate in the data frame
#
# return The centered and rotated data frame
#
# examples
# # Generate sample data
# data <- data.frame(x = c(1, 2, 3, 4, 5), y = c(2, 3, 4, 5, 6))
#
# # Center and rotate the data by 45 degrees around point (3, 4)
# p <- center_and_rotate(data, center_x = 3, center_y = 4, angle = 45)
center_and_rotate <- function(data, center_x, center_y, angle, x="x", y="y"){
  
  data[, x] <- data[, x] - center_x
  data[, y] <- data[, y] - center_y
  angle_rad <- pi * angle / 180
  
  # rotate
  rotation_mat <- matrix(c(cos(angle_rad),
                           -sin(angle_rad),
                           sin(angle_rad),
                           cos(angle_rad)),
                         ncol=2, byrow = T)
  
  rotated_mat <- as.matrix(data[, c(x,y)]) %*% rotation_mat
  colnames(rotated_mat) <- c(x,y)
  data[,c(x,y)] <- rotated_mat[,c(x,y)]
  
  return(data)
}

###########################################################
#
# The get_neighbor_class() function
#
###########################################################
# Find the class of the nearest point (internal function see visium_hull).
#
# This function takes in a centered data frame of coordinates and returns the class of the
# nearest point lying on the right of the center. On x The point is expected to be between 0 and
# step_x * 2 + (step_x * 2) * delta. On y it is expected to be between -step_y/2 and step_y/2. 
# The neighbor class is return.
#
# param data_cr A data frame of coordinates to search for the nearest point in the grid
# param step_x The step size in the x direction for the grid
# param step_y The step size in the y direction for the grid
# param delta The tolerance to use when searching for a point in the grid
#
# return The class of the nearest point in the grid, or NULL if no such point is found
get_neighbor_class <- function(data_cr,
                               step_x=3,
                               step_y=5,
                               delta=0.3,
                               k="k"){
  
  # We are looking for a point on the right (data_cr$x > 0)
  # That it rather close (data_cr$x < step_x * 2 + (step_x * 2) * delta)
  test_x <- data_cr$x > 0 & data_cr$x < step_x * 2 + (step_x * 2) * delta
  
  # We are looking for a point whose y value is enclosed around O
  #  -step_y/2 < y < step_y/2  
  test_y <- data_cr$y > 0 - (step_y/2)  & data_cr$y < 0 + (step_y/2)
  
  # The point that meet that criteria
  hit <- data_cr[test_x & test_y, , drop=FALSE]
  
  # The class of that point
  if(nrow(hit) > 0){
    neighbor_class <- hit[, k]
    neighbor_name <- rownames(hit)
  }else{
    neighbor_class <- NULL
    neighbor_name <- NULL
  }
  return(list(neighbor_class=neighbor_class, 
              neighbor_name=neighbor_name))
}

###########################################################
#
# merging_neighborhood()
#
###########################################################
# name merging_neighborhood
# title Internal function. Whether a barcode/spot is in the known
# or new neighborhood.
merging_neighborhood <- function(neighborhood=NULL, spot_1=NULL, spot_2=NULL){
  
  
  if(is.null(spot_1) | is.null(spot_2))
    return(neighborhood)
  
    test_1 <- unlist(lapply(lapply(neighborhood, "==", spot_1), any))
    test_2 <- unlist(lapply(lapply(neighborhood, "==", spot_2), any))
  
    
    if(which(test_1) != which(test_2)){
      
      print_msg("Merging two neighborhoods...", msg_type = "DEBUG")
      neighborhood[[which(test_1)]] <- unique(c(neighborhood[[which(test_1)]], 
                                          neighborhood[[which(test_2)]]))
      neighborhood <- neighborhood[-which(test_2)]
      names(neighborhood) <- 1:length(neighborhood)
    }

  return(neighborhood)
}

###########################################################
#
# The Visium_hull() function
#
###########################################################

# name visium_hull
# title Internal function. Compute location of segments used to create an orthogonal 
# hull around points related to a particular class.
# data a data.frame.
# param x The column name storing the x coord
# param y The column name storing the y coord
# param k The column name storing the classes (0 not part of the class of interest, 1 part of the class of interest)
# param size_y The size of the square (y axis)
# param size_x The size of the square (x axis)
# param step_y The distance between two points on the y axis.
# param step_x The distance between two points on the x axis.
# param delta Add more or less flexibility to search for neighbor points
# keywords hull, spatial transcriptomics, visium
# return A dataframe with coordinates x, y, xend, yend (see geom_segments).
# seealso display_visium_hull
# examples
# # See display_visium_hull()
visium_hull <- function(data,
                        x="x",
                        y="y",
                        k="k",
                        size_x=2.3,
                        size_y=2,
                        step_x=2.6,
                        step_y=2.4,
                        delta=0.3){
  
  
  print_msg("Creating a dataframe to store output", msg_type="INFO")
  
  df_coord <- data.frame(matrix(NA, ncol=4, nrow=6))
  rownames(df_coord) <- sapply("region_name",
                               paste0, "_",
                               c("north_west", "north_east",
                                 "south_east", "south_west",
                                 "west", "east"))
  
  colnames(df_coord) <- c("x1", "x2", "y1", "y2")
  
  print_msg("Creating a list of neighborhoods.", msg_type="INFO")
  
  neighborhoods <- as.list(rownames(data[data$k==1,]))
  
  print_msg("Looping over the points.", msg_type="INFO")
  
  for(p in 1:nrow(data)){
    
    pts_name <- rownames(data)[p]
    
    x_p <- data[p, x]
    y_p <- data[p, y]
    
    pt_class <- data[p, k]
    
    p_name <- paste0(pts_name, "||", x_p, "_", y_p)
    
    
    if(pt_class == 1){
      
      ## West
      data_cr <- center_and_rotate(data, center_x = x_p, center_y = y_p, angle=180, x=x, y=y)
      neighbor_info <- get_neighbor_class(data_cr, 
                                           step_x=step_x,
                                           step_y=step_y,
                                           delta=delta,
                                           k=k)
      
      neighbor_class <- neighbor_info$neighbor_class
      neighbor_name <- neighbor_info$neighbor_name
      
      if(isTRUE(neighbor_class == 0) | is.null(neighbor_class)){

        df_coord[paste0(p_name, "west"),] <- c("x1" = x_p - size_x,
                                               "x2" = x_p - size_x,
                                               "y1" = y_p -  size_y,
                                               "y2" = y_p +  size_y
        )
      }else if(neighbor_class == 1){
        neighborhoods <- merging_neighborhood(neighborhoods, 
                                             pts_name, 
                                             neighbor_name)
      }
      
      
      ## Check_east
      data_cr <- center_and_rotate(data, center_x = x_p, center_y = y_p, angle=0, x=x, y=y)
      neighbor_info <- get_neighbor_class(data_cr, 
                                           step_x=step_x,
                                           step_y=step_y,
                                           delta=delta,
                                           k=k)

      
      neighbor_class <- neighbor_info$neighbor_class
      neighbor_name <- neighbor_info$neighbor_name
      
      if(isTRUE(neighbor_class == 0) | is.null(neighbor_class)){
        
        df_coord[paste0(p_name, "east"),] <- c("x1" = x_p + size_x,
                                               "x2" = x_p + size_x,
                                               "y1" = y_p -  size_y,
                                               "y2" = y_p +  size_y
        )
      }else if(neighbor_class == 1){
        neighborhoods <- merging_neighborhood(neighborhoods, 
                                             pts_name, 
                                             neighbor_name)
      }
      
      ## North_east
      data_cr <- center_and_rotate(data, center_x = x_p, center_y = y_p, angle=60, x=x, y=y)
      neighbor_info <- get_neighbor_class(data_cr, 
                                           step_x=step_x,
                                           step_y=step_y,
                                           delta=delta,
                                           k=k)
      
      neighbor_class <- neighbor_info$neighbor_class
      neighbor_name <- neighbor_info$neighbor_name

      if(isTRUE(neighbor_class == 0) | is.null(neighbor_class)){
        
        df_coord[paste0(p_name, "north_east"),] <- c("x1" = x_p,
                                                     "x2" = x_p + size_x,
                                                     "y1" = y_p +  size_y,
                                                     "y2" = y_p +  size_y
        )
      }else if(neighbor_class == 1){
        neighborhoods <- merging_neighborhood(neighborhoods, 
                                             pts_name, 
                                             neighbor_name)
      }
      
      ## North_west
      data_cr <- center_and_rotate(data, center_x = x_p, center_y = y_p, angle=120, x=x, y=y)
      neighbor_info <- get_neighbor_class(data_cr, 
                                           step_x=step_x,
                                           step_y=step_y,
                                           delta=delta,
                                           k=k)
      
      neighbor_class <- neighbor_info$neighbor_class
      neighbor_name <- neighbor_info$neighbor_name
      
      if(isTRUE(neighbor_class == 0) | is.null(neighbor_class)){
        
        df_coord[paste0(p_name, "north_west"),] <- c("x1" = x_p -  size_x,
                                                     "x2" = x_p ,
                                                     "y1" = y_p +  size_y,
                                                     "y2" = y_p +  size_y
        )
      }else if(neighbor_class == 1){
        neighborhoods <- merging_neighborhood(neighborhoods, 
                                             pts_name, 
                                             neighbor_name)
      }
      
      ## South_west
      data_cr <- center_and_rotate(data, center_x = x_p, center_y = y_p, angle=240, x=x, y=y)
      neighbor_info <- get_neighbor_class(data_cr, 
                                          step_x=step_x,
                                          step_y=step_y,
                                          delta=delta,
                                          k=k)
      
      neighbor_class <- neighbor_info$neighbor_class
      neighbor_name <- neighbor_info$neighbor_name
      
      if(isTRUE(neighbor_class == 0) | is.null(neighbor_class)){
        
        df_coord[paste0(p_name, "south_west"),] <- c("x1" = x_p - size_x,
                                                     "x2" = x_p,
                                                     "y1" = y_p - size_y,
                                                     "y2" = y_p - size_y
        )
      }else if(neighbor_class == 1){
        neighborhoods <- merging_neighborhood(neighborhoods, 
                                             pts_name, 
                                             neighbor_name)
      }
      
      ## South_east
      data_cr <- center_and_rotate(data, center_x = x_p, center_y = y_p, angle=300, x=x, y=y)
      neighbor_info <- get_neighbor_class(data_cr, 
                                          step_x=step_x,
                                          step_y=step_y,
                                          delta=delta,
                                          k=k)
      
      neighbor_class <- neighbor_info$neighbor_class
      neighbor_name <- neighbor_info$neighbor_name
      
      if(isTRUE(neighbor_class == 0) | is.null(neighbor_class)){
        
        df_coord[paste0(p_name, "south_east"),] <- c("x1" = x_p,
                                                     "x2" = x_p +  size_x,
                                                     "y1" = y_p - size_y,
                                                     "y2" = y_p - size_y
        )
      }else if(neighbor_class == 1){
        neighborhoods <- merging_neighborhood(neighborhoods, 
                                             pts_name, 
                                             neighbor_name)
      }
    }
  }
  
  return(list(coord=stats::na.omit(df_coord), neighborhoods=neighborhoods))
}


###########################################################
#
# The display_hull() function
#
###########################################################

#' @name display_hull
#' @title Draw a hull around a region of interest of Visium data.
#' @param data A Seurat spatial object (Visium technology) or a data.frame (spot coordinates).
#' @param ident A binary vector whose length correspond to the number
#' of columns/spots in the seurat object. 1 corresponds to the spots 
#' of the class of interest for which a hull is to be displayed
#' @param size_y The size of the square (y axis)
#' @param size_x The size of the square (x axis)
#' @param step_y The distance between two points on the y axis.
#' @param step_x The distance between two points on the x axis.
#' @param delta Add more or less flexibility to search for neighbor points
#' @param color A color for the hull.
#' @param expand	For hull_type 'ggforce'. A numeric or unit vector of length one, specifying the expansion amount. 
#' @param radius	For hull_type 'ggforce'. As expand but specifying the corner radius.
#' @param hull_type The method for drawing the hull: 'ggforce' (see ggforce::geom_mark_hull) or 'wall'. 
#' If wall, segment wont be connected but this will frequently allows to obtain a better delimited hull.  
#' @param concavity A measure of the concavity of the hull. 1 is very concave while 
#' it approaches convex as it grows. Defaults to 2.
#' @param size The line width.
#' @keywords hull, spatial transcriptomics, visium
#' @return A ggplot object (with segments).
#' @export display_hull
#' @examples 
#' library(Seurat)
#' load_example_dataset("7870305/files/lymph_node_tiny_2")
#' load_example_dataset("7870305/files/lymph_node_tiny_clusters_2")
#' lymph_node_tiny_2 <- Seurat::AddModuleScore(lymph_node_tiny_2, features = lymph_node_tiny_clusters_2@gene_clusters, nbin = 10)
#' p <- Seurat::SpatialDimPlot(lymph_node_tiny_2, pt.size.factor = 4)
#' p
#' hull <- display_hull(lymph_node_tiny_2, 
#'         ident=ifelse(Seurat::Idents(lymph_node_tiny_2) %in% c(7, 8), 1, 0),
#'         delta=1, size_x=3.4, size_y=3)
#' p + hull
#' p <- plot_spatial(lymph_node_tiny_2, metadata = "Cluster3", pt_size = 5)
#' hull <- display_hull(lymph_node_tiny_2, 
#'                      ident=ifelse(Seurat::  Idents(lymph_node_tiny_2) %in% c(7, 8), 1, 0),
#'                      delta=1, size_x=3.4, size_y=3, color="black")
#' p + hull
display_hull <- function(data=NULL,
                         ident=NULL,
                         size_x=2.4,
                         size_y=2,
                         step_x=2.6,
                         step_y=2.4,
                         delta=0.3,
                         color="white",
                         hull_type=c("walls", "ggforce"),
                         expand = unit(0.05, "mm"),
                         radius = unit(0.5, "mm"),
                         concavity=1,
                         size=0.7
                         ){
  
  hull_type <- match.arg(hull_type)
  
  if(inherits(data, "Seurat")){
    if(length(ident) != length(colnames(data)))
      print_msg("The 'ident' vector should heave same length as the number of spots.",
                msg_type = "STOP")
    
    if(! "Spatial" %in% names(data))
      print_msg("Please provide a Spatial Seurat object.",
                msg_type = "STOP")
  
    coord_st_data <- getFlippedTissueCoordinates(data, 
                                                 as_data_frame=TRUE)
  }else if(inherits(data, "data.frame")){
    if(length(ident) != length(rownames(data)))
      print_msg("The 'ident' vector should heave same length as the number of spots.",
                msg_type = "STOP")
    coord_st_data <- data[,1:2]
    colnames(coord_st_data) <- c("x", "y")
  }
  
  coord_st_data$k <- ident
  
  path <- visium_hull(coord_st_data, 
                      size_x=size_x, 
                      size_y=size_y, 
                      delta=delta)
  

  if(hull_type == "walls"){
    p <- geom_segment(data=path$coord,
                      mapping=aes(x=x1,
                                  y=y1,
                                  xend=x2,
                                  yend=y2),
                      inherit.aes = FALSE,
                      color=color,
                      size=size)

  }else{

    layers_list <- list()
    
    u1 <- path$coord[, c("x1", "y1")]
    u2 <- path$coord[, c("x2", "y2")]
    colnames(u2) <- c("x","y")
    colnames(u1) <- c("x","y")
    s <- rbind(u1, u2)
    
    for(i in 1:length(path$neighborhoods)){
      p <- ggforce::geom_mark_hull(data=s[gsub("\\|\\|.*", "", rownames(s)) %in% path$neighborhoods[[i]],], 
                                   mapping=aes(x=x,
                                               y=y), 
                                   inherit.aes = FALSE,  
                                   expand = expand, 
                                   radius = radius, 
                                   concavity = concavity,
                                   size=size) 
      layers_list[[i]] <- p
    }
    
    p <- layers_list

  }

  return(p)            
}


