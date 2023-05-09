#---------------------------------------
# The Complex9 dataset
#---------------------------------------

#' The Complex9 dataset
#' A set of 2D points corresponding to 9 shapes.
#' 
#' Obtained from Salvador, S. and Chan, P., "Determining the Number of Clusters/Segments 
#' in Hierarchical clustering/Segmentation Algorithm”, ICTAI 2004,576-584."
#' Downloaded from:
#' https://www2.cs.uh.edu/~ml_kdd/restored/Complex&Diamond/complex&diamond.htm
#' 
#' Dataset information:
#'  ·        Number of Instances: 3031
#'  ·        Number of Attributes: 2 numeric (X,Y coordinate) and the class
#'  ·        class: total 9 classes(from 0 to 8)
#'  ·        Missing Attribute Values: None
#'  ·        Format:  x, y, (z) class-label
#'
#' \itemize{
#'   \item x: x coordinates
#'   \item y: y coordinates
#'   \item z: class-label
#' }
#'
#' @docType data
#' @keywords datasets
#' @name complex9
#' @usage data(complex9)
#' @format A data frame with 3031 rows and 3 variables
#' 
'complex9'

#---------------------------------------
# The complex9Noisy dataset
#---------------------------------------

#' The complex9Noisy dataset
#' A noisy version of the complex9 dataset. A set of 2D points corresponding 
#' to 9 shapes embedded in a noisy environment. A uniform noise was added to
#' the dataset.
#' 
#' Modified based on Salvador, S. and Chan, P., "Determining the Number of 
#' Clusters/Segments in Hierarchical clustering/Segmentation Algorithm”, 
#' ICTAI 2004,576-584." Downloaded from:
#' https://www2.cs.uh.edu/~ml_kdd/restored/Complex&Diamond/complex&diamond.htm
#' 
#' Dataset information:
#'  ·        Number of Instances: 3637 (vs 3031 in the original dataset)
#'  ·        Number of Attributes: 2 numeric (x,y coordinate) and the class
#'  ·        class: total 10 classes(from 0 to 9). Last one in added noise.
#'  ·        Missing Attribute Values: None
#'  ·        Format:  x, y, (z) class-label
#'
#' \itemize{
#'   \item x: x coordinates
#'   \item y: y coordinates
#'   \item z: class-label
#' }
#'
#' @docType data
#' @keywords datasets
#' @name complex9Noisy
#' @usage data(complex9Noisy)
#' @format A data frame with 3637 rows and 3 variables
#' 
'complex9Noisy'


#---------------------------------------
# The lymph_node_selected_spot dataset
#---------------------------------------

#' @title The lymph_node_selected_spot dataset
#' @description
#'  A vector with selected spots from 10X Visium Lymph node dataset.
#'  See https://tinyurl.com/y9hv7vj3
#'  and https://tinyurl.com/4bwh759j
#' 
#' 
#' Dataset information:
#'  ·        Number of spots: 601
#'
#'
#' @docType data
#' @keywords datasets
#' @name lymph_node_selected_spot
#' @usage data(lymph_node_selected_spot)
#' @format A vector containing 601 Visium spots.
#' 
#' @note
#' # This dataset was produced using loupe browser (10X genomics)
#' 
'lymph_node_selected_spot'

# -------------------------------------------------------------------------
# The lymph_node_selected_spot dataset (#2)
# -------------------------------------------------------------------------
#' @title The lymph_node_selected_spot dataset (#2)
#' @description
#'  A vector with selected spots from 10X Visium Lymph node dataset.
#'  See https://tinyurl.com/y9hv7vj3
#'  and https://tinyurl.com/4bwh759j
#' 
#' 
#' Dataset information:
#'  ·        Number of spots: 443
#'
#'
#' @docType data
#' @keywords datasets
#' @name lymph_node_selected_spot_2
#' @usage data(lymph_node_selected_spot)
#' @format A vector containing 601 Visium spots.
#' 
#' @note
#' # This dataset was produced using loupe browser (10X genomics)
#' 
'lymph_node_selected_spot_2'