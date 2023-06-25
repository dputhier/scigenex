#' @title Given two spot groups, compute for each spot: (i) the distance of the
#' closest spot in the other class.
#' @param seurat_obj A Seurat object.
#' @param ident The classes of the spots (0/1). A named vector (name = spot name; value = classes). 
#' @param colours Colors for the spot classes.
#' @param colours_stratum A set of colors for the stratum.
#' @param pt_size The size of the points.
#' @param breaks Either a numeric vector of two or more unique cut points or a single number (greater than or equal to 2) giving the number of stratums in which spot are to be cut.
#' @param labels Labels for the levels of the stratum. By default, labels are constructed using "(a,b]" interval notation. If labels = FALSE, simple integer codes are returned instead of a factor.
#' @param diagnostic_plot Whether to produce a diagnostic diagram. Highly recommanded to visually inspect the results.
#' @export stratify_seurat
#' @importFrom ggplot2 aes element_blank element_rect facet_wrap geom_col geom_text scale_color_manual theme theme_minimal
#' @importFrom reshape2 melt
#' @importFrom scales hue_pal
#' @examples 
#' library(Seurat)
#' load_example_dataset("7870305/files/lymph_node_tiny_2")
#' Seurat::SpatialDimPlot(lymph_node_tiny_2)
#' identity <- Idents(lymph_node_tiny_2)
#' classes <- ifelse(Idents(lymph_node_tiny_2) == 7, 1, 0)
#' names(classes) <- names(identity)
#' h <- display_hull(lymph_node_tiny_2, ident=classes, color="black", delta=1, size_x=3.4, size_y=3)
#' plot_spatial(lymph_node_tiny_2, metadata = "seurat_clusters", pt_size = 6) + h
#' strat <- stratify_seurat(lymph_node_tiny_2, ident=classes)
#' strat[[1]]
#' strat[[2]] + coord_equal()
#' strat[[3]] + coord_equal()
#' strat[[4]] + coord_equal() 
stratify_seurat <- function(seurat_obj, 
                        ident=NULL,
                        colours=NULL,
                        colours_stratum=NULL,
                        pt_size = 4,
                        breaks=5,
                        labels = NULL,
                        diagnostic_plot=TRUE){
  
  if(!inherits(seurat_obj, "Seurat")){
     print_msg("Please provide a Seurat object.", msg_type="STOP")
  }
  
  ident_obj <- Idents(seurat_obj)
  
  if(length(ident) != length(ident_obj)){
    print_msg("The number of spots ('ident' argument) should be the same as in the Seurat object.", msg_type="STOP")
  }
  
  if(length(ident) != length(ident_obj)){
    print_msg("The number of spots ('ident' argument) should be the same as in the Seurat object.", msg_type="STOP")
  }
  
  if(is.null(names(ident))){
    print_msg("The 'ident' argument should be a named vector", msg_type="STOP")
  }else{
    
    ni <- names(ident)
    ident <- as.character(ident)
    names(ident) <- ni
    
    if(any(!names(ident) %in% colnames(seurat_obj))){
      print_msg("The 'ident' argument should be have names/spots defined in the Seurat object.", msg_type="STOP")
    }
    
    if(any(sort(unique(ident)) != c("0", "1"))){
      print_msg("Please use 0/1 to label classes.", msg_type="STOP")
    }
  }

  coord_spot <- getFlippedTissueCoordinates(seurat_obj, as_data_frame = TRUE)
  coord_spot <- coord_spot[names(ident), ]

  spot_0 <- names(ident[ident == "0"])
  spot_1 <- names(ident[ident == "1"])
  
  dist_to_point <- pracma::distmat(as.matrix(coord_spot[spot_0, ]), 
                                   as.matrix(coord_spot[spot_1, ]))
  
  value_min_dist_to_spot_0 <- apply(dist_to_point, 1, min)
  names(value_min_dist_to_spot_0) <- spot_0
  
  name_min_dist_to_spot_0 <- apply(dist_to_point, 
                                           1, 
                                           function(x,y) y[which(x==min(x))][1], 
                                           colnames(dist_to_point))
  names(name_min_dist_to_spot_0) <- spot_0
  
  value_min_dist_to_spot_1 <- apply(dist_to_point, 2, min)
  names(value_min_dist_to_spot_1) <- spot_1
  
  name_min_dist_to_spot_1 <- apply(dist_to_point, 
                                           2, 
                                           function(x,y) y[which(x==min(x))][1], 
                                           rownames(dist_to_point)) 
  names(name_min_dist_to_spot_1) <- spot_1
  
  coord_spot$closest_x  <- NA
  coord_spot$closest_y  <- NA
  coord_spot$closest_name  <- NA
  coord_spot$closest_dist  <- NA
  coord_spot$classes  <- NA
  
  coord_spot[spot_0, ]$closest_x <-  coord_spot[name_min_dist_to_spot_0, "x"]
  coord_spot[spot_1, ]$closest_x <-  coord_spot[name_min_dist_to_spot_1, "x"]

  coord_spot[spot_0, ]$closest_y <-  coord_spot[name_min_dist_to_spot_0, "y"]
  coord_spot[spot_1, ]$closest_y <-  coord_spot[name_min_dist_to_spot_1, "y"]
  
  coord_spot[spot_0, ]$closest_name <-  name_min_dist_to_spot_0[spot_0]
  coord_spot[spot_1, ]$closest_name <-  name_min_dist_to_spot_1[spot_1]
  
  coord_spot[spot_0, ]$closest_dist <-  value_min_dist_to_spot_0[spot_0]
  coord_spot[spot_1, ]$closest_dist <-  -value_min_dist_to_spot_1[spot_1]
  
  coord_spot[spot_0, ]$classes <- 0
  coord_spot[spot_1, ]$classes <- 1  
  
  seurat_obj$classes <-  as.factor(setNames(coord_spot$classes, 
                                            rownames(coord_spot)))

  if(is.null(colours)){
    colours <- discrete_palette(n=2, palette = "ggplot")
  }else{
    if(length(colours) != length(levels(seurat_obj$classes))){
      print_msg("Please provide the right number of colors.", msg_type = "STOP")
    }
  }

  p1 <- plot_spatial(seurat_obj, 
                     metadata =  "classes", 
                     colours = colours,
                     pt_size=pt_size) + 
        ggplot2::labs(fill="Classes") +
        ggplot2::geom_segment(data=coord_spot, 
                            inherit.aes = FALSE,
                            ggplot2::aes(x=x, 
                                         y=y, 
                                         xend=closest_x, 
                                         yend=closest_y, 
                                         color=classes), 
                            linewidth=0.3, show.legend = FALSE) 
    
  print_msg("Preparing diagnostic plot #1") 


  coord_spot$stratum  <- cut(coord_spot$closest_dist, breaks=breaks, labels = NULL)
  
  seurat_obj$stratum <-  as.factor(setNames(coord_spot$stratum, 
                                            rownames(coord_spot)))
  
  if(is.null(colours_stratum)){
    colours_stratum <- discrete_palette(n=length(levels(seurat_obj$stratum)), 
                                        palette = "ggplot")
  }else{
    if(length(colours_stratum) != length(levels(seurat_obj$stratum))){
      print_msg("Please provide the right number of colors for stratums.", msg_type = "STOP")
    }
  }
  
  p2 <- plot_spatial(seurat_obj, 
                     metadata =  "stratum", 
                     colours = colours,
                     pt_size=pt_size) +
        ggplot2::scale_fill_manual(values=unname(colours_stratum)) + 
        ggplot2::scale_colour_manual(values=unname(colours_stratum))
  
  print_msg("Preparing diagnostic plot #2") 
  


  p3 <- ggplot2::ggplot(data=coord_spot, 
                  ggplot2::aes(x=stratum, fill=stratum, col=stratum)) +
    ggplot2::geom_bar() +
    ggplot2::theme_bw() +
    ggplot2::coord_flip() + 
    ggplot2::scale_fill_manual(values=unname(colours_stratum)) + 
    ggplot2::scale_colour_manual(values=unname(colours_stratum))
  
  print_msg("Preparing diagnostic plot #3") 
  
  for(i in 2:ncol(coord_spot)){
    seurat_obj[[paste0("dist2border_", colnames(coord_spot)[i])]] <- setNames(coord_spot[ ,i], rownames(coord_spot))
  }
  
  return(list(seurat_obj, p1, p2, p3))
}

#' @title Display the results of stratify_seurat
#' @description
#' Display the results of stratify_seurat
#' @param seurat_obj A Seurat object.
#' @param gene_name The name of the gene to plot.
#' @param metadata Provide the name of a metadata that will be used instead of genes (i.e. from 
#' meta.data) slot of a seurat object.
#' @param intensity_slot The assay slot to use for the gene expression values.
#'        Must be one of "sct", "counts", or "data". Default is "sct".
#' @param colours_stratum A set of colors for the stratum.
#' @param polar Whter to use polar coordinates.
#' @export plot_stratum
#' @importFrom ggplot2 aes element_blank element_rect facet_wrap geom_col geom_text scale_color_manual theme theme_minimal
#' @importFrom reshape2 melt
#' @examples
#' library(Seurat)
#' load_example_dataset("7870305/files/lymph_node_tiny_2")
#' Seurat::SpatialDimPlot(lymph_node_tiny_2)
#' identity <- Idents(lymph_node_tiny_2)
#' classes <- ifelse(Idents(lymph_node_tiny_2) == 7, 1, 0)
#' names(classes) <- names(identity)
#' h <- display_hull(lymph_node_tiny_2, ident=classes, color="black", delta=1, size_x=3.4, size_y=3)
#' plot_spatial(lymph_node_tiny_2, metadata = "seurat_clusters", pt_size = 6) + h
#' strats <- stratify_seurat(lymph_node_tiny_2, ident=classes,
#'                           colours_stratum = rev(discrete_palette(5, "De1")))
#'  p <- plot_stratum(strats[[1]], gene_name="CCL19", polar = TRUE, 
#'                    colours_stratum = rev(discrete_palette(5, "De1")))
#' p + strats[[3]] + coord_equal()  + NoLegend()
plot_stratum <- function(seurat_obj=NULL,
                         gene_name=NULL,
                         metadata=NULL,
                         intensity_slot=c("data", "counts", "sct"),
                         colours_stratum=NULL,
                         polar=FALSE){
  
  intensity_slot <- match.arg(intensity_slot)
  
  if(is.null(colours_stratum)){
    colours_stratum <- discrete_palette(n=length(levels(seurat_obj$stratum)), 
                                        palette = "ggplot")
  }else{
    if(length(colours_stratum) != length(levels(seurat_obj$stratum))){
      print_msg("Please provide the right number of colors for stratums.", msg_type = "STOP")
    }
  }
  
  if(is.null(seurat_obj))
    print_msg("Please provide a seurat object...", msg_type = "STOP")
  
  if(is.null(gene_name) & is.null(metadata))
    print_msg("Please provide a value for gene_name or metadata...", msg_type = "STOP")
  
  if(!"dist2border_stratum" %in% names(seurat_obj@meta.data)){
    print_msg("Please use the stratify_seurat() function first.", msg_type = "STOP")
  }
  
  if(!is.null(metadata)){
    
    if(!metadata %in% colnames(seurat_obj@meta.data))
      print_msg("The metadata was not found in the object", msg_type = "STOP")
    intensities <- t(seurat_obj@meta.data[ , metadata, drop=FALSE])
    
  }else{
    if(intensity_slot=="sct"){
      slot_intensity <- seurat_obj@assays$SC
    }else if(intensity_slot=="counts"){
      slot_intensity <- seurat_obj@assays$Spatial@counts
    }else if(intensity_slot=="data"){
      slot_intensity <- seurat_obj@assays$Spatial@data
    }
    
    if(is.null(slot_intensity))
      print_msg("Slot is empty.", msg_type = "STOP")
    
    if(!gene_name %in% rownames(slot_intensity))
      print_msg("The gene_name was not found in the object", msg_type = "STOP")
    intensities <- slot_intensity[gene_name, , drop=FALSE]
    intensities <- as.matrix(intensities)
  }
  
  intensities <- reshape2::melt(intensities)
  colnames(intensities) <- c("gene", "spot", "value")
  intensities$stratum <- seurat_obj$dist2border_stratum[intensities$spot] 

  p <- ggplot2::ggplot(data=intensities, ggplot2::aes(x=stratum, y=value, fill=stratum)) +
    ggplot2::geom_col() +
    ggplot2::scale_fill_manual(values=unname(colours_stratum)) + 
    ggplot2::scale_colour_manual(values=unname(colours_stratum)) +
    ggplot2::theme_bw() 
  
  if(polar)
    p <- p + ggplot2::coord_polar(theta = "y")

  p
}


