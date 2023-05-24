
# -------------------------------------------------------------------------
# plot_spatial function ---------------------------------------------------
# -------------------------------------------------------------------------

#' @title XY scatter plot of Visium data with hexagonal spots.
#' @description
#'  Given a Seurat Spatial object, this function creates a scatter plot of the 
#' spatial expression of a gene across spots, where the X and Y coordinates 
#' represent the spatial location of spots and the color represents the 
#' expression level of the gene. This function as been tested with visium 
#' data at the moment. Defaut shape is hexagon.
#'
#' @param seurat_obj A Seurat object containing spatial expression data.
#' @param gene_name The name of the gene to plot.
#' @param metadata Provide the name of a metadata that will be used instead of genes (i.e. from 
#' meta.data) slot of a seurat object.
#' @param intensity_slot The assay slot to use for the gene expression values.
#'        Must be one of "sct", "counts", or "data". Default is "sct".
#' @param title The title of the plot. Default is an empty string.
#' @param size_title The size of the title.
#' @param face_title Font face for the title. Possible values are “plain”, “italic”, “bold” and “bold.italic”.
#' @param legend Whether to display a legend for the color scale. Default is FALSE.
#' @param barwidth A numeric or a grid::unit() object specifying the width of the colourbar. 
#' Default value is legend.key.width or legend.key.size in theme() or theme.
#' @param barheight A numeric or a grid::unit() object specifying the height of the colourbar. 
#' Default value is legend.key.height or legend.key.size in theme() or theme.
#' @param axis Whether to display a axis for the color scale. Default is FALSE.
#' @param pt_size The size of the points in the plot. Default is 2.1.
#' @param pt_shape The shape of the points in the plot. Default is 16 (a circle).
#' @param pt_star A boolean. Whether to use ggstar shapes.
#' @param stroke The thickness of margin of points.
#' @param colours A vector of colors.
#' @param coord_flip Whether to flip coordinates.
#' @return A ggplot2 object containing the scatter plot.
#'
#' @importFrom ggplot2 ggplot geom_point scale_color_gradientn theme_void
#'              ggtitle element_text margin guide_colourbar coord_flip
#' @importFrom Seurat NoLegend NoAxes
#' @importFrom ggstar geom_star
#'
#' @examples
#' library(Seurat)
#' load_example_dataset("7870305/files/lymph_node_tiny_2")
#' plot_spatial(seurat_obj = lymph_node_tiny_2, gene_name = "CCL22", intensity_slot="data", pt_size=6)
#' plot_spatial(seurat_obj = lymph_node_tiny_2, metadata = "nCount_Spatial", pt_size=6)
#' hull <- display_hull(lymph_node_tiny_2, 
#'         ident=ifelse(Seurat::Idents(lymph_node_tiny_2) %in% 7, 1, 0),
#'         delta=1, size_x=3.4, size_y=3, color="black")
#' p <- plot_spatial(seurat_obj = lymph_node_tiny_2, gene_name = "VPREB3", intensity_slot="data", pt_size=6)
#' p + hull
#' @export
plot_spatial <- function(seurat_obj=NULL,
                         gene_name=NULL,
                         metadata=NULL,
                         intensity_slot=c("data", "counts", "sct"),
                         title="",
                         size_title=10,
                         face_title=c("plain", "italic", "bold", "bold.italic"),
                         legend=TRUE,
                         barwidth=1,
                         barheight=3,
                         axis=TRUE,
                         pt_size=3.6,
                         pt_shape=6,
                         pt_star=TRUE,
                         stroke=0,
                         colours=colors_for_gradient("Ju1"),
                         coord_flip=T){

  intensity_slot <- match.arg(intensity_slot)
  face_title <- match.arg(face_title)

  if(is.null(seurat_obj))
    print_msg("Please provide a seurat object...", msg_type = "STOP")
  
  if(is.null(gene_name) & is.null(metadata))
    print_msg("Please provide a value for gene_name or metadata...", msg_type = "STOP")
  
  print_msg("Getting x/y coordinates", msg_type = "DEBUG")
  
  if(coord_flip){
    print_msg("Using flipped coordinates.", msg_type = "DEBUG")
    xy_coord <- getFlippedTissueCoordinates(seurat_obj, 
                                            as_data_frame = TRUE)
  }else{
    print_msg("Coordinates won't be flipped.", msg_type = "DEBUG")
    xy_coord <- GetTissueCoordinates(seurat_obj)
    colnames(xy_coord) <- c("x", "y")
  }

  print_msg(paste0("Title is: ", title), msg_type = "DEBUG")
  print_msg("Extracting expression values", msg_type = "DEBUG")
  
  if(!is.null(metadata)){
    
    if(!metadata %in% colnames(seurat_obj@meta.data))
      print_msg("The metadata was not found in the object", msg_type = "STOP")
    intensities <- seurat_obj@meta.data[,metadata]
    
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
    intensities <- as.vector(slot_intensity[gene_name, ])
  }
  
  
  print_msg("Creating a ggplot diagram.", msg_type = "DEBUG")
  

  if(!is.factor(intensities)){
    
    print_msg("Feature is not a factor.", msg_type = "INFO")
    
    if(pt_star){
      the_scale <- scale_fill_gradientn(colours=colours, guide = guide_colourbar(barwidth=barwidth, barheight=barheight))
    }else{
      the_scale <- scale_color_gradientn(colours=colours, guide = guide_colourbar(barwidth=barwidth, barheight=barheight))
    }
    
  }else{
    
    print_msg("Feature is a factor.", msg_type = "INFO")
    
    if(length(colours) < length(levels(intensities))){
      print_msg("Not enough colors supplied. Creating a ggplot-like palette")
      colours <- discrete_palette(n=length(levels(intensities)), palette="ggplot")
    }
      
    if(pt_star){
      the_scale <- scale_fill_manual(values = setNames(colours[1:length(levels(intensities))],
                                                       as.character(levels(intensities))))
    }else{
      the_scale <- scale_color_manual(values = setNames(colours[1:length(levels(intensities))],
                                                        as.character(levels(intensities))))
    }
    
  }
  
  df <- cbind(xy_coord, intensities)
  colnames(df) <- c("x", "y", "intensity")
  
  
  if(pt_star){
    p <- ggplot(df, aes(x=x, y=y, fill=intensity)) + 
      ggstar::geom_star(starshape=pt_shape, 
                               starstroke=stroke,
                               size=pt_size) +
      the_scale
  } else{
    p <- ggplot(df, aes(x=x, y=y, color=intensity))+ 
      geom_point(size=pt_size, shape = pt_shape) +
      the_scale
  }
  
    p <- p + theme_void() +
              ggtitle(title) +
              theme(legend.title = element_text(size=6),
                    legend.text = element_text(size=6),
                    legend.margin = margin(0,0,0,0),
                    plot.title = element_text(face = face_title, 
                                              size=size_title,
                                              ))
  if(!legend)
    p <- p + NoLegend()
  
  if(!axis)
    p <- p + Seurat::NoAxes()
  
  return(p)
}


# -------------------------------------------------------------------------
# plot_spatial panel function  --------------------------------------------
# -------------------------------------------------------------------------


#' @title A planel of XY scatter plots of Visium data with hexagonal spots.
#'
#' @description
#'  This function creates a panel of scatter plots for the spatial expression of a list of genes across spots, 
#' where the X and Y coordinates represent the spatial location of spots and the color represents the 
#' expression level of the gene.
#'
#' @param seurat_obj A Seurat object containing spatial expression data.
#' @param genes A vector of gene names to plot.
#' @param metadata Provide a vector of metadata that will be used instead of genes (i.e. from 
#' meta.data) slot of a seurat object.
#' @param panel_names A vector of panel names to use for each gene plot.  
#' @param size_title The size of the titles.
#' @param face_title Font face for the title. Possible values are “plain”, “italic”, “bold” and “bold.italic”.
#' @param ncol_layout Number of columns to use for the panel layout. Default is the ceiling 
#' of the number of genes divided by 2.
#' @param intensity_slot The assay slot to use for the gene expression values.
#'        Must be one of "sct", "counts", or "data". Default is "sct".
#' @param title The title of the plot. Default is an empty string.
#' @param legend Whether to display a legend for the color scale. Default is FALSE.
#' @param barwidth A numeric or a grid::unit() object specifying the width of the colourbar. 
#' Default value is legend.key.width or legend.key.size in theme() or theme.
#' @param barheight A numeric or a grid::unit() object specifying the height of the colourbar. 
#' Default value is legend.key.height or legend.key.size in theme() or theme.
#' @param axis Whether to display a axis for the color scale. Default is FALSE.
#' @param pt_size The size of the points in the plot. Default is 2.1.
#' @param pt_shape The shape of the points in the plot. Default is 16 (a circle).
#' @param pt_star A boolean. Whether to use ggstar shapes.
#' @param guides  A string specifying how guides should be treated in the layout. See patchwork::plot_layout().
#' @param stroke The thickness of margin of points.
#' @param coord_flip Whether to flip coordinates.
#' @param colours A vector of colors.
#' @importFrom ggplot2 ggplot theme_void
#' @importFrom patchwork plot_layout
#'
#' @return A ggplot2 object containing the panel of scatter plots.
#'
#' @examples
#' library(Seurat)
#' load_example_dataset("7870305/files/lymph_node_tiny_2")
#' load_example_dataset("7870305/files/lymph_node_tiny_clusters_2")
#' lymph_node_tiny_2 <- Seurat::AddModuleScore(lymph_node_tiny_2, features = lymph_node_tiny_clusters_2@gene_clusters, nbin = 15)
#' for(i in 1:nclust(lymph_node_tiny_clusters_2)){ # Normalizing module scores
#'     tmp <- lymph_node_tiny_2[[paste0("Cluster", i, sep="")]] 
#'     max_tmp <- max(tmp)
#'     min_tmp <- min(tmp)
#'     lymph_node_tiny_2[[paste0("Cluster", i, sep="")]]  <- (tmp[,1] - min(tmp))/(max_tmp - min_tmp)
#' }
#' plot_spatial_panel(lymph_node_tiny_2, metadata=paste0("Cluster", 1:4), ncol_layout=2,
#'                    guides='collect', pt_size=2.2, coord_flip=TRUE)
#'                    
#' plot_spatial_panel(lymph_node_tiny_2, gene=c('VPREB3', 'IGHG1', 'PRDX4', 
#'                                              'LTB', 'CCL20', 'LYVE1', 
#'                                              'IL7R', 'RGS9', 'MAPT'), 
#'                    ncol_layout=3,
#'                    pt_size=2, coord_flip=TRUE, 
#'                    panel_names=LETTERS[1:9], 
#'                    size_title = 10)
#' @export
plot_spatial_panel <- function(seurat_obj=NULL,
                               genes=NULL,
                               metadata=NULL,
                               intensity_slot=c("data", "counts", "sct"),
                               title="",
                               size_title=10,
                               face_title=c("plain", "italic", "bold", "bold.italic"),
                               barwidth=1,
                               barheight=3,
                               axis=TRUE,
                               panel_names=NULL,
                               ncol_layout=NULL,
                               legend=TRUE,
                               guides=NULL,
                               pt_size=3.6,
                               pt_shape=6,
                               pt_star=TRUE,
                               stroke=0,
                               coord_flip=TRUE,
                               colours=colors_for_gradient("Ju1")
){

  intensity_slot <- match.arg(intensity_slot)
  face_title <- match.arg(face_title)
  
  if(is.null(panel_names)){
    if(is.null(genes)){
      panel_names <- LETTERS[1:length(metadata)]
    }else{
      panel_names <- LETTERS[1:length(genes)]
    }
  }
  
  if(is.null(ncol_layout)){
    if(is.null(genes)){
      ncol_layout <- ceiling(length(metadata)/2)
    }else{
      ncol_layout <- ceiling(length(genes)/2)
    }
  }
  
  
  if(is.null(seurat_obj))
    print_msg("Please provide a seurat object...", msg_type = "STOP")
  
  if(is.null(genes) & is.null(metadata))
    print_msg("Please provide a list of genes or metadata...", msg_type = "STOP")
  
  
  if(is.null(genes) & is.null(metadata))
    print_msg("Please provide a list of genes or metadata...", msg_type = "STOP") 
  
  print_msg(paste0("Panel names : ",  panel_names), msg_type = "DEBUG")
  
  plot_panels <- NULL
  
  if(is.null(metadata)){
    if(length(panel_names) != length(genes))
      print_msg("panel_names and genes should have same length.", msg_type = "STOP") 
    
    for(i in 1:length(genes)){
      
      gene_curr <- genes[i]
      panel_curr <- panel_names[i]
      
      print_msg(paste0("Creating diagram for gene: ",  gene_curr), msg_type = "DEBUG")
      
      plot_cur <- plot_spatial(seurat_obj=seurat_obj,
                                    gene_name=gene_curr,
                                    intensity_slot=intensity_slot,
                                    title=panel_curr,
                                    legend=legend,
                                    pt_size=pt_size,
                                    pt_shape=pt_shape,
                                    pt_star=pt_star,
                                    size_title=size_title,
                                    stroke=stroke,
                                    colours=colours,
                                    coord_flip=coord_flip)
      
      if(is.null(plot_panels)){
        plot_panels <- plot_cur
      }else{
        plot_panels <- plot_panels + plot_cur
      }
    }
  }else{
    for(i in 1:length(metadata)){
      
      metadata_curr <- metadata[i]
      panel_curr <- panel_names[i]
      
      print_msg(paste0("Creating diagram for metadata: ",  metadata_curr), msg_type = "DEBUG")
      
      plot_cur <- plot_spatial(seurat_obj=seurat_obj,
                                    metadata=metadata_curr,
                                    title=panel_curr,
                                    legend=legend,
                                    pt_size=pt_size,
                                    pt_shape=pt_shape,
                                    pt_star=pt_star,
                                    size_title=size_title,
                                    stroke=stroke,
                                    colours=colours,
                                    coord_flip=coord_flip)
      
      if(is.null(plot_panels)){
        plot_panels <- plot_cur
      }else{
        plot_panels <- plot_panels + plot_cur
      }
    }
  }
  
  print_msg("Preparing diagram layout", msg_type = "DEBUG")
  
  plot_panels + patchwork::plot_layout(ncol=ncol_layout, guides = guides)
}

