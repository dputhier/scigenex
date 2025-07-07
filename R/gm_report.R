#====================================================================
# Create a report (v2) for a scigenex::clusterSet and Seurat object 
#====================================================================
#' @title Create a report from a ClusterSet and Seurat objects
#' @description
#' Create a report from a ClusterSet experiment and Seurat object. If is_spatial_exp argument is set to TRUE Spatial diagram will be added. 
#' The function can call Gemini IA to try to guess cell type and functions associated to each gene module.
#' NB: When dealing with large dataset (large number of cells/spots), it is advisable not activate interactive heatmap by deleting "module_iheatmap" from the 'section' argument.
#' @param cluster_set A clusterSet object.
#' @param seurat_object Seurat object.
#' @param seurat_assay Which assay should be used in the Seurat object.
#' @param seurat_layer Which layer should be used in the Seurat object.
#' @param is_spatial_exp Whether there are some spatial information in the ClusterSet object.
#' @param annotation_src The sources of functional annotation (currently a vector taken from "BP", "CC", "MF").
#' @param report_title A title for the report.
#' @param report_subtitle A subtitle for the report.
#' @param report_author A character string corresponding to one or several authors.
#' @param report_date A date.
#' @param out_dir A directory where to store the bookdown output.
#' @param experimenters A data.frame providing information about experimenters.
#' @param workflow_params A data.frame indicating some workflow parameters.
#' @param bioc_org_db A gene annotation database as provided by Bioconductor (e.g. org.Hs.eg.db, org.Mm.eg.db, org.Rn.eg.db, o
#' rg.Dm.eg.db, org.Dr.eg.db, org.Sc.sgd.db, org.Ce.eg.db...). If NULL, module analysis related to functional annotation are canceled.
#' @param api_key An API key for Gemini to perform cell type / cell function analysis. If NULL, the 'module_cell_annot_IA' section will be skipped.
#' @param smp_species The species of the sample (free text). E.g 'Homo sapiens'.
#' @param smp_stage Development stage of the sample (free text). E.g. 'adult'.
#' @param smp_organ The sample organ (free text). E.g. 'tonsil'.
#' @param smp_region A region in the organ. E.g. 'whole' (which will merge as 'whole tonsil').
#' @param rmd_dir A path where to find the templates for creating the book.
#' @param add_module_score_params Some parameters for Seurat::add_module_score() function.
#' @param plot_profiles_params Some parameters for plot_profiles() function
#' @param plot_multi_profiles_params Some parameters for plot_profiles_multi() function.
#' @param FeaturePlot_params Some parameters for Seurat::FeaturePlot() function.
#' @param SpatialFeaturePlot_params Some parameters for Seurat::SpatialFeaturePlot() function.
#' @param SpatialDimPlot_params Some parameters for Seurat::SpatialDimPlot() function.
#' @param plot_ggheatmap_params Some parameters for plot_ggheatmap() function.
#' @param subsample_by_ident_params The number of cell to take per cell type when subsampling the data for plotting interactive heatmap.
#' Reduce this number if you have a large dataset (e.g. > 10000 cells/spots). Otherwise the dataset will be too large to be handled by the web browser.
#' @param plot_heatmap_params Some parameters for plot_heatmap() function.
#' @param cnetplot_params Some parameters for enrichplot:::cnetplot.enrichResult() function.
#' @param rm_tmpdir Whether to delete temporary directory.
#' @param section Which section to activate/deactivate.
#' @param quiet Whether to run bookdown::render_book() quietly.
#' @return No return value. A report is generated and written to the specified output directory.
#' @examples
#' library(scigenex)
#' library(Seurat)
#' set_verbosity(3)
#' load_example_dataset('7871581/files/pbmc3k_medium_clusters')
#' load_example_dataset('7871581/files/pbmc3k_medium')
#' gm_report(pbmc3k_medium_clusters[1:2,], 
#'                 pbmc3k_medium, 
#'                 smp_species="Homo_sapiens", 
#'                 annotation_src="CC",
#'                 smp_region="total", 
#'                 smp_organ="blood", 
#'                 smp_stage="adult", 
#'                 bioc_org_db="org.Hs.eg.db",
#'                 api_key=NULL)
#' set_verbosity(3)
#' load_example_dataset('7870305/files/lymph_node_tiny_clusters_2')
#' load_example_dataset('7870305/files/lymph_node_tiny_2')
#' gm_report(lymph_node_tiny_clusters_2[1:2,], 
#'                 lymph_node_tiny_2, 
#'                 smp_species="Homo sapiens", 
#'                 smp_region="total", 
#'                 smp_organ="lymph node", 
#'                 smp_stage="adult", 
#'                 annotation_src="CC",
#'                 bioc_org_db="org.Hs.eg.db",
#'                 api_key=NULL,
#'                 is_spatial_exp=TRUE,
#'                 SpatialFeaturePlot_params=list(pt.size.factor = 3000),
#'                 SpatialDimPlot_params=list(pt.size.factor = 3000)) # Object was created with an older seurat version
#' set_verbosity(3)
#' markers <- Seurat::FindAllMarkers(lymph_node_tiny_2, only.pos = TRUE)
#' cs <- cluster_set_from_seurat(lymph_node_tiny_2, markers, p_val_adj=0.001, assay="Spatial")
#' gm_report(cs[1:2,], 
#'                 lymph_node_tiny_2, 
#'                 smp_species="Homo sapiens", 
#'                 smp_region="total", 
#'                 smp_organ="lymph node", 
#'                 smp_stage="adult", 
#'                 annotation_src="CC",
#'                 is_spatial_exp=TRUE,
#'                 bioc_org_db="org.Hs.eg.db",
#'                 api_key=NULL) # Object was created with an older seurat version
#' @importFrom fs path_home
#' @importFrom bookdown render_book
#' @importFrom xaringanExtra use_panelset
#' @importFrom xaringanExtra style_panelset
#' @importFrom xaringanExtra style_panelset_tabs
#' @importFrom BiocManager install
#' @importFrom knitr opts_chunk
#' @importFrom knitr asis_output
#' @importFrom Seurat VlnPlot
#' @importFrom gemini.R setAPI
#' @importFrom gemini.R gemini_chat
#' @importFrom AnnotationDbi select
#' @importFrom dplyr mutate
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_text
#' @importFrom clusterProfiler enrichGO
#' @importFrom enrichplot pairwise_termsim
#' @importFrom enrichplot dotplot
#' @importFrom enrichplot cnetplot
#' @importFrom Seurat NoLegend
#' @importFrom Seurat SpatialFeaturePlot
#' @importFrom Seurat FeaturePlot
#' @importFrom Seurat AddModuleScore
#' @importFrom dplyr "%>%"
#' @importFrom pander pander
#' @importFrom GOSemSim godata
#' @importFrom htmlwidgets JS
#' @export
gm_report <- function(cluster_set = NULL,
                            seurat_object=NULL,
                            seurat_assay=NULL,
                            seurat_layer=NULL,
                            is_spatial_exp=FALSE,
                            annotation_src=c("BP", "CC", "MF"),
                            report_title = "Gene module report",
                            report_subtitle = "An example experiment",
                            report_author = "Undefined",
                            report_date = format(Sys.time(), '%d %B %Y'),
                            out_dir = file.path(fs::path_home(), "gm_book"),
                            experimenters=data.frame(),
                            workflow_params=data.frame(),
                            bioc_org_db=NULL,
                            api_key=NULL,
                            smp_species=NULL,
                            smp_stage=NULL,
                            smp_organ=NULL,
                            smp_region=NULL,
                            rmd_dir = file.path(system.file(package = "scigenex"), 
                                                "rmarkdown"),
                            add_module_score_params=list(nbin = 16,
                                                         ctrl = 100,
                                                         name = "MOD_",
                                                         slot = "data"),
                            plot_profiles_params=list(to_lin=TRUE, averaged = TRUE),
                            plot_multi_profiles_params=list(legend_name="Gene\nModule"),
                            FeaturePlot_params=list(cols=RColorBrewer::brewer.pal(3, "BuPu")),
                            SpatialFeaturePlot_params=list(pt.size.factor = 1.7),
                            SpatialDimPlot_params=list(pt.size.factor = 1.7),
                            plot_ggheatmap_params=list(use_top_genes=FALSE, 
                                                       hide_gene_name=TRUE,
                                                       xlab = "Cells/Spots",
                                                       ylab="Genes"),
                            subsample_by_ident_params=list(nbcell=200),
                            plot_heatmap_params=list(link="complete", 
                                                     use_top_genes=FALSE,
                                                     interactive=TRUE,
                                                     line_size_horizontal=1,
                                                     label_size=7),
                            cnetplot_params=list(color.params = list(edge = TRUE), 
                                                 showCategory = 6, 
                                                 cex.params = list(category_label = 0.6, 
                                                                   gene_label = 0.7)),
                            rm_tmpdir = TRUE,
                            section=c("exp_info",
                                      "exp_metadata",
                                      "exp_experimenters",
                                      "exp_sample",
                                      "exp_params", 
                                      "exp_stats",
                                      "exp_genes",
                                      "exp_heatmap",
                                      "exp_dimplot",
                                      "exp_metrics",
                                      "exp_pop",
                                      "exp_spatial_dist",
                                      "exp_spatial_dimplot",
                                      "exp_mean_1",
                                      "exp_mean_2",
                                      "module_spatial",
                                      "module_heatmap",
                                      "module_iheatmap",
                                      "module_umap",
                                      "module_violin",
                                      "module_genes",
                                      "module_cell_annot_IA",
                                      "module_term_network",
                                      "module_term_barplot_1",
                                      "module_term_barplot_2",
                                      "module_term_network_circ",
                                      "term_table"),
                            quiet=FALSE) {
  
  if(!all(annotation_src %in% c("BP", "CC", "MF")))
    print_msg("Unknow annotation source.", msg_type="STOP") 
    
  verb_level <- get_verbosity()
  
  if("module_cell_annot_IA" %in% section){
    if(is.null(smp_species) | is.null(smp_stage) | is.null(smp_organ) | is.null(smp_region))
      print_msg("Please provide information about species, stage, organ and region for IA-based annotation.")
  }
  
  print_msg("Checking seurat object.", msg_type="DEBUG") 
  
  if (!inherits(seurat_object, "Seurat"))
    print_msg("Please provide a valid Seurat object.", msg_type="STOP") 
  
  check_format_cluster_set(cluster_set)
  
  
  if(nclust(cluster_set) < 1 ) 
    print_msg("No cluster found in ClusterSet.", msg_type="STOP")
  
  if(is.null(bioc_org_db)){
    print_msg("No annotation database found for organism...", msg_type = "WARNING")
    print_msg("Canceling functional annotation", msg_type = "WARNING")
    section <- setdiff(section, c("module_term_network",
                                  "module_term_barplot_1",
                                  "module_term_barplot_2",
                                  "module_term_network_circ",
                                  "term_table"))
  }else{
    if(!require(bioc_org_db, character.only = TRUE, quietly=TRUE)){
      print_msg("The annotation library (see 'bioc_org_db') was not found. Please install it.", msg_type = "STOP")
    }else{
      library(bioc_org_db, character.only = TRUE)
    }
  }

  if(!is_spatial_exp){
    print_msg("This is not a ST experiment...", msg_type = "INFO")
    print_msg("Canceling ST reporting module.", msg_type = "INFO")
    section <- setdiff(section, c("exp_spatial_dist",
                                  "exp_spatial_dimplot",
                                  "module_spatial"))
  }
  
  if(is.null(api_key)){
    print_msg("No Gemini key provided...", msg_type = "INFO")
    print_msg("Canceling IA-based cell type annotation.", msg_type = "INFO")
    section <- setdiff(section, c("exp_spatial",
                                  "module_cell_annot_IA"))
  }
  
  if(length(Idents(seurat_object)) == 0 )
    print_msg("Seurat object needs to contain all identities (check with Idents()).", msg="STOP")
  
  print_msg("Computing module scores...", msg_type = "INFO")
  add_module_score_used_params <- as.list(formals(Seurat:::AddModuleScore.Seurat))
  
  for(i in names(add_module_score_params)){
    add_module_score_used_params[[i]] <- add_module_score_params[[i]]
  }
  
  add_module_score_used_params$object <- seurat_object
  add_module_score_used_params$features <- lapply(cluster_set@gene_clusters, gsub, pattern = "~[0-9]+$", replacement = "_")
  add_module_score_used_params$name <- "MOD_"
  seurat_object <- do.call(Seurat::AddModuleScore, add_module_score_used_params)
  
  tmp_dir <- tempdir(check = FALSE)
  tmp_dir <- paste0(tmp_dir, gsub(" ", "", format(Sys.time(), "%a-%b-%H-%M-%S-%Y")))
  
  dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)
  
  print_msg("Created temporary directory:", msg_type = "DEBUG")
  print_msg(tmp_dir)
  
  print_msg("Copying file in temporary directory.", msg_type = "DEBUG")
  
  out <- file.copy(rmd_dir, tmp_dir, recursive = TRUE)
  
  tmp_dir <- file.path(tmp_dir, "rmarkdown")
  
  for(tpfile in list.files(tmp_dir)){
    print_msg(paste0("The temporary dir contains file: ", tpfile))
  }
   
  print_msg("Computing top genes.", msg_type = "DEBUG")
  cluster_set <- top_genes(cluster_set)
  
  print_msg("Computing centers.", msg_type = "DEBUG")
  cluster_set <- compute_centers(cluster_set)
  

  sem_sim <- list()
  
  if(any(c("module_term_network",
             "module_term_barplot_1",
             "module_term_barplot_2",
             "module_term_network_circ",
             "term_table") %in% section)){
      for(ann in annotation_src){
        print_msg(paste0("Computing semantic similarity (", ann ,") for cnet_plot."), msg_type="INFO") 
        suppressMessages(sem_sim[[ann]] <- GOSemSim::godata(annoDb=eval(bioc_org_db),
                                                   ont=ann))
      }
  }
  
  
  module_rmd <- file.path(tmp_dir, "module.Rmd")
  
  print_msg("Preparing parameters for the report.")
  
  print_msg("Looping through parameterised Reports.")
  
  n <- 1
  
  all_knited_files <- vector()
  
  for(n in 1:nclust(cluster_set)){
    
    print_msg(paste0("Preparing Rmd files for modules ", n), msg_type = "DEBUG")
    fig_path <- paste0("figure_", n, "_")
    
    cur_rmd <- gsub(".Rmd$" , paste0("_", sprintf("%04d", n), ".Rmd"), module_rmd)
    
    print_msg(paste0("file : ", cur_rmd), msg_type = "DEBUG")
    
    code_rmd  <- readLines(module_rmd)
    code_rmd  <- gsub("FIG_PATH", fig_path, x = code_rmd)
    code_rmd  <- gsub("MODULE_NUMBER", n, x = code_rmd)
    
    if(!is.null(api_key)){
      print_msg("Setting Gemini API key.", msg_type='INFO')
      code_rmd  <- gsub("APIKEY", api_key, x = code_rmd)
    }
    
    writeLines(code_rmd, con = cur_rmd)
    
    print_msg(paste0("Preparation of Rmd files for objects ", n, " finished."), msg_type = "DEBUG")
    
    all_knited_files[n] <- basename(cur_rmd)
    
  }
  
  
  print_msg("Deleting sample.Rmd.")
  
  unlink(module_rmd)
  
  print_msg("preparing _bookdown.yml")
  
  all_knited_files_merged <- paste(sapply(all_knited_files, shQuote), collapse = ", ")
  
  code_yml  <- readLines(file.path(tmp_dir, "_bookdown.yml"), warn = FALSE)
  
  code_yml  <- gsub(
    "rmd_files:",
    paste0(
      "rmd_files: [ 'index.Rmd', ",
      all_knited_files_merged,
      ", 'footer.Rmd']"
    ),
    x = code_yml
  )
  
  writeLines(code_yml, con = file.path(tmp_dir, "_bookdown.yml"))
  
  options(knitr.duplicate.label = "allow")
  
  bookdown::render_book(tmp_dir, quiet=quiet)
  set_verbosity(verb_level)
  
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  for(tpfile in list.files(file.path(tmp_dir, "_book")))
    file.copy(file.path(tmp_dir, "_book", tpfile), out_dir, recursive=TRUE)
  
  dir.create(file.path(out_dir, "rmarkdown"), recursive = TRUE, showWarnings = FALSE)  
  
  rmd_input <- list.files(file.path(tmp_dir))
  rmd_input <-   rmd_input[rmd_input != "_book"]              
  
  for(tpfile in rmd_input)
    file.copy(file.path(tmp_dir, tpfile), file.path(out_dir, "rmarkdown"), recursive=TRUE)
  
  print_msg(paste0("Results have been copied to: ", out_dir))
  
  if(rm_tmpdir)
    unlink(tmp_dir, recursive=TRUE)
  
  set_verbosity(verb_level)
}
