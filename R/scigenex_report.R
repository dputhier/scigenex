#===========================================================
# Create a report (v2) for scigenex 
#===========================================================
#' @title Create a report from a ClusterSet experiment.
#' @description
#' Create a report from a ClusterSet experiment and Seurat object. If Seurat contains spatial information set 'spatial' argument to TRUE and selection the 
#' corresponding assay and layer. The function call Gemini IA to try to guess cell type and functions associated to each cluster.
#' @param cluster_set A clusterSet object.
#' @param seurat_object Seurat object.
#' @param seurat_assay Which assay should be used in the Seurat object.
#' @param seurat_layer Which layer should be used in the Seurat object.
#' @param spatial_experiment Whether there are some spatial information in the ClusterSet object.
#' @param report_title A title for the report.
#' @param report_subtitle A subtitle for the report.
#' @param report_author A character string corresponding to one or severe al authors.
#' @param report_date A date.
#' @param out_dir A directory where to store the bookdown output.
#' @param experimenters Name of the experimenter.
#' @param experimenter_labs Laboratory of the experimenter.
#' @param api_key An API key for Gemini to perform cell type / cell function analysis. Otherwise deactivate section 'module_cell_annot_IA'.
#' @param smp_species The species of the sample (free text). E.g 'Homo sapiens'.
#' @param smp_stage Development stage of the sample (free text). E.g. 'adult'.
#' @param smp_organ The sample organ (free text). E.g. 'tonsil'.
#' @param smp_region A region in the organ. E.g. 'whole' (which will merge as 'whole tonsil').
#' @param rmd_dir A path where to find the templates for creating the book.
#' @param rm_tmpdir Whether to delete temporary directory.
#' @param quiet Whether to run bookdown::render_book() quietly.
#' @examples
#' ## TODO:  'org.Hs.eg.db' / enrich_go "Hsapiens" / library(enrichplot) / GOSemSim::godata / patchwork /xaringanExtra / "/Users/puthier/Documents/git/project_dev/scigenex"
#' library(scigenex)
#' library(Seurat)
#' library(clusterProfiler)
#' library(enrichplot)
#' load_example_dataset('7871581/files/pbmc3k_medium_clusters')
#' load_example_dataset('7871581/files/pbmc3k_medium')
#' scigenex_report(pbmc3k_medium_clusters, 
#'                 pbmc3k_medium, 
#'                 smp_species="Homo_sapiens", 
#'                 smp_region="total", 
#'                 smp_organ="blood", 
#'                 smp_stage="adult", rmd_dir="/Users/puthier/Documents/git/project_dev/scigenex/inst/rmarkdown", 
#'                 api_key="AIzaSyDj2dA0w4LoXi6LhXgf62vEOkuOUnpONbY")
#' @export scigenex_report
#' @importFrom fs path_home
#' @importFrom Seurat AddModuleScore
#' @importFrom bookdown render_book
scigenex_report <- function(cluster_set = NULL,
                            seurat_object=NULL,
                            seurat_assay=NULL,
                            seurat_layer=NULL,
                            spatial_experiment=FALSE,
                            report_title = "Spatial transcriptomics report",
                            report_subtitle = "An example experiment",
                            report_author = "Undefined",
                            report_date = format(Sys.time(), '%d %B %Y'),
                            out_dir = file.path(fs::path_home(), "scigenex_book"),
                            experimenters = NULL,
                            experimenter_labs=NULL,
                            api_key=NULL,
                            smp_species=NULL,
                            smp_stage=NULL,
                            smp_organ=NULL,
                            smp_region=NULL,
                            rmd_dir = file.path(system.file(package = "scigenex"), 
                                                                      "rmarkdown"),
                            add_module_score_params=list(pool = NULL,
                                                         nbin = 16,
                                                         ctrl = 100,
                                                         k = FALSE,
                                                         assay = NULL,
                                                         name = "MOD_",
                                                         seed = 1,
                                                         search = FALSE,
                                                         slot = "data"),
                            rm_tmpdir = TRUE,
                            section=c("exp_info", 
                                      "module_heatmap",
                                      "module_iheatmap",
                                      "module_umap",
                                      "module_violin",
                                      "module_cell_annot_IA",
                                      "module_term_network",
                                      "module_term_barplot_1",
                                      "module_term_barplot_2",
                                      "module_term_network_circ",
                                      "term_table"),
                            quiet=FALSE) {
  
  verb_level <- get_verbosity()
  
  if("module_cell_annot_IA" %in% section){
    if(is.null(smp_species) | is.null(smp_stage) | is.null(smp_organ) | is.null(smp_region))
      print_msg("Please provide information about species, stage, organ and region for IA-based annotation.")
  }

  if("module_cell_annot_IA" %in% section){
    if(is.null(api_key))
      print_msg("Please provide a Gemini API key.", msg_type='STOP')
    
    print_msg("Setting Gemini API key.", msg_type='INFO')
    gemini.R::setAPI(api_key)
  }
    
  check_format_cluster_set(cluster_set)
  
  if(nclust(cluster_set) < 1 ) 
    print_msg("No cluster found in ClusterSet.", msg_type="STOP")
  
  print_msg("Checking seurat object.", msg_type="DEBUG") 
  
  if(!is.null(seurat_object)){

    if (!inherits(seurat_object, "Seurat"))
      print_msg("Please provide a valid Seurat object.", msg="STOP") 
    
    if(length(Idents(seurat_object)) == 0 )
      print_msg("Seurat object needs to contain cll identities (check with Idents()).", msg="STOP")
    
    add_module_score_params$object <- seurat_object
    add_module_score_params$features <- cluster_set@gene_clusters
    add_module_score_params$name <- "MOD_"
    seurat_object <- do.call(Seurat::AddModuleScore, add_module_score_params)
  }

  tmp_dir <- tempdir(check = FALSE)
  tmp_dir <- paste0(tmp_dir, format(Sys.time(), "%a-%b-%e-%H-%M-%S-%Y"))
  dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)
  
  print_msg("Created temporary directory:", msg_type = "DEBUG")
  print_msg(tmp_dir)
  
  print_msg("Copying file in temporary directory.", msg_type = "DEBUG")
  
  out <- file.copy(rmd_dir, tmp_dir, recursive = TRUE)
  
  tmp_dir <- file.path(tmp_dir, "rmarkdown")

  # Compute top genes
  cluster_set <- top_genes(cluster_set)
  
  if(any(c("module_term_network",
         "module_term_barplot_1",
         "module_term_barplot_2",
         "module_term_network_circ",
         "term_table") %in% section)){
    print_msg("Computing semantic similarity for cnet_plot.", msg_type="INFO") 
    suppressMessages(sem_sim <- GOSemSim::godata('org.Hs.eg.db', ont="BP"))
  }

  
  module_rmd <- file.path(tmp_dir, "module.rmd")
  
  print_msg("Preparing parameters for the report.")
  
  print_msg("Looping through parameterised Reports.")
  
  n <- 1
  all_knited_files <- vector()
  
  for(n in 1:nclust(cluster_set)){
    
    print_msg(paste0("Preparing rmd files for modules", n), msg_type = "DEBUG")
    fig_path <- paste0("figure_", n, "_")
    
    cur_rmd <- gsub(".rmd$" , paste0("_", sprintf("%04d", n), ".rmd"), module_rmd)
    
    print_msg(paste0("file : ", cur_rmd), msg_type = "DEBUG")
    
    code_rmd  <- readLines(module_rmd)
    code_rmd  <- gsub("FIG_PATH", fig_path, x = code_rmd)
    code_rmd  <- gsub("MODULE_NUMBER", n, x = code_rmd)
    writeLines(code_rmd, con = cur_rmd)
    
    print_msg(paste0("Preparation of rmd files for objects", n, "finished."), msg_type = "DEBUG")
    
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
  code_yml  <- readLines(file.path(tmp_dir, "_bookdown.yml"), warn = FALSE)
  
  options(knitr.duplicate.label = "allow")
  
  bookdown::render_book(tmp_dir, quiet=quiet)
  set_verbosity(verb_level)
  
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  file.copy(file.path(tmp_dir, "_book"), out_dir, recursive=TRUE)
  
  print_msg(paste0("Results have been copied to:", out_dir))
  
  if(rm_tmpdir)
    unlink(tmp_dir, recursive=TRUE)
  
  set_verbosity(verb_level)
}
