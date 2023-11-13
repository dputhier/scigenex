#' @title Create a report from a ClusterSet and Seurat object.
#' @description Create a report from a ClusterSet and Seurat object.
#' @param clusterset_object The ClusterSet object.
#' @param seurat_object The Seurat object.
#' @param file_path A file path where to store the report (html extension).
#' @param force Whether to force erase output file.
#' @param report_title A title for the report.
#' @param report_author Names of report authors.
#' @param heatmap_colors Colors for the heatmap gradient.
#' @param heatmap_color_ident Colors for cell classes.
#' @param coord_flip Whether to flip spatial coordinates.
#' @param pt_size Size of the spots.
#' @param go_info Info about GO stats. Should be from the list: "ONTOLOGY", "ID", "Description", 
#' "GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue", "geneID", "Count"
#' @param verbosity The verbosity level (0-3, defaut to 1.)
#' @param pandoc_heap_size Fix the maximum memory used by pandoc. See https://pandoc.org/MANUAL.html.
#' @examples 
#' # Load datasets
#' # load_example_dataset("8028126/files/pbmc3k_medium_clusters_enr")
#' # load_example_dataset('7871581/files/pbmc3k_medium')
#' # cluster_set_report(clusterset_object = pbmc3k_medium_clusters_enr[1:2,],
#' #                    seurat_object=pbmc3k_medium)
#' # load_example_dataset("7870305/files/lymph_node_tiny_clusters_2")
#' # load_example_dataset("7870305/files/lymph_node_tiny_2")
#' # cluster_set_report(clusterset_object = lymph_node_tiny_clusters_2[1:2,],
#'                    seurat_object=lymph_node_tiny_2,
#'                    pt_size=6)
#' @export cluster_set_report
cluster_set_report <- function(clusterset_object=NULL,
                               seurat_object=NULL,
                               file_path = NULL,
                               force=FALSE,
                               report_title="Scigenex report",
                               report_author="Unknown",
                               heatmap_colors = colors_for_gradient("Ju1"),
                               heatmap_color_ident=NULL,
                               coord_flip=TRUE,
                               pt_size=2.75,
                               go_info=c("ID", "Description", "GeneRatio",  "pvalue", "qvalue", "Count"),
                               verbosity=1,
                               pandoc_heap_size ="512m") {

  check_format_cluster_set(clusterset_object)
  
  print_msg("Preparing report.")

  set_verbosity(verbosity)

  if(! all(go_info %in% c("ONTOLOGY", "ID", "Description", 
                          "GeneRatio", "BgRatio", "pvalue", 
                          "p.adjust", "qvalue", "geneID", "Count"))){
    print_msg("Please check go_info argument.")
  }
  
  if(!is.null(heatmap_color_ident)){
    heatmap_color_ident <- dput(heatmap_color_ident)
  }else{
    heatmap_color_ident <- "NULL"
  }
    

  # Dealing with file path
  
  if(!is.null(file_path)){
    
    file_dir_name <- dirname(file_path)
    file_basename <- basename(file_path)
    
    if(file.exists(file_path)){
      if(!force)
        print_msg("Output file already exists.", msg_type = "STOP")
      
      if(!grepl("\\.html$", file_basename))
        print_msg("Output file should have .html extension.", msg_type = "STOP")
      
    }else{
      
      try_dir <- try(dir.create(file_dir_name, recursive = TRUE, showWarnings = FALSE), 
                     silent = TRUE)
      if(class(try_dir) == "try-error")
        print_msg("Unable to create directory.", msg_type = "STOP")
    }
    
  }else{
    file_path <- paste0(tempfile(pattern="scigenex_report"), ".html")
    dir.create(dirname(file_path), showWarnings = FALSE, recursive = TRUE)
  }
  
  data_path <- gsub("\\.html", ".Rdata", file_path)
  
  print_msg('Dumping dataset')
  save(clusterset_object, 
       seurat_object,
       file=data_path)
  
  rmd_code <- 
    r"(---
title: {{ title }}
author: {{ report_author }}
date: '`r Sys.Date()`'
output:
  html_document:
    fig_caption: yes
    highlight: zenburn
    theme: cerulean
    toc: no
    toc_depth: 3
    pandoc_args: [
      '+RTS', '-K{{ pandoc_heap_size }}',
      '-RTS'
    ]
params:
  year: 
  region: Europe
  printcode: TRUE
  data: file.csv
---

<style type="text/css">
   .row {display: flex;}
   .column {flex: 50%;}
</style>


```{r, echo=FALSE, result='hide', message=FALSE, warning=FALSE }
set_verbosity(0)
load("{{data_path}}")
clusterset_object <- top_genes(clusterset_object)

# Computing addModuleScore

if("Spatial" %in% names(seurat_object@assays)){
  
  seurat_object <- Seurat::AddModuleScore(seurat_object, features = clusterset_object@gene_clusters, ctrl=20)
  
  for(i in 1:nclust(clusterset_object)){ # Normalizing module scores
    tmp <- seurat_object[[paste0("Cluster", i, sep="")]] 
    max_tmp <- max(tmp)
    min_tmp <- min(tmp)
    seurat_object[[paste0("Cluster", i, sep="")]]  <- (tmp[,1] - min(tmp))/(max_tmp - min_tmp)
  }
  
  for(i in 1:nclust(clusterset_object)){ # Normalizing module scores
    tmp <- seurat_object[[paste0("Cluster", i, sep="")]] 
    max_tmp <- max(tmp)
    min_tmp <- min(tmp)
    seurat_object[[paste0("Cluster_", i, sep="")]]  <- (tmp[,1] - min(tmp))/(max_tmp - min_tmp)
  }
} 
```
  
## Cluster summaries {.tabset .tabset-pills}
  
  
### All clusters
  

<div class = "row">
<div class = "column">

```{r, echo=FALSE, result='hide', message=FALSE, warning=FALSE }
DimPlot(seurat_object, reduction = "umap", label = TRUE)
```

</div>
<div class = "column">

```{r, echo=FALSE, result='hide', message=FALSE, warning=FALSE }
if("Spatial" %in% names(seurat_object@assays)){
  SpatialDimPlot(seurat_object, label = TRUE, label.size = 3, pt.size.factor = 1.4)
}
```

</div>
</div>

<div class = "row">
<div class = "column">

```{r, echo=FALSE, message=FALSE, warning=FALSE}
plot_ggheatmap(clusterset_object, 
               use_top_genes = TRUE,
               colors = {{heatmap_colors}},
               color_ident = {{heatmap_color_ident}},  
               ident = Seurat::Idents(seurat_object)) + ggplot2::ggtitle("All gene modules")

```

</div>
<div class = "column">

```{r, echo=FALSE, message=FALSE, warning=FALSE}

```

</div>
</div>


<div class = "row">
<div class = "column">

```{r, echo=FALSE, message=FALSE, warning=FALSE}
df_tmp <- data.frame(gene_cluster(clusterset_object))
DT::datatable(data.frame(Gene=rownames(df_tmp), 
                         Module=gene_cluster(clusterset_object)),
              rownames= FALSE)
```

</div>
<div class = "column">
</div>
</div>
)"

clust_code <- r"(

### Module {{clust_number}}

<div class = "row">
<div class = "column">

```{r, echo=FALSE, message=FALSE, warning=FALSE}
clust_sub <- clusterset_object[{{clust_number}},]
gene <- clust_sub@gene_clusters[[1]]
data_sub <- clust_sub@data
variance <- round(apply(data_sub, 1, var), 2)
sd <- round(apply(data_sub, 1, sd), 2)
sum_counts <- round(apply(data_sub, 1, sum), 2)
nb_pos_cell <- apply(data_sub > 0, 1, sum)  
DT::datatable(data.frame(Gene=gene,
                         Variance=variance,
                         Sd=sd,
                         Sum_counts=sum_counts,
                         Nb_pos_cell=nb_pos_cell), 
              rownames= FALSE)
```
  
</div>

<div class = "column">

```{r echo=FALSE, message=FALSE, warning=FALSE, result='hide'  }
if("Spatial" %in% names(seurat_object@assays)){
 plot_spatial(seurat_obj = seurat_object, 
             metadata = "Cluster_{{clust_number}}", 
             pt_size={{pt_size}}, coord_flip = {{coord_flip}})
}
```
</div>
</div>

<div class = "row">
<div class = "column">

```{r, echo=FALSE, message=FALSE, warning=FALSE}
plot_heatmap(clusterset_object[{{clust_number}}, ], 
               use_top_genes = TRUE,
               colors = {{heatmap_colors}},  
               cell_clusters = Seurat::Idents(seurat_object))
```
  
</div>

<div class = "column">

```{r echo=FALSE, message=FALSE, warning=FALSE, result='hide'  }

```
</div>
</div>


<div class = "row">
<div class = "column">

```{r, echo=FALSE, message=FALSE, warning=FALSE } 

```
</div>
<div class = "column">

```{r, echo=FALSE, message=FALSE, warning=FALSE }

```

</div>
</div>

<div class = "row">
<div class = "column">

```{r, echo=FALSE, message=FALSE, warning=FALSE }
if(length(clusterset_object@gene_cluster_annotations) > 0){
  DT::datatable(clusterset_object@gene_cluster_annotations[[{{clust_number}}]]@result[, {{go_info}}],
              rownames= FALSE)
}
```
</div>
<div class = "column">

```{r, echo=FALSE, message=FALSE, warning=FALSE }

```

</div>
</div>


)"
  
  rmd_code <- jinjar::render(rmd_code, 
                           title=report_title,
                           report_author=report_author,
                           data_path=data_path,
                           heatmap_colors=capture.output(dput(heatmap_colors))[[1]],
                           heatmap_color_ident=heatmap_color_ident,
                           pandoc_heap_size=pandoc_heap_size
  )


  
  capt_heatmap_colors <- capture.output(dput(heatmap_colors))[[1]]
  capt_go_info <- capture.output(dput(go_info))
  capt_go_info <- paste0(capt_go_info, collapse = "")
  coord_flip <- ifelse(coord_flip, "TRUE", "FALSE")
  

  
  
  for(clust_number in 1:nclust(clusterset_object)){

    glue_code <- jinjar::render(clust_code,
                                 heatmap_colors=capt_heatmap_colors,
                                 go_info=capt_go_info,
                                 pt_size=pt_size,
                                 coord_flip=coord_flip,
                                 heatmap_color_ident=heatmap_color_ident,
                                clust_number=clust_number
    )
    
    
    rmd_code <- paste0(rmd_code, 
                       glue_code,
                        collapse = "")
  }

  
  rmd_file <- gsub(pattern = ".html$", ".Rmd", file_path)
  
  print_msg(paste0("Writing report into ", rmd_file))
  writeLines(rmd_code, con = rmd_file)
  
  print_msg(paste0("Rmd file written at : ", rmd_file))
  print_msg("Rendering report.")
  print_msg(paste0("Creating report at : ", file_path))
  rmarkdown::render(rmd_file, quiet = TRUE)

}

