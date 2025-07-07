# Changelog

## v2.0.8

* Various fix for gm_report()
* Fix only_pos=TRUE in Spatial vignette.

## v2.0.7

* Various improvement including better support for SeuratFindAllMarkers() in gm_report()

## v2.0.6

* Improved vignette 'report.Rmd'.

## v2.0.5

* Refactored module.Rmd file called by gm_report().
* Added a vignette about gm_report().

## v2.0.4

* Fi gm_report() (AI prompt).

## v2.0.1-3

* Fi gm_report() (Added MF and CC as enrichment analyses).

## v2.0.0


* The clusterset_report() is deprecated. The gm_report() should be used instead.
* The @data slot in ClusterSet object is now a dcgMatrix (i.e. sparse matrix) !
* The cluster are no more computed by default. use compute_centers().
* The 'which_slot' argument of select_genes() has been replaced by 'layer'.
* Select_genes() function now has an 'assay' argument. 
* Suppressed some message in vignettes.
* Add a rename argument to filter_cluster_size(), filter_nb_supporting_cells(), filter_by_dot_prod(), filter_cluster_sd().
* In plot_profiles() argument color_cell_type has been renamed color_cluster.
* Fix 'Functions or methods with usage in Rd file  but not in code' using @noRd (instead of @keyword internal).

## v1.6.1

* The top_genes() has new argument 'distance_method' to allow overridding the default distance method (which can be
unknown when object is obtained using cluster_set_from_matrix()).
* The top_genes() has a new argument 'fast' to allow fast computation of Pearson-based distances.
* Added plot_multi_profiles().

## v1.6.0

* The FDR in `select_genes()` is now computed differently using BH. This may impact your results. 
* The FDR is no more a purcentage but a ratio. The fdr default value in `select_genes()` has been 
changed accordingly. 
* The 'cluster' argument has been deleted from top_genes(). The method has been refactored.
* The reorder_genes() has been refactored. Added argument 'decreasing'. 

## v1.5.4

* Add write_cname and file_suffix arguments to write_clust(). Changed default behaviour.
* Fix #150.

## v1.5.2

* Fix #147 and #145.
* Added cmp_to_a_list() function.
* Added coord_flip argument to plot_ggheatmap()
* Added top_by_grep()
* Added as_list argument to grep_clust().
* Added single_file to write_clust().
* Added to_lin to plot_profile.
* Deleted dot_prod from cluster_stats() and added CV.
* The row_sum is now normalized by cluster size in cluster_stats().
* Added top_by_intersect().

## v1.5.1

* Fix gene_clustering() and mcl_install(). Added 'mcl_dir' argument to `gene_clustering()` (in case automated installation in difficult, e.g., when running in a proprietary environnement).
* Fix `top_by_go()`. Added argument 'as_list'.
* Documentation has been splitted in two part (scRNA-seq and ST).

## v1.4.13

* Added support for Louvain and walktrap: call_walktrap_clusterset(), call_louvain_clusterset()
* construct_new_graph() was renamed do_closest_neighbor_graph()
* keep_dbf_graph() was renamed do_reciprocal_neighbor_graph()
* code refactoring (mainly comments).
* Deleted code relics (ident argument as list) from plot_ggheatmap().
* Added the subsample_by_ident() method for a ClusterSet object. handy to be used with plot_ggheatmap or plot_heatmap (to put emphasis on small populations).
* Added grep_clust() to search for genes in a ClusterSet object using a REGEXP.
* Added reload_pac() (for development).

## v1.4.10

* Fix NAMESPACE (xlxs dependency).
* Rebuilt doc in html format.

## v1.4.9

* Remove xlxs dependency in favor with WriteXLS dependency (to avoid rJava dependency).

## v1.4.8

* Fix an issue with temp dir.

## v1.4.7

* Change k_g argument to 's' in `gene_clustering()` function for clarity purpose.

* Added dependency to DT library (required by `cluster_set_report()`).

* Updated html doc.

* Fix an issue with pandoc in `clusterset_report()`.

* Improved importfrom in `clusterset_report()`.

## v1.4.6

* Prefix version with a v...

## v1.4.5

* keep_nn argument is deprecated (`gene_clustering()`).

## v1.4.4

* Fix #134

## v1.4.3

* Remove DynmicTreeCut dependency.

* Fix bug #137.

## v1.4.2

* Added color for cell cluster in `plot_ggheatmap()`.

* Updated doc.

## v1.4.1

* Minor bug fix.

* Updated doc

## v1.4.0

*   Update of the documentation
*   Added support for spatial transcriptomics: `getFlippedTissueCoordinates()`, `plot_spatial()`, `plot_spatial_panel()`.
*   Added import from matrix and Seurat object: `cluster_set_from_matrix()`, `cluster_set_from_seurat()`.
*   Added `cluster_stats()` and `plot_cluster_stats()` to perform statistics on clusters.
*   Added export: `cluster_set_to_xls()`.
*   Added `compare_genesets()` and `plot_cmp_genesets()` to compare gene lists with various statistics.
*   Added `plot_clust_enrichments()` to compare functional enrichment over clusters.
*   Added `plot_ggheatmap()` for statically display heatmaps.
*   Added `plot_markers_to_clusters()` to map markers onto clusters.
*   Added `top_by_go()` to select top genes by GO terms.

## v1.3.1

*   Update of the documentation
*   Improve get_genes function : gene names are ordered by gene clusters
*   Fix viz_enrich function : viz_enrich function is now working when using a specific ontology db as ontology parameter in enrich_go function


## v1.3.0

*   Add a clustering method using hierchical clustering and DynamicTreeCut R package
*   Add a cell_clusters slot in the ClusterSet object
*   Create a new documentation
*   Modification of the av_dot_prod_min filter : use the max of the median values of the dot product
*   Add get_genes and get_cells functions
*   Improve the computational speed of DBFMCL function
*   Replace the name of the slot "cluster" by "gene_patterns"
*   Implementation of new test (testthat)


## v1.2.21

*   Update the package used for enrichment analysis to use ClusterProfiler
*   Update visualization of enrichment analysis results : add dotplot and barplot provided by ClusterProfiler
*   Add more parameters to control the density plot of distances
*   Fix issues from Github

## v1.2.2

*   Replace the C code for DBF function by a R implementation to reduce the computational time.
*   Split the main R script into sub-script related to each function of scigenex package.

## v1.1.2

*   Create plot_heatmap function to generate an interactive heatmap of the gene clusters
*   Modify NAMESPACE and DESCRIPTION file and add import_package.R files to automatically update dependencies using roxygen2::roxygenise()
*   Create top_genes function to find the top n genes from each signature. This function is also usefull for the plot_heatmap function.
*   Add a slot in the ClusterSet object containing the filename of the output file including simulated distances and cutting threshold.
*   Modify the input of plot_dist function to a ClusterSet object.
*   Update enric_analysis function : run the functional enrichment analysis for each gene cluster.
*   Remove R implementation of MCL. The package used is no more maintained and the released version needs to be fixed.


## v1.0.2

*   Fix fprint_selected function in C code byt changing type of len varibale to size_t
*   Set default optional_output parameter of DBF and DBFMCL function to TRUE

## v1.0.1

### Bug Fixes

### API/CLI Changes

*   Scigenex now depends on R 4.0.0
*   No more reference to Biobase in the code.

### Code changes

### New Features


## v0.1.7

### Bug Fixes

*   Fix call to MCL.

### API/CLI Changes

### Code changes


### New Features


## v0.1.6

### Bug Fixes

### API/CLI Changes

*   Now provides a support for R MCL library (not recommended).


### Code changes


### New Features



## v0.1.4

### Bug Fixes


### API/CLI Changes

*   The *DBFMCLresult* object has been replaced by a more versatile object: *ClusterSet*.

*   The *plot_dbf()* method has been replaced by *plot_clust()*.

*   The *write_dbf()* method has been replaced by *write_clust()*.


### Code changes

*   Several code changes which were not clearly traced at the moment... Sry.

### New Features

