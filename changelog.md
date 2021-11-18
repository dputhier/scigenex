
# Changelog

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

