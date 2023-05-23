set_verbosity(0)
load_example_dataset('7871581/files/pbmc3k_medium_clusters')

test_that("Check cluster_set_to_xls", {
  
  tp_dir <- file.path(tempdir(), create_rand_str())
  dir.create(tp_dir, showWarnings = F, recursive = TRUE)
  f_path <- file.path(tp_dir, "test.xls")
  cluster_set_to_xls(pbmc3k_medium_clusters, f_path)
  testthat::expect_true(file.exists(f_path))
  testthat::expect_error(cluster_set_to_xls(pbmc3k_medium_clusters, f_path))
 
  tp_dir <- file.path(tempdir(), create_rand_str())
  dir.create(tp_dir, showWarnings = F, recursive = TRUE)
  f_path <- file.path(tp_dir, "test.xls")
  cluster_set_to_xls(pbmc3k_medium_clusters, f_path, single_tab = TRUE)
  testthat::expect_true(file.exists(f_path))
  testthat::expect_error(cluster_set_to_xls(pbmc3k_medium_clusters, f_path)) 
})

