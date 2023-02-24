# Load seurat object
data(pbmc_small, package = "SeuratObject")
set_verbosity(0)
res <- find_gene_clusters(
  data = pbmc_small,
  distance_method = "pearson",
  inflation = 2,
  k = 20,
  row_sum = -Inf,
  highest = 0.3,
  min_nb_supporting_cell = 0,
  fdr = 1e-8
)


test_that("Checking if load_seurat stops when\
          seurat_obj is not a Seurat object", {
  expect_error(load_seurat(
    object = res,
    seurat_obj = "not a seurat object",
    dimplot_obj = DimPlot(pbmc_small)
  ))
})


test_that("Checking if load_seurat stops when\
          dimplot_obj is not a patchwork object", {
  expect_error(load_seurat(
    object = res,
    seurat_obj = pbmc_small,
    dimplot_obj = "not a patchwork object"
  ))
})


test_that("Checking if load_seurat stops when\
          object is not a ClusterSet object", {
  expect_error(load_seurat(
    object = "not a ClusterSet object",
    seurat_obj = pbmc_small,
    dimplot_obj = DimPlot(pbmc_small)
  ))
})


test_that("Checking cells metadata provided by load_seurat()", {
  res_cells_meta <- load_seurat(
    object = res,
    seurat_obj = pbmc_small,
    dimplot_obj = DimPlot(pbmc_small)
  )

  expect_equal(res_cells_meta@cells_metadata$cell_colors_from_seurat, c(
    "#F8766D", "#F8766D", "#F8766D", "#F8766D", "#F8766D", "#F8766D",
    "#F8766D", "#F8766D", "#F8766D", "#F8766D", "#619CFF", "#619CFF",
    "#619CFF", "#619CFF", "#619CFF", "#619CFF", "#619CFF", "#619CFF",
    "#619CFF", "#619CFF", "#00BA38", "#00BA38", "#00BA38", "#00BA38",
    "#00BA38", "#00BA38", "#00BA38", "#00BA38", "#00BA38", "#00BA38",
    "#F8766D", "#619CFF", "#F8766D", "#F8766D", "#F8766D", "#F8766D",
    "#F8766D", "#F8766D", "#F8766D", "#F8766D", "#F8766D", "#619CFF",
    "#F8766D", "#F8766D", "#619CFF", "#619CFF", "#F8766D", "#F8766D",
    "#F8766D", "#F8766D", "#00BA38", "#00BA38", "#619CFF", "#00BA38",
    "#00BA38", "#00BA38", "#00BA38", "#00BA38", "#00BA38", "#F8766D",
    "#619CFF", "#00BA38", "#619CFF", "#00BA38", "#619CFF", "#00BA38",
    "#00BA38", "#00BA38", "#00BA38", "#619CFF", "#F8766D", "#F8766D",
    "#F8766D", "#F8766D", "#F8766D", "#F8766D", "#F8766D", "#F8766D",
    "#F8766D", "#00BA38"
  ))

  expect_equal(sort(res_cells_meta@cells_metadata$cell_legend_from_seurat), c(
    rep(1, 36), rep(2, 25), rep(3, 19)
  ))

  expect_equal(res_cells_meta@cells_metadata$cell_order_from_seurat, c(
    1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 31, 33, 34, 35, 36, 37, 38, 39, 40, 41, 43,
    44, 47, 48, 49, 50, 60, 71, 72, 73, 74, 75, 76, 77, 78, 79, 21, 22, 23, 24,
    25, 26, 27, 28, 29, 30, 51, 52, 54, 55, 56, 57, 58, 59, 62, 64, 66, 67, 68,
    69, 80, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 32, 42, 45, 46, 53, 61, 63,
    65, 70
  ))

  expect_that(res_cells_meta, is_a("ClusterSet"))
})
