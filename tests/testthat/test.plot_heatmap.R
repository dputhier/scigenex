# Create matrix containing 3 signatures
set_verbosity(0)

m <- create_4_rnd_clust()

res_scigenex <- find_gene_clusters(
  data = m,
  distance_method = "pearson",
  inflation = 2,
  k = 75,
  row_sum = -Inf,
  highest = 0.3,
  min_nb_supporting_cell = 0,
  fdr = 1e-8
)

test_that("Checking plot_heatmap()", {
  set.seed(123)
  
  htmp <- plot_heatmap(
    object = res_scigenex,
    cell_clusters = NULL,
    use_top_genes = FALSE,
    show_dendro = TRUE,
    line_size_vertical = 2,
    line_size_horizontal = 2,
    name = "Expression heatmap",
    show_legend = TRUE,
    xlab = "Samples",
    ylab = "Genes"
  )

  htmp_matrix <- htmp@plots@listData$`Exp. level`@data
  htmp_matrix_wo_na <- htmp_matrix[!is.na(htmp_matrix)]

  # Check heatmap format
  expect_that(htmp, is_a("IheatmapHorizontal"))

  # Check dim of heatmap matrix
  expect_equal(ncol(htmp_matrix), 20)
  expect_equal(nrow(htmp_matrix), 349)
  expect_equal(dim(htmp_matrix), c(349, 20))

  # Check NA rows added by line_size_horizontal
  na_rows <- which(apply(htmp_matrix, 1, function(x) all(is.na(x))))
  expect_equal(length(na_rows), 6)
  expect_equal(as.vector(na_rows), c(65, 66, 143, 144, 244, 245))

  # Check NA rows added by line_size_vertical
  na_cols <- which(apply(htmp_matrix, 2, function(x) all(is.na(x))))
  expect_equal(length(na_cols), 0)
  expect_equal(as.vector(na_cols), as.integer())

  # Check some measurements of heatmap matrix
  expect_equal(round(sum(htmp_matrix_wo_na), 2), 218.69)
  expect_equal(round(mean(htmp_matrix_wo_na), 2), 0.03)
  expect_equal(round(median(htmp_matrix_wo_na), 2), 0.12)
  expect_equal(round(sd(htmp_matrix_wo_na), 2), 0.87)

  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["0%"]), 2), -1)
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["25%"]), 2), -1)
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["50%"]), 2), 0.12)
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["75%"]), 2), 1)
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["100%"]), 2), 1)

  # Checking order of columns
  expect_equal(colnames(htmp_matrix), c("sample14", "sample11", "sample12", "sample13", "sample15", 
                                        "sample20", "sample17", "sample18", "sample16", "sample19", "sample1", 
                                        "sample5", "sample10", "sample4", "sample2", "sample3", "sample7", 
                                        "sample8", "sample6", "sample9"))

  # Checking order of rows
  expect_equal(rownames(htmp_matrix), c(
    "gene89", "gene20", "gene22", "gene99",
    "gene96", "gene85", "gene95", "gene58", "gene14", "gene35",
    "gene90", "gene11", "gene68", "gene46", "gene28", "gene71",
    "gene87", "gene17", "gene31", "gene25", "gene63", "gene10",
    "gene83", "gene82", "gene94", "gene42", "gene93", "gene79",
    "gene34", "gene74", "gene56", "gene1", "gene4", "gene13",
    "gene64", "gene66", "gene51", "gene52", "gene9", "gene86",
    "gene57", "gene2", "gene7", "gene39", "gene29", "gene36",
    "gene73", "gene78", "gene61", "gene19", "gene88", "gene45",
    "gene84", "gene12", "gene16", "gene59", "gene80", "gene67",
    "gene55", "gene70", "gene32", "gene37", "gene76", "gene27",
    " ", " ", "gene395", "gene2911", "gene304", "gene360",
    "gene920", "gene400", "gene384", "gene310", "gene323", "gene383",
    "gene344", "gene3609", "gene392", "gene391", "gene328", "gene317",
    "gene373", "gene337", "gene367", "gene345", "gene371", "gene342",
    "gene372", "gene355", "gene339", "gene370", "gene396", "gene351",
    "gene385", "gene301", "gene341", "gene303", "gene332", "gene311",
    "gene394", "gene358", "gene336", "gene378", "gene340", "gene365",
    "gene335", "gene362", "gene366", "gene397", "gene353", "gene319",
    "gene389", "gene327", "gene350", "gene368", "gene331", "gene306",
    "gene313", "gene349", "gene324", "gene329", "gene322", "gene379",
    "gene326", "gene382", "gene320", "gene330", "gene334", "gene363",
    "gene387", "gene398", "gene315", "gene318", "gene369", "gene364",
    "gene386", "gene316", "gene314", "gene312", "gene390", "gene393",
    " ", " ", "gene139", "gene102", "gene169", "gene123",
    "gene156", "gene125", "gene170", "gene110", "gene196", "gene131",
    "gene172", "gene149", "gene152", "gene120", "gene185", "gene164",
    "gene195", "gene134", "gene182", "gene115", "gene118", "gene130",
    "gene154", "gene181", "gene138", "gene179", "gene136", "gene151",
    "gene113", "gene174", "gene148", "gene199", "gene200", "gene167",
    "gene183", "gene119", "gene140", "gene198", "gene142", "gene128",
    "gene171", "gene141", "gene124", "gene175", "gene178", "gene158",
    "gene162", "gene116", "gene121", "gene163", "gene194", "gene105",
    "gene103", "gene177", "gene191", "gene197", "gene173", "gene157",
    "gene127", "gene114", "gene109", "gene104", "gene101", "gene188",
    "gene129", "gene107", "gene147", "gene160", "gene168", "gene137",
    "gene193", "gene176", "gene161", "gene122", "gene184", "gene145",
    "gene126", "gene143", "gene112", "gene166", "gene133", "gene159",
    "gene108", "gene153", "gene111", "gene192", "gene106", "gene132",
    "gene144", "gene155", "gene180", "gene186", "gene165", "gene135",
    "gene190", "gene117", "gene150", "gene146", "gene189", " ",
    " ", "gene2157", "gene2036", "gene847", "gene276", "gene231",
    "gene1931", "gene252", "gene274", "gene264", "gene265", "gene210",
    "gene218", "gene228", "gene3991", "gene209", "gene213", "gene261",
    "gene280", "gene243", "gene211", "gene296", "gene201", "gene224",
    "gene260", "gene271", "gene295", "gene221", "gene205", "gene208",
    "gene291", "gene202", "gene234", "gene220", "gene227", "gene257",
    "gene238", "gene288", "gene268", "gene292", "gene270", "gene246",
    "gene254", "gene279", "gene226", "gene281", "gene297", "gene244",
    "gene242", "gene247", "gene284", "gene294", "gene237", "gene269",
    "gene253", "gene215", "gene232", "gene216", "gene240", "gene207",
    "gene203", "gene258", "gene255", "gene239", "gene262", "gene233",
    "gene287", "gene219", "gene259", "gene212", "gene250", "gene235",
    "gene263", "gene300", "gene283", "gene229", "gene293", "gene217",
    "gene266", "gene223", "gene299", "gene290", "gene222", "gene267",
    "gene236", "gene298", "gene251", "gene206", "gene241", "gene277",
    "gene248", "gene214", "gene275", "gene282", "gene278", "gene286",
    "gene289", "gene225", "gene256", "gene230", "gene249", "gene272",
    "gene245", "gene285", "gene273"
  ))

  # Checking title of the heatmap
  htmp_title <- htmp@annotations@listData$col_title2@data
  expect_equal(htmp_title, "Expression heatmap")

  # Checking x lab of the heatmap
  htmp_xlab <- htmp@annotations@listData$col_title@data
  expect_equal(htmp_xlab, "Samples")

  # Checking y lab of the heatmap
  htmp_ylab <- htmp@annotations@listData$row_title@data
  expect_equal(htmp_ylab, "Genes")

  # Checking colorbar name of the heatmap
  htmp_colors_name <- htmp@colorbars@listData$`Exp. level`@title
  expect_equal(htmp_colors_name, "Exp. level")

  # Checking if colorbar is shown of the heatmap
  htmp_show_colorbar <- htmp@plots@listData$`Exp. level`@show_colorbar
  expect_true(htmp_show_colorbar)

  # Checking if dendrogram is showed
  htmp_show_dendro <- names(htmp@shapes@listData)
  expect_equal(htmp_show_dendro, "col_dendro")
})


# ==============================================================================
# ==============================================================================


test_that("Checking plot_heatmap() using top genes", {
  set.seed(123)
  
  res_scigenex_top <- top_genes(res_scigenex, top = 20)
  htmp <- plot_heatmap(
    object = res_scigenex_top,
    cell_clusters = NULL,
    use_top_genes = TRUE, # Use top genes
    show_dendro = FALSE, # Set to FALSE
    line_size_vertical = 2,
    line_size_horizontal = 2,
    name = NULL, # Set to NULL
    show_legend = FALSE, # Set to FALSE
    xlab = NULL, # Set to NULL
    ylab = NULL
  ) # Set to NULL

  htmp_matrix <- htmp@plots@listData$`Exp. level`@data
  htmp_matrix_wo_na <- htmp_matrix[!is.na(htmp_matrix)]

  # Check heatmpa format
  expect_that(htmp, is_a("IheatmapHorizontal"))

  # Check dim of heatmap matrix
  expect_equal(ncol(htmp_matrix), 20)
  expect_equal(nrow(htmp_matrix), 86)
  expect_equal(dim(htmp_matrix), c(86, 20))

  # Check NA rows added by line_size_horizontal
  na_rows <- which(apply(htmp_matrix, 1, function(x) all(is.na(x))))
  expect_equal(length(na_rows), 6)
  expect_equal(as.vector(na_rows), c(21, 22, 43, 44, 65, 66))

  # Check NA rows added by line_size_horizontal
  na_cols <- which(apply(htmp_matrix, 2, function(x) all(is.na(x))))
  expect_equal(length(na_cols), 0)
  expect_equal(as.vector(na_cols), as.integer())

  # Check some measurements of heatmap matrix
  expect_equal(round(sum(htmp_matrix_wo_na), 2), 64.5)
  expect_equal(round(mean(htmp_matrix_wo_na), 2), 0.04)
  expect_equal(round(median(htmp_matrix_wo_na), 2), 0.24)
  expect_equal(round(sd(htmp_matrix_wo_na), 2), 0.89)

  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["0%"]), 2), -1)
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["25%"]), 2), -1)
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["50%"]), 2), 0.24)
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["75%"]), 2), 1)
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["100%"]), 2), 1)

  # Checking order of columns
  expect_equal(colnames(htmp_matrix), c("sample14", "sample11", "sample12", "sample13", "sample15", 
                                        "sample20", "sample17", "sample18", "sample16", "sample19", "sample1", 
                                        "sample5", "sample10", "sample4", "sample2", "sample3", "sample7", 
                                        "sample8", "sample6", "sample9"))

  # Checking order of rows
  expect_equal(rownames(htmp_matrix), c(
    "gene7", "gene19", "gene57", "gene36",
    "gene78", "gene39", "gene67", "gene88", "gene16", "gene84",
    "gene29", "gene37", "gene27", "gene70", "gene59", "gene12",
    "gene80", "gene32", "gene55", "gene76", " ", " ",
    "gene306", "gene329", "gene397", "gene320", "gene382", "gene369",
    "gene379", "gene398", "gene387", "gene315", "gene364", "gene390",
    "gene334", "gene363", "gene326", "gene314", "gene386", "gene316",
    "gene393", "gene312", " ", " ", "gene111", "gene166",
    "gene133", "gene153", "gene126", "gene108", "gene132", "gene180",
    "gene106", "gene144", "gene186", "gene192", "gene135", "gene165",
    "gene150", "gene117", "gene155", "gene190", "gene146", "gene189",
    " ", " ", "gene206", "gene275", "gene277", "gene289",
    "gene222", "gene241", "gene214", "gene223", "gene236", "gene282",
    "gene225", "gene256", "gene278", "gene249", "gene286", "gene230",
    "gene245", "gene272", "gene273", "gene285"
  ))

  # Checking title of the heatmap
  htmp_slots <- names(htmp@annotations@listData)
  expect_false("col_title2" %in% htmp_slots)

  # Checking x lab of the heatmap
  expect_false("col_title" %in% htmp_slots)

  # Checking y lab of the heatmap
  expect_false("row_title" %in% htmp_slots)

  # Checking colorbar name of the heatmap
  htmp_colors_name <- htmp@colorbars@listData$`Exp. level`@title
  expect_equal(htmp_colors_name, "Exp. level")

  # Checking if colorbar is showed
  htmp_show_colorbar <- htmp@plots@listData$`Exp. level`@show_colorbar
  expect_false(htmp_show_colorbar)

  # Checking if dendrogram is showed
  htmp_show_dendro <- names(htmp@shapes@listData)
  expect_equal(htmp_show_dendro, NULL)
})


# ==============================================================================
# ==============================================================================


test_that("Checking plot_heatmap() using cell clusters", {
  cell_clusters <- c(rep(1, 10), rep(2, 10))
  names(cell_clusters) <- colnames(res_scigenex@data)
  htmp <- plot_heatmap(
    object = res_scigenex,
    cell_clusters = cell_clusters,
    use_top_genes = FALSE,
    show_dendro = TRUE,
    line_size_vertical = 2,
    line_size_horizontal = 2,
    name = "Expression heatmap",
    show_legend = TRUE,
    xlab = "Samples",
    ylab = "Genes"
  )

  htmp_matrix <- htmp@plots@listData$`Exp. level`@data
  htmp_matrix_wo_na <- htmp_matrix[!is.na(htmp_matrix)]

  # Check heatmpa format
  expect_that(htmp, is_a("IheatmapHorizontal"))

  # Check dim of heatmap matrix
  expect_equal(ncol(htmp_matrix), 22)
  expect_equal(nrow(htmp_matrix), 349)
  expect_equal(dim(htmp_matrix), c(349, 22))

  # Check NA rows added by line_size_horizontal
  na_rows <- which(apply(htmp_matrix, 1, function(x) all(is.na(x))))
  expect_equal(length(na_rows), 6)
  expect_equal(as.vector(na_rows), c(65, 66, 143, 144, 244, 245))

  # Check NA rows added by line_size_horizontal
  na_cols <- which(apply(htmp_matrix, 2, function(x) all(is.na(x))))
  expect_equal(length(na_cols), 2)
  expect_equal(as.vector(na_cols), c(11, 12))

  # Check some measurements of heatmap matrix
  expect_equal(round(sum(htmp_matrix_wo_na), 2), 218.69)
  expect_equal(round(mean(htmp_matrix_wo_na), 2), 0.03)
  expect_equal(round(median(htmp_matrix_wo_na), 2), 0.12)
  expect_equal(round(sd(htmp_matrix_wo_na), 2), 0.87)

  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["0%"]), 2), -1)
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["25%"]), 2), -1)
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["50%"]), 2), 0.12)
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["75%"]), 2), 1)
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["100%"]), 2), 1)

  # Checking order of columns
  expect_equal(colnames(htmp_matrix), c("sample1", "sample2", "sample3", "sample4", "sample5", "sample6", 
                                        "sample7", "sample8", "sample9", "sample10", " ", "  ", "sample11", 
                                        "sample12", "sample13", "sample14", "sample15", "sample16", "sample17", 
                                        "sample18", "sample19", "sample20"))

  # Checking order of rows
  expect_equal(rownames(htmp_matrix), c(
    "gene89", "gene20", "gene22", "gene99",
    "gene96", "gene85", "gene95", "gene58", "gene14", "gene35",
    "gene90", "gene11", "gene68", "gene46", "gene28", "gene71",
    "gene87", "gene17", "gene31", "gene25", "gene63", "gene10",
    "gene83", "gene82", "gene94", "gene42", "gene93", "gene79",
    "gene34", "gene74", "gene56", "gene1", "gene4", "gene13",
    "gene64", "gene66", "gene51", "gene52", "gene9", "gene86",
    "gene57", "gene2", "gene7", "gene39", "gene29", "gene36",
    "gene73", "gene78", "gene61", "gene19", "gene88", "gene45",
    "gene84", "gene12", "gene16", "gene59", "gene80", "gene67",
    "gene55", "gene70", "gene32", "gene37", "gene76", "gene27",
    " ", " ", "gene395", "gene2911", "gene304", "gene360",
    "gene920", "gene400", "gene384", "gene310", "gene323", "gene383",
    "gene344", "gene3609", "gene392", "gene391", "gene328", "gene317",
    "gene373", "gene337", "gene367", "gene345", "gene371", "gene342",
    "gene372", "gene355", "gene339", "gene370", "gene396", "gene351",
    "gene385", "gene301", "gene341", "gene303", "gene332", "gene311",
    "gene394", "gene358", "gene336", "gene378", "gene340", "gene365",
    "gene335", "gene362", "gene366", "gene397", "gene353", "gene319",
    "gene389", "gene327", "gene350", "gene368", "gene331", "gene306",
    "gene313", "gene349", "gene324", "gene329", "gene322", "gene379",
    "gene326", "gene382", "gene320", "gene330", "gene334", "gene363",
    "gene387", "gene398", "gene315", "gene318", "gene369", "gene364",
    "gene386", "gene316", "gene314", "gene312", "gene390", "gene393",
    " ", " ", "gene139", "gene102", "gene169", "gene123",
    "gene156", "gene125", "gene170", "gene110", "gene196", "gene131",
    "gene172", "gene149", "gene152", "gene120", "gene185", "gene164",
    "gene195", "gene134", "gene182", "gene115", "gene118", "gene130",
    "gene154", "gene181", "gene138", "gene179", "gene136", "gene151",
    "gene113", "gene174", "gene148", "gene199", "gene200", "gene167",
    "gene183", "gene119", "gene140", "gene198", "gene142", "gene128",
    "gene171", "gene141", "gene124", "gene175", "gene178", "gene158",
    "gene162", "gene116", "gene121", "gene163", "gene194", "gene105",
    "gene103", "gene177", "gene191", "gene197", "gene173", "gene157",
    "gene127", "gene114", "gene109", "gene104", "gene101", "gene188",
    "gene129", "gene107", "gene147", "gene160", "gene168", "gene137",
    "gene193", "gene176", "gene161", "gene122", "gene184", "gene145",
    "gene126", "gene143", "gene112", "gene166", "gene133", "gene159",
    "gene108", "gene153", "gene111", "gene192", "gene106", "gene132",
    "gene144", "gene155", "gene180", "gene186", "gene165", "gene135",
    "gene190", "gene117", "gene150", "gene146", "gene189", " ",
    " ", "gene2157", "gene2036", "gene847", "gene276", "gene231",
    "gene1931", "gene252", "gene274", "gene264", "gene265", "gene210",
    "gene218", "gene228", "gene3991", "gene209", "gene213", "gene261",
    "gene280", "gene243", "gene211", "gene296", "gene201", "gene224",
    "gene260", "gene271", "gene295", "gene221", "gene205", "gene208",
    "gene291", "gene202", "gene234", "gene220", "gene227", "gene257",
    "gene238", "gene288", "gene268", "gene292", "gene270", "gene246",
    "gene254", "gene279", "gene226", "gene281", "gene297", "gene244",
    "gene242", "gene247", "gene284", "gene294", "gene237", "gene269",
    "gene253", "gene215", "gene232", "gene216", "gene240", "gene207",
    "gene203", "gene258", "gene255", "gene239", "gene262", "gene233",
    "gene287", "gene219", "gene259", "gene212", "gene250", "gene235",
    "gene263", "gene300", "gene283", "gene229", "gene293", "gene217",
    "gene266", "gene223", "gene299", "gene290", "gene222", "gene267",
    "gene236", "gene298", "gene251", "gene206", "gene241", "gene277",
    "gene248", "gene214", "gene275", "gene282", "gene278", "gene286",
    "gene289", "gene225", "gene256", "gene230", "gene249", "gene272",
    "gene245", "gene285", "gene273"
  ))

  # Checking title of the heatmap
  htmp_title <- htmp@annotations@listData$col_title2@data
  expect_equal(htmp_title, "Expression heatmap")

  # Checking x lab of the heatmap
  htmp_xlab <- htmp@annotations@listData$col_title@data
  expect_equal(htmp_xlab, "Samples")

  # Checking y lab of the heatmap
  htmp_ylab <- htmp@annotations@listData$row_title@data
  expect_equal(htmp_ylab, "Genes")

  # Checking colorbar name of the heatmap
  htmp_colors_name <- htmp@colorbars@listData$`Exp. level`@title
  expect_equal(htmp_colors_name, "Exp. level")

  # Checking if colorbar is shown of the heatmap
  htmp_show_colorbar <- htmp@plots@listData$`Exp. level`@show_colorbar
  expect_true(htmp_show_colorbar)

  # Checking if dendrogram is showed
  htmp_show_dendro <- names(htmp@shapes@listData)
  expect_equal(htmp_show_dendro, NULL)
})


# ==============================================================================
# ==============================================================================


test_that("Checking plot_heatmap() using top genes and cell clusters", {
  cell_clusters <- c(rep(1, 10), rep(2, 10))
  names(cell_clusters) <- colnames(res_scigenex@data)

  res_scigenex_top <- top_genes(res_scigenex, top = 20)

  htmp <- plot_heatmap(
    object = res_scigenex_top,
    cell_clusters = cell_clusters,
    use_top_genes = TRUE,
    show_dendro = TRUE,
    line_size_vertical = 2,
    line_size_horizontal = 2,
    name = "Expression heatmap",
    show_legend = TRUE,
    xlab = "Samples",
    ylab = "Genes"
  )

  htmp_matrix <- htmp@plots@listData$`Exp. level`@data
  htmp_matrix_wo_na <- htmp_matrix[!is.na(htmp_matrix)]

  # Check heatmpa format
  expect_that(htmp, is_a("IheatmapHorizontal"))

  # Check dim of heatmap matrix
  expect_equal(ncol(htmp_matrix), 22)
  expect_equal(nrow(htmp_matrix), 86)
  expect_equal(dim(htmp_matrix), c(86, 22))

  # Check NA rows added by line_size_horizontal
  na_rows <- which(apply(htmp_matrix, 1, function(x) all(is.na(x))))
  expect_equal(length(na_rows), 6)
  expect_equal(as.vector(na_rows), c(21, 22, 43, 44, 65, 66))

  # Check NA rows added by line_size_horizontal
  na_cols <- which(apply(htmp_matrix, 2, function(x) all(is.na(x))))
  expect_equal(length(na_cols), 2)
  expect_equal(as.vector(na_cols), c(11, 12))

  # Check some measurements of heatmap matrix
  expect_equal(round(sum(htmp_matrix_wo_na), 2), 64.5)
  expect_equal(round(mean(htmp_matrix_wo_na), 2), 0.04)
  expect_equal(round(median(htmp_matrix_wo_na), 2), 0.24)
  expect_equal(round(sd(htmp_matrix_wo_na), 2), 0.89)

  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["0%"]), 2), -1)
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["25%"]), 2), -1)
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["50%"]), 2), 0.24)
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["75%"]), 2), 1)
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["100%"]), 2), 1)

  # Checking order of columns
  expect_equal(colnames(htmp_matrix), c("sample1", "sample2", "sample3", "sample4", "sample5", "sample6", 
                                        "sample7", "sample8", "sample9", "sample10", " ", " ", "sample11", 
                                        "sample12", "sample13", "sample14", "sample15", "sample16", "sample17", 
                                        "sample18", "sample19", "sample20"))

  # Checking order of rows
  expect_equal(rownames(htmp_matrix), c(
    "gene7", "gene19", "gene57", "gene36",
    "gene78", "gene39", "gene67", "gene88", "gene16", "gene84",
    "gene29", "gene37", "gene27", "gene70", "gene59", "gene12",
    "gene80", "gene32", "gene55", "gene76", " ", " ",
    "gene306", "gene329", "gene397", "gene320", "gene382", "gene369",
    "gene379", "gene398", "gene387", "gene315", "gene364", "gene390",
    "gene334", "gene363", "gene326", "gene314", "gene386", "gene316",
    "gene393", "gene312", " ", " ", "gene111", "gene166",
    "gene133", "gene153", "gene126", "gene108", "gene132", "gene180",
    "gene106", "gene144", "gene186", "gene192", "gene135", "gene165",
    "gene150", "gene117", "gene155", "gene190", "gene146", "gene189",
    " ", " ", "gene206", "gene275", "gene277", "gene289",
    "gene222", "gene241", "gene214", "gene223", "gene236", "gene282",
    "gene225", "gene256", "gene278", "gene249", "gene286", "gene230",
    "gene245", "gene272", "gene273", "gene285"
  ))

  # Checking title of the heatmap
  htmp_title <- htmp@annotations@listData$col_title2@data
  expect_equal(htmp_title, "Expression heatmap")

  # Checking x lab of the heatmap
  htmp_xlab <- htmp@annotations@listData$col_title@data
  expect_equal(htmp_xlab, "Samples")

  # Checking y lab of the heatmap
  htmp_ylab <- htmp@annotations@listData$row_title@data
  expect_equal(htmp_ylab, "Genes")

  # Checking colorbar name of the heatmap
  htmp_colors_name <- htmp@colorbars@listData$`Exp. level`@title
  expect_equal(htmp_colors_name, "Exp. level")

  # Checking if colorbar is shown of the heatmap
  htmp_show_colorbar <- htmp@plots@listData$`Exp. level`@show_colorbar
  expect_true(htmp_show_colorbar)

  # Checking if dendrogram is showed
  htmp_show_dendro <- names(htmp@shapes@listData)
  expect_equal(htmp_show_dendro, NULL)
})


# ==============================================================================
# ==============================================================================


test_that("Checking plot_heatmap() increasing size of the blank lines", {
  cell_clusters <- c(rep(1, 10), rep(2, 10))
  names(cell_clusters) <- colnames(res_scigenex@data)

  res_scigenex_top <- top_genes(res_scigenex, top = 20)

  htmp <- plot_heatmap(
    object = res_scigenex_top,
    cell_clusters = cell_clusters,
    use_top_genes = TRUE,
    show_dendro = TRUE,
    line_size_vertical = 10,
    line_size_horizontal = 10,
    name = "Expression heatmap",
    show_legend = TRUE,
    xlab = "Samples",
    ylab = "Genes"
  )

  htmp_matrix <- htmp@plots@listData$`Exp. level`@data
  htmp_matrix_wo_na <- htmp_matrix[!is.na(htmp_matrix)]

  # Check heatmpa format
  expect_that(htmp, is_a("IheatmapHorizontal"))

  # Check dim of heatmap matrix
  expect_equal(ncol(htmp_matrix), 30)
  expect_equal(nrow(htmp_matrix), 110)
  expect_equal(dim(htmp_matrix), c(110, 30))

  # Check NA rows added by line_size_horizontal
  na_rows <- which(apply(htmp_matrix, 1, function(x) all(is.na(x))))
  expect_equal(length(na_rows), 30)
  expect_equal(as.vector(na_rows), c(
    seq(21, 30),
    seq(51, 60),
    seq(81, 90)
  ))

  # Check NA rows added by line_size_horizontal
  na_cols <- which(apply(htmp_matrix, 2, function(x) all(is.na(x))))
  expect_equal(length(na_cols), 10)
  expect_equal(as.vector(na_cols), c(seq(11, 20)))

  # Check some measurements of heatmap matrix
  expect_equal(round(sum(htmp_matrix_wo_na), 2), 64.5)
  expect_equal(round(mean(htmp_matrix_wo_na), 2), 0.04)
  expect_equal(round(median(htmp_matrix_wo_na), 2), 0.24)
  expect_equal(round(sd(htmp_matrix_wo_na), 2), 0.89)

  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["0%"]), 2), -1)
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["25%"]), 2), -1)
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["50%"]), 2), 0.24)
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["75%"]), 2), 1)
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["100%"]), 2), 1)

  # Checking order of columns
  expect_equal(colnames(htmp_matrix), c("sample1", "sample2", "sample3", "sample4", "sample5", "sample6", 
                                        "sample7", "sample8", "sample9", "sample10", " ", "  ", "   ", 
                                        "    ", "     ", "      ", "       ", "        ", "         ", 
                                        "          ", "sample11", "sample12", "sample13", "sample14", 
                                        "sample15", "sample16", "sample17", "sample18", "sample19", "sample20"))

  # Checking order of rows
  expect_equal(rownames(htmp_matrix), c(
    "gene7", "gene19",
    "gene57", "gene36", "gene78", "gene39", "gene67", "gene88",
    "gene16", "gene84", "gene29", "gene37", "gene27", "gene70",
    "gene59", "gene12", "gene80", "gene32", "gene55", "gene76",
    " ", " ", " ", " ", " ", " ",
    " ", " ", " ", " ", "gene306", "gene329",
    "gene397", "gene320", "gene382", "gene369", "gene379", "gene398",
    "gene387", "gene315", "gene364", "gene390", "gene334", "gene363",
    "gene326", "gene314", "gene386", "gene316", "gene393", "gene312",
    " ", " ", " ", " ", " ", " ",
    " ", " ", " ", " ", "gene111", "gene166",
    "gene133", "gene153", "gene126", "gene108", "gene132", "gene180",
    "gene106", "gene144", "gene186", "gene192", "gene135", "gene165",
    "gene150", "gene117", "gene155", "gene190", "gene146", "gene189",
    " ", " ", " ", " ", " ", " ",
    " ", " ", " ", " ", "gene206", "gene275",
    "gene277", "gene289", "gene222", "gene241", "gene214", "gene223",
    "gene236", "gene282", "gene225", "gene256", "gene278", "gene249",
    "gene286", "gene230", "gene245", "gene272", "gene273", "gene285"
  ))

  # Checking title of the heatmap
  htmp_title <- htmp@annotations@listData$col_title2@data
  expect_equal(htmp_title, "Expression heatmap")

  # Checking x lab of the heatmap
  htmp_xlab <- htmp@annotations@listData$col_title@data
  expect_equal(htmp_xlab, "Samples")

  # Checking y lab of the heatmap
  htmp_ylab <- htmp@annotations@listData$row_title@data
  expect_equal(htmp_ylab, "Genes")

  # Checking colorbar name of the heatmap
  htmp_colors_name <- htmp@colorbars@listData$`Exp. level`@title
  expect_equal(htmp_colors_name, "Exp. level")

  # Checking if colorbar is shown of the heatmap
  htmp_show_colorbar <- htmp@plots@listData$`Exp. level`@show_colorbar
  expect_true(htmp_show_colorbar)

  # Checking if dendrogram is showed
  htmp_show_dendro <- names(htmp@shapes@listData)
  expect_equal(htmp_show_dendro, NULL)
})


# ==============================================================================
# ==============================================================================


test_that("Checking plot_heatmap() displaying gene clusters 1 and 4", {
  cell_clusters <- c(rep(1, 10), rep(2, 10))
  names(cell_clusters) <- colnames(res_scigenex@data)

  res_scigenex_top <- top_genes(res_scigenex, top = 20)

  htmp <- plot_heatmap(
    object = res_scigenex_top[c(1, 4), ],
    cell_clusters = cell_clusters,
    use_top_genes = TRUE,
    show_dendro = TRUE,
    line_size_vertical = 2,
    line_size_horizontal = 2,
    name = "Expression heatmap",
    show_legend = TRUE,
    xlab = "Samples",
    ylab = "Genes"
  )

  htmp_matrix <- htmp@plots@listData$`Exp. level`@data
  htmp_matrix_wo_na <- htmp_matrix[!is.na(htmp_matrix)]

  # Check heatmpa format
  expect_that(htmp, is_a("IheatmapHorizontal"))

  # Check dim of heatmap matrix
  expect_equal(ncol(htmp_matrix), 22)
  expect_equal(nrow(htmp_matrix), 42)
  expect_equal(dim(htmp_matrix), c(42, 22))

  # Check NA rows added by line_size_horizontal
  na_rows <- which(apply(htmp_matrix, 1, function(x) all(is.na(x))))
  expect_equal(length(na_rows), 2)
  expect_equal(as.vector(na_rows), c(21, 22))

  # Check NA rows added by line_size_horizontal
  na_cols <- which(apply(htmp_matrix, 2, function(x) all(is.na(x))))
  expect_equal(length(na_cols), 2)
  expect_equal(as.vector(na_cols), c(11, 12))

  # Check some measurements of heatmap matrix
  expect_equal(round(sum(htmp_matrix_wo_na), 2), 77.15)
  expect_equal(round(mean(htmp_matrix_wo_na), 2), 0.1)
  expect_equal(round(median(htmp_matrix_wo_na), 2), 0.36)
  expect_equal(round(sd(htmp_matrix_wo_na), 2), 0.85)

  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["0%"]), 2), -1)
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["25%"]), 2), -1)
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["50%"]), 2), 0.36)
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["75%"]), 2), 1)
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["100%"]), 2), 1)

  # Checking order of columns
  expect_equal(colnames(htmp_matrix), c("sample1", "sample2", "sample3", "sample4", "sample5", "sample6", 
                                        "sample7", "sample8", "sample9", "sample10", " ", " ", "sample11", 
                                        "sample12", "sample13", "sample14", "sample15", "sample16", "sample17", 
                                        "sample18", "sample19", "sample20"))

  # Checking order of rows
  expect_equal(rownames(htmp_matrix), c(
    "gene7", "gene19", "gene57", "gene36",
    "gene78", "gene39", "gene67", "gene88", "gene16", "gene84",
    "gene29", "gene37", "gene27", "gene70", "gene59", "gene12",
    "gene80", "gene32", "gene55", "gene76", " ", " ",
    "gene206", "gene275", "gene277", "gene289", "gene222", "gene241",
    "gene214", "gene223", "gene236", "gene282", "gene225", "gene256",
    "gene278", "gene249", "gene286", "gene230", "gene245", "gene272",
    "gene273", "gene285"
  ))

  # Checking title of the heatmap
  htmp_title <- htmp@annotations@listData$col_title2@data
  expect_equal(htmp_title, "Expression heatmap")

  # Checking x lab of the heatmap
  htmp_xlab <- htmp@annotations@listData$col_title@data
  expect_equal(htmp_xlab, "Samples")

  # Checking y lab of the heatmap
  htmp_ylab <- htmp@annotations@listData$row_title@data
  expect_equal(htmp_ylab, "Genes")

  # Checking colorbar name of the heatmap
  htmp_colors_name <- htmp@colorbars@listData$`Exp. level`@title
  expect_equal(htmp_colors_name, "Exp. level")

  # Checking if colorbar is shown of the heatmap
  htmp_show_colorbar <- htmp@plots@listData$`Exp. level`@show_colorbar
  expect_true(htmp_show_colorbar)

  # Checking if dendrogram is showed
  htmp_show_dendro <- names(htmp@shapes@listData)
  expect_equal(htmp_show_dendro, NULL)
})


# ==============================================================================
# ==============================================================================


test_that("Checking plot_heatmap() displaying only gene cluster 3", {
  
  cell_clusters <- c(rep(1, 10), rep(2, 10))
  names(cell_clusters) <- colnames(res_scigenex@data)
  res_scigenex_top <- top_genes(res_scigenex, top = 20)

  htmp <- plot_heatmap(
    object = res_scigenex_top[3,],
    cell_clusters = cell_clusters,
    use_top_genes = TRUE,
    show_dendro = TRUE,
    line_size_vertical = 2,
    line_size_horizontal = 2,
    name = "Expression heatmap",
    show_legend = TRUE,
    xlab = "Samples",
    ylab = "Genes"
  )

  htmp_matrix <- htmp@plots@listData$`Exp. level`@data
  htmp_matrix_wo_na <- htmp_matrix[!is.na(htmp_matrix)]

  # Check heatmpa format
  expect_that(htmp, is_a("IheatmapHorizontal"))

  # Check dim of heatmap matrix
  expect_equal(ncol(htmp_matrix), 22)
  expect_equal(nrow(htmp_matrix), 20)
  expect_equal(dim(htmp_matrix), c(20, 22))

  # Check NA rows added by line_size_horizontal
  na_rows <- which(apply(htmp_matrix, 1, function(x) all(is.na(x))))
  expect_equal(length(na_rows), 0)
  expect_equal(as.vector(na_rows), integer(0))

  # Check NA rows added by line_size_horizontal
  na_cols <- which(apply(htmp_matrix, 2, function(x) all(is.na(x))))
  expect_equal(length(na_cols), 2)
  expect_equal(as.vector(na_cols), c(11, 12))

  # Check some measurements of heatmap matrix
  expect_equal(round(sum(htmp_matrix_wo_na), 2), -10.54)
  expect_equal(round(mean(htmp_matrix_wo_na), 2), -0.03)
  expect_equal(round(median(htmp_matrix_wo_na), 2), -0.14)
  expect_equal(round(sd(htmp_matrix_wo_na), 2), 0.88)

  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["0%"]), 2), -1)
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["25%"]), 2), -1)
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["50%"]), 2), -0.14)
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["75%"]), 2), 1)
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["100%"]), 2), 1)

  # Checking order of columns
  expect_equal(colnames(htmp_matrix), c(
    "sample1", "sample2", "sample3", "sample4", "sample5",
    "sample6", "sample7", "sample8", "sample9", "sample10", " ",
    " ", "sample11", "sample12", "sample13", "sample14", "sample15",
    "sample16", "sample17", "sample18", "sample19", "sample20"
  ))

  # Checking order of rows
  expect_equal(rownames(htmp_matrix), c(
    "gene306", "gene329", "gene397", "gene320",
    "gene382", "gene369", "gene379", "gene398", "gene387", "gene315",
    "gene364", "gene390", "gene334", "gene363", "gene326", "gene314",
    "gene386", "gene316", "gene393", "gene312"
  ))

  # Checking title of the heatmap
  htmp_title <- htmp@annotations@listData$col_title2@data
  expect_equal(htmp_title, "Expression heatmap")

  # Checking x lab of the heatmap
  htmp_xlab <- htmp@annotations@listData$col_title@data
  expect_equal(htmp_xlab, "Samples")

  # Checking y lab of the heatmap
  htmp_ylab <- htmp@annotations@listData$row_title@data
  expect_equal(htmp_ylab, "Genes")

  # Checking colorbar name of the heatmap
  htmp_colors_name <- htmp@colorbars@listData$`Exp. level`@title
  expect_equal(htmp_colors_name, "Exp. level")

  # Checking if colorbar is shown of the heatmap
  htmp_show_colorbar <- htmp@plots@listData$`Exp. level`@show_colorbar
  expect_true(htmp_show_colorbar)

  # Checking if dendrogram is showed
  htmp_show_dendro <- names(htmp@shapes@listData)
  expect_equal(htmp_show_dendro, NULL)
})


testthat::test_that("Check plot_heatmap and subsetting", {
  
  set_verbosity(0)
  data(pbmc_small, package = "SeuratObject")
  # Compute the signatures using find_gene_clusters()
  ident <- Seurat::Idents(pbmc_small)
  
  clust_set <- find_gene_clusters(pbmc_small, k=50, no_dknn_filter=TRUE)
  
  expect_error(plot_heatmap(clust_set, cell_clusters=Seurat::Idents(pbmc_small), colors_cell_clusters = "green"))
  p <- plot_heatmap(clust_set[, names(ident[ident==2])], cell_clusters=Seurat::Idents(pbmc_small), colors_cell_clusters = "red")
  testthat::expect_error(print(p), NA)
  expect_error(plot_heatmap(clust_set[, names(ident[ident==3])], cell_clusters=Seurat::Idents(pbmc_small)))
  p <- plot_heatmap(clust_set, cell_clusters=Seurat::Idents(pbmc_small))
  testthat::expect_error(print(p), NA)
  p <- plot_heatmap(clust_set[,colnames(clust_set@data)[1:80]], cell_clusters=Seurat::Idents(pbmc_small))
  testthat::expect_error(print(p), NA)
  p <- plot_heatmap(clust_set[,colnames(clust_set@data)[1:10]], cell_clusters=Seurat::Idents(pbmc_small))
  testthat::expect_error(print(p), NA)
  p <- plot_heatmap(clust_set[,colnames(clust_set@data)[1:20]], cell_clusters=Seurat::Idents(pbmc_small))
  testthat::expect_error(print(p), NA)
  p <- plot_heatmap(clust_set[1,colnames(clust_set@data)[1:10]], cell_clusters=Seurat::Idents(pbmc_small))
  testthat::expect_error(print(p), NA)
  p <- plot_heatmap(clust_set[2 ,colnames(clust_set@data)[1:10]], cell_clusters=Seurat::Idents(pbmc_small))
  testthat::expect_error(print(p), NA)
  p <- plot_heatmap(clust_set[2 ,colnames(clust_set@data)[1:10]], cell_clusters=Seurat::Idents(pbmc_small))
  testthat::expect_error(print(p), NA)
  p <- plot_heatmap(clust_set[1,colnames(clust_set@data)[1:10]], cell_clusters=Seurat::Idents(pbmc_small))
  testthat::expect_error(print(p), NA)
})

testthat::test_that("Check plot_heatmap, subsetting and clustering", {
  
  set_verbosity(0)
  data(pbmc_small, package = "SeuratObject")
  # Compute the signatures using find_gene_clusters()
  ident <- Seurat::Idents(pbmc_small)
  
  clust_set <- find_gene_clusters(pbmc_small, k=50, no_dknn_filter=TRUE)
  
  expect_error(plot_heatmap(clust_set, cell_clusters=Seurat::Idents(pbmc_small), colors_cell_clusters = "green"))
  plot_heatmap(clust_set[, names(ident[ident==2])], cell_clusters=Seurat::Idents(pbmc_small), colors_cell_clusters = "red")
  
  expect_error(plot_heatmap(clust_set[, names(ident[ident==3])], cell_clusters=Seurat::Idents(pbmc_small)))
  p <- plot_heatmap(clust_set, cell_clusters=Seurat::Idents(pbmc_small))
  testthat::expect_error(print(p), NA)
  p <- plot_heatmap(clust_set[,colnames(clust_set@data)[1:80]], cell_clusters=Seurat::Idents(pbmc_small))
  testthat::expect_error(print(p), NA)
  p <- plot_heatmap(clust_set[,colnames(clust_set@data)[1:10]], cell_clusters=Seurat::Idents(pbmc_small))
  testthat::expect_error(print(p), NA)
  p <- plot_heatmap(clust_set[,colnames(clust_set@data)[1:20]], cell_clusters=Seurat::Idents(pbmc_small))
  testthat::expect_error(print(p), NA)
  p <- plot_heatmap(clust_set[1,colnames(clust_set@data)[1:10]], cell_clusters=Seurat::Idents(pbmc_small))
  testthat::expect_error(print(p), NA)
  p <- plot_heatmap(clust_set[2 ,colnames(clust_set@data)[1:10]], cell_clusters=Seurat::Idents(pbmc_small))
  testthat::expect_error(print(p), NA)
  p <- plot_heatmap(clust_set[2 ,colnames(clust_set@data)[1:10]], cell_clusters=Seurat::Idents(pbmc_small))
  testthat::expect_error(print(p), NA)
  
})


testthat::test_that("Test non interactive heatmaps.", {
  
    set_verbosity(0)
  data(pbmc_small, package = "SeuratObject")
  # Compute the signatures using find_gene_clusters()
  ident <- Seurat::Idents(pbmc_small)
  
  clust_set <- find_gene_clusters(pbmc_small, k=50, no_dknn_filter=TRUE)
  
  n_gene_clusters <- 5 
  n_cell_clusters <- 3
  p_int <- plot_heatmap(clust_set, cell_clusters=Seurat::Idents(pbmc_small), 
                        colors_cell_clusters =  c("green", "red", "yellow"))
  p_noi <- plot_heatmap(clust_set, cell_clusters=Seurat::Idents(pbmc_small), 
                        colors_cell_clusters =  c("green", "red", "yellow"),
                        interactive = FALSE)
  p_int_mat <- htmp_matrix <- p_int@plots@listData$`Exp. level`@data
  p_noi_rown <- p_noi[["gtable"]][["grobs"]][[2]][["label"]]
  testthat::expect_equal(p_noi_rown, rev(rownames(p_int_mat)))
  
  p_noi_colnames <- names(p_noi[["gtable"]][["grobs"]][[3]][["gp"]][["fill"]][,1])
  expect_true(all(colnames(p_int_mat) == p_noi_colnames))
  expect_true(length(grep(" ", p_noi_colnames)) == (n_cell_clusters-1) * 3)
  expect_true(length(grep(" ", p_noi_rown)) == (n_gene_clusters-1) * 3)
  
  n_gene_clusters <- 1 
  n_cell_clusters <- 3
  p_noi <- plot_heatmap(clust_set[1,], cell_clusters=Seurat::Idents(pbmc_small), 
                        colors_cell_clusters =  c("green", "red", "yellow"),
                        interactive = FALSE)
  p_noi_rown <- p_noi[["gtable"]][["grobs"]][[2]][["label"]]
  expect_true(length(grep(" ", p_noi_colnames)) == (n_cell_clusters-1) * 3)
  expect_true(length(grep(" ", p_noi_rown)) == (n_gene_clusters-1) * 3)
  
  n_gene_clusters <- 2 
  n_cell_clusters <- 3
  p_noi <- plot_heatmap(clust_set[c(1, 3), colnames(pbmc_small)], cell_clusters=Seurat::Idents(pbmc_small), 
                        colors_cell_clusters =  c("green", "red", "yellow"),
                        interactive = FALSE)
  p_noi_rown <- p_noi[["gtable"]][["grobs"]][[2]][["label"]]
  expect_true(length(grep(" ", p_noi_colnames)) == (n_cell_clusters-1) * 3)
  expect_true(length(grep(" ", p_noi_rown)) == (n_gene_clusters-1) * 3)
  
})
