# Set verbosity to 0
set_verbosity(0)

#Create matrix containing 3 signatures
m <- create_4_rnd_clust()

## Select informative genes
res <- select_genes(data=m,
                    distance_method="kendall",
                    k=75,
                    row_sum=-Inf,
                    highest=0.3,
                    dist_threads = 6,
                    fdr = 1e-8)

## Cluster genes
res <- gene_clustering(object = res,
                       inflation = 1.2,
                       keep_nn = FALSE,
                       k = 5,
                       threads = 6)




test_that("Checking plot_heatmap()", {
  set.seed(123)
  
  htmp <- plot_heatmap(
    object = res,
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
  expect_equal(nrow(htmp_matrix), 365)
  expect_equal(dim(htmp_matrix), c(365, 20))
  
  # Check NA rows added by line_size_horizontal
  na_rows <- which(apply(htmp_matrix, 1, function(x) all(is.na(x))))
  expect_equal(length(na_rows), 6)
  expect_equal(as.vector(na_rows), c(68, 69, 151, 152, 241, 242))
  
  # Check NA rows added by line_size_vertical
  na_cols <- which(apply(htmp_matrix, 2, function(x) all(is.na(x))))
  expect_equal(length(na_cols), 0)
  expect_equal(as.vector(na_cols), as.integer())
  
  # Check some measurements of heatmap matrix
  expect_equal(round(sum(htmp_matrix_wo_na), 2), 188.56)
  expect_equal(round(mean(htmp_matrix_wo_na), 2), 0.03)
  expect_equal(round(median(htmp_matrix_wo_na), 2), 0.1)
  expect_equal(round(sd(htmp_matrix_wo_na), 2), 0.86)
  
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["0%"]), 2), -1)
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["25%"]), 2), -1)
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["50%"]), 2), 0.1)
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["75%"]), 2), 1)
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["100%"]), 2), 1)
  
  # Checking order of columns
  expect_equal(colnames(htmp_matrix), c(
    "sample14", "sample11", "sample12", "sample13", "sample17", "sample18",
    "sample16", "sample19", "sample15", "sample20", "sample1", "sample5",
    "sample10", "sample4", "sample2", "sample3", "sample7", "sample6",
    "sample8", "sample9"
  ))
  
  # Checking order of rows
  expect_equal(rownames(htmp_matrix), c(
    "gene98", "gene52", "gene99", "gene43", "gene3", "gene95",
    "gene54", "gene85", "gene35", "gene93", "gene82", "gene20",
    "gene11", "gene63", "gene46", "gene56", "gene13", "gene14",
    "gene4", "gene71", "gene1", "gene42", "gene80", "gene10",
    "gene32", "gene96", "gene66", "gene9", "gene94", "gene70",
    "gene45", "gene78", "gene34", "gene29", "gene28", "gene25",
    "gene36", "gene83", "gene17", "gene79", "gene51", "gene7",
    "gene12", "gene59", "gene57", "gene67", "gene2", "gene73",
    "gene88", "gene16", "gene39", "gene64", "gene21", "gene55",
    "gene90", "gene74", "gene61", "gene76", "gene86", "gene6",
    "gene37", "gene31", "gene58", "gene87", "gene27", "gene84",
    "gene19", " ", " ", "gene2490", "gene1423", "gene760",
    "gene396", "gene375", "gene325", "gene342", "gene310", "gene302",
    "gene2911", "gene370", "gene360", "gene376", "gene328", "gene344",
    "gene345", "gene317", "gene392", "gene3609", "gene400", "gene372",
    "gene335", "gene355", "gene353", "gene351", "gene341", "gene367",
    "gene394", "gene301", "gene365", "gene327", "gene350", "gene371",
    "gene389", "gene331", "gene330", "gene359", "gene336", "gene313",
    "gene340", "gene378", "gene306", "gene339", "gene373", "gene337",
    "gene303", "gene349", "gene352", "gene304", "gene329", "gene382",
    "gene387", "gene316", "gene358", "gene397", "gene311", "gene369",
    "gene393", "gene332", "gene379", "gene366", "gene398", "gene390",
    "gene319", "gene386", "gene322", "gene362", "gene385", "gene364",
    "gene320", "gene324", "gene391", "gene314", "gene363", "gene368",
    "gene318", "gene323", "gene326", "gene312", "gene334", "gene315",
    " ", " ", "gene1931", "gene292", "gene268", "gene252",
    "gene234", "gene280", "gene213", "gene299", "gene287", "gene254",
    "gene276", "gene228", "gene246", "gene297", "gene293", "gene255",
    "gene275", "gene201", "gene247", "gene296", "gene259", "gene266",
    "gene294", "gene241", "gene281", "gene283", "gene240", "gene238",
    "gene261", "gene212", "gene220", "gene203", "gene224", "gene3991",
    "gene227", "gene249", "gene289", "gene257", "gene207", "gene288",
    "gene300", "gene236", "gene223", "gene263", "gene267", "gene282",
    "gene248", "gene277", "gene258", "gene269", "gene290", "gene262",
    "gene217", "gene285", "gene273", "gene245", "gene253", "gene298",
    "gene208", "gene242", "gene237", "gene256", "gene225", "gene235",
    "gene291", "gene233", "gene286", "gene232", "gene272", "gene229",
    "gene206", "gene251", "gene244", "gene214", "gene230", "gene219",
    "gene226", "gene278", "gene250", "gene239", "gene2930", "gene202",
    "gene284", "gene279", "gene215", "gene222", "gene211", "gene270",
    " ", " ", "gene3924", "gene3779", "gene3610", "gene2902",
    "gene2553", "gene2421", "gene2662", "gene2223", "gene2112", "gene1915",
    "gene1651", "gene3488", "gene3485", "gene3442", "gene2725", "gene2573",
    "gene2184", "gene3528", "gene1503", "gene1856", "gene1616", "gene838",
    "gene1922", "gene156", "gene187", "gene174", "gene170", "gene200",
    "gene181", "gene172", "gene169", "gene185", "gene134", "gene123",
    "gene110", "gene130", "gene102", "gene128", "gene179", "gene177",
    "gene152", "gene125", "gene149", "gene140", "gene196", "gene131",
    "gene157", "gene129", "gene113", "gene175", "gene147", "gene182",
    "gene167", "gene104", "gene184", "gene195", "gene199", "gene109",
    "gene193", "gene183", "gene115", "gene163", "gene158", "gene154",
    "gene151", "gene138", "gene139", "gene161", "gene132", "gene164",
    "gene136", "gene168", "gene162", "gene165", "gene118", "gene124",
    "gene178", "gene114", "gene171", "gene191", "gene121", "gene190",
    "gene137", "gene120", "gene127", "gene197", "gene166", "gene148",
    "gene116", "gene160", "gene153", "gene194", "gene112", "gene133",
    "gene145", "gene111", "gene122", "gene108", "gene159", "gene135",
    "gene107", "gene106", "gene180", "gene173", "gene186", "gene103",
    "gene119", "gene101", "gene188", "gene142", "gene105", "gene198",
    "gene117", "gene143", "gene146", "gene144", "gene150", "gene141",
    "gene176", "gene192", "gene189", "gene155", "gene126"
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
  
  res_scigenex_top <- top_genes(res, top = 20)
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
  expect_equal(round(sum(htmp_matrix_wo_na), 2), 52.05)
  expect_equal(round(mean(htmp_matrix_wo_na), 2), 0.03)
  expect_equal(round(median(htmp_matrix_wo_na), 2), 0.15)
  expect_equal(round(sd(htmp_matrix_wo_na), 2), 0.88)
  
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["0%"]), 2), -1)
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["25%"]), 2), -1)
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["50%"]), 2), 0.15)
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["75%"]), 2), 1)
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["100%"]), 2), 1)
  
  # Checking order of columns
  expect_equal(colnames(htmp_matrix), c(
    "sample14", "sample11", "sample12", "sample13", "sample17", "sample18",
    "sample16", "sample19", "sample15", "sample20", "sample1", "sample5",
    "sample10", "sample4", "sample2", "sample3", "sample7", "sample6",
    "sample8", "sample9"
  ))
  
  # Checking order of rows
  expect_equal(rownames(htmp_matrix), c(
    "gene2", "gene70", "gene12", "gene29", "gene31", "gene80",
    "gene61", "gene55", "gene27", "gene57", "gene78", "gene7",
    "gene88", "gene76", "gene32", "gene64", "gene51", "gene19",
    "gene84", "gene37", " ", " ", "gene326", "gene363",
    "gene349", "gene369", "gene339", "gene397", "gene382", "gene387",
    "gene350", "gene316", "gene368", "gene393", "gene398", "gene334",
    "gene314", "gene320", "gene390", "gene312", "gene364", "gene315",
    " ", " ", "gene288", "gene232", "gene202", "gene212",
    "gene229", "gene269", "gene248", "gene219", "gene277", "gene285",
    "gene283", "gene282", "gene281", "gene256", "gene273", "gene207",
    "gene206", "gene223", "gene245", "gene286", " ", " ",
    "gene115", "gene119", "gene155", "gene163", "gene158", "gene198",
    "gene109", "gene122", "gene120", "gene126", "gene121", "gene176",
    "gene141", "gene143", "gene171", "gene116", "gene137", "gene108",
    "gene191", "gene160"
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
  names(cell_clusters) <- colnames(res@data)
  htmp <- plot_heatmap(
    object = res,
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
  expect_equal(nrow(htmp_matrix), 365)
  expect_equal(dim(htmp_matrix), c(365, 22))
  
  # Check NA rows added by line_size_horizontal
  na_rows <- which(apply(htmp_matrix, 1, function(x) all(is.na(x))))
  expect_equal(length(na_rows), 6)
  expect_equal(as.vector(na_rows), c(68, 69, 151, 152, 241, 242))
  
  # Check NA rows added by line_size_horizontal
  na_cols <- which(apply(htmp_matrix, 2, function(x) all(is.na(x))))
  expect_equal(length(na_cols), 2)
  expect_equal(as.vector(na_cols), c(11, 12))
  
  # Check some measurements of heatmap matrix
  expect_equal(round(sum(htmp_matrix_wo_na), 2), 188.56)
  expect_equal(round(mean(htmp_matrix_wo_na), 2), 0.03)
  expect_equal(round(median(htmp_matrix_wo_na), 2), 0.1)
  expect_equal(round(sd(htmp_matrix_wo_na), 2), 0.86)
  
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["0%"]), 2), -1)
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["25%"]), 2), -1)
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["50%"]), 2), 0.1)
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["75%"]), 2), 1)
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["100%"]), 2), 1)
  
  # Checking order of columns
  expect_equal(colnames(htmp_matrix), c(
    "sample1", "sample2", "sample3", "sample4", "sample5", "sample6",
    "sample7", "sample8", "sample9", "sample10", " ", "  ",
    "sample11", "sample12", "sample13", "sample14", "sample15", "sample16",
    "sample17", "sample18", "sample19", "sample20"
  ))
  
  # Checking order of rows
  expect_equal(rownames(htmp_matrix), c(
    "gene98", "gene52", "gene99", "gene43", "gene3", "gene95",
    "gene54", "gene85", "gene35", "gene93", "gene82", "gene20",
    "gene11", "gene63", "gene46", "gene56", "gene13", "gene14",
    "gene4", "gene71", "gene1", "gene42", "gene80", "gene10",
    "gene32", "gene96", "gene66", "gene9", "gene94", "gene70",
    "gene45", "gene78", "gene34", "gene29", "gene28", "gene25",
    "gene36", "gene83", "gene17", "gene79", "gene51", "gene7",
    "gene12", "gene59", "gene57", "gene67", "gene2", "gene73",
    "gene88", "gene16", "gene39", "gene64", "gene21", "gene55",
    "gene90", "gene74", "gene61", "gene76", "gene86", "gene6",
    "gene37", "gene31", "gene58", "gene87", "gene27", "gene84",
    "gene19", " ", " ", "gene2490", "gene1423", "gene760",
    "gene396", "gene375", "gene325", "gene342", "gene310", "gene302",
    "gene2911", "gene370", "gene360", "gene376", "gene328", "gene344",
    "gene345", "gene317", "gene392", "gene3609", "gene400", "gene372",
    "gene335", "gene355", "gene353", "gene351", "gene341", "gene367",
    "gene394", "gene301", "gene365", "gene327", "gene350", "gene371",
    "gene389", "gene331", "gene330", "gene359", "gene336", "gene313",
    "gene340", "gene378", "gene306", "gene339", "gene373", "gene337",
    "gene303", "gene349", "gene352", "gene304", "gene329", "gene382",
    "gene387", "gene316", "gene358", "gene397", "gene311", "gene369",
    "gene393", "gene332", "gene379", "gene366", "gene398", "gene390",
    "gene319", "gene386", "gene322", "gene362", "gene385", "gene364",
    "gene320", "gene324", "gene391", "gene314", "gene363", "gene368",
    "gene318", "gene323", "gene326", "gene312", "gene334", "gene315",
    " ", " ", "gene1931", "gene292", "gene268", "gene252",
    "gene234", "gene280", "gene213", "gene299", "gene287", "gene254",
    "gene276", "gene228", "gene246", "gene297", "gene293", "gene255",
    "gene275", "gene201", "gene247", "gene296", "gene259", "gene266",
    "gene294", "gene241", "gene281", "gene283", "gene240", "gene238",
    "gene261", "gene212", "gene220", "gene203", "gene224", "gene3991",
    "gene227", "gene249", "gene289", "gene257", "gene207", "gene288",
    "gene300", "gene236", "gene223", "gene263", "gene267", "gene282",
    "gene248", "gene277", "gene258", "gene269", "gene290", "gene262",
    "gene217", "gene285", "gene273", "gene245", "gene253", "gene298",
    "gene208", "gene242", "gene237", "gene256", "gene225", "gene235",
    "gene291", "gene233", "gene286", "gene232", "gene272", "gene229",
    "gene206", "gene251", "gene244", "gene214", "gene230", "gene219",
    "gene226", "gene278", "gene250", "gene239", "gene2930", "gene202",
    "gene284", "gene279", "gene215", "gene222", "gene211", "gene270",
    " ", " ", "gene3924", "gene3779", "gene3610", "gene2902",
    "gene2553", "gene2421", "gene2662", "gene2223", "gene2112", "gene1915",
    "gene1651", "gene3488", "gene3485", "gene3442", "gene2725", "gene2573",
    "gene2184", "gene3528", "gene1503", "gene1856", "gene1616", "gene838",
    "gene1922", "gene156", "gene187", "gene174", "gene170", "gene200",
    "gene181", "gene172", "gene169", "gene185", "gene134", "gene123",
    "gene110", "gene130", "gene102", "gene128", "gene179", "gene177",
    "gene152", "gene125", "gene149", "gene140", "gene196", "gene131",
    "gene157", "gene129", "gene113", "gene175", "gene147", "gene182",
    "gene167", "gene104", "gene184", "gene195", "gene199", "gene109",
    "gene193", "gene183", "gene115", "gene163", "gene158", "gene154",
    "gene151", "gene138", "gene139", "gene161", "gene132", "gene164",
    "gene136", "gene168", "gene162", "gene165", "gene118", "gene124",
    "gene178", "gene114", "gene171", "gene191", "gene121", "gene190",
    "gene137", "gene120", "gene127", "gene197", "gene166", "gene148",
    "gene116", "gene160", "gene153", "gene194", "gene112", "gene133",
    "gene145", "gene111", "gene122", "gene108", "gene159", "gene135",
    "gene107", "gene106", "gene180", "gene173", "gene186", "gene103",
    "gene119", "gene101", "gene188", "gene142", "gene105", "gene198",
    "gene117", "gene143", "gene146", "gene144", "gene150", "gene141",
    "gene176", "gene192", "gene189", "gene155", "gene126"
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
  names(cell_clusters) <- colnames(res@data)
  
  res_scigenex_top <- top_genes(res, top = 20)
  
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
  expect_equal(round(sum(htmp_matrix_wo_na), 2), 52.05)
  expect_equal(round(mean(htmp_matrix_wo_na), 2), 0.03)
  expect_equal(round(median(htmp_matrix_wo_na), 2), 0.15)
  expect_equal(round(sd(htmp_matrix_wo_na), 2), 0.88)
  
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["0%"]), 2), -1)
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["25%"]), 2), -1)
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["50%"]), 2), 0.15)
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["75%"]), 2), 1)
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["100%"]), 2), 1)
  
  # Checking order of columns
  expect_equal(colnames(htmp_matrix), c("sample1", "sample2", "sample3", "sample4", "sample5", "sample6", 
                                        "sample7", "sample8", "sample9", "sample10", " ", "  ", "sample11", 
                                        "sample12", "sample13", "sample14", "sample15", "sample16", "sample17", 
                                        "sample18", "sample19", "sample20"))
  
  # Checking order of rows
  expect_equal(rownames(htmp_matrix), c(
    "gene2", "gene70", "gene12", "gene29", "gene31", "gene80",
    "gene61", "gene55", "gene27", "gene57", "gene78", "gene7",
    "gene88", "gene76", "gene32", "gene64", "gene51", "gene19",
    "gene84", "gene37", " ", " ", "gene326", "gene363",
    "gene349", "gene369", "gene339", "gene397", "gene382", "gene387",
    "gene350", "gene316", "gene368", "gene393", "gene398", "gene334",
    "gene314", "gene320", "gene390", "gene312", "gene364", "gene315",
    " ", " ", "gene288", "gene232", "gene202", "gene212",
    "gene229", "gene269", "gene248", "gene219", "gene277", "gene285",
    "gene283", "gene282", "gene281", "gene256", "gene273", "gene207",
    "gene206", "gene223", "gene245", "gene286", " ", " ",
    "gene115", "gene119", "gene155", "gene163", "gene158", "gene198",
    "gene109", "gene122", "gene120", "gene126", "gene121", "gene176",
    "gene141", "gene143", "gene171", "gene116", "gene137", "gene108",
    "gene191", "gene160"
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
  names(cell_clusters) <- colnames(res@data)
  
  res_scigenex_top <- top_genes(res, top = 20)
  
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
  expect_equal(round(sum(htmp_matrix_wo_na), 2), 52.05)
  expect_equal(round(mean(htmp_matrix_wo_na), 2), 0.03)
  expect_equal(round(median(htmp_matrix_wo_na), 2), 0.15)
  expect_equal(round(sd(htmp_matrix_wo_na), 2), 0.88)
  
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["0%"]), 2), -1)
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["25%"]), 2), -1)
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["50%"]), 2), 0.15)
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
    "gene2", "gene70", "gene12", "gene29", "gene31", "gene80",
    "gene61", "gene55", "gene27", "gene57", "gene78", "gene7",
    "gene88", "gene76", "gene32", "gene64", "gene51", "gene19",
    "gene84", "gene37", " ", " ", " ", " ",
    " ", " ", " ", " ", " ", " ",
    "gene326", "gene363", "gene349", "gene369", "gene339", "gene397",
    "gene382", "gene387", "gene350", "gene316", "gene368", "gene393",
    "gene398", "gene334", "gene314", "gene320", "gene390", "gene312",
    "gene364", "gene315", " ", " ", " ", " ",
    " ", " ", " ", " ", " ", " ",
    "gene288", "gene232", "gene202", "gene212", "gene229", "gene269",
    "gene248", "gene219", "gene277", "gene285", "gene283", "gene282",
    "gene281", "gene256", "gene273", "gene207", "gene206", "gene223",
    "gene245", "gene286", " ", " ", " ", " ",
    " ", " ", " ", " ", " ", " ",
    "gene115", "gene119", "gene155", "gene163", "gene158", "gene198",
    "gene109", "gene122", "gene120", "gene126", "gene121", "gene176",
    "gene141", "gene143", "gene171", "gene116", "gene137", "gene108",
    "gene191", "gene160"
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
  names(cell_clusters) <- colnames(res@data)
  
  res_scigenex_top <- top_genes(res, top = 20)
  
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
  expect_equal(round(sum(htmp_matrix_wo_na), 2), -1.26)
  expect_equal(round(mean(htmp_matrix_wo_na), 2), 0)
  expect_equal(round(median(htmp_matrix_wo_na), 2), 0.06)
  expect_equal(round(sd(htmp_matrix_wo_na), 2), 0.91)
  
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["0%"]), 2), -1)
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["25%"]), 2), -1)
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["50%"]), 2), 0.06)
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["75%"]), 2), 1)
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["100%"]), 2), 1)
  
  # Checking order of columns
  expect_equal(colnames(htmp_matrix), c("sample1", "sample2", "sample3", "sample4", "sample5", "sample6", 
                                        "sample7", "sample8", "sample9", "sample10", " ", "  ", "sample11", 
                                        "sample12", "sample13", "sample14", "sample15", "sample16", "sample17", 
                                        "sample18", "sample19", "sample20"))
  
  # Checking order of rows
  expect_equal(rownames(htmp_matrix), c(
    "gene2", "gene70", "gene12", "gene29", "gene31", "gene80",
    "gene61", "gene55", "gene27", "gene57", "gene78", "gene7",
    "gene88", "gene76", "gene32", "gene64", "gene51", "gene19",
    "gene84", "gene37", " ", " ", "gene115", "gene119",
    "gene155", "gene163", "gene158", "gene198", "gene109", "gene122",
    "gene120", "gene126", "gene121", "gene176", "gene141", "gene143",
    "gene171", "gene116", "gene137", "gene108", "gene191", "gene160"
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
  names(cell_clusters) <- colnames(res@data)
  res_scigenex_top <- top_genes(res, top = 20)
  
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
  expect_equal(round(sum(htmp_matrix_wo_na), 2), -3.42)
  expect_equal(round(mean(htmp_matrix_wo_na), 2), -0.01)
  expect_equal(round(median(htmp_matrix_wo_na), 2), -0.06)
  expect_equal(round(sd(htmp_matrix_wo_na), 2), 0.86)
  
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["0%"]), 2), -1)
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["25%"]), 2), -1)
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["50%"]), 2), -0.06)
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["75%"]), 2), 1)
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["100%"]), 2), 1)
  
  # Checking order of columns
  expect_equal(colnames(htmp_matrix), c(
    "sample1", "sample2", "sample3", "sample4", "sample5",
    "sample6", "sample7", "sample8", "sample9", "sample10", " ",
    "  ", "sample11", "sample12", "sample13", "sample14", "sample15",
    "sample16", "sample17", "sample18", "sample19", "sample20"
  ))
  
  # Checking order of rows
  expect_equal(rownames(htmp_matrix), c(
    "gene326", "gene363", "gene349", "gene369", "gene339", "gene397",
    "gene382", "gene387", "gene350", "gene316", "gene368", "gene393",
    "gene398", "gene334", "gene314", "gene320", "gene390", "gene312",
    "gene364", "gene315"
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




testthat::test_that("Check plot_heatmap, subsetting and clustering", {
  # Set verbosity to 0
  set_verbosity(0)
  
  # Load seurat object
  data(pbmc_small, package = "SeuratObject")
  
  ## Select informative genes
  clust_set <- select_genes(data=pbmc_small,
                            distance_method="pearson",
                            k=10,
                            row_sum=-Inf,
                            highest=0.95,
                            fdr = 1e-6)
  
  ## Cluster genes
  clust_set <- gene_clustering(object = clust_set,
                               inflation = 1.2,
                               keep_nn = FALSE,
                               k = 5,
                               threads = 1)
  
  ident <- Seurat::Idents(pbmc_small)
  
  # Check color arg
  expect_error(plot_heatmap(clust_set, cell_clusters=Seurat::Idents(pbmc_small), colors_cell_clusters = "green"))
  p <- plot_heatmap(clust_set[, names(ident[ident==2])], cell_clusters=Seurat::Idents(pbmc_small), colors_cell_clusters = "red")
  testthat::expect_error(print(p), NA)
  expect_equal(unique(p@colorbars@listData$Ident.@colors), "red")
  
  # Check if seurat identities is considered
  expect_error(plot_heatmap(clust_set[, names(ident[ident==3])], cell_clusters=Seurat::Idents(pbmc_small)))
  p <- plot_heatmap(clust_set, cell_clusters=Seurat::Idents(pbmc_small))
  testthat::expect_error(print(p), NA)
  expect_equal(p@colorbars@listData$Ident.@ticktext, c("1", "2", "3"))
  
  p <- plot_heatmap(clust_set[,colnames(clust_set@data)[1:80]], cell_clusters=Seurat::Idents(pbmc_small))
  testthat::expect_error(print(p), NA)
  expect_equal(dim(p@plots@listData$`Exp. level`@data), c(178, 86)) # 80 + 6 NA columns
  
  p <- plot_heatmap(clust_set[,colnames(clust_set@data)[1:10]], cell_clusters=Seurat::Idents(pbmc_small))
  testthat::expect_error(print(p), NA)
  expect_equal(dim(p@plots@listData$`Exp. level`@data), c(178, 10))
  
  p <- plot_heatmap(clust_set[,colnames(clust_set@data)[1:20]], cell_clusters=Seurat::Idents(pbmc_small))
  testthat::expect_error(print(p), NA)
  expect_equal(dim(p@plots@listData$`Exp. level`@data), c(178, 23)) # 20 + 3 NA columns
  
  p <- plot_heatmap(clust_set[1,colnames(clust_set@data)[1:10]], cell_clusters=Seurat::Idents(pbmc_small))
  testthat::expect_error(print(p), NA)
  expect_equal(dim(p@plots@listData$`Exp. level`@data), c(91, 10))
  expect_equal(dim(p@plots@listData$`Exp. level`@data), c(length(clust_set@gene_clusters$`1`), 10))
  
  p <- plot_heatmap(clust_set[2 ,colnames(clust_set@data)[1:10]], cell_clusters=Seurat::Idents(pbmc_small))
  testthat::expect_error(print(p), NA)
  expect_equal(dim(p@plots@listData$`Exp. level`@data), c(52, 10))
  expect_equal(dim(p@plots@listData$`Exp. level`@data), c(length(clust_set@gene_clusters$`2`), 10))
})




testthat::test_that("Test non interactive heatmaps.", {
  # Set verbosity to 0
  set_verbosity(0)
  
  # Load seurat object
  data(pbmc_small, package = "SeuratObject")
  
  ## Select informative genes
  clust_set <- select_genes(data=pbmc_small,
                            distance_method="pearson",
                            k=10,
                            row_sum=-Inf,
                            highest=0.95,
                            fdr = 1e-6)
  
  ## Cluster genes
  clust_set <- gene_clustering(object = clust_set,
                               inflation = 1.2,
                               keep_nn = FALSE,
                               k = 5,
                               threads = 1)
                               
  ident <- Seurat::Idents(pbmc_small)
  
  n_gene_clusters <- length(clust_set@gene_clusters) 
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
