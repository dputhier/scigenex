# Set verbosity to 0
set_verbosity(0)

#Create matrix containing 3 signatures
m <- create_4_rnd_clust()




test_that("Cheking get_genes is providing the right list of genes - pearson", {
  ## Select informative genes
  res <- select_genes(data=m,
                      distance_method="kendall",
                      k=75,
                      row_sum=-Inf,
                      highest=0.3,
                      fdr = 1e-8)
  
  ## Cluster genes
  res <- gene_clustering(object = res,
                         inflation = 1.2,
                         keep_nn = FALSE,
                         k = 5,
                         threads = 1)
  
  # ========================================
  # Top 20
  res_20 <- suppressWarnings(top_genes(res, cluster = "all", top = 20))
  
  # Test top genes list in all cluster
  gene_name_to_check <- c(
    "gene160", "gene191", "gene108", "gene137", "gene116", "gene171",
    "gene143", "gene141", "gene176", "gene121", "gene126", "gene120",
    "gene122", "gene109", "gene198", "gene158", "gene163", "gene155",
    "gene119", "gene115", "gene286", "gene245", "gene223", "gene206",
    "gene207", "gene273", "gene256", "gene281", "gene282", "gene283",
    "gene285", "gene277", "gene219", "gene248", "gene269", "gene229",
    "gene212", "gene202", "gene232", "gene288", "gene315", "gene364",
    "gene312", "gene390", "gene320", "gene314", "gene334", "gene398",
    "gene393", "gene368", "gene316", "gene350", "gene387", "gene382",
    "gene397", "gene339", "gene369", "gene349", "gene363", "gene326",
    "gene37", "gene84", "gene19", "gene51", "gene64", "gene32",
    "gene76", "gene88", "gene7", "gene78", "gene57", "gene27",
    "gene55", "gene61", "gene80", "gene31", "gene29", "gene12",
    "gene70", "gene2"
  )
  expect_equal(unlist(res_20@top_genes, use.names = FALSE), gene_name_to_check)
  expect_that(res@top_genes, is_a("list"))
  
  
  # Test top genes list in cluster 1
  gene_name_to_check <- c(
    "gene160", "gene191", "gene108", "gene137", "gene116", "gene171",
    "gene143", "gene141", "gene176", "gene121", "gene126", "gene120",
    "gene122", "gene109", "gene198", "gene158", "gene163", "gene155",
    "gene119", "gene115"
  )
  expect_equal(res_20@top_genes$`1`, gene_name_to_check)
  expect_equal(length(res_20@top_genes$`1`), 20)
  
  
  # Test top genes list in cluster 2 and 3
  gene_name_to_check <- c(
    "gene286", "gene245", "gene223", "gene206", "gene207", "gene273",
    "gene256", "gene281", "gene282", "gene283", "gene285", "gene277",
    "gene219", "gene248", "gene269", "gene229", "gene212", "gene202",
    "gene232", "gene288", "gene315", "gene364", "gene312", "gene390",
    "gene320", "gene314", "gene334", "gene398", "gene393", "gene368",
    "gene316", "gene350", "gene387", "gene382", "gene397", "gene339",
    "gene369", "gene349", "gene363", "gene326"
  )
  expect_equal(c(res_20@top_genes$`2`, res_20@top_genes$`3`),
               gene_name_to_check)
  expect_equal(length(c(res_20@top_genes$`2`, res_20@top_genes$`3`)), 40)
  
  
  # ========================================
  # Top 10
  res_10 <- top_genes(res, cluster = "all", top = 10)
  
  # Test top genes list in all cluster
  gene_name_to_check <- c(
    "gene160", "gene191", "gene108", "gene137", "gene116", "gene171",
    "gene143", "gene141", "gene176", "gene121", "gene286", "gene245",
    "gene223", "gene206", "gene207", "gene273", "gene256", "gene281",
    "gene282", "gene283", "gene315", "gene364", "gene312", "gene390",
    "gene320", "gene314", "gene334", "gene398", "gene393", "gene368",
    "gene37", "gene84", "gene19", "gene51", "gene64", "gene32",
    "gene76", "gene88", "gene7", "gene78"
  )
  expect_equal(unlist(res_10@top_genes, use.names = FALSE), gene_name_to_check)
  expect_that(res@top_genes, is_a("list"))
  expect_equal(length(unlist(res_10@top_genes, use.names = FALSE)), 40)
  
  
  # Test top genes list in cluster 1
  gene_name_to_check <- c(
    "gene160", "gene191", "gene108", "gene137", "gene116", "gene171",
    "gene143", "gene141", "gene176", "gene121"
  )
  expect_equal(res_10@top_genes$`1`, gene_name_to_check)
  expect_equal(length(res_10@top_genes$`1`), 10)
  
  
  # Test top genes list in cluster 2 and 3
  gene_name_to_check <- c(
    "gene286", "gene245", "gene223", "gene206", "gene207", "gene273",
    "gene256", "gene281", "gene282", "gene283", "gene315", "gene364",
    "gene312", "gene390", "gene320", "gene314", "gene334", "gene398",
    "gene393", "gene368"
  )
  expect_equal(
    c(res_10@top_genes$`2`, res_10@top_genes$`3`),
    gene_name_to_check
  )
  expect_equal(length(c(res_10@top_genes$`2`, res_10@top_genes$`3`)), 20)
  
  
  
  # ========================================
  # Top 5
  res_5 <- top_genes(res, cluster = "all", top = 5)
  
  # Test top genes list in all cluster
  gene_name_to_check <- c(
    "gene160", "gene191", "gene108", "gene137", "gene116", "gene286",
    "gene245", "gene223", "gene206", "gene207", "gene315", "gene364",
    "gene312", "gene390", "gene320", "gene37", "gene84", "gene19",
    "gene51", "gene64"
  )
  expect_equal(unlist(res_5@top_genes, use.names = FALSE), gene_name_to_check)
  expect_that(res@top_genes, is_a("list"))
  expect_equal(length(unlist(res_5@top_genes, use.names = FALSE)), 20)
  
  
  # Test top genes list in cluster 1
  gene_name_to_check <- c(
    "gene160", "gene191", "gene108", "gene137", "gene116"
  )
  expect_equal(res_5@top_genes$`1`, gene_name_to_check)
  
  
  # Test top genes list in cluster 2 and 3
  gene_name_to_check <- c(
    "gene286", "gene245", "gene223", "gene206", "gene207", "gene315",
    "gene364", "gene312", "gene390", "gene320"
  )
  expect_equal(c(res_5@top_genes$`2`, res_5@top_genes$`3`), gene_name_to_check)
  
  
  
  # ========================================
  # Top 100
  res_100 <- suppressWarnings(top_genes(res, cluster = "all", top = 100))
  
  # Test top genes list in all cluster
  gene_name_to_check <- c(
    "gene160", "gene191", "gene108", "gene137", "gene116", "gene171",
    "gene143", "gene141", "gene176", "gene121", "gene126", "gene120",
    "gene122", "gene109", "gene198", "gene158", "gene163", "gene155",
    "gene119", "gene115", "gene148", "gene154", "gene197", "gene124",
    "gene188", "gene159", "gene161", "gene112", "gene150", "gene144",
    "gene192", "gene118", "gene184", "gene146", "gene175", "gene147",
    "gene199", "gene173", "gene114", "gene103", "gene127", "gene183",
    "gene194", "gene113", "gene101", "gene181", "gene180", "gene153",
    "gene165", "gene135", "gene111", "gene166", "gene142", "gene186",
    "gene131", "gene104", "gene136", "gene168", "gene134", "gene152",
    "gene149", "gene179", "gene151", "gene178", "gene139", "gene133",
    "gene162", "gene138", "gene107", "gene189", "gene200", "gene140",
    "gene132", "gene105", "gene193", "gene195", "gene145", "gene123",
    "gene185", "gene164", "gene177", "gene157", "gene106", "gene190",
    "gene110", "gene174", "gene102", "gene169", "gene117", "gene167",
    "gene125", "gene172", "gene129", "gene130", "gene196", "gene182",
    "gene128", "gene187", "gene170", "gene156", "gene286", "gene245",
    "gene223", "gene206", "gene207", "gene273", "gene256", "gene281",
    "gene282", "gene283", "gene285", "gene277", "gene219", "gene248",
    "gene269", "gene229", "gene212", "gene202", "gene232", "gene288",
    "gene275", "gene236", "gene289", "gene230", "gene292", "gene203",
    "gene244", "gene276", "gene299", "gene272", "gene241", "gene261",
    "gene222", "gene259", "gene253", "gene258", "gene217", "gene279",
    "gene290", "gene270", "gene284", "gene278", "gene296", "gene233",
    "gene225", "gene267", "gene235", "gene238", "gene255", "gene287",
    "gene300", "gene266", "gene234", "gene239", "gene293", "gene215",
    "gene251", "gene291", "gene262", "gene257", "gene226", "gene242",
    "gene254", "gene214", "gene240", "gene250", "gene298", "gene249",
    "gene247", "gene263", "gene227", "gene220", "gene208", "gene3991",
    "gene294", "gene201", "gene237", "gene213", "gene211", "gene224",
    "gene228", "gene2930", "gene297", "gene268", "gene1931", "gene252",
    "gene280", "gene246", "gene315", "gene364", "gene312", "gene390",
    "gene320", "gene314", "gene334", "gene398", "gene393", "gene368",
    "gene316", "gene350", "gene387", "gene382", "gene397", "gene339",
    "gene369", "gene349", "gene363", "gene326", "gene318", "gene386",
    "gene341", "gene366", "gene324", "gene306", "gene378", "gene362",
    "gene329", "gene331", "gene303", "gene365", "gene340", "gene322",
    "gene353", "gene379", "gene330", "gene313", "gene335", "gene327",
    "gene345", "gene394", "gene311", "gene389", "gene360", "gene337",
    "gene301", "gene351", "gene372", "gene385", "gene336", "gene358",
    "gene400", "gene332", "gene323", "gene355", "gene2911", "gene342",
    "gene367", "gene319", "gene373", "gene396", "gene371", "gene310",
    "gene352", "gene359", "gene376", "gene392", "gene317", "gene370",
    "gene391", "gene328", "gene3609", "gene375", "gene760", "gene344",
    "gene304", "gene325", "gene1423", "gene302", "gene2490", "gene37",
    "gene84", "gene19", "gene51", "gene64", "gene32", "gene76",
    "gene88", "gene7", "gene78", "gene57", "gene27", "gene55",
    "gene61", "gene80", "gene31", "gene29", "gene12", "gene70",
    "gene2", "gene59", "gene56", "gene36", "gene67", "gene79",
    "gene16", "gene34", "gene13", "gene45", "gene1", "gene73",
    "gene93", "gene28", "gene39", "gene4", "gene87", "gene86",
    "gene94", "gene9", "gene6", "gene11", "gene25", "gene82",
    "gene10", "gene83", "gene71", "gene74", "gene96", "gene58",
    "gene46", "gene42", "gene63", "gene85", "gene95", "gene43",
    "gene20", "gene66", "gene52", "gene90", "gene14", "gene17",
    "gene35", "gene54", "gene21", "gene98", "gene99", "gene3"
  )
  expect_equal(unlist(res_100@top_genes, use.names = FALSE), gene_name_to_check)
  expect_that(res@top_genes, is_a("list"))
  
  
  # Test top genes list in cluster 1
  gene_name_to_check <- c(
    "gene160", "gene191", "gene108", "gene137", "gene116", "gene171",
    "gene143", "gene141", "gene176", "gene121", "gene126", "gene120",
    "gene122", "gene109", "gene198", "gene158", "gene163", "gene155",
    "gene119", "gene115", "gene148", "gene154", "gene197", "gene124",
    "gene188", "gene159", "gene161", "gene112", "gene150", "gene144",
    "gene192", "gene118", "gene184", "gene146", "gene175", "gene147",
    "gene199", "gene173", "gene114", "gene103", "gene127", "gene183",
    "gene194", "gene113", "gene101", "gene181", "gene180", "gene153",
    "gene165", "gene135", "gene111", "gene166", "gene142", "gene186",
    "gene131", "gene104", "gene136", "gene168", "gene134", "gene152",
    "gene149", "gene179", "gene151", "gene178", "gene139", "gene133",
    "gene162", "gene138", "gene107", "gene189", "gene200", "gene140",
    "gene132", "gene105", "gene193", "gene195", "gene145", "gene123",
    "gene185", "gene164", "gene177", "gene157", "gene106", "gene190",
    "gene110", "gene174", "gene102", "gene169", "gene117", "gene167",
    "gene125", "gene172", "gene129", "gene130", "gene196", "gene182",
    "gene128", "gene187", "gene170", "gene156"
  )
  expect_equal(res_100@top_genes$`1`, gene_name_to_check)
  
  
  # Test top genes list in cluster 2 and 3
  gene_name_to_check <- c(
    "gene286", "gene245", "gene223", "gene206", "gene207", "gene273",
    "gene256", "gene281", "gene282", "gene283", "gene285", "gene277",
    "gene219", "gene248", "gene269", "gene229", "gene212", "gene202",
    "gene232", "gene288", "gene275", "gene236", "gene289", "gene230",
    "gene292", "gene203", "gene244", "gene276", "gene299", "gene272",
    "gene241", "gene261", "gene222", "gene259", "gene253", "gene258",
    "gene217", "gene279", "gene290", "gene270", "gene284", "gene278",
    "gene296", "gene233", "gene225", "gene267", "gene235", "gene238",
    "gene255", "gene287", "gene300", "gene266", "gene234", "gene239",
    "gene293", "gene215", "gene251", "gene291", "gene262", "gene257",
    "gene226", "gene242", "gene254", "gene214", "gene240", "gene250",
    "gene298", "gene249", "gene247", "gene263", "gene227", "gene220",
    "gene208", "gene3991", "gene294", "gene201", "gene237", "gene213",
    "gene211", "gene224", "gene228", "gene2930", "gene297", "gene268",
    "gene1931", "gene252", "gene280", "gene246", "gene315", "gene364",
    "gene312", "gene390", "gene320", "gene314", "gene334", "gene398",
    "gene393", "gene368", "gene316", "gene350", "gene387", "gene382",
    "gene397", "gene339", "gene369", "gene349", "gene363", "gene326",
    "gene318", "gene386", "gene341", "gene366", "gene324", "gene306",
    "gene378", "gene362", "gene329", "gene331", "gene303", "gene365",
    "gene340", "gene322", "gene353", "gene379", "gene330", "gene313",
    "gene335", "gene327", "gene345", "gene394", "gene311", "gene389",
    "gene360", "gene337", "gene301", "gene351", "gene372", "gene385",
    "gene336", "gene358", "gene400", "gene332", "gene323", "gene355",
    "gene2911", "gene342", "gene367", "gene319", "gene373", "gene396",
    "gene371", "gene310", "gene352", "gene359", "gene376", "gene392",
    "gene317", "gene370", "gene391", "gene328", "gene3609", "gene375",
    "gene760", "gene344", "gene304", "gene325", "gene1423", "gene302",
    "gene2490"
  )
  expect_equal(
    c(res_100@top_genes$`2`, res_100@top_genes$`3`),
    gene_name_to_check
  )
})





test_that("Cheking get_genes is providing the right list of genes - cosine", {
  
  ## Select informative genes
  res <- select_genes(data=m,
                      distance_method="cosine",
                      k=75,
                      row_sum=-Inf,
                      highest=0.3,
                      fdr = 1e-8)
  
  ## Cluster genes
  res <- gene_clustering(object = res,
                         inflation = 1.2,
                         keep_nn = FALSE,
                         k = 5,
                         threads = 1)
  
  # ========================================
  # Top 20
  res_20 <- suppressWarnings(top_genes(res, cluster = "all", top = 20))
  
  # Test top genes list in all cluster
  gene_name_to_check <- c(
    "gene189", "gene155", "gene192", "gene165", "gene190", "gene150",
    "gene146", "gene180", "gene186", "gene166", "gene132", "gene159",
    "gene112", "gene153", "gene143", "gene122", "gene176", "gene161",
    "gene184", "gene117", "gene316", "gene393", "gene386", "gene379",
    "gene330", "gene398", "gene368", "gene382", "gene315", "gene390",
    "gene334", "gene364", "gene312", "gene329", "gene349", "gene389",
    "gene397", "gene318", "gene365", "gene314", "gene230", "gene286",
    "gene285", "gene249", "gene256", "gene272", "gene245", "gene282",
    "gene223", "gene214", "gene222", "gene289", "gene298", "gene267",
    "gene293", "gene241", "gene206", "gene259", "gene225", "gene251",
    "gene76", "gene67", "gene16", "gene59", "gene88", "gene12",
    "gene32", "gene55", "gene80", "gene73", "gene9", "gene39",
    "gene19", "gene2", "gene56", "gene36", "gene74", "gene84",
    "gene6", "gene29"
  )
  expect_equal(unlist(res_20@top_genes, use.names = FALSE), gene_name_to_check)
  expect_that(res@top_genes, is_a("list"))
})



test_that("Cheking get_genes is providing the right list of genes - euclidean", {
  
  ## Select informative genes
  res <- select_genes(data=m,
                      distance_method="euclidean",
                      k=75,
                      row_sum=-Inf,
                      highest=0.3,
                      fdr = 1e-8)
  
  ## Cluster genes
  res <- gene_clustering(object = res,
                         inflation = 1.2,
                         keep_nn = FALSE,
                         k = 5,
                         threads = 1)
  
  # ========================================
  # Top 20
  res_20 <- suppressWarnings(top_genes(res, cluster = "all", top = 20))
  
  # Test top genes list in all cluster
  gene_name_to_check <- c(
    "gene568", "gene1938", "gene717", "gene2616", "gene2037", "gene3169",
    "gene1430", "gene464", "gene864", "gene3478", "gene3224", "gene1015",
    "gene1873", "gene793", "gene2041", "gene3122", "gene490", "gene824",
    "gene1225", "gene796"
  )
  expect_equal(unlist(res_20@top_genes, use.names = FALSE), gene_name_to_check)
  expect_that(res@top_genes, is_a("list"))
})




test_that("Cheking get_genes is providing the right list of genes - kendall", {
  
  ## Select informative genes
  res <- select_genes(data=m,
                      distance_method="kendall",
                      k=75,
                      row_sum=-Inf,
                      highest=0.3,
                      fdr = 1e-8)
  
  ## Cluster genes
  res <- gene_clustering(object = res,
                         inflation = 1.2,
                         keep_nn = FALSE,
                         k = 5,
                         threads = 1)
  
  # ========================================
  # Top 20
  res_20 <- suppressWarnings(top_genes(res, cluster = "all", top = 20))
  
  # Test top genes list in all cluster
  gene_name_to_check <- c(
    "gene160", "gene191", "gene108", "gene137", "gene116", "gene171",
    "gene143", "gene141", "gene176", "gene121", "gene126", "gene120",
    "gene122", "gene109", "gene198", "gene158", "gene163", "gene155",
    "gene119", "gene115", "gene286", "gene245", "gene223", "gene206",
    "gene207", "gene273", "gene256", "gene281", "gene282", "gene283",
    "gene285", "gene277", "gene219", "gene248", "gene269", "gene229",
    "gene212", "gene202", "gene232", "gene288", "gene315", "gene364",
    "gene312", "gene390", "gene320", "gene314", "gene334", "gene398",
    "gene393", "gene368", "gene316", "gene350", "gene387", "gene382",
    "gene397", "gene339", "gene369", "gene349", "gene363", "gene326",
    "gene37", "gene84", "gene19", "gene51", "gene64", "gene32",
    "gene76", "gene88", "gene7", "gene78", "gene57", "gene27",
    "gene55", "gene61", "gene80", "gene31", "gene29", "gene12",
    "gene70", "gene2"
  )
  expect_equal(unlist(res_20@top_genes, use.names = FALSE), gene_name_to_check)
  expect_that(res@top_genes, is_a("list"))
})
