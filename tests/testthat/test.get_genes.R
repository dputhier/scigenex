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
                    fdr = 1e-8)

## Cluster genes
res <- gene_clustering(object = res,
                       inflation = 1.2,
                       keep_nn = FALSE,
                       k = 5,
                       threads = 1)




test_that("Checking if get_genes stop if no slot top_genes in ClusterSet obj", {
  expect_error(get_genes(res, cluster = "all", top = TRUE))
})

test_that("Checking get_genes is providing the right list of genes", {
  
  # =======================================================
  # Test gene list in all cluster
  gene_names <- get_genes(res, cluster = "all")
  gene_name_to_check <- c(
    "gene126", "gene155", "gene189", "gene192", "gene176", "gene141",
    "gene150", "gene144", "gene146", "gene143", "gene117", "gene198",
    "gene105", "gene142", "gene188", "gene101", "gene119", "gene103",
    "gene186", "gene173", "gene180", "gene106", "gene107", "gene135",
    "gene159", "gene108", "gene122", "gene111", "gene145", "gene133",
    "gene112", "gene194", "gene153", "gene160", "gene116", "gene148",
    "gene166", "gene197", "gene127", "gene120", "gene137", "gene190",
    "gene121", "gene191", "gene171", "gene114", "gene178", "gene124",
    "gene118", "gene165", "gene162", "gene168", "gene136", "gene164",
    "gene132", "gene161", "gene139", "gene138", "gene151", "gene154",
    "gene158", "gene163", "gene115", "gene183", "gene193", "gene109",
    "gene199", "gene195", "gene184", "gene104", "gene167", "gene182",
    "gene147", "gene175", "gene113", "gene129", "gene157", "gene131",
    "gene196", "gene140", "gene149", "gene125", "gene152", "gene177",
    "gene179", "gene128", "gene102", "gene130", "gene110", "gene123",
    "gene134", "gene185", "gene169", "gene172", "gene181", "gene200",
    "gene170", "gene174", "gene187", "gene156", "gene1922", "gene838",
    "gene1616", "gene1856", "gene1503", "gene3528", "gene2184", "gene2573",
    "gene2725", "gene3442", "gene3485", "gene3488", "gene1651", "gene1915",
    "gene2112", "gene2223", "gene2662", "gene2421", "gene2553", "gene2902",
    "gene3610", "gene3779", "gene3924", "gene270", "gene211", "gene222",
    "gene215", "gene279", "gene284", "gene202", "gene2930", "gene239",
    "gene250", "gene278", "gene226", "gene219", "gene230", "gene214",
    "gene244", "gene251", "gene206", "gene229", "gene272", "gene232",
    "gene286", "gene233", "gene291", "gene235", "gene225", "gene256",
    "gene237", "gene242", "gene208", "gene298", "gene253", "gene245",
    "gene273", "gene285", "gene217", "gene262", "gene290", "gene269",
    "gene258", "gene277", "gene248", "gene282", "gene267", "gene263",
    "gene223", "gene236", "gene300", "gene288", "gene207", "gene257",
    "gene289", "gene249", "gene227", "gene3991", "gene224", "gene203",
    "gene220", "gene212", "gene261", "gene238", "gene240", "gene283",
    "gene281", "gene241", "gene294", "gene266", "gene259", "gene296",
    "gene247", "gene201", "gene275", "gene255", "gene293", "gene297",
    "gene246", "gene228", "gene276", "gene254", "gene287", "gene299",
    "gene213", "gene280", "gene234", "gene252", "gene268", "gene292",
    "gene1931", "gene315", "gene334", "gene312", "gene326", "gene323",
    "gene318", "gene368", "gene363", "gene314", "gene391", "gene324",
    "gene320", "gene364", "gene385", "gene362", "gene322", "gene386",
    "gene319", "gene390", "gene398", "gene366", "gene379", "gene332",
    "gene393", "gene369", "gene311", "gene397", "gene358", "gene316",
    "gene387", "gene382", "gene329", "gene304", "gene352", "gene349",
    "gene303", "gene337", "gene373", "gene339", "gene306", "gene378",
    "gene340", "gene313", "gene336", "gene359", "gene330", "gene331",
    "gene389", "gene371", "gene350", "gene327", "gene365", "gene301",
    "gene394", "gene367", "gene341", "gene351", "gene353", "gene355",
    "gene335", "gene372", "gene400", "gene3609", "gene392", "gene317",
    "gene345", "gene344", "gene328", "gene376", "gene360", "gene370",
    "gene2911", "gene302", "gene310", "gene342", "gene325", "gene375",
    "gene396", "gene760", "gene1423", "gene2490", "gene19", "gene84",
    "gene27", "gene87", "gene58", "gene31", "gene37", "gene6",
    "gene86", "gene76", "gene61", "gene74", "gene90", "gene55",
    "gene21", "gene64", "gene39", "gene16", "gene88", "gene73",
    "gene2", "gene67", "gene57", "gene59", "gene12", "gene7",
    "gene51", "gene79", "gene17", "gene83", "gene36", "gene25",
    "gene28", "gene29", "gene34", "gene78", "gene45", "gene70",
    "gene94", "gene9", "gene66", "gene96", "gene32", "gene10",
    "gene80", "gene42", "gene1", "gene71", "gene4", "gene14",
    "gene13", "gene56", "gene46", "gene63", "gene11", "gene20",
    "gene82", "gene93", "gene35", "gene85", "gene54", "gene95",
    "gene3", "gene43", "gene99", "gene52", "gene98"
  )
  
  
  expect_equal(gene_names, gene_name_to_check)
  expect_equal(gene_names, rownames(res@data))
  expect_equal(gene_names, unlist(res@gene_clusters, use.names = FALSE))
  expect_that(gene_names, is_a("character"))
  
  
  # =======================================================
  # Test gene list in cluster 1
  gene_names <- get_genes(res, cluster = 1)
  gene_name_to_check <- c(
    "gene126", "gene155", "gene189", "gene192", "gene176", "gene141",
    "gene150", "gene144", "gene146", "gene143", "gene117", "gene198",
    "gene105", "gene142", "gene188", "gene101", "gene119", "gene103",
    "gene186", "gene173", "gene180", "gene106", "gene107", "gene135",
    "gene159", "gene108", "gene122", "gene111", "gene145", "gene133",
    "gene112", "gene194", "gene153", "gene160", "gene116", "gene148",
    "gene166", "gene197", "gene127", "gene120", "gene137", "gene190",
    "gene121", "gene191", "gene171", "gene114", "gene178", "gene124",
    "gene118", "gene165", "gene162", "gene168", "gene136", "gene164",
    "gene132", "gene161", "gene139", "gene138", "gene151", "gene154",
    "gene158", "gene163", "gene115", "gene183", "gene193", "gene109",
    "gene199", "gene195", "gene184", "gene104", "gene167", "gene182",
    "gene147", "gene175", "gene113", "gene129", "gene157", "gene131",
    "gene196", "gene140", "gene149", "gene125", "gene152", "gene177",
    "gene179", "gene128", "gene102", "gene130", "gene110", "gene123",
    "gene134", "gene185", "gene169", "gene172", "gene181", "gene200",
    "gene170", "gene174", "gene187", "gene156", "gene1922", "gene838",
    "gene1616", "gene1856", "gene1503", "gene3528", "gene2184", "gene2573",
    "gene2725", "gene3442", "gene3485", "gene3488", "gene1651", "gene1915",
    "gene2112", "gene2223", "gene2662", "gene2421", "gene2553", "gene2902",
    "gene3610", "gene3779", "gene3924"
  )
  expect_equal(gene_names, gene_name_to_check)
  expect_equal(gene_names, res@gene_clusters$`1`)
  expect_that(gene_names, is_a("character"))
  
  
  # =======================================================
  # Test gene list in cluster 2
  gene_names <- get_genes(res, cluster = 2)
  gene_name_to_check <- c(
    "gene270", "gene211", "gene222", "gene215", "gene279", "gene284",
    "gene202", "gene2930", "gene239", "gene250", "gene278", "gene226",
    "gene219", "gene230", "gene214", "gene244", "gene251", "gene206",
    "gene229", "gene272", "gene232", "gene286", "gene233", "gene291",
    "gene235", "gene225", "gene256", "gene237", "gene242", "gene208",
    "gene298", "gene253", "gene245", "gene273", "gene285", "gene217",
    "gene262", "gene290", "gene269", "gene258", "gene277", "gene248",
    "gene282", "gene267", "gene263", "gene223", "gene236", "gene300",
    "gene288", "gene207", "gene257", "gene289", "gene249", "gene227",
    "gene3991", "gene224", "gene203", "gene220", "gene212", "gene261",
    "gene238", "gene240", "gene283", "gene281", "gene241", "gene294",
    "gene266", "gene259", "gene296", "gene247", "gene201", "gene275",
    "gene255", "gene293", "gene297", "gene246", "gene228", "gene276",
    "gene254", "gene287", "gene299", "gene213", "gene280", "gene234",
    "gene252", "gene268", "gene292", "gene1931"
  )
  expect_equal(gene_names, gene_name_to_check)
  expect_equal(gene_names, res@gene_clusters$`2`)
  expect_that(gene_names, is_a("character"))
  
  
  # =======================================================
  # Test gene list in cluster 3
  gene_names <- get_genes(res, cluster = 3)
  gene_name_to_check <- c(
    "gene315", "gene334", "gene312", "gene326", "gene323", "gene318",
    "gene368", "gene363", "gene314", "gene391", "gene324", "gene320",
    "gene364", "gene385", "gene362", "gene322", "gene386", "gene319",
    "gene390", "gene398", "gene366", "gene379", "gene332", "gene393",
    "gene369", "gene311", "gene397", "gene358", "gene316", "gene387",
    "gene382", "gene329", "gene304", "gene352", "gene349", "gene303",
    "gene337", "gene373", "gene339", "gene306", "gene378", "gene340",
    "gene313", "gene336", "gene359", "gene330", "gene331", "gene389",
    "gene371", "gene350", "gene327", "gene365", "gene301", "gene394",
    "gene367", "gene341", "gene351", "gene353", "gene355", "gene335",
    "gene372", "gene400", "gene3609", "gene392", "gene317", "gene345",
    "gene344", "gene328", "gene376", "gene360", "gene370", "gene2911",
    "gene302", "gene310", "gene342", "gene325", "gene375", "gene396",
    "gene760", "gene1423", "gene2490"
  )
  expect_equal(gene_names, gene_name_to_check)
  expect_equal(gene_names, res@gene_clusters$`3`)
  expect_that(gene_names, is_a("character"))
  
  
  # =======================================================
  # Test gene list in cluster 4
  gene_names <- get_genes(res, cluster = 4)
  gene_name_to_check <- c(
    "gene19", "gene84", "gene27", "gene87", "gene58", "gene31",
    "gene37", "gene6", "gene86", "gene76", "gene61", "gene74",
    "gene90", "gene55", "gene21", "gene64", "gene39", "gene16",
    "gene88", "gene73", "gene2", "gene67", "gene57", "gene59",
    "gene12", "gene7", "gene51", "gene79", "gene17", "gene83",
    "gene36", "gene25", "gene28", "gene29", "gene34", "gene78",
    "gene45", "gene70", "gene94", "gene9", "gene66", "gene96",
    "gene32", "gene10", "gene80", "gene42", "gene1", "gene71",
    "gene4", "gene14", "gene13", "gene56", "gene46", "gene63",
    "gene11", "gene20", "gene82", "gene93", "gene35", "gene85",
    "gene54", "gene95", "gene3", "gene43", "gene99", "gene52",
    "gene98"
  )
  expect_equal(gene_names, gene_name_to_check)
  expect_equal(gene_names, res@gene_clusters$`4`)
  expect_that(gene_names, is_a("character"))
  
  
  # =======================================================
  # Test gene list in cluster 3 and 4
  gene_names <- get_genes(res, cluster = c(3, 4))
  gene_name_to_check <- c(
    "gene315", "gene334", "gene312", "gene326", "gene323", "gene318",
    "gene368", "gene363", "gene314", "gene391", "gene324", "gene320",
    "gene364", "gene385", "gene362", "gene322", "gene386", "gene319",
    "gene390", "gene398", "gene366", "gene379", "gene332", "gene393",
    "gene369", "gene311", "gene397", "gene358", "gene316", "gene387",
    "gene382", "gene329", "gene304", "gene352", "gene349", "gene303",
    "gene337", "gene373", "gene339", "gene306", "gene378", "gene340",
    "gene313", "gene336", "gene359", "gene330", "gene331", "gene389",
    "gene371", "gene350", "gene327", "gene365", "gene301", "gene394",
    "gene367", "gene341", "gene351", "gene353", "gene355", "gene335",
    "gene372", "gene400", "gene3609", "gene392", "gene317", "gene345",
    "gene344", "gene328", "gene376", "gene360", "gene370", "gene2911",
    "gene302", "gene310", "gene342", "gene325", "gene375", "gene396",
    "gene760", "gene1423", "gene2490", "gene19", "gene84", "gene27",
    "gene87", "gene58", "gene31", "gene37", "gene6", "gene86",
    "gene76", "gene61", "gene74", "gene90", "gene55", "gene21",
    "gene64", "gene39", "gene16", "gene88", "gene73", "gene2",
    "gene67", "gene57", "gene59", "gene12", "gene7", "gene51",
    "gene79", "gene17", "gene83", "gene36", "gene25", "gene28",
    "gene29", "gene34", "gene78", "gene45", "gene70", "gene94",
    "gene9", "gene66", "gene96", "gene32", "gene10", "gene80",
    "gene42", "gene1", "gene71", "gene4", "gene14", "gene13",
    "gene56", "gene46", "gene63", "gene11", "gene20", "gene82",
    "gene93", "gene35", "gene85", "gene54", "gene95", "gene3",
    "gene43", "gene99", "gene52", "gene98"
  )
  expect_equal(gene_names, gene_name_to_check)
  expect_equal(gene_names, c(res@gene_clusters$`3`, res@gene_clusters$`4`))
  expect_that(gene_names, is_a("character"))
})








test_that("Checking get_genes using top genes.", {
  res <- top_genes(res, top = 20, cluster = "all")
  
  # =======================================================
  # Test top gene list in all cluster
  gene_names <- get_genes(res, cluster = "all", top = T)
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
  expect_equal(gene_names, gene_name_to_check)
  expect_equal(sort(gene_names), sort(unlist(res@top_genes, use.names = F)))
  expect_that(gene_names, is_a("character"))
  
  
  # =======================================================
  # Test top gene list in cluster 1
  gene_names <- get_genes(res, cluster = 1, top = T)
  gene_name_to_check <- c(
    "gene160", "gene191", "gene108", "gene137", "gene116", "gene171",
    "gene143", "gene141", "gene176", "gene121", "gene126", "gene120",
    "gene122", "gene109", "gene198", "gene158", "gene163", "gene155",
    "gene119", "gene115"
  )
  expect_equal(gene_names, gene_name_to_check)
  expect_equal(sort(gene_names), sort(res@top_genes$`1`))
  expect_that(gene_names, is_a("character"))
  
  
  # =======================================================
  # Test top gene list in cluster 2 and 3
  gene_names <- get_genes(res, cluster = 2:3, top = T)
  gene_name_to_check <- c(
    "gene286", "gene245", "gene223", "gene206", "gene207", "gene273",
    "gene256", "gene281", "gene282", "gene283", "gene285", "gene277",
    "gene219", "gene248", "gene269", "gene229", "gene212", "gene202",
    "gene232", "gene288", "gene315", "gene364", "gene312", "gene390",
    "gene320", "gene314", "gene334", "gene398", "gene393", "gene368",
    "gene316", "gene350", "gene387", "gene382", "gene397", "gene339",
    "gene369", "gene349", "gene363", "gene326"
  )
  expect_equal(gene_names, gene_name_to_check)
  expect_equal(sort(gene_names), sort(c(
    res@top_genes$`2`,
    res@top_genes$`3`
  )))
  expect_that(gene_names, is_a("character"))
})
