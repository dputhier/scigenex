set_verbosity(0)

# Create matrix containing 4 signatures
m <- create_4_rnd_clust()

## A rather stringent version
res <- find_gene_clusters(
  data = m,
  distance_method = "pearson",
  inflation = 2,
  k = 75,
  row_sum = -Inf,
  highest = 0.3,
  min_nb_supporting_cell = 0,
  fdr = 1e-8
)

test_that("Checking if get_genes stop if no slot top_genes in ClusterSet obj", {
  expect_error(get_genes(res, cluster = "all", top = TRUE))
})

test_that("Checking get_genes is providing the right list of genes", {

  # =======================================================
  # Test gene list in all cluster
  gene_names <- get_genes(res, cluster = "all")
  gene_name_to_check <- c(
    "gene273", "gene285", "gene245", "gene272", "gene249", "gene230",
    "gene256", "gene225", "gene289", "gene286", "gene278", "gene282",
    "gene275", "gene214", "gene248", "gene277", "gene241", "gene206",
    "gene251", "gene298", "gene236", "gene267", "gene222", "gene290",
    "gene299", "gene223", "gene266", "gene217", "gene293", "gene229",
    "gene283", "gene300", "gene263", "gene235", "gene250", "gene212",
    "gene259", "gene219", "gene287", "gene233", "gene262", "gene239",
    "gene255", "gene258", "gene203", "gene207", "gene240", "gene216",
    "gene232", "gene215", "gene253", "gene269", "gene237", "gene294",
    "gene284", "gene247", "gene242", "gene244", "gene297", "gene281",
    "gene226", "gene279", "gene254", "gene246", "gene270", "gene292",
    "gene268", "gene288", "gene238", "gene257", "gene227", "gene220",
    "gene234", "gene202", "gene291", "gene208", "gene205", "gene221",
    "gene295", "gene271", "gene260", "gene224", "gene201", "gene296",
    "gene211", "gene243", "gene280", "gene261", "gene213", "gene209",
    "gene3991", "gene228", "gene218", "gene210", "gene265", "gene264",
    "gene274", "gene252", "gene1931", "gene231", "gene276", "gene847",
    "gene2036", "gene2157", "gene189", "gene146", "gene150", "gene117",
    "gene190", "gene135", "gene165", "gene186", "gene180", "gene155",
    "gene144", "gene132", "gene106", "gene192", "gene111", "gene153",
    "gene108", "gene159", "gene133", "gene166", "gene112", "gene143",
    "gene126", "gene145", "gene184", "gene122", "gene161", "gene176",
    "gene193", "gene137", "gene168", "gene160", "gene147", "gene107",
    "gene129", "gene188", "gene101", "gene104", "gene109", "gene114",
    "gene127", "gene157", "gene173", "gene197", "gene191", "gene177",
    "gene103", "gene105", "gene194", "gene163", "gene121", "gene116",
    "gene162", "gene158", "gene178", "gene175", "gene124", "gene141",
    "gene171", "gene128", "gene142", "gene198", "gene140", "gene119",
    "gene183", "gene167", "gene200", "gene199", "gene148", "gene174",
    "gene113", "gene151", "gene136", "gene179", "gene138", "gene181",
    "gene154", "gene130", "gene118", "gene115", "gene182", "gene134",
    "gene195", "gene164", "gene185", "gene120", "gene152", "gene149",
    "gene172", "gene131", "gene196", "gene110", "gene170", "gene125",
    "gene156", "gene123", "gene169", "gene102", "gene139", "gene393",
    "gene390", "gene312", "gene314", "gene316", "gene386", "gene364",
    "gene369", "gene318", "gene315", "gene398", "gene387", "gene363",
    "gene334", "gene330", "gene320", "gene382", "gene326", "gene379",
    "gene322", "gene329", "gene324", "gene349", "gene313", "gene306",
    "gene331", "gene368", "gene350", "gene327", "gene389", "gene319",
    "gene353", "gene397", "gene366", "gene362", "gene335", "gene365",
    "gene340", "gene378", "gene336", "gene358", "gene394", "gene311",
    "gene332", "gene303", "gene341", "gene301", "gene385", "gene351",
    "gene396", "gene370", "gene339", "gene355", "gene372", "gene342",
    "gene371", "gene345", "gene367", "gene337", "gene373", "gene317",
    "gene328", "gene391", "gene392", "gene3609", "gene344", "gene383",
    "gene323", "gene310", "gene384", "gene400", "gene920", "gene360",
    "gene304", "gene2911", "gene395", "gene27", "gene76", "gene37",
    "gene32", "gene70", "gene55", "gene67", "gene80", "gene59", "gene16",
    "gene12", "gene84", "gene45", "gene88", "gene19", "gene61", "gene78",
    "gene73", "gene36", "gene29", "gene39", "gene7", "gene2", "gene57",
    "gene86", "gene9", "gene52", "gene51", "gene66", "gene64", "gene13",
    "gene4", "gene1", "gene56", "gene74", "gene34", "gene79", "gene93",
    "gene42", "gene94", "gene82", "gene83", "gene10", "gene63", "gene25",
    "gene31", "gene17", "gene87", "gene71", "gene28", "gene46", "gene68",
    "gene11", "gene90", "gene35", "gene14", "gene58", "gene95", "gene85",
    "gene96", "gene99", "gene22", "gene20", "gene89"
  )


  expect_equal(gene_names, gene_name_to_check)
  expect_equal(gene_names, rownames(res@data))
  expect_equal(gene_names, unlist(res@gene_clusters, use.names = FALSE))
  expect_that(gene_names, is_a("character"))


  # =======================================================
  # Test gene list in cluster 1
  gene_names <- get_genes(res, cluster = 1)
  gene_name_to_check <- c(
    "gene273", "gene285", "gene245", "gene272", "gene249", "gene230",
    "gene256", "gene225", "gene289", "gene286", "gene278", "gene282",
    "gene275", "gene214", "gene248", "gene277", "gene241", "gene206",
    "gene251", "gene298", "gene236", "gene267", "gene222", "gene290",
    "gene299", "gene223", "gene266", "gene217", "gene293", "gene229",
    "gene283", "gene300", "gene263", "gene235", "gene250", "gene212",
    "gene259", "gene219", "gene287", "gene233", "gene262", "gene239",
    "gene255", "gene258", "gene203", "gene207", "gene240", "gene216",
    "gene232", "gene215", "gene253", "gene269", "gene237", "gene294",
    "gene284", "gene247", "gene242", "gene244", "gene297", "gene281",
    "gene226", "gene279", "gene254", "gene246", "gene270", "gene292",
    "gene268", "gene288", "gene238", "gene257", "gene227", "gene220",
    "gene234", "gene202", "gene291", "gene208", "gene205", "gene221",
    "gene295", "gene271", "gene260", "gene224", "gene201", "gene296",
    "gene211", "gene243", "gene280", "gene261", "gene213", "gene209",
    "gene3991", "gene228", "gene218", "gene210", "gene265", "gene264",
    "gene274", "gene252", "gene1931", "gene231", "gene276", "gene847",
    "gene2036", "gene2157"
  )
  expect_equal(gene_names, gene_name_to_check)
  expect_equal(gene_names, res@gene_clusters$`1`)
  expect_that(gene_names, is_a("character"))


  # =======================================================
  # Test gene list in cluster 2
  gene_names <- get_genes(res, cluster = 2)
  gene_name_to_check <- c(
    "gene189", "gene146", "gene150", "gene117", "gene190", "gene135",
    "gene165", "gene186", "gene180", "gene155", "gene144", "gene132",
    "gene106", "gene192", "gene111", "gene153", "gene108", "gene159",
    "gene133", "gene166", "gene112", "gene143", "gene126", "gene145",
    "gene184", "gene122", "gene161", "gene176", "gene193", "gene137",
    "gene168", "gene160", "gene147", "gene107", "gene129", "gene188",
    "gene101", "gene104", "gene109", "gene114", "gene127", "gene157",
    "gene173", "gene197", "gene191", "gene177", "gene103", "gene105",
    "gene194", "gene163", "gene121", "gene116", "gene162", "gene158",
    "gene178", "gene175", "gene124", "gene141", "gene171", "gene128",
    "gene142", "gene198", "gene140", "gene119", "gene183", "gene167",
    "gene200", "gene199", "gene148", "gene174", "gene113", "gene151",
    "gene136", "gene179", "gene138", "gene181", "gene154", "gene130",
    "gene118", "gene115", "gene182", "gene134", "gene195", "gene164",
    "gene185", "gene120", "gene152", "gene149", "gene172", "gene131",
    "gene196", "gene110", "gene170", "gene125", "gene156", "gene123",
    "gene169", "gene102", "gene139"
  )
  expect_equal(gene_names, gene_name_to_check)
  expect_equal(gene_names, res@gene_clusters$`2`)
  expect_that(gene_names, is_a("character"))


  # =======================================================
  # Test gene list in cluster 3
  gene_names <- get_genes(res, cluster = 3)
  gene_name_to_check <- c(
    "gene393", "gene390", "gene312", "gene314", "gene316", "gene386",
    "gene364", "gene369", "gene318", "gene315", "gene398", "gene387",
    "gene363", "gene334", "gene330", "gene320", "gene382", "gene326",
    "gene379", "gene322", "gene329", "gene324", "gene349", "gene313",
    "gene306", "gene331", "gene368", "gene350", "gene327", "gene389",
    "gene319", "gene353", "gene397", "gene366", "gene362", "gene335",
    "gene365", "gene340", "gene378", "gene336", "gene358", "gene394",
    "gene311", "gene332", "gene303", "gene341", "gene301", "gene385",
    "gene351", "gene396", "gene370", "gene339", "gene355", "gene372",
    "gene342", "gene371", "gene345", "gene367", "gene337", "gene373",
    "gene317", "gene328", "gene391", "gene392", "gene3609", "gene344",
    "gene383", "gene323", "gene310", "gene384", "gene400", "gene920",
    "gene360", "gene304", "gene2911", "gene395"
  )
  expect_equal(gene_names, gene_name_to_check)
  expect_equal(gene_names, res@gene_clusters$`3`)
  expect_that(gene_names, is_a("character"))


  # =======================================================
  # Test gene list in cluster 4
  gene_names <- get_genes(res, cluster = 4)
  gene_name_to_check <- c(
    "gene27", "gene76", "gene37", "gene32", "gene70", "gene55", "gene67",
    "gene80", "gene59", "gene16", "gene12", "gene84", "gene45", "gene88",
    "gene19", "gene61", "gene78", "gene73", "gene36", "gene29", "gene39",
    "gene7", "gene2", "gene57", "gene86", "gene9", "gene52", "gene51",
    "gene66", "gene64", "gene13", "gene4", "gene1", "gene56", "gene74",
    "gene34", "gene79", "gene93", "gene42", "gene94", "gene82", "gene83",
    "gene10", "gene63", "gene25", "gene31", "gene17", "gene87", "gene71",
    "gene28", "gene46", "gene68", "gene11", "gene90", "gene35", "gene14",
    "gene58", "gene95", "gene85", "gene96", "gene99", "gene22", "gene20",
    "gene89"
  )
  expect_equal(gene_names, gene_name_to_check)
  expect_equal(gene_names, res@gene_clusters$`4`)
  expect_that(gene_names, is_a("character"))


  # =======================================================
  # Test gene list in cluster 3 and 4
  gene_names <- get_genes(res, cluster = c(3, 4))
  gene_name_to_check <- c(
    "gene393", "gene390", "gene312", "gene314", "gene316", "gene386",
    "gene364", "gene369", "gene318", "gene315", "gene398", "gene387",
    "gene363", "gene334", "gene330", "gene320", "gene382", "gene326",
    "gene379", "gene322", "gene329", "gene324", "gene349", "gene313",
    "gene306", "gene331", "gene368", "gene350", "gene327", "gene389",
    "gene319", "gene353", "gene397", "gene366", "gene362", "gene335",
    "gene365", "gene340", "gene378", "gene336", "gene358", "gene394",
    "gene311", "gene332", "gene303", "gene341", "gene301", "gene385",
    "gene351", "gene396", "gene370", "gene339", "gene355", "gene372",
    "gene342", "gene371", "gene345", "gene367", "gene337", "gene373",
    "gene317", "gene328", "gene391", "gene392", "gene3609", "gene344",
    "gene383", "gene323", "gene310", "gene384", "gene400", "gene920",
    "gene360", "gene304", "gene2911", "gene395",
    "gene27", "gene76", "gene37", "gene32", "gene70", "gene55", "gene67",
    "gene80", "gene59", "gene16", "gene12", "gene84", "gene45", "gene88",
    "gene19", "gene61", "gene78", "gene73", "gene36", "gene29", "gene39",
    "gene7", "gene2", "gene57", "gene86", "gene9", "gene52", "gene51",
    "gene66", "gene64", "gene13", "gene4", "gene1", "gene56", "gene74",
    "gene34", "gene79", "gene93", "gene42", "gene94", "gene82", "gene83",
    "gene10", "gene63", "gene25", "gene31", "gene17", "gene87", "gene71",
    "gene28", "gene46", "gene68", "gene11", "gene90", "gene35", "gene14",
    "gene58", "gene95", "gene85", "gene96", "gene99", "gene22", "gene20",
    "gene89"
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
    "gene285", "gene273", "gene272", "gene245", "gene230", "gene286",
    "gene249", "gene278", "gene256", "gene225", "gene282", "gene236",
    "gene223", "gene214", "gene241", "gene222", "gene289", "gene277",
    "gene275", "gene206", "gene189", "gene146", "gene190", "gene155",
    "gene117", "gene150", "gene165", "gene135", "gene192", "gene186",
    "gene144", "gene106", "gene180", "gene132", "gene108", "gene126",
    "gene153", "gene133", "gene166", "gene111", "gene312", "gene393",
    "gene316", "gene386", "gene314", "gene326", "gene363", "gene334",
    "gene390", "gene364", "gene315", "gene387", "gene398", "gene379",
    "gene369", "gene382", "gene320", "gene397", "gene329", "gene306",
    "gene76", "gene55", "gene32", "gene80", "gene12", "gene59",
    "gene70", "gene27", "gene37", "gene29", "gene84", "gene16",
    "gene88", "gene67", "gene39", "gene78", "gene36", "gene57",
    "gene19", "gene7"
  )
  expect_equal(gene_names, gene_name_to_check)
  expect_equal(sort(gene_names), sort(unlist(res@top_genes, use.names = F)))
  expect_that(gene_names, is_a("character"))


  # =======================================================
  # Test top gene list in cluster 1
  gene_names <- get_genes(res, cluster = 1, top = T)
  gene_name_to_check <- c(
    "gene285", "gene273", "gene272", "gene245", "gene230", "gene286",
    "gene249", "gene278", "gene256", "gene225", "gene282", "gene236",
    "gene223", "gene214", "gene241", "gene222", "gene289", "gene277",
    "gene275", "gene206"
  )
  expect_equal(gene_names, gene_name_to_check)
  expect_equal(sort(gene_names), sort(res@top_genes$`1`))
  expect_that(gene_names, is_a("character"))


  # =======================================================
  # Test top gene list in cluster 2 and 3
  gene_names <- get_genes(res, cluster = 2:3, top = T)
  gene_name_to_check <- c(
    "gene189", "gene146", "gene190", "gene155", "gene117", "gene150",
    "gene165", "gene135", "gene192", "gene186", "gene144", "gene106",
    "gene180", "gene132", "gene108", "gene126", "gene153", "gene133",
    "gene166", "gene111", "gene312", "gene393", "gene316", "gene386",
    "gene314", "gene326", "gene363", "gene334", "gene390", "gene364",
    "gene315", "gene387", "gene398", "gene379", "gene369", "gene382",
    "gene320", "gene397", "gene329", "gene306"
  )
  expect_equal(gene_names, gene_name_to_check)
  expect_equal(sort(gene_names), sort(c(
    res@top_genes$`2`,
    res@top_genes$`3`
  )))
  expect_that(gene_names, is_a("character"))
})

# =======================================================
# Remove output files
file.remove("test.dbf_out.txt")
file.remove("test.mcl_out.txt")
