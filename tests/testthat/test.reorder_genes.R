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




# writeLines(paste0('"', paste(res@gene_clusters$`1`, collapse = '","'), '"'), "genes1.txt")
# writeLines(paste0('"', paste(res@gene_clusters$`2`, collapse = '","'), '"'), "genes2.txt")
# writeLines(paste0('"', paste(res@gene_clusters$`3`, collapse = '","'), '"'), "genes3.txt")
# writeLines(paste0('"', paste(res@gene_clusters$`4`, collapse = '","'), '"'), "genes4.txt")

test_that("Checking reorder_genes function using order_by='gene_names'...", {
  res <- reorder_genes(res, order_by = "gene_names")
  expect_equal(res@gene_clusters$`1`, c(
    "gene1931", "gene201", "gene202", "gene203", "gene2036", "gene205",
    "gene206", "gene207", "gene208", "gene209", "gene210", "gene211",
    "gene212", "gene213", "gene214", "gene215", "gene2157", "gene216",
    "gene217", "gene218", "gene219", "gene220", "gene221", "gene222",
    "gene223", "gene224", "gene225", "gene226", "gene227", "gene228",
    "gene229", "gene230", "gene231", "gene232", "gene233", "gene234",
    "gene235", "gene236", "gene237", "gene238", "gene239", "gene240",
    "gene241", "gene242", "gene243", "gene244", "gene245", "gene246",
    "gene247", "gene248", "gene249", "gene250", "gene251", "gene252",
    "gene253", "gene254", "gene255", "gene256", "gene257", "gene258",
    "gene259", "gene260", "gene261", "gene262", "gene263", "gene264",
    "gene265", "gene266", "gene267", "gene268", "gene269", "gene270",
    "gene271", "gene272", "gene273", "gene274", "gene275", "gene276",
    "gene277", "gene278", "gene279", "gene280", "gene281", "gene282",
    "gene283", "gene284", "gene285", "gene286", "gene287", "gene288",
    "gene289", "gene290", "gene291", "gene292", "gene293", "gene294",
    "gene295", "gene296", "gene297", "gene298", "gene299", "gene300",
    "gene3991", "gene847"
  ))

  expect_equal(res@gene_clusters$`2`, c(
    "gene101", "gene102", "gene103", "gene104", "gene105", "gene106",
    "gene107", "gene108", "gene109", "gene110", "gene111", "gene112",
    "gene113", "gene114", "gene115", "gene116", "gene117", "gene118",
    "gene119", "gene120", "gene121", "gene122", "gene123", "gene124",
    "gene125", "gene126", "gene127", "gene128", "gene129", "gene130",
    "gene131", "gene132", "gene133", "gene134", "gene135", "gene136",
    "gene137", "gene138", "gene139", "gene140", "gene141", "gene142",
    "gene143", "gene144", "gene145", "gene146", "gene147", "gene148",
    "gene149", "gene150", "gene151", "gene152", "gene153", "gene154",
    "gene155", "gene156", "gene157", "gene158", "gene159", "gene160",
    "gene161", "gene162", "gene163", "gene164", "gene165", "gene166",
    "gene167", "gene168", "gene169", "gene170", "gene171", "gene172",
    "gene173", "gene174", "gene175", "gene176", "gene177", "gene178",
    "gene179", "gene180", "gene181", "gene182", "gene183", "gene184",
    "gene185", "gene186", "gene188", "gene189", "gene190", "gene191",
    "gene192", "gene193", "gene194", "gene195", "gene196", "gene197",
    "gene198", "gene199", "gene200"
  ))

  expect_equal(res@gene_clusters$`3`, c(
    "gene2911", "gene301", "gene303", "gene304", "gene306", "gene310",
    "gene311", "gene312", "gene313", "gene314", "gene315", "gene316",
    "gene317", "gene318", "gene319", "gene320", "gene322", "gene323",
    "gene324", "gene326", "gene327", "gene328", "gene329", "gene330",
    "gene331", "gene332", "gene334", "gene335", "gene336", "gene337",
    "gene339", "gene340", "gene341", "gene342", "gene344", "gene345",
    "gene349", "gene350", "gene351", "gene353", "gene355", "gene358",
    "gene360", "gene3609", "gene362", "gene363", "gene364", "gene365",
    "gene366", "gene367", "gene368", "gene369", "gene370", "gene371",
    "gene372", "gene373", "gene378", "gene379", "gene382", "gene383",
    "gene384", "gene385", "gene386", "gene387", "gene389", "gene390",
    "gene391", "gene392", "gene393", "gene394", "gene395", "gene396",
    "gene397", "gene398", "gene400", "gene920"
  ))

  expect_equal(res@gene_clusters$`4`, c(
    "gene1", "gene10", "gene11", "gene12", "gene13", "gene14", "gene16",
    "gene17", "gene19", "gene2", "gene20", "gene22", "gene25", "gene27",
    "gene28", "gene29", "gene31", "gene32", "gene34", "gene35", "gene36",
    "gene37", "gene39", "gene4", "gene42", "gene45", "gene46", "gene51",
    "gene52", "gene55", "gene56", "gene57", "gene58", "gene59", "gene61",
    "gene63", "gene64", "gene66", "gene67", "gene68", "gene7", "gene70",
    "gene71", "gene73", "gene74", "gene76", "gene78", "gene79", "gene80",
    "gene82", "gene83", "gene84", "gene85", "gene86", "gene87", "gene88",
    "gene89", "gene9", "gene90", "gene93", "gene94", "gene95", "gene96",
    "gene99"
  ))
})




test_that("Checking reorder_genes function using order_by='hclust'...", {
  set.seed(123)
  res <- reorder_genes(res, order_by = "hclust")

  expect_equal(res@gene_clusters$`1`, c(
    "gene847", "gene209", "gene2036", "gene278", "gene202", "gene3991",
    "gene2157", "gene243", "gene1931", "gene231", "gene218", "gene276",
    "gene264", "gene271", "gene205", "gene246", "gene297", "gene293",
    "gene255", "gene281", "gene245", "gene233", "gene288", "gene260",
    "gene257", "gene213", "gene280", "gene229", "gene232", "gene295",
    "gene210", "gene207", "gene242", "gene208", "gene275", "gene272",
    "gene253", "gene277", "gene258", "gene269", "gene283", "gene227",
    "gene273", "gene282", "gene300", "gene285", "gene223", "gene286",
    "gene206", "gene225", "gene235", "gene290", "gene249", "gene289",
    "gene248", "gene296", "gene214", "gene267", "gene244", "gene298",
    "gene259", "gene241", "gene266", "gene236", "gene294", "gene251",
    "gene230", "gene219", "gene254", "gene279", "gene270", "gene211",
    "gene240", "gene261", "gene215", "gene284", "gene212", "gene268",
    "gene224", "gene262", "gene226", "gene217", "gene263", "gene239",
    "gene291", "gene252", "gene216", "gene299", "gene287", "gene247",
    "gene221", "gene203", "gene237", "gene234", "gene238", "gene256",
    "gene222", "gene292", "gene265", "gene250", "gene201", "gene274",
    "gene220", "gene228"
  ))

  expect_equal(res@gene_clusters$`2`, c(
    "gene169", "gene110", "gene172", "gene174", "gene196", "gene151",
    "gene125", "gene136", "gene164", "gene138", "gene139", "gene199",
    "gene195", "gene135", "gene105", "gene152", "gene130", "gene102",
    "gene156", "gene120", "gene147", "gene123", "gene134", "gene185",
    "gene113", "gene173", "gene103", "gene181", "gene154", "gene140",
    "gene170", "gene182", "gene104", "gene167", "gene131", "gene127",
    "gene157", "gene178", "gene116", "gene148", "gene122", "gene171",
    "gene163", "gene115", "gene149", "gene193", "gene101", "gene142",
    "gene184", "gene188", "gene145", "gene111", "gene108", "gene179",
    "gene161", "gene176", "gene137", "gene112", "gene194", "gene200",
    "gene118", "gene180", "gene162", "gene165", "gene114", "gene124",
    "gene129", "gene177", "gene128", "gene132", "gene117", "gene190",
    "gene109", "gene175", "gene168", "gene144", "gene106", "gene133",
    "gene198", "gene119", "gene166", "gene197", "gene107", "gene183",
    "gene126", "gene159", "gene141", "gene189", "gene153", "gene186",
    "gene160", "gene192", "gene143", "gene146", "gene150", "gene155",
    "gene158", "gene191", "gene121"
  ))

  expect_equal(res@gene_clusters$`3`, c(
    "gene3609", "gene920", "gene328", "gene310", "gene395", "gene387",
    "gene304", "gene378", "gene306", "gene384", "gene337", "gene369",
    "gene363", "gene397", "gene383", "gene373", "gene358", "gene339",
    "gene301", "gene365", "gene350", "gene327", "gene336", "gene392",
    "gene313", "gene315", "gene379", "gene323", "gene353", "gene335",
    "gene370", "gene344", "gene314", "gene326", "gene317", "gene372",
    "gene366", "gene329", "gene400", "gene396", "gene341", "gene320",
    "gene385", "gene318", "gene303", "gene391", "gene386", "gene319",
    "gene334", "gene322", "gene324", "gene311", "gene349", "gene340",
    "gene355", "gene368", "gene332", "gene342", "gene312", "gene390",
    "gene364", "gene345", "gene398", "gene331", "gene393", "gene316",
    "gene330", "gene351", "gene382", "gene362", "gene394", "gene389",
    "gene367", "gene360", "gene371", "gene2911"
  ))

  expect_equal(res@gene_clusters$`4`, c(
    "gene63", "gene68", "gene52", "gene61", "gene90", "gene7", "gene51",
    "gene27", "gene37", "gene28", "gene95", "gene10", "gene64", "gene85",
    "gene74", "gene71", "gene89", "gene34", "gene55", "gene29", "gene32",
    "gene80", "gene93", "gene20", "gene78", "gene42", "gene45", "gene94",
    "gene88", "gene39", "gene1", "gene22", "gene46", "gene99", "gene57",
    "gene9", "gene96", "gene56", "gene16", "gene76", "gene59", "gene66",
    "gene67", "gene2", "gene25", "gene17", "gene83", "gene12", "gene79",
    "gene4", "gene14", "gene13", "gene70", "gene36", "gene82", "gene35",
    "gene31", "gene86", "gene87", "gene58", "gene73", "gene11", "gene84",
    "gene19"
  ))
})




test_that("Checking reorder_genes function using order_by='correlation'...", {
  res <- reorder_genes(res, order_by = "correlation")

  expect_equal(res@gene_clusters$`1`, c(
    "gene285", "gene273", "gene272", "gene245", "gene230", "gene286",
    "gene249", "gene278", "gene256", "gene225", "gene282", "gene236",
    "gene223", "gene214", "gene241", "gene222", "gene289", "gene277",
    "gene275", "gene206", "gene298", "gene251", "gene266", "gene290",
    "gene248", "gene267", "gene235", "gene300", "gene293", "gene219",
    "gene233", "gene217", "gene259", "gene229", "gene212", "gene283",
    "gene258", "gene287", "gene239", "gene255", "gene250", "gene207",
    "gene253", "gene299", "gene232", "gene237", "gene297", "gene203",
    "gene263", "gene216", "gene244", "gene269", "gene247", "gene284",
    "gene262", "gene281", "gene294", "gene242", "gene279", "gene215",
    "gene292", "gene288", "gene254", "gene240", "gene226", "gene246",
    "gene268", "gene205", "gene291", "gene257", "gene234", "gene270",
    "gene238", "gene227", "gene221", "gene202", "gene208", "gene209",
    "gene271", "gene261", "gene260", "gene220", "gene295", "gene296",
    "gene280", "gene213", "gene201", "gene210", "gene211", "gene274",
    "gene243", "gene218", "gene224", "gene3991", "gene265", "gene231",
    "gene252", "gene228", "gene276", "gene264", "gene2157", "gene1931",
    "gene847", "gene2036"
  ))

  expect_equal(res@gene_clusters$`2`, c(
    "gene189", "gene146", "gene190", "gene155", "gene117", "gene150",
    "gene165", "gene135", "gene192", "gene186", "gene144", "gene106",
    "gene180", "gene132", "gene108", "gene126", "gene153", "gene133",
    "gene166", "gene111", "gene112", "gene159", "gene145", "gene168",
    "gene143", "gene122", "gene176", "gene109", "gene184", "gene137",
    "gene161", "gene160", "gene147", "gene104", "gene101", "gene188",
    "gene114", "gene105", "gene193", "gene107", "gene127", "gene129",
    "gene191", "gene124", "gene197", "gene177", "gene173", "gene157",
    "gene171", "gene194", "gene116", "gene103", "gene121", "gene162",
    "gene128", "gene163", "gene158", "gene140", "gene175", "gene119",
    "gene142", "gene141", "gene183", "gene198", "gene167", "gene138",
    "gene178", "gene200", "gene148", "gene199", "gene136", "gene151",
    "gene179", "gene118", "gene174", "gene113", "gene115", "gene130",
    "gene154", "gene195", "gene182", "gene134", "gene181", "gene164",
    "gene185", "gene152", "gene120", "gene149", "gene172", "gene110",
    "gene131", "gene125", "gene102", "gene196", "gene170", "gene169",
    "gene123", "gene156", "gene139"
  ))

  expect_equal(res@gene_clusters$`3`, c(
    "gene312", "gene393", "gene316", "gene386", "gene314", "gene326",
    "gene363", "gene334", "gene390", "gene364", "gene315", "gene387",
    "gene398", "gene379", "gene369", "gene382", "gene320", "gene397",
    "gene329", "gene306", "gene330", "gene368", "gene350", "gene365",
    "gene324", "gene322", "gene318", "gene349", "gene340", "gene362",
    "gene353", "gene327", "gene366", "gene389", "gene378", "gene331",
    "gene319", "gene303", "gene385", "gene335", "gene367", "gene341",
    "gene336", "gene311", "gene339", "gene358", "gene313", "gene332",
    "gene345", "gene396", "gene301", "gene351", "gene371", "gene394",
    "gene337", "gene373", "gene323", "gene3609", "gene342", "gene355",
    "gene372", "gene391", "gene317", "gene328", "gene400", "gene360",
    "gene2911", "gene395", "gene370", "gene392", "gene344", "gene310",
    "gene384", "gene383", "gene304", "gene920"
  ))

  expect_equal(res@gene_clusters$`4`, c(
    "gene76", "gene55", "gene32", "gene80", "gene12", "gene59", "gene70",
    "gene27", "gene37", "gene29", "gene84", "gene16", "gene88", "gene67",
    "gene39", "gene78", "gene36", "gene57", "gene19", "gene7", "gene2",
    "gene13", "gene45", "gene1", "gene34", "gene51", "gene73", "gene56",
    "gene9", "gene64", "gene93", "gene86", "gene61", "gene79", "gene74",
    "gene94", "gene31", "gene4", "gene42", "gene71", "gene82", "gene66",
    "gene25", "gene52", "gene83", "gene10", "gene96", "gene46", "gene11",
    "gene14", "gene58", "gene20", "gene63", "gene22", "gene87", "gene35",
    "gene89", "gene95", "gene17", "gene28", "gene90", "gene85", "gene99",
    "gene68"
  ))
})
