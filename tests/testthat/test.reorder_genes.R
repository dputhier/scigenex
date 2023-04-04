# Set verbosity to 0
set_verbosity(0)

data("scigenex_test_I1.2")
res <- scigenex_test_I1.2

test_that("Checking if reorder_genes stops if needed", {
  expect_error(reorder_genes(res, order_by = "not_ok"))
})


test_that("Checking reorder_genes function using order_by='gene_names'...", {
  res <- reorder_genes(res, order_by = "gene_names")
  expect_equal(res@gene_clusters$`1`, c(
    "gene101", "gene102", "gene103", "gene104", "gene105", "gene106",
    "gene107", "gene108", "gene109", "gene110", "gene111", "gene112",
    "gene113", "gene114", "gene115", "gene116", "gene117", "gene118",
    "gene119", "gene120", "gene121", "gene122", "gene123", "gene124",
    "gene125", "gene126", "gene127", "gene128", "gene129", "gene130",
    "gene131", "gene132", "gene133", "gene134", "gene135", "gene136",
    "gene137", "gene138", "gene139", "gene140", "gene141", "gene142",
    "gene143", "gene144", "gene145", "gene146", "gene147", "gene148",
    "gene149", "gene150", "gene1503", "gene151", "gene152", "gene153",
    "gene154", "gene155", "gene156", "gene157", "gene158", "gene159",
    "gene160", "gene161", "gene1616", "gene162", "gene163", "gene164",
    "gene165", "gene1651", "gene166", "gene167", "gene168", "gene169",
    "gene170", "gene171", "gene172", "gene173", "gene174", "gene175",
    "gene176", "gene177", "gene178", "gene179", "gene180", "gene181",
    "gene182", "gene183", "gene184", "gene185", "gene1856", "gene186",
    "gene187", "gene188", "gene189", "gene190", "gene191", "gene1915",
    "gene192", "gene1922", "gene193", "gene194", "gene195", "gene196",
    "gene197", "gene198", "gene199", "gene200", "gene2112", "gene2184",
    "gene2223", "gene2421", "gene2553", "gene2573", "gene2662", "gene2725",
    "gene2902", "gene3442", "gene3485", "gene3488", "gene3528", "gene3610",
    "gene3779", "gene3924", "gene838"
  ))
  
  expect_equal(res@gene_clusters$`2`, c(
    "gene1931", "gene201", "gene202", "gene203", "gene206", "gene207",
    "gene208", "gene211", "gene212", "gene213", "gene214", "gene215",
    "gene217", "gene219", "gene220", "gene222", "gene223", "gene224",
    "gene225", "gene226", "gene227", "gene228", "gene229", "gene230",
    "gene232", "gene233", "gene234", "gene235", "gene236", "gene237",
    "gene238", "gene239", "gene240", "gene241", "gene242", "gene244",
    "gene245", "gene246", "gene247", "gene248", "gene249", "gene250",
    "gene251", "gene252", "gene253", "gene254", "gene255", "gene256",
    "gene257", "gene258", "gene259", "gene261", "gene262", "gene263",
    "gene266", "gene267", "gene268", "gene269", "gene270", "gene272",
    "gene273", "gene275", "gene276", "gene277", "gene278", "gene279",
    "gene280", "gene281", "gene282", "gene283", "gene284", "gene285",
    "gene286", "gene287", "gene288", "gene289", "gene290", "gene291",
    "gene292", "gene293", "gene2930", "gene294", "gene296", "gene297",
    "gene298", "gene299", "gene300", "gene3991"
  ))
  
  expect_equal(res@gene_clusters$`3`, c(
    "gene1423", "gene2490", "gene2911", "gene301", "gene302", "gene303",
    "gene304", "gene306", "gene310", "gene311", "gene312", "gene313",
    "gene314", "gene315", "gene316", "gene317", "gene318", "gene319",
    "gene320", "gene322", "gene323", "gene324", "gene325", "gene326",
    "gene327", "gene328", "gene329", "gene330", "gene331", "gene332",
    "gene334", "gene335", "gene336", "gene337", "gene339", "gene340",
    "gene341", "gene342", "gene344", "gene345", "gene349", "gene350",
    "gene351", "gene352", "gene353", "gene355", "gene358", "gene359",
    "gene360", "gene3609", "gene362", "gene363", "gene364", "gene365",
    "gene366", "gene367", "gene368", "gene369", "gene370", "gene371",
    "gene372", "gene373", "gene375", "gene376", "gene378", "gene379",
    "gene382", "gene385", "gene386", "gene387", "gene389", "gene390",
    "gene391", "gene392", "gene393", "gene394", "gene396", "gene397",
    "gene398", "gene400", "gene760"
  ))
  
  expect_equal(res@gene_clusters$`4`, c(
    "gene1", "gene10", "gene11", "gene12", "gene13", "gene14",
    "gene16", "gene17", "gene19", "gene2", "gene20", "gene21",
    "gene25", "gene27", "gene28", "gene29", "gene3", "gene31",
    "gene32", "gene34", "gene35", "gene36", "gene37", "gene39",
    "gene4", "gene42", "gene43", "gene45", "gene46", "gene51",
    "gene52", "gene54", "gene55", "gene56", "gene57", "gene58",
    "gene59", "gene6", "gene61", "gene63", "gene64", "gene66",
    "gene67", "gene7", "gene70", "gene71", "gene73", "gene74",
    "gene76", "gene78", "gene79", "gene80", "gene82", "gene83",
    "gene84", "gene85", "gene86", "gene87", "gene88", "gene9",
    "gene90", "gene93", "gene94", "gene95", "gene96", "gene98",
    "gene99"
  ))
})




test_that("Checking reorder_genes function using order_by='hclust'...", {
  set.seed(123)
  res <- reorder_genes(res, order_by = "hclust")
  
  expect_equal(res@gene_clusters$`1`, c(
    "gene1616", "gene1651", "gene3528", "gene2902", "gene2223", "gene2662",
    "gene838", "gene1856", "gene1503", "gene3610", "gene3442", "gene2421",
    "gene2573", "gene3485", "gene3488", "gene2553", "gene1922", "gene2112",
    "gene2184", "gene1915", "gene3779", "gene2725", "gene3924", "gene187",
    "gene169", "gene110", "gene172", "gene196", "gene174", "gene151",
    "gene125", "gene136", "gene164", "gene139", "gene138", "gene199",
    "gene195", "gene105", "gene135", "gene152", "gene102", "gene130",
    "gene156", "gene120", "gene147", "gene123", "gene134", "gene185",
    "gene113", "gene103", "gene173", "gene154", "gene181", "gene140",
    "gene170", "gene182", "gene104", "gene167", "gene131", "gene127",
    "gene178", "gene157", "gene116", "gene148", "gene122", "gene171",
    "gene163", "gene115", "gene149", "gene101", "gene193", "gene142",
    "gene188", "gene184", "gene145", "gene108", "gene111", "gene179",
    "gene161", "gene176", "gene137", "gene112", "gene194", "gene200",
    "gene118", "gene180", "gene162", "gene165", "gene114", "gene124",
    "gene129", "gene177", "gene128", "gene132", "gene117", "gene190",
    "gene109", "gene175", "gene168", "gene144", "gene106", "gene133",
    "gene198", "gene119", "gene166", "gene197", "gene107", "gene183",
    "gene126", "gene141", "gene159", "gene189", "gene153", "gene186",
    "gene160", "gene192", "gene143", "gene146", "gene155", "gene150",
    "gene158", "gene121", "gene191"
  ))
  
  expect_equal(res@gene_clusters$`2`, c(
    "gene1931", "gene276", "gene220", "gene228", "gene252", "gene277",
    "gene269", "gene258", "gene254", "gene279", "gene270", "gene211",
    "gene240", "gene261", "gene215", "gene284", "gene212", "gene268",
    "gene247", "gene287", "gene299", "gene234", "gene238", "gene222",
    "gene256", "gene248", "gene296", "gene250", "gene201", "gene227",
    "gene283", "gene224", "gene226", "gene262", "gene278", "gene237",
    "gene203", "gene217", "gene263", "gene239", "gene291", "gene292",
    "gene288", "gene207", "gene257", "gene213", "gene246", "gene297",
    "gene255", "gene293", "gene280", "gene229", "gene232", "gene281",
    "gene233", "gene245", "gene242", "gene208", "gene275", "gene272",
    "gene253", "gene273", "gene282", "gene300", "gene285", "gene223",
    "gene206", "gene286", "gene235", "gene225", "gene290", "gene289",
    "gene249", "gene214", "gene267", "gene244", "gene298", "gene259",
    "gene241", "gene266", "gene236", "gene294", "gene251", "gene219",
    "gene230", "gene3991", "gene202", "gene2930"
  ))
  
  expect_equal(res@gene_clusters$`3`, c(
    "gene376", "gene310", "gene302", "gene317", "gene325", "gene352",
    "gene373", "gene358", "gene339", "gene337", "gene369", "gene363",
    "gene397", "gene306", "gene378", "gene341", "gene320", "gene385",
    "gene301", "gene365", "gene350", "gene327", "gene336", "gene392",
    "gene313", "gene315", "gene323", "gene379", "gene359", "gene372",
    "gene353", "gene335", "gene370", "gene344", "gene326", "gene314",
    "gene318", "gene391", "gene303", "gene386", "gene319", "gene334",
    "gene322", "gene324", "gene311", "gene349", "gene340", "gene355",
    "gene368", "gene332", "gene342", "gene312", "gene364", "gene390",
    "gene375", "gene345", "gene398", "gene331", "gene393", "gene316",
    "gene330", "gene351", "gene362", "gene382", "gene394", "gene389",
    "gene367", "gene396", "gene366", "gene329", "gene400", "gene328",
    "gene387", "gene304", "gene360", "gene371", "gene2911", "gene760",
    "gene2490", "gene3609", "gene1423"
  ))
  
  expect_equal(res@gene_clusters$`4`, c(
    "gene43", "gene61", "gene90", "gene63", "gene3", "gene52",
    "gene28", "gene98", "gene4", "gene14", "gene7", "gene51",
    "gene27", "gene37", "gene6", "gene54", "gene31", "gene87",
    "gene86", "gene58", "gene73", "gene11", "gene19", "gene84",
    "gene82", "gene35", "gene74", "gene21", "gene95", "gene1",
    "gene71", "gene20", "gene78", "gene42", "gene10", "gene64",
    "gene85", "gene46", "gene99", "gene17", "gene25", "gene83",
    "gene12", "gene79", "gene45", "gene94", "gene39", "gene88",
    "gene57", "gene9", "gene96", "gene56", "gene16", "gene76",
    "gene59", "gene66", "gene2", "gene67", "gene13", "gene36",
    "gene70", "gene93", "gene34", "gene55", "gene29", "gene32",
    "gene80"
  ))
})




test_that("Checking reorder_genes function using order_by='correlation'...", {
  res <- reorder_genes(res, order_by = "correlation")
  
  expect_equal(res@gene_clusters$`1`, c(
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
    "gene128", "gene187", "gene170", "gene156", "gene838", "gene1616",
    "gene1856", "gene3528", "gene3485", "gene2725", "gene2573", "gene3488",
    "gene2112", "gene1503", "gene2662", "gene1651", "gene2223", "gene2421",
    "gene2553", "gene2902", "gene3442", "gene1915", "gene3779", "gene2184",
    "gene3610", "gene3924", "gene1922"
  ))
  
  expect_equal(res@gene_clusters$`2`, c(
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
    "gene1931", "gene252", "gene280", "gene246"
  ))
  
  expect_equal(res@gene_clusters$`3`, c(
    "gene315", "gene364", "gene312", "gene390", "gene320", "gene314",
    "gene334", "gene398", "gene393", "gene368", "gene316", "gene350",
    "gene387", "gene382", "gene397", "gene339", "gene369", "gene349",
    "gene363", "gene326", "gene318", "gene386", "gene341", "gene366",
    "gene324", "gene306", "gene378", "gene362", "gene329", "gene331",
    "gene303", "gene365", "gene340", "gene322", "gene353", "gene379",
    "gene330", "gene313", "gene335", "gene327", "gene345", "gene394",
    "gene311", "gene389", "gene360", "gene337", "gene301", "gene351",
    "gene372", "gene385", "gene336", "gene358", "gene400", "gene332",
    "gene323", "gene355", "gene2911", "gene342", "gene367", "gene319",
    "gene373", "gene396", "gene371", "gene310", "gene352", "gene359",
    "gene376", "gene392", "gene317", "gene370", "gene391", "gene328",
    "gene3609", "gene375", "gene760", "gene344", "gene304", "gene325",
    "gene1423", "gene302", "gene2490"
  ))
  
  expect_equal(res@gene_clusters$`4`, c(
    "gene37", "gene84", "gene19", "gene51", "gene64", "gene32",
    "gene76", "gene88", "gene7", "gene78", "gene57", "gene27",
    "gene55", "gene61", "gene80", "gene31", "gene29", "gene12",
    "gene70", "gene2", "gene59", "gene56", "gene36", "gene67",
    "gene79", "gene16", "gene34", "gene13", "gene45", "gene1",
    "gene73", "gene93", "gene28", "gene39", "gene4", "gene87",
    "gene86", "gene94", "gene9", "gene6", "gene11", "gene25",
    "gene82", "gene10", "gene83", "gene71", "gene74", "gene96",
    "gene58", "gene46", "gene42", "gene63", "gene85", "gene95",
    "gene43", "gene20", "gene66", "gene52", "gene90", "gene14",
    "gene17", "gene35", "gene54", "gene21", "gene98", "gene99",
    "gene3"
  ))
})
