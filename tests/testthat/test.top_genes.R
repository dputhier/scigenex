# Create matrix containing 3 signatures

set_verbosity(0)

m <- create_3_rnd_clust()


test_that("Cheking get_genes is providing the right list of genes - pearson", {
  
  res <- find_gene_clusters(
    data = m,
    name = "test",
    distance_method = "pearson",
    inflation = 2,
    k = 25,
    fdr = 10
  )
  
  # ========================================
  # Top 20
  res_20 <- suppressWarnings(top_genes(res, cluster = "all", top = 20))
  
  # Test top genes list in all cluster
  gene_name_to_check <- c(
    "gene76", "gene37", "gene59", "gene80", "gene55", "gene67",
    "gene84", "gene32", "gene27", "gene78", "gene12", "gene88",
    "gene9", "gene29", "gene2", "gene19", "gene39", "gene70",
    "gene73", "gene1", "gene146", "gene155", "gene150", "gene189",
    "gene190", "gene144", "gene165", "gene117", "gene135", "gene108",
    "gene126", "gene133", "gene145", "gene111", "gene180", "gene106",
    "gene153", "gene186", "gene112", "gene143", "gene109", "gene192",
    "gene104", "gene147", "gene166", "gene167", "gene122", "gene132",
    "gene175", "gene197", "gene129", "gene174", "gene140", "gene172",
    "gene156", "gene123"
  )
  expect_equal(unlist(res_20@top_genes, use.names = FALSE), gene_name_to_check)
  expect_that(res@top_genes, is_a("list"))
  
  
  # Test top genes list in cluster 1
  gene_name_to_check <- c(
    "gene76", "gene37", "gene59", "gene80", "gene55", "gene67",
    "gene84", "gene32", "gene27", "gene78", "gene12", "gene88",
    "gene9", "gene29", "gene2", "gene19", "gene39", "gene70",
    "gene73", "gene1"
  )
  expect_equal(res_20@top_genes$`1`, gene_name_to_check)
  
  
  # Test top genes list in cluster 2 and 3
  gene_name_to_check <- c(
    "gene146", "gene155", "gene150", "gene189", "gene190", "gene144",
    "gene165", "gene117", "gene135", "gene108", "gene126", "gene133",
    "gene145", "gene111", "gene180", "gene106", "gene153", "gene186",
    "gene112", "gene143", "gene109", "gene192", "gene104", "gene147",
    "gene166", "gene167", "gene122", "gene132", "gene175", "gene197",
    "gene129", "gene174", "gene140", "gene172", "gene156", "gene123"
  )
  expect_equal(c(res_20@top_genes$`2`, res_20@top_genes$`3`),
               gene_name_to_check)
  
  
  # ========================================
  # Top 10
  res_10 <- top_genes(res, cluster = "all", top = 10)
  
  # Test top genes list in all cluster
  gene_name_to_check <- c(
    "gene76", "gene37", "gene59", "gene80", "gene55", "gene67",
    "gene84", "gene32", "gene27", "gene78", "gene146", "gene155",
    "gene150", "gene189", "gene190", "gene144", "gene165", "gene117",
    "gene135", "gene108", "gene109", "gene192", "gene104", "gene147",
    "gene166", "gene167", "gene122", "gene132", "gene175", "gene197"
  )
  expect_equal(unlist(res_10@top_genes, use.names = FALSE), gene_name_to_check)
  expect_that(res@top_genes, is_a("list"))
  
  
  # Test top genes list in cluster 1
  gene_name_to_check <- c(
    "gene76", "gene37", "gene59", "gene80", "gene55", "gene67",
    "gene84", "gene32", "gene27", "gene78"
  )
  expect_equal(res_10@top_genes$`1`, gene_name_to_check)
  
  
  # Test top genes list in cluster 2 and 3
  gene_name_to_check <- c(
    "gene146", "gene155", "gene150", "gene189", "gene190", "gene144",
    "gene165", "gene117", "gene135", "gene108", "gene109", "gene192",
    "gene104", "gene147", "gene166", "gene167", "gene122", "gene132",
    "gene175", "gene197"
  )
  expect_equal(
    c(res_10@top_genes$`2`, res_10@top_genes$`3`),
    gene_name_to_check
  )
  
  
  
  # ========================================
  # Top 5
  res_5 <- top_genes(res, cluster = "all", top = 5)
  
  # Test top genes list in all cluster
  gene_name_to_check <- c(
    "gene76", "gene37", "gene59", "gene80", "gene55", "gene146",
    "gene155", "gene150", "gene189", "gene190", "gene109", "gene192",
    "gene104", "gene147", "gene166"
  )
  expect_equal(unlist(res_5@top_genes, use.names = FALSE), gene_name_to_check)
  expect_that(res@top_genes, is_a("list"))
  
  
  # Test top genes list in cluster 1
  gene_name_to_check <- c(
    "gene76", "gene37", "gene59", "gene80", "gene55"
  )
  expect_equal(res_5@top_genes$`1`, gene_name_to_check)
  
  
  # Test top genes list in cluster 2 and 3
  gene_name_to_check <- c(
    "gene146", "gene155", "gene150", "gene189", "gene190", "gene109",
    "gene192", "gene104", "gene147", "gene166"
  )
  expect_equal(c(res_5@top_genes$`2`, res_5@top_genes$`3`), gene_name_to_check)
  
  
  
  # ========================================
  # Top 100
  res_100 <- suppressWarnings(top_genes(res, cluster = "all", top = 100))
  
  # Test top genes list in all cluster
  gene_name_to_check <- c(
    "gene76", "gene37", "gene59", "gene80", "gene55", "gene67",
    "gene84", "gene32", "gene27", "gene78", "gene12", "gene88",
    "gene9", "gene29", "gene2", "gene19", "gene39", "gene70",
    "gene73", "gene1", "gene16", "gene7", "gene58", "gene57",
    "gene34", "gene86", "gene36", "gene74", "gene71", "gene51",
    "gene64", "gene93", "gene13", "gene45", "gene79", "gene66",
    "gene96", "gene94", "gene4", "gene56", "gene42", "gene25",
    "gene61", "gene89", "gene46", "gene28", "gene31", "gene20",
    "gene6", "gene82", "gene63", "gene52", "gene10", "gene11",
    "gene81", "gene21", "gene23", "gene95", "gene87", "gene47",
    "gene68", "gene43", "gene41", "gene40", "gene38", "gene26",
    "gene91", "gene146", "gene155", "gene150", "gene189", "gene190",
    "gene144", "gene165", "gene117", "gene135", "gene108", "gene126",
    "gene133", "gene145", "gene111", "gene180", "gene106", "gene153",
    "gene186", "gene112", "gene143", "gene184", "gene137", "gene176",
    "gene168", "gene161", "gene159", "gene188", "gene194", "gene160",
    "gene105", "gene162", "gene114", "gene118", "gene101", "gene124",
    "gene142", "gene107", "gene191", "gene103", "gene177", "gene200",
    "gene173", "gene136", "gene128", "gene183", "gene198", "gene141",
    "gene195", "gene185", "gene152", "gene134", "gene109", "gene192",
    "gene104", "gene147", "gene166", "gene167", "gene122", "gene132",
    "gene175", "gene197", "gene129", "gene174", "gene140", "gene172",
    "gene156", "gene123"
  )
  expect_equal(unlist(res_100@top_genes, use.names = FALSE), gene_name_to_check)
  expect_that(res@top_genes, is_a("list"))
  
  
  # Test top genes list in cluster 1
  gene_name_to_check <- c(
    "gene76", "gene37", "gene59", "gene80", "gene55", "gene67",
    "gene84", "gene32", "gene27", "gene78", "gene12", "gene88",
    "gene9", "gene29", "gene2", "gene19", "gene39", "gene70",
    "gene73", "gene1", "gene16", "gene7", "gene58", "gene57",
    "gene34", "gene86", "gene36", "gene74", "gene71", "gene51",
    "gene64", "gene93", "gene13", "gene45", "gene79", "gene66",
    "gene96", "gene94", "gene4", "gene56", "gene42", "gene25",
    "gene61", "gene89", "gene46", "gene28", "gene31", "gene20",
    "gene6", "gene82", "gene63", "gene52", "gene10", "gene11",
    "gene81", "gene21", "gene23", "gene95", "gene87", "gene47",
    "gene68", "gene43", "gene41", "gene40", "gene38", "gene26",
    "gene91"
  )
  expect_equal(res_100@top_genes$`1`, gene_name_to_check)
  
  
  # Test top genes list in cluster 2 and 3
  gene_name_to_check <- c(
    "gene146", "gene155", "gene150", "gene189", "gene190", "gene144",
    "gene165", "gene117", "gene135", "gene108", "gene126", "gene133",
    "gene145", "gene111", "gene180", "gene106", "gene153", "gene186",
    "gene112", "gene143", "gene184", "gene137", "gene176", "gene168",
    "gene161", "gene159", "gene188", "gene194", "gene160", "gene105",
    "gene162", "gene114", "gene118", "gene101", "gene124", "gene142",
    "gene107", "gene191", "gene103", "gene177", "gene200", "gene173",
    "gene136", "gene128", "gene183", "gene198", "gene141", "gene195",
    "gene185", "gene152", "gene134", "gene109", "gene192", "gene104",
    "gene147", "gene166", "gene167", "gene122", "gene132", "gene175",
    "gene197", "gene129", "gene174", "gene140", "gene172", "gene156",
    "gene123"
  )
  expect_equal(
    c(res_100@top_genes$`2`, res_100@top_genes$`3`),
    gene_name_to_check
  )
})





test_that("Cheking get_genes is providing the right list of genes - cosine", {
  
  res <- find_gene_clusters(
    data = m,
    name = "test",
    distance_method = "cosine",
    inflation = 2,
    k = 25,
    fdr = 10
  )
  
  # ========================================
  # Top 20
  res_20 <- suppressWarnings(top_genes(res, cluster = "all", top = 20))
  
  # Test top genes list in all cluster
  gene_name_to_check <- c(
    "gene189", "gene146", "gene165", "gene155", "gene150", "gene186",
    "gene192", "gene153", "gene180", "gene112", "gene190", "gene160",
    "gene143", "gene111", "gene166", "gene176", "gene132", "gene144",
    "gene159", "gene108", "gene76", "gene80", "gene9", "gene67",
    "gene55", "gene59", "gene32", "gene29", "gene16", "gene12",
    "gene25", "gene78", "gene39", "gene2", "gene93", "gene36",
    "gene66", "gene89", "gene73", "gene74", "gene19", "gene84",
    "gene88", "gene58", "gene86", "gene31", "gene87", "gene37",
    "gene7", "gene45", "gene23", "gene21", "gene11", "gene14",
    "gene4", "gene51", "gene91", "gene2897", "gene1369", "gene1093",
    "gene2625", "gene1509", "gene2476", "gene1721", "gene3705", "gene2184",
    "gene3474", "gene3387", "gene1942", "gene424", "gene667", "gene2351",
    "gene3964", "gene1582", "gene582", "gene2005", "gene3108", "gene1927",
    "gene944", "gene762", "gene3042", "gene1924", "gene2164", "gene2839",
    "gene2160", "gene1809", "gene813", "gene3", "gene69", "gene97",
    "gene44", "gene17", "gene33", "gene50", "gene52", "gene49",
    "gene10", "gene53"
  )
  expect_equal(unlist(res_20@top_genes, use.names = FALSE), gene_name_to_check)
  expect_that(res@top_genes, is_a("list"))
})



test_that("Cheking get_genes is providing the right list of genes - euclidean", {
  
  res <- find_gene_clusters(
    data = m,
    name = "test",
    distance_method = "euclidean",
    inflation = 2,
    k = 25,
    fdr = 10, 
    min_cluster_size = 5
  )
  
  # ========================================
  # Top 20
  res_20 <- suppressWarnings(top_genes(res, cluster = "all", top = 20))
  
  # Test top genes list in all cluster
  gene_name_to_check <- c(
    "gene754", "gene908", "gene1123", "gene3438", "gene463", "gene2445",
    "gene1129", "gene325", "gene2489", "gene3362", "gene2624", "gene1291",
    "gene1771", "gene2877", "gene3674", "gene3225", "gene2653", "gene3735",
    "gene1205", "gene323", "gene2079", "gene3610", "gene3788", "gene1222",
    "gene2263"
  )
  expect_equal(unlist(res_20@top_genes, use.names = FALSE), gene_name_to_check)
  expect_that(res@top_genes, is_a("list"))
})




test_that("Cheking get_genes is providing the right list of genes - kendall", {
  
  res <- find_gene_clusters(
    data = m,
    name = "test",
    distance_method = "kendall",
    inflation = 2,
    k = 25,
    fdr = 10
  )
  
  # ========================================
  # Top 20
  res_20 <- suppressWarnings(top_genes(res, cluster = "all", top = 20))
  
  # Test top genes list in all cluster
  gene_name_to_check <- c(
    "gene94", "gene56", "gene9", "gene81", "gene25", "gene22",
    "gene36", "gene99", "gene34", "gene17", "gene13", "gene29",
    "gene85", "gene89", "gene32", "gene45", "gene46", "gene39",
    "gene2", "gene55", "gene176", "gene116", "gene126", "gene191",
    "gene160", "gene141", "gene109", "gene171", "gene155", "gene153",
    "gene143", "gene148", "gene121", "gene137", "gene192", "gene101",
    "gene158", "gene159", "gene188", "gene108", "gene84", "gene21",
    "gene87", "gene19", "gene11", "gene37", "gene7", "gene58",
    "gene63", "gene27", "gene23", "gene31", "gene86", "gene61",
    "gene6", "gene33", "gene41", "gene75", "gene88", "gene90",
    "gene105", "gene195", "gene135", "gene162", "gene132", "gene117",
    "gene180", "gene178", "gene139", "gene130", "gene138", "gene200",
    "gene165", "gene157", "gene102", "gene172", "gene169", "gene154",
    "gene181", "gene179", "gene1415", "gene2621", "gene1970", "gene3194",
    "gene2420", "gene1663", "gene727", "gene1828", "gene1714", "gene1702",
    "gene3792", "gene3052", "gene2680", "gene2581", "gene673", "gene2798",
    "gene1960", "gene1263", "gene779", "gene742", "gene1808", "gene1406",
    "gene879", "gene1953", "gene1079", "gene427", "gene129", "gene104",
    "gene190", "gene128", "gene177", "gene174", "gene106", "gene182",
    "gene193", "gene110", "gene134"
  )
  expect_equal(unlist(res_20@top_genes, use.names = FALSE), gene_name_to_check)
  expect_that(res@top_genes, is_a("list"))
})
