# Create matrix containing 3 signatures

set_verbosity(0)

m <- create_3_rnd_clust()
res <- find_gene_clusters(
  data = m,
  name = "test",
  distance_method = "pearson",
  inflation = 2,
  k = 25,
  fdr = 10
)

test_that("Cheking get_genes is providing the right list of genes", {

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


# Remove output files
file.remove("test.dbf_out.txt")
file.remove("test.mcl_out.txt")
