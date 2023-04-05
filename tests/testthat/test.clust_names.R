data("scigenex_test_I1.2")
data("scigenex_test_I4")

test_that("Check clust_names() is working.", {
  a <- unname(gene_cluster(scigenex_test_I1.2))
  b <- unname(clust_names(scigenex_test_I1.2))
  expect_true(all(a==b))
})