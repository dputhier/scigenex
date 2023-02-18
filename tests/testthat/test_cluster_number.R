test_that("Cheking get_genes is providing the right number of clusters", {
  #Create matrix containing 3 signatures
  set.seed(123)
  m <- create_3_rnd_clust()
  res <- find_gene_clusters(data=m,
                            name = "test",
                            distance_method="pearson",
                            inflation = 2,
                            k=25,
                            fdr = 10)
  expect(res@cluster_number, 3)
})