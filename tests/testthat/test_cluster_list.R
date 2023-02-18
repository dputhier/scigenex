test_that("Cheking get_genes is providing the right list of clusters", {
  #Create matrix containing 3 signatures
  m <- create_3_rnd_clust()
  res <- find_gene_clusters(data=m,
                            name = "test",
                            distance_method="pearson",
                            inflation = 2,
                            k=25,
                            fdr = 10)
  expect(length(res@cluster_list), 3)
})