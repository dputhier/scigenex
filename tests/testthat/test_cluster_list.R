test_that("Cheking get_genes is providing the right list of clusters", {
  #Create matrix containing 3 signatures
  set.seed(123)
  m <- matrix(rnorm(40000), nc=20)
  m[1:100,1:10] <- m[1:100,1:10] + 4
  m[101:200,11:20] <- m[101:200,11:20] + 3
  m[201:300,5:15] <- m[201:300,5:15] + -2
  res <- find_gene_clusters(data=m,
                            name = "test",
                            distance_method="pearson",
                            inflation = 2,
                            k=25,
                            fdr = 10)
  expect(length(res@cluster_list), 3)
})