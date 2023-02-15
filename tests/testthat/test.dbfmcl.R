test_that("Cheking DBFMCL is providing the right number of genes", {
  #Create matrix containing 3 signatures
  set.seed(123)
  m <- matrix(rnorm(80000), nc=20)
  m[1:100, 1:10] <- m[1:100, 1:10] + 3
  m[101:200, 11:20] <- m[101:200, 11:20] + 3
  m[301:400, 5:15] <- m[201:300, 5:15] - 3
  m[201:300, c(1,5,10,15,20)] <- m[201:300, c(1,5,10,15,20)] - 3
  res <- find_gene_clusters(data=m,
                                distance_method="pearson",
                                inflation = 2,
                                k=75,
                                row_sum=-Inf,
                                highest=0.3,
                                min_nb_supporting_cell = 0,
                                fdr = 1e-8)
  
  #Test number of cluster
  expect_equal(length(res@size), 4)
  
  #Test number of genes
  expect_equal(sum(res@size), 436)
  
  #Test number of genes in each gene cluster
  expect_equal(res@size[1], 114)
  expect_equal(res@size[2], 112)
  expect_equal(res@size[3], 106)
  
  #Remove output files
  file.remove("test.dbf_out.txt")
  file.remove("test.mcl_out.txt")
})

# test_that("Checking DBF is writting the output in the correct directory...", {
#   set.seed(123)
#   #Create matrix
#   m <- matrix(rnorm(80000), nc=20)
#   m[1:100, 1:10] <- m[1:100, 1:10] + 3
#   m[101:200, 11:20] <- m[101:200, 11:20] + 3
#   m[301:400, 5:15] <- m[201:300, 5:15] - 3
#   m[201:300, c(1,5,10,15,20)] <- m[201:300, c(1,5,10,15,20)] - 3
#   res <- find_gene_clusters(data=m,
#                             distance_method="pearson",
#                             inflation = 2,
#                             k=75,
#                             row_sum=-Inf,
#                             highest=0.3,
#                             min_nb_supporting_cell = 0,
#                             fdr = 1e-8)
#   #Test results
#   out <- file.exists("tmp_output/test_.dbf_out.txt")
#   expect_equal(out, TRUE)
#   #Remove temporary files
#   file.remove("tmp_output/test_.dbf_out.txt")
#   file.remove("tmp_output/test_.mcl_out.txt")
#   file.remove("tmp_output")
# })
