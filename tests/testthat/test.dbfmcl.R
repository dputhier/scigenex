test_that("Cheking DBFMCL is providing the right number of genes", {
  set.seed(123)
  #Create matrix containing 3 signatures
  m <- matrix(rnorm(80000), nc=20)
  m[1:100,1:10] <- m[1:100,1:10] + 4
  m[101:200,11:20] <- m[101:200,11:20] + 3
  m[201:300,5:15] <- m[201:300,5:15] + -2
  res <- DBFMCL(data=m,
                distance_method="pearson",
                av_dot_prod_min = 0,
                inflation = 2,
                k=50,
                fdr = 10,
                mcl_cmd_line = TRUE)
  #plot_clust(res, ceil = 10, floor = -10)
  
  #Test number of cluster
  expect_equal(length(res@size), 3)
  
  #Test number of genes
  expect_gt(sum(res@size), 360)
  expect_lt(sum(res@size), 390)
  
  #Test number of genes in each gene cluster
  expect_gt(res@size[1], 120)
  expect_lt(res@size[1], 150)
  
  expect_gt(res@size[2], 120)
  expect_lt(res@size[2], 150)
  
  expect_gt(res@size[3], 90)
  expect_lt(res@size[3], 120)
})





test_that("Checking DBF is writting the output in the correct directory...", {
  set.seed(123)
  #Create matrix
  m <- matrix(rnorm(80000), nc=20)
  m[1:100,1:10] <- m[1:100,1:10] + 4
  m[101:200,11:20] <- m[101:200,11:20] + 3
  m[201:300,5:15] <- m[201:300,5:15] + -2
    #Create an output directory in the current folder
    dir.create(path = "tmp_output")
  #Run DBF
  res <- DBFMCL(data = m,
                output_path = "tmp_output",
                name = "test_",
                distance_method = "pearson",
                k = 50,
                random = 3,
                fdr = 10,
                memory_used = 1024,
                silent = FALSE,
                av_dot_prod_min = 0,
                inflation = 2,
                mcl_cmd_line = TRUE)
  #Test results
  out <- file.exists("tmp_output/test_.dbf_out.txt")
  expect_equal(out, TRUE)
  #Remove temporary files
  file.remove("tmp_output/test_.dbf_out.txt")
  file.remove("tmp_output")
})
