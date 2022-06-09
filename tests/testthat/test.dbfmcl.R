test_that("Cheking DBFMCL is providing the right number of genes", {
  #Create matrix containing 3 signatures
  set.seed(123)
  m <- matrix(rnorm(80000), nc=20)
  m[1:100,1:10] <- m[1:100,1:10] + 4
  m[101:200,11:20] <- m[101:200,11:20] + 3
  m[201:300,5:15] <- m[201:300,5:15] + -2
  res <- DBFMCL(data=m,
                name = "test",
                distance_method="pearson",
                av_dot_prod_min = 0,
                inflation = 2,
                k=50,
                fdr = 10,
                seed = 123)
  
  #Test number of cluster
  expect_equal(length(res@size), 3)
  
  #Test number of genes
  expect_equal(sum(res@size), 372)
  
  #Test number of genes in each gene cluster
  expect_equal(res@size[1], 134)
  expect_equal(res@size[2], 134)
  expect_equal(res@size[3], 104)
  
  #Remove output files
  file.remove("test.dbf_out.txt")
  file.remove("test.mcl_out.txt")
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
                fdr = 10,
                silent = FALSE,
                av_dot_prod_min = 0,
                inflation = 2,
                seed = 123)
  #Test results
  out <- file.exists("tmp_output/test_.dbf_out.txt")
  expect_equal(out, TRUE)
  #Remove temporary files
  file.remove("tmp_output/test_.dbf_out.txt")
  file.remove("tmp_output/test_.mcl_out.txt")
  file.remove("tmp_output")
})
