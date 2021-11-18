test_that("Cheking DBF is providing the right number of genes", {
  #Create matrix
  set.seed(123)
  m <- matrix(rnorm(80000), nc=20)
  m[1:100,1:10] <- m[1:100,1:10] + 4
  m[101:200,11:20] <- m[101:200,11:20] + 3
  m[201:300,5:15] <- m[201:300,5:15] + -2
  #Run DBF
  res <- DBF(data = m,
             distance_method = "pearson",
             k = 50,
             fdr = 10,
             silent = FALSE,
             seed = 123)
  #Test results
  expect_equal(res@size, 370)
  #Remove output files
  file.remove("exprs.dbf_out.txt")
})


test_that("Checking DBF is writting the output in the correct directory...", {
  #Create matrix
  set.seed(123)
  m <- matrix(rnorm(80000), nc=20)
  m[1:100,1:10] <- m[1:100,1:10] + 4
  m[101:200,11:20] <- m[101:200,11:20] + 3
  m[201:300,5:15] <- m[201:300,5:15] + -2
  #Create an output directory in the current folder
  unlink("tmp_output", recursive = TRUE)
  dir.create(path = "tmp_output")
  #Run DBF
  res <- DBF(data = m,
             output_path = "tmp_output",
             name = "test_",
             distance_method = "pearson",
             k = 50,
             fdr = 10,
             silent = FALSE,
             seed = 123)
  #Test results
  out <- file.exists("tmp_output/test_.dbf_out.txt")
  expect_equal(out, TRUE)
  #Remove temporary files
  file.remove("tmp_output/test_.dbf_out.txt")
  unlink("tmp_output", recursive = TRUE)
})
