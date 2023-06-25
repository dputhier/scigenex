set_verbosity(0)
set.seed(123)
set_1 <- list(letters[1:10], letters[11:20], letters[21:26], LETTERS[1:26])
x <- sample(c(letters, LETTERS))
set_2 <- list(x[1:5], x[6:20], LETTERS[1:10])

test_that("Check diagrams are working (plot_cmp_genesets).", {
  expect_error(plot_cmp_genesets(set_1 = NULL, set_2 = NULL))
  expect_true(ggplot2::is.ggplot(plot_cmp_genesets(set_1 = set_1, set_2 = set_2, stat = "jaccard")))
  expect_true(ggplot2::is.ggplot(plot_cmp_genesets(set_1 = set_2, set_2 = set_1, stat = "jaccard")))
  expect_true(ggplot2::is.ggplot(plot_cmp_genesets(set_1 = set_2, set_2 = set_1, stat = "hypergeom")))
  expect_true(ggplot2::is.ggplot(plot_cmp_genesets(set_1 = set_2, set_2 = set_1, stat = "diff_set_1")))
  expect_true(ggplot2::is.ggplot(plot_cmp_genesets(set_1 = set_2, set_2 = set_1, stat = "diff_set_2")))
  expect_true(ggplot2::is.ggplot(plot_cmp_genesets(set_1 = set_2, set_2 = set_1, stat = "size_set_1")))
  expect_true(ggplot2::is.ggplot(plot_cmp_genesets(set_1 = set_2, set_2 = set_1, stat = "size_set_2")))
  expect_true(ggplot2::is.ggplot(plot_cmp_genesets(set_1 = set_2, set_2 = set_1, layout="square")))
})

library(testthat)
set_verbosity(0)

test_that("compare_genesets returns a message when input is empty", {
  expect_error(compare_genesets(set_1 = NULL, set_2 = NULL, stat = "jaccard"))
})


test_that("compare_genesets returns a message when input is not a list", {
  expect_error(compare_genesets(set_1 = "not_a_list", set_2 = "not_a_list", stat = "jaccard"))
})


set_1 <- list(c("gene1", "gene2", "gene3"), c("gene1", "gene2"))
set_2 <- list(c("gene1", "gene2", "gene4", "gene5"), c("gene2", "gene3"))
test_that("compare_genesets returns correct output when size_set_1 is used", {
  expect_equal(as.vector(compare_genesets(set_1, set_2, stat = "size_set_1")), 
               c(3, 2, 3, 2))
})


test_that("compare_genesets returns correct output when size_set_2 is used", {
  expect_equal(as.vector(compare_genesets(set_1, set_2, stat = "size_set_2")), 
               c(4, 4, 2, 2))
})

set_1 <- list(c("gene1", "gene2", "gene3", "gene5", "gene6"), c("gene1", "gene2", "gene4"))
set_2 <- list(c("gene1", "gene2", "gene4", "gene5"), c("gene2", "gene3"))
test_that("compare_genesets returns correct output when diff_set_1 is used", {
  expect_equal(as.vector(compare_genesets(set_1, set_2, stat = "diff_set_1")), 
               c(2, 0, 3, 2))
})


test_that("compare_genesets returns correct output when diff_set_2 is used", {
  expect_equal(as.vector(compare_genesets(set_1, set_2, stat = "diff_set_2")), 
               c(1, 1, 0, 1))
})

test_that("compare_genesets returns correct output when method='union' is used", {
  expect_equal(as.vector(compare_genesets(set_1, set_2, stat = "union")), 
               c(6, 4, 5, 4))
})

test_that("compare_genesets returns correct output when method='intersection' is used", {
  expect_equal(as.vector(compare_genesets(set_1, set_2, stat = "intersection")), 
               c(3, 3, 2, 1))
})
i <- as.vector(compare_genesets(set_1, set_2, stat = "intersection"))
u <- as.vector(compare_genesets(set_1, set_2, stat = "union"))
j <- i/u
test_that("compare_genesets returns correct output when method='intersection' is used", {
  expect_equal(as.vector(compare_genesets(set_1, set_2, stat = "jaccard")), j) 
})

test_that("compare_genesets returns correct output when method='hypergeom' is used", {
  expect_equal(round(as.vector(compare_genesets(set_1, set_2, stat = "hypergeom")), 2), 
               c(1, 0.2, 0.67, 0.8))
})

