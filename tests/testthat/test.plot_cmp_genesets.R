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

