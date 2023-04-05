# Set verbosity to 0
set_verbosity(0)

data("scigenex_test_I1.2")
res <- scigenex_test_I1.2

test_that("Check 'which_clust' is working.", {
  expect_true(all(which_clust(res, c('gene27', 'gene336', 'gene187')) == c(4,3,1)))
  expect_true(all(is.na(which_clust(res, c('gene27', 'gene336', 'bla'))) == c(F,F,T)))
})
