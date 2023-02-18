test_that("Just some check for create_3_rnd_clust()", {
  m <- create_3_rnd_clust()
  expect_equal(ncol(m), 20)
  expect_equal(nrow(m), 4000)
  expect_equal(round(sum(m), 2), 4904.76)
})
