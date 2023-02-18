test_that("Just some check about create_4_rnd_clust()", {
  m <- create_4_rnd_clust()
  expect_equal(ncol(m), 20)
  expect_equal(nrow(m), 4000)
  expect_equal(round(sum(m), 2), 2390.88)
})
