# Create matrix containing 4 signatures
m <- create_4_rnd_clust()

set_verbosity(0)

## A rather stringent version
res <- DBF(
  data = m,
  name = "test",
  distance_method = "pearson",
  k = 75,
  row_sum = -Inf,
  highest = 0.3,
  fdr = 1e-8,
  output_path = tempdir()
)

test_that("Checking results obtained with plot_dist()", {
  p_dist <- plot_dist(res)
  
  # Check class
  expect_equal(class(p_dist), c("gg", "ggplot"))
  
  #Check graph elements
  expect_equal(p_dist$labels$y, "Count")
  expect_equal(p_dist$labels$x, "Distance with KNN")
  expect_equal(p_dist$labels$fill, "Type")
  
  p_dist_infos <- ggplot_build(p_dist)
  expect_equal(unique(p_dist_infos$data[[1]]$fill), c("#36949D", "#FB8500"))
  expect_equal(p_dist_infos$plot$layers[[1]]$aes_params$alpha, 0.5)
  expect_equal(p_dist_infos$plot$layers[[1]]$computed_stat_params$bins, 150)
  
  # Check data
  expect_equal(round(mean(p_dist$data$DKNN), 4), 0.528)
  expect_equal(p_dist$data[p_dist$data$Type == "Observed", "DKNN"],
               unname(res@dbf_output$dknn))
  expect_equal(p_dist$data[p_dist$data$Type == "Simulated", "DKNN"],
               unname(res@dbf_output$simulated_dknn))
})



test_that("Checking if plot_dist() stops when object argument is not a\
          ClusterSet object", {
  expect_error(plot_dist(object = "Not a ClusterSet object"))
})
