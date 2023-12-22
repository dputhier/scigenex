library(Seurat)
library(ggplot2)

# Set verbosity to 0
set_verbosity(0)

load_example_dataset("7870305/files/lymph_node_tiny_2")
identity <- Idents(lymph_node_tiny_2)
classes <- ifelse(Idents(lymph_node_tiny_2) == 7, 1, 0)
names(classes) <- names(identity)
strat <- stratify_seurat(lymph_node_tiny_2, ident=classes)

test_that("Check stratify_seurat() is working.", {

  testthat::expect_true(length(strat) == 4)
  testthat::expect_true(length(strat[[1]]$dist2border_classes) == 442)
  testthat::expect_true(all(unname(table(strat[[1]]$dist2border_classes)) == c(395, 47)))
})



p <- plot_stratum(strat[[1]], gene_name="CCL19", polar = TRUE, 
                  colours_stratum = rev(discrete_palette(5, "De1")))

test_that("Check plot_stratum() is working.", {
  testthat::expect_true(ggplot2::is.ggplot(p))
})