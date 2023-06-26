library(testthat)
library(org.Hs.eg.db)

# Set verbosit  y to 0
set_verbosity(0)

# Load seurat object
data(pbmc_small, package = "SeuratObject")

## Select informative genes
res <- select_genes(data=pbmc_small,
                    distance_method="pearson",
                    k=10,
                    row_sum=-Inf,
                    noise_level=0.95,
                    fdr = 1e-6)

## Cluster genes
res <- gene_clustering(object = res,
                       inflation = 1.2,
                       keep_nn = FALSE,
                       k = 5,
                       threads = 1)
set_verbosity(2)


test_that("Check if enrich_go stops when species argument is invalid", {
  expect_error(enrich_go(res, species = "Not working"))
})

res_enrich <- enrich_go(res[1,])
test_that("Check enrich_go results with all ontologies", {
  set_verbosity(0)
  
  # Check results from gene cluster 1
  expect_equal(length(res_enrich@gene_cluster_annotations), 1)
  expect_equal(res_enrich@gene_cluster_annotations$`1`@pvalueCutoff, 0.05)
  expect_equal(res_enrich@gene_cluster_annotations$`1`@pAdjustMethod, "BH")
  expect_equal(res_enrich@gene_cluster_annotations$`1`@qvalueCutoff, 0.2)
  expect_equal(res_enrich@gene_cluster_annotations$`1`@organism, "Homo sapiens")
  expect_equal(res_enrich@gene_cluster_annotations$`1`@ontology, "GOALL")
  
  
  expect_equal(res_enrich@gene_cluster_annotations$`1`@keytype, "ENTREZID")
  expect_true(length(res_enrich@gene_cluster_annotations$`1`@geneSets) > 1000)
  expect_true(res_enrich@gene_cluster_annotations$`1`@readable)
  
  
  
  expect_true("MHC class II protein complex assembly" %in% res_enrich@gene_cluster_annotations$`1`@result$Description[1:10])
  expect_true(all(c("BP", "CC") %in% res_enrich@gene_cluster_annotations$`1`@result$ONTOLOGY))
  expect_true(length(grep("CST3", res_enrich@gene_cluster_annotations$`1`@result$geneID)) > 0)
  expect_true(length(grep("HLA-DMB", res_enrich@gene_cluster_annotations$`1`@result$geneID)) > 0)
  
  
  # #########################
  # Test viz_enrich()
  # #########################
  plot_res_enrich <- viz_enrich(res_enrich)
  
  expect_equal(length(plot_res_enrich), 2)
  expect_equal(class(plot_res_enrich[[1]]), c("gg", "ggplot"))
  expect_equal(class(plot_res_enrich[[2]]), c("gg", "ggplot"))
  
  # Barplot
  expect_equal(plot_res_enrich[[1]]$labels$title, "Gene cluster 1")
  expect_equal(plot_res_enrich[[1]]$labels$fill, "p.adjust")
  expect_equal(unique(plot_res_enrich[[1]]$data$ONTOLOGY), c("BP", "CC", "MF"))

  expect_equal(plot_res_enrich[[1]]$labels$title, "Gene cluster 1")
  expect_equal(plot_res_enrich[[1]]$labels$fill, "p.adjust")
  expect_equal(unique(plot_res_enrich[[1]]$data$ONTOLOGY), c("BP", "CC", "MF"))
  
  expect_equal(plot_res_enrich[[1]]$data$ONTOLOGY, c(rep("BP", 20),
                                                     rep("CC", 20),
                                                     rep("MF", 20)))
  expect_equal(plot_res_enrich[[1]]$data$ID[1:20],
               res_enrich@gene_cluster_annotations$`1`@result$ID[1:20])
  expect_equal(as.vector(plot_res_enrich[[1]]$data$Description[1:20]),
               res_enrich@gene_cluster_annotations$`1`@result$Description[1:20])
  expect_equal(plot_res_enrich[[1]]$data$pvalue[1:20],
               res_enrich@gene_cluster_annotations$`1`@result$pvalue[1:20])
  expect_equal(plot_res_enrich[[1]]$data$p.adjust[1:20],
               res_enrich@gene_cluster_annotations$`1`@result$p.adjust[1:20])
  expect_equal(plot_res_enrich[[1]]$data$qvalue[1:20],
               res_enrich@gene_cluster_annotations$`1`@result$qvalue[1:20])
  expect_equal(plot_res_enrich[[1]]$data$geneID[1:20],
               res_enrich@gene_cluster_annotations$`1`@result$geneID[1:20])
  expect_equal(plot_res_enrich[[1]]$data$Count[1:20],
               res_enrich@gene_cluster_annotations$`1`@result$Count[1:20])
  
 
})


