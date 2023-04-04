# Set verbosity to 0
set_verbosity(0)

# Load seurat object
data(pbmc_small, package = "SeuratObject")

## Select informative genes
res <- select_genes(data=pbmc_small,
                    distance_method="pearson",
                    k=10,
                    row_sum=-Inf,
                    highest=0.95,
                    fdr = 1e-6)

## Cluster genes
res <- gene_clustering(object = res,
                       inflation = 1.2,
                       keep_nn = FALSE,
                       k = 5,
                       threads = 1)
set_verbosity(2)
msg <- R.utils::captureOutput(res_enrich <- enrich_go(res, 
                                             species = "Hsapiens"))


test_that("Check if enrich_go stops when species argument is invalid", {
  expect_error(enrich_go(res, species = "Not working"))
})


test_that("Check message printed by enrich_go()", {

  
  # I'm commenting some of the following message 
  # as the number of EntrezId will change from one
  # release to another one.
  expect_equal(msg[1], "|-- INFO :  Species used : Homo sapiens ")
  
  expect_equal(msg[2], "|-- INFO :  Enrichment analysis for cluster 1 ")
  expect_equal(msg[3], "|-- DEBUG :  Cluster 1")
  #expect_equal(msg[4], "\t|--Number of gene names converted to EntrezId : 20")
  #expect_equal(msg[5], "\t|--Number of gene names not converted to EntrezId : 3")
  expect_equal(msg[6], "\t|-- Calling enrichGO. ")
  expect_equal(msg[7], "|-- INFO :  Enrichment analysis for cluster 2 ")
  expect_equal(msg[8], "|-- DEBUG :  Cluster 2")
  #expect_equal(msg[9], "\t|--Number of gene names converted to EntrezId : 5")
  #expect_equal(msg[10], "\t|--Number of gene names not converted to EntrezId : 0")
  
  expect_equal(msg[11], "\t|-- Calling enrichGO. ")
  expect_equal(msg[12], "|-- INFO :  Enrichment analysis for cluster 3 ")
  expect_equal(msg[13], "|-- DEBUG :  Cluster 3")
  #expect_equal(msg[14], "\t|--Number of gene names converted to EntrezId : 3")
  #expect_equal(msg[15], "\t|--Number of gene names not converted to EntrezId : 0")
  expect_equal(msg[16], "\t|-- Calling enrichGO. ")
})

test_that("Check if enrich_go stops when species argument is invalid", {
  expect_error(enrich_go(res, species = "Not working"))
})

test_that("Check enrich_go results with all ontologies", {
  set_verbosity(0)
  
  # Check results from gene cluster 1
  expect_equal(length(res_enrich@gene_cluster_annotations), 3)
  expect_equal(res_enrich@gene_cluster_annotations$`1`@pvalueCutoff, 0.05)
  expect_equal(res_enrich@gene_cluster_annotations$`1`@pAdjustMethod, "BH")
  expect_equal(res_enrich@gene_cluster_annotations$`1`@qvalueCutoff, 0.2)
  expect_equal(res_enrich@gene_cluster_annotations$`1`@organism, "Homo sapiens")
  expect_equal(res_enrich@gene_cluster_annotations$`1`@ontology, "GOALL")
  
  # I'm commenting some of the following message 
  # as the number of EntrezId will change from one
  # release to another one.
  #expect_equal(res_enrich@gene_cluster_annotations$`1`@gene, c(
  #  "760", "2815", NA, "171558", "55287", "2791",
  #  "8848", "340205", "5196", "1191", "10398", "911",
  #  "6678", "10462", "3674", "10857", "928", "84519",
  #  "81027", "29780", "4900"
  #))
  expect_equal(res_enrich@gene_cluster_annotations$`1`@keytype, "ENTREZID")
  #expect_equal(unname(res_enrich@gene_cluster_annotations$`1`@gene2Symbol), c(
  #  "CA2", "GP9", NA, "PTCRA", "TMEM40", "GNG11",
  #  "TSC22D1", "TREML1", "PF4", "CLU", "MYL9", "CD1C",
  #  "SPARC", "CLEC10A", "ITGA2B", "PGRMC1", "CD9", "ACRBP",
  #  "TUBB1", "PARVB", "NRGN"
  #))
  # expect_true(length(res_enrich@gene_cluster_annotations$`1`@geneSets) > 1219)
  
  # Adding a more flexible test
  expect_true(length(res_enrich@gene_cluster_annotations$`1`@geneSets) > 1000)
  expect_true(res_enrich@gene_cluster_annotations$`1`@readable)
  
  
  #expect_equal(
  #  res_enrich@gene_cluster_annotations$`1`@result$Description[1:10],
  #  c(
  #    "platelet activation",
  #    "platelet degranulation",
  #    "blood coagulation",
  #    "hemostasis",
  #    "coagulation",
  #    "platelet aggregation",
  #    "regulation of megakaryocyte differentiation",
  #    "homotypic cell-cell adhesion",
  #    "regulation of myeloid cell differentiation",
  #    "megakaryocyte differentiation"
  #  )
  #)
  
  # Adding a more flexible test
  
  expect_true("MHC class II protein complex assembly" %in% res_enrich@gene_cluster_annotations$`1`@result$Description[1:10])
  
  #expect_equal(res_enrich@gene_cluster_annotations$`1`@result$ONTOLOGY, c(
  #  "BP", "BP", "BP", "BP", "BP", "BP", "BP", "BP", "BP", "BP",
  #  "BP", "BP", "CC", "CC", "CC", "CC"
  #))
  
  # Adding a more flexible test
  expect_true(all(c("BP", "CC") %in% res_enrich@gene_cluster_annotations$`1`@result$ONTOLOGY))
  
  #expect_equal(res_enrich@gene_cluster_annotations$`1`@result$GeneRatio, c(
  #  "6/19", "5/19", "6/19", "6/19", "6/19", "3/19", "3/19", "3/19", "4/19",
  #  "3/19", "4/19", "4/19", "6/20", "3/20", "5/20", "3/20"
  #))
  
  #expect_equal(round(median(
  #  res_enrich@gene_cluster_annotations$`1`@result$pvalue
  #), 4), 0)
  #expect_equal(round(mean(
  #  res_enrich@gene_cluster_annotations$`1`@result$pvalue
  #), 4), 1e-04)
  #expect_equal(round(sum(
  #  res_enrich@gene_cluster_annotations$`1`@result$pvalue
  #), 4), 0.0019)
  #expect_equal(round(sd(
  #  res_enrich@gene_cluster_annotations$`1`@result$pvalue
  #), 4), 2e-04)
  
  #expect_equal(round(median(
  #  res_enrich@gene_cluster_annotations$`1`@result$p.adjust
  #), 4), 7e-04)
  #expect_equal(round(mean(
  #  res_enrich@gene_cluster_annotations$`1`@result$p.adjust
  #), 4), 0.0051)
  #expect_equal(round(sum(
  #  res_enrich@gene_cluster_annotations$`1`@result$p.adjust
  #), 4), 0.0808)
  #expect_equal(round(sd(
  #  res_enrich@gene_cluster_annotations$`1`@result$p.adjust
  #), 4), 0.0091)
  
  
  # expect_equal(res_enrich@gene_cluster_annotations$`1`@result$geneID, c(
  #   "GP9/TREML1/PF4/MYL9/ITGA2B/CD9", "PF4/CLU/SPARC/ITGA2B/CD9",
  #   "GP9/TREML1/PF4/MYL9/ITGA2B/CD9", "GP9/TREML1/PF4/MYL9/ITGA2B/CD9",
  #   "GP9/TREML1/PF4/MYL9/ITGA2B/CD9", "MYL9/ITGA2B/CD9",
  #   "PF4/MYL9/ITGA2B", "MYL9/ITGA2B/CD9",
  #   "CA2/PF4/MYL9/ITGA2B", "PF4/MYL9/ITGA2B",
  #   "CA2/PF4/MYL9/ITGA2B", "CA2/PF4/MYL9/ITGA2B",
  #   "TREML1/PF4/CLU/SPARC/ITGA2B/CD9", "SPARC/ITGA2B/CD9",
  #   "SPARC/ITGA2B/PGRMC1/CD9/ACRBP", "PF4/CLU/SPARC"
  # ))
  
  # Adding a more flexible test
  expect_true(length(grep("CST3", res_enrich@gene_cluster_annotations$`1`@result$geneID)) > 0)
  expect_true(length(grep("HLA-DMB", res_enrich@gene_cluster_annotations$`1`@result$geneID)) > 0)
  
  #expect_equal(
  #  unique(res_enrich@gene_cluster_annotations$`1`@result$Count),
  #  c(6, 5, 3, 4)
  #)
  
  #expect_equal(sum(res_enrich@gene_cluster_annotations$`1`@result$Count), 70)
  
  
  # Check results from gene cluster 2
  expect_equal(res_enrich@gene_cluster_annotations$`2`@pvalueCutoff, 0.05)
  expect_equal(res_enrich@gene_cluster_annotations$`2`@pAdjustMethod, "BH")
  expect_equal(res_enrich@gene_cluster_annotations$`2`@qvalueCutoff, 0.2)
  expect_equal(res_enrich@gene_cluster_annotations$`2`@organism, "Homo sapiens")
  expect_equal(res_enrich@gene_cluster_annotations$`2`@ontology, "GOALL")
  
  #expect_equal(res_enrich@gene_cluster_annotations$`2`@gene, c(
  #  "2205", "27309", "649446", "282969", "5996"
  #))
  
  # Adding a more flexible test
  expect_true(all(c("3002", "6351", "83888", "5551", "10578") %in% res_enrich@gene_cluster_annotations$`2`@gene))
  
  expect_equal(res_enrich@gene_cluster_annotations$`2`@keytype, "ENTREZID")
  #expect_equal(unname(res_enrich@gene_cluster_annotations$`2`@gene2Symbol), c(
  #  "FCER1A", "ZNF330", "DLGAP1-AS1", "FUOM", "RGS1"
  #))
  
  # Adding a more flexible test
  expect_true(all(c("GZMB", "GNLY", "GZMA", "IL7R", "FCGR3A") %in% unname(res_enrich@gene_cluster_annotations$`2`@gene2Symbol)))
  
  # expect_equal(length(res_enrich@gene_cluster_annotations$`2`@geneSets), 130)
  
  # Adding a more flexible test 
  expect_true(length(res_enrich@gene_cluster_annotations$`2`@geneSets) > 100)
  
  expect_true(res_enrich@gene_cluster_annotations$`2`@readable)
  
  
  # expect_equal(
  #   res_enrich@gene_cluster_annotations$`2`@result$Description[1:10],
  #   c(
  #     "fucose metabolic process",
  #     "fucosylation",
  #     paste(
  #       "racemase and epimerase activity,",
  #       "acting on carbohydrates and derivatives"
  #     ),
  #     "racemase and epimerase activity",
  #     "immunoglobulin binding",
  #     "G-protein alpha-subunit binding",
  #     "monosaccharide binding",
  #     NA,
  #     NA,
  #     NA
  #   )
  # )
  
  # Adding a more flexible test 
  expect_true(grep("cytolysis",
                   res_enrich@gene_cluster_annotations$`2`@result$Description[1:10])
              > 0)
  
  #expect_equal(res_enrich@gene_cluster_annotations$`2`@result$ONTOLOGY, c(
  #  "BP", "BP", "MF", "MF", "MF", "MF", "MF"
  #))
  
  expect_true(all(c("BP", "MF") %in% res_enrich@gene_cluster_annotations$`2`@result$ONTOLOGY))
  
  #expect_equal(res_enrich@gene_cluster_annotations$`2`@result$GeneRatio, c(
  #  "1/4", "1/4", "1/4", "1/4", "1/4", "1/4", "1/4"
  #))
  
  # expect_equal(round(median(
  #   res_enrich@gene_cluster_annotations$`2`@result$pvalue
  # ), 4), 0.0037)
  # expect_equal(round(mean(
  #   res_enrich@gene_cluster_annotations$`2`@result$pvalue
  # ), 4), 0.0055)
  # expect_equal(round(sum(
  #   res_enrich@gene_cluster_annotations$`2`@result$pvalue
  # ), 4), 0.0383)
  # expect_equal(round(sd(
  #   res_enrich@gene_cluster_annotations$`2`@result$pvalue
  # ), 4), 0.0041)
  
  # expect_equal(round(median(
  #   res_enrich@gene_cluster_annotations$`2`@result$p.adjust
  # ), 4), 0.015)
  # expect_equal(round(mean(
  #   res_enrich@gene_cluster_annotations$`2`@result$p.adjust
  # ), 4), 0.0193)
  # expect_equal(round(sum(
  #   res_enrich@gene_cluster_annotations$`2`@result$p.adjust
  # ), 4), 0.1351)
  # expect_equal(round(sd(
  #   res_enrich@gene_cluster_annotations$`2`@result$p.adjust
  # ), 4), 0.0064)
  
  # expect_equal(res_enrich@gene_cluster_annotations$`2`@result$geneID, c(
  #   "FUOM", "FUOM", "FUOM", "FUOM", "FCER1A", "RGS1", "FUOM"
  # ))
  # Adding a more flexible tests
  expect_true(all(c("CCR7/CCL5", "GZMB/GZMH/GZMA/GZMM") %in% res_enrich@gene_cluster_annotations$`2`@result$geneID))
  expect_true(all(res_enrich@gene_cluster_annotations$`2`@result$Count > 0))
  expect_true(sum(res_enrich@gene_cluster_annotations$`2`@result$Count) > 0)
  
  
  # Check results from gene cluster 3
  expect_equal(res_enrich@gene_cluster_annotations$`3`@pvalueCutoff, 0.05)
  expect_equal(res_enrich@gene_cluster_annotations$`3`@pAdjustMethod, "BH")
  expect_equal(res_enrich@gene_cluster_annotations$`3`@qvalueCutoff, 0.2)
  expect_equal(res_enrich@gene_cluster_annotations$`3`@organism, "Homo sapiens")
  expect_equal(res_enrich@gene_cluster_annotations$`3`@ontology, "GOALL")
  # expect_equal(res_enrich@gene_cluster_annotations$`3`@gene, c(
  #   "4345", "27240", "100423062"
  # ))
  # Adding a more flexible tests
  expect_true(length(res_enrich@gene_cluster_annotations$`3`@gene) > 0)
  expect_equal(res_enrich@gene_cluster_annotations$`3`@keytype, "ENTREZID")
  expect_true(all(c("PF4", "PPBP", "GP9") %in% unname(res_enrich@gene_cluster_annotations$`3`@gene2Symbol)))
  
  # expect_equal(length(res_enrich@gene_cluster_annotations$`3`@geneSets), 315)
  
  # Adding a more flexible tests
  expect_true(length(res_enrich@gene_cluster_annotations$`3`@geneSets) > 200)
  expect_true(res_enrich@gene_cluster_annotations$`3`@readable)
  
  expect_true(length(grep("platelet activation",  res_enrich@gene_cluster_annotations$`3`@result$Description[1:10])) > 0)
  
  # expect_equal(res_enrich@gene_cluster_annotations$`3`@result$ONTOLOGY, c(
  #   "BP", "BP", "BP", "BP", "BP", "BP", "BP", "BP", "BP", "BP",
  #   "BP", "BP", "BP", "BP", "BP", "BP", "BP", "BP", "BP", "BP",
  #   "BP", "BP", "BP", "BP", "BP", "BP", "BP", "BP", "CC", "MF",
  #   "MF", "MF", "MF", "MF", "MF"
  # ))
  # Adding a more flexible tests
  expect_true("BP" %in% res_enrich@gene_cluster_annotations$`3`@result$ONTOLOGY)
  
  # expect_equal(res_enrich@gene_cluster_annotations$`3`@result$GeneRatio, c(
  #   "1/3", "1/3", "1/3", "1/3", "1/3", "1/3", "1/3", "1/3", "1/3", "1/3",
  #   "1/3", "1/3", "1/3", "1/3", "1/3", "1/3", "1/3", "1/3", "1/3", "1/3",
  #   "1/3", "1/3", "1/3", "1/3", "1/3", "1/3", "1/3", "1/3", "1/3", "1/3",
  #   "1/3", "1/3", "1/3", "1/3", "1/3"
  # ))
  
  # expect_equal(round(median(
  #   res_enrich@gene_cluster_annotations$`3`@result$pvalue
  # ), 4), 0.007)
  # expect_equal(round(mean(
  #   res_enrich@gene_cluster_annotations$`3`@result$pvalue
  # ), 4), 0.0077)
  # expect_equal(round(sum(
  #   res_enrich@gene_cluster_annotations$`3`@result$pvalue
  # ), 4), 0.2699)
  # expect_equal(round(sd(
  #   res_enrich@gene_cluster_annotations$`3`@result$pvalue
  # ), 4), 0.0049)
  # 
  # expect_equal(round(median(
  #   res_enrich@gene_cluster_annotations$`3`@result$p.adjust
  # ), 4), 0.0412)
  # expect_equal(round(mean(
  #   res_enrich@gene_cluster_annotations$`3`@result$p.adjust
  # ), 4), 0.0371)
  # expect_equal(round(sum(
  #   res_enrich@gene_cluster_annotations$`3`@result$p.adjust
  # ), 4), 1.297)
  # expect_equal(round(sd(
  #   res_enrich@gene_cluster_annotations$`3`@result$p.adjust
  # ), 4), 0.0102)
  
  
  # expect_equal(res_enrich@gene_cluster_annotations$`3`@result$geneID, c(
  #   "CD200", "CD200", "CD200", "CD200", "CD200", "CD200",
  #   "CD200", "CD200", "CD200", "SIT1", "CD200", "CD200",
  #   "CD200", "CD200", "CD200", "CD200", "CD200", "CD200",
  #   "SIT1", "CD200", "CD200", "CD200", "CD200", "CD200",
  #   "CD200", "CD200", "SIT1", "CD200", "IGLL5", "CD200",
  #   "SIT1", "CD200", "CD200", "IGLL5", "IGLL5"
  # ))
  
  # Adding a more flexible tests
  expect_true(all(c("CD9/MYL9/FERMT3", "CLU/GPX1") %in% res_enrich@gene_cluster_annotations$`3`@result$geneID))
  expect_true(sum(res_enrich@gene_cluster_annotations$`3`@result$Count) > 20)
  
  
  # #########################
  # Test viz_enrich()
  # #########################
  plot_res_enrich <- viz_enrich(res_enrich)
  
  expect_equal(length(plot_res_enrich), 6)
  expect_equal(class(plot_res_enrich[[1]]), c("gg", "ggplot"))
  expect_equal(class(plot_res_enrich[[2]]), c("gg", "ggplot"))
  expect_equal(class(plot_res_enrich[[3]]), c("gg", "ggplot"))
  expect_equal(class(plot_res_enrich[[4]]), c("gg", "ggplot"))
  expect_equal(class(plot_res_enrich[[5]]), c("gg", "ggplot"))
  expect_equal(class(plot_res_enrich[[6]]), c("gg", "ggplot"))
  
  # Barplot
  expect_equal(plot_res_enrich[[1]]$labels$title, "Gene cluster 1")
  expect_equal(plot_res_enrich[[1]]$labels$fill, "p.adjust")
  expect_equal(unique(plot_res_enrich[[1]]$data$ONTOLOGY), c("BP", "CC", "MF"))
  
  #expect_equal(round(sum(plot_res_enrich[[1]]$data$GeneRatio), 4), 3.6395)
  #expect_equal(round(mean(plot_res_enrich[[1]]$data$GeneRatio), 4), 0.2275)
  #expect_equal(round(median(plot_res_enrich[[1]]$data$GeneRatio), 4), 0.2105)
  #expect_equal(round(sd(plot_res_enrich[[1]]$data$GeneRatio), 4), 0.0685)
  
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
  
  # Dotplot
  expect_equal(plot_res_enrich[[4]]$labels$title, "Gene cluster 1")
  expect_equal(plot_res_enrich[[4]]$labels$size, "Count")
  expect_equal(plot_res_enrich[[4]]$labels$colour, "p.adjust")
  expect_equal(plot_res_enrich[[4]]$labels$x, "GeneRatio")
  
  expect_equal(unique(plot_res_enrich[[4]]$data$ONTOLOGY), c("BP", "CC", "MF"))
  
  #expect_equal(round(sum(plot_res_enrich[[4]]$data$GeneRatio), 4), 3.6395)
  #expect_equal(round(mean(plot_res_enrich[[4]]$data$GeneRatio), 4), 0.2275)
  #expect_equal(round(median(plot_res_enrich[[4]]$data$GeneRatio), 4), 0.2105)
  #expect_equal(round(sd(plot_res_enrich[[4]]$data$GeneRatio), 4), 0.0685)
  
  expect_equal(plot_res_enrich[[4]]$data$ONTOLOGY[1:20],
               res_enrich@gene_cluster_annotations$`1`@result$ONTOLOGY[1:20])
  expect_equal(plot_res_enrich[[4]]$data$ID[1:20],
               res_enrich@gene_cluster_annotations$`1`@result$ID[1:20])
  expect_equal(as.vector(plot_res_enrich[[4]]$data$Description[1:20]),
               res_enrich@gene_cluster_annotations$`1`@result$Description[1:20])
  expect_equal(plot_res_enrich[[4]]$data$pvalue[1:20],
               res_enrich@gene_cluster_annotations$`1`@result$pvalue[1:20])
  expect_equal(plot_res_enrich[[4]]$data$p.adjust[1:20],
               res_enrich@gene_cluster_annotations$`1`@result$p.adjust[1:20])
  expect_equal(plot_res_enrich[[4]]$data$qvalue[1:20],
               res_enrich@gene_cluster_annotations$`1`@result$qvalue[1:20])
  expect_equal(plot_res_enrich[[4]]$data$geneID[1:20],
               res_enrich@gene_cluster_annotations$`1`@result$geneID[1:20])
  expect_equal(plot_res_enrich[[4]]$data$Count[1:20],
               res_enrich@gene_cluster_annotations$`1`@result$Count[1:20])
  
  
  # Test arguments clusters set to 1 and 3
  plot_res_enrich <- viz_enrich(res_enrich, clusters = c(1,3))
  expect_equal(length(plot_res_enrich), 4)
  expect_equal(plot_res_enrich[[1]]$labels$title, "Gene cluster 1")
  expect_equal(plot_res_enrich[[2]]$labels$title, "Gene cluster 3")
  expect_equal(plot_res_enrich[[3]]$labels$title, "Gene cluster 1")
  expect_equal(plot_res_enrich[[4]]$labels$title, "Gene cluster 3")
  
  plot_res_enrich <- viz_enrich(res_enrich, clusters = 2)
  expect_equal(length(plot_res_enrich), 2)
  expect_equal(plot_res_enrich[[1]]$labels$title, "Gene cluster 2")
  expect_equal(plot_res_enrich[[2]]$labels$title, "Gene cluster 2")
  
  # Test arguments type set to barplot
  plot_res_enrich <- viz_enrich(res_enrich, clusters = 1, type = "barplot")
  expect_equal(plot_res_enrich[[1]]$labels$fill, "p.adjust")
  
  # Test arguments type set to dotplot
  plot_res_enrich <- viz_enrich(res_enrich, clusters = 1, type = "dotplot")
  expect_equal(plot_res_enrich[[1]]$labels$size, "Count")
  expect_equal(plot_res_enrich[[1]]$labels$colour, "p.adjust")
  expect_equal(plot_res_enrich[[1]]$labels$x, "GeneRatio")
  
  # Test arguments nb_terms set to 10
  plot_res_enrich <- viz_enrich(res_enrich,
                                clusters = 1,
                                type = "dotplot",
                                nb_terms = 10)
  # expect_equal(nrow(
  #   plot_res_enrich[[1]]$data[plot_res_enrich[[1]]$data$ONTOLOGY == "BP", ]
  #   ),
  #   10
  # )
  
  # expect_equal(nrow(
  #   plot_res_enrich[[1]]$data[plot_res_enrich[[1]]$data$ONTOLOGY == "CC", ]
  # ),
  # 4
  # )
})


test_that("Check enrich_go results with BP ontology", {
  set_verbosity(0)
  
  # Use Biological Process (BP)
  res_enrich <- enrich_go(res, species = "Hsapiens", ontology = "BP")
  expect_equal(res_enrich@gene_cluster_annotations$`1`@ontology, "BP")
  expect_equal(res_enrich@gene_cluster_annotations$`2`@ontology, "BP")
  expect_equal(res_enrich@gene_cluster_annotations$`3`@ontology, "BP")
  
  # ## For cluster 1
  # expect_equal(sort(unique(
  #   res_enrich@gene_cluster_annotations$`1`@result$Count
  # )), seq(1, 6))
  # expect_equal(sum(res_enrich@gene_cluster_annotations$`1`@result$Count), 553)
  
  # expect_equal(
  #   length(
  #     res_enrich@gene_cluster_annotations$`1`@result$Description
  #   ),
  #   465
  # )
  
  # expect_equal(
  #   res_enrich@gene_cluster_annotations$`1`@result$Description[1:10],
  #   c(
  #     "platelet activation",
  #     "platelet degranulation",
  #     "blood coagulation",
  #     "hemostasis",
  #     "coagulation",
  #     "platelet aggregation",
  #     "regulation of megakaryocyte differentiation",
  #     "homotypic cell-cell adhesion",
  #     "regulation of myeloid cell differentiation",
  #     "megakaryocyte differentiation"
  #   )
  # )
  
  # Adding a more flexible test
  expect_true(length(grep("neutrophil", 
                          res_enrich@gene_cluster_annotations$`1`@result$Description[1:10])) > 0) 
  ## For cluster 2
  expect_true(sum(res_enrich@gene_cluster_annotations$`2`@result$Count) > 0)
  
  # expect_equal(
  #   length(
  #     res_enrich@gene_cluster_annotations$`2`@result$Description
  #   ),
  #   12
  # )
  
  # expect_true(
  #   res_enrich@gene_cluster_annotations$`2`@result$Description[1:10],
  #   c(
  #     "fucose metabolic process",
  #     "fucosylation",
  #     paste(
  #       "adenylate cyclase-inhibiting G protein-coupled",
  #       "receptor signaling pathway"
  #     ),
  #     "regulation of G protein-coupled receptor signaling pathway",
  #     "Fc-epsilon receptor signaling pathway",
  #     paste(
  #       "adenylate cyclase-modulating G protein-coupled",
  #       "receptor signaling pathway"
  #     ),
  #     "Fc receptor signaling pathway",
  #     "hexose metabolic process",
  #     "glycosylation",
  #     "monosaccharide metabolic process"
  #   )
  # )
  
  expect_true(length(grep("T cell", 
                          res_enrich@gene_cluster_annotations$`2`@result$Description[1:10])) > 0) 
  
  
  ## For cluster 3
  expect_true(sum(res_enrich@gene_cluster_annotations$`3`@result$Count) > 50)
  
  expect_true(
    length(
      res_enrich@gene_cluster_annotations$`3`@result$Description
    ) > 50
  )
  expect_true(
    length(grep("platelet" , res_enrich@gene_cluster_annotations$`3`@result$Description[1:10])) > 0)
})
