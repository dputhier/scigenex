library(testthat)

testthat::test_that("Test the indexing operator of the ClusterSet", {
  
  clusters <- 1:5
  cluster_size <- c(12,10,8,6,4)
  m <- matrix(rnorm(20000), nc=500, nr=40)
  colnames(m) <- paste0("S_", 1:500)
  gn <- paste0("G_", 1:40)
  cls <- unlist(mapply(rep, clusters, cluster_size))
  gene_cluster <- split(gn, cls) 
  rownames(m) <- paste0("G_", 1:40)
  cls <- unlist(mapply(rep, clusters, cluster_size))
  gene_clusters_metadata <- list()
  gene_clusters_metadata$cluster_id <- as.numeric(names(gene_cluster))
  gene_clusters_metadata$number <- length(clusters)
  gene_clusters_metadata$size <- cluster_size
  names(gene_clusters_metadata$size) <- clusters
  gene_cluster_annotations <- list()
  top_genes <- lapply(gene_cluster, "[", 1:3)
  dbf_output <- list()
  dbf_output$dknn <- runif(1000)
  dbf_output$simulated_dknn <- runif(1000)
  dbf_output$critical_distance <- 0.5
  dbf_output$fdr <- 0.05
  dbf_output$center <- do.call(rbind, 
                               lapply(split.data.frame(as.data.frame(m), 
                                                       as.factor(cls)), 
                                      apply, 2, mean))
  
  colnames(dbf_output$center) <- colnames(m)
  rownames(dbf_output$center) <- names(gene_cluster)
  
  res <- new("ClusterSet",
              data = m,
              gene_clusters = split(gn, cls),
              top_genes = top_genes,
              gene_clusters_metadata = gene_clusters_metadata,
              gene_cluster_annotations = gene_cluster_annotations,
              cells_metadata = data.frame(cells_barcode=colnames(m), 
                                          row.names = colnames(m)),
              dbf_output = dbf_output,
              parameters = list(filename="oP0H7t6Sc0", distance_method="pearson", 
                                k=50, inflation=2, noise_level=0.3, fdr=0.05, 
                                min_nb_supporting_cell=0, min_pct_gene_expressed=40, 
                                min_cluster_size=10, row_sum=-Inf, seed=123)
            )
  testthat::expect_equal(dim(res[1:3,4:5]), c(30, 2))
  testthat::expect_equal(dim(res[c(2,4),1:5]), c(16, 5))
  testthat::expect_equal(length(res[c(2,4),1:5]@gene_clusters), 2)
  testthat::expect_equal(length(res[c(2,4),1:5]@gene_clusters), 2)
  testthat::expect_equal(length(res[c(1,2,4),]@gene_clusters), 3)
  testthat::expect_equal(dim(res[logical(), logical()]), c(0,0))
  testthat::expect_equal(dim(res[c("1", "2"), ]), c(22, 500))
  testthat::expect_equal(res[c("1", "2"), ]@gene_clusters_metadata$number, 2)
  testthat::expect_equal(res[2:4, ]@gene_clusters_metadata$number, 3)
  testthat::expect_equal(res[c("1", "2", "3"), ]@gene_clusters_metadata$number, 3)
  testthat::expect_equal(res[logical(), ]@gene_clusters_metadata$number, 0)
  testthat::expect_equal(dim(res[, ]), c(40, 500))
  testthat::expect_equal(dim(res[, 10:20]), c(40, 11))
  testthat::expect_equal(dim(res[1:3, 10:20]@cells_metadata), c(11, 1))
  testthat::expect_equal(dim(res[1:3, 10:20]@cells_metadata), c(11, 1))
  testthat::expect_equal(dim(res[1, 1]@cells_metadata), c(1, 1))
  testthat::expect_equal(dim(res[1, 1]), c(12, 1))
  testthat::expect_equal(length(res[1:2, 1]@top_genes), 2)
  testthat::expect_equal(length(res[1:2, c()]@top_genes), 2)
  
  testthat::expect_identical(res[c(1, 3), ], res[c(1, 3), seq_len(ncol(res))])
  testthat::expect_identical(res[, c("S_408", "S_409")], res[, c(408, 409)])
  testthat::expect_identical(res[c(1, 3), c("S_3", "S_4")], res[c(1, 3), 3:4])
  testthat::expect_identical(res[-c(1, 2), ], res[c(3:nclust(res)), ])
  testthat::expect_identical(res[as.character(2:5), c("S_3", "S_4")], res[2:5, 3:4])
  testthat::expect_identical(res[as.character(2:5), ], res[2:5, ])
  testthat::expect_identical(res[,paste0("S_",2:5)], res[,2:5])
  testthat::expect_identical(res[c("2", "4"), ]@gene_clusters, res[c(2, 4), ]@gene_clusters)
  testthat::expect_identical(res[c("1", "3"), ]@gene_clusters_metadata$number, res[c(1, 3), ]@gene_clusters_metadata$number)
  testthat::expect_identical(res[3, ]@gene_clusters_metadata, res["3", ]@gene_clusters_metadata)
  testthat::expect_error(res["a", "b"])
  })

test_that("Test messages when looking at ClusterSet object", {
  obj <- new("ClusterSet")
  msg_clusterset <- capture.output(obj)
  
  expect_equal(msg_clusterset[1], "\t\tAn object of class ClusterSet")
  expect_equal(msg_clusterset[2], "\t\tName: ")
  expect_equal(msg_clusterset[3], "\t\tMemory used:  2520 ")
  expect_equal(msg_clusterset[4], "\t\tNumber of cells:  0 ")
  expect_equal(msg_clusterset[5], "\t\tNumber of informative genes:  0 ")
  expect_equal(msg_clusterset[6], "\t\tNumber of gene clusters:  ")
  expect_equal(msg_clusterset[7], "\t\tThis object contains the following informations:")
  expect_equal(msg_clusterset[8], "\t\t\t -  data ")
  expect_equal(msg_clusterset[9], "\t\t\t -  gene_clusters ")
  expect_equal(msg_clusterset[10], "\t\t\t -  top_genes ")
  expect_equal(msg_clusterset[11], "\t\t\t -  gene_clusters_metadata ")
  expect_equal(msg_clusterset[12], "\t\t\t -  gene_cluster_annotations ")
  expect_equal(msg_clusterset[13], "\t\t\t -  cells_metadata ")
  expect_equal(msg_clusterset[14], "\t\t\t -  dbf_output ")
  expect_equal(msg_clusterset[15], "\t\t\t -  parameters ")
  
})
set_verbosity(0)
load_example_dataset("7871581/files/pbmc3k_medium_clusters")

test_that("Check clust_names() is working.", {
  a <- unique(unname(gene_cluster(pbmc3k_medium_clusters)))
  b <- clust_names(pbmc3k_medium_clusters)
  expect_true(all(as.character(a)==b))
})

test_that("Check clust_names() is working.", {
  a <- show_methods()
  expect_true(all(c("[", "clust_names", "clust_size", "cluster_stats", "col_names", 
                "dim", "enrich_go", "gene_cluster", "nclust", "rename_clust", 
                "row_names", "show", "top_genes", "viz_enrich", "which_clust")  %in% a))
})# Set verbosity to 0
set_verbosity(0)
library(Seurat)

load_example_dataset("7871581/files/pbmc3k_medium_clusters")
load_example_dataset("7871581/files/pbmc3k_medium")
gn <- pbmc3k_medium_clusters

test_that("Check 'rename' is working.", {
  expect_true(nclust(rename_clust(gn, 1:nclust(gn))) == nclust(gn))
  expect_true(all(names(rename_clust(gn, nclust(gn):1)@gene_clusters) == 15:1))
  rn <- rename_clust(gn, letters[1:nclust(gn)])
  expect_true(ggplot2::is.ggplot(plot_profiles(rn, ident = Seurat::Idents(pbmc3k_medium))))
})

set_verbosity(0)
library(Seurat)
load_example_dataset("7871581/files/pbmc3k_medium_clusters")
res <- pbmc3k_medium_clusters


test_that("Check 'which_clust' is working.", {
  expect_true(all(which_clust(res, c('APOBEC3A', 'RPL11', 'PF4')) == c(3,1,2)))
  expect_true(all(is.na(which_clust(res, c('APOBEC3A', 'RPL11', 'bla'))) == c(F,F,T)))
})
set_verbosity(0)
load_example_dataset("7871581/files/pbmc3k_medium_clusters")

test_that("Check clust_names() is working.", {
  a <- unique(unname(gene_cluster(pbmc3k_medium_clusters)))
  b <- clust_names(pbmc3k_medium_clusters)
  expect_true(all(as.character(a)==b))
})

# Set verbosity to 0
set_verbosity(0)

load_example_dataset("7871581/files/pbmc3k_medium_clusters")

res <- pbmc3k_medium_clusters

test_that("Checking ncol()", {
  expect_equal(ncol(res), 361)
})

test_that("Checking nrow()", {
  expect_equal(nrow(res), 291)
})

test_that("Checking dim()", {
  expect_equal(dim(res), c(291, 361))
})


test_that("Checking colnames()", {
  expect_equal(col_names(res)[1:10], c("GATCTACTGGTGAG-1", "ACAGTGACTCACCC-1", "AGACGTACAGAGGC-1", 
                                       "GACGTAACCTGTGA-1", "TATACAGATCCAGA-1", "CGGATAACAGCTCA-1", "TTACCATGGTTGAC-1", 
                                       "TCCCACGATCATTC-1", "ATAGCGTGCAGATC-1", "TGTAGGTGTGCTGA-1"))
  expect_equal(tail(col_names(res)), c("GGCATATGGGGAGT-1", "TTACGTACGTTCAG-1", "GGAACACTTCAGAC-1", 
                                      "ATCATCTGACACCA-1", "ACGAACTGGCTATG-1", "TAACACCTTGTTTC-1"))
})

test_that("Checking rownames()", {
  expect_equal(row_names(res), c("MALAT1", "RPL9", "RPL23A", "RPL21", "RPS25", "RPS19", "RPS3", 
                                 "RPS12", "RPL18A", "RPL11", "RPL32", "RPLP2", "RPS18", "RPL3", 
                                 "RPL29", "RPL14", "RPS3A", "RPS15A", "RPS23", "RPL27A", "RPS16", 
                                 "RPL12", "EEF1A1", "RPL31", "RPS4X", "RPS8", "RPL28", "RPS29", 
                                 "RPS10", "RPS13", "RPS11", "RPS26", "RPL35A", "TPT1", "BTG1", 
                                 "RPL15", "RPL18", "RPS20", "RPS28", "RPL36", "RPL5", "EEF1B2", 
                                 "RPL7A", "RPL17", "NACA", "RPL27", "RPL26", "RPS9", "RPLP1", 
                                 "RPSAP58", "RPL6", "PTCRA", "ACRBP", "SDPR", "GNG11", "PF4", 
                                 "HIST1H2AC", "TMEM40", "TUBB1", "SPARC", "PPBP", "GP9", "MMD", 
                                 "AP001189.4", "CLU", "RGS18", "MPP1", "CA2", "F13A1", "TPM1", 
                                 "MYL9", "TREML1", "NGFRAP1", "PTGS1", "GAS2L1", "GRAP2", "ODC1", 
                                 "SNCA", "TPM4", "CTSA", "ILK", "PDLIM1", "ACTN1", "FERMT3", "PLA2G12A", 
                                 "HIST1H2BK", "SNN", "APP", "YWHAH", "ASAH1", "RAP1B", "MFSD1", 
                                 "GSN", "MARCH2", "H2AFJ", "TPTEP1", "PARVB", "MYL12A", "CD151", 
                                 "PVALB", "TYMP", "COTL1", "CTSS", "CFD", "S100A11", "LGALS2", 
                                 "GABARAP", "APOBEC3A", "TNFSF13B", "CPVL", "MNDA", "CD14", "IGSF6", 
                                 "CSTB", "MAFB", "AP2S1", "WARS", "IFI6", "LILRB2", "LILRB4", 
                                 "C1orf162", "SLC7A7", "IFNGR2", "SOD2", "TGFBI", "BLVRA", "CAPG", 
                                 "BID", "CTSZ", "NAAA", "FCGR2A", "CD300LF", "S100A10", "GAPDH", 
                                 "ALDH2", "EIF4EBP1", "GRINA", "NUP214", "CD302", "ATG3", "CAMK1", 
                                 "CYBB", "CNPY3", "GCA", "POU2F2", "CCL5", "GZMA", "CST7", "CTSW", 
                                 "HLA-C", "FGFBP2", "GZMB", "GNLY", "PRF1", "SPON2", "CCL4", "PTPRCAP", 
                                 "CLIC3", "IL2RB", "HOPX", "CD247", "AKR1C3", "IGFBP7", "KLRF1", 
                                 "CD99", "CD7", "XCL1", "EFHD2", "PLEK", "HEMGN", "EGFL7", "AC137932.6", 
                                 "hsa-mir-1199", "SEC14L5", "HEXIM2", "GATA2", "SH3BGRL2", "KIAA0513", 
                                 "PCP2", "PYGL", "DPY19L1", "RAB27B", "CLIC4", "LYL1", "FAM110A", 
                                 "RILP", "DAB2", "KIFC3", "BTK", "SPOCD1", "ARHGAP21", "ITGB3", 
                                 "ALOX12", "SCFD2", "ATP9A", "FAM212A", "ABCC3", "RP11-359I18.5", 
                                 "TGFB1I1", "CTNNAL1", "ENDOD1", "NEXN", "ABHD16A", "ARG2", "GLA", 
                                 "WIPI1", "PPP1R14A", "CLEC1B", "CMTM5", "GP1BA", "RP11-367G6.3", 
                                 "HIST1H2BD", "FNTB", "SENCR", "SEPT4", "MLH3", "SCGB1C1", "HIST2H2BE", 
                                 "C1orf198", "H1F0", "GLUL", "LGALSL", "SMOX", "SSX2IP", "AC147651.3", 
                                 "C19orf33", "FAM212B", "SERPINE2", "ZNF778", "TMCC2", "TRIM58", 
                                 "FN3K", "NFE2", "MGLL", "MOB3B", "HLA-DPA1", "HLA-DRB1", "HLA-DRA", 
                                 "HLA-DRB5", "HLA-DPB1", "CD74", "HLA-DQA1", "HLA-DQB1", "FAM26F", 
                                 "LY86", "HLA-DMB", "SLC25A6", "CD37", "RNASET2", "C2orf88", "CLDN5", 
                                 "HIST1H2BJ", "FAM63A", "SPHK1", "SLC40A1", "TUBA1C", "TMEM91", 
                                 "XPNPEP1", "PRKAR2B", "TJP2", "AC092295.7", "ENKUR", "LCN2", 
                                 "GFI1B", "ANKRD9", "PRUNE", "SMIM5", "PTPN18", "CD3D", "IL32", 
                                 "CD3E", "LDHB", "IL7R", "AES", "CD27", "HES4", "FCGR3A", "RHOC", 
                                 "IFITM2", "ABI3", "CLIC1", "RP11-879F14.2", "MEST", "ZHX1-C8ORF76", 
                                 "MTURN", "CAPZA2", "MS4A1", "LINC00926", "CD79B", "ISG20", "PKIG"
  ))

})

test_that("Checking clust_size()", {
  expect_equal(clust_size(res), c(`1` = 51L, `2` = 49L, `3` = 45L, `4` = 24L, `5` = 20L, `6` = 18L, 
                                  `7` = 14L, `8` = 14L, `9` = 14L, `10` = 12L, `11` = 7L, `12` = 7L, 
                                  `13` = 6L, `14` = 5L, `15` = 5L))

})

test_that("Checking gene_cluster()", {
  expect_equal(unname(gene_cluster(res)), c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
                                            1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
                                            1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
                                            1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 
                                            2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 
                                            2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 
                                            2L, 2L, 2L, 2L, 2L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 
                                            3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 
                                            3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 
                                            3L, 3L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 
                                            4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 5L, 5L, 5L, 5L, 5L, 5L, 
                                            5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 6L, 6L, 
                                            6L, 6L, 6L, 6L, 6L, 6L, 6L, 6L, 6L, 6L, 6L, 6L, 6L, 6L, 6L, 6L, 
                                            7L, 7L, 7L, 7L, 7L, 7L, 7L, 7L, 7L, 7L, 7L, 7L, 7L, 7L, 8L, 8L, 
                                            8L, 8L, 8L, 8L, 8L, 8L, 8L, 8L, 8L, 8L, 8L, 8L, 9L, 9L, 9L, 9L, 
                                            9L, 9L, 9L, 9L, 9L, 9L, 9L, 9L, 9L, 9L, 10L, 10L, 10L, 10L, 10L, 
                                            10L, 10L, 10L, 10L, 10L, 10L, 10L, 11L, 11L, 11L, 11L, 11L, 11L, 
                                            11L, 12L, 12L, 12L, 12L, 12L, 12L, 12L, 13L, 13L, 13L, 13L, 13L, 
                                            13L, 14L, 14L, 14L, 14L, 14L, 15L, 15L, 15L, 15L, 15L))

})

# Set verbosity to 0
set_verbosity(0)

load_example_dataset("7871581/files/pbmc3k_medium_clusters")

res <- pbmc3k_medium_clusters

test_that("Checking if get_genes stop if no slot top_genes in ClusterSet obj", {
  expect_error(get_genes(res, cluster = "all", top = TRUE))
})

test_that("Checking get_genes is providing the right list of genes", {
  
  # =======================================================
  # Test gene list in all cluster
  gene_names <- get_genes(res, cluster = "all")
  gene_name_to_check <- c("MALAT1", "RPL9", "RPL23A", "RPL21", "RPS25", "RPS19", "RPS3", 
                          "RPS12", "RPL18A", "RPL11", "RPL32", "RPLP2", "RPS18", "RPL3", 
                          "RPL29", "RPL14", "RPS3A", "RPS15A", "RPS23", "RPL27A", "RPS16", 
                          "RPL12", "EEF1A1", "RPL31", "RPS4X", "RPS8", "RPL28", "RPS29", 
                          "RPS10", "RPS13", "RPS11", "RPS26", "RPL35A", "TPT1", "BTG1", 
                          "RPL15", "RPL18", "RPS20", "RPS28", "RPL36", "RPL5", "EEF1B2", 
                          "RPL7A", "RPL17", "NACA", "RPL27", "RPL26", "RPS9", "RPLP1", 
                          "RPSAP58", "RPL6", "PTCRA", "ACRBP", "SDPR", "GNG11", "PF4", 
                          "HIST1H2AC", "TMEM40", "TUBB1", "SPARC", "PPBP", "GP9", "MMD", 
                          "AP001189.4", "CLU", "RGS18", "MPP1", "CA2", "F13A1", "TPM1", 
                          "MYL9", "TREML1", "NGFRAP1", "PTGS1", "GAS2L1", "GRAP2", "ODC1", 
                          "SNCA", "TPM4", "CTSA", "ILK", "PDLIM1", "ACTN1", "FERMT3", "PLA2G12A", 
                          "HIST1H2BK", "SNN", "APP", "YWHAH", "ASAH1", "RAP1B", "MFSD1", 
                          "GSN", "MARCH2", "H2AFJ", "TPTEP1", "PARVB", "MYL12A", "CD151", 
                          "PVALB", "TYMP", "COTL1", "CTSS", "CFD", "S100A11", "LGALS2", 
                          "GABARAP", "APOBEC3A", "TNFSF13B", "CPVL", "MNDA", "CD14", "IGSF6", 
                          "CSTB", "MAFB", "AP2S1", "WARS", "IFI6", "LILRB2", "LILRB4", 
                          "C1orf162", "SLC7A7", "IFNGR2", "SOD2", "TGFBI", "BLVRA", "CAPG", 
                          "BID", "CTSZ", "NAAA", "FCGR2A", "CD300LF", "S100A10", "GAPDH", 
                          "ALDH2", "EIF4EBP1", "GRINA", "NUP214", "CD302", "ATG3", "CAMK1", 
                          "CYBB", "CNPY3", "GCA", "POU2F2", "CCL5", "GZMA", "CST7", "CTSW", 
                          "HLA-C", "FGFBP2", "GZMB", "GNLY", "PRF1", "SPON2", "CCL4", "PTPRCAP", 
                          "CLIC3", "IL2RB", "HOPX", "CD247", "AKR1C3", "IGFBP7", "KLRF1", 
                          "CD99", "CD7", "XCL1", "EFHD2", "PLEK", "HEMGN", "EGFL7", "AC137932.6", 
                          "hsa-mir-1199", "SEC14L5", "HEXIM2", "GATA2", "SH3BGRL2", "KIAA0513", 
                          "PCP2", "PYGL", "DPY19L1", "RAB27B", "CLIC4", "LYL1", "FAM110A", 
                          "RILP", "DAB2", "KIFC3", "BTK", "SPOCD1", "ARHGAP21", "ITGB3", 
                          "ALOX12", "SCFD2", "ATP9A", "FAM212A", "ABCC3", "RP11-359I18.5", 
                          "TGFB1I1", "CTNNAL1", "ENDOD1", "NEXN", "ABHD16A", "ARG2", "GLA", 
                          "WIPI1", "PPP1R14A", "CLEC1B", "CMTM5", "GP1BA", "RP11-367G6.3", 
                          "HIST1H2BD", "FNTB", "SENCR", "SEPT4", "MLH3", "SCGB1C1", "HIST2H2BE", 
                          "C1orf198", "H1F0", "GLUL", "LGALSL", "SMOX", "SSX2IP", "AC147651.3", 
                          "C19orf33", "FAM212B", "SERPINE2", "ZNF778", "TMCC2", "TRIM58", 
                          "FN3K", "NFE2", "MGLL", "MOB3B", "HLA-DPA1", "HLA-DRB1", "HLA-DRA", 
                          "HLA-DRB5", "HLA-DPB1", "CD74", "HLA-DQA1", "HLA-DQB1", "FAM26F", 
                          "LY86", "HLA-DMB", "SLC25A6", "CD37", "RNASET2", "C2orf88", "CLDN5", 
                          "HIST1H2BJ", "FAM63A", "SPHK1", "SLC40A1", "TUBA1C", "TMEM91", 
                          "XPNPEP1", "PRKAR2B", "TJP2", "AC092295.7", "ENKUR", "LCN2", 
                          "GFI1B", "ANKRD9", "PRUNE", "SMIM5", "PTPN18", "CD3D", "IL32", 
                          "CD3E", "LDHB", "IL7R", "AES", "CD27", "HES4", "FCGR3A", "RHOC", 
                          "IFITM2", "ABI3", "CLIC1", "RP11-879F14.2", "MEST", "ZHX1-C8ORF76", 
                          "MTURN", "CAPZA2", "MS4A1", "LINC00926", "CD79B", "ISG20", "PKIG"
  )
  
  
  expect_equal(gene_names, gene_name_to_check)
  expect_equal(gene_names, rownames(res@data))
  expect_equal(gene_names, unlist(res@gene_clusters, use.names = FALSE))

  # =======================================================
  # Test gene list in cluster 1
  gene_names <- get_genes(res, cluster = 1)
  gene_name_to_check <- c("MALAT1", "RPL9", "RPL23A", "RPL21", "RPS25", "RPS19", "RPS3", 
                          "RPS12", "RPL18A", "RPL11", "RPL32", "RPLP2", "RPS18", "RPL3", 
                          "RPL29", "RPL14", "RPS3A", "RPS15A", "RPS23", "RPL27A", "RPS16", 
                          "RPL12", "EEF1A1", "RPL31", "RPS4X", "RPS8", "RPL28", "RPS29", 
                          "RPS10", "RPS13", "RPS11", "RPS26", "RPL35A", "TPT1", "BTG1", 
                          "RPL15", "RPL18", "RPS20", "RPS28", "RPL36", "RPL5", "EEF1B2", 
                          "RPL7A", "RPL17", "NACA", "RPL27", "RPL26", "RPS9", "RPLP1", 
                          "RPSAP58", "RPL6")
  expect_equal(gene_names, gene_name_to_check)
  expect_equal(gene_names, res@gene_clusters$`1`)
  
  # =======================================================
  # Test gene list in cluster 2
  gene_names <- get_genes(res, cluster = 2)
  gene_name_to_check <- c("PTCRA", "ACRBP", "SDPR", "GNG11", "PF4", "HIST1H2AC", "TMEM40", 
                          "TUBB1", "SPARC", "PPBP", "GP9", "MMD", "AP001189.4", "CLU", 
                          "RGS18", "MPP1", "CA2", "F13A1", "TPM1", "MYL9", "TREML1", "NGFRAP1", 
                          "PTGS1", "GAS2L1", "GRAP2", "ODC1", "SNCA", "TPM4", "CTSA", "ILK", 
                          "PDLIM1", "ACTN1", "FERMT3", "PLA2G12A", "HIST1H2BK", "SNN", 
                          "APP", "YWHAH", "ASAH1", "RAP1B", "MFSD1", "GSN", "MARCH2", "H2AFJ", 
                          "TPTEP1", "PARVB", "MYL12A", "CD151", "PVALB")
  expect_equal(gene_names, gene_name_to_check)
  expect_equal(gene_names, res@gene_clusters$`2`)
  
  # =======================================================
  # Test gene list in cluster 3
  gene_names <- get_genes(res, cluster = 3)
  gene_name_to_check <- c("TYMP", "COTL1", "CTSS", "CFD", "S100A11", "LGALS2", "GABARAP", 
                          "APOBEC3A", "TNFSF13B", "CPVL", "MNDA", "CD14", "IGSF6", "CSTB", 
                          "MAFB", "AP2S1", "WARS", "IFI6", "LILRB2", "LILRB4", "C1orf162", 
                          "SLC7A7", "IFNGR2", "SOD2", "TGFBI", "BLVRA", "CAPG", "BID", 
                          "CTSZ", "NAAA", "FCGR2A", "CD300LF", "S100A10", "GAPDH", "ALDH2", 
                          "EIF4EBP1", "GRINA", "NUP214", "CD302", "ATG3", "CAMK1", "CYBB", 
                          "CNPY3", "GCA", "POU2F2")
  expect_equal(gene_names, gene_name_to_check)
  expect_equal(gene_names, res@gene_clusters$`3`)
})








test_that("Checking get_genes using top genes.", {
  res <- top_genes(res, top = 20, cluster = "all")
  
  # =======================================================
  # Test top gene list in all cluster
  gene_names <- get_genes(res, cluster = "all", top = T)
  gene_name_to_check <- c("RPL11", "RPL32", "RPS12", "RPL18A", "RPS18", "RPS23", "RPS16", 
                          "RPL9", "RPL23A", "RPS3", "RPL12", "RPS15A", "RPLP2", "RPL28", 
                          "RPS19", "RPS25", "EEF1A1", "RPS8", "RPS4X", "RPL29", "SDPR", 
                          "GNG11", "PF4", "SPARC", "PPBP", "TUBB1", "HIST1H2AC", "GP9", 
                          "CLU", "AP001189.4", "RGS18", "PTCRA", "TREML1", "F13A1", "CA2", 
                          "TMEM40", "MMD", "MYL9", "PTGS1", "MPP1", "TYMP", "COTL1", "CTSS", 
                          "S100A11", "CFD", "LGALS2", "GABARAP", "TNFSF13B", "APOBEC3A", 
                          "AP2S1", "IGSF6", "CSTB", "CPVL", "MAFB", "IFI6", "SLC7A7", "TGFBI", 
                          "C1orf162", "CAPG", "BLVRA", "GZMB", "PRF1", "GNLY", "FGFBP2", 
                          "GZMA", "CST7", "CCL4", "SPON2", "CTSW", "CLIC3", "CD247", "IL2RB", 
                          "HOPX", "HLA-C", "CD7", "AKR1C3", "KLRF1", "XCL1", "IGFBP7", 
                          "CCL5", "AC137932.6", "HEMGN", "hsa-mir-1199", "EGFL7", "SEC14L5", 
                          "GATA2", "DPY19L1", "HEXIM2", "PCP2", "KIFC3", "RAB27B", "SH3BGRL2", 
                          "KIAA0513", "DAB2", "RILP", "CLIC4", "PYGL", "FAM110A", "BTK", 
                          "LYL1", "ARHGAP21", "SPOCD1", "ATP9A", "FAM212A", "CTNNAL1", 
                          "ABCC3", "SCFD2", "ARG2", "ITGB3", "ALOX12", "ENDOD1", "NEXN", 
                          "RP11-359I18.5", "TGFB1I1", "ABHD16A", "GLA", "WIPI1", "PPP1R14A", 
                          "CLEC1B", "GP1BA", "CMTM5", "SCGB1C1", "HIST1H2BD", "RP11-367G6.3", 
                          "FNTB", "SEPT4", "SENCR", "C1orf198", "MLH3", "HIST2H2BE", "GLUL", 
                          "H1F0", "LGALSL", "SMOX", "SSX2IP", "FAM212B", "AC147651.3", 
                          "TMCC2", "ZNF778", "SERPINE2", "TRIM58", "MOB3B", "NFE2", "C19orf33", 
                          "MGLL", "FN3K", "HLA-DRB1", "HLA-DRA", "HLA-DPA1", "HLA-DRB5", 
                          "HLA-DPB1", "CD74", "HLA-DQA1", "HLA-DQB1", "LY86", "HLA-DMB", 
                          "FAM26F", "CD37", "RNASET2", "SLC25A6", "HIST1H2BJ", "CLDN5", 
                          "C2orf88", "FAM63A", "TUBA1C", "SLC40A1", "SPHK1", "TMEM91", 
                          "PRKAR2B", "TJP2", "AC092295.7", "XPNPEP1", "ENKUR", "LCN2", 
                          "GFI1B", "ANKRD9", "PRUNE", "SMIM5", "PTPN18", "CD3D", "IL32", 
                          "LDHB", "CD3E", "IL7R", "CD27", "AES", "FCGR3A", "RHOC", "IFITM2", 
                          "ABI3", "HES4", "CLIC1", "RP11-879F14.2", "ZHX1-C8ORF76", "MTURN", 
                          "MEST", "CAPZA2", "LINC00926", "CD79B", "MS4A1", "ISG20", "PKIG"
  )
  expect_equal(gene_names, gene_name_to_check)
  expect_equal(sort(gene_names), sort(unlist(res@top_genes, use.names = F)))
  
  
  # =======================================================
  # Test top gene list in cluster 1
  gene_names <- get_genes(res, cluster = 1, top = T)
  gene_name_to_check <- c("RPL11", "RPL32", "RPS12", "RPL18A", "RPS18", "RPS23", "RPS16", 
                          "RPL9", "RPL23A", "RPS3", "RPL12", "RPS15A", "RPLP2", "RPL28", 
                          "RPS19", "RPS25", "EEF1A1", "RPS8", "RPS4X", "RPL29")
  expect_equal(gene_names, gene_name_to_check)
  expect_equal(sort(gene_names), sort(res@top_genes$`1`))
  
  # =======================================================
  # Test top gene list in cluster 2 and 3
  gene_names <- get_genes(res, cluster = 2:3, top = T)
  gene_name_to_check <- c("SDPR", "GNG11", "PF4", "SPARC", "PPBP", "TUBB1", "HIST1H2AC", 
                          "GP9", "CLU", "AP001189.4", "RGS18", "PTCRA", "TREML1", "F13A1", 
                          "CA2", "TMEM40", "MMD", "MYL9", "PTGS1", "MPP1", "TYMP", "COTL1", 
                          "CTSS", "S100A11", "CFD", "LGALS2", "GABARAP", "TNFSF13B", "APOBEC3A", 
                          "AP2S1", "IGSF6", "CSTB", "CPVL", "MAFB", "IFI6", "SLC7A7", "TGFBI", 
                          "C1orf162", "CAPG", "BLVRA")
  expect_equal(gene_names, gene_name_to_check)
  expect_equal(sort(gene_names), sort(c(
    res@top_genes$`2`,
    res@top_genes$`3`
  )))
})

set_verbosity(0)
load_example_dataset('7871581/files/pbmc3k_medium_clusters')

set_verbosity(0)
load_example_dataset('7871581/files/pbmc3k_medium_clusters')

test_that("Check cluster_set_to_xls", {
  
  tp_dir <- file.path(tempdir(), create_rand_str())
  dir.create(tp_dir, showWarnings = F, recursive = TRUE)
  f_path <- file.path(tp_dir, "test.xls")
  unlink(f_path)
  cluster_set_to_xls(pbmc3k_medium_clusters, f_path)
  testthat::expect_true(file.exists(f_path))
  testthat::expect_error(cluster_set_to_xls(pbmc3k_medium_clusters, f_path))
  
  tp_dir <- file.path(tempdir(), create_rand_str())
  dir.create(tp_dir, showWarnings = F, recursive = TRUE)
  f_path <- file.path(tp_dir, "test.xls")
  unlink(f_path)
  cluster_set_to_xls(pbmc3k_medium_clusters, f_path)
  testthat::expect_true(file.exists(f_path))
  testthat::expect_error(cluster_set_to_xls(pbmc3k_medium_clusters, f_path))
})

test_that("Checking gene_cluster() is working...", { 
  
  # Set verbosity to 0
  set_verbosity(0)
  
  load_example_dataset("7871581/files/pbmc3k_medium_clusters")
  res <- pbmc3k_medium_clusters
  
  expect_equal(length(gene_cluster(res)), 291)
  expect_equal(paste0(table(gene_cluster(res)), collapse = " "), "51 49 45 24 20 18 14 14 14 12 7 7 6 5 5")
  expect_equal(as.vector(table(gene_cluster(res, cluster = 1))), 51)
  expect_equal(as.vector(table(gene_cluster(res, cluster = 2))), 49)
  expect_equal(as.vector(table(gene_cluster(res, cluster = 3))), 45)
  expect_equal(as.vector(table(gene_cluster(res, cluster = 4))), 24)
  expect_equal(as.vector(table(gene_cluster(res, cluster = c(1,4)))), c(51, 24))
  expect_equal(as.vector(table(gene_cluster(res, cluster = c(4,4)))), 24)
  expect_error(gene_cluster(res, cluster = c(-1,8)))
  expect_error(gene_cluster(res, cluster = c(0,9)))
  expect_error(gene_cluster(res, cluster = c(0:8)))
  expect_error(gene_cluster(res, cluster = c(1,40)))
  expect_equal(paste0(res@gene_clusters[[1]], collapse = " "), 
               paste0(names(gene_cluster(res, cluster = 1)), collapse = " "))
  expect_equal(paste0(unlist(res@gene_clusters[1:2]), collapse = " "), 
               paste0(names(gene_cluster(res, cluster = 1:2)), collapse = " "))
  expect_equal(paste0(unlist(res@gene_clusters[1:4]), collapse = " "), 
               paste0(names(gene_cluster(res, cluster = 1:4)), collapse = " "))
})
