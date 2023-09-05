library(testthat)
library(iheatmapr)
library(ggplot2)
library(Seurat)

# Set verbosity to 0
set_verbosity(0)
load_example_dataset("7871581/files/pbmc3k_medium")
load_example_dataset("7871581/files/pbmc3k_medium_clusters")

res <- pbmc3k_medium_clusters

test_that("Checking plot_heatmap()", {
  set.seed(123)
  
  htmp <- plot_heatmap(
    object = res,
    cell_clusters = NULL,
    use_top_genes = FALSE,
    show_dendro = TRUE,
    line_size_vertical = 2,
    line_size_horizontal = 2,
    name = "Expression heatmap",
    show_legend = TRUE,
    xlab = "Samples",
    ylab = "Genes"
  )
  
  htmp_matrix <- htmp@plots@listData$`Exp. level`@data
  htmp_matrix_wo_na <- htmp_matrix[!is.na(htmp_matrix)]
  
  # Check heatmap format

  # Check dim of heatmap matrix
  expect_equal(ncol(htmp_matrix), 361)
  expect_equal(nrow(htmp_matrix), 319)
  expect_equal(dim(htmp_matrix), c(319, 361))
  
  # Check NA rows added by line_size_horizontal
  na_rows <- which(apply(htmp_matrix, 1, function(x) all(is.na(x))))
  expect_equal(length(na_rows), 28)
  expect_equal(as.vector(na_rows), c(6L, 7L, 13L, 14L, 21L, 22L, 30L, 31L, 39L, 40L, 53L, 54L, 69L, 
                                     70L, 85L, 86L, 101L, 102L, 121L, 122L, 143L, 144L, 169L, 170L, 
                                     216L, 217L, 267L, 268L))
  
  # Check NA rows added by line_size_vertical
  na_cols <- which(apply(htmp_matrix, 2, function(x) all(is.na(x))))
  expect_equal(length(na_cols), 0)
  expect_equal(as.vector(na_cols), as.integer())
  
  # Check some measurements of heatmap matrix
  expect_equal(round(sum(htmp_matrix_wo_na), 2), -3631.17)
  expect_equal(round(mean(htmp_matrix_wo_na), 2), -0.03)
  expect_equal(round(median(htmp_matrix_wo_na), 2), -0.07)
  expect_equal(round(sd(htmp_matrix_wo_na), 2), 0.49)
  
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["0%"]), 2), -1)
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["25%"]), 2), -0.25)
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["50%"]), 2), -0.07)
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["75%"]), 2), -0.02)
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["100%"]), 2), 1)
  
  # Checking order of columns
  expect_equal(head(colnames(htmp_matrix),10), c("TAACTCACGTATCG-1", "CTAATGCTTGTGGT-1", "TTCGAGGATAGAAG-1", 
                                                 "TTCTTACTCTGGAT-1", "ATAAGTTGGTACGT-1", "TTACCATGAATCGC-1", "GGCGCATGCCTAAG-1", 
                                                 "ATACGGACCTACTT-1", "GGAACACTCACTTT-1", "AAGGTCACGGTTAC-1"))
  
  # Checking order of columns
  expect_equal(tail(colnames(htmp_matrix),10), c("TTCGTATGAAAAGC-1", "CCAAGAACTACTGG-1", "TCGAATCTCTGGTA-1", 
                                                 "GTCATACTGCGATT-1", "GGCAATACGGCATT-1", "GCTATACTAGCGTT-1", "CAGCTCTGCAAGCT-1", 
                                                 "TCAATCACACTCTT-1", "ATCCTAACGCTACA-1", "ACGATTCTACGGGA-1"))
  
  # Checking order of rows
  expect_equal(rownames(htmp_matrix), c("PKIG", "ISG20", "CD79B", "LINC00926", "MS4A1", " ", " ", "CAPZA2", 
                                        "MTURN", "ZHX1-C8ORF76", "MEST", "RP11-879F14.2", " ", " ", "CLIC1", 
                                        "ABI3", "IFITM2", "RHOC", "FCGR3A", "HES4", " ", " ", "CD27", 
                                        "AES", "IL7R", "LDHB", "CD3E", "IL32", "CD3D", " ", " ", "PTPN18", 
                                        "SMIM5", "PRUNE", "ANKRD9", "GFI1B", "LCN2", "ENKUR", " ", " ", 
                                        "AC092295.7", "TJP2", "PRKAR2B", "XPNPEP1", "TMEM91", "TUBA1C", 
                                        "SLC40A1", "SPHK1", "FAM63A", "HIST1H2BJ", "CLDN5", "C2orf88", 
                                        " ", " ", "RNASET2", "CD37", "SLC25A6", "HLA-DMB", "LY86", "FAM26F", 
                                        "HLA-DQB1", "HLA-DQA1", "CD74", "HLA-DPB1", "HLA-DRB5", "HLA-DRA", 
                                        "HLA-DRB1", "HLA-DPA1", " ", " ", "MOB3B", "MGLL", "NFE2", "FN3K", 
                                        "TRIM58", "TMCC2", "ZNF778", "SERPINE2", "FAM212B", "C19orf33", 
                                        "AC147651.3", "SSX2IP", "SMOX", "LGALSL", " ", " ", "GLUL", "H1F0", 
                                        "C1orf198", "HIST2H2BE", "SCGB1C1", "MLH3", "SEPT4", "SENCR", 
                                        "FNTB", "HIST1H2BD", "RP11-367G6.3", "GP1BA", "CMTM5", "CLEC1B", 
                                        " ", " ", "PPP1R14A", "WIPI1", "GLA", "ARG2", "ABHD16A", "NEXN", 
                                        "ENDOD1", "CTNNAL1", "TGFB1I1", "RP11-359I18.5", "ABCC3", "FAM212A", 
                                        "ATP9A", "SCFD2", "ALOX12", "ITGB3", "ARHGAP21", "SPOCD1", " ", 
                                        " ", "BTK", "KIFC3", "DAB2", "RILP", "FAM110A", "LYL1", "CLIC4", 
                                        "RAB27B", "DPY19L1", "PYGL", "PCP2", "KIAA0513", "SH3BGRL2", 
                                        "GATA2", "HEXIM2", "SEC14L5", "hsa-mir-1199", "AC137932.6", "EGFL7", 
                                        "HEMGN", " ", " ", "PLEK", "EFHD2", "XCL1", "CD7", "CD99", "KLRF1", 
                                        "IGFBP7", "AKR1C3", "CD247", "HOPX", "IL2RB", "CLIC3", "PTPRCAP", 
                                        "CCL4", "SPON2", "PRF1", "GNLY", "GZMB", "FGFBP2", "HLA-C", "CTSW", 
                                        "CST7", "GZMA", "CCL5", " ", " ", "POU2F2", "GCA", "CNPY3", "CYBB", 
                                        "CAMK1", "ATG3", "CD302", "NUP214", "GRINA", "EIF4EBP1", "ALDH2", 
                                        "GAPDH", "S100A10", "CD300LF", "FCGR2A", "NAAA", "CTSZ", "BID", 
                                        "CAPG", "BLVRA", "TGFBI", "SOD2", "IFNGR2", "SLC7A7", "C1orf162", 
                                        "LILRB4", "LILRB2", "IFI6", "WARS", "AP2S1", "MAFB", "CSTB", 
                                        "IGSF6", "CD14", "MNDA", "CPVL", "TNFSF13B", "APOBEC3A", "GABARAP", 
                                        "LGALS2", "S100A11", "CFD", "CTSS", "COTL1", "TYMP", " ", " ", 
                                        "PVALB", "CD151", "MYL12A", "PARVB", "TPTEP1", "H2AFJ", "MARCH2", 
                                        "GSN", "MFSD1", "RAP1B", "ASAH1", "YWHAH", "APP", "SNN", "HIST1H2BK", 
                                        "PLA2G12A", "FERMT3", "ACTN1", "PDLIM1", "ILK", "CTSA", "TPM4", 
                                        "SNCA", "ODC1", "GRAP2", "GAS2L1", "PTGS1", "NGFRAP1", "TREML1", 
                                        "MYL9", "TPM1", "F13A1", "CA2", "MPP1", "RGS18", "CLU", "AP001189.4", 
                                        "MMD", "GP9", "PPBP", "SPARC", "TUBB1", "TMEM40", "HIST1H2AC", 
                                        "PF4", "GNG11", "SDPR", "ACRBP", "PTCRA", " ", " ", "RPL6", "RPSAP58", 
                                        "RPLP1", "RPS9", "RPL26", "RPL27", "NACA", "RPL17", "RPL7A", 
                                        "EEF1B2", "RPL5", "RPL36", "RPS28", "RPS20", "RPL18", "RPL15", 
                                        "BTG1", "TPT1", "RPL35A", "RPS26", "RPS11", "RPS13", "RPS10", 
                                        "RPS29", "RPL28", "RPS8", "RPS4X", "RPL31", "EEF1A1", "RPL12", 
                                        "RPS16", "RPL27A", "RPS23", "RPS15A", "RPS3A", "RPL14", "RPL29", 
                                        "RPL3", "RPS18", "RPLP2", "RPL32", "RPL11", "RPL18A", "RPS12", 
                                        "RPS3", "RPS19", "RPS25", "RPL21", "RPL23A", "RPL9", "MALAT1"
  ))
  
  # Checking title of the heatmap
  htmp_title <- htmp@annotations@listData$col_title2@data
  expect_equal(htmp_title, "Expression heatmap")
  
  # Checking x lab of the heatmap
  htmp_xlab <- htmp@annotations@listData$col_title@data
  expect_equal(htmp_xlab, "Samples")
  
  # Checking y lab of the heatmap
  htmp_ylab <- htmp@annotations@listData$row_title@data
  expect_equal(htmp_ylab, "Genes")
  
  # Checking colorbar name of the heatmap
  htmp_colors_name <- htmp@colorbars@listData$`Exp. level`@title
  expect_equal(htmp_colors_name, "Exp. level")
  
  # Checking if colorbar is shown of the heatmap
  htmp_show_colorbar <- htmp@plots@listData$`Exp. level`@show_colorbar
  expect_true(htmp_show_colorbar)
  
  # Checking if dendrogram is showed
  htmp_show_dendro <- names(htmp@shapes@listData)
  expect_equal(htmp_show_dendro, "col_dendro")
})

test_that("Checking plot_heatmap() using top genes", {
  set.seed(123)
  
  res_scigenex_top <- top_genes(res, top = 20)
  htmp <- plot_heatmap(
    object = res_scigenex_top,
    cell_clusters = NULL,
    use_top_genes = TRUE, # Use top genes
    show_dendro = FALSE, # Set to FALSE
    line_size_vertical = 2,
    line_size_horizontal = 2,
    name = NULL, # Set to NULL
    show_legend = FALSE, # Set to FALSE
    xlab = NULL, # Set to NULL
    ylab = NULL
  ) # Set to NULL
  
  htmp_matrix <- htmp@plots@listData$`Exp. level`@data
  htmp_matrix_wo_na <- htmp_matrix[!is.na(htmp_matrix)]
  
  # Check heatmap format

  # Check dim of heatmap matrix
  expect_equal(ncol(htmp_matrix), 361)
  expect_equal(nrow(htmp_matrix), 230)
  expect_equal(dim(htmp_matrix), c(230, 361))
  
  # Check NA rows added by line_size_horizontal
  na_rows <- which(apply(htmp_matrix, 1, function(x) all(is.na(x))))
  expect_equal(length(na_rows), 28)
  expect_equal(as.vector(na_rows), c(6L, 7L, 13L, 14L, 21L, 22L, 30L, 31L, 39L, 40L, 53L, 54L, 69L, 
                                     70L, 85L, 86L, 101L, 102L, 121L, 122L, 143L, 144L, 165L, 166L, 
                                     187L, 188L, 209L, 210L))
  
  # Check NA rows added by line_size_horizontal
  na_cols <- which(apply(htmp_matrix, 2, function(x) all(is.na(x))))
  expect_equal(length(na_cols), 0)
  expect_equal(as.vector(na_cols), as.integer())
  
  # Check some measurements of heatmap matrix
  expect_equal(round(sum(htmp_matrix_wo_na), 2), -3394.67)
  expect_equal(round(mean(htmp_matrix_wo_na), 2), -0.05)
  expect_equal(round(median(htmp_matrix_wo_na), 2), -0.05)
  expect_equal(round(sd(htmp_matrix_wo_na), 2), 0.48)
  
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["0%"]), 2), -1)
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["25%"]), 2), -0.22)
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["50%"]), 2), -0.05)
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["75%"]), 2), -0.02)
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["100%"]), 2), 1)
  
  # Checking order of columns
  expect_equal(head(colnames(htmp_matrix),10), c("TAACTCACGTATCG-1", "CTAATGCTTGTGGT-1", "TTCGAGGATAGAAG-1", 
                                                 "TTCTTACTCTGGAT-1", "ATAAGTTGGTACGT-1", "TTACCATGAATCGC-1", "GGCGCATGCCTAAG-1", 
                                                 "ATACGGACCTACTT-1", "GGAACACTCACTTT-1", "AAGGTCACGGTTAC-1"))
  
  expect_equal(tail(colnames(htmp_matrix),10), c("TTCGTATGAAAAGC-1", "CCAAGAACTACTGG-1", "TCGAATCTCTGGTA-1", 
                                                 "GTCATACTGCGATT-1", "GGCAATACGGCATT-1", "GCTATACTAGCGTT-1", "CAGCTCTGCAAGCT-1", 
                                                 "TCAATCACACTCTT-1", "ATCCTAACGCTACA-1", "ACGATTCTACGGGA-1"))  
  # Checking order of rows
  expect_equal(rownames(htmp_matrix), c("PKIG", "ISG20", "MS4A1", "CD79B", "LINC00926", " ", " ", "CAPZA2", 
                                        "MEST", "MTURN", "ZHX1-C8ORF76", "RP11-879F14.2", " ", " ", "CLIC1", 
                                        "HES4", "ABI3", "IFITM2", "RHOC", "FCGR3A", " ", " ", "AES", 
                                        "CD27", "IL7R", "CD3E", "LDHB", "IL32", "CD3D", " ", " ", "PTPN18", 
                                        "SMIM5", "PRUNE", "ANKRD9", "GFI1B", "LCN2", "ENKUR", " ", " ", 
                                        "XPNPEP1", "AC092295.7", "TJP2", "PRKAR2B", "TMEM91", "SPHK1", 
                                        "SLC40A1", "TUBA1C", "FAM63A", "C2orf88", "CLDN5", "HIST1H2BJ", 
                                        " ", " ", "SLC25A6", "RNASET2", "CD37", "FAM26F", "HLA-DMB", 
                                        "LY86", "HLA-DQB1", "HLA-DQA1", "CD74", "HLA-DPB1", "HLA-DRB5", 
                                        "HLA-DPA1", "HLA-DRA", "HLA-DRB1", " ", " ", "FN3K", "MGLL", 
                                        "C19orf33", "NFE2", "MOB3B", "TRIM58", "SERPINE2", "ZNF778", 
                                        "TMCC2", "AC147651.3", "FAM212B", "SSX2IP", "SMOX", "LGALSL", 
                                        " ", " ", "H1F0", "GLUL", "HIST2H2BE", "MLH3", "C1orf198", "SENCR", 
                                        "SEPT4", "FNTB", "RP11-367G6.3", "HIST1H2BD", "SCGB1C1", "CMTM5", 
                                        "GP1BA", "CLEC1B", " ", " ", "PPP1R14A", "WIPI1", "GLA", "ABHD16A", 
                                        "TGFB1I1", "RP11-359I18.5", "NEXN", "ENDOD1", "ALOX12", "ITGB3", 
                                        "ARG2", "SCFD2", "ABCC3", "CTNNAL1", "FAM212A", "ATP9A", "SPOCD1", 
                                        "ARHGAP21", " ", " ", "LYL1", "BTK", "FAM110A", "PYGL", "CLIC4", 
                                        "RILP", "DAB2", "KIAA0513", "SH3BGRL2", "RAB27B", "KIFC3", "PCP2", 
                                        "HEXIM2", "DPY19L1", "GATA2", "SEC14L5", "EGFL7", "hsa-mir-1199", 
                                        "HEMGN", "AC137932.6", " ", " ", "CCL5", "IGFBP7", "XCL1", "KLRF1", 
                                        "AKR1C3", "CD7", "HLA-C", "HOPX", "IL2RB", "CD247", "CLIC3", 
                                        "CTSW", "SPON2", "CCL4", "CST7", "GZMA", "FGFBP2", "GNLY", "PRF1", 
                                        "GZMB", " ", " ", "BLVRA", "CAPG", "C1orf162", "TGFBI", "SLC7A7", 
                                        "IFI6", "MAFB", "CPVL", "CSTB", "IGSF6", "AP2S1", "APOBEC3A", 
                                        "TNFSF13B", "GABARAP", "LGALS2", "CFD", "S100A11", "CTSS", "COTL1", 
                                        "TYMP", " ", " ", "MPP1", "PTGS1", "MYL9", "MMD", "TMEM40", "CA2", 
                                        "F13A1", "TREML1", "PTCRA", "RGS18", "AP001189.4", "CLU", "GP9", 
                                        "HIST1H2AC", "TUBB1", "PPBP", "SPARC", "PF4", "GNG11", "SDPR", 
                                        " ", " ", "RPL29", "RPS4X", "RPS8", "EEF1A1", "RPS25", "RPS19", 
                                        "RPL28", "RPLP2", "RPS15A", "RPL12", "RPS3", "RPL23A", "RPL9", 
                                        "RPS16", "RPS23", "RPS18", "RPL18A", "RPS12", "RPL32", "RPL11"
  ))
  
  # Checking title of the heatmap
  htmp_slots <- names(htmp@annotations@listData)
  expect_false("col_title2" %in% htmp_slots)
  
  # Checking x lab of the heatmap
  expect_false("col_title" %in% htmp_slots)
  
  # Checking y lab of the heatmap
  expect_false("row_title" %in% htmp_slots)
  
  # Checking colorbar name of the heatmap
  htmp_colors_name <- htmp@colorbars@listData$`Exp. level`@title
  expect_equal(htmp_colors_name, "Exp. level")
  
  # Checking if colorbar is showed
  htmp_show_colorbar <- htmp@plots@listData$`Exp. level`@show_colorbar
  expect_false(htmp_show_colorbar)
  
  # Checking if dendrogram is showed
  htmp_show_dendro <- names(htmp@shapes@listData)
  expect_equal(htmp_show_dendro, NULL)
})

test_that("Checking plot_heatmap() using cell clusters", {
  
  set.seed(123)
  cell_clusters <- Idents(pbmc3k_medium)
  htmp <- plot_heatmap(
    object = res,
    cell_clusters = cell_clusters,
    use_top_genes = FALSE,
    show_dendro = TRUE,
    line_size_vertical = 2,
    line_size_horizontal = 2,
    name = "Expression heatmap",
    show_legend = TRUE,
    xlab = "Samples",
    ylab = "Genes"
  )
  
  htmp_matrix <- htmp@plots@listData$`Exp. level`@data
  htmp_matrix_wo_na <- htmp_matrix[!is.na(htmp_matrix)]
  
  # Check heatmap format

  # Check dim of heatmap matrix
  expect_equal(ncol(htmp_matrix), 377)
  expect_equal(nrow(htmp_matrix), 319)
  expect_equal(dim(htmp_matrix), c(319, 377))
  
  # Check NA rows added by line_size_horizontal
  na_rows <- which(apply(htmp_matrix, 1, function(x) all(is.na(x))))
  expect_equal(length(na_rows), 28)
  expect_equal(as.vector(na_rows), c(6L, 7L, 13L, 14L, 21L, 22L, 30L, 31L, 39L, 40L, 53L, 54L, 69L, 
                                     70L, 85L, 86L, 101L, 102L, 121L, 122L, 143L, 144L, 169L, 170L, 
                                     216L, 217L, 267L, 268L))
  
  # Check NA rows added by line_size_vertical
  na_cols <- which(apply(htmp_matrix, 2, function(x) all(is.na(x))))
  expect_equal(length(na_cols), 16)
  expect_equal(as.vector(na_cols), c(48L, 49L, 97L, 98L, 147L, 148L, 198L, 199L, 246L, 247L, 295L, 
                                     296L, 334L, 335L, 362L, 363L))
  
  # Check some measurements of heatmap matrix
  expect_equal(round(sum(htmp_matrix_wo_na), 2), -3631.17)
  expect_equal(round(mean(htmp_matrix_wo_na), 2), -0.03)
  expect_equal(round(median(htmp_matrix_wo_na), 2), -0.07)
  expect_equal(round(sd(htmp_matrix_wo_na), 2), 0.49)
  
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["0%"]), 2), -1)
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["25%"]), 2), -0.25)
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["50%"]), 2), -0.07)
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["75%"]), 2), -0.02)
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["100%"]), 2), 1)
  
  # Checking order of columns
  expect_equal(head(colnames(htmp_matrix),10), c("GATCTACTGGTGAG-1", "ACAGTGACTCACCC-1", "AGACGTACAGAGGC-1", 
                                                 "GACGTAACCTGTGA-1", "TATACAGATCCAGA-1", "CGGATAACAGCTCA-1", "TTACCATGGTTGAC-1", 
                                                 "TCCCACGATCATTC-1", "ATAGCGTGCAGATC-1", "TGTAGGTGTGCTGA-1"))
  
  # Checking order of columns
  expect_equal(tail(colnames(htmp_matrix),10), c("ACCTGAGATATCGG-1", "GACGCTCTCTCTCG-1", "GCGCATCTGGTTAC-1", 
                                                 "AGTCTTACTTCGGA-1", "GGCATATGGGGAGT-1", "TTACGTACGTTCAG-1", "GGAACACTTCAGAC-1", 
                                                 "ATCATCTGACACCA-1", "ACGAACTGGCTATG-1", "TAACACCTTGTTTC-1"))
  
  # Checking order of rows
  expect_equal(rownames(htmp_matrix), c("PKIG", "ISG20", "CD79B", "LINC00926", "MS4A1", " ", " ", "CAPZA2", 
                                        "MTURN", "ZHX1-C8ORF76", "MEST", "RP11-879F14.2", " ", " ", "CLIC1", 
                                        "ABI3", "IFITM2", "RHOC", "FCGR3A", "HES4", " ", " ", "CD27", 
                                        "AES", "IL7R", "LDHB", "CD3E", "IL32", "CD3D", " ", " ", "PTPN18", 
                                        "SMIM5", "PRUNE", "ANKRD9", "GFI1B", "LCN2", "ENKUR", " ", " ", 
                                        "AC092295.7", "TJP2", "PRKAR2B", "XPNPEP1", "TMEM91", "TUBA1C", 
                                        "SLC40A1", "SPHK1", "FAM63A", "HIST1H2BJ", "CLDN5", "C2orf88", 
                                        " ", " ", "RNASET2", "CD37", "SLC25A6", "HLA-DMB", "LY86", "FAM26F", 
                                        "HLA-DQB1", "HLA-DQA1", "CD74", "HLA-DPB1", "HLA-DRB5", "HLA-DRA", 
                                        "HLA-DRB1", "HLA-DPA1", " ", " ", "MOB3B", "MGLL", "NFE2", "FN3K", 
                                        "TRIM58", "TMCC2", "ZNF778", "SERPINE2", "FAM212B", "C19orf33", 
                                        "AC147651.3", "SSX2IP", "SMOX", "LGALSL", " ", " ", "GLUL", "H1F0", 
                                        "C1orf198", "HIST2H2BE", "SCGB1C1", "MLH3", "SEPT4", "SENCR", 
                                        "FNTB", "HIST1H2BD", "RP11-367G6.3", "GP1BA", "CMTM5", "CLEC1B", 
                                        " ", " ", "PPP1R14A", "WIPI1", "GLA", "ARG2", "ABHD16A", "NEXN", 
                                        "ENDOD1", "CTNNAL1", "TGFB1I1", "RP11-359I18.5", "ABCC3", "FAM212A", 
                                        "ATP9A", "SCFD2", "ALOX12", "ITGB3", "ARHGAP21", "SPOCD1", " ", 
                                        " ", "BTK", "KIFC3", "DAB2", "RILP", "FAM110A", "LYL1", "CLIC4", 
                                        "RAB27B", "DPY19L1", "PYGL", "PCP2", "KIAA0513", "SH3BGRL2", 
                                        "GATA2", "HEXIM2", "SEC14L5", "hsa-mir-1199", "AC137932.6", "EGFL7", 
                                        "HEMGN", " ", " ", "PLEK", "EFHD2", "XCL1", "CD7", "CD99", "KLRF1", 
                                        "IGFBP7", "AKR1C3", "CD247", "HOPX", "IL2RB", "CLIC3", "PTPRCAP", 
                                        "CCL4", "SPON2", "PRF1", "GNLY", "GZMB", "FGFBP2", "HLA-C", "CTSW", 
                                        "CST7", "GZMA", "CCL5", " ", " ", "POU2F2", "GCA", "CNPY3", "CYBB", 
                                        "CAMK1", "ATG3", "CD302", "NUP214", "GRINA", "EIF4EBP1", "ALDH2", 
                                        "GAPDH", "S100A10", "CD300LF", "FCGR2A", "NAAA", "CTSZ", "BID", 
                                        "CAPG", "BLVRA", "TGFBI", "SOD2", "IFNGR2", "SLC7A7", "C1orf162", 
                                        "LILRB4", "LILRB2", "IFI6", "WARS", "AP2S1", "MAFB", "CSTB", 
                                        "IGSF6", "CD14", "MNDA", "CPVL", "TNFSF13B", "APOBEC3A", "GABARAP", 
                                        "LGALS2", "S100A11", "CFD", "CTSS", "COTL1", "TYMP", " ", " ", 
                                        "PVALB", "CD151", "MYL12A", "PARVB", "TPTEP1", "H2AFJ", "MARCH2", 
                                        "GSN", "MFSD1", "RAP1B", "ASAH1", "YWHAH", "APP", "SNN", "HIST1H2BK", 
                                        "PLA2G12A", "FERMT3", "ACTN1", "PDLIM1", "ILK", "CTSA", "TPM4", 
                                        "SNCA", "ODC1", "GRAP2", "GAS2L1", "PTGS1", "NGFRAP1", "TREML1", 
                                        "MYL9", "TPM1", "F13A1", "CA2", "MPP1", "RGS18", "CLU", "AP001189.4", 
                                        "MMD", "GP9", "PPBP", "SPARC", "TUBB1", "TMEM40", "HIST1H2AC", 
                                        "PF4", "GNG11", "SDPR", "ACRBP", "PTCRA", " ", " ", "RPL6", "RPSAP58", 
                                        "RPLP1", "RPS9", "RPL26", "RPL27", "NACA", "RPL17", "RPL7A", 
                                        "EEF1B2", "RPL5", "RPL36", "RPS28", "RPS20", "RPL18", "RPL15", 
                                        "BTG1", "TPT1", "RPL35A", "RPS26", "RPS11", "RPS13", "RPS10", 
                                        "RPS29", "RPL28", "RPS8", "RPS4X", "RPL31", "EEF1A1", "RPL12", 
                                        "RPS16", "RPL27A", "RPS23", "RPS15A", "RPS3A", "RPL14", "RPL29", 
                                        "RPL3", "RPS18", "RPLP2", "RPL32", "RPL11", "RPL18A", "RPS12", 
                                        "RPS3", "RPS19", "RPS25", "RPL21", "RPL23A", "RPL9", "MALAT1"
  ))
  
  # Checking title of the heatmap
  htmp_title <- htmp@annotations@listData$col_title2@data
  expect_equal(htmp_title, "Expression heatmap")
  
  # Checking x lab of the heatmap
  htmp_xlab <- htmp@annotations@listData$col_title@data
  expect_equal(htmp_xlab, "Samples")
  
  # Checking y lab of the heatmap
  htmp_ylab <- htmp@annotations@listData$row_title@data
  expect_equal(htmp_ylab, "Genes")
  
  # Checking colorbar name of the heatmap
  htmp_colors_name <- htmp@colorbars@listData$`Exp. level`@title
  expect_equal(htmp_colors_name, "Exp. level")
  
  # Checking if colorbar is shown of the heatmap
  htmp_show_colorbar <- htmp@plots@listData$`Exp. level`@show_colorbar
  expect_true(htmp_show_colorbar)
  
  # Checking if dendrogram is showed
  htmp_show_dendro <- names(htmp@shapes@listData)
  expect_equal(htmp_show_dendro, NULL)
})

test_that("Checking plot_heatmap() using top genes and cell clusters", {
  
  set.seed(123)
  cell_clusters <- Idents(pbmc3k_medium)
  res <- top_genes(res, top = 20)
  
  htmp <- plot_heatmap(
    object = res,
    cell_clusters = cell_clusters,
    use_top_genes = TRUE,
    show_dendro = TRUE,
    line_size_vertical = 2,
    line_size_horizontal = 2,
    name = "Expression heatmap",
    show_legend = TRUE,
    xlab = "Samples",
    ylab = "Genes"
  )
  
  htmp_matrix <- htmp@plots@listData$`Exp. level`@data
  htmp_matrix_wo_na <- htmp_matrix[!is.na(htmp_matrix)]
  
  # Check heatmap format

  # Check dim of heatmap matrix
  expect_equal(ncol(htmp_matrix), 377)
  expect_equal(nrow(htmp_matrix), 230)
  expect_equal(dim(htmp_matrix), c(230, 377))
  
  # Check NA rows added by line_size_horizontal
  na_rows <- which(apply(htmp_matrix, 1, function(x) all(is.na(x))))
  expect_equal(length(na_rows), 28)
  expect_equal(as.vector(na_rows), c(6L, 7L, 13L, 14L, 21L, 22L, 30L, 31L, 39L, 40L, 53L, 54L, 69L, 
                                     70L, 85L, 86L, 101L, 102L, 121L, 122L, 143L, 144L, 165L, 166L, 
                                     187L, 188L, 209L, 210L))
  
  # Check NA rows added by line_size_vertical
  na_cols <- which(apply(htmp_matrix, 2, function(x) all(is.na(x))))
  expect_equal(length(na_cols), 16)
  expect_equal(as.vector(na_cols), c(48L, 49L, 97L, 98L, 147L, 148L, 198L, 199L, 246L, 247L, 295L, 
                                     296L, 334L, 335L, 362L, 363L))
  
  # Check some measurements of heatmap matrix
  expect_equal(round(sum(htmp_matrix_wo_na), 2), -3394.67)
  expect_equal(round(mean(htmp_matrix_wo_na), 2), -0.05)
  expect_equal(round(median(htmp_matrix_wo_na), 2), -0.05)
  expect_equal(round(sd(htmp_matrix_wo_na), 2), 0.48)
  
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["0%"]), 2), -1)
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["25%"]), 2), -0.22)
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["50%"]), 2), -0.05)
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["75%"]), 2), -0.02)
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["100%"]), 2), 1)
  
  # Checking order of columns
  expect_equal(head(colnames(htmp_matrix),10), c("GATCTACTGGTGAG-1", "ACAGTGACTCACCC-1", "AGACGTACAGAGGC-1", 
                                                 "GACGTAACCTGTGA-1", "TATACAGATCCAGA-1", "CGGATAACAGCTCA-1", "TTACCATGGTTGAC-1", 
                                                 "TCCCACGATCATTC-1", "ATAGCGTGCAGATC-1", "TGTAGGTGTGCTGA-1"))
  
  # Checking order of columns
  expect_equal(tail(colnames(htmp_matrix),10), c("ACCTGAGATATCGG-1", "GACGCTCTCTCTCG-1", "GCGCATCTGGTTAC-1", 
                                                 "AGTCTTACTTCGGA-1", "GGCATATGGGGAGT-1", "TTACGTACGTTCAG-1", "GGAACACTTCAGAC-1", 
                                                 "ATCATCTGACACCA-1", "ACGAACTGGCTATG-1", "TAACACCTTGTTTC-1"))
  
  # Checking order of rows
  expect_equal(rownames(htmp_matrix), c("PKIG", "ISG20", "MS4A1", "CD79B", "LINC00926", " ", " ", "CAPZA2", 
                                        "MEST", "MTURN", "ZHX1-C8ORF76", "RP11-879F14.2", " ", " ", "CLIC1", 
                                        "HES4", "ABI3", "IFITM2", "RHOC", "FCGR3A", " ", " ", "AES", 
                                        "CD27", "IL7R", "CD3E", "LDHB", "IL32", "CD3D", " ", " ", "PTPN18", 
                                        "SMIM5", "PRUNE", "ANKRD9", "GFI1B", "LCN2", "ENKUR", " ", " ", 
                                        "XPNPEP1", "AC092295.7", "TJP2", "PRKAR2B", "TMEM91", "SPHK1", 
                                        "SLC40A1", "TUBA1C", "FAM63A", "C2orf88", "CLDN5", "HIST1H2BJ", 
                                        " ", " ", "SLC25A6", "RNASET2", "CD37", "FAM26F", "HLA-DMB", 
                                        "LY86", "HLA-DQB1", "HLA-DQA1", "CD74", "HLA-DPB1", "HLA-DRB5", 
                                        "HLA-DPA1", "HLA-DRA", "HLA-DRB1", " ", " ", "FN3K", "MGLL", 
                                        "C19orf33", "NFE2", "MOB3B", "TRIM58", "SERPINE2", "ZNF778", 
                                        "TMCC2", "AC147651.3", "FAM212B", "SSX2IP", "SMOX", "LGALSL", 
                                        " ", " ", "H1F0", "GLUL", "HIST2H2BE", "MLH3", "C1orf198", "SENCR", 
                                        "SEPT4", "FNTB", "RP11-367G6.3", "HIST1H2BD", "SCGB1C1", "CMTM5", 
                                        "GP1BA", "CLEC1B", " ", " ", "PPP1R14A", "WIPI1", "GLA", "ABHD16A", 
                                        "TGFB1I1", "RP11-359I18.5", "NEXN", "ENDOD1", "ALOX12", "ITGB3", 
                                        "ARG2", "SCFD2", "ABCC3", "CTNNAL1", "FAM212A", "ATP9A", "SPOCD1", 
                                        "ARHGAP21", " ", " ", "LYL1", "BTK", "FAM110A", "PYGL", "CLIC4", 
                                        "RILP", "DAB2", "KIAA0513", "SH3BGRL2", "RAB27B", "KIFC3", "PCP2", 
                                        "HEXIM2", "DPY19L1", "GATA2", "SEC14L5", "EGFL7", "hsa-mir-1199", 
                                        "HEMGN", "AC137932.6", " ", " ", "CCL5", "IGFBP7", "XCL1", "KLRF1", 
                                        "AKR1C3", "CD7", "HLA-C", "HOPX", "IL2RB", "CD247", "CLIC3", 
                                        "CTSW", "SPON2", "CCL4", "CST7", "GZMA", "FGFBP2", "GNLY", "PRF1", 
                                        "GZMB", " ", " ", "BLVRA", "CAPG", "C1orf162", "TGFBI", "SLC7A7", 
                                        "IFI6", "MAFB", "CPVL", "CSTB", "IGSF6", "AP2S1", "APOBEC3A", 
                                        "TNFSF13B", "GABARAP", "LGALS2", "CFD", "S100A11", "CTSS", "COTL1", 
                                        "TYMP", " ", " ", "MPP1", "PTGS1", "MYL9", "MMD", "TMEM40", "CA2", 
                                        "F13A1", "TREML1", "PTCRA", "RGS18", "AP001189.4", "CLU", "GP9", 
                                        "HIST1H2AC", "TUBB1", "PPBP", "SPARC", "PF4", "GNG11", "SDPR", 
                                        " ", " ", "RPL29", "RPS4X", "RPS8", "EEF1A1", "RPS25", "RPS19", 
                                        "RPL28", "RPLP2", "RPS15A", "RPL12", "RPS3", "RPL23A", "RPL9", 
                                        "RPS16", "RPS23", "RPS18", "RPL18A", "RPS12", "RPL32", "RPL11"
  ))
  
  # Checking title of the heatmap
  htmp_title <- htmp@annotations@listData$col_title2@data
  expect_equal(htmp_title, "Expression heatmap")
  
  # Checking x lab of the heatmap
  htmp_xlab <- htmp@annotations@listData$col_title@data
  expect_equal(htmp_xlab, "Samples")
  
  # Checking y lab of the heatmap
  htmp_ylab <- htmp@annotations@listData$row_title@data
  expect_equal(htmp_ylab, "Genes")
  
  # Checking colorbar name of the heatmap
  htmp_colors_name <- htmp@colorbars@listData$`Exp. level`@title
  expect_equal(htmp_colors_name, "Exp. level")
  
  # Checking if colorbar is shown of the heatmap
  htmp_show_colorbar <- htmp@plots@listData$`Exp. level`@show_colorbar
  expect_true(htmp_show_colorbar)
  
  # Checking if dendrogram is showed
  htmp_show_dendro <- names(htmp@shapes@listData)
  expect_equal(htmp_show_dendro, NULL)
})


# ==============================================================================
# ==============================================================================


test_that("Checking plot_heatmap() displaying gene clusters 1 and 4", {
  set.seed(123)
  cell_clusters <- Idents(pbmc3k_medium)
  res <- top_genes(res, top = 20)
  
  htmp <- plot_heatmap(
    object = res[c(1,4), ],
    cell_clusters = cell_clusters,
    use_top_genes = TRUE,
    show_dendro = TRUE,
    line_size_vertical = 2,
    line_size_horizontal = 2,
    name = "Expression heatmap",
    show_legend = TRUE,
    xlab = "Samples",
    ylab = "Genes"
  )
  
  htmp_matrix <- htmp@plots@listData$`Exp. level`@data
  htmp_matrix_wo_na <- htmp_matrix[!is.na(htmp_matrix)]
  
  # Check heatmap format

  # Check dim of heatmap matrix
  expect_equal(ncol(htmp_matrix), 377)
  expect_equal(nrow(htmp_matrix), 42)
  expect_equal(dim(htmp_matrix), c(42, 377))
  
  # Check NA rows added by line_size_horizontal
  na_rows <- which(apply(htmp_matrix, 1, function(x) all(is.na(x))))
  expect_equal(length(na_rows), 2)
  expect_equal(as.vector(na_rows), c(21, 22))
  
  # Check NA rows added by line_size_vertical
  na_cols <- which(apply(htmp_matrix, 2, function(x) all(is.na(x))))
  expect_equal(length(na_cols), 16)
  expect_equal(as.vector(na_cols), c(48L, 49L, 97L, 98L, 147L, 148L, 198L, 199L, 246L, 247L, 295L, 
                                     296L, 334L, 335L, 362L, 363L))
  
  # Check some measurements of heatmap matrix
  expect_equal(round(sum(htmp_matrix_wo_na), 2), -509.92)
  expect_equal(round(mean(htmp_matrix_wo_na), 2), -0.04)
  expect_equal(round(median(htmp_matrix_wo_na), 2), -0.11)
  expect_equal(round(sd(htmp_matrix_wo_na), 2), 0.57)
  
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["0%"]), 2), -1)
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["25%"]), 2), -0.49)
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["50%"]), 2), -0.11)
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["75%"]), 2), 0.4)
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["100%"]), 2), 1)
  
  # Checking order of columns
  expect_equal(head(colnames(htmp_matrix),10), c("GATCTACTGGTGAG-1", "ACAGTGACTCACCC-1", "AGACGTACAGAGGC-1", 
                                                 "GACGTAACCTGTGA-1", "TATACAGATCCAGA-1", "CGGATAACAGCTCA-1", "TTACCATGGTTGAC-1", 
                                                 "TCCCACGATCATTC-1", "ATAGCGTGCAGATC-1", "TGTAGGTGTGCTGA-1"))
  
  # Checking order of columns
  expect_equal(tail(colnames(htmp_matrix),10), c("ACCTGAGATATCGG-1", "GACGCTCTCTCTCG-1", "GCGCATCTGGTTAC-1", 
                                                 "AGTCTTACTTCGGA-1", "GGCATATGGGGAGT-1", "TTACGTACGTTCAG-1", "GGAACACTTCAGAC-1", 
                                                 "ATCATCTGACACCA-1", "ACGAACTGGCTATG-1", "TAACACCTTGTTTC-1"))
  
  # Checking order of rows
  expect_equal(rownames(htmp_matrix), c("CCL5", "IGFBP7", "XCL1", "KLRF1", "AKR1C3", "CD7", "HLA-C", 
                                        "HOPX", "IL2RB", "CD247", "CLIC3", "CTSW", "SPON2", "CCL4", "CST7", 
                                        "GZMA", "FGFBP2", "GNLY", "PRF1", "GZMB", " ", " ", "RPL29", 
                                        "RPS4X", "RPS8", "EEF1A1", "RPS25", "RPS19", "RPL28", "RPLP2", 
                                        "RPS15A", "RPL12", "RPS3", "RPL23A", "RPL9", "RPS16", "RPS23", 
                                        "RPS18", "RPL18A", "RPS12", "RPL32", "RPL11"))
  
  # Checking title of the heatmap
  htmp_title <- htmp@annotations@listData$col_title2@data
  expect_equal(htmp_title, "Expression heatmap")
  
  # Checking x lab of the heatmap
  htmp_xlab <- htmp@annotations@listData$col_title@data
  expect_equal(htmp_xlab, "Samples")
  
  # Checking y lab of the heatmap
  htmp_ylab <- htmp@annotations@listData$row_title@data
  expect_equal(htmp_ylab, "Genes")
  
  # Checking colorbar name of the heatmap
  htmp_colors_name <- htmp@colorbars@listData$`Exp. level`@title
  expect_equal(htmp_colors_name, "Exp. level")
  
  # Checking if colorbar is shown of the heatmap
  htmp_show_colorbar <- htmp@plots@listData$`Exp. level`@show_colorbar
  expect_true(htmp_show_colorbar)
  
  # Checking if dendrogram is showed
  htmp_show_dendro <- names(htmp@shapes@listData)
  expect_equal(htmp_show_dendro, NULL)
})


# ==============================================================================
# ==============================================================================


test_that("Checking plot_heatmap() displaying only gene cluster 3", {
  
  set.seed(123)
  cell_clusters <- Idents(pbmc3k_medium)
  res <- top_genes(res, top = 20)
  
  htmp <- plot_heatmap(
    object = res[3, ],
    cell_clusters = cell_clusters,
    use_top_genes = TRUE,
    show_dendro = TRUE,
    line_size_vertical = 2,
    line_size_horizontal = 2,
    name = "Expression heatmap",
    show_legend = TRUE,
    xlab = "Samples",
    ylab = "Genes"
  )
  
  htmp_matrix <- htmp@plots@listData$`Exp. level`@data
  htmp_matrix_wo_na <- htmp_matrix[!is.na(htmp_matrix)]
  
  # Check heatmap format

  # Check dim of heatmap matrix
  expect_equal(ncol(htmp_matrix), 377)
  expect_equal(nrow(htmp_matrix), 20)
  expect_equal(dim(htmp_matrix), c(20, 377))
  
  # Check NA rows added by line_size_horizontal
  na_rows <- which(apply(htmp_matrix, 1, function(x) all(is.na(x))))
  expect_equal(length(na_rows), 0)
  expect_equal(as.vector(na_rows), integer(0))
  
  # Check NA rows added by line_size_vertical
  na_cols <- which(apply(htmp_matrix, 2, function(x) all(is.na(x))))
  expect_equal(length(na_cols), 16)
  expect_equal(as.vector(na_cols), c(48L, 49L, 97L, 98L, 147L, 148L, 198L, 199L, 246L, 247L, 295L, 
                                     296L, 334L, 335L, 362L, 363L))
  
  # Check some measurements of heatmap matrix
  expect_equal(round(sum(htmp_matrix_wo_na), 2), -663.16)
  expect_equal(round(mean(htmp_matrix_wo_na), 2), -0.09)
  expect_equal(round(median(htmp_matrix_wo_na), 2), -0.27)
  expect_equal(round(sd(htmp_matrix_wo_na), 2), 0.72)
  
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["0%"]), 2), -1)
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["25%"]), 2), -0.59)
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["50%"]), 2), -0.27)
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["75%"]), 2), 0.78)
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["100%"]), 2), 1)
  
  # Checking order of columns
  expect_equal(head(colnames(htmp_matrix),10), c("GATCTACTGGTGAG-1", "ACAGTGACTCACCC-1", "AGACGTACAGAGGC-1", 
                                                 "GACGTAACCTGTGA-1", "TATACAGATCCAGA-1", "CGGATAACAGCTCA-1", "TTACCATGGTTGAC-1", 
                                                 "TCCCACGATCATTC-1", "ATAGCGTGCAGATC-1", "TGTAGGTGTGCTGA-1"))
  
  # Checking order of columns
  expect_equal(tail(colnames(htmp_matrix),10), c("ACCTGAGATATCGG-1", "GACGCTCTCTCTCG-1", "GCGCATCTGGTTAC-1", 
                                                 "AGTCTTACTTCGGA-1", "GGCATATGGGGAGT-1", "TTACGTACGTTCAG-1", "GGAACACTTCAGAC-1", 
                                                 "ATCATCTGACACCA-1", "ACGAACTGGCTATG-1", "TAACACCTTGTTTC-1"))
  
  # Checking order of rows
  expect_equal(rownames(htmp_matrix), c("BLVRA", "CAPG", "C1orf162", "TGFBI", "SLC7A7", "IFI6", "MAFB", 
                                        "CPVL", "CSTB", "IGSF6", "AP2S1", "APOBEC3A", "TNFSF13B", "GABARAP", 
                                        "LGALS2", "CFD", "S100A11", "CTSS", "COTL1", "TYMP"))
  
  # Checking title of the heatmap
  htmp_title <- htmp@annotations@listData$col_title2@data
  expect_equal(htmp_title, "Expression heatmap")
  
  # Checking x lab of the heatmap
  htmp_xlab <- htmp@annotations@listData$col_title@data
  expect_equal(htmp_xlab, "Samples")
  
  # Checking y lab of the heatmap
  htmp_ylab <- htmp@annotations@listData$row_title@data
  expect_equal(htmp_ylab, "Genes")
  
  # Checking colorbar name of the heatmap
  htmp_colors_name <- htmp@colorbars@listData$`Exp. level`@title
  expect_equal(htmp_colors_name, "Exp. level")
  
  # Checking if colorbar is shown of the heatmap
  htmp_show_colorbar <- htmp@plots@listData$`Exp. level`@show_colorbar
  expect_true(htmp_show_colorbar)
  
  # Checking if dendrogram is showed
  htmp_show_dendro <- names(htmp@shapes@listData)
  expect_equal(htmp_show_dendro, NULL)
})




testthat::test_that("Check plot_heatmap, subsetting and clustering", {
  # Set verbosity to 0
  set_verbosity(0)
  
  # Load seurat object
  data(pbmc_small, package = "SeuratObject")
  
  ## Select informative genes
  clust_set <- select_genes(data=pbmc_small,
                            distance_method="pearson",
                            k=10,
                            row_sum=-Inf,
                            noise_level=0.95,
                            fdr = 1e-6)
  
  ## Cluster genes
  clust_set <- gene_clustering(object = clust_set,
                               inflation = 1.2,
                               keep_nn = FALSE,
                               k = 5,
                               threads = 1)
  
  ident <- Seurat::Idents(pbmc_small)
  
  # Check color arg
  expect_error(plot_heatmap(clust_set, cell_clusters=Seurat::Idents(pbmc_small), colors_cell_clusters = "green"))
  p <- plot_heatmap(clust_set[, names(ident[ident==2])], cell_clusters=Seurat::Idents(pbmc_small), colors_cell_clusters = "red")
  testthat::expect_error(print(p), NA)
  expect_equal(unique(p@colorbars@listData$Ident.@colors), "red")
  
  # Check if seurat identities is considered
  
  p <- plot_heatmap(clust_set[,colnames(clust_set@data)[1:80]], cell_clusters=Seurat::Idents(pbmc_small))
  testthat::expect_error(print(p), NA)
  expect_equal(dim(p@plots@listData$`Exp. level`@data), c(178, 86)) # 80 + 6 NA columns
  
  p <- plot_heatmap(clust_set[,colnames(clust_set@data)[1:10]], cell_clusters=Seurat::Idents(pbmc_small))
  testthat::expect_error(print(p), NA)
  expect_equal(dim(p@plots@listData$`Exp. level`@data), c(178, 10))
  
  p <- plot_heatmap(clust_set[,colnames(clust_set@data)[1:20]], cell_clusters=Seurat::Idents(pbmc_small))
  testthat::expect_error(print(p), NA)
  expect_equal(dim(p@plots@listData$`Exp. level`@data), c(178, 23)) # 20 + 3 NA columns
  
  p <- plot_heatmap(clust_set[1,colnames(clust_set@data)[1:10]], cell_clusters=Seurat::Idents(pbmc_small))
  testthat::expect_error(print(p), NA)
  expect_equal(dim(p@plots@listData$`Exp. level`@data), c(91, 10))
  expect_equal(dim(p@plots@listData$`Exp. level`@data), c(length(clust_set@gene_clusters$`1`), 10))
  
  p <- plot_heatmap(clust_set[2 ,colnames(clust_set@data)[1:10]], cell_clusters=Seurat::Idents(pbmc_small))
  testthat::expect_error(print(p), NA)
  expect_equal(dim(p@plots@listData$`Exp. level`@data), c(52, 10))
  expect_equal(dim(p@plots@listData$`Exp. level`@data), c(length(clust_set@gene_clusters$`2`), 10))
})




testthat::test_that("Test non interactive heatmaps.", {
  # Set verbosity to 0
  set_verbosity(0)
  
  # Load seurat object
  data(pbmc_small, package = "SeuratObject")
  
  ## Select informative genes
  clust_set <- select_genes(data=pbmc_small,
                            distance_method="pearson",
                            k=10,
                            row_sum=-Inf,
                            noise_level=0.95,
                            fdr = 1e-6)
  
  ## Cluster genes
  clust_set <- gene_clustering(object = clust_set,
                               inflation = 1.2,
                               keep_nn = FALSE,
                               k = 5,
                               threads = 1)
  
  ident <- Seurat::Idents(pbmc_small)
  
  n_gene_clusters <- length(clust_set@gene_clusters) 
  n_cell_clusters <- 3
  p_int <- plot_heatmap(clust_set, cell_clusters=Seurat::Idents(pbmc_small), 
                        colors_cell_clusters =  c("green", "red", "yellow"))
  p_noi <- plot_heatmap(clust_set, cell_clusters=Seurat::Idents(pbmc_small), 
                        colors_cell_clusters =  c("green", "red", "yellow"),
                        interactive = FALSE)
  p_int_mat <- htmp_matrix <- p_int@plots@listData$`Exp. level`@data
  p_noi_rown <- p_noi[["gtable"]][["grobs"]][[2]][["label"]]
  testthat::expect_equal(p_noi_rown, rev(rownames(p_int_mat)))
  
  p_noi_colnames <- names(p_noi[["gtable"]][["grobs"]][[3]][["gp"]][["fill"]][,1])
  expect_true(all(colnames(p_int_mat) == p_noi_colnames))
  expect_true(length(grep(" ", p_noi_colnames)) == (n_cell_clusters-1) * 3)
  expect_true(length(grep(" ", p_noi_rown)) == (n_gene_clusters-1) * 3)
  
  n_gene_clusters <- 1 
  n_cell_clusters <- 3
  p_noi <- plot_heatmap(clust_set[1,], cell_clusters=Seurat::Idents(pbmc_small), 
                        colors_cell_clusters =  c("green", "red", "yellow"),
                        interactive = FALSE)
  p_noi_rown <- p_noi[["gtable"]][["grobs"]][[2]][["label"]]
  expect_true(length(grep(" ", p_noi_colnames)) == (n_cell_clusters-1) * 3)
  expect_true(length(grep(" ", p_noi_rown)) == (n_gene_clusters-1) * 3)
  
  n_gene_clusters <- 2 
  n_cell_clusters <- 3
  p_noi <- plot_heatmap(clust_set[c(1, 3), colnames(pbmc_small)], cell_clusters=Seurat::Idents(pbmc_small), 
                        colors_cell_clusters =  c("green", "red", "yellow"),
                        interactive = FALSE)
  p_noi_rown <- p_noi[["gtable"]][["grobs"]][[2]][["label"]]
  expect_true(length(grep(" ", p_noi_colnames)) == (n_cell_clusters-1) * 3)
  expect_true(length(grep(" ", p_noi_rown)) == (n_gene_clusters-1) * 3)
})



# Create matrix containing 4 signatures
m <- create_4_rnd_clust()

set_verbosity(0)

## A rather stringent version
res <- select_genes(data = m,
                    distance_method = "pearson",
                    k = 75,
                    row_sum = -Inf,
                    noise_level = 0.3,
                    fdr = 1e-8)

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
# Set verbosity to 0
set_verbosity(0)
load_example_dataset("7871581/files/pbmc3k_medium")
load_example_dataset("7871581/files/pbmc3k_medium_clusters")

res <- pbmc3k_medium_clusters

test_that("Checking plot_heatmap()", {
  set.seed(123)
  
  htmp <- plot_heatmap(
    object = res,
    cell_clusters = NULL,
    use_top_genes = FALSE,
    show_dendro = TRUE,
    line_size_vertical = 2,
    line_size_horizontal = 2,
    name = "Expression heatmap",
    show_legend = TRUE,
    xlab = "Samples",
    ylab = "Genes"
  )
  
  htmp_matrix <- htmp@plots@listData$`Exp. level`@data
  htmp_matrix_wo_na <- htmp_matrix[!is.na(htmp_matrix)]
  
  # Check heatmap format

  # Check dim of heatmap matrix
  expect_equal(ncol(htmp_matrix), 361)
  expect_equal(nrow(htmp_matrix), 319)
  expect_equal(dim(htmp_matrix), c(319, 361))
  
  # Check NA rows added by line_size_horizontal
  na_rows <- which(apply(htmp_matrix, 1, function(x) all(is.na(x))))
  expect_equal(length(na_rows), 28)
  expect_equal(as.vector(na_rows), c(6L, 7L, 13L, 14L, 21L, 22L, 30L, 31L, 39L, 40L, 53L, 54L, 69L, 
                                     70L, 85L, 86L, 101L, 102L, 121L, 122L, 143L, 144L, 169L, 170L, 
                                     216L, 217L, 267L, 268L))
  
  # Check NA rows added by line_size_vertical
  na_cols <- which(apply(htmp_matrix, 2, function(x) all(is.na(x))))
  expect_equal(length(na_cols), 0)
  expect_equal(as.vector(na_cols), as.integer())
  
  # Check some measurements of heatmap matrix
  expect_equal(round(sum(htmp_matrix_wo_na), 2), -3631.17)
  expect_equal(round(mean(htmp_matrix_wo_na), 2), -0.03)
  expect_equal(round(median(htmp_matrix_wo_na), 2), -0.07)
  expect_equal(round(sd(htmp_matrix_wo_na), 2), 0.49)
  
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["0%"]), 2), -1)
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["25%"]), 2), -0.25)
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["50%"]), 2), -0.07)
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["75%"]), 2), -0.02)
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["100%"]), 2), 1)
  
  # Checking order of columns
  expect_equal(head(colnames(htmp_matrix),10), c("TAACTCACGTATCG-1", "CTAATGCTTGTGGT-1", "TTCGAGGATAGAAG-1", 
                                                 "TTCTTACTCTGGAT-1", "ATAAGTTGGTACGT-1", "TTACCATGAATCGC-1", "GGCGCATGCCTAAG-1", 
                                                 "ATACGGACCTACTT-1", "GGAACACTCACTTT-1", "AAGGTCACGGTTAC-1"))
  
  # Checking order of columns
  expect_equal(tail(colnames(htmp_matrix),10), c("TTCGTATGAAAAGC-1", "CCAAGAACTACTGG-1", "TCGAATCTCTGGTA-1", 
                                                 "GTCATACTGCGATT-1", "GGCAATACGGCATT-1", "GCTATACTAGCGTT-1", "CAGCTCTGCAAGCT-1", 
                                                 "TCAATCACACTCTT-1", "ATCCTAACGCTACA-1", "ACGATTCTACGGGA-1"))
  
  # Checking order of rows
  expect_equal(rownames(htmp_matrix), c("PKIG", "ISG20", "CD79B", "LINC00926", "MS4A1", " ", " ", "CAPZA2", 
                                        "MTURN", "ZHX1-C8ORF76", "MEST", "RP11-879F14.2", " ", " ", "CLIC1", 
                                        "ABI3", "IFITM2", "RHOC", "FCGR3A", "HES4", " ", " ", "CD27", 
                                        "AES", "IL7R", "LDHB", "CD3E", "IL32", "CD3D", " ", " ", "PTPN18", 
                                        "SMIM5", "PRUNE", "ANKRD9", "GFI1B", "LCN2", "ENKUR", " ", " ", 
                                        "AC092295.7", "TJP2", "PRKAR2B", "XPNPEP1", "TMEM91", "TUBA1C", 
                                        "SLC40A1", "SPHK1", "FAM63A", "HIST1H2BJ", "CLDN5", "C2orf88", 
                                        " ", " ", "RNASET2", "CD37", "SLC25A6", "HLA-DMB", "LY86", "FAM26F", 
                                        "HLA-DQB1", "HLA-DQA1", "CD74", "HLA-DPB1", "HLA-DRB5", "HLA-DRA", 
                                        "HLA-DRB1", "HLA-DPA1", " ", " ", "MOB3B", "MGLL", "NFE2", "FN3K", 
                                        "TRIM58", "TMCC2", "ZNF778", "SERPINE2", "FAM212B", "C19orf33", 
                                        "AC147651.3", "SSX2IP", "SMOX", "LGALSL", " ", " ", "GLUL", "H1F0", 
                                        "C1orf198", "HIST2H2BE", "SCGB1C1", "MLH3", "SEPT4", "SENCR", 
                                        "FNTB", "HIST1H2BD", "RP11-367G6.3", "GP1BA", "CMTM5", "CLEC1B", 
                                        " ", " ", "PPP1R14A", "WIPI1", "GLA", "ARG2", "ABHD16A", "NEXN", 
                                        "ENDOD1", "CTNNAL1", "TGFB1I1", "RP11-359I18.5", "ABCC3", "FAM212A", 
                                        "ATP9A", "SCFD2", "ALOX12", "ITGB3", "ARHGAP21", "SPOCD1", " ", 
                                        " ", "BTK", "KIFC3", "DAB2", "RILP", "FAM110A", "LYL1", "CLIC4", 
                                        "RAB27B", "DPY19L1", "PYGL", "PCP2", "KIAA0513", "SH3BGRL2", 
                                        "GATA2", "HEXIM2", "SEC14L5", "hsa-mir-1199", "AC137932.6", "EGFL7", 
                                        "HEMGN", " ", " ", "PLEK", "EFHD2", "XCL1", "CD7", "CD99", "KLRF1", 
                                        "IGFBP7", "AKR1C3", "CD247", "HOPX", "IL2RB", "CLIC3", "PTPRCAP", 
                                        "CCL4", "SPON2", "PRF1", "GNLY", "GZMB", "FGFBP2", "HLA-C", "CTSW", 
                                        "CST7", "GZMA", "CCL5", " ", " ", "POU2F2", "GCA", "CNPY3", "CYBB", 
                                        "CAMK1", "ATG3", "CD302", "NUP214", "GRINA", "EIF4EBP1", "ALDH2", 
                                        "GAPDH", "S100A10", "CD300LF", "FCGR2A", "NAAA", "CTSZ", "BID", 
                                        "CAPG", "BLVRA", "TGFBI", "SOD2", "IFNGR2", "SLC7A7", "C1orf162", 
                                        "LILRB4", "LILRB2", "IFI6", "WARS", "AP2S1", "MAFB", "CSTB", 
                                        "IGSF6", "CD14", "MNDA", "CPVL", "TNFSF13B", "APOBEC3A", "GABARAP", 
                                        "LGALS2", "S100A11", "CFD", "CTSS", "COTL1", "TYMP", " ", " ", 
                                        "PVALB", "CD151", "MYL12A", "PARVB", "TPTEP1", "H2AFJ", "MARCH2", 
                                        "GSN", "MFSD1", "RAP1B", "ASAH1", "YWHAH", "APP", "SNN", "HIST1H2BK", 
                                        "PLA2G12A", "FERMT3", "ACTN1", "PDLIM1", "ILK", "CTSA", "TPM4", 
                                        "SNCA", "ODC1", "GRAP2", "GAS2L1", "PTGS1", "NGFRAP1", "TREML1", 
                                        "MYL9", "TPM1", "F13A1", "CA2", "MPP1", "RGS18", "CLU", "AP001189.4", 
                                        "MMD", "GP9", "PPBP", "SPARC", "TUBB1", "TMEM40", "HIST1H2AC", 
                                        "PF4", "GNG11", "SDPR", "ACRBP", "PTCRA", " ", " ", "RPL6", "RPSAP58", 
                                        "RPLP1", "RPS9", "RPL26", "RPL27", "NACA", "RPL17", "RPL7A", 
                                        "EEF1B2", "RPL5", "RPL36", "RPS28", "RPS20", "RPL18", "RPL15", 
                                        "BTG1", "TPT1", "RPL35A", "RPS26", "RPS11", "RPS13", "RPS10", 
                                        "RPS29", "RPL28", "RPS8", "RPS4X", "RPL31", "EEF1A1", "RPL12", 
                                        "RPS16", "RPL27A", "RPS23", "RPS15A", "RPS3A", "RPL14", "RPL29", 
                                        "RPL3", "RPS18", "RPLP2", "RPL32", "RPL11", "RPL18A", "RPS12", 
                                        "RPS3", "RPS19", "RPS25", "RPL21", "RPL23A", "RPL9", "MALAT1"
  ))
  
  # Checking title of the heatmap
  htmp_title <- htmp@annotations@listData$col_title2@data
  expect_equal(htmp_title, "Expression heatmap")
  
  # Checking x lab of the heatmap
  htmp_xlab <- htmp@annotations@listData$col_title@data
  expect_equal(htmp_xlab, "Samples")
  
  # Checking y lab of the heatmap
  htmp_ylab <- htmp@annotations@listData$row_title@data
  expect_equal(htmp_ylab, "Genes")
  
  # Checking colorbar name of the heatmap
  htmp_colors_name <- htmp@colorbars@listData$`Exp. level`@title
  expect_equal(htmp_colors_name, "Exp. level")
  
  # Checking if colorbar is shown of the heatmap
  htmp_show_colorbar <- htmp@plots@listData$`Exp. level`@show_colorbar
  expect_true(htmp_show_colorbar)
  
  # Checking if dendrogram is showed
  htmp_show_dendro <- names(htmp@shapes@listData)
  expect_equal(htmp_show_dendro, "col_dendro")
})

test_that("Checking plot_heatmap() using top genes", {
  set.seed(123)
  
  res_scigenex_top <- top_genes(res, top = 20)
  htmp <- plot_heatmap(
    object = res_scigenex_top,
    cell_clusters = NULL,
    use_top_genes = TRUE, # Use top genes
    show_dendro = FALSE, # Set to FALSE
    line_size_vertical = 2,
    line_size_horizontal = 2,
    name = NULL, # Set to NULL
    show_legend = FALSE, # Set to FALSE
    xlab = NULL, # Set to NULL
    ylab = NULL
  ) # Set to NULL
  
  htmp_matrix <- htmp@plots@listData$`Exp. level`@data
  htmp_matrix_wo_na <- htmp_matrix[!is.na(htmp_matrix)]
  
  # Check heatmap format

  # Check dim of heatmap matrix
  expect_equal(ncol(htmp_matrix), 361)
  expect_equal(nrow(htmp_matrix), 230)
  expect_equal(dim(htmp_matrix), c(230, 361))
  
  # Check NA rows added by line_size_horizontal
  na_rows <- which(apply(htmp_matrix, 1, function(x) all(is.na(x))))
  expect_equal(length(na_rows), 28)
  expect_equal(as.vector(na_rows), c(6L, 7L, 13L, 14L, 21L, 22L, 30L, 31L, 39L, 40L, 53L, 54L, 69L, 
                                     70L, 85L, 86L, 101L, 102L, 121L, 122L, 143L, 144L, 165L, 166L, 
                                     187L, 188L, 209L, 210L))
  
  # Check NA rows added by line_size_horizontal
  na_cols <- which(apply(htmp_matrix, 2, function(x) all(is.na(x))))
  expect_equal(length(na_cols), 0)
  expect_equal(as.vector(na_cols), as.integer())
  
  # Check some measurements of heatmap matrix
  expect_equal(round(sum(htmp_matrix_wo_na), 2), -3394.67)
  expect_equal(round(mean(htmp_matrix_wo_na), 2), -0.05)
  expect_equal(round(median(htmp_matrix_wo_na), 2), -0.05)
  expect_equal(round(sd(htmp_matrix_wo_na), 2), 0.48)
  
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["0%"]), 2), -1)
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["25%"]), 2), -0.22)
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["50%"]), 2), -0.05)
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["75%"]), 2), -0.02)
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["100%"]), 2), 1)
  
  # Checking order of columns
  expect_equal(head(colnames(htmp_matrix),10), c("TAACTCACGTATCG-1", "CTAATGCTTGTGGT-1", "TTCGAGGATAGAAG-1", 
                                                 "TTCTTACTCTGGAT-1", "ATAAGTTGGTACGT-1", "TTACCATGAATCGC-1", "GGCGCATGCCTAAG-1", 
                                                 "ATACGGACCTACTT-1", "GGAACACTCACTTT-1", "AAGGTCACGGTTAC-1"))
  
  expect_equal(tail(colnames(htmp_matrix),10), c("TTCGTATGAAAAGC-1", "CCAAGAACTACTGG-1", "TCGAATCTCTGGTA-1", 
                                                 "GTCATACTGCGATT-1", "GGCAATACGGCATT-1", "GCTATACTAGCGTT-1", "CAGCTCTGCAAGCT-1", 
                                                 "TCAATCACACTCTT-1", "ATCCTAACGCTACA-1", "ACGATTCTACGGGA-1"))  
  # Checking order of rows
  expect_equal(rownames(htmp_matrix), c("PKIG", "ISG20", "MS4A1", "CD79B", "LINC00926", " ", " ", "CAPZA2", 
                                        "MEST", "MTURN", "ZHX1-C8ORF76", "RP11-879F14.2", " ", " ", "CLIC1", 
                                        "HES4", "ABI3", "IFITM2", "RHOC", "FCGR3A", " ", " ", "AES", 
                                        "CD27", "IL7R", "CD3E", "LDHB", "IL32", "CD3D", " ", " ", "PTPN18", 
                                        "SMIM5", "PRUNE", "ANKRD9", "GFI1B", "LCN2", "ENKUR", " ", " ", 
                                        "XPNPEP1", "AC092295.7", "TJP2", "PRKAR2B", "TMEM91", "SPHK1", 
                                        "SLC40A1", "TUBA1C", "FAM63A", "C2orf88", "CLDN5", "HIST1H2BJ", 
                                        " ", " ", "SLC25A6", "RNASET2", "CD37", "FAM26F", "HLA-DMB", 
                                        "LY86", "HLA-DQB1", "HLA-DQA1", "CD74", "HLA-DPB1", "HLA-DRB5", 
                                        "HLA-DPA1", "HLA-DRA", "HLA-DRB1", " ", " ", "FN3K", "MGLL", 
                                        "C19orf33", "NFE2", "MOB3B", "TRIM58", "SERPINE2", "ZNF778", 
                                        "TMCC2", "AC147651.3", "FAM212B", "SSX2IP", "SMOX", "LGALSL", 
                                        " ", " ", "H1F0", "GLUL", "HIST2H2BE", "MLH3", "C1orf198", "SENCR", 
                                        "SEPT4", "FNTB", "RP11-367G6.3", "HIST1H2BD", "SCGB1C1", "CMTM5", 
                                        "GP1BA", "CLEC1B", " ", " ", "PPP1R14A", "WIPI1", "GLA", "ABHD16A", 
                                        "TGFB1I1", "RP11-359I18.5", "NEXN", "ENDOD1", "ALOX12", "ITGB3", 
                                        "ARG2", "SCFD2", "ABCC3", "CTNNAL1", "FAM212A", "ATP9A", "SPOCD1", 
                                        "ARHGAP21", " ", " ", "LYL1", "BTK", "FAM110A", "PYGL", "CLIC4", 
                                        "RILP", "DAB2", "KIAA0513", "SH3BGRL2", "RAB27B", "KIFC3", "PCP2", 
                                        "HEXIM2", "DPY19L1", "GATA2", "SEC14L5", "EGFL7", "hsa-mir-1199", 
                                        "HEMGN", "AC137932.6", " ", " ", "CCL5", "IGFBP7", "XCL1", "KLRF1", 
                                        "AKR1C3", "CD7", "HLA-C", "HOPX", "IL2RB", "CD247", "CLIC3", 
                                        "CTSW", "SPON2", "CCL4", "CST7", "GZMA", "FGFBP2", "GNLY", "PRF1", 
                                        "GZMB", " ", " ", "BLVRA", "CAPG", "C1orf162", "TGFBI", "SLC7A7", 
                                        "IFI6", "MAFB", "CPVL", "CSTB", "IGSF6", "AP2S1", "APOBEC3A", 
                                        "TNFSF13B", "GABARAP", "LGALS2", "CFD", "S100A11", "CTSS", "COTL1", 
                                        "TYMP", " ", " ", "MPP1", "PTGS1", "MYL9", "MMD", "TMEM40", "CA2", 
                                        "F13A1", "TREML1", "PTCRA", "RGS18", "AP001189.4", "CLU", "GP9", 
                                        "HIST1H2AC", "TUBB1", "PPBP", "SPARC", "PF4", "GNG11", "SDPR", 
                                        " ", " ", "RPL29", "RPS4X", "RPS8", "EEF1A1", "RPS25", "RPS19", 
                                        "RPL28", "RPLP2", "RPS15A", "RPL12", "RPS3", "RPL23A", "RPL9", 
                                        "RPS16", "RPS23", "RPS18", "RPL18A", "RPS12", "RPL32", "RPL11"
  ))
  
  # Checking title of the heatmap
  htmp_slots <- names(htmp@annotations@listData)
  expect_false("col_title2" %in% htmp_slots)
  
  # Checking x lab of the heatmap
  expect_false("col_title" %in% htmp_slots)
  
  # Checking y lab of the heatmap
  expect_false("row_title" %in% htmp_slots)
  
  # Checking colorbar name of the heatmap
  htmp_colors_name <- htmp@colorbars@listData$`Exp. level`@title
  expect_equal(htmp_colors_name, "Exp. level")
  
  # Checking if colorbar is showed
  htmp_show_colorbar <- htmp@plots@listData$`Exp. level`@show_colorbar
  expect_false(htmp_show_colorbar)
  
  # Checking if dendrogram is showed
  htmp_show_dendro <- names(htmp@shapes@listData)
  expect_equal(htmp_show_dendro, NULL)
})

test_that("Checking plot_heatmap() using cell clusters", {
  
  set.seed(123)
  cell_clusters <- Idents(pbmc3k_medium)
  htmp <- plot_heatmap(
    object = res,
    cell_clusters = cell_clusters,
    use_top_genes = FALSE,
    show_dendro = TRUE,
    line_size_vertical = 2,
    line_size_horizontal = 2,
    name = "Expression heatmap",
    show_legend = TRUE,
    xlab = "Samples",
    ylab = "Genes"
  )
  
  htmp_matrix <- htmp@plots@listData$`Exp. level`@data
  htmp_matrix_wo_na <- htmp_matrix[!is.na(htmp_matrix)]
  
  # Check heatmap format

  # Check dim of heatmap matrix
  expect_equal(ncol(htmp_matrix), 377)
  expect_equal(nrow(htmp_matrix), 319)
  expect_equal(dim(htmp_matrix), c(319, 377))
  
  # Check NA rows added by line_size_horizontal
  na_rows <- which(apply(htmp_matrix, 1, function(x) all(is.na(x))))
  expect_equal(length(na_rows), 28)
  expect_equal(as.vector(na_rows), c(6L, 7L, 13L, 14L, 21L, 22L, 30L, 31L, 39L, 40L, 53L, 54L, 69L, 
                                     70L, 85L, 86L, 101L, 102L, 121L, 122L, 143L, 144L, 169L, 170L, 
                                     216L, 217L, 267L, 268L))
  
  # Check NA rows added by line_size_vertical
  na_cols <- which(apply(htmp_matrix, 2, function(x) all(is.na(x))))
  expect_equal(length(na_cols), 16)
  expect_equal(as.vector(na_cols), c(48L, 49L, 97L, 98L, 147L, 148L, 198L, 199L, 246L, 247L, 295L, 
                                     296L, 334L, 335L, 362L, 363L))
  
  # Check some measurements of heatmap matrix
  expect_equal(round(sum(htmp_matrix_wo_na), 2), -3631.17)
  expect_equal(round(mean(htmp_matrix_wo_na), 2), -0.03)
  expect_equal(round(median(htmp_matrix_wo_na), 2), -0.07)
  expect_equal(round(sd(htmp_matrix_wo_na), 2), 0.49)
  
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["0%"]), 2), -1)
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["25%"]), 2), -0.25)
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["50%"]), 2), -0.07)
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["75%"]), 2), -0.02)
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["100%"]), 2), 1)
  
  # Checking order of columns
  expect_equal(head(colnames(htmp_matrix),10), c("GATCTACTGGTGAG-1", "ACAGTGACTCACCC-1", "AGACGTACAGAGGC-1", 
                                                 "GACGTAACCTGTGA-1", "TATACAGATCCAGA-1", "CGGATAACAGCTCA-1", "TTACCATGGTTGAC-1", 
                                                 "TCCCACGATCATTC-1", "ATAGCGTGCAGATC-1", "TGTAGGTGTGCTGA-1"))
  
  # Checking order of columns
  expect_equal(tail(colnames(htmp_matrix),10), c("ACCTGAGATATCGG-1", "GACGCTCTCTCTCG-1", "GCGCATCTGGTTAC-1", 
                                                 "AGTCTTACTTCGGA-1", "GGCATATGGGGAGT-1", "TTACGTACGTTCAG-1", "GGAACACTTCAGAC-1", 
                                                 "ATCATCTGACACCA-1", "ACGAACTGGCTATG-1", "TAACACCTTGTTTC-1"))
  
  # Checking order of rows
  expect_equal(rownames(htmp_matrix), c("PKIG", "ISG20", "CD79B", "LINC00926", "MS4A1", " ", " ", "CAPZA2", 
                                        "MTURN", "ZHX1-C8ORF76", "MEST", "RP11-879F14.2", " ", " ", "CLIC1", 
                                        "ABI3", "IFITM2", "RHOC", "FCGR3A", "HES4", " ", " ", "CD27", 
                                        "AES", "IL7R", "LDHB", "CD3E", "IL32", "CD3D", " ", " ", "PTPN18", 
                                        "SMIM5", "PRUNE", "ANKRD9", "GFI1B", "LCN2", "ENKUR", " ", " ", 
                                        "AC092295.7", "TJP2", "PRKAR2B", "XPNPEP1", "TMEM91", "TUBA1C", 
                                        "SLC40A1", "SPHK1", "FAM63A", "HIST1H2BJ", "CLDN5", "C2orf88", 
                                        " ", " ", "RNASET2", "CD37", "SLC25A6", "HLA-DMB", "LY86", "FAM26F", 
                                        "HLA-DQB1", "HLA-DQA1", "CD74", "HLA-DPB1", "HLA-DRB5", "HLA-DRA", 
                                        "HLA-DRB1", "HLA-DPA1", " ", " ", "MOB3B", "MGLL", "NFE2", "FN3K", 
                                        "TRIM58", "TMCC2", "ZNF778", "SERPINE2", "FAM212B", "C19orf33", 
                                        "AC147651.3", "SSX2IP", "SMOX", "LGALSL", " ", " ", "GLUL", "H1F0", 
                                        "C1orf198", "HIST2H2BE", "SCGB1C1", "MLH3", "SEPT4", "SENCR", 
                                        "FNTB", "HIST1H2BD", "RP11-367G6.3", "GP1BA", "CMTM5", "CLEC1B", 
                                        " ", " ", "PPP1R14A", "WIPI1", "GLA", "ARG2", "ABHD16A", "NEXN", 
                                        "ENDOD1", "CTNNAL1", "TGFB1I1", "RP11-359I18.5", "ABCC3", "FAM212A", 
                                        "ATP9A", "SCFD2", "ALOX12", "ITGB3", "ARHGAP21", "SPOCD1", " ", 
                                        " ", "BTK", "KIFC3", "DAB2", "RILP", "FAM110A", "LYL1", "CLIC4", 
                                        "RAB27B", "DPY19L1", "PYGL", "PCP2", "KIAA0513", "SH3BGRL2", 
                                        "GATA2", "HEXIM2", "SEC14L5", "hsa-mir-1199", "AC137932.6", "EGFL7", 
                                        "HEMGN", " ", " ", "PLEK", "EFHD2", "XCL1", "CD7", "CD99", "KLRF1", 
                                        "IGFBP7", "AKR1C3", "CD247", "HOPX", "IL2RB", "CLIC3", "PTPRCAP", 
                                        "CCL4", "SPON2", "PRF1", "GNLY", "GZMB", "FGFBP2", "HLA-C", "CTSW", 
                                        "CST7", "GZMA", "CCL5", " ", " ", "POU2F2", "GCA", "CNPY3", "CYBB", 
                                        "CAMK1", "ATG3", "CD302", "NUP214", "GRINA", "EIF4EBP1", "ALDH2", 
                                        "GAPDH", "S100A10", "CD300LF", "FCGR2A", "NAAA", "CTSZ", "BID", 
                                        "CAPG", "BLVRA", "TGFBI", "SOD2", "IFNGR2", "SLC7A7", "C1orf162", 
                                        "LILRB4", "LILRB2", "IFI6", "WARS", "AP2S1", "MAFB", "CSTB", 
                                        "IGSF6", "CD14", "MNDA", "CPVL", "TNFSF13B", "APOBEC3A", "GABARAP", 
                                        "LGALS2", "S100A11", "CFD", "CTSS", "COTL1", "TYMP", " ", " ", 
                                        "PVALB", "CD151", "MYL12A", "PARVB", "TPTEP1", "H2AFJ", "MARCH2", 
                                        "GSN", "MFSD1", "RAP1B", "ASAH1", "YWHAH", "APP", "SNN", "HIST1H2BK", 
                                        "PLA2G12A", "FERMT3", "ACTN1", "PDLIM1", "ILK", "CTSA", "TPM4", 
                                        "SNCA", "ODC1", "GRAP2", "GAS2L1", "PTGS1", "NGFRAP1", "TREML1", 
                                        "MYL9", "TPM1", "F13A1", "CA2", "MPP1", "RGS18", "CLU", "AP001189.4", 
                                        "MMD", "GP9", "PPBP", "SPARC", "TUBB1", "TMEM40", "HIST1H2AC", 
                                        "PF4", "GNG11", "SDPR", "ACRBP", "PTCRA", " ", " ", "RPL6", "RPSAP58", 
                                        "RPLP1", "RPS9", "RPL26", "RPL27", "NACA", "RPL17", "RPL7A", 
                                        "EEF1B2", "RPL5", "RPL36", "RPS28", "RPS20", "RPL18", "RPL15", 
                                        "BTG1", "TPT1", "RPL35A", "RPS26", "RPS11", "RPS13", "RPS10", 
                                        "RPS29", "RPL28", "RPS8", "RPS4X", "RPL31", "EEF1A1", "RPL12", 
                                        "RPS16", "RPL27A", "RPS23", "RPS15A", "RPS3A", "RPL14", "RPL29", 
                                        "RPL3", "RPS18", "RPLP2", "RPL32", "RPL11", "RPL18A", "RPS12", 
                                        "RPS3", "RPS19", "RPS25", "RPL21", "RPL23A", "RPL9", "MALAT1"
  ))
  
  # Checking title of the heatmap
  htmp_title <- htmp@annotations@listData$col_title2@data
  expect_equal(htmp_title, "Expression heatmap")
  
  # Checking x lab of the heatmap
  htmp_xlab <- htmp@annotations@listData$col_title@data
  expect_equal(htmp_xlab, "Samples")
  
  # Checking y lab of the heatmap
  htmp_ylab <- htmp@annotations@listData$row_title@data
  expect_equal(htmp_ylab, "Genes")
  
  # Checking colorbar name of the heatmap
  htmp_colors_name <- htmp@colorbars@listData$`Exp. level`@title
  expect_equal(htmp_colors_name, "Exp. level")
  
  # Checking if colorbar is shown of the heatmap
  htmp_show_colorbar <- htmp@plots@listData$`Exp. level`@show_colorbar
  expect_true(htmp_show_colorbar)
  
  # Checking if dendrogram is showed
  htmp_show_dendro <- names(htmp@shapes@listData)
  expect_equal(htmp_show_dendro, NULL)
})



test_that("Checking plot_heatmap() using top genes and cell clusters", {
  
  set.seed(123)
  cell_clusters <- Idents(pbmc3k_medium)
  res <- top_genes(res, top = 20)
  
  htmp <- plot_heatmap(
    object = res,
    cell_clusters = cell_clusters,
    use_top_genes = TRUE,
    show_dendro = TRUE,
    line_size_vertical = 2,
    line_size_horizontal = 2,
    name = "Expression heatmap",
    show_legend = TRUE,
    xlab = "Samples",
    ylab = "Genes"
  )
  
  htmp_matrix <- htmp@plots@listData$`Exp. level`@data
  htmp_matrix_wo_na <- htmp_matrix[!is.na(htmp_matrix)]
  
  # Check heatmap format

  # Check dim of heatmap matrix
  expect_equal(ncol(htmp_matrix), 377)
  expect_equal(nrow(htmp_matrix), 230)
  expect_equal(dim(htmp_matrix), c(230, 377))
  
  # Check NA rows added by line_size_horizontal
  na_rows <- which(apply(htmp_matrix, 1, function(x) all(is.na(x))))
  expect_equal(length(na_rows), 28)
  expect_equal(as.vector(na_rows), c(6L, 7L, 13L, 14L, 21L, 22L, 30L, 31L, 39L, 40L, 53L, 54L, 69L, 
                                     70L, 85L, 86L, 101L, 102L, 121L, 122L, 143L, 144L, 165L, 166L, 
                                     187L, 188L, 209L, 210L))
  
  # Check NA rows added by line_size_vertical
  na_cols <- which(apply(htmp_matrix, 2, function(x) all(is.na(x))))
  expect_equal(length(na_cols), 16)
  expect_equal(as.vector(na_cols), c(48L, 49L, 97L, 98L, 147L, 148L, 198L, 199L, 246L, 247L, 295L, 
                                     296L, 334L, 335L, 362L, 363L))
  
  # Check some measurements of heatmap matrix
  expect_equal(round(sum(htmp_matrix_wo_na), 2), -3394.67)
  expect_equal(round(mean(htmp_matrix_wo_na), 2), -0.05)
  expect_equal(round(median(htmp_matrix_wo_na), 2), -0.05)
  expect_equal(round(sd(htmp_matrix_wo_na), 2), 0.48)
  
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["0%"]), 2), -1)
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["25%"]), 2), -0.22)
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["50%"]), 2), -0.05)
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["75%"]), 2), -0.02)
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["100%"]), 2), 1)
  
  # Checking order of columns
  expect_equal(head(colnames(htmp_matrix),10), c("GATCTACTGGTGAG-1", "ACAGTGACTCACCC-1", "AGACGTACAGAGGC-1", 
                                                 "GACGTAACCTGTGA-1", "TATACAGATCCAGA-1", "CGGATAACAGCTCA-1", "TTACCATGGTTGAC-1", 
                                                 "TCCCACGATCATTC-1", "ATAGCGTGCAGATC-1", "TGTAGGTGTGCTGA-1"))
  
  # Checking order of columns
  expect_equal(tail(colnames(htmp_matrix),10), c("ACCTGAGATATCGG-1", "GACGCTCTCTCTCG-1", "GCGCATCTGGTTAC-1", 
                                                 "AGTCTTACTTCGGA-1", "GGCATATGGGGAGT-1", "TTACGTACGTTCAG-1", "GGAACACTTCAGAC-1", 
                                                 "ATCATCTGACACCA-1", "ACGAACTGGCTATG-1", "TAACACCTTGTTTC-1"))
  
  # Checking order of rows
  expect_equal(rownames(htmp_matrix), c("PKIG", "ISG20", "MS4A1", "CD79B", "LINC00926", " ", " ", "CAPZA2", 
                                        "MEST", "MTURN", "ZHX1-C8ORF76", "RP11-879F14.2", " ", " ", "CLIC1", 
                                        "HES4", "ABI3", "IFITM2", "RHOC", "FCGR3A", " ", " ", "AES", 
                                        "CD27", "IL7R", "CD3E", "LDHB", "IL32", "CD3D", " ", " ", "PTPN18", 
                                        "SMIM5", "PRUNE", "ANKRD9", "GFI1B", "LCN2", "ENKUR", " ", " ", 
                                        "XPNPEP1", "AC092295.7", "TJP2", "PRKAR2B", "TMEM91", "SPHK1", 
                                        "SLC40A1", "TUBA1C", "FAM63A", "C2orf88", "CLDN5", "HIST1H2BJ", 
                                        " ", " ", "SLC25A6", "RNASET2", "CD37", "FAM26F", "HLA-DMB", 
                                        "LY86", "HLA-DQB1", "HLA-DQA1", "CD74", "HLA-DPB1", "HLA-DRB5", 
                                        "HLA-DPA1", "HLA-DRA", "HLA-DRB1", " ", " ", "FN3K", "MGLL", 
                                        "C19orf33", "NFE2", "MOB3B", "TRIM58", "SERPINE2", "ZNF778", 
                                        "TMCC2", "AC147651.3", "FAM212B", "SSX2IP", "SMOX", "LGALSL", 
                                        " ", " ", "H1F0", "GLUL", "HIST2H2BE", "MLH3", "C1orf198", "SENCR", 
                                        "SEPT4", "FNTB", "RP11-367G6.3", "HIST1H2BD", "SCGB1C1", "CMTM5", 
                                        "GP1BA", "CLEC1B", " ", " ", "PPP1R14A", "WIPI1", "GLA", "ABHD16A", 
                                        "TGFB1I1", "RP11-359I18.5", "NEXN", "ENDOD1", "ALOX12", "ITGB3", 
                                        "ARG2", "SCFD2", "ABCC3", "CTNNAL1", "FAM212A", "ATP9A", "SPOCD1", 
                                        "ARHGAP21", " ", " ", "LYL1", "BTK", "FAM110A", "PYGL", "CLIC4", 
                                        "RILP", "DAB2", "KIAA0513", "SH3BGRL2", "RAB27B", "KIFC3", "PCP2", 
                                        "HEXIM2", "DPY19L1", "GATA2", "SEC14L5", "EGFL7", "hsa-mir-1199", 
                                        "HEMGN", "AC137932.6", " ", " ", "CCL5", "IGFBP7", "XCL1", "KLRF1", 
                                        "AKR1C3", "CD7", "HLA-C", "HOPX", "IL2RB", "CD247", "CLIC3", 
                                        "CTSW", "SPON2", "CCL4", "CST7", "GZMA", "FGFBP2", "GNLY", "PRF1", 
                                        "GZMB", " ", " ", "BLVRA", "CAPG", "C1orf162", "TGFBI", "SLC7A7", 
                                        "IFI6", "MAFB", "CPVL", "CSTB", "IGSF6", "AP2S1", "APOBEC3A", 
                                        "TNFSF13B", "GABARAP", "LGALS2", "CFD", "S100A11", "CTSS", "COTL1", 
                                        "TYMP", " ", " ", "MPP1", "PTGS1", "MYL9", "MMD", "TMEM40", "CA2", 
                                        "F13A1", "TREML1", "PTCRA", "RGS18", "AP001189.4", "CLU", "GP9", 
                                        "HIST1H2AC", "TUBB1", "PPBP", "SPARC", "PF4", "GNG11", "SDPR", 
                                        " ", " ", "RPL29", "RPS4X", "RPS8", "EEF1A1", "RPS25", "RPS19", 
                                        "RPL28", "RPLP2", "RPS15A", "RPL12", "RPS3", "RPL23A", "RPL9", 
                                        "RPS16", "RPS23", "RPS18", "RPL18A", "RPS12", "RPL32", "RPL11"
  ))
  
  # Checking title of the heatmap
  htmp_title <- htmp@annotations@listData$col_title2@data
  expect_equal(htmp_title, "Expression heatmap")
  
  # Checking x lab of the heatmap
  htmp_xlab <- htmp@annotations@listData$col_title@data
  expect_equal(htmp_xlab, "Samples")
  
  # Checking y lab of the heatmap
  htmp_ylab <- htmp@annotations@listData$row_title@data
  expect_equal(htmp_ylab, "Genes")
  
  # Checking colorbar name of the heatmap
  htmp_colors_name <- htmp@colorbars@listData$`Exp. level`@title
  expect_equal(htmp_colors_name, "Exp. level")
  
  # Checking if colorbar is shown of the heatmap
  htmp_show_colorbar <- htmp@plots@listData$`Exp. level`@show_colorbar
  expect_true(htmp_show_colorbar)
  
  # Checking if dendrogram is showed
  htmp_show_dendro <- names(htmp@shapes@listData)
  expect_equal(htmp_show_dendro, NULL)
})


# ==============================================================================
# ==============================================================================


test_that("Checking plot_heatmap() displaying gene clusters 1 and 4", {
  set.seed(123)
  cell_clusters <- Idents(pbmc3k_medium)
  res <- top_genes(res, top = 20)
  
  htmp <- plot_heatmap(
    object = res[c(1,4), ],
    cell_clusters = cell_clusters,
    use_top_genes = TRUE,
    show_dendro = TRUE,
    line_size_vertical = 2,
    line_size_horizontal = 2,
    name = "Expression heatmap",
    show_legend = TRUE,
    xlab = "Samples",
    ylab = "Genes"
  )
  
  htmp_matrix <- htmp@plots@listData$`Exp. level`@data
  htmp_matrix_wo_na <- htmp_matrix[!is.na(htmp_matrix)]
  
  # Check heatmap format

  # Check dim of heatmap matrix
  expect_equal(ncol(htmp_matrix), 377)
  expect_equal(nrow(htmp_matrix), 42)
  expect_equal(dim(htmp_matrix), c(42, 377))
  
  # Check NA rows added by line_size_horizontal
  na_rows <- which(apply(htmp_matrix, 1, function(x) all(is.na(x))))
  expect_equal(length(na_rows), 2)
  expect_equal(as.vector(na_rows), c(21, 22))
  
  # Check NA rows added by line_size_vertical
  na_cols <- which(apply(htmp_matrix, 2, function(x) all(is.na(x))))
  expect_equal(length(na_cols), 16)
  expect_equal(as.vector(na_cols), c(48L, 49L, 97L, 98L, 147L, 148L, 198L, 199L, 246L, 247L, 295L, 
                                     296L, 334L, 335L, 362L, 363L))
  
  # Check some measurements of heatmap matrix
  expect_equal(round(sum(htmp_matrix_wo_na), 2), -509.92)
  expect_equal(round(mean(htmp_matrix_wo_na), 2), -0.04)
  expect_equal(round(median(htmp_matrix_wo_na), 2), -0.11)
  expect_equal(round(sd(htmp_matrix_wo_na), 2), 0.57)
  
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["0%"]), 2), -1)
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["25%"]), 2), -0.49)
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["50%"]), 2), -0.11)
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["75%"]), 2), 0.4)
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["100%"]), 2), 1)
  
  # Checking order of columns
  expect_equal(head(colnames(htmp_matrix),10), c("GATCTACTGGTGAG-1", "ACAGTGACTCACCC-1", "AGACGTACAGAGGC-1", 
                                                 "GACGTAACCTGTGA-1", "TATACAGATCCAGA-1", "CGGATAACAGCTCA-1", "TTACCATGGTTGAC-1", 
                                                 "TCCCACGATCATTC-1", "ATAGCGTGCAGATC-1", "TGTAGGTGTGCTGA-1"))
  
  # Checking order of columns
  expect_equal(tail(colnames(htmp_matrix),10), c("ACCTGAGATATCGG-1", "GACGCTCTCTCTCG-1", "GCGCATCTGGTTAC-1", 
                                                 "AGTCTTACTTCGGA-1", "GGCATATGGGGAGT-1", "TTACGTACGTTCAG-1", "GGAACACTTCAGAC-1", 
                                                 "ATCATCTGACACCA-1", "ACGAACTGGCTATG-1", "TAACACCTTGTTTC-1"))
  
  # Checking order of rows
  expect_equal(rownames(htmp_matrix), c("CCL5", "IGFBP7", "XCL1", "KLRF1", "AKR1C3", "CD7", "HLA-C", 
                                        "HOPX", "IL2RB", "CD247", "CLIC3", "CTSW", "SPON2", "CCL4", "CST7", 
                                        "GZMA", "FGFBP2", "GNLY", "PRF1", "GZMB", " ", " ", "RPL29", 
                                        "RPS4X", "RPS8", "EEF1A1", "RPS25", "RPS19", "RPL28", "RPLP2", 
                                        "RPS15A", "RPL12", "RPS3", "RPL23A", "RPL9", "RPS16", "RPS23", 
                                        "RPS18", "RPL18A", "RPS12", "RPL32", "RPL11"))
  
  # Checking title of the heatmap
  htmp_title <- htmp@annotations@listData$col_title2@data
  expect_equal(htmp_title, "Expression heatmap")
  
  # Checking x lab of the heatmap
  htmp_xlab <- htmp@annotations@listData$col_title@data
  expect_equal(htmp_xlab, "Samples")
  
  # Checking y lab of the heatmap
  htmp_ylab <- htmp@annotations@listData$row_title@data
  expect_equal(htmp_ylab, "Genes")
  
  # Checking colorbar name of the heatmap
  htmp_colors_name <- htmp@colorbars@listData$`Exp. level`@title
  expect_equal(htmp_colors_name, "Exp. level")
  
  # Checking if colorbar is shown of the heatmap
  htmp_show_colorbar <- htmp@plots@listData$`Exp. level`@show_colorbar
  expect_true(htmp_show_colorbar)
  
  # Checking if dendrogram is showed
  htmp_show_dendro <- names(htmp@shapes@listData)
  expect_equal(htmp_show_dendro, NULL)
})


# ==============================================================================
# ==============================================================================


test_that("Checking plot_heatmap() displaying only gene cluster 3", {
  
  set.seed(123)
  cell_clusters <- Idents(pbmc3k_medium)
  res <- top_genes(res, top = 20)
  
  htmp <- plot_heatmap(
    object = res[3, ],
    cell_clusters = cell_clusters,
    use_top_genes = TRUE,
    show_dendro = TRUE,
    line_size_vertical = 2,
    line_size_horizontal = 2,
    name = "Expression heatmap",
    show_legend = TRUE,
    xlab = "Samples",
    ylab = "Genes"
  )
  
  htmp_matrix <- htmp@plots@listData$`Exp. level`@data
  htmp_matrix_wo_na <- htmp_matrix[!is.na(htmp_matrix)]
  
  # Check heatmap format

  # Check dim of heatmap matrix
  expect_equal(ncol(htmp_matrix), 377)
  expect_equal(nrow(htmp_matrix), 20)
  expect_equal(dim(htmp_matrix), c(20, 377))
  
  # Check NA rows added by line_size_horizontal
  na_rows <- which(apply(htmp_matrix, 1, function(x) all(is.na(x))))
  expect_equal(length(na_rows), 0)
  expect_equal(as.vector(na_rows), integer(0))
  
  # Check NA rows added by line_size_vertical
  na_cols <- which(apply(htmp_matrix, 2, function(x) all(is.na(x))))
  expect_equal(length(na_cols), 16)
  expect_equal(as.vector(na_cols), c(48L, 49L, 97L, 98L, 147L, 148L, 198L, 199L, 246L, 247L, 295L, 
                                     296L, 334L, 335L, 362L, 363L))
  
  # Check some measurements of heatmap matrix
  expect_equal(round(sum(htmp_matrix_wo_na), 2), -663.16)
  expect_equal(round(mean(htmp_matrix_wo_na), 2), -0.09)
  expect_equal(round(median(htmp_matrix_wo_na), 2), -0.27)
  expect_equal(round(sd(htmp_matrix_wo_na), 2), 0.72)
  
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["0%"]), 2), -1)
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["25%"]), 2), -0.59)
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["50%"]), 2), -0.27)
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["75%"]), 2), 0.78)
  expect_equal(round(sum(quantile(htmp_matrix_wo_na)["100%"]), 2), 1)
  
  # Checking order of columns
  expect_equal(head(colnames(htmp_matrix),10), c("GATCTACTGGTGAG-1", "ACAGTGACTCACCC-1", "AGACGTACAGAGGC-1", 
                                                 "GACGTAACCTGTGA-1", "TATACAGATCCAGA-1", "CGGATAACAGCTCA-1", "TTACCATGGTTGAC-1", 
                                                 "TCCCACGATCATTC-1", "ATAGCGTGCAGATC-1", "TGTAGGTGTGCTGA-1"))
  
  # Checking order of columns
  expect_equal(tail(colnames(htmp_matrix),10), c("ACCTGAGATATCGG-1", "GACGCTCTCTCTCG-1", "GCGCATCTGGTTAC-1", 
                                                 "AGTCTTACTTCGGA-1", "GGCATATGGGGAGT-1", "TTACGTACGTTCAG-1", "GGAACACTTCAGAC-1", 
                                                 "ATCATCTGACACCA-1", "ACGAACTGGCTATG-1", "TAACACCTTGTTTC-1"))
  
  # Checking order of rows
  expect_equal(rownames(htmp_matrix), c("BLVRA", "CAPG", "C1orf162", "TGFBI", "SLC7A7", "IFI6", "MAFB", 
                                        "CPVL", "CSTB", "IGSF6", "AP2S1", "APOBEC3A", "TNFSF13B", "GABARAP", 
                                        "LGALS2", "CFD", "S100A11", "CTSS", "COTL1", "TYMP"))
  
  # Checking title of the heatmap
  htmp_title <- htmp@annotations@listData$col_title2@data
  expect_equal(htmp_title, "Expression heatmap")
  
  # Checking x lab of the heatmap
  htmp_xlab <- htmp@annotations@listData$col_title@data
  expect_equal(htmp_xlab, "Samples")
  
  # Checking y lab of the heatmap
  htmp_ylab <- htmp@annotations@listData$row_title@data
  expect_equal(htmp_ylab, "Genes")
  
  # Checking colorbar name of the heatmap
  htmp_colors_name <- htmp@colorbars@listData$`Exp. level`@title
  expect_equal(htmp_colors_name, "Exp. level")
  
  # Checking if colorbar is shown of the heatmap
  htmp_show_colorbar <- htmp@plots@listData$`Exp. level`@show_colorbar
  expect_true(htmp_show_colorbar)
  
  # Checking if dendrogram is showed
  htmp_show_dendro <- names(htmp@shapes@listData)
  expect_equal(htmp_show_dendro, NULL)
})




testthat::test_that("Check plot_heatmap, subsetting and clustering", {
  # Set verbosity to 0
  set_verbosity(0)
  
  # Load seurat object
  data(pbmc_small, package = "SeuratObject")
  
  ## Select informative genes
  clust_set <- select_genes(data=pbmc_small,
                            distance_method="pearson",
                            k=10,
                            row_sum=-Inf,
                            noise_level=0.95,
                            fdr = 1e-6)
  
  ## Cluster genes
  clust_set <- gene_clustering(object = clust_set,
                               inflation = 1.2,
                               keep_nn = FALSE,
                               k = 5,
                               threads = 1)
  
  ident <- Seurat::Idents(pbmc_small)
  
  # Check color arg
  expect_error(plot_heatmap(clust_set, cell_clusters=Seurat::Idents(pbmc_small), colors_cell_clusters = "green"))
  p <- plot_heatmap(clust_set[, names(ident[ident==2])], cell_clusters=Seurat::Idents(pbmc_small), colors_cell_clusters = "red")
  testthat::expect_error(print(p), NA)
  expect_equal(unique(p@colorbars@listData$Ident.@colors), "red")
  
  # Check if seurat identities is considered
  
  p <- plot_heatmap(clust_set[,colnames(clust_set@data)[1:80]], cell_clusters=Seurat::Idents(pbmc_small))
  testthat::expect_error(print(p), NA)
  expect_equal(dim(p@plots@listData$`Exp. level`@data), c(178, 86)) # 80 + 6 NA columns
  
  p <- plot_heatmap(clust_set[,colnames(clust_set@data)[1:10]], cell_clusters=Seurat::Idents(pbmc_small))
  testthat::expect_error(print(p), NA)
  expect_equal(dim(p@plots@listData$`Exp. level`@data), c(178, 10))
  
  p <- plot_heatmap(clust_set[,colnames(clust_set@data)[1:20]], cell_clusters=Seurat::Idents(pbmc_small))
  testthat::expect_error(print(p), NA)
  expect_equal(dim(p@plots@listData$`Exp. level`@data), c(178, 23)) # 20 + 3 NA columns
  
  p <- plot_heatmap(clust_set[1,colnames(clust_set@data)[1:10]], cell_clusters=Seurat::Idents(pbmc_small))
  testthat::expect_error(print(p), NA)
  expect_equal(dim(p@plots@listData$`Exp. level`@data), c(91, 10))
  expect_equal(dim(p@plots@listData$`Exp. level`@data), c(length(clust_set@gene_clusters$`1`), 10))
  
  p <- plot_heatmap(clust_set[2 ,colnames(clust_set@data)[1:10]], cell_clusters=Seurat::Idents(pbmc_small))
  testthat::expect_error(print(p), NA)
  expect_equal(dim(p@plots@listData$`Exp. level`@data), c(52, 10))
  expect_equal(dim(p@plots@listData$`Exp. level`@data), c(length(clust_set@gene_clusters$`2`), 10))
})




testthat::test_that("Test non interactive heatmaps.", {
  # Set verbosity to 0
  set_verbosity(0)
  
  # Load seurat object
  data(pbmc_small, package = "SeuratObject")
  
  ## Select informative genes
  clust_set <- select_genes(data=pbmc_small,
                            distance_method="pearson",
                            k=10,
                            row_sum=-Inf,
                            noise_level=0.95,
                            fdr = 1e-6)
  
  ## Cluster genes
  clust_set <- gene_clustering(object = clust_set,
                               inflation = 1.2,
                               keep_nn = FALSE,
                               k = 5,
                               threads = 1)
  
  ident <- Seurat::Idents(pbmc_small)
  
  n_gene_clusters <- length(clust_set@gene_clusters) 
  n_cell_clusters <- 3
  p_int <- plot_heatmap(clust_set, cell_clusters=Seurat::Idents(pbmc_small), 
                        colors_cell_clusters =  c("green", "red", "yellow"))
  p_noi <- plot_heatmap(clust_set, cell_clusters=Seurat::Idents(pbmc_small), 
                        colors_cell_clusters =  c("green", "red", "yellow"),
                        interactive = FALSE)
  p_int_mat <- htmp_matrix <- p_int@plots@listData$`Exp. level`@data
  p_noi_rown <- p_noi[["gtable"]][["grobs"]][[2]][["label"]]
  testthat::expect_equal(p_noi_rown, rev(rownames(p_int_mat)))
  
  p_noi_colnames <- names(p_noi[["gtable"]][["grobs"]][[3]][["gp"]][["fill"]][,1])
  expect_true(all(colnames(p_int_mat) == p_noi_colnames))
  expect_true(length(grep(" ", p_noi_colnames)) == (n_cell_clusters-1) * 3)
  expect_true(length(grep(" ", p_noi_rown)) == (n_gene_clusters-1) * 3)
  
  n_gene_clusters <- 1 
  n_cell_clusters <- 3
  p_noi <- plot_heatmap(clust_set[1,], cell_clusters=Seurat::Idents(pbmc_small), 
                        colors_cell_clusters =  c("green", "red", "yellow"),
                        interactive = FALSE)
  p_noi_rown <- p_noi[["gtable"]][["grobs"]][[2]][["label"]]
  expect_true(length(grep(" ", p_noi_colnames)) == (n_cell_clusters-1) * 3)
  expect_true(length(grep(" ", p_noi_rown)) == (n_gene_clusters-1) * 3)
  
  n_gene_clusters <- 2 
  n_cell_clusters <- 3
  p_noi <- plot_heatmap(clust_set[c(1, 3), colnames(pbmc_small)], cell_clusters=Seurat::Idents(pbmc_small), 
                        colors_cell_clusters =  c("green", "red", "yellow"),
                        interactive = FALSE)
  p_noi_rown <- p_noi[["gtable"]][["grobs"]][[2]][["label"]]
  expect_true(length(grep(" ", p_noi_colnames)) == (n_cell_clusters-1) * 3)
  expect_true(length(grep(" ", p_noi_rown)) == (n_gene_clusters-1) * 3)
})



