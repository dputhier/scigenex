# Set verbosity to 0
set_verbosity(0)
library(Seurat)
library(testthat)
load_example_dataset("7871581/files/pbmc3k_medium_clusters")

res <- pbmc3k_medium_clusters

test_that("Cheking top_gene()", {

  # ========================================
  # Top 20
  res_20 <- top_genes(res,top = 20)
  
  # Test top genes list in all cluster
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
  expect_equal(unlist(res_20@top_genes, use.names = FALSE), gene_name_to_check)

  
  # Test top genes list in cluster 1
  gene_name_to_check <- c("RPL11", "RPL32", "RPS12", "RPL18A", "RPS18", "RPS23", "RPS16", 
                          "RPL9", "RPL23A", "RPS3", "RPL12", "RPS15A", "RPLP2", "RPL28", 
                          "RPS19", "RPS25", "EEF1A1", "RPS8", "RPS4X", "RPL29")
  expect_equal(res_20@top_genes$`1`, gene_name_to_check)
  expect_equal(length(res_20@top_genes$`1`), 20)
  
  
  # Test top genes list in cluster 2 and 3
  gene_name_to_check <- c("SDPR", "GNG11", "PF4", "SPARC", "PPBP", "TUBB1", "HIST1H2AC", 
                          "GP9", "CLU", "AP001189.4", "RGS18", "PTCRA", "TREML1", "F13A1", 
                          "CA2", "TMEM40", "MMD", "MYL9", "PTGS1", "MPP1", "TYMP", "COTL1", 
                          "CTSS", "S100A11", "CFD", "LGALS2", "GABARAP", "TNFSF13B", "APOBEC3A", 
                          "AP2S1", "IGSF6", "CSTB", "CPVL", "MAFB", "IFI6", "SLC7A7", "TGFBI", 
                          "C1orf162", "CAPG", "BLVRA")
  expect_equal(c(res_20@top_genes$`2`, res_20@top_genes$`3`),
               gene_name_to_check)
  expect_equal(length(c(res_20@top_genes$`2`, res_20@top_genes$`3`)), 40)
  
  
  # ========================================
  # Top 10
  res_10 <- top_genes(res, top = 10)
  
  # Test top genes list in all cluster
  gene_name_to_check <- c("RPL11", "RPL32", "RPS12", "RPL18A", "RPS18", "RPS23", "RPS16", 
                          "RPL9", "RPL23A", "RPS3", "SDPR", "GNG11", "PF4", "SPARC", "PPBP", 
                          "TUBB1", "HIST1H2AC", "GP9", "CLU", "AP001189.4", "TYMP", "COTL1", 
                          "CTSS", "S100A11", "CFD", "LGALS2", "GABARAP", "TNFSF13B", "APOBEC3A", 
                          "AP2S1", "GZMB", "PRF1", "GNLY", "FGFBP2", "GZMA", "CST7", "CCL4", 
                          "SPON2", "CTSW", "CLIC3", "AC137932.6", "HEMGN", "hsa-mir-1199", 
                          "EGFL7", "SEC14L5", "GATA2", "DPY19L1", "HEXIM2", "PCP2", "KIFC3", 
                          "ARHGAP21", "SPOCD1", "ATP9A", "FAM212A", "CTNNAL1", "ABCC3", 
                          "SCFD2", "ARG2", "ITGB3", "ALOX12", "CLEC1B", "GP1BA", "CMTM5", 
                          "SCGB1C1", "HIST1H2BD", "RP11-367G6.3", "FNTB", "SEPT4", "SENCR", 
                          "C1orf198", "LGALSL", "SMOX", "SSX2IP", "FAM212B", "AC147651.3", 
                          "TMCC2", "ZNF778", "SERPINE2", "TRIM58", "MOB3B", "HLA-DRB1", 
                          "HLA-DRA", "HLA-DPA1", "HLA-DRB5", "HLA-DPB1", "CD74", "HLA-DQA1", 
                          "HLA-DQB1", "LY86", "HLA-DMB", "HIST1H2BJ", "CLDN5", "C2orf88", 
                          "FAM63A", "TUBA1C", "SLC40A1", "SPHK1", "TMEM91", "PRKAR2B", 
                          "TJP2", "ENKUR", "LCN2", "GFI1B", "ANKRD9", "PRUNE", "SMIM5", 
                          "PTPN18", "CD3D", "IL32", "LDHB", "CD3E", "IL7R", "CD27", "AES", 
                          "FCGR3A", "RHOC", "IFITM2", "ABI3", "HES4", "CLIC1", "RP11-879F14.2", 
                          "ZHX1-C8ORF76", "MTURN", "MEST", "CAPZA2", "LINC00926", "CD79B", 
                          "MS4A1", "ISG20", "PKIG")
  expect_equal(unlist(res_10@top_genes, use.names = FALSE), gene_name_to_check)
  expect_equal(length(unlist(res_10@top_genes, use.names = FALSE)), 130)
  
  
  # Test top genes list in cluster 1
  gene_name_to_check <- c("RPL11", "RPL32", "RPS12", "RPL18A", "RPS18", "RPS23", "RPS16", 
                          "RPL9", "RPL23A", "RPS3")
  expect_equal(res_10@top_genes$`1`, gene_name_to_check)
  expect_equal(length(res_10@top_genes$`1`), 10)
  
  
  # Test top genes list in cluster 2 and 3
  gene_name_to_check <- c("SDPR", "GNG11", "PF4", "SPARC", "PPBP", "TUBB1", "HIST1H2AC", 
                          "GP9", "CLU", "AP001189.4", "TYMP", "COTL1", "CTSS", "S100A11", 
                          "CFD", "LGALS2", "GABARAP", "TNFSF13B", "APOBEC3A", "AP2S1")
  expect_equal(
    c(res_10@top_genes$`2`, res_10@top_genes$`3`),
    gene_name_to_check
  )
  expect_equal(length(c(res_10@top_genes$`2`, res_10@top_genes$`3`)), 20)
  
  # ========================================
  # Top 5
  res_5 <- top_genes(res, top = 5)
  
  # Test top genes list in all cluster
  gene_name_to_check <- c("RPL11", "RPL32", "RPS12", "RPL18A", "RPS18", "SDPR", "GNG11", 
                          "PF4", "SPARC", "PPBP", "TYMP", "COTL1", "CTSS", "S100A11", "CFD", 
                          "GZMB", "PRF1", "GNLY", "FGFBP2", "GZMA", "AC137932.6", "HEMGN", 
                          "hsa-mir-1199", "EGFL7", "SEC14L5", "ARHGAP21", "SPOCD1", "ATP9A", 
                          "FAM212A", "CTNNAL1", "CLEC1B", "GP1BA", "CMTM5", "SCGB1C1", 
                          "HIST1H2BD", "LGALSL", "SMOX", "SSX2IP", "FAM212B", "AC147651.3", 
                          "HLA-DRB1", "HLA-DRA", "HLA-DPA1", "HLA-DRB5", "HLA-DPB1", "HIST1H2BJ", 
                          "CLDN5", "C2orf88", "FAM63A", "TUBA1C", "ENKUR", "LCN2", "GFI1B", 
                          "ANKRD9", "PRUNE", "CD3D", "IL32", "LDHB", "CD3E", "IL7R", "FCGR3A", 
                          "RHOC", "IFITM2", "ABI3", "HES4", "RP11-879F14.2", "ZHX1-C8ORF76", 
                          "MTURN", "MEST", "CAPZA2", "LINC00926", "CD79B", "MS4A1", "ISG20", 
                          "PKIG")
  expect_equal(unlist(res_5@top_genes, use.names = FALSE), gene_name_to_check)
  expect_equal(length(unlist(res_5@top_genes, use.names = FALSE)), 75)
  
  
  # Test top genes list in cluster 1
  gene_name_to_check <- c("RPL11", "RPL32", "RPS12", "RPL18A", "RPS18")
  expect_equal(res_5@top_genes$`1`, gene_name_to_check)
  
  
  # Test top genes list in cluster 2 and 3
  gene_name_to_check <- c("SDPR", "GNG11", "PF4", "SPARC", "PPBP", "TYMP", "COTL1", "CTSS", 
                          "S100A11", "CFD")
  expect_equal(c(res_5@top_genes$`2`, res_5@top_genes$`3`), gene_name_to_check)
  
  
  
  # ========================================
  # Top 100
  res_100 <- suppressWarnings(top_genes(res, top = 100))
  
  # Test top genes list in all cluster
  gene_name_to_check <- c("RPL11", "RPL32", "RPS12", "RPL18A", "RPS18", "RPS23", "RPS16", 
                          "RPL9", "RPL23A", "RPS3", "RPL12", "RPS15A", "RPLP2", "RPL28", 
                          "RPS19", "RPS25", "EEF1A1", "RPS8", "RPS4X", "RPL29", "RPL27A", 
                          "RPL18", "RPL14", "RPL3", "RPL35A", "RPL21", "RPLP1", "RPS28", 
                          "RPL6", "RPS10", "RPL26", "RPL7A", "RPL5", "RPL27", "RPL15", 
                          "RPS13", "RPL17", "RPL31", "RPS3A", "RPS29", "RPL36", "MALAT1", 
                          "RPS20", "TPT1", "RPS26", "RPS9", "EEF1B2", "RPS11", "NACA", 
                          "RPSAP58", "BTG1", "SDPR", "GNG11", "PF4", "SPARC", "PPBP", "TUBB1", 
                          "HIST1H2AC", "GP9", "CLU", "AP001189.4", "RGS18", "PTCRA", "TREML1", 
                          "F13A1", "CA2", "TMEM40", "MMD", "MYL9", "PTGS1", "MPP1", "ACTN1", 
                          "ACRBP", "TPM4", "TPM1", "SNCA", "NGFRAP1", "GAS2L1", "GRAP2", 
                          "CTSA", "ILK", "ODC1", "SNN", "PDLIM1", "FERMT3", "HIST1H2BK", 
                          "MARCH2", "PARVB", "PLA2G12A", "MFSD1", "PVALB", "CD151", "TPTEP1", 
                          "GSN", "YWHAH", "RAP1B", "ASAH1", "APP", "MYL12A", "H2AFJ", "TYMP", 
                          "COTL1", "CTSS", "S100A11", "CFD", "LGALS2", "GABARAP", "TNFSF13B", 
                          "APOBEC3A", "AP2S1", "IGSF6", "CSTB", "CPVL", "MAFB", "IFI6", 
                          "SLC7A7", "TGFBI", "C1orf162", "CAPG", "BLVRA", "BID", "MNDA", 
                          "LILRB2", "S100A10", "CD14", "WARS", "NAAA", "NUP214", "ATG3", 
                          "CNPY3", "CYBB", "GAPDH", "IFNGR2", "ALDH2", "EIF4EBP1", "CD300LF", 
                          "GCA", "GRINA", "LILRB4", "SOD2", "CTSZ", "CD302", "CAMK1", "POU2F2", 
                          "FCGR2A", "GZMB", "PRF1", "GNLY", "FGFBP2", "GZMA", "CST7", "CCL4", 
                          "SPON2", "CTSW", "CLIC3", "CD247", "IL2RB", "HOPX", "HLA-C", 
                          "CD7", "AKR1C3", "KLRF1", "XCL1", "IGFBP7", "CCL5", "CD99", "PTPRCAP", 
                          "PLEK", "EFHD2", "AC137932.6", "HEMGN", "hsa-mir-1199", "EGFL7", 
                          "SEC14L5", "GATA2", "DPY19L1", "HEXIM2", "PCP2", "KIFC3", "RAB27B", 
                          "SH3BGRL2", "KIAA0513", "DAB2", "RILP", "CLIC4", "PYGL", "FAM110A", 
                          "BTK", "LYL1", "ARHGAP21", "SPOCD1", "ATP9A", "FAM212A", "CTNNAL1", 
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
  expect_equal(unlist(res_100@top_genes, use.names = FALSE), gene_name_to_check)
  
  # Test top genes list in cluster 1
  gene_name_to_check <- c("RPL11", "RPL32", "RPS12", "RPL18A", "RPS18", "RPS23", "RPS16", 
                          "RPL9", "RPL23A", "RPS3", "RPL12", "RPS15A", "RPLP2", "RPL28", 
                          "RPS19", "RPS25", "EEF1A1", "RPS8", "RPS4X", "RPL29", "RPL27A", 
                          "RPL18", "RPL14", "RPL3", "RPL35A", "RPL21", "RPLP1", "RPS28", 
                          "RPL6", "RPS10", "RPL26", "RPL7A", "RPL5", "RPL27", "RPL15", 
                          "RPS13", "RPL17", "RPL31", "RPS3A", "RPS29", "RPL36", "MALAT1", 
                          "RPS20", "TPT1", "RPS26", "RPS9", "EEF1B2", "RPS11", "NACA", 
                          "RPSAP58", "BTG1")
  expect_equal(res_100@top_genes$`1`, gene_name_to_check)
  
  
  # Test top genes list in cluster 2 and 3
  gene_name_to_check <- c("SDPR", "GNG11", "PF4", "SPARC", "PPBP", "TUBB1", "HIST1H2AC", 
                         "GP9", "CLU", "AP001189.4", "RGS18", "PTCRA", "TREML1", "F13A1", 
                         "CA2", "TMEM40", "MMD", "MYL9", "PTGS1", "MPP1", "ACTN1", "ACRBP", 
                         "TPM4", "TPM1", "SNCA", "NGFRAP1", "GAS2L1", "GRAP2", "CTSA", 
                         "ILK", "ODC1", "SNN", "PDLIM1", "FERMT3", "HIST1H2BK", "MARCH2", 
                         "PARVB", "PLA2G12A", "MFSD1", "PVALB", "CD151", "TPTEP1", "GSN", 
                         "YWHAH", "RAP1B", "ASAH1", "APP", "MYL12A", "H2AFJ", "TYMP", 
                         "COTL1", "CTSS", "S100A11", "CFD", "LGALS2", "GABARAP", "TNFSF13B", 
                         "APOBEC3A", "AP2S1", "IGSF6", "CSTB", "CPVL", "MAFB", "IFI6", 
                         "SLC7A7", "TGFBI", "C1orf162", "CAPG", "BLVRA", "BID", "MNDA", 
                         "LILRB2", "S100A10", "CD14", "WARS", "NAAA", "NUP214", "ATG3", 
                         "CNPY3", "CYBB", "GAPDH", "IFNGR2", "ALDH2", "EIF4EBP1", "CD300LF", 
                         "GCA", "GRINA", "LILRB4", "SOD2", "CTSZ", "CD302", "CAMK1", "POU2F2", 
                         "FCGR2A")
  expect_equal(
    c(res_100@top_genes$`2`, res_100@top_genes$`3`),
    gene_name_to_check
  )
  
  # Test 'fast' arg.
  b <- top_genes(pbmc3k_medium_clusters, fast=F)
  a <- top_genes(pbmc3k_medium_clusters, fast=T)
  
  expect_equal(all(unlist(mapply("==", a@top_genes, b@top_genes))), TRUE)
  
  
})


