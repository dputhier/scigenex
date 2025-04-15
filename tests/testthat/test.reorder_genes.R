# Set verbosity to 0
set_verbosity(0)

load_example_dataset("7871581/files/pbmc3k_medium_clusters")
res <- pbmc3k_medium_clusters

test_that("Checking if reorder_genes stops if needed", {
  expect_error(reorder_genes(res, order_by = "not_ok"))
})


test_that("Checking reorder_genes function using order_by='gene_names'...", {
  res <- reorder_genes(res, order_by = "gene_names")
  expect_equal(res@gene_clusters$`1`, c("BTG1", "EEF1A1", "EEF1B2", "MALAT1", "NACA", "RPL11", "RPL12", 
                                        "RPL14", "RPL15", "RPL17", "RPL18", "RPL18A", "RPL21", "RPL23A", 
                                        "RPL26", "RPL27", "RPL27A", "RPL28", "RPL29", "RPL3", "RPL31", 
                                        "RPL32", "RPL35A", "RPL36", "RPL5", "RPL6", "RPL7A", "RPL9", 
                                        "RPLP1", "RPLP2", "RPS10", "RPS11", "RPS12", "RPS13", "RPS15A", 
                                        "RPS16", "RPS18", "RPS19", "RPS20", "RPS23", "RPS25", "RPS26", 
                                        "RPS28", "RPS29", "RPS3", "RPS3A", "RPS4X", "RPS8", "RPS9", "RPSAP58", 
                                        "TPT1"))
  
  expect_equal(res@gene_clusters$`2`, c("ACRBP", "ACTN1", "AP001189.4", "APP", "ASAH1", "CA2", "CD151", 
                                        "CLU", "CTSA", "F13A1", "FERMT3", "GAS2L1", "GNG11", "GP9", "GRAP2", 
                                        "GSN", "H2AFJ", "HIST1H2AC", "HIST1H2BK", "ILK", "MARCH2", "MFSD1", 
                                        "MMD", "MPP1", "MYL12A", "MYL9", "NGFRAP1", "ODC1", "PARVB", 
                                        "PDLIM1", "PF4", "PLA2G12A", "PPBP", "PTCRA", "PTGS1", "PVALB", 
                                        "RAP1B", "RGS18", "SDPR", "SNCA", "SNN", "SPARC", "TMEM40", "TPM1", 
                                        "TPM4", "TPTEP1", "TREML1", "TUBB1", "YWHAH"))
  
  expect_equal(res@gene_clusters$`3`, c("ALDH2", "AP2S1", "APOBEC3A", "ATG3", "BID", "BLVRA", "C1orf162", 
                                        "CAMK1", "CAPG", "CD14", "CD300LF", "CD302", "CFD", "CNPY3", 
                                        "COTL1", "CPVL", "CSTB", "CTSS", "CTSZ", "CYBB", "EIF4EBP1", 
                                        "FCGR2A", "GABARAP", "GAPDH", "GCA", "GRINA", "IFI6", "IFNGR2", 
                                        "IGSF6", "LGALS2", "LILRB2", "LILRB4", "MAFB", "MNDA", "NAAA", 
                                        "NUP214", "POU2F2", "S100A10", "S100A11", "SLC7A7", "SOD2", "TGFBI", 
                                        "TNFSF13B", "TYMP", "WARS"))
  
  expect_equal(res@gene_clusters$`4`, c("AKR1C3", "CCL4", "CCL5", "CD247", "CD7", "CD99", "CLIC3", 
                                        "CST7", "CTSW", "EFHD2", "FGFBP2", "GNLY", "GZMA", "GZMB", "HLA-C", 
                                        "HOPX", "IGFBP7", "IL2RB", "KLRF1", "PLEK", "PRF1", "PTPRCAP", 
                                        "SPON2", "XCL1"))
  
  res <- reorder_genes(res, order_by = "gene_names", decreasing=TRUE)
  expect_equal(res@gene_clusters$`1`, rev(c("BTG1", "EEF1A1", "EEF1B2", "MALAT1", "NACA", "RPL11", "RPL12", 
                                        "RPL14", "RPL15", "RPL17", "RPL18", "RPL18A", "RPL21", "RPL23A", 
                                        "RPL26", "RPL27", "RPL27A", "RPL28", "RPL29", "RPL3", "RPL31", 
                                        "RPL32", "RPL35A", "RPL36", "RPL5", "RPL6", "RPL7A", "RPL9", 
                                        "RPLP1", "RPLP2", "RPS10", "RPS11", "RPS12", "RPS13", "RPS15A", 
                                        "RPS16", "RPS18", "RPS19", "RPS20", "RPS23", "RPS25", "RPS26", 
                                        "RPS28", "RPS29", "RPS3", "RPS3A", "RPS4X", "RPS8", "RPS9", "RPSAP58", 
                                        "TPT1")))
  
  expect_equal(res@gene_clusters$`2`, rev(c("ACRBP", "ACTN1", "AP001189.4", "APP", "ASAH1", "CA2", "CD151", 
                                        "CLU", "CTSA", "F13A1", "FERMT3", "GAS2L1", "GNG11", "GP9", "GRAP2", 
                                        "GSN", "H2AFJ", "HIST1H2AC", "HIST1H2BK", "ILK", "MARCH2", "MFSD1", 
                                        "MMD", "MPP1", "MYL12A", "MYL9", "NGFRAP1", "ODC1", "PARVB", 
                                        "PDLIM1", "PF4", "PLA2G12A", "PPBP", "PTCRA", "PTGS1", "PVALB", 
                                        "RAP1B", "RGS18", "SDPR", "SNCA", "SNN", "SPARC", "TMEM40", "TPM1", 
                                        "TPM4", "TPTEP1", "TREML1", "TUBB1", "YWHAH")))
  
  expect_equal(res@gene_clusters$`3`, rev(c("ALDH2", "AP2S1", "APOBEC3A", "ATG3", "BID", "BLVRA", "C1orf162", 
                                        "CAMK1", "CAPG", "CD14", "CD300LF", "CD302", "CFD", "CNPY3", 
                                        "COTL1", "CPVL", "CSTB", "CTSS", "CTSZ", "CYBB", "EIF4EBP1", 
                                        "FCGR2A", "GABARAP", "GAPDH", "GCA", "GRINA", "IFI6", "IFNGR2", 
                                        "IGSF6", "LGALS2", "LILRB2", "LILRB4", "MAFB", "MNDA", "NAAA", 
                                        "NUP214", "POU2F2", "S100A10", "S100A11", "SLC7A7", "SOD2", "TGFBI", 
                                        "TNFSF13B", "TYMP", "WARS")))
  
  expect_equal(res@gene_clusters$`4`, rev(c("AKR1C3", "CCL4", "CCL5", "CD247", "CD7", "CD99", "CLIC3", 
                                        "CST7", "CTSW", "EFHD2", "FGFBP2", "GNLY", "GZMA", "GZMB", "HLA-C", 
                                        "HOPX", "IGFBP7", "IL2RB", "KLRF1", "PLEK", "PRF1", "PTPRCAP", 
                                        "SPON2", "XCL1")))
})



test_that("Checking reorder_genes function using order_by='hclust'...", {
  set.seed(123)
  res <- reorder_genes(res, order_by = "hclust")
  
  expect_equal(res@gene_clusters$`1`, c("RPSAP58", "BTG1", "EEF1B2", "RPS26", "RPS29", "RPS11", "RPS20", 
                                        "RPL36", "NACA", "RPL5", "RPS28", "RPS13", "RPL17", "RPL27", 
                                        "RPS10", "RPS3A", "RPL31", "RPS9", "TPT1", "RPL7A", "RPL35A", 
                                        "RPS25", "RPL29", "RPL18", "RPL9", "MALAT1", "RPL21", "RPL27A", 
                                        "RPL3", "RPL14", "RPS16", "RPLP2", "RPL26", "RPLP1", "RPS19", 
                                        "RPL12", "RPL28", "RPS12", "RPL18A", "RPS18", "RPL11", "RPL32", 
                                        "RPS4X", "RPS3", "EEF1A1", "RPS8", "RPL15", "RPS15A", "RPL6", 
                                        "RPL23A", "RPS23"))
  
  expect_equal(res@gene_clusters$`2`, c("APP", "H2AFJ", "YWHAH", "ASAH1", "FERMT3", "RAP1B", "MYL12A", 
                                        "TPTEP1", "CD151", "PVALB", "GSN", "MFSD1", "PLA2G12A", "MARCH2", 
                                        "PARVB", "HIST1H2BK", "PDLIM1", "SNN", "ODC1", "ILK", "CTSA", 
                                        "GAS2L1", "SNCA", "GRAP2", "TPM4", "TPM1", "MYL9", "MPP1", "PTGS1", 
                                        "NGFRAP1", "ACTN1", "MMD", "TMEM40", "PTCRA", "ACRBP", "F13A1", 
                                        "CA2", "TREML1", "RGS18", "CLU", "GP9", "AP001189.4", "TUBB1", 
                                        "HIST1H2AC", "PF4", "PPBP", "SPARC", "SDPR", "GNG11"))
  
  expect_equal(res@gene_clusters$`3`, c("CAMK1", "FCGR2A", "CD300LF", "CD302", "GRINA", "POU2F2", "LILRB2", 
                                        "LILRB4", "CYBB", "CTSZ", "EIF4EBP1", "GCA", "CD14", "ALDH2", 
                                        "WARS", "NUP214", "SOD2", "SLC7A7", "IFNGR2", "CPVL", "MNDA", 
                                        "NAAA", "TGFBI", "MAFB", "IGSF6", "ATG3", "APOBEC3A", "BLVRA", 
                                        "BID", "CAPG", "CNPY3", "TNFSF13B", "LGALS2", "CFD", "IFI6", 
                                        "C1orf162", "AP2S1", "CSTB", "S100A10", "COTL1", "GAPDH", "GABARAP", 
                                        "TYMP", "CTSS", "S100A11"))
  
  expect_equal(res@gene_clusters$`4`, c("EFHD2", "PLEK", "AKR1C3", "IGFBP7", "KLRF1", "XCL1", "CCL5", 
                                        "CD99", "HLA-C", "PTPRCAP", "IL2RB", "CD247", "CD7", "CLIC3", 
                                        "HOPX", "SPON2", "FGFBP2", "PRF1", "GZMB", "GNLY", "CCL4", "CTSW", 
                                        "GZMA", "CST7"))
})




test_that("Checking reorder_genes function using order_by='correlation'...", {
  res <- reorder_genes(res, order_by = "correlation")
  
  expect_equal(res@gene_clusters$`1`, c("RPL11", "RPL32", "RPS12", "RPL18A", "RPS18", "RPS23", "RPS16", 
                                        "RPL9", "RPL23A", "RPS3", "RPL12", "RPS15A", "RPLP2", "RPL28", 
                                        "RPS19", "RPS25", "EEF1A1", "RPS8", "RPS4X", "RPL29", "RPL27A", 
                                        "RPL18", "RPL14", "RPL3", "RPL35A", "RPL21", "RPLP1", "RPS28", 
                                        "RPL6", "RPS10", "RPL26", "RPL7A", "RPL5", "RPL27", "RPL15", 
                                        "RPS13", "RPL17", "RPL31", "RPS3A", "RPS29", "RPL36", "MALAT1", 
                                        "RPS20", "TPT1", "RPS26", "RPS9", "EEF1B2", "RPS11", "NACA", 
                                        "RPSAP58", "BTG1"))
  
  expect_equal(res@gene_clusters$`2`, c("SDPR", "GNG11", "PF4", "SPARC", "PPBP", "TUBB1", "HIST1H2AC", 
                                        "GP9", "CLU", "AP001189.4", "RGS18", "PTCRA", "TREML1", "F13A1", 
                                        "CA2", "TMEM40", "MMD", "MYL9", "PTGS1", "MPP1", "ACTN1", "ACRBP", 
                                        "TPM4", "TPM1", "SNCA", "NGFRAP1", "GAS2L1", "GRAP2", "CTSA", 
                                        "ILK", "ODC1", "SNN", "PDLIM1", "FERMT3", "HIST1H2BK", "MARCH2", 
                                        "PARVB", "PLA2G12A", "MFSD1", "PVALB", "CD151", "TPTEP1", "GSN", 
                                        "YWHAH", "RAP1B", "ASAH1", "APP", "MYL12A", "H2AFJ"))
  
  expect_equal(res@gene_clusters$`3`, c("TYMP", "COTL1", "CTSS", "S100A11", "CFD", "LGALS2", "GABARAP", 
                                        "TNFSF13B", "APOBEC3A", "AP2S1", "IGSF6", "CSTB", "CPVL", "MAFB", 
                                        "IFI6", "SLC7A7", "TGFBI", "C1orf162", "CAPG", "BLVRA", "BID", 
                                        "MNDA", "LILRB2", "S100A10", "CD14", "WARS", "NAAA", "NUP214", 
                                        "ATG3", "CNPY3", "CYBB", "GAPDH", "IFNGR2", "ALDH2", "EIF4EBP1", 
                                        "CD300LF", "GCA", "GRINA", "LILRB4", "SOD2", "CTSZ", "CD302", 
                                        "CAMK1", "POU2F2", "FCGR2A"))
  
  expect_equal(res@gene_clusters$`4`, c("GZMB", "PRF1", "GNLY", "FGFBP2", "GZMA", "CST7", "CCL4", "SPON2", 
                                        "CTSW", "CLIC3", "CD247", "IL2RB", "HOPX", "HLA-C", "CD7", "AKR1C3", 
                                        "KLRF1", "XCL1", "IGFBP7", "CCL5", "CD99", "PTPRCAP", "PLEK", 
                                        "EFHD2"))
})

