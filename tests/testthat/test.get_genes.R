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
