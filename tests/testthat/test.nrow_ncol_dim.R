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
  expect_that(col_names(res), is_a("character"))
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
  expect_that(row_names(res), is_a("character"))
})

test_that("Checking clust_size()", {
  expect_equal(clust_size(res), c(`1` = 51L, `2` = 49L, `3` = 45L, `4` = 24L, `5` = 20L, `6` = 18L, 
                                  `7` = 14L, `8` = 14L, `9` = 14L, `10` = 12L, `11` = 7L, `12` = 7L, 
                                  `13` = 6L, `14` = 5L, `15` = 5L))
  expect_that(clust_size(res), is_a("integer"))
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
  expect_that(gene_cluster(res), is_a("integer"))
})

