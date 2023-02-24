test_that(paste("Checking get_data_for_scigenex stops if data provided is not",
                "a dataframe, a matrix or a Seurat object"), {
                  set_verbosity(0)
                  data_for_scigenex <- seq(1,20)
                  expect_error(get_data_for_scigenex(data_for_scigenex))
                }
)

test_that(paste("Checking get_data_for_scigenex stops if data provided NULL"), {
  set_verbosity(0)
  data_for_scigenex <- NULL
  expect_error(get_data_for_scigenex(data_for_scigenex))
})

test_that(paste("Checking matrix obtained with get_data_for_scigenex",
                "using a seurat object as data argument"), {
                  
                  set_verbosity(0)
                  data(pbmc_small, package = "SeuratObject")
                  data_for_scigenex <- pbmc_small
                  expr_matrix <- get_data_for_scigenex(data_for_scigenex)

                  expect_equal(class(expr_matrix)[1], "dgCMatrix")
                  expect_equal(round(mean(as.matrix(expr_matrix)), 6), 1.110883)
                  expect_equal(round(median(as.matrix(expr_matrix)), 6), 0)
                  expect_equal(round(sd(as.matrix(expr_matrix)), 6), 2.031546)
                  expect_equal(round(sum(as.matrix(expr_matrix^2)), 2), 98642.66)
                  expect_equal(length(as.matrix(expr_matrix)), 18400)

                  expect_equal(round(sum(quantile(as.matrix(expr_matrix))["0%"]), 6), 0)
                  expect_equal(round(sum(quantile(as.matrix(expr_matrix))["25%"]), 6), 0)
                  expect_equal(round(sum(quantile(as.matrix(expr_matrix))["50%"]), 6), 0)
                  expect_equal(round(sum(quantile(as.matrix(expr_matrix))["75%"]), 6), 0)
                  expect_equal(round(sum(quantile(as.matrix(expr_matrix))["100%"]), 6), 8.092033)

                  expect_equal(round(sum(colMeans(expr_matrix)),5), 88.87063)
                  expect_equal(round(sum(rowMeans(expr_matrix)),4), 255.5031)

                  expect_equal(ncol(expr_matrix), 80)
                  expect_equal(colnames(expr_matrix), c(
                    "ATGCCAGAACGACT", "CATGGCCTGTGCAT", "GAACCTGATGAACC", "TGACTGGATTCTCA",
                    "AGTCAGACTGCACA", "TCTGATACACGTGT", "TGGTATCTAAACAG", "GCAGCTCTGTTTCT",
                    "GATATAACACGCAT", "AATGTTGACAGTCA", "AGGTCATGAGTGTC", "AGAGATGATCTCGC",
                    "GGGTAACTCTAGTG", "CATGAGACACGGGA", "TACGCCACTCCGAA", "CTAAACCTGTGCAT",
                    "GTAAGCACTCATTC", "TTGGTACTGAATCC", "CATCATACGGAGCA", "TACATCACGCTAAC",
                    "TTACCATGAATCGC", "ATAGGAGAAACAGA", "GCGCACGACTTTAC", "ACTCGCACGAAAGT",
                    "ATTACCTGCCTTAT", "CCCAACTGCAATCG", "AAATTCGAATCACG", "CCATCCGATTCGCC",
                    "TCCACTCTGAGCTT", "CATCAGGATGCACA", "CTAAACCTCTGACA", "GATAGAGAAGGGTG",
                    "CTAACGGAACCGAT", "AGATATACCCGTAA", "TACTCTGAATCGAC", "GCGCATCTTGCTCC",
                    "GTTGACGATATCGG", "ACAGGTACTGGTGT", "GGCATATGCTTATC", "CATTACACCAACTG",
                    "TAGGGACTGAACTC", "GCTCCATGAGAAGT", "TACAATGATGCTAG", "CTTCATGACCGAAT",
                    "CTGCCAACAGGAGC", "TTGCATTGAGCTAC", "AAGCAAGAGCTTAG", "CGGCACGAACTCAG",
                    "GGTGGAGATTACTC", "GGCCGATGTACTCT", "CGTAGCCTGTATGC", "TGAGCTGAATGCTG",
                    "CCTATAACGAGACG", "ATAAGTTGGTACGT", "AAGCGACTTTGACG", "ACCAGTGAATACCG",
                    "ATTGCACTTGCTTT", "CTAGGTGATGGTTG", "GCACTAGACCTTTA", "CATGCGCTAGTCAC",
                    "TTGAGGACTACGCA", "ATACCACTCTAAGC", "CATATAGACTAAGC", "TTTAGCTGTACTCT",
                    "GACATTCTCCACCT", "ACGTGATGCCATGA", "ATTGTAGATTCCCG", "GATAGAGATCACGA",
                    "AATGCGTGGACGGA", "GCGTAAACACGGTT", "ATTCAGCTCATTGG", "GGCATATGGGGAGT",
                    "ATCATCTGACACCA", "GTCATACTTCGCCT", "TTACGTACGTTCAG", "GAGTTGTGGTAGCT",
                    "GACGCTCTCTCTCG", "AGTCTTACTTCGGA", "GGAACACTTCAGAC", "CTTGATTGATCTTC"
                  ))

                  expect_equal(nrow(expr_matrix), 230)
                  expect_equal(rownames(expr_matrix), c(
                    "MS4A1", "CD79B", "CD79A", "HLA-DRA", "TCL1A", "HLA-DQB1",
                    "HVCN1", "HLA-DMB", "LTB", "LINC00926", "FCER2", "SP100",
                    "NCF1", "PPP3CC", "EAF2", "PPAPDC1B", "CD19", "KIAA0125",
                    "CYB561A3", "CD180", "RP11-693J15.5", "FAM96A", "CXCR4", "STX10",
                    "SNHG7", "NT5C", "BANK1", "IGLL5", "CD200", "FCRLA",
                    "CD3D", "NOSIP", "SAFB2", "CD2", "IL7R", "PIK3IP1",
                    "MPHOSPH6", "KHDRBS1", "MAL", "CCR7", "THYN1", "TAF7",
                    "LDHB", "TMEM123", "CCDC104", "EPC1", "EIF4A2", "CD3E",
                    "TMUB1", "BLOC1S4", "ACSM3", "TMEM204", "SRSF7", "ACAP1",
                    "TNFAIP8", "CD7", "TAGAP", "DNAJB1", "ASNSD1", "S1PR4",
                    "CTSW", "GZMK", "NKG7", "IL32", "DNAJC2", "LYAR",
                    "CST7", "LCK", "CCL5", "HNRNPH1", "SSR2", "DLGAP1-AS1",
                    "GIMAP1", "MMADHC", "ZNF76", "CD8A", "PTPN22", "GYPC",
                    "HNRNPF", "RPL7L1", "KLRG1", "CRBN", "SATB1", "SIT1",
                    "PMPCB", "NRBP1", "TCF7", "HNRNPA3", "S100A8", "S100A9",
                    "LYZ", "CD14", "FCN1", "TYROBP", "ASGR1", "NFKBIA",
                    "TYMP", "CTSS", "TSPO", "RBP7", "CTSB", "LGALS1",
                    "FPR1", "VSTM1", "BLVRA", "MPEG1", "BID", "SMCO4",
                    "CFD", "LINC00936", "LGALS2", "MS4A6A", "FCGRT", "LGALS3",
                    "NUP214", "SCO2", "IL17RA", "IFI6", "HLA-DPA1", "FCER1A",
                    "CLEC10A", "HLA-DMA", "RGS1", "HLA-DPB1", "HLA-DQA1", "RNF130",
                    "HLA-DRB5", "HLA-DRB1", "CST3", "IL1B", "POP7", "HLA-DQA2",
                    "CD1C", "GSTP1", "EIF3G", "VPS28", "LY86", "ZFP36L1",
                    "ZNF330", "ANXA2", "GRN", "CFP", "HSP90AA1", "FUOM",
                    "LST1", "AIF1", "PSAP", "YWHAB", "MYO1G", "SAT1",
                    "RGS2", "SERPINA1", "IFITM3", "FCGR3A", "LILRA3", "S100A11",
                    "FCER1G", "TNFRSF1B", "IFITM2", "WARS", "IFI30", "MS4A7",
                    "C5AR1", "HCK", "COTL1", "LGALS9", "CD68", "RP11-290F20.3",
                    "RHOC", "CARD16", "LRRC25", "COPS6", "ADAR", "PPBP",
                    "GPX1", "TPM4", "PF4", "SDPR", "NRGN", "SPARC",
                    "GNG11", "CLU", "HIST1H2AC", "NCOA4", "GP9", "FERMT3",
                    "ODC1", "CD9", "RUFY1", "TUBB1", "TALDO1", "TREML1",
                    "NGFRAP1", "PGRMC1", "CA2", "ITGA2B", "MYL9", "TMEM40",
                    "PARVB", "PTCRA", "ACRBP", "TSC22D1", "VDAC3", "GZMB",
                    "GZMA", "GNLY", "FGFBP2", "AKR1C3", "CCL4", "PRF1",
                    "GZMH", "XBP1", "GZMM", "PTGDR", "IGFBP7", "TTC38",
                    "KLRD1", "ARHGDIA", "IL2RB", "CLIC3", "PPP1R18", "CD247",
                    "ALOX5AP", "XCL2", "C12orf75", "RARRES3", "PCMT1", "LAMP1",
                    "SPON2", "S100B"
                  ))

                  expect_that(expr_matrix, is_a("dgCMatrix"))
                })


test_that(paste("Checking matrix obtained with get_data_for_scigenex",
                "using a dataframe as data argument"), {
                  set_verbosity(0)
                  data_for_scigenex <- as.data.frame(create_4_rnd_clust())
                  expr_matrix <- get_data_for_scigenex(data_for_scigenex)
                  
                  expect_equal(round(mean(expr_matrix), 6), 0.029886)
                  expect_equal(round(median(expr_matrix), 6), 0.00703)
                  expect_equal(round(sd(expr_matrix), 6), 1.169623)
                  expect_equal(round(sum(expr_matrix^2), 1), 109511.6)
                  expect_equal(length(expr_matrix), 80000)
                  
                  expect_equal(round(sum(quantile(expr_matrix)["0%"]), 6), -5.806817)
                  expect_equal(round(sum(quantile(expr_matrix)["25%"]), 6), -0.695268)
                  expect_equal(round(sum(quantile(expr_matrix)["50%"]), 6), 0.00703)
                  expect_equal(round(sum(quantile(expr_matrix)["75%"]), 6), 0.718151)
                  expect_equal(round(sum(quantile(expr_matrix)["100%"]), 6), 7.073935)
                  
                  expect_equal(round(sum(colMeans(expr_matrix)),5), 0.59772)
                  expect_equal(round(sum(rowMeans(expr_matrix)),4), 119.544)
                  
                  expect_equal(ncol(expr_matrix), 20)
                  expect_equal(colnames(expr_matrix), paste0("V", seq(1,20)))
                  
                  expect_equal(nrow(expr_matrix), 4000) 
                  
                  expect_that(expr_matrix, is_a("matrix"))
                }
)


test_that(paste("Checking matrix obtained with get_data_for_scigenex",
                "using a matrix as data argument"), {
                  
                  set_verbosity(0)

                  data_for_scigenex <- create_4_rnd_clust()
                  expr_matrix <- get_data_for_scigenex(data_for_scigenex)
                  
                  expect_equal(round(mean(expr_matrix), 6), 0.029886)
                  expect_equal(round(median(expr_matrix), 6), 0.00703)
                  expect_equal(round(sd(expr_matrix), 6), 1.169623)
                  expect_equal(round(sum(expr_matrix^2), 1), 109511.6)
                  expect_equal(length(expr_matrix), 80000)
                  
                  expect_equal(round(sum(quantile(expr_matrix)["0%"]), 6), -5.806817)
                  expect_equal(round(sum(quantile(expr_matrix)["25%"]), 6), -0.695268)
                  expect_equal(round(sum(quantile(expr_matrix)["50%"]), 6), 0.00703)
                  expect_equal(round(sum(quantile(expr_matrix)["75%"]), 6), 0.718151)
                  expect_equal(round(sum(quantile(expr_matrix)["100%"]), 6), 7.073935)
                  
                  expect_equal(round(sum(colMeans(expr_matrix)),5), 0.59772)
                  expect_equal(round(sum(rowMeans(expr_matrix)),4), 119.544)
                  
                  expect_equal(ncol(expr_matrix), 20)
                  expect_equal(colnames(expr_matrix), paste0("sample", seq(1,20)))
                  
                  expect_equal(nrow(expr_matrix), 4000) 
                  
                  expect_that(expr_matrix, is_a("matrix"))
                }
)

test_that(paste("Checking matrix obtained with get_data_for_scigenex",
                "using a data.frame as data argument"), {
                  
                  set_verbosity(0)

                  data_for_scigenex <- as.data.frame(create_4_rnd_clust())
                  expr_matrix <- get_data_for_scigenex(data_for_scigenex)
                  
                  expect_equal(round(mean(expr_matrix), 6), 0.029886)
                  expect_equal(round(median(expr_matrix), 6), 0.00703)
                  expect_equal(round(sd(expr_matrix), 6), 1.169623)
                  expect_equal(round(sum(expr_matrix^2), 1), 109511.6)
                  expect_equal(length(expr_matrix), 80000)
                  
                  expect_equal(round(sum(quantile(expr_matrix)["0%"]), 6), -5.806817)
                  expect_equal(round(sum(quantile(expr_matrix)["25%"]), 6), -0.695268)
                  expect_equal(round(sum(quantile(expr_matrix)["50%"]), 6), 0.00703)
                  expect_equal(round(sum(quantile(expr_matrix)["75%"]), 6), 0.718151)
                  expect_equal(round(sum(quantile(expr_matrix)["100%"]), 6), 7.073935)
                  
                  expect_equal(round(sum(colMeans(expr_matrix)),5), 0.59772)
                  expect_equal(round(sum(rowMeans(expr_matrix)),4), 119.544)
                  
                  expect_equal(ncol(expr_matrix), 20)
                  expect_equal(colnames(expr_matrix), paste0("V", seq(1,20)))
                  
                  expect_equal(nrow(expr_matrix), 4000) 
                  
                  expect_that(expr_matrix, is_a("matrix"))
                }
)


test_that(paste("Checking SCT with pmbc_samll"), {
                  
                  set_verbosity(0)
                  data(pbmc_small, package = "SeuratObject")
                  data_for_scigenex <- pbmc_small
                  
                  expect_error(get_data_for_scigenex(data_for_scigenex, which_slot ="sct"))
                  
                  data_for_scigenex <- suppressWarnings(Seurat::SCTransform(pbmc_small, verbose = FALSE))
                  expr_matrix <- get_data_for_scigenex(data_for_scigenex, which_slot ="sct")
                  
                  expect_equal(round(mean(as.matrix(expr_matrix[!is.na(expr_matrix)])), 3), 0.282)
                  expect_equal(ncol(expr_matrix), 80)
                  expect_equal(nrow(expr_matrix), 220) 
                  
                  expect_that(expr_matrix, is_a("dgCMatrix"))
                }
)
  