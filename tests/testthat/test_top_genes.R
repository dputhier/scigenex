test_that("Cheking get_genes is providing the right list of genes", {
  #Create matrix containing 3 signatures
  set.seed(123)
  m <- matrix(rnorm(40000), nc=20)
  m[1:100,1:10] <- m[1:100,1:10] + 4
  m[101:200,11:20] <- m[101:200,11:20] + 3
  m[201:300,5:15] <- m[201:300,5:15] + -2
  res <- find_gene_clusters(data=m,
                            name = "test",
                            distance_method="pearson",
                            av_dot_prod_min = 0,
                            inflation = 2,
                            k=25,
                            fdr = 10)
  
  genes_top <- top_genes(res, cluster = "all", top = 20)
  #Test 20 top genes list in all cluster
  gene_name_to_check <- matrix(c("gene58", "gene37", "gene81", "gene27", "gene94", "gene54", "gene4", "gene65", "gene11", "gene35", "gene84", "gene14", 
                                 "gene28", "gene82", "gene20", "gene23", "gene48", "gene55", "gene12", "gene79", "gene166", "gene186", "gene192", "gene183", 
                                 "gene117", "gene155", "gene180", "gene114", "gene122", "gene163", "gene121", "gene150", "gene168", "gene181", "gene200", 
                                 "gene171", "gene104", "gene113", "gene165", "gene189", "gene201", "gene279", "gene212", "gene289", "gene266", "gene249", 
                                 "gene293", "gene272", "gene242", "gene277", "gene213", "gene210", "gene290", "gene233", "gene273", "gene278", "gene215", 
                                 "gene228", "gene248", "gene298"), nrow = 3, ncol = 20, byrow = T)
  colnames(gene_name_to_check) <- paste0("gene_top_", 1:20)
  rownames(gene_name_to_check) <- paste0("cluster_", 1:3)
  expect_equal(genes_top@top_genes, gene_name_to_check)
  
  genes_top <- top_genes(res, cluster = 1, top = 20)
  #Test 20 top genes list in cluster 1
  gene_name_to_check <- matrix(c("gene58", "gene37", "gene81", "gene27", "gene94", "gene54", "gene4", "gene65", "gene11", "gene35", "gene84", "gene14", 
                          "gene28", "gene82", "gene20", "gene23", "gene48", "gene55", "gene12", "gene79"), nrow = 1, ncol = 20, byrow = T)
  colnames(gene_name_to_check) <- paste0("gene_top_", 1:20)
  rownames(gene_name_to_check) <- paste0("cluster_", 1)
  expect_equal(genes_top@top_genes, gene_name_to_check)
  
  genes_top <- top_genes(res, cluster = 2:3, top = 20)
  #Test 20 top genes list in cluster 2 and 3
  gene_name_to_check <- matrix(c("gene166", "gene186", "gene192", "gene183", "gene117", "gene155", "gene180", "gene114", "gene122", "gene163", "gene121", 
                          "gene150", "gene168", "gene181", "gene200", "gene171", "gene104", "gene113", "gene165", "gene189", "gene201", "gene279", 
                          "gene212", "gene289", "gene266", "gene249", "gene293", "gene272", "gene242", "gene277", "gene213", "gene210", "gene290", 
                          "gene233", "gene273", "gene278", "gene215", "gene228", "gene248", "gene298"), nrow = 2, ncol = 20, byrow = T)
  colnames(gene_name_to_check) <- paste0("gene_top_", 1:20)
  rownames(gene_name_to_check) <- paste0("cluster_", 2:3)
  expect_equal(genes_top@top_genes, gene_name_to_check)
  
  genes_top <- top_genes(res, cluster = "all", top = 10)
  #Test 20 top genes list in all cluster
  gene_name_to_check <- matrix(c("gene58", "gene37", "gene81", "gene27", "gene94", "gene54", "gene4", "gene65", "gene11", "gene35",
                                 "gene166", "gene186", "gene192", "gene183", "gene117", "gene155", "gene180", "gene114", "gene122", "gene163",
                                 "gene201", "gene279", "gene212", "gene289", "gene266", "gene249", "gene293", "gene272", "gene242", "gene277"), nrow = 3, ncol = 10, byrow = T)
  colnames(gene_name_to_check) <- paste0("gene_top_", 1:10)
  rownames(gene_name_to_check) <- paste0("cluster_", 1:3)
  expect_equal(genes_top@top_genes, gene_name_to_check)
  
  genes_top <- top_genes(res, cluster = 1, top = 10)
  #Test 20 top genes list in cluster 1
  gene_name_to_check <- matrix(c("gene58", "gene37", "gene81", "gene27", "gene94", "gene54", "gene4", "gene65", "gene11", "gene35"), nrow = 1, ncol = 10, byrow = T)
  colnames(gene_name_to_check) <- paste0("gene_top_", 1:10)
  rownames(gene_name_to_check) <- paste0("cluster_", 1)
  expect_equal(genes_top@top_genes, gene_name_to_check)
  
  genes_top <- top_genes(res, cluster = 2:3, top = 10)
  #Test 20 top genes list in cluster 2 and 3
  gene_name_to_check <- matrix(c("gene166", "gene186", "gene192", "gene183", "gene117", "gene155", "gene180", "gene114", "gene122", "gene163",
                                 "gene201", "gene279", "gene212", "gene289", "gene266", "gene249", "gene293", "gene272", "gene242", "gene277"), nrow = 2, ncol = 10, byrow = T)
  colnames(gene_name_to_check) <- paste0("gene_top_", 1:10)
  rownames(gene_name_to_check) <- paste0("cluster_", 2:3)
  expect_equal(genes_top@top_genes, gene_name_to_check)
  
  genes_top <- top_genes(res, cluster = "all", top = 5)
  #Test 20 top genes list in all cluster
  gene_name_to_check <- matrix(c("gene58", "gene37", "gene81", "gene27", "gene94",
                                 "gene166", "gene186", "gene192", "gene183", "gene117", 
                                 "gene201", "gene279", "gene212", "gene289", "gene266"), nrow = 3, ncol = 5, byrow = T)
  colnames(gene_name_to_check) <- paste0("gene_top_", 1:5)
  rownames(gene_name_to_check) <- paste0("cluster_", 1:3)
  expect_equal(genes_top@top_genes, gene_name_to_check)
  
  genes_top <- top_genes(res, cluster = 1, top = 5)
  #Test 20 top genes list in cluster 1
  gene_name_to_check <- matrix(c("gene58", "gene37", "gene81", "gene27", "gene94"), nrow = 1, ncol = 5, byrow = T)
  colnames(gene_name_to_check) <- paste0("gene_top_", 1:5)
  rownames(gene_name_to_check) <- paste0("cluster_", 1)
  expect_equal(genes_top@top_genes, gene_name_to_check)
  
  genes_top <- top_genes(res, cluster = 2:3, top = 5)
  #Test 20 top genes list in cluster 2 and 3
  gene_name_to_check <- matrix(c("gene166", "gene186", "gene192", "gene183", "gene117",
                                 "gene201", "gene279", "gene212", "gene289", "gene266"), nrow = 2, ncol = 5, byrow = T)
  colnames(gene_name_to_check) <- paste0("gene_top_", 1:5)
  rownames(gene_name_to_check) <- paste0("cluster_", 2:3)
  expect_equal(genes_top@top_genes, gene_name_to_check)
  
  #Remove output files
  file.remove("test.dbf_out.txt")
  file.remove("test.mcl_out.txt")
})
