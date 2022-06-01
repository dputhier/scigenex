test_that("Cheking get_genes is providing the right list of genes", {
  #Create matrix containing 3 signatures
  set.seed(123)
  m <- matrix(rnorm(40000), nc=20)
  m[1:100,1:10] <- m[1:100,1:10] + 4
  m[101:200,11:20] <- m[101:200,11:20] + 3
  m[201:300,5:15] <- m[201:300,5:15] + -2
  res <- DBFMCL(data=m,
                distance_method="pearson",
                av_dot_prod_min = 0,
                inflation = 2,
                k=25,
                fdr = 10)
  
  gene_names <- get_genes(res, cluster = "all")
  #Test gene list in all cluster
  gene_name_to_check <- c("gene94",   "gene58",   "gene37",   "gene81",   "gene4",    "gene54",   "gene27",   "gene84",   "gene20",   "gene14",   "gene17",   "gene39",   "gene67",   "gene48",   "gene16",   "gene98",  
                          "gene33",   "gene13",   "gene59",   "gene24",   "gene6",    "gene92",   "gene91",   "gene9",    "gene56",   "gene51",   "gene11",   "gene65",   "gene82",   "gene89",   "gene35",   "gene74",  
                          "gene50",   "gene79",   "gene41",   "gene88",   "gene8",    "gene23",   "gene44",   "gene100",  "gene75",   "gene2",    "gene15",   "gene47",   "gene40",   "gene28",   "gene12",   "gene32",  
                          "gene60",   "gene90",   "gene49",   "gene30",   "gene95",   "gene25",   "gene99",   "gene55",   "gene5",    "gene85",   "gene29",   "gene36",   "gene80",   "gene61",   "gene93",   "gene70",  
                          "gene73",   "gene31",   "gene19",   "gene7",    "gene76",   "gene42",   "gene34",   "gene46",   "gene64",   "gene63",   "gene69",   "gene78",   "gene71",   "gene96",   "gene62",   "gene52",  
                          "gene1",    "gene45",   "gene3",    "gene21",   "gene77",   "gene18",   "gene10",   "gene26",   "gene22",   "gene83",   "gene38",   "gene68",   "gene57",   "gene43",   "gene66",   "gene97",  
                          "gene53",   "gene86",   "gene87",   "gene72",   "gene740",  "gene637",  "gene186",  "gene166",  "gene183",  "gene180",  "gene192",  "gene155",  "gene117",  "gene114",  "gene122",  "gene163", 
                          "gene171",  "gene121",  "gene150",  "gene190",  "gene200",  "gene135",  "gene104",  "gene138",  "gene167",  "gene151",  "gene193",  "gene143",  "gene132",  "gene110",  "gene182",  "gene136", 
                          "gene168",  "gene181",  "gene113",  "gene125",  "gene133",  "gene106",  "gene172",  "gene184",  "gene175",  "gene116",  "gene173",  "gene170",  "gene187",  "gene152",  "gene160",  "gene107", 
                          "gene142",  "gene120",  "gene123",  "gene158",  "gene197",  "gene162",  "gene165",  "gene129",  "gene144",  "gene146",  "gene103",  "gene105",  "gene126",  "gene189",  "gene147",  "gene109", 
                          "gene195",  "gene159",  "gene131",  "gene185",  "gene149",  "gene174",  "gene139",  "gene119",  "gene140",  "gene124",  "gene157",  "gene115",  "gene178",  "gene127",  "gene176",  "gene179", 
                          "gene101",  "gene128",  "gene130",  "gene134",  "gene194",  "gene118",  "gene198",  "gene199",  "gene102",  "gene169",  "gene145",  "gene108",  "gene191",  "gene111",  "gene137",  "gene148", 
                          "gene153",  "gene164",  "gene112",  "gene188",  "gene156",  "gene141",  "gene1566", "gene154",  "gene177",  "gene240",  "gene161",  "gene196",  "gene279",  "gene201",  "gene212",  "gene249", 
                          "gene289",  "gene266",  "gene293",  "gene242",  "gene228",  "gene278",  "gene213",  "gene272",  "gene273",  "gene215",  "gene214",  "gene239",  "gene243",  "gene223",  "gene271",  "gene265", 
                          "gene205",  "gene220",  "gene222",  "gene261",  "gene297",  "gene227",  "gene210",  "gene277",  "gene263",  "gene219",  "gene238",  "gene252",  "gene211",  "gene233",  "gene300",  "gene248", 
                          "gene267",  "gene250",  "gene218",  "gene232",  "gene291",  "gene260",  "gene206",  "gene298",  "gene244",  "gene251",  "gene259",  "gene236",  "gene290",  "gene258",  "gene257",  "gene245", 
                          "gene269",  "gene276",  "gene296",  "gene246",  "gene280",  "gene292",  "gene295",  "gene230",  "gene216",  "gene275",  "gene207",  "gene231",  "gene235",  "gene466",  "gene256",  "gene255", 
                          "gene288",  "gene285",  "gene286",  "gene283",  "gene241",  "gene299",  "gene282")
  expect_equal(gene_names, gene_name_to_check)
  
  gene_names <- get_genes(res, cluster = 1)
  #Test gene list in cluster 1
  gene_name_to_check <- c("gene94",   "gene58",   "gene37",   "gene81",   "gene4",    "gene54",   "gene27",   "gene84",   "gene20",   "gene14",   "gene17",   "gene39",   "gene67",   "gene48",   "gene16",   "gene98",  
                          "gene33",   "gene13",   "gene59",   "gene24",   "gene6",    "gene92",   "gene91",   "gene9",    "gene56",   "gene51",   "gene11",   "gene65",   "gene82",   "gene89",   "gene35",   "gene74",  
                          "gene50",   "gene79",   "gene41",   "gene88",   "gene8",    "gene23",   "gene44",   "gene100",  "gene75",   "gene2",    "gene15",   "gene47",   "gene40",   "gene28",   "gene12",   "gene32",  
                          "gene60",   "gene90",   "gene49",   "gene30",   "gene95",   "gene25",   "gene99",   "gene55",   "gene5",    "gene85",   "gene29",   "gene36",   "gene80",   "gene61",   "gene93",   "gene70",  
                          "gene73",   "gene31",   "gene19",   "gene7",    "gene76",   "gene42",   "gene34",   "gene46",   "gene64",   "gene63",   "gene69",   "gene78",   "gene71",   "gene96",   "gene62",   "gene52",  
                          "gene1",    "gene45",   "gene3",    "gene21",   "gene77",   "gene18",   "gene10",   "gene26",   "gene22",   "gene83",   "gene38",   "gene68",   "gene57",   "gene43",   "gene66",   "gene97",  
                          "gene53",   "gene86",   "gene87",   "gene72",   "gene740",  "gene637")
  expect_equal(gene_names, gene_name_to_check)
  
  gene_names <- get_genes(res, cluster = 2)
  #Test gene list in cluster 2
  gene_name_to_check <- c("gene186",  "gene166",  "gene183",  "gene180",  "gene192",  "gene155",  "gene117",  "gene114",  "gene122",  "gene163", 
                          "gene171",  "gene121",  "gene150",  "gene190",  "gene200",  "gene135",  "gene104",  "gene138",  "gene167",  "gene151",  "gene193",  "gene143",  "gene132",  "gene110",  "gene182",  "gene136", 
                          "gene168",  "gene181",  "gene113",  "gene125",  "gene133",  "gene106",  "gene172",  "gene184",  "gene175",  "gene116",  "gene173",  "gene170",  "gene187",  "gene152",  "gene160",  "gene107", 
                          "gene142",  "gene120",  "gene123",  "gene158",  "gene197",  "gene162",  "gene165",  "gene129",  "gene144",  "gene146",  "gene103",  "gene105",  "gene126",  "gene189",  "gene147",  "gene109", 
                          "gene195",  "gene159",  "gene131",  "gene185",  "gene149",  "gene174",  "gene139",  "gene119",  "gene140",  "gene124",  "gene157",  "gene115",  "gene178",  "gene127",  "gene176",  "gene179", 
                          "gene101",  "gene128",  "gene130",  "gene134",  "gene194",  "gene118",  "gene198",  "gene199",  "gene102",  "gene169",  "gene145",  "gene108",  "gene191",  "gene111",  "gene137",  "gene148", 
                          "gene153",  "gene164",  "gene112",  "gene188",  "gene156",  "gene141",  "gene1566", "gene154",  "gene177",  "gene240",  "gene161",  "gene196")
  expect_equal(gene_names, gene_name_to_check)
  
  gene_names <- get_genes(res, cluster = 3)
  #Test gene list in cluster 3
  gene_name_to_check <- c("gene279",  "gene201",  "gene212",  "gene249", 
                          "gene289",  "gene266",  "gene293",  "gene242",  "gene228",  "gene278",  "gene213",  "gene272",  "gene273",  "gene215",  "gene214",  "gene239",  "gene243",  "gene223",  "gene271",  "gene265", 
                          "gene205",  "gene220",  "gene222",  "gene261",  "gene297",  "gene227",  "gene210",  "gene277",  "gene263",  "gene219",  "gene238",  "gene252",  "gene211",  "gene233",  "gene300",  "gene248", 
                          "gene267",  "gene250",  "gene218",  "gene232",  "gene291",  "gene260",  "gene206",  "gene298",  "gene244",  "gene251",  "gene259",  "gene236",  "gene290",  "gene258",  "gene257",  "gene245", 
                          "gene269",  "gene276",  "gene296",  "gene246",  "gene280",  "gene292",  "gene295",  "gene230",  "gene216",  "gene275",  "gene207",  "gene231",  "gene235",  "gene466",  "gene256",  "gene255", 
                          "gene288",  "gene285",  "gene286",  "gene283",  "gene241",  "gene299",  "gene282")
  expect_equal(gene_names, gene_name_to_check)
})
