test_that("Check colors_for_gradient()", {
  expect_equal(colors_for_gradient(palette = "Seurat_Like"),
               c("#5D50A3", "#9FD7A4", "#FBFDBA", "#FEB163", "#A80B44"))
  
  expect_equal(colors_for_gradient(palette = "Ju1"),
               c("#A9D6E5", "#2166AC", "#000000", "#B2182B", "#FFCA3A"))
  
  expect_equal(colors_for_gradient(palette = "De1"),
               c("#d73027","#fc8d59","#fee090","#e0f3f8","#91bfdb","#253494"))
  
  expect_equal(colors_for_gradient(palette = "De2"),
               c("#FFF7FB","#ECE2F0","#D0D1E6","#A6BDDB","#67A9CF","#3690C0","#02818A","#016450"))
  
  expect_equal(colors_for_gradient(palette = "De3"),
               c("#1A1835","#15464E","#2B6F39","#757B33","#C17A70","#D490C6","#C3C1F2","#CFEBEF"))
  
  expect_equal(colors_for_gradient(palette = "De4"),
               c("#0000FF","#00FFFF","#80FF80","#FFFF00","#FF0000"))
  
  expect_equal(colors_for_gradient(palette = "De5"),
               c("#0000AA","#0000FF","#00FFFF","#80FF80","#FFFF00","#FF0000","#AA0000"))
  
  expect_equal(colors_for_gradient(palette = "De6"),
               c("#4575b4","#74add1","#abd9e9","#e0f3f8","#fee090","#fdae61","#f46d43","#d73027"))
  
  expect_equal(colors_for_gradient(palette = "De7"),
               c("#67001f","#b2182b","#d6604d","#f4a582","#fddbc7","#f7f7f7","#d1e5f0","#92c5de","#4393c3","#2166ac","#053061"))
  
  expect_equal(colors_for_gradient(palette = "De8"),
               c("#2b83ba","#abdda4","#fdae61","#d7191c"))
  
  expect_equal(colors_for_gradient(palette = "De9"),
               c("#0000BF","#0000FF","#0080FF","#00FFFF","#40FFBF","#80FF80","#BFFF40","#FFFF00","#FF8000","#FF0000","#BF0000"))
  
  expect_equal(colors_for_gradient(palette = "Je1"),
               c("#27408B", "#3A5FCD", "#3288BD", "#66C2A5","#ABDDA4", "#E6F598","#FEE08B", "#FDAE61","#F46D43","#D53E4F","#8B2323"))
})


test_that("Check discrete_palette()", {
  expect_equal(length(discrete_palette(n = 10)), 10)
  expect_equal(length(discrete_palette(n = 20)), 20)
  expect_equal(length(discrete_palette(n = 30)), 30)
  expect_equal(length(discrete_palette(n = 40)), 40)
  expect_equal(length(discrete_palette(n = 50)), 50)
  
  expect_equal(unname(discrete_palette(n = 10)), c(
    "#9F1717", "#517416", "#9D1C70", "#FFCA3A", "#9579B9", "#DA7316",
    "#1882C0", "#EF9292", "#B6E36A", "#EA8AC9"
  ))
  
  expect_equal(unname(discrete_palette(n = 20)), c(
    "#9F1717", "#B77009", "#637712", "#235485", "#872772", "#E96B63",
    "#F6AF3D", "#90C927", "#6799D6", "#D162B5", "#DA311D", "#E79808",
    "#689B26", "#426BAE", "#C52894", "#EF9D8B", "#FBD078", "#ACDE7E",
    "#97AAD8", "#EA8AC9"
  ))
  
  expect_equal(unname(discrete_palette(n = 30)), c(
    "#9F1717", "#AC5611", "#C08602", "#687911", "#226269", "#3F4A7E",
    "#802B73", "#C4426B", "#EA7A56", "#F3A73E", "#DAC933", "#7BC254",
    "#59A4DF", "#9C75B8", "#E159AE", "#DB2225", "#DA6A17", "#ECA303",
    "#8CA116", "#328A8D", "#4F64A8", "#A53896", "#DE5993", "#F0A089",
    "#F6C17A", "#EADC73", "#A9DC85", "#81BDE8", "#B197C9", "#EA8AC9"
  ))
  
  expect_equal(unname(discrete_palette(n = 40)), c(
    "#9F1717", "#A94612", "#B66E0A", "#BB8B01", "#6B7910", "#33684B",
    "#1B5787", "#4D457A", "#7D2D73", "#AE2D6D", "#E36167", "#EA8250",
    "#F1A33F", "#FFCA3A", "#ADC92C", "#72BE70", "#52A9E4", "#8484C3",
    "#B86BB6", "#E0529D", "#DC2936", "#DA4C1A", "#DF8010", "#EEA901",
    "#9FA312", "#52944F", "#1882C0", "#5560A5", "#954097", "#D22C93",
    "#E77892", "#F0A288", "#F3BA7C", "#FDD477", "#D7DE70", "#A8DB88",
    "#7FC6E3", "#99A8D7", "#BF93C9", "#EA8AC9"
  ))
  
  expect_equal(unname(discrete_palette(n = 50)), c(
    "#9F1717", "#A73C13", "#B0600F", "#BC7B05", "#AC8704", "#6D7A10",
    "#3D6C39", "#1A5F79", "#2D5082", "#554278", "#7B2E73", "#A1206F",
    "#CB496A", "#E96E60", "#EB864C", "#F1A13F", "#FBBF3B", "#D4C932",
    "#93C927", "#6DBB81", "#4EABE7", "#768FCD", "#9E75B8", "#C865B5",
    "#E04E92", "#DC2D41", "#DA3B1C", "#DA6817", "#E38C0C", "#EFAC00",
    "#AAA40F", "#669A2B", "#378B84", "#2879B9", "#595EA4", "#8C4498",
    "#C02B94", "#DC4E93", "#ED8B92", "#F0A387", "#F2B67C", "#F9CA79",
    "#F4DB75", "#CCE06D", "#A7DB8A", "#87CAD3", "#8CB4E1", "#A79CCC",
    "#C891C9", "#EA8AC9"
  ))
})
