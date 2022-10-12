library(BioNAR)
library(testthat)
test_that("GO annotation is added",{
    file <- system.file("extdata", "PPI_Presynaptic.gml", package = "BioNAR")
    gg <- igraph::read.graph(file, format="gml")
    ggGO <- annotateGOont(gg,orgDB = org.Hs.eg.db, keytype = "ENTREZID")
    expect_equal(length(grep('GO',vertex_attr_names(ggGO))),6)
})
