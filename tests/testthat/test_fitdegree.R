library(BioNAR)
library(testthat)
file <- system.file("extdata", "PPI_Presynaptic.gml", package = "BioNAR")
gg <- igraph::read.graph(file, format="gml")
louvain4<-induced_subgraph(gg,V(gg)[V(gg)$louvain==4])

test_that('Test log seq',{
    expect_equal(BioNAR:::lseqBy(1,100),c(1,10,100))
    yTICKS <- round(BioNAR:::lseqBy(1e-2, 1), 4)
    yLABELS <- BioNAR:::changeSciNotation(yTICKS)
    str<-c("10^{\n    phantom() - 2\n}",
           "10^{\n    phantom() - 1\n}",
           "10^{\n    0\n}" )
    expect_equal(as.character(yLABELS),str)
})

test_that('FitDegree',{
    set.seed(100)
    pFit <- fitDegree( as.vector(igraph::degree(graph=louvain4)),
                       threads=1, Nsim=5)
    expect_equal(pFit@alpha,2.61454,tolerance = 0.001)
})

# test_that('FitDegree plot',{
#     set.seed(100)
#     el<-as.vector(igraph::degree(graph=gg))
#     vdiffr::expect_doppelganger("PowerFitPlot",
#                                 fitDegree(el,Nsim = 5,threads=1,
#                                           plot=TRUE,showErr = FALSE))
# })
