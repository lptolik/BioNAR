library(BioNAR)
library(testthat)
library(synaptome.db)
#context("Testing network creation")
file <- system.file("extdata", "PPI_Presynaptic.gml", package = "BioNAR")
gg <- igraph::read.graph(file, format="gml")
louvain4<-induced_subgraph(gg,V(gg)[V(gg)$louvain==4])

test_that('LCC finder',{
    g1 <- make_star(10, mode="undirected") %du% make_ring(7) %du% make_ring(5)
    lcc<-findLCC(g1)
    expect_equal(vcount(lcc),10)
    expect_equal(ecount(lcc),9)
    expect_equal(calcSparsness(gg),0.003688761,tolerance = 1e-5)
    expect_equal(calcSparsness(louvain4),0.03520164,tolerance = 1e-5)
})

test_that('annotations',{
    res<-BioNAR:::findTERM(V(gg)$TopOntoOVG,'FDR')
    expect_false(any(res))
    res<-BioNAR:::findTERM(V(gg)$TopOntoOVG,'FTD')
    expect_true(any(res))
    expect_equal(length(which(res)),84)
    gga<-gg
    E(gga)$att<-1:ecount(gg)
    ggr <- addEdgeAtts(gga, louvain4)
    expect_equal(E(ggr)$att[1:4],c("870",  "888",  "889","880"))
})

test_that('Build network',{
    f<-data.frame(A=c('A', 'A', 'B', 'D'), B=c('B', 'C', 'C', 'E'))
    g<-buildNetwork(f)
    expect_equal(vcount(g),3)
    expect_equal(ecount(g),3)

    # cid<-match('Presynaptic', getCompartments()$Name)
    # geneTable<-getAllGenes4Compartment(cid)
    # g<-buildFromSynaptomeByEntrez(geneTable$HumanEntrez)
    # expect_equal(vcount(g),2275)
    # expect_equal(ecount(g),9190)
    #
    # g<-buildFromSynaptomeGeneTable(geneTable)
    # expect_equal(vcount(g),2274)
    # expect_equal(ecount(g),9189)

})
