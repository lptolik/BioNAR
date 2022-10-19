library(BioNAR)
library(testthat)
file <- system.file("extdata", "PPI_Presynaptic.gml", package = "BioNAR")
gg <- igraph::read.graph(file, format="gml")
louvainG<-induced_subgraph(gg,V(gg)[V(gg)$louvain%in%c(4,8,10,12)])

test_that('Scale',{
    expect_equal(BioNAR:::scale(1:11),
                 c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0))
})

test_that('Karate Bridgenes',{
    data(karate, package='igraphdata')
    set.seed(100)
    g <- calcClustering(karate, 'louvain')
    cnmat <- makeConsensusMatrix(g, N=10, alg = 'louvain', type = 2, mask = 10)
    br<-getBridgeness(g, alg = 'louvain', cnmat)
    expect_equal(dim(br),c(34,2))
    expect_equal(br$BRIDGENESS.louvain[br$ID=='John A'],0.2407267,
                 tolerance = 0.01)
})

test_that('Presynaptic Bridgenes',{
    set.seed(100)
    cnmat <- makeConsensusMatrix(louvainG, N=10, alg = 'louvain',
                                 type = 2, mask = 10)
    br<-getBridgeness(louvainG, alg = 'louvain', cnmat)
    expect_equal(dim(br),c(212,3))
    expect_equal(br$BRIDGENESS.louvain[br$GENE.NAME == 'ACTN2'],0.2385256,
                 tolerance = 0.01)
    expect_error(getBridgeness(louvainG, alg = 'lec',cnmat),
                 '.*calcClustering.*')
})

test_that('Norm Modularity',{
    set.seed(100)
    nm<-normModularity(gg, alg='louvain',Nint=10)
    expect_equal(nm,0.01347063,tolerance = 0.001)
})
