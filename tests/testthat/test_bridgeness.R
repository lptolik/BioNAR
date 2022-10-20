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
    agg<-calcBridgeness(g, alg = 'louvain', cnmat)
    expect_true(any(grepl('louvain',vertex_attr_names(agg))))
    expect_true(any(grepl('BRIDGENESS.louvain',vertex_attr_names(agg))))
    idx<-match(br$ID,V(agg)$name)
    expect_false(any(is.na(idx)))
    expect_equal(br$BRIDGENESS.louvain,V(agg)$BRIDGENESS.louvain[idx])
    expect_error(getBridgeness(louvainG, alg = 'lec',cnmat),
                 '.*calcClustering.*')
    expect_error(plotBridgeness(agg, alg = 'lec',VIPs=c("Mr Hi","John A")),
                 '.*SL.*')
    agg <- calcCentrality(agg)
    expect_error(plotBridgeness(agg, alg = 'lec',VIPs=c("Mr Hi","John A")),
                 '.*BRIDGENESS.lec.*')
    vdiffr::expect_doppelganger("KarateBridgenessPlot",
                                plotBridgeness(agg,alg = 'louvain',
                                               VIPs=c("Mr Hi","John A"),
                                               Xatt='SL',
                                               Xlab = "SL",
                                               Ylab = "Bridgeness",
                                               MainDivSize = 0.8,
                                               xmin = 0,
                                               xmax = 1,
                                               ymin = 0,
                                               ymax = 1,
                                               baseColor="royalblue2",
                                               SPColor="royalblue2"))
})

test_that('Presynaptic Bridgenes',{
    set.seed(100)
    cnmat <- makeConsensusMatrix(louvainG, N=10, alg = 'louvain',
                                 type = 2, mask = 10)
    br<-getBridgeness(louvainG, alg = 'louvain', cnmat)
    expect_equal(dim(br),c(212,3))
    expect_equal(br$BRIDGENESS.louvain[br$GENE.NAME == 'ACTN2'],0.2385256,
                 tolerance = 0.01)
    agg<-calcBridgeness(louvainG, alg = 'louvain', cnmat)
    expect_true(any(grepl('louvain',vertex_attr_names(agg))))
    expect_true(any(grepl('BRIDGENESS.louvain',vertex_attr_names(agg))))
    idx<-match(br$ID,V(agg)$name)
    expect_false(any(is.na(idx)))
    expect_equal(br$BRIDGENESS.louvain,V(agg)$BRIDGENESS.louvain[idx])
    agg <- calcCentrality(agg)
    vdiffr::expect_doppelganger("PresynBridgenessPlot",
                                plotBridgeness(agg,alg = 'louvain',
                                               VIPs=c('8495','22999','8927',
                                                      '8573','26059','8497',
                                                      '27445','8499'),
                                               Xatt='SL',
                                               Xlab = "SL",
                                               Ylab = "Bridgeness",
                                               MainDivSize = 0.8,
                                               xmin = 0,
                                               xmax = 1,
                                               ymin = 0,
                                               ymax = 1,
                                               baseColor="royalblue2",
                                               SPColor="royalblue2"))
})

test_that('Norm Modularity',{
    set.seed(100)
    nm<-normModularity(gg, alg='louvain',Nint=10)
    expect_equal(nm,0.01347063,tolerance = 0.001)
})

test_that('Perturbation entropy',{
    set.seed(100)
    e<- getEntropy(louvainG)
    expect_equal(dim(e),c(vcount(louvainG),5))
    expect_equal(e$DEGREE[e$ENTREZ.ID==2891],4)
    expect_equal(e$UP[e$ENTREZ.ID==2891],0.7948626,tolerance = 0.001)
    expect_equal(e$DOWN[e$ENTREZ.ID==2891],0.8231825,tolerance = 0.001)
    ove<-BioNAR:::getEntropyOverExpressed(e)
    expect_equal(dim(ove),c(2,5))
    expect_equal(as.vector(ove$ENTREZ.ID),c('10342','81488'))
    g<- calcEntropy(louvainG)
    expect_false(any(is.na(match(c('SR_UP', 'SR_DOWN'),vertex_attr_names(g)))))
    set.seed(100)
    ent <- getEntropyRate(louvainG)
    vdiffr::expect_doppelganger("EntropyPlot",
                                plotEntropy(e, subTIT = "Entropy",
                                            SRo = ent$SRo, maxSr = ent$maxSr))

})

