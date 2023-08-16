library(BioNAR)
library(testthat)
library(BiocParallel)
file <- system.file("extdata", "PPI_Presynaptic.gml", package = "BioNAR")
gg <- igraph::read.graph(file, format="gml")
louvain4<-induced_subgraph(gg,V(gg)[V(gg)$louvain==4])
if(.Platform$OS.type=="windows"){
  register(SerialParam(),default = TRUE)
}
mnSP<-c(3.939, 5, 5.091, 5.667, 6.545, 6.576, 6.515, 5.485, 5.03, 5.242, 5.697,
        6.848, 4.879, 5.242, 5.848, 7.273, 9.212, 5.182, 5.333, 4.03, 4.909,
        5.667, 6.303, 7.212, 6.333, 8.091, 5.909, 6.394, 4.939, 5.727, 5.758,
        4.788, 5, 3.97)
sdSP<-c(1.694, 1.871, 2.199, 1.882, 2.063, 2, 2.078, 1.752, 1.447, 2.122,
        1.879, 1.822, 1.781, 1.458, 2.063, 2.254, 2.595, 1.828, 2.189, 1.51,
        1.99, 1.963, 2.229, 2.147, 1.831, 2.052, 2.006, 2.474, 1.731, 1.989,
        2.194, 1.474, 2.236, 1.912)
data(karate,package='igraphdata')
data(macaque,package='igraphdata')

test_that('calcCentrality',{
    gc<-calcCentrality(louvain4)
    idx<-match(c("DEG", "BET", "CC", "SL", "mnSP", "PR", "sdSP"),
               vertex_attr_names(gc))
    expect_false(any(is.na(idx)))
})

test_that('calcDirectedCentrality',{
    gc<-calcCentrality(macaque)
    idx<-match(c("DEG", "iDEG", "oDEG", "BET", "dBET", "CC", "SL",
                 "mnSP", "PR", "dPR", "sdSP"),
               vertex_attr_names(gc))
    expect_false(any(is.na(idx)))
})

test_that('calcCentralitySerial',{
    gc<-calcCentrality(louvain4,BPparam=SerialParam())
    idx<-match(c("DEG", "BET", "CC", "SL", "mnSP", "PR", "sdSP"),
               vertex_attr_names(gc))
    expect_false(any(is.na(idx)))
})

test_that('calcDirectedCentralitySerial',{
    gc<-calcCentrality(macaque,BPparam=SerialParam())
    idx<-match(c("DEG", "iDEG", "oDEG", "BET", "dBET", "CC", "SL",
                 "mnSP", "PR", "dPR", "sdSP"),
               vertex_attr_names(gc))
    expect_false(any(is.na(idx)))
})

test_that('Apply matrix',{
    g1 <- make_star(10, mode="undirected")
    V(g1)$name <- letters[1:10]
    m<-cbind(ID=letters[1:10],capital=LETTERS[1:10])
    g2<-applpMatrixToGraph(g1,m)
    expect_equal(V(g2)$capital,LETTERS[1:10])
    m<-data.frame(ID=letters[3:10],partial=LETTERS[1:8],digits=1:8)
    g3<-applpMatrixToGraph(g2,m)
    expect_equal(length(which(is.na(V(g3)$partial))),2)
    expect_equal(V(g3)$digits[3],1)
})

test_that('Errors and warnings',{
    g1 <- make_star(10, mode="undirected")
    V(g1)$name <- letters[1:10]
    expect_error(applpMatrixToGraph(g1,data.frame(ID=1,Val=LETTERS[1:10])),
                 '.*unique.')
})

test_that('SP centrality',{
    cm<-getCentralityMatrix(karate)
    expect_equal(cm$mnSP,mnSP,tolerance = 0.01)
    expect_equal(cm$sdSP,sdSP,tolerance = 0.01)
})

test_that('Random centrality',{
    cm<-getCentralityMatrix(karate)
    set.seed(100)
    m<-getRandomGraphCentrality(gg=karate,N=1,type='pa',
                                BPparam=SerialParam(RNGseed = 100))[[1]]
    expect_equal(m[1,2],4,ignore_attr = TRUE)
    set.seed(100)
    pFit <- fitDegree( as.vector(igraph::degree(graph=karate)),
                       Nsim=10, plot=FALSE,threads=1)
    pwr <- slot(pFit,'alpha')
    lpa<-getRandomGraphCentrality(gg=karate,N=5,type='pa',
                power=pwr,weights = NULL,BPparam=SerialParam(RNGseed = 100))
    iDlpa<-calcCentralityInternalDistances(lpa)
    eDlpa<-calcCentralityExternalDistances(cm,lpa)
    sigPA<-evalCentralitySignificance(iDlpa,eDlpa)
    expect_equal(sapply(sigPA,function(.x).x$pval),
                 c(0.000666000666000666, 0.000666000666000666,
                   0.000333000333000333, 0.000666000666000666,
                   0.000666000666000666, 0.0193140193140193,
                   0.000666000666000666),
                 tolerance = 1e-5,ignore_attr = TRUE)
    lgnp<-getRandomGraphCentrality(gg=karate,N=5,type='gnp',
                                   BPparam=SerialParam(RNGseed = 100))
    iDlgnp<-calcCentralityInternalDistances(lgnp)
    eDlgnp<-calcCentralityExternalDistances(cm,lgnp)
    sigGNP<-evalCentralitySignificance(iDlgnp,eDlgnp)
    expect_equal(sapply(sigGNP,function(.x).x$pval),
                 c(0.000666000666000666, 0.000666000666000666,
                   0.000666000666000666, 0.654678654678655,
                   0.000666000666000666, 0.000666000666000666,
                   0.000666000666000666),
                 tolerance = 1e-5,ignore_attr = TRUE)
    lcgnp<-getRandomGraphCentrality(gg=karate,N=5,type='cgnp',
                                    BPparam=SerialParam(RNGseed = 100))
    iDlcgnp<-calcCentralityInternalDistances(lcgnp)
    eDlcgnp<-calcCentralityExternalDistances(cm,lcgnp)
    sigCGNP<-evalCentralitySignificance(iDlcgnp,eDlcgnp)
    expect_equal(sapply(sigCGNP,function(.x).x$pval),
                 c(0.154845154845155, 0.000666000666000666,
                   0.654678654678655, 0.350649350649351,
                   0.000666000666000666, 0.0193140193140193,
                   0.000666000666000666),
                 tolerance = 1e-5,ignore_attr = TRUE)

})


test_that('Cluster subgraph',{
    gc3<-getClusterSubgraphByID(4,gg,as.numeric(vertex_attr(gg, 'louvain')))
    expect_true(isomorphic(louvain4,gc3))
})

test_that('Layouts',{
    alg<-'louvain'
    set.seed(100)
    mem<-calcMembership(karate,alg = alg)
    set.seed(100)
    lay<-layoutByCluster(karate,mem)
    expect_equal(lay[1,],c(-10.1634700,9.7602631),tolerance = 0.001)
    set.seed(100)
    remem<-calcReclusterMatrix(karate,mem,alg,10)
    expect_equal(unlist(remem[34,c(2,3)]),c(3,4),ignore_attr = TRUE)
    set.seed(100)
    lay<-layoutByRecluster(karate,remem)
    expect_equal(lay[1,],c(6.0216997,14.6233495),tolerance = 0.001)
    cg<-getCommunityGraph(karate,mem$membership)
    expect_true(isomorphic(cg,graph_from_literal(A-B-C)))
})
