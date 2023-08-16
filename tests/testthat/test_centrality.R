library(BioNAR)
library(testthat)
file <- system.file("extdata", "PPI_Presynaptic.gml", package = "BioNAR")
gg <- igraph::read.graph(file, format="gml")
louvain4<-induced_subgraph(gg,V(gg)[V(gg)$louvain==4])
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

test_that('Random centrality',{
    set.seed(100)
    cm<-getCentralityMatrix(karate)
    m<-getRandomGraphCentrality(karate,'pa',threads=1)
    expect_equal(m[1,2],33,ignore_attr = TRUE)
    set.seed(100)
    pFit <- fitDegree( as.vector(igraph::degree(graph=karate)),
                       Nsim=10, plot=FALSE,threads=1)
    pwr <- slot(pFit,'alpha')
    set.seed(100)
    lpa<-lapply(1:5,getRandomGraphCentrality,gg=karate,type='pa',
                power=pwr,weights = NULL)
    iDlpa<-calcCentralityInternalDistances(lpa)
    eDlpa<-calcCentralityExternalDistances(cm,lpa)
    sigPA<-evalCentralitySignificance(iDlpa,eDlpa)
    expect_equal(sapply(sigPA,function(.x).x$pval),
                 c(0.0606060606,0.6546786547,0.0036630037,0.6546786547,
                   0.0006660007,0.6546786547,0.0006660007),
                 tolerance = 1e-5,ignore_attr = TRUE)
    set.seed(100)
    lgnp<-lapply(1:5,getRandomGraphCentrality,gg=karate,type='gnp')
    iDlgnp<-calcCentralityInternalDistances(lgnp)
    eDlgnp<-calcCentralityExternalDistances(cm,lgnp)
    sigGNP<-evalCentralitySignificance(iDlgnp,eDlgnp)
    expect_equal(sapply(sigGNP,function(.x).x$pval),
                 c(0.0006660007,0.0006660007,0.0006660007,
                   0.9190809191,0.0006660007,0.0006660007,0.0006660007),
                 tolerance = 1e-5,ignore_attr = TRUE)
    set.seed(100)
    lcgnp<-lapply(1:5,getRandomGraphCentrality,gg=karate,type='cgnp')
    iDlcgnp<-calcCentralityInternalDistances(lcgnp)
    eDlcgnp<-calcCentralityExternalDistances(cm,lcgnp)
    sigCGNP<-evalCentralitySignificance(iDlcgnp,eDlcgnp)
    expect_equal(sapply(sigCGNP,function(.x).x$pval),
                 c(0.0506160506,0.0006660007,0.0006660007,0.9190809191,
                   0.0006660007,0.0006660007,0.0006660007),
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
