library(BioNAR)
library(testthat)
data(karate,package='igraphdata')
alg<-'louvain'
#context("Testing robustness calculations")
test_that("proper algorithm selected",{
    expect_error(getRobustness(karate,'alg',matrix(0,ncol = 1,nrow = 1)),"Membership.*")
    expect_error(getClustering(karate,'alg'),".*lec*")
})

test_that('Membership robustness',{
    set.seed(100)
    gg<-calcClustering(karate, alg = alg)
    conmat<-makeConsensusMatrix(gg, N=10, mask = 10, alg = alg, type = 2)
    expect_equal(dim(conmat),c(34,34))
    clrob<-getRobustness(gg, alg = alg, conmat)
    expect_equal(dim(clrob)[2],5)
})

test_that('clustering',{
    m<-calcMembership(karate, 'lec')
    expect_equal(dim(m),c(34,2))
    expect_equal(max(table(m$membership)),12)
    set.seed(100)
    g1<-calcAllClustering(karate)
    df<-clusteringSummary(g1)
    expect_equal(dim(df),c(9,11))
})

test_that('spinglass on disconnected graph',{
    g1 <- make_star(10, mode="undirected") %du% make_ring(7) %du% make_ring(5)
    expect_warning(m<-calcMembership(g1, 'sgG1'),'.*NULL.*')
    expect_equal(dim(m)[1],0)
    expect_warning(g2<-calcClustering(g1, 'sgG1'),'.*NULL.*')
    expect_identical(g2,g1)

})
