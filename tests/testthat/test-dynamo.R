library(BioNAR)
library(testthat)
#karate <- make_graph("Zachary")
data(karate, package='igraphdata')
upgrade_graph(karate)
test_that("dynamo karate weight", {
    d<-getDYNAMO(karate,attr='weight')
    expect_equal(dim(d),c(34,34))
    expect_true(inherits(d,'dgCMatrix'))
})

test_that("dynamo karate NULL", {
    d<-getDYNAMO(karate,attr=NULL)
    expect_equal(dim(d),c(34,34))
    expect_true(inherits(d,'dgCMatrix'))
})

test_that("dynamo NA weights fail", {
    g<-karate
    E(g)$weight[1]<-NA
    expect_error(getDYNAMO(g,attr='weight'),".+NAs.+")
})

test_that("dynamo NA weights fail", {
    g<-karate
    E(g)$weight<-as.character(E(g)$weight)
    E(g)$weight[1]<-".NA!"
    expect_error(getDYNAMO(g,attr='weight'),".+numeric+")
})
