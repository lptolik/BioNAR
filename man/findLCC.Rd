% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/buildNetwork.R
\name{findLCC}
\alias{findLCC}
\title{Find Largest Connected Component of the graph}
\usage{
findLCC(GG)
}
\arguments{
\item{GG}{igraph object to analyze}
}
\value{
igraph representation LCC
}
\description{
Find Largest Connected Component of the graph
}
\examples{
g1 <- make_star(10, mode="undirected") \%du\% make_ring(7) \%du\% make_ring(5)
lcc<-findLCC(g1)
summary(lcc)
}
