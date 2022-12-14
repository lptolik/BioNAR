% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sampling.R
\name{getClusterSubgraphByID}
\alias{getClusterSubgraphByID}
\title{Return induced subgraph for cluster}
\usage{
getClusterSubgraphByID(clID, gg, mem)
}
\arguments{
\item{clID}{cluster ID to extracte}

\item{gg}{graph to analyze}

\item{mem}{membership vector}
}
\value{
induced subgraph as igraph object
}
\description{
Function reads in a graph \code{gg}, vertex cluster membership vector
\code{mem}, and returns an induced subgraph given a cluster membership
number 'clID'.
}
\examples{
data(karate,package='igraphdata')
alg<-'louvain'
c<-getClustering(karate,alg = alg)
gc3<-getClusterSubgraphByID(3,karate,membership(c))
#plot(gc3,vertex.label=V(gc3)$name)
}
