% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sampling.R
\name{calcReclusterMatrix}
\alias{calcReclusterMatrix}
\title{Hierarchial graph clustering}
\usage{
calcReclusterMatrix(gg, mem, alg, CnMAX = 10, keepSplit = FALSE)
}
\arguments{
\item{gg}{graph to cluster}

\item{mem}{data.frame with previous level clustering results}

\item{alg}{algorithm to apply}

\item{CnMAX}{maximus size of the cluster in \code{mem} that will not be
processed}
}
\value{
remembership matrix, that contains vertex ID membership and
result of reclustering
}
\description{
Function takes graph \code{gg} and its membership matrix \code{mem}
as returned \code{calcMembership} and apply clustering algorithm \code{alg}
to all clusters larger than \code{CnMAX}
}
\examples{
data(karate,package='igraphdata')
alg<-'louvain'
mem<-calcMembership(karate,alg = alg)
remem<-calcReclusterMatrix(karate,mem,alg,10)
}