% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/runCentrality.R
\name{getRandomGraphCentrality}
\alias{getRandomGraphCentrality}
\title{Generate a random graph  that mimic somehow properties of the input graph
and calls \link{\code{getCentralityMatrix}} to
calculate all available centrality measires. There are four different}
\usage{
getRandomGraphCentrality(gg, type = c("gnp", "pa", "cgnp", "rw"), ...)
}
\arguments{
\item{gg}{template graph to mimic}

\item{type:}{\itemize{
\item gnp -- G(n,p) Erdos-Renyi model
\item pa --  Barabasi-Albert model
\item cgnp -- new random graph from a given graph by randomly a
dding/removing edges
\item rw -- new random graph from a given graph by rewiring 25\% of
edges preserving the degree distribution
}}
}
\value{
matrix of random graph vertices centrality measure.
}
\description{
Generate a random graph  that mimic somehow properties of the input graph
and calls \link{\code{getCentralityMatrix}} to
calculate all available centrality measires. There are four different
}
\examples{
data(karate,package='igraphdata')
m<-getRandomGraphCentrality(karate,'pa')
}