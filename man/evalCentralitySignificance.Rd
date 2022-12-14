% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/runCentrality.R
\name{evalCentralitySignificance}
\alias{evalCentralitySignificance}
\title{Compare distance distributions of internal and external distances}
\usage{
evalCentralitySignificance(dmi, dme)
}
\arguments{
\item{dmi}{distribution of internal distances between random graph
centralities}

\item{dme}{distribution of external distances between random and
original graph centralities}
}
\value{
list of lists for each centrality value in the input matrix three
element list is created where \code{ks} contains Kolmogorov-Smirnov test
result from class \code{ks.test}; \code{pval} contains Kolmogorov-Smirnov
test pvalue;
and \code{dt} contains input distribution.
}
\description{
Function to compare two distance distributions using the Kolmogorov-Smirnov
test. Where the first distance distribution is generated internally and
calculates the distance between random graph centralities. The second
distance distribution is generated externally, and measures the distance
between random and the original graph centralities.
}
\examples{
data(karate,package='igraphdata')
m<-getCentralityMatrix(karate)
gnp<-list()
for(i in 1:10){
    gnp[[i]]<-getRandomGraphCentrality(karate,type = 'gnp')
}
gnpIDist<-calcCentralityInternalDistances(gnp)
gnpEDist<-calcCentralityExternalDistances(m,gnp)

simSig<-evalCentralitySignificance(gnpIDist,gnpEDist)
sapply(simSig,function(.x).x$ks$p.value)
}
\seealso{
ks.test
}
