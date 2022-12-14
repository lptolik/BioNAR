% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/annotateNetwork.R
\name{annotateSCHanno}
\alias{annotateSCHanno}
\title{Add SCHanno synaptic functional groups}
\usage{
annotateSCHanno(gg, anno)
}
\arguments{
\item{gg}{igraph object to annotate}

\item{anno}{annotation matrix in Pairs form}
}
\value{
annotated igraph object
}
\description{
The function loads an annotation data matrix of functional groups for
schizopherina risk genes (1) called anno, which contains three columns; the
first containing gene Entrez IDs, the second gene functional group ID terms,
the third gene functional group description terms. The function then performs
a many-to-one mapping of each matrix row to a network vertex using matching
Entrez IDs, filling the \code{SCHanno} vertices attribute.
}
\details{
References:
\enumerate{
\item Lips E, Cornelisse L, Toonen R, Min J, Hultman C, the International
Schizophernia Consortium, Holmans P, Donovan M, Purcell S, Smit A, Verhage M,
Sullivan P, Visscher P, D P: Functional gene group analysis identifies
synaptic gene groups as risk factor for schizophrenia.
Molecular Psychiatry 2012,17:996–1006.
}
}
\examples{
library(synaptome.db)
cid<-match('Presynaptic', getCompartments()$Name)
t<-getAllGenes4Compartment(cid)
gg<-buildFromSynaptomeByEntrez(t$HumanEntrez)
afile<-system.file("extdata", "SCH_flatfile.csv", package = "BioNAR")
dis    <- read.table(afile, sep="\t", skip=1, header=FALSE,
strip.white=TRUE, quote="")
agg<-annotateSCHanno(gg, dis)
}
\seealso{
getAnnotationVertexList
}
