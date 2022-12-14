% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/enrichment.R
\name{clusterORA}
\alias{clusterORA}
\title{Calculate annotation enrichment for clusters in the graph}
\usage{
clusterORA(g, alg, name, vid = "name", alpha = 0.1, col = COLLAPSE)
}
\arguments{
\item{g}{graph to get annotation from}

\item{alg}{cluster algorithm and membership attribute name}

\item{name}{annotation attribute name}

\item{vid}{attribute to be used as a vertex ID}

\item{alpha}{probability threshold}

\item{col}{list separation character in attribute, by
default is \code{;}}
}
\value{
A table with overrepresentation results.
Each row corresponds to a tested annotation in particular cluster.
The columns are the following:
\itemize{
\item pathway – name of the enriched term as in 'names(pathway)';
\item pval – an enrichment p-value from hypergeometric test;
\item padj – a BH-adjusted p-value;
\item overlap – size of the overlap;
\item size – size of the gene set;
\item leadingEdge – vector with overlapping genes.
\item cl – cluster ID
}
}
\description{
Calculate the cluster enrichment of a graph given a clustering algorithm
'alg' and vertex annotation attribute 'name'. Function generates an
enrichment table, one row for each cluster, containing: the cluster ID, the
cluster size, overlap of annotation terms in cluster, p.value of enrichment
using the Hypergeometric test, adjusted p.value Bonferroni correction (BH).
}
\examples{
options("show.error.messages"=TRUE)
file <- system.file("extdata", "PPI_Presynaptic.gml", package = "BioNAR")
g <- igraph::read.graph(file, format="gml")
anL<-getAnnotationVertexList(g, 'TopOntoOVGHDOID')
res<-clusterORA(g, alg='louvain', name='TopOntoOVGHDOID', vid='name')
andf<-unique(data.frame(ID=get.vertex.attribute(g, 'TopOntoOVGHDOID'),
Term=get.vertex.attribute(g, 'TopOntoOVG')))
rr<-merge(andf, res, by.y='pathway', by.x='ID')
rr[order(rr$cl), ]
}
