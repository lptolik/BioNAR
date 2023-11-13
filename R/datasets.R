
#' Table of protein protein interactions for presynaptic compartment
#'
#' Protein-protein interactions (PPIS) for presynaptic compartment, extracted
#' from Synaptome.db, in a csv form. Columns A and B correspond to Entrez IDs
#' for interacting proteins A and B (node names); column We contains the edge
#' weights, if available.
#'
#'
#' @name PPI_Presynaptic.csv
#' @keywords file
#' @docType data
#' @seealso \code{\link{buildNetwork}}
NULL

#' PPI graph for presynaptic compartment
#'
#' Protein-protein interactions (PPIS) for presynaptic compartment, extracted
#' from Synaptome.db, and saved in a graph format. Graph contains node
#' attributes, such as names (Entrez IDs), Gene Names, disease association
#' (TopOntoOVG, TopOntoOVGHDOID), annotation with schizophrenia-related genes
#' (Schanno (v/c), function annotation from GO (GOBPID, GOBP, GOMFID, GOMF,
#' GOCCID, GOCC), centrality measures (DEG - degree, BET - betweenness, CC -
#' clustering coefficient,  SL - semilocal centrality, mnSP - mean shortest
#' path, PR - page rank, sdSP - standard deviation of the shortest path), and
#' clustering memberships for 8 clustering algorithms (lec, wt, fc, infomap,
#' louvain, sgG1, sgG2, sgG5)
#'
#'
#' @name PPI_Presynaptic.gml
#' @keywords file
#' @docType data
#'
NULL


#' Presynaptic genes specific functional annotation
#'
#' Presynaptic genes functional annotation derived from Boyken at al. (2013)
#' <doi:10.1016/j.neuron.2013.02.027>. The table has columns: the first
#' containing  functional group ID terms, the second - gene functional group
#' description terms, third - gene Human Entrez Ids; in csv format
#'
#' @name PresynAn.csv
#' @keywords file
#' @docType data
#' @seealso \code{\link{annotatePresynaptic}}
NULL

#' Schizopherina related synaptic gene functional annotation.
#'
#' Annotation, manually curated from an external file: Lips et al., (2012)
#' doi:10.1038/mp.2011.117.The table has columns: the first
#' containing gene Human Entrez IDs, the second gene functional group ID terms,
#' the third gene functional group description terms; in csv format
#'
#' @name SCH_flatfile.csv
#' @keywords file
#' @docType data
#' @seealso \code{\link{annotateSCHanno}}
NULL

#' Barabasi's Diseasome Network
#'
#' In the paper Goh.t al. (2007) doi:10.1073/pnas.0701361104 Barabasi with
#' colleagues published Diseasome: a network of disorders and disease genes
#' linked by known disorder–gene associations. We extract definition of the
#' genes, disorders and interactions from papers supplementary materials and
#' store it as \code{\link[igraph]{graph}} object.
#'
#' Diseasesome is a bipartite graph that have nodes of two types \code{gene}
#' and \code{disease} and links are allowed only between nodes of different
#' types. It could be projected to Human Disease Network (HDN) and Disease
#' Gene Network (DGN).
#'
#' @format
#'   A bipartite graph as \code{\link[igraph]{graph}} object.
#'
#'   Vertex attributes: \sQuote{name} for the node ID, \sQuote{Name} for the
#'   human readable node name, \sQuote{Disorder.class},
#'   \sQuote{Type} for the human readable node type,
#'   \sQuote{label} and \sQuote{shape} for plotting the graph,
#'   \sQuote{type} the node type for bipartite \code{\link[igraph]{graph}}
#'   representation.
#'
#'
#'
#' @keywords diseasome
#' @keywords graphs
#' @source Goh, K.-I. et al. The human disease network. Proc. Natl. Acad. Sci.
#' U.S.A. 104, 8685–8690 (2007).
#' https://pnas.org/doi/full/10.1073/pnas.0701361104
#' @docType data
#'
"diseasome"

#' Annotation from Gene Ontology Biological Process (GO_BP)
#'
#' Annotation, downloaded from Gene Ontology for Biological Proceess domain.
#' The table has columns: the first containing gene gene functional group ID
#' terms, the second gene functional group description terms,
#' the third - Human gene Entrez IDs; in csv format
#'
#' @name flatfile.go.BP.csv
#' @keywords file
#' @docType data
#' @seealso \code{\link{annotateGoBP}}
NULL

#' Annotation from Gene Ontology Cellular Compartment (GO_CC)
#'
#' Annotation, downloaded from Gene Ontology for Cellular Compartment domain.
#' The table has columns: the first containing gene gene functional group ID
#' terms, the second gene functional group description terms,
#' the third - Human gene Entrez IDs; in csv format
#'
#' @name flatfile.go.CC.csv
#' @keywords file
#' @docType data
#' @seealso \code{\link{annotateGoCC}}
NULL

#' Annotation from Gene Ontology Molecular Function (GO_MF)
#'
#' Annotation, downloaded from Gene Ontology for Molecular Function domain.
#' The table has columns: the first containing gene gene functional group ID
#' terms, the second gene functional group description terms,
#' the third - Human gene Entrez IDs; in csv format
#'
#' @name flatfile.go.MF.csv
#' @keywords file
#' @docType data
#' @seealso \code{\link{annotateGoMF}}
NULL

#' Human Gene Disease Associations (GDA)
#'
#' Annotation derived from Human Disease Ontology database (HDO). The table
#' contains three columns; the first containing gene Entrez IDs, the second gene
#' Human Disease Ontology (HDO) ID terms, the third gene HDO description terms;
#' in csv format
#'
#'
#' @name flatfile_human_gene2HDO.csv
#' @keywords file
#' @docType data
#' @seealso \code{\link{annotateTopOntoOVG}}
NULL
