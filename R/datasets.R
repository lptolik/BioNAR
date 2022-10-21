
#' Table of protein protein interactions for presynaptic compartment
#'
#' Protein-protein interactions (PPIS) for presynaptic compartment, extracted
#' from Synaptome.db, in a csv form. Columns A and B correspond to Entrez IDs
#' for interacting proteins A and B (node names); column We contains the edge
#' weights, if available.
#'
#'
#' @name PPI_Presynaptic.csv
#' @keywords PPI_Presynaptic.csv
#' @docType data
#' @details
#'     details for PPI_Presynaptic.csv
#'
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
#' @keywords PPI_Presynaptic.gml
#' @docType data
#' @details
#'     details for PPI_Presynaptic.gml
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
#' @keywords PresynAn.csv
#' @docType data
#' @details
#'     details for PresynAn.csv
#'
NULL

#' Schizopherina related synaptic gene functional annotation.
#'
#' Annotation, manually curated from an external file: Lips et al., (2012)
#' doi:10.1038/mp.2011.117.The table has columns: the first
#' containing gene Human Entrez IDs, the second gene functional group ID terms,
#' the third gene functional group description terms; in csv format
#'
#' @name SCH_flatfile.csv
#' @keywords SCH_flatfile.csv
#' @docType data
#' @details
#'     details for SCH_flatfile.csv
#'
NULL

#' Title for diseasome.rda
#'
#' description for diseasome.rda
#'
#' @name diseasome.rda
#' @keywords diseasome.rda
#' @docType data
#' @details
#'     details for diseasome.rda
#'
NULL

#' Annotation from Gene Ontology Biological Process (GO_BP)
#'
#' Annotation, downloaded from Gene Ontology for Biological Proceess domain.
#' The table has columns: the first containing gene gene functional group ID
#' terms, the second gene functional group description terms,
#' the third - Human gene Entrez IDs; in csv format
#'
#' @name flatfile.go.BP.csv
#' @keywords flatfile.go.BP.csv
#' @docType data
#' @details
#'     details for flatfile.go.BP.csv
#'
NULL

#' Annotation from Gene Ontology Cellular Compartment (GO_CC)
#'
#' Annotation, downloaded from Gene Ontology for Cellular Compartment domain.
#' The table has columns: the first containing gene gene functional group ID
#' terms, the second gene functional group description terms,
#' the third - Human gene Entrez IDs; in csv format
#'
#' @name flatfile.go.CC.csv
#' @keywords flatfile.go.CC.csv
#' @docType data
#' @details
#'     details for flatfile.go.CC.csv
#'
NULL

#' Annotation from Gene Ontology Molecular Function (GO_CC)
#'
#' Annotation, downloaded from Gene Ontology for Molecular Function domain.
#' The table has columns: the first containing gene gene functional group ID
#' terms, the second gene functional group description terms,
#' the third - Human gene Entrez IDs; in csv format
#'
#' @name flatfile.go.MF.csv
#' @keywords flatfile.go.MF.csv
#' @docType data
#' @details
#'     details for flatfile.go.MF.csv
#'
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
#' @keywords flatfile_human_gene2HDO.csv
#' @docType data
#' @details
#'     details for flatfile_human_gene2HDO.csv
#'
NULL
