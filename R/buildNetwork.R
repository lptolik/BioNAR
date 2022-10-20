#---Find Largest CC
#' Find Largest Connected Component of the graph
#'
#' @param GG igraph object to analyze
#'
#' @return igraph representation largest CC
#' @export
#' @import igraph
#'
#' @examples
#' g1 <- make_star(10, mode="undirected") %du% make_ring(7) %du% make_ring(5)
#' lcc<-findLCC(g1)
#' summary(lcc)
findLCC <- function(GG){

    dec <- decompose.graph(GG)
    d <- 1
    CC <- length(V(dec[[1]]))
    for( i in seq_along(dec) ){
    if(length(V(dec[[i]])) > CC){
        d <- i
        CC <- length(V(dec[[i]]))
    }
    }

    GG  <- decompose.graph(GG)[[d]]
    return(GG)

}

#' Find term in the Annotation List
#'
#' @param eatt annotation list
#' @param TERMS vector of terms
#'
#' @return logical vector
#' @noRd
findTERM <- function(eatt, TERMS){

    eatt  <- as.vector(eatt)
    found <- rep(FALSE, length(eatt))

    for( t in seq_along(TERMS) ){
    temp  <- grepl(TERMS[t], eatt)
    found <- as.logical(found) | as.logical(temp)
    }

    return(found)

}

#
# filterPMIDs <- function(TERMS=NULL){
#
#   filterIDs <- c()
#
#   if( !is.null(TERMS ) && length(TERMS) != 0 ){
#
#     pmids <- read.delim("/afs/inf.ed.ac.uk/user/c/cmclean5/ownCloud/
#Synaptic_proteome/anaysis_17_05_2019/mined_PPIs/pmid_keywords.csv", sep="\t",
#header=TRUE)
#
#     indX <- list()
#
#     for( i in 1:length(TERMS) ){
#       N = length(names(indX))
#       indX[[N+1]] <- grepl(TERMS[i], pmids[,4])
#       names(indX)[N+1] <- sprintf("%s_title",TERMS[i])
#       N = length(names(indX))
#       indX[[N+1]] <- grepl(TERMS[i], pmids[,5])
#       names(indX)[N+1] <- sprintf("%s_keywords",TERMS[i])
#     }
#
#     exc <- indX[[1]]
#     for( i in 2:length(names(indX)) ){
#       exc <- exc | indX[[i]]
#     }
#
#     filterIDs <- as.vector(pmids[exc,1])
#
#   }
#
#   return(filterIDs)
#
# }
#

#--- add attributes to igraph edges from it raw file
#' Copy edge attributes from one graph to another
#'
#' @param GG igraph object, source of attributes
#' @param gg igraph object, attributes recipient
#'
#' @return annotated version of \code{gg} igraph object
#' @export
#'
#' @examples
#' file <- system.file("extdata", "PPI_Presynaptic.gml", package = "BioNAR")
#' GG <- igraph::read.graph(file, format="gml")
#' gg<-findLCC(GG)
#' gg <- addEdgeAtts(GG, gg)
#' edge_attr_names(gg)
addEdgeAtts <- function(GG, gg){

    ATTS <- names(edge.attributes(GG))

    if( !is.null(ATTS) ){

    ed <- get.edgelist(gg)
    M <- length(E(gg))
    ED <- get.edgelist(GG)

    VALUES <- list()

    for( a in seq_along(ATTS) ){
        VALUES[[a]] <- get.edge.attribute(GG, ATTS[a], E(GG))
        names(VALUES)[a] <- ATTS[a]
    }

    # cat("\n")
    # cat("scanning edges...")
    RES <-  matrix("", nrow=M, ncol=length(ATTS))

    for( e in seq_len(M) ){

        indx  <-  (ed[e, 1] == ED[, 1] & ed[e, 2] == ED[, 2]) |
        (ed[e, 1] == ED[, 2] & ed[e, 2] == ED[, 1])

        for( a in seq_along(ATTS) ){

        res <- VALUES[[a]][indx]

        if( res != "" ){
            res <- unique(res)
            if( length(res) == 1 ){
            RES[e, a] <- res
            } else {
            RES[e, a] <- paste(as.character(res), collapse=';')
            }
        }
        }
    }

    # cat("done.\n")

    for( a in seq_along(ATTS) ){
        gg <- set.edge.attribute(gg, ATTS[a], E(gg), as.character(RES[, a]))
    }

    }

    return(gg)

}

#' Build network from data.table
#'
#' Wrapper for \code{\link[igraph]{graph.data.frame}} function which will always
#' return the largest connect component for a given network \code{ff}.
#' The function will also #' annotated the edges in \code{ff} with PubMed data
#' if it exits.
#'
#'
#' @param ff network structure data.frame with first two columns defining the
#' network edge nodes
#' @param kw pmid keyword annotation data.frame. If `NA`
#' no annotation will be added
#'
#' @return igraph object of the largest connected component
#' @export
#' @import igraph
#'
#' @examples
#' f<-data.frame(A=c('A', 'A', 'B', 'D'), B=c('B', 'C', 'C', 'E'))
#' gg<-buildNetwork(f)
#' V(gg)$name
buildNetwork<-function(ff, kw=NA){
    #--- build raw graph
    GG <- graph.data.frame(ff[, seq_len(2)], directed=FALSE)
    if( !is.na(kw) ){
    GG <- set.edge.attribute(GG, "METHOD", E(GG), as.character(ff[, 3]))
    GG <- set.edge.attribute(GG, "TYPE", E(GG), as.character(ff[, 7]))

    PMIDS <- ifelse(!grepl("unassigned", ff[, 4]),
                    sprintf("PMID:%s", ff[, 4]), ff[, 4])
    GG <- set.edge.attribute(GG, "PUBMED", E(GG), PMIDS)

    YEARS <- kw[match(gsub("PMID:", "", E(GG)$PUBMED), kw[, 1]), 3]
    YEARS <- ifelse(is.na(YEARS), "na", YEARS)
    GG <- set.edge.attribute(GG, "YEAR", E(GG), YEARS)
    #---

    }
    #--- build igraph, removing multiple edges and loops
    gg <- simplify(GG, remove.multiple=TRUE, remove.loops=TRUE)
    #---Find Largest CC
    gg  <- findLCC(gg)
    gg <- addEdgeAtts(GG, gg)
}

#' Utility function to create network from
#' \code{\link[synaptome.db]{synaptome.db}} data
#'
#' @param entrez vector of EntrezIDs for network vertices
#'
#' @return  largest connected component of network defined by the gene list
#' @export
#' @import synaptome.db
#'
#' @examples
#' library(synaptome.db)
#' cid<-match('Presynaptic', getCompartments()$Name)
#' geneTable<-getAllGenes4Compartment(cid)
#' gg<-buildFromSynaptomeByEntrez(geneTable$HumanEntrez)
buildFromSynaptomeByEntrez<-function(entrez){
    geneTable<-findGenesByEntrez(entrez)
    gg<-buildFromSynaptomeGeneTable(geneTable)
    return(gg)
}

#' Utility function to create network from
#' \code{\link[synaptome.db]{synaptome.db}} data
#'
#' @param geneTable data.frame described in
#'        \code{\link[synaptome.db]{getGenesByID}}
#'
#' @return largest connected component of network defined by the gene table
#' @export
#'
#' @examples
#' library(synaptome.db)
#' cid<-match('Presynaptic', getCompartments()$Name)
#' geneTable<-getAllGenes4Compartment(cid)
#' gg<-buildFromSynaptomeGeneTable(geneTable)
buildFromSynaptomeGeneTable<-function(geneTable){
    p<-getPPIbyIDs(geneTable$GeneID, type = 'limited')
    aidx<-match(p$A, geneTable$GeneID)
    bidx<-match(p$B, geneTable$GeneID)
    gg<-buildNetwork(data.frame(A=geneTable$HumanEntrez[aidx],
                                B=geneTable$HumanEntrez[bidx]))
    return(gg)
}

#' Calculate sparsness of the graph.
#'
#' Sparsness is defined as ratio of graph edge numbers to the edge numbers of
#' the full graph of the same size.
#'
#' @param gg graph to evaluate
#'
#' @return sparsness value
#' @export
#'
#' @examples
#' library(synaptome.db)
#' cid<-match('Presynaptic', getCompartments()$Name)
#' t<-getAllGenes4Compartment(cid)
#' gg<-buildFromSynaptomeByEntrez(t$HumanEntrez)
#' calcSparsness(gg)
calcSparsness<-function(gg){
    N<-igraph::vcount(gg)
    E<-igraph::ecount(gg)
    sp<-2.0*E/(N*(N-1))
    return(sp)
}
