#---Find Largest CC
#' Find Largest Connected Component of the graph
#'
#' @param GG igraph object to analyze
#'
#' @return igraph representation LCC
#' @export
#' @import igraph
#'
#' @examples
#' g1 <- make_star(10, mode="undirected") %du% make_ring(7) %du% make_ring(5)
#' lcc<-findLCC(g1)
#' summary(lcc)
findLCC <- function(GG){

    dec <- decompose(GG,mode = 'weak')
    vc<-vapply(dec,vcount, FUN.VALUE = 0)
    resG<-dec[[which.max(vc)]]
    return(resG)

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
#' GG <- igraph::read_graph(file, format="gml")
#' gg<-findLCC(GG)
#' gg <- addEdgeAtts(GG, gg)
#' edge_attr_names(gg)
addEdgeAtts <- function(GG, gg){

    ATTS <- names(edge.attributes(GG))

    if( !is.null(ATTS) ){

    ed <- as_edgelist(gg)
    M <- length(E(gg))
    ED <- as_edgelist(GG)

    VALUES <- list()

    for( a in seq_along(ATTS) ){
        VALUES[[a]] <- edge_attr(GG, ATTS[a], E(GG))
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
        gg <- set_edge_attr(gg, ATTS[a], E(gg), as.character(RES[, a]))
    }

    }

    return(gg)

}

#' Build network from data.table
#'
#' Wrapper for \code{\link[igraph]{graph_from_data_frame}} function which will
#' always return the largest connect component for a given network \code{ff}.
#' The function will also annotated the edges in \code{ff} with PubMed data
#' from \code{kw} if provided.
#'
#'
#' @param ff network structure data.frame with first two columns defining the
#' network edge nodes
#' @param kw pmid keyword annotation data.frame. If `NA`
#' no annotation will be added
#' @param LCC if TRUE only largest connected component is returned
#' @param simplify if TRUE loops and multiple edges will be removed
#'
#' @return igraph object of the largest connected component
#' @export
#' @import igraph
#'
#' @examples
#' f<-data.frame(A=c('A', 'A', 'B', 'D'), B=c('B', 'C', 'C', 'E'))
#' gg<-buildNetwork(f)
#' V(gg)$name
buildNetwork<-function(ff, kw=NA,LCC=TRUE,simplify=TRUE){
    #--- build raw graph
    GG <- graph_from_data_frame(ff[, seq_len(2)], directed=FALSE)
    if( !is.na(kw) ){
    GG <- set_edge_attr(GG, "METHOD", E(GG), as.character(ff[, 3]))
    GG <- set_edge_attr(GG, "TYPE", E(GG), as.character(ff[, 7]))

    PMIDS <- ifelse(!grepl("unassigned", ff[, 4]),
                    sprintf("PMID:%s", ff[, 4]), ff[, 4])
    GG <- set_edge_attr(GG, "PUBMED", E(GG), PMIDS)

    YEARS <- kw[match(gsub("PMID:", "", E(GG)$PUBMED), kw[, 1]), 3]
    YEARS <- ifelse(is.na(YEARS), "na", YEARS)
    GG <- set_edge_attr(GG, "YEAR", E(GG), YEARS)
    #---

    }
    #--- build igraph, removing multiple edges and loops
    gg <- simplify(GG, remove.multiple=TRUE, remove.loops=TRUE)
    #---Find Largest CC
    gg  <- findLCC(gg)
    gg <- addEdgeAtts(GG, gg)
}


#' Calculate sparsness of the graph.
#'
#' For a simple unweighted, undirected graph G(N,E).
#' Network sparseness is defined as the ratio of the actual number of graph
#' edges (E) to the maximum number of edges possible in a graph with same number
#' of vertices (N):  E/binom(N,2)
#'
#' @param gg graph to evaluate
#'
#' @return sparsness value
#' @export
#'
#' @examples
#' file <- system.file("extdata", "PPI_Presynaptic.csv", package = "BioNAR")
#' tbl <- read.csv(file, sep="\t")
#' gg <- buildNetwork(tbl)
#' calcSparsness(gg)
calcSparsness<-function(gg){
    N<-igraph::vcount(gg)
    E<-igraph::ecount(gg)
    sp<-2.0*E/(N*(N-1))
    return(sp)
}
