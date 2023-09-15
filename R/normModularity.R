##---------------------------------------------------------------------
## Ref: Takemoto, K. & Kihara, K. Modular organization of cancer s
##  ignaling networks is associated with patient survivability.
##  Biosystems 113, 149–154 (2013).
##
## To allow the comparison of network modularity with networks of different
## size and connectivity.
##
## Used the normalized network modularity value
## Qm based on the previous studies by Parter et al., 2007, Takemoto, 2012,
## Takemoto, 2013, Takemoto and Borjigin, 2011, which was defined as:
## Qm = (Qreal-Qrand)/(Qmax-Qrand)
## Where Qreal is the network modularity of a real-world signaling network and,
## Qrand is the average network modularity value obtained from 10,000
## randomized networks constructed from its real-world network. Qmax was
## estimated as: 1 − 1/M, where M is the number of modules in the real network.
##
## Randomized networks were generated from a real-world network using the
## edge-rewiring algorithm (Maslov and Sneppen, 2002).
##---------------------------------------------------------------------



#' Calculates the normalised network modularity value.
#'
#' Function to compare network Modularity of input network with networks of
#' different size and connectivity.
#'
#' Used the normalised network modularity value
#' Qm based on the previous studies by Parter et al., 2007, Takemoto, 2012,
#' Takemoto, 2013, Takemoto and Borjigin, 2011, which was defined as:
#' \deqn{Q_m = \frac{Q_{real}-Q_{rand}}{Q_{max}-Q_{rand}}}{Qm = (Qreal-Qrand)/(Qmax-Qrand)}
#' Where \eqn{Q_{real}}{Qreal} is the network modularity of a real-world signalling
#' network and, \eqn{Q_{rand}}{Qrand} is the average network modularity value obtained
#' from 10,000 randomised networks constructed from its real-world network.
#' \eqn{Q_{max}}{Qmax} was estimated as: 1 - 1/M, where M is the number of
#' modules in the real network.
#'
#' Randomised networks were generated from a real-world network using the
#' edge-rewiring algorithm (Maslov and Sneppen, 2002).
#'
#'
#' @references Takemoto, K. & Kihara, K. Modular organization of cancer
#'             signaling networks is associated with patient survivability.
#'             Biosystems 113, 149–154 (2013).
#'
#' @param gg graph object to analyze
#' @param alg clustering algorithm
#' @param Nint number of iterations
#' @param weights The weights of the edges. It must be a positive numeric
#'        vector, NULL or NA. If it is NULL and the input graph has a ‘weight’
#'        edge attribute, then that attribute will be used. If NULL and no such
#'        attribute is present, then the edges will have equal weights. Set
#'        this to NA if the graph was a ‘weight’ edge attribute, but you don't
#'        want to use it for community detection. A larger edge weight means a
#'        stronger connection for this function. The weights value is ignored
#'        for the \code{spectral} clustering.
#'
#' @return normalized modularity value
#' @export
#'
#' @examples
#' file <- system.file("extdata", "PPI_Presynaptic.csv", package = "BioNAR")
#' tbl <- read.csv(file, sep="\t")
#' gg <- buildNetwork(tbl)
#'
#' nm<-normModularity(gg, alg='louvain',Nint=10)
normModularity <- function(gg,
                           alg = c('lec', 'wt', 'fc', 'infomap',
                                   'louvain', 'sgG1', 'sgG2', 'sgG5'),
                           Nint = 1000,weights=NULL) {
    cl <- getClustering(gg, alg,weights=weights)
    Qobs <- max(cl$modularity)

    ##--- get max number of communities
    Kmax <- max(as.numeric(cl$membership))

    ##--- define Qmax
    Qmax <- 1 - 1 / Kmax

    ##--- generate random graph from our graph using rewring function,
    ##    preserving original graphs degree distribution.
    ##    Using 1000 rewing studies to get random modularity
    ##    for graph with cl clustering
    Qrnd <- 0

    for (i in seq_len(Nint)) {
        ggRnd <- igraph::rewire(graph = gg,
                                 with = keeping_degseq(loops = FALSE,
                                                       niter = 100))
        clusteRnd <- getClustering(ggRnd, alg,weights=weights)
        Qrnd <- Qrnd + max(as.numeric(clusteRnd$modularity))
        rm(ggRnd, clusteRnd)
    }

    ##--- random modularity for graph given cl clustering
    Qrnd <- Qrnd / Nint

    ##--- normalised modularity for graph given cl clustering
    Qnorm <- (Qobs - Qrnd) / (Qmax - Qrnd)
    return(Qnorm)
}
