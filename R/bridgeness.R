#' Calculate bridginess from consensus matrix
#'
#' Bridginess takes into account a vertices shared community membership
#' together with its local neighbourhood. It was proposed in
#' Nepusz et al., 2008 <doi:10.1103/PhysRevE.77.016107>.
#'
#' Algorithm assume that clustering was already calculated and its
#' membership is stored in the appropriate vertex attribute. If
#' \link{\code{alg}} attribute is not found error will be issued.
#'
#' @param gg igraph object
#' @param alg clustering algorithm
#' @param conmat consensus matrix calculated with that algorithm
#'
#' @return data.frame with first column contains vertex ID, if GeneName
#'         attribute assigned to the vertices its value will be stored as a
#'         second column, the last column contains bridginess values for the
#          selected clustering algorithm.
#' @export
#' @importFrom igraph get.edgelist get.vertex.attribute
#' @examples
#' library(BioNAR)
#' data(karate,package='igraphdata')
#' gg <- calcClustering(karate, 'louvain')
#' conmat <- makeConsensusMatrix(gg, N=10,alg = 'louvain', type = 2, mask = 10)
#' br<-getBridgeness(gg,alg = 'louvain',conmat)
getBridgeness <- function(gg, alg, conmat) {
    #---number of vertices/genes
    N    <- length(V(gg))
    #---number of edges/PPIs
    M    <- length(E(gg))
    if (!alg %in% names(vertex.attributes(gg))) {
        stop(
            'Clustering membership attribute "',
            alg,
            '" not found.\n',
            'Please run calcClustering and makeConsensusMatrix before ',
            'bridginess calculations.\n'
        )
    }
    #---container column names
    if ("GeneName" %in% names(vertex.attributes(gg))) {
        CN   <- c('ID',
                    'GENE.NAME',
                    sprintf("BRIDGENESS.%s", alg))
        FROM <- 3
    } else{
        CN   <- c('ID',  sprintf("BRIDGENESS.%s", alg))
        FROM <- 2
    }
    #---container to store Bridgeness for algorithm 'alg'
    meas <- matrix(0, nrow = N, ncol = length(CN))
    colnames(meas) <- CN
    meas[, 1] <- as.character(V(gg)$name)
    if ("GeneName" %in% names(vertex.attributes(gg))) {
        meas[, 2] <- as.character(V(gg)$GeneName)
    }
    ##get consensus matrix indices for each edge in edgelist
    indA <- match(get.edgelist(gg)[, 1], rownames(conmat))
    indB <- match(get.edgelist(gg)[, 2], rownames(conmat))
    dat  <- data.frame(indA, indB)
    ##get community assigned to each vertex in edgelist from the algorithm 'alg'
    elA <- get.vertex.attribute(gg, alg,
                                V(gg))[match(get.edgelist(gg)[, 1],
                                                V(gg)$name)]
    elB <- get.vertex.attribute(gg, alg,
                                V(gg))[match(get.edgelist(gg)[, 2],
                                                V(gg)$name)]
    ##for each edge record the community assigned to each vertex and it's
    ##consensus matrix value
    ed      <- matrix(ncol = 6, nrow = length(E(gg)))
    ed[, 1]  <- igraph::get.edgelist(gg)[, 1]
    ed[, 2]  <- igraph::get.edgelist(gg)[, 2]
    ed[, 3]  <- elA
    ed[, 4]  <- elB
    ed[, 5]  <-
        apply(dat, 1, function(x, mat)
            mat[x[1], x[2]], mat = conmat)
    ed[, 6]  <- (as.numeric(elA) - as.numeric(elB))
    ##maximum number of communities found by clustering algorithm
    Cmax  <-
        max(as.numeric(igraph::get.vertex.attribute(gg, alg, V(gg))))
    ##loop over each vertex in the graph
    for (i in seq_along(V(gg))) {
        ##get edges belonging to the i'th veretx
        ind <-
            which(ed[, 1] == V(gg)$name[i] | ed[, 2] == V(gg)$name[i])
        ##get community belonging to the i'th vertex
        c <- igraph::get.vertex.attribute(gg, alg, V(gg))[i]
        ##reorder edge communities, so ed[,3] equals current community no: 'c'
        for (k in seq_along(ind)) {
            if (ed[ind[k], 6] != 0 && ed[ind[k], 4] == c) {
                ed[ind[k], 4] <- ed[ind[k], 3]
                ed[ind[k], 3] <- c
            }
        }
        ##number of communities i'th vertex is connected too (via it's edges)
        cc <- unique(ed[ind, 4])
        ##use sum of consensus values to calculate the likelihood of i'th
        ##vertex beloning to to k'th community.
        prob <- vector(length = length(cc))
        for (k in seq_along(cc)) {
            prob[k] <-
                sum(as.numeric(ed[which(ed[ind, 4] == cc[k]), 5])) / length(ind)
        }
        ##normalise
        prob <- prob / sum(prob)
        ##calculate bridgeness of i'th vertex
        ##Fuzzy communities and the concept of bridgeness in complex networks,
        ##T. Nepusz, arXiv, 2007
        b    <- sum((prob - 1 / Cmax) * (prob - 1 / Cmax))
        Kzero <- Cmax - length(cc)
        b <- b + sum(rep((1 / (Cmax * Cmax)), times = Kzero))
        ##store values
        ##BRIDGENESS.
        meas[i, (FROM)]  <- 1 - sqrt(Cmax / (Cmax - 1) * b)
    }
    res <- as.data.frame(meas)
    res[, FROM] <- as.numeric(res[, FROM])
    return(res)
}

#### Code for plot ####
scale <- function(x, VALUE = NULL) {
    x <- as.numeric(as.vector(x))
    
    xmin <- min(x, na.rm = TRUE)
    xmax <- max(x, na.rm = TRUE)
    
    if (is.null(VALUE)) {
        x  <- x - xmin
        x  <- ifelse(!is.na(x), x / (xmax - xmin), NA)
        
        return(x)
    }
    
    value <- as.numeric(as.vector(value)[1])
    value <- value - xmin
    value <- ifelse(!is.na(value), value / (xmax - xmin), NA)
    return(value)
}
