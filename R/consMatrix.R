
#' Build a consensus matrix from list of resampled clustering matrices outputted
#' from the function \code{\link{sampleGraphClust}}
#'
#' @details
#'
#' Function build a consensus matrix from list of membership matrices, which
#' are a three column matrix: the first column contains the vertex IDs of
#' input network; the second column the vertex IDs of the subsampled network,
#' or -1 if the vertex has been masked; the third column the cluster membership
#' of subsampled network, or -1 if vertex has been masked.
#' The randomised resampled membership matrices could be obtained from the
#' function \code{\link{sampleGraphClust}}.
#'
#' @param lcc list of membership matrices obtained from the
#' \code{\link{sampleGraphClust}}
#'
#' @return consensus matrix of Nvert X Nvert
#' @export
#' @keywords internal
buildConsensusMatrix <- function(lcc) {
    N <- NULL
    I <- NULL
    M <- NULL
    C <- NULL
    initM <- TRUE

    NJobs <- 0
    max_com <- 0
    min_com <- 500
    for (i in seq_along(lcc)) {
        tb <- lcc[[i]]
        ## make sure node id == -1 if node com == -1
        indx <- tb[, 3] == -1
        tb[indx, 2] <- -1

        if (initM) {
            N <- dim(tb)[1]
            temp <- matrix(0, nrow = N, ncol = N)
            I <- temp
            M <- temp
            initM <- FALSE
            rm(temp)
        }
        study <- calculateConsensusMat(data = tb)
        I <- I + study$I
        M <- M + study$M
        NJobs <- NJobs + 1

        rm(tb, study)
    }
    if (!is.null(N)) {
        C <- do.call(cbind, lapply(seq_len(N),
                                   function(s)
                                       matrixDiv(M[, s], I[, s])))
    }
    return(C)
}

## divide columns of matrix X by Y
matrixDiv <- function(x, y) {
    N <- length(x)
    res <- rep(0, N)
    indx <- y != 0
    if (sum(indx) != 0) {
        res[indx] <- x[indx] / y[indx]
    }

    return(res)
}

#' Function to make random resampling consensus matrix in memory
#'
#' @details
#'
#' Function to assess the robustness of network clustering. A randomisation
#' study is performed apply the same clustering algorithm to N perturbed
#' networks, and which returns the consensus matrix where each vertex pair is
#' assigned the probability of belong to the same cluster. The inputted network
#' is perturbed by randomly removing a \code{mask} percentage of edges
#' (\code{type=1}) or vertices (\code{type=2}) from the network before
#' clustering.
#'
#'
#'
#' @param gg graph to perturb
#' @param N number of perturbation steps
#' @param mask percentage of elements to perturbe
#' @param alg clustering alg.
#' @param type edges (1) or nodes (2)  to mask
#' @param weights The weights of the edges. It must be a positive numeric
#'        vector, NULL or NA. If it is NULL and the input graph has a ‘weight’
#'        edge attribute, then that attribute will be used. If NULL and no such
#'        attribute is present, then the edges will have equal weights. Set
#'        this to NA if the graph was a ‘weight’ edge attribute, but you don't
#'        want to use it for community detection. A larger edge weight means a
#'        stronger connection for this function. The weights value is ignored
#'        for the \code{spectral} clustering.
#' @param reclust logical to decide wether to invoke reclustering via
#'        \code{\link{recluster}}
#' @param Cnmax maximum size of the cluster in \code{mem} that will not be
#'        processed if reclustering is invoked
#'
#' @return consensus matrix of Nvert X Nvert
#' @export
#' @family {Robustness functions}
#' @examples
#' karate <- make_graph("Zachary")
#' # We need vertex ID in the 'name' attribute of the vertex
#' V(karate)$name<-c(LETTERS,letters)[1:vcount(karate)]
#' alg<-'louvain'
#' gg<-calcClustering(karate, alg = alg)
#' conmat<-makeConsensusMatrix(gg, N=100, mask = 10, alg = alg, type = 2)
#' dim(conmat)
makeConsensusMatrix <- function(gg,
                                N = 500,
                                mask = 20,
                                alg,
                                type,
                                weights=NULL,
                                reclust = FALSE,
                                Cnmax = 10) {
    lcc <- lapply(seq_len(N), function(.x)
        sampleGraphClust(
            gg = gg,
            mask = mask,
            alg = alg,
            type = type,
            weights=weights,
            reclust = reclust,
            Cnmax = Cnmax
        ))
    mm <- buildConsensusMatrix(lcc)
    colnames(mm) <- V(gg)$name
    rownames(mm) <- V(gg)$name
    return(mm)
}

## calculate the identity matrix "I" and tally matrix "M"
## given each node pairs community assignment
calculateConsensusMat <- function(data = NULL) {
    I <- NULL
    M <- NULL

    if (!is.null(data)) {
        N <- dim(data)[1]
        tempI <- matrix(0, nrow = N, ncol = N)
        tempM <- matrix(0, nrow = N, ncol = N)

        for (i in seq_len(N)) {
            comi <- as.numeric(data[i, 3])
            keyi <- data[i, 2]
            jj <- seq(i, N, 1)
            if (keyi != '-1') {
                comj <- as.numeric(data[jj, 3])
                keyj <- data[jj, 2]

                ## I
                indxJ <- jj[keyj != '-1']
                Nindx <- length(indxJ)

                tempI[i, indxJ] <- as.numeric(rep(1, Nindx))
                tempI[i, i] <- as.numeric(0.5)

                ## M
                indxC <- jj[comj != -1 & comi == comj]
                Nindx <- length(indxC)

                tempM[i, indxC] <- as.numeric(rep(1, Nindx))
                tempM[i, i] <- as.numeric(0.5)

            }
        }

        M <- tempM + t(tempM)
        I <- tempI + t(tempI)

        rm(tempM, tempI)

    }

    return(list(M = M, I = I))

}
