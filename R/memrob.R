
#' Calculate cluster robustness from consensus matrix
#'
#' This function takes as argument a network (\code{gg}), the name of a
#' clustering algorithm (\code{alg}) which can be found in the network, and
#' a consensus matrix (\code{conmat}) generated from the clustering network.
#' The function uses the consensus matrix to generate a measure of cluster
#' robustness \eqn{C_{rob}} (\code{Crob}) for  each cluster (\code{C}) using the R function
#' \code{\link[clusterCons]{clrob}}.
#' Briefly, this is done by summing elements of the consensus matrix that
#' are found in the same cluster, and dividing this by the total number of
#' entries in the matrix:
#' \deqn{C_{rob}=\frac{2}{C_n(C_n-1)} \sum_{i,j \in I_C \atop i \le j} conmat_{i,j}}{Crob=(2/(C_n(C_n-1)))\Sigma conmat_ij; i<j, i,j \in I_C}
#' where \eqn{I_C} -- indices of vertices of the cluster \eqn{C},
#' \eqn{C_n} is the number of nodes found inside the cluster \eqn{C}.
#'
#' @param gg igroph object
#' @param alg clustering algorithm
#' @param conmat consensus matrix
#'
#' @return data.frame that for each cluster \code{C} shows
#' * its size \eqn{C_n} (\code{Cn}),
#' * robustness \eqn{C_{rob}} (\code{Crob}) and
#' * robustness scaled to range between 0 and 1 (\code{CrobScaled}).
#'
#' @export
#' @family {Robustness functions}
#'
#' @examples
#' karate <- make_graph("Zachary")
#' # We need vertex ID in the 'name' attribute of the vertex
#' V(karate)$name<-c(LETTERS,letters)[1:vcount(karate)]
#' alg<-'louvain'
#' gg<-calcClustering(karate, alg = alg)
#' conmat<-makeConsensusMatrix(gg, N=100, mask = 10, alg = alg, type = 2)
#' clrob<-getRobustness(gg, alg = alg, conmat)
#' clrob
getRobustness <- function(gg, alg, conmat) {
    if (!alg %in% igraph::vertex_attr_names(gg)) {
        stop(
            "Membership for the ",
            alg,
            "algorithm should be stored in the graph.\n",
            "See calcClustering.\n"
        )
    }
    rm <-
        data.frame(cm = as.numeric(factor(igraph::vertex_attr(gg, alg,
                                                                V(gg)))))
    cm           <- data.frame(conmat)

    names(cm)    <- rownames(rm)

    rownames(cm) <- rownames(rm)

    cm           <- as.matrix(cm)

    kk     <- max(as.numeric(as.vector(rm$cm)))

    ##--- make the consensus matrix object for clusterCons so you can use
    ##--- its functions
    # out <- new(
    #     'consmatrix',
    #     cm = cm,
    #     rm = rm,
    #     k = kk,
    #     a = alg
    # )


    ##--- get cluster robustness values from clusterCons
    #cr <- .clrob(out)
    cr <- .clrob(rm,cm)


    ##--- the scaled cluster robustness values
    crScales <- cr$rob
    crScales <-
        (crScales - min(crScales)) / (max(crScales) - min(crScales))
    oo <- data.frame(
        alg = as.character(rep(alg, length(rownames(
            cr
        )))),
        C = as.numeric(rownames(cr)),
        Cn = as.numeric(table(rm[, 1])),
        Crob = as.numeric(cr$rob),
        CrobScaled = as.numeric(crScales)
    )

}

#Function adopted from https://github.com/biomedicalinformaticsgroup/clusterCons
#due to exclusion of the clusterCons package from CRAN
#.clrob<-function (x, rm = data.frame())
.clrob<-function (cmref,cmr, rm = data.frame())
{
    # if (class(x) == "consmatrix") {
    #     cmref <- x@rm
    # }
    # else {
    #     if (length(rm) == 0) {
    #         stop("You need to specify a reference matrix for a merge consensus matrix")
    #     }
    #     else {
    #         cmref <- rm
    #     }
    # }
    # cmr <- x@cm
    clnum <- length(unique(cmref$cm))
    cl_rob = data.frame(matrix(0, clnum, 1))
    row.names(cl_rob) <- seq(1:clnum)
    names(cl_rob) <- c("rob")
    for (k in 1:clnum) {
        sm = as.matrix(cmr[row.names(cmref)[cmref$cm == k],
                           row.names(cmref)[cmref$cm == k]])
        cl_sum <- 0
        cl_size <- dim(sm)[1]
        for (i in 1:cl_size) {
            for (j in 1:cl_size) {
                if (i < j) {
                    cl_sum = cl_sum + sm[i, j]
                }
            }
        }
        meas_sum = 1/(cl_size * (cl_size - 1)/2)
        curr_cl_rob = cl_sum * meas_sum
        if (is.na(curr_cl_rob)) {
            cl_rob[k, 1] = 0
        }
        else {
            cl_rob[k, 1] = curr_cl_rob
        }
    }
    return(cl_rob)
}
