
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
#' @import clusterCons
#' @family {Robustness functions}
#'
#' @examples
#' data(karate, package='igraphdata')
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
        data.frame(cm = as.numeric(igraph::get.vertex.attribute(gg, alg,
                                                                V(gg))))
    cm           <- data.frame(conmat)

    names(cm)    <- rownames(rm)

    rownames(cm) <- rownames(rm)

    cm           <- as.matrix(cm)

    kk     <- max(as.numeric(as.vector(rm$cm)))

    ##--- make the consensus matrix object for clusterCons so you can use
    ##--- its functions
    out <- new(
        'consmatrix',
        cm = cm,
        rm = rm,
        k = kk,
        a = alg
    )


    ##--- get cluster robustness values from clusterCons
    cr <- clusterCons::clrob(out)


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
