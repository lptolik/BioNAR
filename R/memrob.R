

#' function to get the member robustness from the consensus matrices #BUG
#' fixed 14/03/12 TIS
#'
#' @param x consensus matrix
#' @param rm reference matrix
#'
#' @return Returns a list of memroblist class objects, one for each cluster,
#' and the full membership robustness matrix as a memrobmatrix class object.
#' @import clusterCons
#'
#' @noRd
#' @seealso clusterCons::memrob
memrob <- function(x, rm = data.frame()) {
    if (is(x, 'consmatrix')) {
        cmref <- x@rm
    } else{
        if (length(rm) == 0) {
            stop('You need to specify a reference matrix',
                    'for a merge consensus matrix')
        }else{
            cmref <- rm
        }
    }
     consensus <- x@cm
    lvlsCM <- levels(as.factor(cmref$cm))
    mem_rob <- matrix(0,
                        dim(consensus)[1],
                        length(lvlsCM),
                        dimnames = list(row.names(consensus),
                                        seq_along(lvlsCM)))
    for (k in seq_along(lvlsCM)) {
        for (i in seq_len(dim(consensus)[1])) {
            Ik <- row.names(cmref)[cmref$cm == k] #where k is the cluster number
            ind <- Ik[Ik != row.names(consensus)[i]]
            sigma <- apply(as.matrix(consensus[ind, i]), 2, sum)
            ei <- row.names(consensus)[i]
            Nk <- summary(as.factor(cmref$cm), maxsum = 10000)[k]
            if (sum(ei == Ik) == 1) {
                mik <- (1 / (Nk - 1)) * sigma
            }else{
                mik <- (1 / Nk) * sigma
            }
            mem_rob[i, k] <- mik
        }
    }
    mem_rob_list <- list()
    for (i in seq_len(dim(mem_rob)[2])) {
        cl_mem_rob <- (mem_rob[(cmref$cm == i), i])
        current_list <- data.frame(sort(cl_mem_rob, dec = TRUE))
        names(current_list) <- 'mem_rob'
        current_mem_rob_list <- new('memroblist', mrl = current_list)
        mem_rob_list[paste('cluster', i, sep = '')] <-
            current_mem_rob_list
    }
    mem_rob_list['resultmatrix'] <- new('memrobmatrix', mrm = mem_rob)
    mem_rob_list['algo'] <- x@a
    mem_rob_list['type'] <- class(x)
    return(mem_rob_list)
}

#' Calculate cluster robustness from consensus matrix
#'
#' @param gg igroph object
#' @param alg clustering algorithm
#' @param conmat consensus matrix
#'
#' @return data.frame that for each cluster \code{C} shows its size \code{Cn},
#' robustness \code{Crob} and robustness scaled to range [0,1]
#' \code{CrobScaled}.
#' @export
#'
#' @examples
#' data(karate,package='igraphdata')
#' alg<-'louvain'
#' gg<-calcClustering(karate,alg = alg)
#' conmat<-makeConsensusMatrix(gg,N=100,mask = 10,alg = alg,type = 2)
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
    crScales <- (crScales - min(crScales)) / (max(crScales) - min(crScales))
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
