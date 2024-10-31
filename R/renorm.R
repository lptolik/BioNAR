#' Calculates Magnetic Laplacian for directed network
#'
#' Algorithm is based upon following publication:
#'
#' ref: https://arxiv.org/pdf/1606.08266
#' ref: https://doi.org/10.1016/j.acha.2017.01.004
#' Magnetic Eigenmaps for the Visualization of Directed Networks
#' Micha ̈el Fanuel, Carlos M. Ala ́ız, A ́ngela Fern ́andez, Johan A. K. Suykens
#' KU Leuven, Department of Electrical Engineering (ESAT), Kasteelpark Arenberg 10, B-3001 Leuven, Belgium
#'
#' @param gg igraph object
#' @param g an electric charge parameter, which value should be in [0,0.5].
#'
#' @return directed graph Magnetic Laplacian matrix
#' @export
#'
#' @examples
#' data(foodwebs,package='igraphdata')
#' g<-foodwebs[['ChesLower']]
#' upgrade_graph(g)
#' H <- getHmatrix(g)
getHmatrix <- function(gg, g=0.5){



    ## 0 <= g < 0.5
    ## g = 1/3, 1/4, 2/5

    ## Adjacency Matrix
    adj = as.matrix(as_adjacency_matrix(gg))
    W   = .5 * (adj + t(adj))
    D   = diag(rowSums(W))

    ## Initialize the Hermitian matrix
    N <- nrow(adj)
    H <- matrix(exp(1i*2*pi*g*0), N, N)

    ## Build the Hermitian matrix
    for (ii in 1:N) {
        rindx <- adj[ii, ]
        cindx <- which(rindx > 0)
        for (jj in cindx) {
            i <- ii
            j <- jj

            adj_ij <- adj[i, j]
            adj_ji <- adj[j, i]

            if (adj_ij > 0 && adj_ji == 0) {
                H[i,j] = exp(1i*2*pi*g*-1)
                H[j,i] = exp(1i*2*pi*g*1)
            }

            if (adj_ij == 0 && adj_ji > 0) {
                H[i,j] = exp(1i*2*pi*g*1)
                H[j,i] = exp(1i*2*pi*g*-1)
            }

        }
    }

    L  = D - H*W

}


#' Prepare Laplacian matrix from graph
#'
#'
#' @param gg igraph object
#' @param el an electric charge parameter for \code{\link{getHmatrix}}
#' @param weights an optional weights for weighted Laplacian matrix
#' @param type type of Laplacian to calculate for directed graph
#'
#' @return
#' @importFrom igraph laplacian_matrix
#' @seealso getHmatrix
#' @export
#'
#' @examples
getLaplacian <- function(gg,el=0.5,weights = NULL,type = c('igraph','magnetic')){
    type <- match.arg(type)
    L = switch(type,
        igraph = igraph::laplacian_matrix(gg,weights = weights),
        magnetic = getHmatrix(gg,g=el)
    )

    return(L)

}

#' Calculate eigenvalues and eigenvectors of Laplacian matrix
#'
#' @param x
#' @param only.values
#' @param inv.vec
#'
#' @return
#' @export
#'
#' @examples
getEigen <- function(x, only.values=TRUE, inv.vec=FALSE){

    ## Eigenvalues of L
    Eigen = eigen(x, only.values=only.values, symmetric=TRUE)
    E     = Eigen$values
    if( only.values ){
        V    = NULL
        Vinv = NULL
    } else {
        V    = Eigen$vectors
        Vinv = NULL
        if( inv.vec ){ Vinv=solve(V); }
    }

    return(list(E=E, V=V, Vinv=Vinv))
}

data_cdf_probs <- function(x){
    ## x ==> graph node degree frequency
    ## as.vector(degree(gg))
    ## Use for plotting
    ## https://github.com/csgillespie/poweRlaw/blob/main/R/dist_data_cdf_methods.R
    p = x/sum(x)
    p = rev(cumsum(rev(p)))
    p
}

meta.graph <- function(rho,complex=FALSE){

    N    = dim(rho)[1]
    adj  = NULL
    if( complex ){
        re  = diag(Re(rho)) ## Probability of being at node i at time t
        im  = Im(rho)       ## Probability flow between node i to j at time t
        adj = matrix(0,N,N)

        for( i in 1:N ){
            edges   = im[i,]/min(re[i],re)
            adj[i,] = ifelse(edges-1>0 & sign(edges),1,0)
        }

    } else {
        mn   = diag(rho)
        adj  = sapply(1:N, function(i) meta.binary(x=rho[i,],
                                                   mn=mn,
                                                   mn_ii=mn[i]) )
    }

    return(adj)
}

get.supernodes <- function(adj){

    ## build network from rho...
    ## remember to add node names to columns/row of rho.
    cc = igraph::graph_from_adjacency_matrix(adj,
                                             mode=ifelse(isSymmetric(adj),
                                                         "undirected",
                                                         "directed"))

    ## function clusters gives a bit more information than decompose
    clusters  = igraph::clusters(cc)

}

coarse.grain.graph <- function(gg, supernodes){

    ## group nodes into supernodes
    gg2 = igraph::contract(gg, supernodes)

    ## remove self-loops and multiple edges from
    gg2 = igraph::simplify(gg2, remove.loops=T, remove.multiple=T)

}


## Laplacian Renormalisation Group in real-space
real.LRG <- function(e, v, vinv=NULL, L=NULL, t, gg, complex=FALSE,
                     method=c("eigen", "balanced", "square"),
                     expm_method=c("Higham08.b"), tol=1e-5, order=1){

    method <- match.arg(method)

    # method = c("approx", "eigen", "balanced")
    # expm_methods = c("Higham08.b", "Higham08",
    #          "AlMohy-Hi09", "Ward77", "PadeRBS",
    #.         "Pade", "Taylor", "PadeO", "TaylorO",
    #          "R_Eigen", "R_Pade", "R_Ward77",
    #.         "hybrid_Eigen_Ward")

    ## gg   == original graph
    ## e    == eigenvalues of laplacian
    ## v    == eigenvectors of laplacian
    ## t    == perform coase-graining at tau
    ## herm == is eigenvalues/vectors from complex laplacian

    gn = V(gg)$name

    ## Calculate rho(tau)
    rho = switch(method,
                 "approx"={ cal.rho.approx(L=L, t=t, order=order)},
                 "eigen"={cal.rho.eigen(e=e, v=v, vinv=vinv, t=t)},
                 "balanced"={cal.rho.balanced(L=L, t=t, method=expm_method, tol=tol, order=order)},
                 "square"={cal.rho.square(L=L, t=t, order=order)},
                 cal.rho.eigen(e=e, v=v, vinv=vinv, t=t)
    )
    #rho = cal.rho(e=e, v=v, vinv=vinv, t=t)

    ## build meta-graph
    adj  = meta.graph(rho=rho, complex=complex)

    colnames(adj) = gn
    rownames(adj) = gn

    sn = get.supernodes(adj=adj)

    ## record node mapping between levels
    mapping = cbind(gn, sn$membership)

    gg2 = coarse.grain.graph(gg=gg, supernodes=sn$membership)

    return(list(rho=rho, meta.adj=adj, mapping=mapping,
                supernodes=sn$csize, gg=gg2))

}

