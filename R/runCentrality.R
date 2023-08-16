##---Semi-local Centrality (Cl)
##   Identifying influential nodes in complex networks,
##   D. Chen et al., Physica A, 2012
Semilocal <- function(gg) {
    N    <- vcount(gg)
    meas <- matrix(0, nrow = N, ncol = 3)
    for (i in seq_len(N)) {
        ids <- as.character(V(gg)[i]$name)
        neig <- igraph::neighbors(gg, v = ids, mode = "all")
        if (length(neig) > 0) {
            for (w in seq_along(neig)) {
                vnames <- as.character(V(gg)$name[neig[w]])
                neig <- c(neig,igraph::neighbors(gg,
                                                 v = vnames,
                                                 mode = "all"))
            }
            neig <- unique(neig)
            meas[i, 1] <- length(neig) - 1
        }
    }
    for (i in seq_len(N)) {
        ids <- as.character(V(gg)[i]$name)
        neig <- igraph::neighbors(gg, v = ids, mode = "all")
        meas[i, 2] <- sum(as.numeric(meas[neig, 1]))
    }
    for (i in seq_len(N)) {
        ids <- as.character(V(gg)[i]$name)
        neig <- igraph::neighbors(gg, v = ids, mode = "all")
        meas[i, 3] <- sum(as.numeric(meas[neig, 2]))
    }
    return(as.numeric(meas[, 3]))
}
fSemilocal <- function(gg) {
    N    <- vcount(gg)
    meas <- matrix(0, nrow = N, ncol = 3)
    meas[, 1] <- ego_size(gg, order = 2, mode = 'all') - 1
    neigSum <- function(i, graph, vec) {
        neig <- igraph::neighbors(graph, v = i, mode = "all")
        return(sum(vec[neig]))
    }
    meas[, 2] <-
        vapply(seq_len(N),
               neigSum,
               c(sum = 0),
               graph = gg,
               vec = meas[, 1])
    meas[, 3] <-
        vapply(seq_len(N),
               neigSum,
               c(sum = 0),
               graph = gg,
               vec = meas[, 2])
    return(as.numeric(meas[, 3]))
}
##calculate the mean and sd of the shortest paths for each gene
calShorestPaths <- function(gg,distL = NULL,BPparam=bpparam()) {
    N    <- vcount(gg)
    getSPstat<-function(i,gg,distL){
        sp <- as.numeric(igraph::shortest.paths(gg, i,mode='all'),weights=distL)
        sp <- sp[-i]
        sp <- sp[!sp == Inf]
        res<-c(i,round(min(sp),3),round(mean(sp), 3),round(stats::sd(sp), 3))
        return(res)
    }
    measL<-bplapply(seq_len(N),getSPstat,gg,distL)
    meas<-do.call(rbind,measL)
    meas<-meas[order(meas[,1]),]
    return(meas[,-1])
}
#
# filterZeroDegree <- function( DIR, SUB ){
#
#   ft <- list()
#
#   for( d in 1:length(DIR) ){
#     tt <- read.table(sprintf("%s/random_%s_permuteDEGREE.csv",DIR[d],SUB),
#     sep="\t",header=T)
#
#     Nc <- length(tt[1,])
#     tt <- tt[,2:Nc]
#
#     tmp <- apply(tt, 2, function(x) {ifelse(x == 0, NA, 1)})
#
#     ft[[d]] <- tmp
#
#     names(ft)[d] <- sprintf("RANDOM%d",d)
#
#   }
#
#   return(ft)
#
# }
#
#
# readRandomDataFiles <- function( DIR, SUB, CENT, C, FILTER , COLN ){
#
#   ranDF <- list()
#
#   for( d in 1:length(DIR) ){
#
#     if( C == 6 ){
#       tt <- read.table(sprintf("%s/random_MEAN_%s_permute%s.csv",DIR[d],
#     SUB,CENT[C]),sep="\t",header=T)
#     } else {
#       tt <- read.table(sprintf("%s/random_%s_permute%s.csv",DIR[d],SUB,
#     CENT[C]),sep="\t",header=T)
#     }
#
#
#     Nc <- length(tt[1,])
#
#     tt <- tt[,2:Nc]
#
#     if( !is.null(FILTER) ){
#       tt <- tt*FILTER[[d]]
#     }
#
#     ##test
#     tt  <- tt[,1]
#
#     tmp <- as.numeric(unlist(tt))
#
#     tmp <- tmp[!is.na(tmp)]
#
#     tmp <- data.frame(x=as.numeric(tmp),group=COLN[d])
#
#     ranDF[[d]] <- tmp
#
#     names(ranDF)[d] <- sprintf("RANDOM%d",d)
#
#     rm(tt,Nc,tmp)
#
#   }
#
#   return(ranDF)
#
# }
#
#
formatLogLogPlot <- function(X, GROUP) {
    X <- as.vector(X)
    mm <- stats::ecdf(X)
    df <- data.frame(x = sort(X),
                     y = 1 - mm(sort(X)),
                     group = GROUP)
    return(df)
}
filterLogLog <- function(df, xMAX, yMIN) {
    if (!is.null(xMAX)) {
        df <- df[df$x <= xMAX,]
    }
    if (!is.null(yMIN)) {
        df <- df[df$y >= yMIN,]
    }
    indx <- which(df$y == 0)
    if (length(indx) != 0) {
        df <- df[-indx, ]
    }
    return(df)
}
##Calculate the Median absolute difference
MAD <- function(X) {
    X <- as.numeric(X)
    Xmd <- stats::median(X)
    MAD <- stats::median(abs(X - Xmd))
    return(MAD)
}

#' Calculate centrality measures for graph nodes.
#'
#' @details
#' The edge attribute \code{weights} treated differently by different functions
#' calculating centrality measures. For example,
#' \code{\link[igraph]{betweenness}} use \code{weights} as an edge length,
#' while in \code{\link[igraph]{page.rank}} "an edge with a larger weight is
#' more likely to be selected by the surfer", which infer the opposite meaning.
#' Taking into account that all methods in \code{\link{getClustering}} treat
#' edge \code{weights} in the same way as \code{\link[igraph]{page.rank}}, we
#' calculate the \code{distance}=1/\code{weights} as edge weights for
#' \code{BET}, \code{dBET}, \code{mnSP}, and \code{sdSP} values. So we treat
#' \code{weights} in the package consistently as the strength and closeness of
#' vertices, rather the distance between them.
#'
#' This function is able to use parallel environment, when avaliable. Parameter
#' \code{BPparam} specifies an optional \code{\link[BiocParallel]{BiocParallelParam}}
#'        instance defining the parallel back-end to be used during evaluation
#'        to be used by \code{\link[BiocParallel]{bplapply}} function.
#'
#' @param gg igraph object
#' @param weights Possibly a numeric vector giving edge weights. If this is
#'        NULL and the graph has a weight edge attribute, then the attribute
#'        is used. If this is NA then no weights are used (even if the graph
#'        has a weight attribute).
#' @param BPparam An optional \code{\link[BiocParallel]{BiocParallelParam}}
#'        instance defining the parallel back-end to be used during evaluation
#'        to be used by \code{\link[BiocParallel]{bplapply}} function.
#'
#' @return data.frame with following columns:
#' * ID   - vertex ID
#' * DEG  - degree
#' * iDEG - in-degree (directed graph only)
#' * oDEG - out-degree (directed graph only)
#' * BET  - betweenness for undirected graph
#' * dBET - betweenness when directionality is taken into account
#'          (directed graph only)
#' * CC   - clustering coefficient
#' * SL   - semilocal centrality
#' * mnSP - mean shortest path
#' * PR   - page rank  for undirected graph
#' * dPR  - page rank when directionality is taken into account
#'          (directed graph only)
#' * sdSP - standard deviation of the shortest path
#' @export
#' @import BiocParallel
#'
#' @family {Parallel Functions}
#'
#' @examples
#' library(synaptome.db)
#' cid<-match('Presynaptic',getCompartments()$Name)
#' t<-getAllGenes4Compartment(cid)
#' gg<-buildFromSynaptomeByEntrez(t$HumanEntrez)
#' system.time(m<-getCentralityMatrix(gg))
#' system.time(m0<-getCentralityMatrix(gg,BPparam=BiocParallel::SerialParam()))
#' identical(m,m0)
getCentralityMatrix <- function(gg,weights = NULL,BPparam=bpparam()) {
    tmp <- makeCentralityMatrix(gg,weights = weights,BPparam=BPparam)
    return(tmp)
}

getWDistance<-function(weights){
    if(any(weights==0)){
        delta <- 1e-2/max(weights)
        distL <- 1/(delta+weights)
    }else{
        distL <- 1/weights
    }
    return(distL)
}

makeCentralityMatrix <- function(gg,weights = NULL,BPparam=bpparam()) {
    if (is.null(weights) && "weight" %in% edge_attr_names(gg)) {
        weights  <- E(gg)$weight
        distL<-getWDistance(weights)
    }
    if (!is.null(weights) && any(!is.na(weights))) {
        weights <- as.numeric(weights)
        distL<-getWDistance(weights)
    }
    else {
        distL <- NA
        weights <- NA
    }
    N  <- vcount(gg)
    if(is.directed(gg)){
        CN  <- c("ID", "DEG", "iDEG", "oDEG", "BET", "dBET", "CC", "SL",
                 "PR", "dPR")
    }else{
        CN  <- c("ID", "DEG", "BET", "CC", "SL", "PR")
    }
    calcCentVec<-function(type=c("ID", "DEG", "iDEG", "oDEG",
                                "BET", "dBET", "CC", "SL",
                                "PR", "dPR"),gg,weights,distL){
        type <- match.arg(type)
        cl <- try(switch(
            type,
            ID = igraph::V(gg)$name,
            DEG = igraph::degree(graph = gg,mode = 'total'),
            iDEG = igraph::degree(graph = gg,mode = 'in'),
            oDEG = igraph::degree(graph = gg,mode = 'out'),
            BET = igraph::betweenness(gg,directed = FALSE,weights = distL),
            dBET = igraph::betweenness(gg,directed = TRUE,weights = distL),
            PR = igraph::page.rank(
                graph = gg,
                vids = igraph::V(gg),
                directed = FALSE,
                weights = weights,
                options = igraph::igraph.arpack.default
            )$vector,
            dPR = igraph::page.rank(
                graph = gg,
                vids = igraph::V(gg),
                directed = TRUE,
                weights = weights,
                options = igraph::igraph.arpack.default
            )$vector,
            CC = igraph::transitivity(gg, "local"),
            SL = BioNAR:::fSemilocal(gg)
        ),silent = TRUE)
        if (inherits(cl, "try-error")) {
            warning('Centrality of type "',
                    type,
                    '" failed. NULL is returned')
            cl<-NULL
        }
        l<-list()
        l[[type]]<-cl
        return(l)

    }
    toStop <- FALSE
    if(!bpisup(BPparam)){
        bpstart(BPparam)
        toStop <- TRUE
    }
    centL<-bplapply(CN,calcCentVec,gg=gg,weights=weights,distL=distL,
                    BPPARAM = BPparam)
    tmp<-as.data.frame(centL)
    tmp<-tmp[,CN]
    res <- calShorestPaths(gg,distL = distL,BPparam=BPparam)
    if(toStop) bpstop(BPparam)
    tmp$mnSP  <- res[, 2]
    tmp$sdSP  <- res[, 3]
    return(tmp)
}

#' Convert character matrix to data.frame
#'
#' All columns that names are not in \code{keep} list will be converted
#' to numbers by \code{\link{as.numeric}}.
#'
#' @param m input matrix, for example from \code{\link{vertex_attr_names}}
#' @param keep vector of column names to keep as characters
#'
#' @return data.frame with numerical columns when needed
#' @noRd
makeDataFrame <- function(m, keep = c('ID')) {
    keep.idx <- which(colnames(m) %in% keep)
    numcol <- dim(m)[2]
    df <- as.data.frame(m)
    for (i in seq_len(numcol)[-keep.idx]) {
        df[, i] <- suppressWarnings(as.numeric(df[, i]))
    }
    return(df)
}

#' Add attributes to the vertex.
#'
#' This function suits more for updating calculated vertex properties rathe
#' than node annotation. For the later case use \code{\link{annotateVertex}}.
#'
#' Unlike \code{\link{annotateVertex}}, which is able to collapse multiple
#' annotation terms, this function assume that vertex ID values are unique
#' in the \code{m} matrix.
#'
#' @param gg igraph object
#' @param m matrix of values to be applied as vertex attributes.
#'     matrix should contains column "ID" to map value to the vertex.
#'
#' @return modified igraph object
#'
#' @seealso annotateVertex
#' @export
#' @examples
#' g1 <- make_star(10, mode="undirected")
#' V(g1)$name <- letters[1:10]
#' m<-cbind(ID=letters[1:10],capital=LETTERS[1:10])
#' g1<-BioNAR::applpMatrixToGraph(g1,m)
#' V(g1)$capital
applpMatrixToGraph <- function(gg, m) {
    ggm <- gg
    measures <- colnames(m)
    id.col <- which(measures == 'ID')
    if(any(table(m[,id.col])>1)){
        stop("Vertex IDs suppose to be unique.")
    }
    meas.col <- which(measures != 'ID')
    for (i in meas.col) {
        #remove previous annotation of that name
        #check does it really needed
        ggm <- removeVertexTerm(ggm, measures[i])
        idx <- match(V(gg)$name, m[, id.col])
        naid <- which(is.na(idx))
        if (length(naid) == 0) {
            ggm <- set.vertex.attribute(
                graph = ggm,
                name = measures[i],
                index = V(ggm),
                value = m[idx, i]
            )
        } else{
            gindex <- which(!is.na(idx))
            ggm <- set.vertex.attribute(
                graph = ggm,
                name = measures[i],
                index = gindex,
                value = m[idx[gindex], i]
            )
        }
    }
    return(ggm)
}
#' Calculate the vertex centrality measures
#'
#' @description
#' Calculate the vertex centrality measures (degree, betweenness, closeness,
#' semi-local, etc....) for each graph vertex and store each result as
#' new vertex attribute in the graph.
#'
#' @details
#' A wrapper function that first calls \code{\link{getCentralityMatrix}}, to
#' calculate all vertex centrality measures, and then
#' \code{\link{applpMatrixToGraph}} to store each centrality result as a new
#' vertex attribute in the graph. The use of \code{weights} explained in
#' details in \code{\link{getCentralityMatrix}}.
#'
#' This function is able to use parallel environment, when avaliable. Parameter
#' \code{BPparam} specifies an optional \code{\link[BiocParallel]{BiocParallelParam}}
#'        instance defining the parallel back-end to be used during evaluation
#'        to be used by \code{\link[BiocParallel]{bplapply}} function.
#'
#' @param gg igraph object
#' @param weights Possibly a numeric vector giving edge weights. If this is
#'        NULL and the graph has a weight edge attribute, then the attribute
#'        is used. If this is NA then no weights are used (even if the graph
#'        has a weight attribute).
#' @param BPparam An optional \code{\link[BiocParallel]{BiocParallelParam}}
#'        instance defining the parallel back-end to be used during evaluation
#'        to be used by \code{\link[BiocParallel]{bplapply}} function.
#'
#' @return modified igraph object
#' @export
#' @seealso [getCentralityMatrix()]
#' @family {Parallel Functions}
#'
#' @examples
#' data(karate,package='igraphdata')
#' ggm<-calcCentrality(karate)
#' V(ggm)$DEG
calcCentrality <- function(gg,weights = NULL,BPparam=bpparam()) {
    m <- makeCentralityMatrix(gg,weights = weights,BPparam=BPparam)
    ggm <- applpMatrixToGraph(gg, m)
    return(ggm)
}

#get centrality measures for random graph
#' Centrality measures for random graphs induced by input one
#'
#' Generate a random graph that mimics the properties of the input graph and
#' calls \code{\link{getCentralityMatrix}} to calculate all available vertex
#' centrality measures. There are four different types of random graph to
#' generate
#'
#'
#' @param gg template graph to mimic
#' @param type type of random graph to generate:
#' * gnp -- G(n,p) Erdos-Renyi model (\code{\link[igraph]{sample_gnp}})
#' * pa --  Barabasi-Albert model (\code{\link[igraph]{sample_pa}})
#' * cgnp -- new random graph from a given graph by randomly a
#' dding/removing edges (\code{\link[igraph]{sample_correlated_gnp}})
#' * rw -- new random graph from a given graph by rewiring 25% of
#' edges preserving the degree distribution
#' \code{\link[igraph]{sample_gnp}},
#' \code{\link[igraph]{sample_correlated_gnp}}, and
#' \code{\link[igraph]{sample_pa}}
#' @param power optional argument of the power of the preferential attachment
#' to be passed to \code{\link[igraph]{sample_pa}}. If \code{power} is
#' \code{NULL} the power of the preferential attachment will be estimated
#' from \code{\link{fitDegree}} function.
#' @param weights Possibly a numeric vector giving edge weights. If this is
#'        NULL and the graph has a weight edge attribute, then the attribute
#'        is used. If this is NA then no weights are used (even if the graph
#'        has a weight attribute).
#' @param BPparam An optional \code{\link[BiocParallel]{BiocParallelParam}}
#'        instance defining the parallel back-end to be used during evaluation
#'        to be used by \code{\link{getCentralityMatrix} function.
#' @param ... other parameters passed to random graph generation functions
#'
#' @return matrix of random graph vertices centrality measure.
#' @export
#' @seealso [getCentralityMatrix()] for explanation of the use of \code{weights}
#'         and \code{BPparam}.
#' @family {Parallel Functions}
#'
#' @examples
#' data(karate,package='igraphdata')
#' m<-getRandomGraphCentrality(karate,'pa',threads=1)
#' # to avoid repetitive costy computation of PowerLaw fit
#' # power parameter could be send explicitly:
#' pFit <- fitDegree( as.vector(igraph::degree(graph=karate)),
#' Nsim=10, plot=FALSE,threads=1)
#' pwr <- slot(pFit,'alpha')
#' m<-getRandomGraphCentrality(karate,'pa',power=pwr)
#' lpa<-lapply(1:5,getRandomGraphCentrality,gg=karate,type='pa',
#'             power=pwr,weights = NULL)
getRandomGraphCentrality <- function(gg,
                                     type = c('gnp', 'pa', 'cgnp', 'rw'),
                                     power = NULL,weights = NULL,
                                     BPparam=bpparam(),
                                     ...) {
    op <- options(warn = -1)
    type <- match.arg(type)
    nv <- vcount(gg)
    ne <- ecount(gg)
    prob <- (2 * ne) / (nv * (nv - 1))
    rg <- switch (
        type,
        gnp = getGNP(gg, ...),
        pa = getPA(gg, pwr = power, ...),
        cgnp = sample_correlated_gnp(gg, corr = 0.75, ...),
        rw = rewire(gg, keeping_degseq(niter = 0.25 * ne))
    )
    V(rg)$name <- V(gg)$name
    m <- makeCentralityMatrix(rg,weights = weights,BPparam=BPparam)
    options(op)
    return(m)
}

#' Generate random graph from reference
#'
#' Function generates random G(n,p) Erdos-Renyi graph
#' (\code{\link[igraph]{sample_gnp}}) with the same number of vertices and
#' edges as in in the reference graph \code{gg}.
#'
#' @param gg reference graph
#' @param ... additional arguments to be passed to
#'          \code{\link[igraph]{sample_gnp}}
#'
#' @return new instance of the random graph.
#' @export
#'
#' @examples
#' data(karate,package='igraphdata')
#' vcount(karate)
#' ecount(karate)
#' rg<- getGNP(karate)
#' vcount(rg)
#' ecount(rg)
getGNP <- function(gg, ...) {
    nv <- vcount(gg)
    ne <- ecount(gg)
    prob <- (2 * ne) / (nv * (nv - 1))
    g <- sample_gnp(nv, p = prob, ...)
    return(g)
}

#' Generate random graph from reference
#'
#' The function generates random Barabasi-Albert graph
#' (\code{\link[igraph]{sample_pa}}) with the same vertex number as in the
#' reference graph \code{gg} and the power specified by parameter \code{pwr}.
#' If pwr is missing, we are trying to estimate pwr from the reference
#' graph \code{gg}.
#'
#' @param gg reference graph
#' @param pwr the power parameter for the \code{\link[igraph]{sample_pa}}
#' @param ... additional parameters to be passed to the
#'          \code{\link[igraph]{sample_pa}}
#'
#' @return new instance of the random graph.
#' @export
#'
#' @examples
#' data(karate,package='igraphdata')
#' vcount(karate)
#' ecount(karate)
#' rg<- getPA(karate,pwr=1.25)
#' vcount(rg)
#' ecount(rg)
getPA <- function(gg, pwr, ...) {
    nv <- vcount(gg)
    args <- list(...)
    if (is.null(pwr)) {
        pFit <- do.call(function(...) {
            fitDegree(as.vector(igraph::degree(graph = gg)),
                      Nsim = 100,
                      plot = FALSE,
                      ...)
        }
        , args[names(args) %in% formalArgs(fitDegree)])
        pwr <- pFit@alpha
    }
    g <- do.call(function(...) {
        sample_pa(nv, power = pwr, directed = FALSE, ...)
    }, args[names(args) %in% formalArgs(sample_pa)])
    return(g)
}
#' Convert centrality matrix into ECDF
#'
#' @param m centrality matrix from \code{\link{getCentralityMatrix}}
#' invocation.
#'
#' @return list of several ecdf objects, corresponding to values in
#' centrality matrix from \code{\link{getCentralityMatrix}} invocation.
#'
#' @seealso getCentralityMatrix
#' @export
#' @examples
#' library(synaptome.db)
#' cid<-match('Presynaptic',getCompartments()$Name)
#' t<-getAllGenes4Compartment(cid)
#' gg<-buildFromSynaptomeByEntrez(t$HumanEntrez)
#' m<-getCentralityMatrix(gg)
#' ecdfL<-getGraphCentralityECDF(m)
getGraphCentralityECDF <- function(m) {
    idx <- which(colnames(m) != 'ID')
    l <- list()
    for (i in 2:8) {
        n <- colnames(m)[i]
        l[[n]] <- stats::ecdf(m[, i])
    }
    return(l)
}

#' Extracts particular measure from matrix and convert for distance
#' calculation by \code{\link{calcCentralityInternalDistances}} and
#' \code{\link{calcCentralityExternalDistances}} functions.
#'
#' @param m matrix of centrality measures as returned by getCentralityMatrix
#' @param nm name of the measure from m
#' @param keepOrder if FALSE valuess will be sorted
#'
#' @return vector of length \code{dim(m)[1]}
#' @seealso calcCentralityInternalDistances
#' @seealso calcCentralityExternalDistances
#' @noRd
getCM <- function(m, nm, keepOrder) {
    v <- as.numeric(m[, which(colnames(m) == nm)])
    if (keepOrder) {
        return(v)
    } else{
        return(sort(v, decreasing = FALSE, na.last = TRUE))
    }
}
#' Function calculates a set of distance metrics between each vertex pair  given
#' a list of vertex centrality matrices
#'
#' @param l list of matrices, for example centrality obtained by invocation
#'         \code{\link{getRandomGraphCentrality}}
#' @param keepOrder if FALSE values will be sorted before distance calculations
#' @param dist methods available from \code{\link{dist}} function
#'
#' @return matrix with seven columns containing distances between all pairs of
#'         \code{l} elements.
#'
#' @seealso getRandomGraphCentrality
#' @seealso getCentralityMatrix
#' @seealso calcCentralityExternalDistances
#' @export
#' @examples
#' data(karate,package='igraphdata')
#' m<-getCentralityMatrix(karate)
#' gnp<-list()
#' for(i in 1:10){
#'     gnp[[i]]<-getRandomGraphCentrality(karate,type = 'gnp')
#' }
#' gnpIDist<-calcCentralityInternalDistances(gnp)
#' summary(gnpIDist)
calcCentralityInternalDistances <-
    function(l, keepOrder = FALSE, dist = 'euclidean') {
        CN  <- c("ID", "DEG", "BET", "CC", "SL", "mnSP", "PR", "sdSP")
        resl <- list()
        for (i in 2:length(CN)) {
            nm <- CN[i]
            res <- do.call(cbind, lapply(l, getCM, nm = nm,
                                         keepOrder = keepOrder))
            if (is.matrix(res)) {
                resl[[nm]] <- as.vector(dist(t(res), method = dist))
            }
        }
        resm <- do.call(cbind, resl)
        return(resm)
    }
#' Function to calculate a distance matrix between a list of permuted vertex
#' centrality matrices and a unperturbed reference matrix.
#'
#' @param m reference matrix, for example centrality obtained by invocation
#'         \code{\link{getCentralityMatrix}}
#' @param l list of permuted matrix, for example centrality obtained by
#'         invocation \code{\link{getRandomGraphCentrality}}
#' @param keepOrder if FALSE valuess will be sorted
#' @param dist methods available from dist function
#'
#' @return matrix with seven columns  containing distances between each element
#'         of \code{l} and reference matrix \code{m}
#'
#' @seealso getRandomGraphCentrality
#' @seealso getCentralityMatrix
#' @seealso calcCentralityInternalDistances
#' @export
#' @examples
#' data(karate,package='igraphdata')
#' m<-getCentralityMatrix(karate)
#' gnp<-list()
#' for(i in 1:10){
#'     gnp[[i]]<-getRandomGraphCentrality(karate,type = 'gnp')
#' }
#' gnpEDist<-calcCentralityExternalDistances(m,gnp)
#' summary(gnpEDist)
calcCentralityExternalDistances <-
    function(m,
             l,
             keepOrder = FALSE,
             dist = 'euclidean') {
        CN  <- c("ID", "DEG", "BET", "CC", "SL", "mnSP", "PR", "sdSP")
        resl <- list()
        for (i in 2:length(CN)) {
            nm <- CN[i]
            rm <- getCM(m, nm = nm, keepOrder = keepOrder)
            res <- do.call(cbind, lapply(l, getCM, nm = nm,
                                         keepOrder = keepOrder))
            if (is.matrix(res)) {
                cmm <- cbind(rm, res)
                cmd <- as.matrix(dist(t(cmm), method = dist))
                resl[[nm]] <- as.vector(cmd[-1, 1])
            }
        }
        resm <- do.call(cbind, resl)
        return(resm)
    }
#' Compare distance distributions of internal and external distances
#'
#' Function to compare two distance distributions using the Kolmogorov-Smirnov
#' test. Where the first distance distribution is generated internally and
#' calculates the distance between random graph centralities. The second
#' distance distribution is generated externally, and measures the distance
#' between random and the original graph centralities.
#'
#' @param dmi distribution of internal distances between random graph
#' centralities
#' @param dme distribution of external distances between random and
#' original graph centralities
#'
#' @return list of lists for each centrality value in the input matrix three
#' element list is created where \code{ks} contains Kolmogorov-Smirnov test
#' result from class \code{ks.test}; \code{pval} contains Kolmogorov-Smirnov
#' test pvalue;
#' and \code{dt} contains input distribution.
#' @export
#' @seealso ks.test
#' @examples
#' data(karate,package='igraphdata')
#' m<-getCentralityMatrix(karate)
#' gnp<-list()
#' for(i in 1:10){
#'     gnp[[i]]<-getRandomGraphCentrality(karate,type = 'gnp')
#' }
#' gnpIDist<-calcCentralityInternalDistances(gnp)
#' gnpEDist<-calcCentralityExternalDistances(m,gnp)
#'
#' simSig<-evalCentralitySignificance(gnpIDist,gnpEDist)
#' sapply(simSig,function(.x).x$ks$p.value)
evalCentralitySignificance <- function(dmi, dme) {
    nmi <- colnames(dmi)
    nme <- colnames(dme)
    nms <- intersect(nmi, nme)
    l <- list()
    for (nm in nms) {
        mi <- dmi[, colnames(dmi) == nm]
        me <- dme[, colnames(dme) == nm]
        ks <- stats::ks.test(mi, me)
        l[[nm]] <- list(
            pval = ks$p.value,
            ks = ks,
            dt = data.frame(val = c(mi, me),
                            cl = factor(c(
                                rep('perm', length(mi)),
                                rep('graph', length(me))
                            )))
        )
    }
    return(l)
}
