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
calShorestPaths <- function(gg) {
    N    <- vcount(gg)
    meas <- matrix(0, nrow = N, ncol = 3)
    for (i in seq_len(N)) {
        sp <- as.numeric(igraph::shortest.paths(gg, i))
        sp <- sp[-i]
        sp <- sp[!sp == Inf]
        meas[i, 1] <- min(sp)
        meas[i, 2] <- round(mean(sp), 3)
        meas[i, 3] <- round(stats::sd(sp), 3)
    }
    return(meas)
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
#' @param gg igraph object
#'
#' @return matrix with following columns:
#' * ID - vertex ID
#' * DEG - degree
#' * BET - betweenness
#' * CC - clustering coefficient
#' * SL - semilocal centrality
#' * mnSP - mean shortest path
#' * PR - page rank
#' * sdSP - standard deviation of the shortest path
#' @export
#'
#' @examples
#' library(synaptome.db)
#' cid<-match('Presynaptic',getCompartments()$Name)
#' t<-getAllGenes4Compartment(cid)
#' gg<-buildFromSynaptomeByEntrez(t$HumanEntrez)
#' m<-getCentralityMatrix(gg)
getCentralityMatrix <- function(gg) {
    tmp <- makeDataFrame(makeCentralityMatrix(gg))
    return(tmp)
}

makeCentralityMatrix <- function(gg) {
    ID <- V(gg)$name
    N  <- length(ID)
    CN  <- c("ID", "DEG", "BET", "CC", "SL", "mnSP", "PR", "sdSP")
    tmp <- matrix(0, nrow = N, ncol = length(CN))
    colnames(tmp) <- CN
    tmp[, 1] <- ID
    tmp[, 2] <- as.vector(igraph::degree(graph = gg))
    tmp[, 3] <- as.character(round(betweenness(gg), 3))
    tmp[, 4] <- as.character(round(transitivity(gg, "local"), 3))
    sl <- fSemilocal(gg)
    tmp[, 5] <- as.character(round(sl, 3))
    res <- calShorestPaths(gg)
    tmp[, 6]  <- as.character(res[, 2])
    tmp[, 7]  <- as.character(round(as.vector(
        page.rank(
            graph = gg,
            vids = V(gg),
            directed = FALSE,
            options = igraph.arpack.default
        )$vector
    ), 6))
    tmp[, 8]  <- as.character(res[, 3])
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
        df[, i] <- as.numeric(df[, i])
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
#' Calculate the vertex centrality measures (degree, betweenness, closeness,
#' semi-local, etc....) for each graph vertex and store each result as
#' new vertex attribute in the graph.
#'
#' A wrapper function that first calls \code{\link{getCentralityMatrix}}, to
#' calculate all vertex centrality measures, and then
#' \code{\link{applpMatrixToGraph}} to store each centrality result as a new
#' vertex attribute in the graph.
#'
#'
#' @param gg igraph object
#'
#' @return modified igraph object
#' @export
#'
#' @examples
#' data(karate,package='igraphdata')
#' ggm<-calcCentrality(karate)
#' V(ggm)$DEG
calcCentrality <- function(gg) {
    m <- makeCentralityMatrix(gg)
    ggm <- applpMatrixToGraph(gg, m)
    return(ggm)
}
#get centrality measures for random graph
#' Centrality measures for random graphs induced by input one
#'
#' Generate a random graph  that mimic somehow properties of the input graph
#' and calls \code{\link{getCentralityMatrix}} to
#' calculate all available centrality measires. There are four different
#'
#' @param gg template graph to mimic
#' @param type type of random graph to generate:
#' * gnp -- G(n,p) Erdos-Renyi model (\code{\link[igraph]{sample_gnp}})
#' * pa --  Barabasi-Albert model (\code{\link[igraph]{sample_pa}})
#' * cgnp -- new random graph from a given graph by randomly a
#' dding/removing edges (\code{\link[igraph]{sample_correlated_gnp}})
#' * rw -- new random graph from a given graph by rewiring 25% of
#' edges preserving the degree distribution
#' @param ... other parameters passed to random graph generation functions
#' \code{\link[igraph]{sample_gnp}},
#' \code{\link[igraph]{sample_correlated_gnp}}, and
#' \code{\link[igraph]{sample_pa}}
#' @param power optional argument of the power of the preferential attachment
#' to be passed to \code{\link[igraph]{sample_pa}}. If \code{power} is
#' \code{NULL} the power of the preferential attachment will be estimated
#' from \code{\link{fitDegree}} function.
#'
#' @return matrix of random graph vertices centrality measure.
#' @export
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
getRandomGraphCentrality <- function(gg,
                                     type = c('gnp', 'pa', 'cgnp', 'rw'),
                                     power = NULL,
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
    m <- makeCentralityMatrix(rg)
    options(op)
    return(m)
}
getGNP <- function(gg, ...) {
    nv <- vcount(gg)
    ne <- ecount(gg)
    prob <- (2 * ne) / (nv * (nv - 1))
    g <- sample_gnp(nv, p = prob, ...)
    return(g)
}

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
#' @return list of sever ecdf objects, corresponding to values in
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
#' calculation by \code{\linc{calcCentralityInternalDistances}} and
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
