maxLSi <- function(XX, BASE = 0) {
    XX  <- as.vector(as.numeric(XX))

    XXo <- XX[XX != 0]
    if (BASE == 2) {
        return (-sum(XXo * log2(XXo)))
    }

    if (BASE == 10) {
        return (-sum(XXo * log10(XXo)))
    }

    if (BASE == 0) {
        return (-sum(XXo * log(XXo)))
    }

}

#' Calculate parameters for Entropy plot and calculations.
#'
#' @param gg igroph object
#'
#' @return list with values of maxSr and SRo
#' @export
#' @import RSpectra
#'
#' @examples
#' data(karate, package='igraphdata')
#' ent <- getEntropyRate(karate)
getEntropyRate <- function(gg) {
    V    <- length(V(gg))
    E    <- length(E(gg))
    ki   <- as.vector(igraph::degree(gg))
    Kbar <- mean(ki)

    #--- get adjacency matrix for graph
    A    <- get.adjacency(gg)

    #--- get leading eigenvalue and vector
    #R     <- eigen(A)
    #Rindx <- which.max(R$values)
    #gamma <- R$values[Rindx]
    #nu    <- R$vectors[,Rindx]
    R <- RSpectra::eigs(A, 1)
    gamma <- R$values[1]
    nu    <- R$vectors[, 1]


    #--- calculate max entropy rate, maxSr
    Pij   <- (A * nu) / (gamma * nu)
    Pi    <- ki / (2 * E)
    maxSr <- sum(as.vector(Pi * apply(Pij, 1, maxLSi, BASE = 0)))


    #--- calculate initial configuration
    Norm <- as.numeric(V * Kbar) #as.numeric(2*E)
    SRo  <- as.numeric(1 / Norm) * sum(ki * log(ki))

    return(list(maxSr = maxSr, SRo = SRo))
}

#' Calculates perturbation entropy and save it as attribute on the graph.
#'
#' @param gg igraph object
#' @param maxSr maxSr value, if NULL \code{getEntropyRate} will be called.
#' @param exVal expression values boundaries.
#' Two columns are expected: \code{xx} and \code{lambda}. If NULL
#' default values \code{c(2,14)} and \code{c(-14,14)} will be used
#' for \code{xx}  and \code{lambda} respectively.
#'
#' @return graph with SR_UP and SR_DOWN attributes storing entropy values for
#' over- and underexpression respectively
#' @export
#'
#' @examples
#' library(synaptome.db)
#' cid<-match('Presynaptic', getCompartments()$Name)
#' t<-getAllGenes4Compartment(cid)
#' gg<-buildFromSynaptomeByEntrez(t$HumanEntrez)
#' gg<-annotateGeneNames(gg)
#' gg<- calcEntropy(gg)
calcEntropy <- function(gg, maxSr = NULL, exVal = NULL) {
    SRprime <- getEntropy(gg, maxSr = maxSr, exVal = exVal)
    SRprime <- SRprime[, -c(2, 3)]
    names(SRprime) <- c('ID', 'SR_UP', 'SR_DOWN')
    gg <- applpMatrixToGraph(gg, SRprime)
    return(gg)
}

#' Calculates perturbation entropy
#'
#' @param gg igraph object
#' @param maxSr maxSr value, if NULL \code{getEntropyRate} will be called.
#' @param exVal expression values boundaries.
#' Two columns are expected: \code{xx} and \code{lambda}.
#' If NULL default values \code{c(2,14)} and \code{c(-14,14)} will
#' be used for \code{xx} and \code{lambda} respectively.
#'
#' @return matrix with EntrexID, GeneNames, Degree, Enttropy UP and Down
#' @export
#'
#' @examples
#' library(synaptome.db)
#' cid<-match('Presynaptic', getCompartments()$Name)
#' t<-getAllGenes4Compartment(cid)
#' gg<-buildFromSynaptomeByEntrez(t$HumanEntrez)
#' gg<-annotateGeneNames(gg)
#' e<- getEntropy(gg)
getEntropy <- function(gg, maxSr = NULL, exVal = NULL) {
    if (!"GeneName" %in% vertex_attr_names(gg)) {
        V(gg)$GeneName <- V(gg)$name
    }
    V    <- vcount(gg)
    E    <- ecount(gg)
    ki   <- as.vector(igraph::degree(gg))
    Kbar <- mean(ki)
    A    <- get.adjacency(gg)
    if (is.null(maxSr)) {
        par <- getEntropyRate(gg)
        maxSr <- par$maxSr
    }
    if (is.null(exVal)) {
        xx    <- vector(length = 2)
        xx[1] <- 2  #active
        xx[2] <- 16 #inactive
        lambda    <- vector(length = 2)
        lambda[1] <- 14   #active
        lambda[2] <- -14  #inactive
        NORM      <- vector(length = 2)
        NORM[1]   <- 0
        NORM[2]   <- 0
    } else{
        if (!('xx' %in% names(exVal) && 'lambda' %in% names(exVal))) {
            stop('exVal should contain two columns: xx and lambda')
        }
        xx <- exVal$xx
        lambda <- exVal$lambda
        NORM <- rep.int(0, length(xx))
    }
    #colnames(SRprime) <- c("ENTREZ.ID","GENE.NAME","DEGREE","UP","DOWN")
    SRprime <- data.frame(
        ENTREZ.ID = V(gg)$name,
        GENE.NAME = V(gg)$GeneName,
        DEGREE = ki,
        UP = rep(NA, V),
        DOWN = rep(NA, V)
    )
    for (v in seq_len(V)) {
        GN     <- as.character(SRprime[v, 2])
        GNindx <- which(V(gg)$GeneName == GN)
        PIprime <- cbind(rep(NA, V), rep(NA, V))
        LSprime <- cbind(rep(NA, V), rep(NA, V))
        NORM <- rep.int(0, length(xx))
        for (s in seq_along(lambda)) {
            X               <- rep(xx[s], V)
            X[GNindx[1]]    <- X[GNindx[1]] + lambda[s]
            NORM[s]         <- X %*% A %*% X
        }
        Nv <- V(gg)$name[neighbors(gg, GNindx, mode = "all")]
        oo <- cbind(ki, !(V(gg)$name %in% Nv))
        for (s in seq_along(lambda)) {
            PIprime[, s] <- ifelse(oo[, 2] == 1,
                                   (1 / NORM[s] * xx[s] * xx[s] * oo[, 1]),
                                   NA)
        }
        for (s in seq_along(lambda)) {
            X   <- (xx[s])
            lam <- (lambda[s])
            DEG <- (oo[GNindx[1], 1])

            PIprime[GNindx[1], s] <- ((X + lam) * DEG * X) / NORM[s]

        }
        for (s in seq_along(lambda)) {
            PIprime[, s] <- ifelse(oo[, 2] == 0,
                                   (1 / NORM[s] * xx[s] *
                                        (xx[s] + lambda[s] +
                                             (oo[, 1] - 1) * xx[s])),
                                   PIprime[, s])
        }
        for (s in seq_along(lambda)) {
            X <- (xx[s])
            LSprime[, s] <- ifelse(oo[, 2] == 1,
                                   (-log(X) + log(X * (oo[, 1]))), NA)
        }
        Ni <- which(0 == oo[, 2])
        for (i in seq_along(Ni)) {
            DEGi <- (oo[Ni[i], 1])
            SUM  <- DEGi - 1
            for (s in seq_along(lambda)) {
                X   <- (xx[s])
                lam <- (lambda[s])
                dem <- X + lam + (DEGi - 1) * X
                pij <- X / dem
                pi1 <- (X + lam) / dem
                LSi <- pij * log(pij)
                LSi <- -SUM * LSi - pi1 * log(pi1)
                LSprime[Ni[i], s]  <- LSi
            }
        }
        SRprime[v, 4] <- sum((PIprime[, 1]) * (LSprime[, 1]))
        SRprime[v, 5] <- sum((PIprime[, 2]) * (LSprime[, 2]))
    }
    SRprime[, 4] <- (SRprime[, 4]) / maxSr
    SRprime[, 5] <- (SRprime[, 5]) / maxSr
    colnames(SRprime) <-
        c("ENTREZ.ID", "GENE.NAME", "DEGREE", "UP", "DOWN")
    return(SRprime)
}

getEntropyOverExpressed <- function(SRprime, perc = 1) {
    #--- Bottom 1% UP, i.e. OVER-EXPRESSED
    V <- dim(SRprime)[1]
    XX  <- as.numeric(SRprime[, 4])
    MIN <- min(XX)
    MAX <- max(XX)
    XX2 <- (XX - MIN) / (MAX - MIN)
    oo2 <- cbind(SRprime[, 2], XX2)
    oo2 <- oo2[order(as.numeric(oo2[, 2])),]
    ii  <- floor(perc / 100 * V)
    GN  <- oo2[seq_len(ii), 1]
    DF3 <- SRprime[match(GN, SRprime[, 2]), c(1, 2, 3, 4)]
    DF3 <- cbind(DF3, rep(paste0(perc, "%"), length(GN)))
    return(DF3)
}

#' Plot entropy values
#'
#' @param SRprime results of \code{\link{getEntropy}} invocation
#' @param subTIT entropy axis label
#' @param maxSr first element of results of
#'            \code{\link{getEntropyRate}} invocation
#' @param SRo second element of results of
#'            \code{\link{getEntropyRate}} invocation
#'
#' @return ggplot2 object with diagram
#' @export
#' @import ggplot2
#'
#' @examples
#' library(synaptome.db)
#' cid<-match('Presynaptic',getCompartments()$Name)
#' t<-getAllGenes4Compartment(cid)
#' gg<-buildFromSynaptomeByEntrez(t$HumanEntrez)
#' gg<-annotateGeneNames(gg)
#' ent <- getEntropyRate(gg)
#' SRprime <- getEntropy(gg, maxSr = NULL)
#' plotEntropy(SRprime, subTIT = "Entropy", SRo = ent$SRo, maxSr = ent$maxSr)
plotEntropy <- function(SRprime,
                        subTIT = 'Entropy',
                        SRo = NULL,
                        maxSr = NULL) {
    colours <- c('lawngreen', 'firebrick2')

    V <- dim(SRprime)[1]
    DF1 <- SRprime[, c(1, 2, 3, 4)]
    names(DF1) <- c("ENTREZ.ID", "GENE.NAME", "DEGREE", "SR")
    DF1$GROUP <- "SR_UP"

    DF2 <- SRprime[, c(1, 2, 3, 5)]
    names(DF2) <- c("ENTREZ.ID", "GENE.NAME", "DEGREE", "SR")
    DF2$GROUP <- "SR_DOWN"

    DF  <- rbind(DF1, DF2)
    colnames(DF) <- c("ENTREZ.ID", "GENE.NAME", "DEGREE", "SR", "GROUP")

    DF <- as.data.frame(DF)
    DF$DEGREE <- as.numeric(DF$DEGREE)
    DF$SR <- as.numeric(DF$SR)
    DF$GROUP <- as.factor(DF$GROUP)
    DF$lDEGREE <- log(DF$DEGREE)

    gplot <-
        ggplot(data = DF, aes_string(x = 'lDEGREE', y = 'SR',
                                     colour = 'GROUP')) +
        geom_point() +
        labs(x = "log(k)", y = "SR", title = subTIT) +
        guides(
            color = guide_legend(override.aes = list(fill = NA, size = 4)),
            fill  = FALSE,
            group = FALSE,
            alpha = FALSE
        ) +
        theme(
            axis.title.x = element_text(face = "bold", size = rel(2)),
            axis.text.x = element_text(face = "bold", size = rel(2)),
            axis.title.y = element_text(face = "bold", size = rel(2)),
            axis.text.y = element_text(face = "bold", size = rel(2)),
            legend.title = element_text(face = "bold", size = rel(1.5)),
            legend.text = element_text(face = "bold", size = rel(1.5)),
            legend.position = "top",
            legend.key = element_blank()
        ) +
        theme(
            panel.grid.major = element_line(colour = "grey40", size = 0.2),
            panel.grid.minor = element_line(colour = "grey40", size = 0.1),
            panel.background = element_rect(fill = "white"),
            panel.border = element_rect(linetype = "solid", fill = NA)
        ) +
        scale_color_manual("",
                           breaks = levels(DF$GROUP),
                           values = c(colours))#+
    if (!is.null(SRo) && !is.null(maxSr)) {
        gplot <- gplot +
            geom_hline(
                yintercept = SRo / maxSr,
                colour = "black",
                size = 2,
                linetype = 2,
                show.legend = FALSE
            )
    }
    # geom_hline(yintercept=SRo/maxSr,colour="grey40",
    #            size=2,linetype=2,show.legend=F)
    return(gplot)
}
