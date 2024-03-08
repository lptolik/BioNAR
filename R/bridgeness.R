#' Calculate bridginess from consensus matrix
#'
#' Bridginess takes into account a vertices shared community membership
#' together with its local neighbourhood. It was proposed in
#' Nepusz et al., 2008 <doi:10.1103/PhysRevE.77.016107>.
#'
#'
#' Function assumes clustering already been performed by the clustering
#' algorithm, and its membership values stored in vertex attributes. If
#' clustering algorithm vertex  \code{alg} attribute is not found an
#' error will be issued.
#'
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
#' karate <- make_graph("Zachary")
#' # We need vertex ID in the 'name' attribute of the vertex
#' V(karate)$name<-c(LETTERS,letters)[1:vcount(karate)]
#' gg <- calcClustering(karate, 'louvain')
#' cnmat <- makeConsensusMatrix(gg, N=10, alg = 'louvain', type = 2, mask = 10)
#' br<-getBridgeness(gg, alg = 'louvain', cnmat)
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
        ##reorder edge communities, so ed[, 3] equals current community no: 'c'
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

#' Helper function that uses \code{\link{getBridgeness}} to calculate
#' graph node bridgeness values for selected algorithm and consensus matrix
#' and save them as a graph attribute \code{BRIDGENESS.<alg>} with \code{<alg>}
#' replaced by the selected algorithm name.
#'
#' @inheritParams getBridgeness
#' @return graph with additional attributes to store Bridgeness value
#' @export
#' @seealso getBridgeness
#'
#' @examples
#' library(BioNAR)
#' karate <- make_graph("Zachary")
#' # We need vertex ID in the 'name' attribute of the vertex
#' V(karate)$name<-c(LETTERS,letters)[1:vcount(karate)]
#' set.seed(100)
#' gg <- calcClustering(karate, 'louvain')
#' cnmat <- makeConsensusMatrix(gg, N=10, alg = 'louvain', type = 2, mask = 10)
#' gg<-calcBridgeness(gg, alg = 'louvain', cnmat)
#' hist(V(gg)$BRIDGENESS.louvain)
calcBridgeness <- function(gg, alg, conmat) {
    br<-getBridgeness(gg, alg = alg, conmat)
    agg<-applpMatrixToGraph(gg,br[,grep('(ID|BRDIDGENESS.+)',names(br))])
    return(agg)
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

#' Plot Bridgeness values
#'
#' Semi-local centrality measure (Chen et al., 2011)
#' lies between 0 and 1 indicating whether protein is important globally or
#' locally. By plotting Bridgeness against semi-local centrality we can
#' categorises the influence each protein found in our network has on the
#' overall network structure:
#' * Region 1, proteins having a 'global' rather than 'local' influence in
#' the network (also been called bottle-neck bridges, connector or kinless
#' hubs (0<Sl<0.5; 0.5<Br<1).
#' * Region 2, proteins having 'global' and 'local' influence (0.5<Sl<1,
#' 0.5<Br<1).
#' * Region 3, proteins centred within the community they belong to, but also
#' communicating with a few other specific communities (0<Sl<0.5; 0.1<Br<0.5).
#' * Region 4, proteins with 'local' impact , primarily within one or two
#' communities (local or party hubs, 0.5<Sl<1, 0<Br<0.5).
#'
#' @param gg igraph object with bridgenes values stored as attributes,
#' after call to \code{\link{calcBridgeness}}
#' @param alg clustering algorithm that was used to calculate bridgeness values
#' @param VIPs list of 'specical' genes to be marked on the plot
#' @param Xatt name of the attribute that stores values to be used as X-axis
#' values. By default \code{SL} for semi-local centrality
#' @param Xlab label for the X-axis
#' @param Ylab label for the Y-axis
#' @param bsize point size for genes
#' @param spsize point size for 'specical' genes
#' @param MainDivSize size of the line for the region separation lines
#' @param xmin low limit for X-axis
#' @param xmax upper limit for X-axis
#' @param ymin low limit for Y-axis
#' @param ymax upper limit for Y-axis
#' @param baseColor basic color for genes
#' @param SPColor colour highlighting any 'specical' genes
#'
#' @return \code{\link[ggplot2]{ggplot}} object with plot
#' @export
#' @importFrom igraph get.vertex.attribute
#' @importFrom ggrepel geom_label_repel
#' @import ggplot2
#'
#' @examples
#' karate <- make_graph("Zachary")
#' # We need vertex ID in the 'name' attribute of the vertex
#' V(karate)$name<-c(LETTERS,letters)[1:vcount(karate)]
#' set.seed(100)
#' gg <- calcClustering(karate, 'louvain')
#' gg <- calcCentrality(gg)
#' cnmat <- makeConsensusMatrix(gg, N=10, alg = 'louvain', type = 2, mask = 10)
#' gg<-calcBridgeness(gg, alg = 'louvain', cnmat)
#' plotBridgeness(gg,alg = 'louvain',VIPs=c("Mr Hi","John A"))
plotBridgeness<-function(gg,alg,VIPs,
                         Xatt='SL',
                         Xlab = "Semilocal Centrality (SL)",
                         Ylab = "Bridgeness (B)",
                         bsize = 3,
                         spsize =7,
                         MainDivSize = 0.8,
                         xmin = 0,
                         xmax = 1,
                         ymin = 0,
                         ymax = 1,
                         baseColor="royalblue2",
                         SPColor="royalblue2"){
    #VIPs=c('8495','22999','8927','8573','26059','8497','27445','8499')
    #VIPs=c('81876','10890','51552','5874','5862','11021','54734','5865','5864',
    #        '9522','192683','10067','10396','9296','527','9114','537','535',
    #       '528','51382','534','51606','523','80331','114569','127262','57084',
    #        '57030','388662','6853','6854','8224','9900','9899','9145','9143',
    #        '6855','132204','6857','127833','6861','529','526','140679','7781',
    #        '81615','6844','6843')
    indx   <- match(V(gg)$name,VIPs)
    group <- ifelse( is.na(indx), 0,1)
    # MainDivSize <- 0.8
    # xmin        <- 0
    # xmax        <- 1
    # ymin        <- 0
    # ymax        <- 1
    # Xlab <- "Semilocal Centrality (SL)"
    # Ylab <- "Bridgeness (B)"
    X    <- as.numeric(get.vertex.attribute(gg,Xatt,V(gg)))
    if(length(X)==0){
        stop('Graph vertices have no numerical attribute "',Xatt,'"\n')
    }
    X    <- scale(X)
    Y   <- as.numeric(get.vertex.attribute(gg,
                                           sprintf("BRIDGENESS.%s", alg),
                                           V(gg)))
    if(length(Y)==0){
        stop('Graph vertices have no numerical attribute "',
             sprintf("BRIDGENESS.%s", alg),'"\n',
             "Check that you've calculated bridginess.\n")
    }
    if('GeneName' %in% vertex_attr_names(gg)){
        lbls <- ifelse(!is.na(indx),V(gg)$GeneName,"")
        name <- V(gg)$GeneName
        dt<-data.frame(X=X,Y=Y,vips=group,entres=V(gg)$name,name=name)
    }else{
        lbls <- ifelse(!is.na(indx),V(gg)$name,"")
        name <- V(gg)$name
        dt<-data.frame(X=X,Y=Y,vips=group,entres=V(gg)$name,name=name)
    }
    dt_vips<-dt[dt$vips==1,]
    dt_res<-dt[dt$vips==0,]
    ##--- baseColor of genes
    #baseColor="royalblue2"

    ##--- SPcolor, colour highlighting any 'specical' genes
    #SPColor="royalblue2"

    ##--- PSDColor, colour of core PSD95 interactor genes
    #PSDColor="magenta"

    g<-ggplot(dt,aes(x=X,y=Y,label=name))+#geom_point()+
        geom_point(data=dt_vips,
                   aes(x=X,y=Y),colour=baseColor,size = spsize,
                   shape=15,show.legend=FALSE)+
        geom_point(data=dt_res,
                   aes(x=X,y=Y, alpha=(X*Y)), size = bsize,
                   shape=16,show.legend=FALSE)+
        geom_label_repel(aes(label=as.vector(lbls)),
                         fontface='bold',color='black',
                         fill='white',box.padding=0.1,
                         point.padding=NA,label.padding=0.15,
                         segment.color='black',
                         force=1,size=rel(3.8),show.legend=FALSE,
                         max.overlaps=300)+
        labs(x=Xlab,y=Ylab,title=sprintf("%s",alg))+
        scale_x_continuous(expand = c(0, 0), limits = c(xmin, xmax)) +
        scale_y_continuous(expand = c(0, 0), limits = c(ymin, ymax))+
        theme(
            axis.title.x=element_text(face="bold",size=rel(2.5)),
            axis.title.y=element_text(face="bold",size=rel(2.5)),
            legend.title=element_text(face="bold",size=rel(1.5)),
            legend.text=element_text(face="bold",size=rel(1.5)),
            legend.key=element_blank())+
        theme(panel.grid.major = element_line(colour="grey40",linewidth=0.2),
              panel.grid.minor = element_line(colour="grey40",linewidth=0.1),
              panel.background = element_rect(fill="white"),
              panel.border = element_rect(linetype="solid",fill=NA))+
        geom_vline(xintercept=0.5,colour="grey40",linewidth=MainDivSize,
                   linetype=2,show.legend=FALSE)+
        geom_hline(yintercept=0.5,colour="grey40",linewidth=MainDivSize,
                   linetype=2,show.legend=FALSE)
    return(g)
}
