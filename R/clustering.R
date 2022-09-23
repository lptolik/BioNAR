
#' Calculate cluster memberships for the graph.
#' 
#' Calculate clustering membership for ten algoritms defined in 
#' \code{\link{getClustering}}
#'
#' @param gg igraph object to cluster
#' @param alg algorithm name
#'
#' @return data.frame with columns \code{names} and \code{membership}
#' @export
#'
#' @seealso getClustering
#' @examples
#' data(karate,package='igraphdata')
#' m<-calcMembership(karate,'lec')
#' head(m)
calcMembership<-function(gg,
                            alg=c('lec','wt','fc','infomap',
                                'louvain','sgG1','sgG2','sgG5','spectral')){
    ids <- V(gg)$name
    cl<-getClustering(gg,alg)
    if(!is.null(cl)){
    cc       <- data.frame(names=cl$names,membership=cl$membership)
    }else{
    cc <- data.frame(names='names',membership=0)[FALSE,]
    }
    return(cc)
}

#' Calculate memberships for all clustering algorithms and store them on the
#' graph vertices.
#'
#' Function apply \link{\code{calcClusttering}} for each name in standard
#' list of available algorithms. In the case when clustering could not been
#' calculated warning will be issued and no attributes added to the graph.
#'
#' @param gg graph for analysis
#'
#' @return new graph object with all membership results stored as a vertex
#'         attribute.
#' @export
#' @seealso calcClusttering
#'
#' @examples
#' g1 <- make_star(10, mode="undirected")
#' V(g1)$name <- letters[1:10]
#' g1<-calcAllClustering(g1)
#' clusteringSummary(g1)
calcAllClustering<-function(gg){
    ids <- V(gg)$name
    cnames<-c('ID','lec','wt','fc','infomap',
            'louvain','sgG1','sgG2','sgG5','spectral')
    l<-list()
    l[[cnames[1]]]<-ids
    for(ai in 2:length(cnames)){
    an<-cnames[ai]
    cm<-calcMembership(gg,an)
    if(dim(cm)[1]>0){
        l[[an]]<-as.character(cm$membership)
        mod<-modularity(gg,cm$membership)
        gg<-set.graph.attribute(gg,an,mod)
    }
    }
    m<-do.call(cbind,l)
    ggm<-applpMatrixToGraph(gg,m)
    return(ggm)
}

#' Calculate memberships for particular clustering algorithms and store 
#' them on the graph vertices.
#' 
#' Results of clustering algorithm application to the same graph could differ
#' between runs due to use of stochastic algorithm. To allow reproducible 
#' downstream analysis clustering results are stored as graph attributes. This
#' function call \code{\link{getClustering}} and store membership as vertex 
#' attribute and modularity as graph attribute with attribute name defined by
#' \code{alg} value.
#'
#' NOTE: \code{\link{getClustering}} verifies algorithm names with 
#' \code{\link[base]{match.arg}} so correct membership will be calculated, but
#' name of the attribute is taken from \code{alg} argument, so it is possible 
#' that vertex attribute name won't exactly match name of the algorithm from
#' \code{link{getClustering}}.
#'
#' @param gg igraph object to cluster
#' @param alg algorithm to apply
#'
#' @seealso getClustering
#' 
#' @return modified igraph object with calculated membership stored as a vertex
#'         attribute and modularity as a graph attribute
#' @export
#'
#' @examples
#' data(karate,package='igraphdata')
#' g<-calcClustering(karate,'louvain')
#' vertex_attr_names(g)
#' graph_attr(g,'louvain')
calcClustering<-function(gg,alg){
    cl<-getClustering(gg,alg)
    if(!is.null(cl)){
    ids <- V(gg)$name
    m      <- matrix(NA, ncol=2, nrow=length(ids))
    colnames(m)<-c('ID',alg)
    m[,1]<-ids
    m[,2]<-as.character(cl$membership)
    ggm<-applpMatrixToGraph(gg,m)
    mod<-modularity(ggm,cl$membership)
    ggm<-set.graph.attribute(ggm,alg,mod)
    return(ggm)
    }else{
    return(gg)
    }
}

#' Get clustering results for the graph.
#' 
#' Wrapper function for calculation of clustering for predefined set of ten 
#' algorithms:
#' * lec -- leading eigenvector community (version of 
#' \code{\link[igraph]{leading.eigenvector.community}});
#' * wt -- walktrap community \code{\link[igraph]{walktrap.community}};
#' * fc -- fastgreedy community \code{\link[igraph]{fastgreedy.community}};
#' * infomap -- infomap community \code{\link[igraph]{cluster_infomap}};
#' * louvain -- cluster_louvain \code{\link[igraph]{cluster_louvain}};
#' * sgG1 -- spin-glass model and simulated annealing clustering (version of 
#' \code{\link[igraph]{spinglass.community}} with spins=500 and gamma=1);
#' * sgG2 -- spin-glass model and simulated annealing clustering (version of 
#' \code{\link[igraph]{spinglass.community}} with spins=500 and gamma=2);
#' * sgG5 -- spin-glass model and simulated annealing clustering (version of 
#' \code{\link[igraph]{spinglass.community}} with spins=500 and gamma=7);
#' * spectral -- spectral modularity clustering 
#' \code{\link[rSpectral]{spectral_igraph_communities}};
#'
#' graph suppose to be undirected. If algorithm failed warning will be issued
#' and function returned NULL.
#' 
#' Algorithm names are verified with \code{\link[base]{match.arg}}.
#' 
#' @md
#' @param gg igraph object to cluster
#' @param alg clustering algorithm name
#'
#' @return \code{\link[igraph]{communities}} object or NULL if algorithm failed.
#' @export
#'
#' @examples
#' data(karate,package='igraphdata')
#' c<-getClustering(karate,'lec')
#' c$modularity
getClustering<-function(gg,
                        alg=c('lec','wt','fc','infomap',
                                'louvain','sgG1','sgG2','sgG5','spectral')){
    alg <- match.arg(alg)
    lec<-function(gg){
    lec     <- igraph::leading.eigenvector.community(gg)
    ll      <- igraph::leading.eigenvector.community(gg, start=membership(lec))
    }
    cl <- try(switch(
        alg,
        lec = lec(gg),
        wt = igraph::walktrap.community(gg),
        fc = igraph::fastgreedy.community(gg),
        infomap = igraph::cluster_infomap(gg),
        louvain = igraph::cluster_louvain(gg),
        sgG1 = igraph::spinglass.community(gg,
                                            spins = as.numeric(500), gamma =
                                                1),
        sgG2 = igraph::spinglass.community(gg,
                                            spins = as.numeric(500), gamma =
                                                2),
        sgG5 = igraph::spinglass.community(gg,
                                            spins = as.numeric(500), gamma =
                                                5),
        spectral = rSpectral::spectral_igraph_communities(gg)
    ))
    if(inherits(cl, "try-error")){
    warning('Clustering calculations for algorithm "',alg,
            '" failed. NULL is returned')
    return(NULL)
    }
    return(cl)
}

#' Matrix of cluster characteristics
#' 
#' Function calculates statistics of clustering membership:
#' * mod -- clustering modularity \code{\link[igraph]{modularity}}
#' * C -- number of clusters
#' * Cn1 -- number of singletones (clusters of size 1)
#' * Cn100 -- number of clusters containing more than 100 nodes
#' * mu -- fraction of intercluster edges
#' * Min. C -- minimum of the cluster size
#' * 1st Qu. C -- first quartile of the cluster size
#' * Median C -- median of the cluster size
#' * Mean C -- average cluster size
#' * 3rd Qu. C  -- third quartile of the cluster size
#' * Max. C -- maximum of the cluster size
#'
#' @param gg graph to analyse
#' @param att vector of attribute names that contains membership data
#'
#' @return matrix of clustering characteristics
#' @export
#' @md
#'
#' @examples
#' data(karate,package='igraphdata')
#' g<-calcAllClustering(karate)
#' clusteringSummary(g)
clusteringSummary<-function(gg,
                            att=c('lec','wt','fc','infomap',
                                    'louvain','sgG1','sgG2','sgG5','spectral')){
    attN<-vertex_attr_names(gg)
    idx<-match(attN,att)
    clusterings<-attN[!is.na(idx)]
    res<-list()
    for(c in clusterings){
    cmem<-as.numeric(vertex_attr(gg,c))
    mod<-modularity(gg,cmem)
    Cn<-table(cmem)
    C <- length(Cn)
    Cn1<-length(which(Cn==1))
    Cn100<-length(which(Cn>=100))
    summary(as.vector(Cn))->s
    names(s)<-paste(names(s),'C')
    sgraphs<-lapply(names(Cn),getClusterSubgraphByID,gg=gg,mem=cmem)
    ug <- disjoint_union(sgraphs)
    mu<- 1-ecount(ug)/ecount(gg)
    r1<-c(mod,C,Cn1,Cn100,mu)
    names(r1)<-c('mod','C','Cn1','Cn100','mu')
    res[[c]]<-c(r1,s)
    }
    return(do.call(rbind,res))
}
