############################
# Code reproduce DYNAMO sensitivity matrix calculations from paper:
# Santolini,M. and Barabasi,A.-L. (2018) Predicting perturbation patterns from
# the topology of biological networks. Proc Natl Acad Sci USA, 169, 201720589.
#
# Structure of algorithm is inspired by code in
# https://github.com/msantolini/dynamo/

#' Convert sparce matrix into triplet \code{data.frame}.
#'
#' For very large graphs handling adjacency-like matrices is difficult due to
#' its sparse nature. This function convert sparse matrix into triplet
#' \code{data.frame} with row and column indices and names, and cell value.
#'
#' @param sparceM sparce matrix to convert into triplet \code{data.frame}
#'
#' @import Matrix
#' @return \code{data.frame} with three colums:
#' * i -- row index;
#' * j -- column index;
#' * x -- cell value;
#' * Rname -- i-th row name;
#' * Cname -- j-th column name.
#' @export
#' @examples
#' data(karate, package='igraphdata')
#' upgrade_graph(karate)
#' Ws <- as_adjacency_matrix(karate,type='both',attr='weight',sparse = TRUE)
#' mdf<-metlMatrix(Ws)
#' head(mdf)
metlMatrix<-function(sparceM){
    m <- as(sparceM, "TsparseMatrix")
    d <- data.frame(
        i = m@i + 1,  # m@i is 0-based, not 1-based like everything else in R
        j = m@j + 1,  # m@j is 0-based, not 1-based like everything else in R
        x = m@x)
    d$Rname <- colnames(m)[d$i]
    d$Cname <- colnames(m)[d$j]
    return(d)
}

#' Calculate DYNAMO sensitivity matrix.
#'
#' This function calculates sensitivity matrix that represents perturbation
#' patterns defined by topology and edge weights of the network. If weights
#' are signed value sensitivity matrix is able to reproduce not only activation
#' but inhibition relationships in the network.
#'
#' Algorithm proposed in:
#'
#' Santolini,M. and Barabasi,A.-L. (2018) Predicting perturbation patterns from
#' the topology of biological networks. Proc Natl Acad Sci USA, 169, 201720589.
#'
#'
#'
#' @param g igraph object
#' @param attr NULL or the name of edge attribute containing numerical weight
#'         values
#' @param vid name of the vertex attribute to be used as row and column names
#' @param alpha parameter characterizing the propagation strength, default
#'         value 0.9 taken from Santolini paper.
#'
#' @return sparce sensitivity matrix defined by the network topology and
#'         edge values
#' @export
#' @import Matrix
#' @importFrom igraph as_adjacency_matrix
#' @importFrom igraph edge_attr
#' @importFrom igraph edge_attr_names
#' @examples
#' data(karate, package='igraphdata')
#' upgrade_graph(karate)
#' d<-getDYNAMO(karate,attr='weight')
#' df<-metlMatrix(d)
#' head(df)
getDYNAMO<-function(g,attr=NULL,vid='name',alpha = .9){
    if(!is.null(attr) ){
        if(!(attr %in% edge_attr_names(g))){
            stop('Attr shold be either NULL or edge attribute name.\n')
        }
        attv<-edge_attr(g,attr)
        if(any(is.na(attv))){
            stop('Attr values shold not contain NAs.\n')
        }
        if(any(!is.numeric(attv))){
            stop('Attr values shold be numeric.\n')
        }
    }
    vcg<- length(V(g))
    Ws <- as_adjacency_matrix(g,type='both',attr=attr,sparse = TRUE)
    D1s <- Diagonal(n=vcg,
                    x=ifelse(sqrt(apply(abs(Ws),1,sum))==0,0,
                             1/sqrt(apply(abs(Ws),1,sum))))
    D2s <- Diagonal(n=vcg,
                    x=ifelse(sqrt(apply(abs(Ws),2,sum))==0,0,
                             1/sqrt(apply(abs(Ws),2,sum))))
    Wps <- D1s %*% Ws %*% D2s
    Is <- Diagonal(n=vcg)
    Msprince_dir_sign <- solve(Is-alpha * Wps) * (1-alpha)
    if(!is.null(vid) && (vid %in% vertex_attr_names(g))){
    nms<-vertex_attr(g,vid)
    colnames(Msprince_dir_sign) <- nms
    rownames(Msprince_dir_sign) <- nms
    }
    return(Msprince_dir_sign)
}


