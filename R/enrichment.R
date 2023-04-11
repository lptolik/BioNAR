#' Calculate annotation enrichment for clusters in the graph
#'
#' Calculate the cluster enrichment of a graph given a clustering algorithm
#' 'alg' and vertex annotation attribute 'name'. Function generates an
#' enrichment table, one row for each cluster, containing: the cluster ID, the
#' cluster size, overlap of annotation terms in cluster, p.value of enrichment
#' using the Hypergeometric test, adjusted p.value Bonferroni correction (BH).
#'
#' @param g graph to get annotation from
#' @param alg cluster algorithm and membership attribute name
#' @param name annotation attribute name
#' @param col list separation character in attribute, by
#' default is \code{;}
#' @param vid attribute to be used as a vertex ID
#' @param alpha probability threshold
#'
#' @return A table with overrepresentation results.
#' Each row corresponds to a tested annotation in particular cluster.
#' The columns are the following:
#' \itemize{
#'   \item pathway – name of the enriched term as in 'names(pathway)';
#'   \item pval – an enrichment p-value from hypergeometric test;
#'   \item padj – a BH-adjusted p-value;
#'   \item overlap – size of the overlap;
#'   \item size – size of the gene set;
#'   \item leadingEdge – vector with overlapping genes.
#'   \item cl – cluster ID
#'   }
#' @export
#' @import data.table
#' @importFrom fgsea fora
#' @examples
#' options("show.error.messages"=TRUE)
#' file <- system.file("extdata", "PPI_Presynaptic.gml", package = "BioNAR")
#' g <- igraph::read.graph(file, format="gml")
#' anL<-getAnnotationVertexList(g, 'TopOntoOVGHDOID')
#' res<-clusterORA(g, alg='louvain', name='TopOntoOVGHDOID', vid='name')
#' andf<-unique(data.frame(ID=get.vertex.attribute(g, 'TopOntoOVGHDOID'),
#' Term=get.vertex.attribute(g, 'TopOntoOVG')))
#' rr<-merge(andf, res, by.y='pathway', by.x='ID')
#' rr[order(rr$cl), ]
clusterORA <- function(g,
                       alg,
                       name,
                       vid = 'name',
                       alpha = 1.0,
                       col = COLLAPSE) {
    anL <- getAnnotationVertexList(g, name)
    cl <- make_clusters(g, as.numeric(get.vertex.attribute(g, alg)))
    vcnt <- vcount(g)
    forafun <- function(.i) {
        gids<-which(membership(cl) == .i)
        fres <- fora(
            anL,
            get.vertex.attribute(g, vid)[gids],
            universe = as.character(get.vertex.attribute(g, vid))
        )
        cn<-length(gids)
        den<-cn/vcnt
        res <- fres[,.(cl=.i,FL=pathway,N=vcnt,F=size,Cn=cn,
                    Mu=overlap,
                    OR=log(overlap*(vcnt-size+overlap-cn)/
                              ((cn-overlap)*(size-overlap))),
                    CIw= sqrt(1/overlap+1/(cn-overlap)+1/(size-overlap)+
                              1/(vcnt-size+overlap-cn)),
                    Fe=(overlap/size)/den,
                    Fc=(overlap/cn)/den,
                    pval,padj,
                    overlapGenes)]
        res$palt <- sapply(1:dim(res)[1], function(.i){fisher(res$Mu[.i],res$F[.i],res$Cn[.i],res$N[.i],alternative = 'less')})

        res<-as.data.frame(res[,.(alg=alg,cl,FL,N,F,Cn,Mu,OR,
                                  CIl=OR-1.96*CIw,CIu=OR+1.96*CIw,
                               Fe,Fc,pval,padj,palt,paltadj=p.adjust(palt,method = 'BH'),
                               overlapGenes)])

        return(res)
    }
    resL <- lapply(seq_along(cl), forafun)
    res <- do.call(rbind, resL)
    res <- res[res$padj < alpha, ]
    #res$overlapGenes<-sapply(res$overlapGenes,paste,collapse = ', ')
    res$overlapGenes <-
        unlist(lapply(res$overlapGenes, paste, collapse = ', '))
    return(res)
}

fisher <- function( mu, P, F, N,alternative = 'less' ){

    mm  <- matrix(c(mu,(P-mu),(F-mu),(N-P-F+mu)),2,2)
    res <- fisher.test(mm,alternative = alternative)
    return(res$p.value)
}

lsum<-function(x){
    return(length(which(x)))
}
#' Calculate summary statistics from enrichment table
#'
#' @param RES enrichment results \code{data.frame}
#' @param ALPHA p-value cut-off
#' @param usePadj logical, wether to use plain or adjusted p-value
#' @param FeMAX max of the FE
#' @param FcMAX max of the FC
#'
#' @return list of \code{data.frame}
#' @export
summaryStats <- function( RES, ALPHA, usePadj=FALSE, FeMAX=0, FcMAX=0 ){

    cand <- data.frame(a=as.character(),b=as.character(),c=as.character(),d=as.character())

    CN <- colnames(RES[[1]])

    ALGi   = which(CN=="alg")[1]
    Pvi    = which(CN=="pval")[1]
    PvALTi = which(CN=="palt")[1]
    if( usePadj ){
        Pvi    = which(CN=="padj")[1]
        PvALTi = which(CN=="paltadj")[1]
    }
    ORi    = which(CN=="OR")[1]
    CIli   = which(CN=="CIl")[1]
    Fei    = which(CN=="Fe")[1]
    Fci    = which(CN=="Fc")[1]

    Ci     = which(CN=="cl")[1]
    Fli    = which(CN=="FL")[1]
    Cni    = which(CN=="Cn")[1]
    Mui    = which(CN=="Mu")[1]

    Ni     = which(CN=="N")[1]
    N      = as.numeric(RES[[1]][1,Ni])

    hh0  = c("alg","FN","CN","FNxCN","Psig","PALTsig","OR>1","ORsig","Psig&ORsig","PALTsig&ORsig","FEsig","Psig&ORsig&FEsig","EnrichedComs(%)","p.value")
    sum1 = matrix("",ncol=length(hh0),nrow=length(names(RES)))
    colnames(sum1) <- hh0

    cmin = seq(0,10,0.1)
    hh   = sprintf("Cmin_%.1f&Cmax_%.1f",cmin,10)
    sum2 = matrix("",nrow=length(names(RES)), ncol=(2+length(hh)) )
    colnames(sum2) = c("Alg","Psig&ORsig",hh)

    steps = 100
    hh2i  = round(FeMAX/steps,3)
    femin = seq(0,ceiling(FeMAX), hh2i)
    hh2   = sprintf("%.1f",femin)
    sum3  = matrix("",nrow=length(names(RES)), ncol=(2+length(hh2)) )
    colnames(sum3) = c("Alg","Psig&ORsig",hh2)

    hh3i  = round(FcMAX/steps,3)
    fcmin = seq(0,ceiling(FcMAX), hh3i)
    hh3   = sprintf("%.1f",fcmin)
    sum4  = matrix("",nrow=length(names(RES)), ncol=(2+length(hh3)) )
    colnames(sum4) = c("Alg","Psig&ORsig",hh3)

    for( i in 1:length(RES) ){

        Ncn = length(RES[[i]][,1])
        P   = as.numeric(RES[[i]][,Pvi])
        Palt= as.numeric(RES[[i]][,PvALTi])
        OR  = as.numeric(RES[[i]][,ORi])
        CI  = as.numeric(RES[[i]][,CIli])
        FE  = as.numeric(RES[[i]][,Fei])
        FC  = as.numeric(RES[[i]][,Fci])

        CNo = as.numeric(RES[[i]][,Cni])

        Cmax = max(as.numeric(RES[[i]][,Ci]))

        sum1[i,1] = as.character(RES[[i]][1,ALGi])

        sum1[i,2] = Ncn / Cmax

        sum1[i,3] = Cmax

        sum1[i,4] = Ncn ##

        sum1[i,5] = lsum(P <= ALPHA)

        sum1[i,6] = lsum(Palt <= ALPHA)

        sum1[i,7] =  lsum(OR > 1)

        sum1[i,8] =  lsum(OR > 1 & CI > 1)

        sum1[i,9] =  lsum(OR > 1 & CI > 1 & P <= ALPHA)

        sum1[i,10] =  lsum(OR > 1 & CI > 1 & Palt <= ALPHA)

        sum1[i,11] = lsum( log2(FE) > 0.5 & log2(FE) < 4.8 )

        sum1[i,12] = lsum(OR > 1 & CI > 1 & P <= ALPHA & log2(FE) > 0.5 & log2(FE) < 4.8 )

        sum2[i,1] =  as.character(RES[[i]][1,ALGi])
        sum2[i,2] =  lsum(OR > 1 & CI > 1 & P <= ALPHA)
        for( j in 1:length(hh) ){
            sum2[i,(j+2)] = lsum(OR > 1 & CI > 1 & P <= ALPHA & CNo > ((cmin[j]*N)/100) & CNo < ((10*N)/100) )
        }

        sum3[i,1] =  as.character(RES[[i]][1,ALGi])
        sum3[i,2] =  lsum(OR > 1 & CI > 1 & P <= ALPHA)
        for( j in 1:length(hh2) ){
            sum3[i,(j+2)] = lsum(OR > 1 & CI > 1 & P <= ALPHA & log2(FE) > femin[j] )
        }

        sum4[i,1] =  as.character(RES[[i]][1,ALGi])
        sum4[i,2] =  lsum(OR > 1 & CI > 1 & P <= ALPHA)
        for( j in 1:length(hh3) ){
            sum4[i,(j+2)] = lsum(OR > 1 & CI > 1 & P <= ALPHA & log2(FC) > fcmin[j] )
        }

        sum1[i,13] = (as.numeric(sum1[i,9]) / as.numeric(sum1[i,4])) * 100

        sum1[i,14] = fisher(as.numeric(sum1[i,12]), as.numeric(sum1[i,9]),
                            as.numeric(sum1[i,11]), as.numeric(sum1[i,4]))

        #store candidate enriched functional communities
        indx <- which(OR > 1 & CI > 1 & P <= ALPHA & log2(FE) > 0)

        n    <- length(indx)
        A    <- rep(sum1[i,1],n)
        B    <- RES[[i]][indx,Fli]
        C    <- RES[[i]][indx,Ci]
        D    <- RES[[i]][indx,Mui]
        cand <- rbind(cand,data.frame(A,B,C,D))
    }

    colnames(cand) <- c("ALG","Fn","C","Mu")

    return(list(SUM=sum1,
               SUM2=sum2,
               SUM3=sum3,
               SUM4=sum4,
               CAN=cand))
}


#' Plot fraction of enriched communities
#'
#' @param xx enrichment statistics
#' @param desc plot subtitle
#' @param anno name of annotation used
#' @param LEGtextSize size of the text
#' @param LEGlineSize width of the line
#' @param type type of the plot
#'
#' @return ggplot object
#' @export
#'
#' @examples
plotRatio <- function(xx,
                      desc="", anno="",
                      LEGtextSize=1.5,
                      LEGlineSize=4,
                      type=c('1','2','3','4')){
}
