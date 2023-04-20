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
#' rr<-merge(andf, res, by.y='FL', by.x='ID')
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
#' @param x enrichment statistics
#' @param desc plot subtitle
#' @param anno name of annotation used
#' @param LEGtextSize size of the text
#' @param LEGlineSize width of the line
#' @param type type of the plot
#'
#' @return ggplot object
#' @export
#' @import viridis
plotRatio <- function(x,
                      desc="", anno="",
                      LEGtextSize=1.5,
                      LEGlineSize=4,
                      type=NULL){
    #type=c('1','2','3','4')){


    #---For p.values
    xx  = x$SUM3 #Fe
    Nxx = length(colnames(xx))
    df  = data.frame()

    #--- labels
    xlab = colnames(xx)[3:Nxx]
    xval = seq(0,(length(xlab)-1),1)
    xlim = c(0,25,50,75,100)

    indx = match(xlim,xval)
    xval = xval[indx]
    xval = factor(xval)
    xlab = xlab[indx]
    brks = as.numeric(xlab)
    #---

    #--- test intervals
    X1 = 7
    X2 = 54
    X3 = 90

    rank <- matrix("",ncol=3,nrow=length(xx[,1]))

    for( i in 1:length(xx[,1]) ){

        size  = rep(1.8,length(xx[i,1]))
        alpha = rep(0.7,length(xx[i,1]))
        col   = rep("grey", length(xx[i,1]))
        xlabs = colnames(xx)[3:Nxx]

        tmp = cbind(xx[i,1],seq(1,(Nxx-2),1),as.numeric(xx[i,3:Nxx])/as.numeric(xx[i,2]), size, alpha, col, xlabs)

        df = rbind(df,tmp)

        zz0 = as.numeric(as.vector(xx[i,2]))
        zz  = as.numeric(as.vector(xx[i,3:Nxx]))

        rank[i,1] = xx[i,1]
        rank[i,2] = ifelse( zz0 == 0, 0, (zz[7]  - zz[54])/zz0)
        rank[i,3] = ifelse( zz0 == 0, 0, (zz[54] - zz[90])/zz0)

    }

    #---
    rank = rank[order(as.numeric(rank[,2]),decreasing=T),]
    colnames(rank) <- c("Alg","FracComsEnriched_log2(FE)>0.5_log2(FE)<4.8","FracComsEnriched_log2(FE)>4.8_log2(FE)<8.0")
    #write.table(rank, sprintf("ranking_%s.csv",desc), sep="\t", row.names=F, col.names=T, quote=F)


    colnames(df) <- c("ALG","X","Y","LSIZE","ALPHA","COL","XLAB")
    df           <- df[order(match(df[,1],rank[,1])),]
    df$ALG       <- factor(df$ALG, levels=rank[,1])
    df           <- as.data.frame(df)

    if( length(which(df$ALG == rank[1,1])) != 0 ){
        df$COL[df$ALG   == rank[1,1]]="royalblue"#lawngreen"
            df$ALPHA[df$ALG == rank[1,1]]= 1.0

    }

    if( length(which(df$ALG == rank[2,1])) != 0 ){
        df$COL[df$ALG   == rank[2,1]]="green2"
            df$ALPHA[df$ALG == rank[2,1]]= 1.0

    }


    if( length(which(df$ALG == rank[3,1])) != 0 ){
        df$COL[df$ALG   == rank[3,1]]="magenta"
            df$ALPHA[df$ALG == rank[3,1]]= 1.0
    }
    #---

    #---Generate plot
    SIZEa=2
    SIZEb=2

    #---legend
    #LEGtextSize=0.75 #ALL Algs
    #LEGtextSize=1.5   #Selected ALgs

    #LEGlineSize=2    #ALL Algs
    #LEGlineSize=4    #Selected ALgs

    colours = df$COL[match(levels(factor(df$ALG)),df$ALG)]
    gplot <- ggplot(df,aes(x=(as.numeric(df$X)),y=as.numeric(as.vector(df$Y)),colour=df$ALG))+
        geom_line(size=as.numeric(as.vector(df$LSIZE)),alpha=as.numeric(as.vector(df$ALPHA)))+
        labs(x="log2(Fe)",y="Fraction of Enriched Communities",title=anno)+
        theme(axis.title.x=element_text(face="bold",size=rel(1.5)),
              axis.title.y=element_text(face="bold",size=rel(1.5)),
              legend.text=element_text(face="bold",size=rel(LEGtextSize)),
              plot.title=element_text(face="bold",size=rel(1.5)),
              legend.position="bottom")+
        scale_color_manual("",breaks=c(levels(factor(df$ALG))),values=c(colours))+
        scale_y_continuous(expand=c(0,0),limits=c(0,1))+
#        scale_x_discrete(expand=c(0,0), limit=xval, labels=xlab)+
#        scale_x_continuous(expand=c(0,0), breaks = 2^brks, labels=xlab)+
        theme(panel.grid.major = element_line(colour = "grey40"),
              panel.grid.minor = element_line(colour="grey40",size=0.1),
              panel.background = element_rect(fill = "white"),
              panel.border = element_rect(linetype="solid",fill=NA))+
        geom_vline(xintercept=(as.numeric(X1)),colour="grey10",size=SIZEb,linetype=2,show.legend=F)+
        geom_vline(xintercept=(as.numeric(X2)),colour="grey10",size=SIZEb,linetype=2,show.legend=F)+
        geom_vline(xintercept=(as.numeric(X3)),colour="grey10",size=SIZEb,linetype=2,show.legend=F)+
        guides(color = guide_legend(override.aes = list(size=LEGlineSize)),
               alpha = 'none',
               size  = 'none')

    gplot2 <- ggplot(df,aes(x=log(as.numeric(df$X)),y=as.numeric(as.vector(df$Y)),colour=df$ALG))+
        geom_line(size=as.numeric(as.vector(df$LSIZE)),alpha=as.numeric(as.vector(df$ALPHA)))+
        labs(x="log(log2(Fe))",y="Fraction of Enriched Communities",title=anno)+
        theme(axis.title.x=element_text(face="bold",size=rel(1.5)),
              axis.title.y=element_text(face="bold",size=rel(1.5)),
              legend.text=element_text(face="bold",size=rel(LEGtextSize)),
              plot.title=element_text(face="bold",size=rel(1.5)),
              legend.position="bottom")+
        scale_color_manual("",breaks=c(levels(factor(df$ALG))),values=c(colours))+
        scale_y_continuous(expand=c(0,0),limits=c(0,1))+
#        scale_x_discrete(expand=c(0,0), limit=xval, labels=xlab)+
        scale_x_continuous(expand=c(0,0), breaks = brks, labels=xlab)+
        theme(panel.grid.major = element_line(colour = "grey40"),
              panel.grid.minor = element_line(colour="grey40",size=0.1),
              panel.background = element_rect(fill = "white"),
              panel.border = element_rect(linetype="solid",fill=NA))+
        geom_vline(xintercept=log(as.numeric(X1)),colour="grey10",size=SIZEb,linetype=2,show.legend=F)+
        geom_vline(xintercept=log(as.numeric(X2)),colour="grey10",size=SIZEb,linetype=2,show.legend=F)+
        geom_vline(xintercept=log(as.numeric(X3)),colour="grey10",size=SIZEb,linetype=2,show.legend=F)+
        guides(color = guide_legend(override.aes = list(size=LEGlineSize)),
               alpha = 'none',
               size  = 'none')

    gplot3 <- ggplot(df,aes(x=(as.numeric(df$X)),y=as.numeric(as.vector(df$Y)),colour=df$ALG))+
        geom_line(size=as.numeric(as.vector(df$LSIZE)),alpha=as.numeric(as.vector(df$ALPHA)))+
        labs(x="log2(Fe)",y="Fraction of Enriched Communities",title=anno)+
        theme(axis.title.x=element_text(face="bold",size=rel(1.5)),
              axis.title.y=element_text(face="bold",size=rel(1.5)),
              legend.text=element_text(face="bold",size=rel(LEGtextSize)),
              plot.title=element_text(face="bold",size=rel(1.5)),
              legend.position="bottom")+
        scale_color_viridis("",discrete = TRUE, option = "D")+
        scale_y_continuous(expand=c(0,0),limits=c(0,1))+
 #       scale_x_discrete(expand=c(0,0), limit=xval, labels=xlab)+
 #       scale_x_continuous(expand=c(0,0), breaks = 2^brks, labels=xlab)+
        theme(panel.grid.major = element_line(colour = "grey40"),
              panel.grid.minor = element_line(colour="grey40",size=0.1),
              panel.background = element_rect(fill = "white"),
              panel.border = element_rect(linetype="solid",fill=NA))+
        geom_vline(xintercept=(as.numeric(X1)),colour="grey10",size=SIZEb,linetype=2,show.legend=F)+
        geom_vline(xintercept=(as.numeric(X2)),colour="grey10",size=SIZEb,linetype=2,show.legend=F)+
        geom_vline(xintercept=(as.numeric(X3)),colour="grey10",size=SIZEb,linetype=2,show.legend=F)+
        guides(color = guide_legend(override.aes = list(size=LEGlineSize)),
               alpha = 'none',
               size  = 'none')

    gplot4 <- ggplot(df,aes(x=log(as.numeric(df$X)),y=as.numeric(as.vector(df$Y)),colour=df$ALG))+
        geom_line(size=as.numeric(as.vector(df$LSIZE)),alpha=as.numeric(as.vector(df$ALPHA)))+
        labs(x="log(log2(Fe))",y="Fraction of Enriched Communities",title=anno)+
        theme(axis.title.x=element_text(face="bold",size=rel(1.5)),
              axis.title.y=element_text(face="bold",size=rel(1.5)),
              legend.text=element_text(face="bold",size=rel(LEGtextSize)),
              plot.title=element_text(face="bold",size=rel(1.5)),
              legend.position="bottom")+
        scale_color_viridis("",discrete = TRUE, option = "D")+
        scale_y_continuous(expand=c(0,0),limits=c(0,1))+
#        scale_x_discrete(expand=c(0,0), limit=xval, labels=xlab)+
        scale_x_continuous(expand=c(0,0), breaks = brks, labels=xlab)+
        theme(panel.grid.major = element_line(colour = "grey40"),
              panel.grid.minor = element_line(colour="grey40",size=0.1),
              panel.background = element_rect(fill = "white"),
              panel.border = element_rect(linetype="solid",fill=NA))+
        geom_vline(xintercept=log(as.numeric(X1)),colour="grey10",size=SIZEb,linetype=2,show.legend=F)+
        geom_vline(xintercept=log(as.numeric(X2)),colour="grey10",size=SIZEb,linetype=2,show.legend=F)+
        geom_vline(xintercept=log(as.numeric(X3)),colour="grey10",size=SIZEb,linetype=2,show.legend=F)+
        guides(color = guide_legend(override.aes = list(size=LEGlineSize)),
               alpha = 'none',
               size  = 'none')

if(is.null(type)){
    return(list(p1=gplot,p2=gplot2,p3=gplot3,p4=gplot4,ranktable=rank))
        }else if(type==1){

        return(gplot)

    }else if(type == 2){

        return(gplot2)
    }else if(type == 3){
        return(gplot3)
    }else if(type == 4){
        return(gplot4)
}
}

