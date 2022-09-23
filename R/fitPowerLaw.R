#library(igraph);
#library(lattice);
#library(methods);
#library(poweRlaw);
#library(scales);
#library(grid);
#library(latex2exp);
#library(stringr);
##---WIDTH and HEIGHT for plots
WIDTH <- 480
HEIGHT <- 480
##set igraph as S4
setClass("poweRlaw")
#' Result of PawerLaw fit
#'
#' @slot fit \code{\link[poweRlaw]{displ-class}} result of power law fit. 
#' @slot p numeric. 
#' @slot alpha numeric degree of power-law. 
#' @slot SDxmin numeric bootstrap sd of Xmin. 
#' @slot SDalpha numeric bootstrap sd of alpha. 
#'
#' @export
setClass(Class="law",representation(
    fit="displ",
    p="numeric",
    alpha="numeric",
    SDxmin="numeric",
    SDalpha="numeric")
)
##log sequence of numbers
lseqBy <- function(from=1, to=100000, by=1, length.out=log10(to/from)+1) {
    tmp <- exp(seq(log(from), log(to), length.out = length.out))
    tmp[seq(1, length(tmp), by)]
}
changeSciNot <- function(n) {
    n <- format(n, scientific = TRUE)
    oo <- strsplit(as.character(n),"e")
    out <- vector(length=length(oo))
    out[1] <- TeX(sprintf("10^{%s}",sub("-0?","-",oo[[1]][2])))
    for( i in 2:length(oo) ){
    if(grepl("-",oo[[i]][2])){
        out[i] <- TeX(sprintf("10^{%s}",sub("-0?","-",oo[[i]][2])))
    } else {
        out[i] <- TeX(sprintf("10^{%s}",sub("\\+0?","",oo[[i]][2])))
    }
    }
    return(out)
}
#' Fit Power Law to degree distribution.
#'
#' @param DEG degree distribution
#' @param Nsim number of bootstrap iterations
#' @param DATAleg legend string for degree data
#' @param WIDTH width of the plot in ptx
#' @param HEIGHT heigth of the plot in ptx
#' @param plot logical, do you want plot to be drawn
#' @param threads number of parallel computational threads
#' @param legpos position of the legend @seealso{legend}
#'
#' @return an object of class \code{\link{law}} with results of fitting
#' @export
#' @import poweRlaw latex2exp methods grid scales
#' @importFrom stringr str_sub
#'
#' @examples
#' ##No: of bootstrap iterations
#' nsim <- 100
#'
#' ##Legend Titles
#' Legend <- "Presynaptic PPI"
#'
#' file <- system.file("extdata", "PPI_Presynaptic.gml", package = "BioNAR")
#' gg <- igraph::read.graph(file,format="gml")
#' pFit <- FitDegree( as.vector(igraph::degree(graph=gg)),
#' DATAleg=Legend,threads=1,Nsim=nsim)
FitDegree <- function(DEG,Nsim=100,  plot=FALSE,
                        DATAleg='Fit power-law', threads=4,
                        WIDTH=480, HEIGHT=480 ,legpos="bottomleft"){
    DEG <- DEG[DEG > 0]
    data <- DEG
    m_pl <- displ$new(data)
    est <- estimate_xmin(m_pl)
    m_pl$setXmin(est)
    suppressMessages(
        gof <- poweRlaw::bootstrap_p(m_pl, no_of_sims = Nsim, threads=threads)
    )
    if(plot){
    op<-options(warn= -1)
    x_lab <- "k"    ##degree
    y_lab <- "P(k)" ## the CDFs P(k) for the PPI network data
    leg_x <- max(data)
    leg_y <- 1
    d <- plot(m_pl,draw=FALSE)
    Xmax <- max(d$x) - max(d$x)*0.5
    yTICKS <- round(lseqBy(min(d$y),1,0.5),4)
    yLABELS <- changeSciNot(yTICKS)
    poweRlaw::plot(m_pl, xlab=sprintf("%s",x_lab), ylab=y_lab,
            panel.first=graphics::grid(col="grey60"),
            pch=22, bg='black', axes = FALSE, cex.lab = 1.5, yaxt='n' )
    graphics::box(col='black')
    graphics::axis(1, cex.axis = 1.5, font = 1.5, family = 'sans')
    graphics::axis(2, cex.axis = 1.5, font = 1.5, family = 'sans', at=yTICKS,
            labels=yLABELS)
    poweRlaw::lines(m_pl, col=2, lwd=3)
    S1 <- round(m_pl$xmin,2)
    S2 <- round(m_pl$pars,2)
    S3 <- round(gof$p,2)
    sdS1 <- round(stats::sd(gof$bootstraps$xmin),0)
    sdS2 <- round(stats::sd(gof$bootstraps$pars),2)
    errS1 <- str_sub(as.character(sdS1),-1,-1)
    errS2 <- str_sub(as.character(sdS2),-1,-1)
    suppressMessages(
        fitl <-
            TeX(sprintf(
                "Power-law $\\alpha = %.2f(%s), $k_{min} = %.0f(%s)",
                S2,
                errS2,
                S1,
                errS1
            ))
    )
    graphics::legend(legpos,c(DATAleg,fitl),lty=c(1,1),lwd=c(4,4),
            col=c('black',2),merge=TRUE, cex = 1.5)
    options(op)
    }
    return(new("law",
                fit=m_pl,
                p=as.numeric(gof$p),
                alpha=as.numeric(est$pars),
                SDxmin=as.numeric(stats::sd(gof$bootstraps$xmin)),
                SDalpha=as.numeric(stats::sd(gof$bootstraps$pars)))
    )
}
