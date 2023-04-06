#require("minpack.lm")
#require(cowplot)

### Sigmoid function ### create a function to generate sigmoid pattern
sigmoid <- function(pars, xx){

    a = as.numeric(pars[1])#lower asymptote,          ideal == 0
    b = as.numeric(pars[2])#upper asymptote,          ideal == 1
    c = as.numeric(pars[3])#gradiant, rate, or slope, ideal == -2
    d = as.numeric(pars[4])#inflextion point,         ideal == median(xx) = 3

    return( a + ((b-a)/(1+exp(-c*(xx-d)))) )
}

residFun <- function(pars, observed, xx) observed - sigmoid(pars,xx)

#' Plot results of the sigmoid fit
#'
#' @param x
#' @param rates
#' @param model
#' @param alg
#' @param pv
#'
#' @return
#' @export
plotSigmoid <- function( x, rates, model, alg="", pv=0 ){
}

addNoise <- function( Y, MN=0, SD=0.05 ){
    return( Y+rnorm(length(Y),mean=MN, sd=SD) )
}


##goodness of fit test, KS
#' Goodnes of fit KS test
#'
#' @param x
#' @param rate
#' @param model
#' @param sigma2
#' @param countDATA
#'
#' @return
#' @export
gofs <- function(x, rate, model, sigma2=NULL, countDATA=TRUE ){

    y    = model$m$lhs()
    yhat = fitted(model)

    R  = length(rate)
    KS = list()

    for( r in 1:R ){

        pp = list(a=0, b=1, c=rate[r], d=round(median(x)) )
        yi = sigmoid(pars=pp, xx=x)

        KS[[r]]      = ks.test(yhat, yi)
        names(KS)[r] = sprintf("rate_%f.1", rate[r])
    }

    return(KS)

}

#' Fit Fe distribution to sigmoid function
#'
#' @param df enrichment \code{data.frame}
#' @param SDv vector of noise SD values
#'
#' @return list of fitted functions tables and plots
#' @export
fitSigmoid<-function(df,SDv){

}
