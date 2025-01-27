#' Calculates Magnetic Laplacian for directed network
#'
#' Algorithm is based upon following publication:
#'
#' ref: https://arxiv.org/pdf/1606.08266
#' ref: https://doi.org/10.1016/j.acha.2017.01.004
#' Magnetic Eigenmaps for the Visualization of Directed Networks
#' Micha ̈el Fanuel, Carlos M. Ala ́ız, A ́ngela Fern ́andez, Johan A. K. Suykens
#' KU Leuven, Department of Electrical Engineering (ESAT), Kasteelpark Arenberg 10, B-3001 Leuven, Belgium
#'
#' @param gg igraph object
#' @param g an electric charge parameter, which value should be in range
#' between 0 and 0.5.
#'
#' @return directed graph Magnetic Laplacian matrix
#' @export
#'
#' @examples
#' data(foodwebs,package='igraphdata')
#' g<-foodwebs[['ChesLower']]
#' upgrade_graph(g)
#' H <- getHmatrix(g)
getHmatrix <- function(gg, g=0.5){



    ## 0 <= g < 0.5
    ## g = 1/3, 1/4, 2/5

    ## Adjacency Matrix
    adj = as.matrix(as_adjacency_matrix(gg))
    W   = .5 * (adj + t(adj))
    D   = diag(rowSums(W))

    ## Initialize the Hermitian matrix
    N <- nrow(adj)
    H <- matrix(exp(1i*2*pi*g*0), N, N)

    ## Build the Hermitian matrix
    for (ii in 1:N) {
        rindx <- adj[ii, ]
        cindx <- which(rindx > 0)
        for (jj in cindx) {
            i <- ii
            j <- jj

            adj_ij <- adj[i, j]
            adj_ji <- adj[j, i]

            if (adj_ij > 0 && adj_ji == 0) {
                H[i,j] = exp(1i*2*pi*g*-1)
                H[j,i] = exp(1i*2*pi*g*1)
            }

            if (adj_ij == 0 && adj_ji > 0) {
                H[i,j] = exp(1i*2*pi*g*1)
                H[j,i] = exp(1i*2*pi*g*-1)
            }

        }
    }

    L  = D - H*W

}


#' Prepare Laplacian matrix from graph
#'
#'
#' @param gg igraph object
#' @param el an electric charge parameter for \code{\link{getHmatrix}}
#' @param weights an optional weights for weighted Laplacian matrix
#' @param type type of Laplacian to calculate for directed graph
#'
#' @return Lapcacian matrix
#' @importFrom igraph laplacian_matrix
#' @seealso getHmatrix
#' @export
getLaplacian <- function(gg,el=0.5,weights = NULL,type = c('igraph','magnetic')){
    type <- match.arg(type)
    L = switch(type,
        igraph = igraph::laplacian_matrix(gg,weights = weights),
        magnetic = getHmatrix(gg,g=el)
    )
    attr(L,'type')<-type
    return(L)

}

#' Calculate eigenvalues and eigenvectors of Laplacian matrix
#'
#' @param x Laplasian matrix from \code{\link{getLaplacian}}
#' @param only.values If TRUE only eigenvalues are returned.
#' @param inv.vec if TRUE inverted eigenvector matrix is returned.
#'
#' @return list with three slots: eigenvalues, eigenvectors and inverted
#' eigenvectors. If \code{only.values} is FALSE, both slots for eigenvectors and inverted
#' eigenvectors are NULL, if \code{inv.vec} is FALSE and \code{only.values} is TRUE, then
#' inverted eigenvectors slot is NULL.
#' @export
getEigen <- function(x, only.values=TRUE, inv.vec=FALSE){

    ## Eigenvalues of L
    Eigen = eigen(x, only.values=only.values, symmetric=TRUE)
    E     = Eigen$values
    if( only.values ){
        V    = NULL
        Vinv = NULL
    } else {
        V    = Eigen$vectors
        Vinv = NULL
        if( inv.vec ){ Vinv=solve(V); }
    }

    return(list(E=E, V=V, Vinv=Vinv))
}


#' Facade function for the LRG analysis.
#'
#' @param gg original graph
#' @param complex should rho(tau) be treated as complex matrix
#' @param method method of rho calculation
#' @param expm_method method to compute exponential of the matrix
#' \code{\link[expm]{expm}}
#' @param tol tolerance used to check if matrix for the exponential is
#' computationally singular \code{\link[expm]{expm}}
#' @param order the Taylor approxymation order for the rho calculation
#' @param n.steps
#' @param max_iter
#' @param dt
#' @param restarts
#' @param cooling_rate
#' @param seed
#' @param comb Specifies how to combine the vertex attributes in the
#' coars-grained graph. Please see \code{\link[igraph]{contract}}
#' for details.
#' @param simplify logical, should coarse-grained graph be simplified at all.
#' Please see \code{\link[igraph]{simplify}} for details.
#' @param remove.loops logical, should self-loops be removed in the
#' coarse-grained graph.
#' Please see \code{\link[igraph]{simplify}} for details.
#' @param remove.multiple logical, should parallel edges be removed in the
#' coarse-grained graph.
#' Please see \code{\link[igraph]{simplify}} for details.
#' @param verbatim turns extensive logging on, do not use in Rmd or within other functions
#'
#' @importFrom tibble tibble
#' @importFrom dplyr bind_rows distinct arrange
#' @importFrom Rbeast beast
#' @return
#' @export
runLRG<-function(gg, complex=FALSE,
                 method=c("eigen", "balanced", "square"),
                 expm_method=c("Higham08.b"), tol=1e-5, order=1,n.steps=100,
                 max_iter=100, dt=1, restarts=25, cooling_rate=0.85, seed=NULL,
                 comb='concat',simplify=TRUE,
                 remove.loops=TRUE, remove.multiple=TRUE,
                 verbatim=FALSE){
    if(verbatim){
        prnt<-1
        cat(format(Sys.time(), "%b %d %X"),'LRG analysis starts.\n')
    }else{
            prnt<-0
        }
    Nv <- vcount(gg)
    L = getLaplacian(gg)
    E = getEigen(L, only.values=FALSE, inv.vec=TRUE)

    Emean = E[[1]]
    Evec  = E[[2]]
    Evinv = E[[3]]
    if(verbatim){cat(format(Sys.time(), "%b %d %X"),'Laplacian calculated.\n')}

    # Create the first tibble with time intervals from 10^-4 to 10^4 with step 0.1
    df1 <- tibble(tau = seq(10^-4, 10^4, by = 0.1))

    # Create the second tibble with time intervals from 10^-4 to 1 with step 10^-2
    df2 <- tibble(tau = seq(10^-4, 1, by = 10^-2))

    # Combine the two tibbles
    df  <- bind_rows(df1, df2)

    # Remove duplicates and reorder the rows by tau in ascending order
    df  <- df %>%
        distinct() %>%  # Remove duplicates
        arrange(tau)    # Reorder by tau (ascending)
    y1      = S_wrapper(e=Emean, t=df$tau, n=Nv)
    y2      = dS_dlogt_wrapper(e=Emean, t=df$tau, n=Nv)
    not.nan = !is.nan(y2)
    y1      = y1[not.nan]
    y2      = y2[not.nan]
    t.set   = df$tau[not.nan]
    ## prepare plot
    tmp.y2 = y2

    #cut.off=0.1#0.05
    cut.off=y2[1]
    tmp.y2[tmp.y2<=cut.off]=cut.off

    b2        = range(tmp.y2)
    b3        = log(seq(b2[1],b2[2], length=4))
    b3        = exp(b3)
    b3        = signif(b3,1)
    y2.breaks = unique(c(0,b3))
    log.y2    = log(tmp.y2)

    a=range(y1)
    b=range(log.y2)
    dat = data.frame(x=t.set, y1=y1, y2=log.y2)
    if(verbatim){cat(format(Sys.time(), "%b %d %X"),'C = -dS(t)/d(log(t)) calculated.\n')}
    X<-10^(seq(min(log10(df$tau)),max(log10(df$tau)),length.out = n.steps+1))#df$tau
    Y<-dC_dt_test_wrapper(e=Emean,t=X,n=Nv)
    tps = turing.points(X,Y, seed = seed,verbatim = verbatim)
    if(verbatim){cat(format(Sys.time(), "%b %d %X"),'BEAST calculated.\n')}
    lrg = scan.tau(gg=gg, e=Emean, v=Evec, vinv=Evinv, L=L, tau=tps$tp$tp.x_axis,
                       complex=complex, method=method, order=order, print=prnt)
    t.upper = tps$tp$tp.max.x_axis[max(lrg$df$Iter)+1]
    t.lower = tps$tp$tp.min.x_axis[max(1,min(lrg$df$Iter)-1)]
    if(verbatim){cat(format(Sys.time(), "%b %d %X"),'Tau bounds found:[',t.lower,',',t.upper,'].\n')}
    lrg.ann = anneal.tau(gg=gg, L=L, e=Emean, v=Evec, vinv=Evinv,seed=seed,
                          t.lower=t.lower, t.upper=t.upper, cooling_rate=cooling_rate,
                          complex=complex, method=method, max_iter=max_iter, dt=dt, print=prnt)
    indx.best = nrow(lrg.ann$df)
    tau.best  = lrg.ann$df$new_t[indx.best]
    if(verbatim){cat(format(Sys.time(), "%b %d %X"),'Best Tau found: [',tau.best,'].\n')}
    pt6 <- ggplot(dat, aes(x))+
        geom_line(aes(y=y1, color="1-S"), linewidth=1)+
        geom_line(aes(y=(y2-b[1])/diff(b), color="C log(N)"),linewidth=1, alpha=.8)+
        ##geom_vline(xintercept=tau.best, linetype="dashed", color="black", linewidth=1)+
        scale_color_manual(values = c("1-S"   = "blue",
                                      "C log(N)" = "red"),
                           breaks = c("1-S","C log(N)"),
                           name = "")+
        labs(x = "log(tau)")+
        annotate("rect", xmin=t.lower, xmax=t.upper,
                 ymin=0, ymax=1, color="grey", alpha = .3)+
        geom_vline(xintercept=tau.best, linetype="dashed", color="black", linewidth=1)+
        scale_x_log10()+
        scale_y_continuous(
            name = "1-S",
            labels=~round(.,2),
            breaks=seq(0,1,.25),
            sec.axis=sec_axis(~.*b2[2],
                              name   = "log(C log(N))",
                              labels = scales::label_scientific(digits = 1),
                              breaks = y2.breaks
            ))+
        annotation_logticks(sides="b", outside=TRUE)+
        coord_cartesian(clip = "off")+
        theme_light()+
        theme(legend.position = "right")
    if(verbatim){cat(format(Sys.time(), "%b %d %X"),'Final plot prepared.\n')}
    R = real.LRG(e=Emean, v=Evec, vinv=Evinv, L=L, t=tau.best, gg=gg, method=method,
                 order=order,complex = complex,tol=tol,
                 comb=comb,simplify=simplify,
                 remove.loops=remove.loops, remove.multiple=remove.multiple)
    if(verbatim){cat(format(Sys.time(), "%b %d %X"),'Coarse-grained graph built.\n')}
    return(list(gg=gg,L=L,Emean=Emean, Evec=Evec, Evinv=Evinv,
                t=tau.best,tps=tps,R=R,finplot=pt6,logplot=lrg.ann$plots[[indx.best]]))

}

#' Laplacian Renormalisation Group in real-space
#'
#' @param e Laplacian eigenvalues from \code{\link{getEigen}}
#' @param v Laplacian eigenvectors from \code{\link{getEigen}}
#' @param vinv Laplacian inverted eigenvectors from \code{\link{getEigen}}
#' @param L Laplacian
#' @param t tau value for coarse-graining
#' @param gg original graph
#' @param complex should rho(tau) be treated as complex matrix
#' @param method method of rho calculation
#' @param expm_method method to compute exponential of the matrix
#' \code{\link[expm]{expm}}
#' @param tol tolerance used to check if matrix for the exponential is
#' computationally singular \code{\link[expm]{expm}}
#' @param order the Taylor approxymation order for the rho calculation
#' @param comb Specifies how to combine the vertex attributes in the
#' coars-grained graph. Please see \code{\link[igraph]{contract}}
#' for details.
#' @param simplify logical, should coarse-grained graph be simplified at all.
#' Please see \code{\link[igraph]{simplify}} for details.
#' @param remove.loops logical, should self-loops be removed in the
#' coarse-grained graph.
#' Please see \code{\link[igraph]{simplify}} for details.
#' @param remove.multiple logical, should parallel edges be removed in the
#' coarse-grained graph.
#' Please see \code{\link[igraph]{simplify}} for details.
#'
#' @return list with the following slots:
#' \itemize{
#'    \item rho - Calculated rho(tau).
#'    \item meta.adj - meta-graph adjacency table.
#'    \item mapping - node mapping between levels.
#'    \item supernodes - sizes of supernodes.
#'    \item gg - coarse-grained graph.
#' }
#'
#' @export
#' @import Deriv
#' @importFrom igraph contract
#' @importFrom igraph simplify
real.LRG <- function(e, v, vinv=NULL, L=NULL, t, gg, complex=FALSE,
                     method=c("eigen", "balanced", "square"),
                     expm_method=c("Higham08.b"), tol=1e-5, order=1,
                     comb='concat',simplify=TRUE,
                     remove.loops=TRUE, remove.multiple=TRUE){

    method <- match.arg(method)

    # method = c("approx", "eigen", "balanced")
    # expm_methods = c("Higham08.b", "Higham08",
    #          "AlMohy-Hi09", "Ward77", "PadeRBS",
    #.         "Pade", "Taylor", "PadeO", "TaylorO",
    #          "R_Eigen", "R_Pade", "R_Ward77",
    #.         "hybrid_Eigen_Ward")

    ## gg   == original graph
    ## e    == eigenvalues of laplacian
    ## v    == eigenvectors of laplacian
    ## t    == perform coase-graining at tau
    ## herm == is eigenvalues/vectors from complex laplacian

    gn = V(gg)$name

    ## Calculate rho(tau)
    rho = switch(method,
                 "approx"={ cal.rho.approx(L=L, t=t, order=order)},
                 "eigen"={cal.rho.eigen(e=e, v=v, vinv=vinv, t=t)},
                 "balanced"={cal.rho.balanced(L=L, t=t, method=expm_method, tol=tol, order=order)},
                 "square"={cal.rho.square(L=L, t=t, order=order)},
                 cal.rho.eigen(e=e, v=v, vinv=vinv, t=t)
    )

    ## build meta-graph
    adj  = meta.graph(rho=rho, complex=complex)

    colnames(adj) = gn
    rownames(adj) = gn

    sn = get.supernodes(adj=adj)

    ## record node mapping between levels
    mapping = cbind(gn, sn$membership)

    gg2 = coarse.grain.graph(gg=gg, supernodes=sn$membership)

    return(list(rho=rho, meta.adj=adj, mapping=mapping,
                supernodes=sn$csize, gg=gg2))

}

#' Scan vector tau values for goodness of fit.
#'
#' Function calculates coarse-grained graph for selected set of tau values
#' and compare degree distribution of the original and coarse-grained graphs.
#'
#' @param gg original graph see [real.LRG()]
#' @param e Laplacian eigenvalues see [real.LRG()]
#' @param v Laplacian eigenvectors see [real.LRG()]
#' @param vinv Laplacian inverted eigenvectors see [real.LRG()]
#' @param L Laplacian see [real.LRG()]
#' @param tau vector of tau values to calculate statistics on
#' @param complex should rho(tau) be treated as complex matrix see [real.LRG()]
#' @param method method of rho calculation see [real.LRG()]
#' @param expm_method method to compute exponential of the matrix
#' see [real.LRG()]
#' @param tol tolerance used to check if matrix for the exponential is
#' computationally singular see [real.LRG()]
#' @param order the Taylor approxymation order for the rho calculation see [real.LRG()]
#' @param print logical, does printing of intermediate results required
#'
#' @import poweRlaw
#' @seealso [real.LRG()]
#' @return list with the following slots:
#' \itemize{
#'    \item df - \code{data.frame} with calculated statistics.
#'    \item plots - list of \code{\link[ggplot2]{ggplot}} objects for succesful
#'              iterations.
#'    \item plots2 - list of \code{\link[ggplot2]{ggplot}} objects for succesful
#'              iterations.
#' }
#'
#' @export
scan.tau <- function(gg, e=NULL, v=NULL, vinv=NULL, L=NULL, tau,
                     complex=FALSE, method=c("eigen", "balanced", "square"),
                     expm_method=c("Higham08.b"), tol=1e-5, order=1, print=FALSE){

    method <- match.arg(method)

    res     = list()
    k       = 1

    ## get the degree distribution for the network
    d = as.numeric(degree(gg))

    ## find alpha & xmin for d using poweRlaw's MLE
    ## create new discrete power-law distribution
    X=poweRlaw::displ$new(d)

    ## the degree values
    x.range = X$internal$values

    ## estimate best xmin and alpha from sequence of xmin values
    X.est    = estimate_xmin(X, seq(1,max(d),1))
    xmin     = X.est$xmin
    zl       = sort(X$dat)
    zl       = zl[zl>=xmin]
    alpha    = X.est$pars

    ## store node probability distribution for original network
    gg.dat  = data.frame(x=X$internal$values,
                         y=data_cdf_probs(X$internal$values))

    plots  = list()
    plots2 = list()


    for( t in 1:length(tau) ){

        ## Perform coarse graining at tau: t
        lrg   = real.LRG(e=e, v=v, vinv=vinv, L=L, t=tau[t], gg=gg,
                         complex=complex, method=method,
                         expm_method=expm_method,
                         tol=tol, order=order)

        ## degree distribution of coarse-grained network
        lrg.d = as.numeric(degree(lrg$gg))

        if( length(lrg.d) > 1 ){

            ## find alpha & xmin for d using poweRlaw's MLE
            ## create new discrete power-law distribution
            cont = tryCatch({
                lrg.X=poweRlaw::displ$new(lrg.d)
                ## If successful, return TRUE
                TRUE
            },
            error = function(e) {
                ## If there's an error, return FALSE
                FALSE
            }
            )

            if( cont ){
                ## estimate best xmin and alpha from sequence of xmin values
                lrg.est    = estimate_xmin(lrg.X, seq(1,max(lrg.d),1))
                lrg.xmin   = lrg.est$xmin
                lrg.alpha  = lrg.est$pars

                ## store node probability distribution for coarse-grain network
                lrg.dat  = data.frame(x=lrg.X$internal$values,
                                      y=data_cdf_probs(lrg.X$internal$values))

                lrg.sn.sz.min = 1
                lrg.sn.sz.max = igraph::vcount(lrg$gg)
                lrg.sn        = table(lrg$supernodes)
                lrg.sn.sz     = as.numeric(names(lrg.sn))
                lrg.sn.sz     = lrg.sn.sz[-c(lrg.sn.sz.min, lrg.sn.sz.max)]
                lrg.d         = as.numeric(degree(lrg$gg))

                if( length(lrg.sn.sz) >= 1 && sum(lrg.d!=1) > 1 ){

                    gg_pdf  = pareto_powerlaw(x.range, xmin, alpha)
                    lrg_pdf = pareto_powerlaw(x.range, lrg.xmin, lrg.alpha)

                    ks_dist = ks_dist(x=zl, xmin=lrg.xmin, alpha=lrg.alpha)

                    ## calculate KL divergence
                    kl_pq = kl.divergence (ks_dist[[2]],ks_dist[[3]])
                    kl_qp = kl.divergence (ks_dist[[3]],ks_dist[[2]])

                    ## calculate JS divergence
                    js    = js.distance(x=cbind(ks_dist[[2]],ks_dist[[3]]),we=rep(0.5,2))

                    ## cross-entropy loss
                    loss = cross.entropy.loss(ks_dist[[2]],ks_dist[[3]])

                    res[[k]]    = c("Iter"=t,
                                    "ks_dist"=ks_dist[[1]], "kl_pq"=kl_pq, "kl_qp"=kl_qp,
                                    "js.div"=js$js.div, "jsd"=js$jsd, "js.norm"=js$js.norm,
                                    "loss"=loss, "lrg.sn"=length(lrg.sn),
                                    "alpha_gg"=alpha, "alpha_lrg"=lrg.alpha,
                                    "xmin_gg"=xmin, "xmin_lrg"=lrg.xmin,
                                    "t"=tau[t])

                    k=k+1

                    ## Print progress
                    if(print){cat(sprintf("Iteration: %d, ks_dist: %f, jsd: %f, loss: %f, lrg.sn: %d, alpha: %f, alpha.lrg: %f, xmin: %d, xmin.lrg: %d, t: %f\n",
                                          t, ks_dist[[1]], js$jsd, loss, length(lrg.sn), alpha, lrg.alpha, xmin, lrg.xmin, tau[t]))}


                    ## Log-Log plot of original and coarse-grain graph's degree distributions
                    plots[[t]] = log_plot(df.x=gg.dat,
                                              df.y=lrg.dat,
                                              x.lab="K", y.lab="LRG.K")


                    ## Log-Log plot of original and coarse-grain graph's degree distributions
                    plots2[[t]] = log_plot(df.x=data.frame(x=x.range,
                                                               y=gg_pdf),
                                               df.y=data.frame(x=x.range,
                                                               y=lrg_pdf),
                                               x.lab="K", y.lab="LRG.K")
                }else{
                    ## Print progress
                    if(print){cat(sprintf("Iteration: %d, length(lrg.sn.sz): %d, sum(lrg.d!=1): %d, lrg.sn: %d, t: %f\n",
                                          t, length(lrg.sn.sz), sum(lrg.d!=1), length(lrg.sn), tau[t]))}

                }

            }
        }
    }

    df = data.frame(do.call(rbind, lapply(res, unlist)))

    return(list(df=df,plots=plots, plots2=plots2))

}

#' Simulated Annealing
#'
#' Function calculates coarse-grained graph for selected set of tau values
#' and minimize the difference between degree distribution of the original
#' and coarse-grained graphs.
#'
#' @param gg original graph see [real.LRG()]
#' @param e Laplacian eigenvalues see [real.LRG()]
#' @param v Laplacian eigenvectors see [real.LRG()]
#' @param vinv Laplacian inverted eigenvectors see [real.LRG()]
#' @param L Laplacian see [real.LRG()]
#' @param t.lower minimal value for tau value search
#' @param t.upper maximal value for tau value search
#' @param complex should rho(tau) be treated as complex matrix see [real.LRG()]
#' @param method method of rho calculation see [real.LRG()]
#' @param expm_method method to compute exponential of the matrix
#' see [real.LRG()]
#' @param tol tolerance used to check if matrix for the exponential is
#' computationally singular see [real.LRG()]
#' @param order the Taylor approxymation order for the rho calculation see [real.LRG()]
#' @param print logical, does printing of intermediate results required
#' @param max_iter maximum of iterations
#' @param dt initial standard deviation for the search
#' @param restarts maximal number of restarts
#' @param cooling_rate cooling rate
#' @param seed random number generator seed
#'
#' @seealso [real.LRG()]
#' @importFrom pracma zeta isempty
#' @return list with the following slots:
#' \itemize{
#'    \item df - \code{data.frame} with calculated statistics.
#'    \item tau - vector of tested tau values.
#'    \item plots - list of \code{\link[ggplot2]{ggplot}} objects for succesful
#'              iterations.
#'    \item plots2 - list of \code{\link[ggplot2]{ggplot}} objects for succesful
#'              iterations.
#'    \item restarts - number of restarts left.
#' }
#' @export
anneal.tau <- function(gg, e=NULL, v=NULL, vinv=NULL, L=NULL, t.lower, t.upper, complex=FALSE,
                        method=c("eigen", "balanced", "square"),
                        expm_method=c("Higham08.b"), tol=1e-5, order=1, print=FALSE,
                        max_iter=100, dt=1, restarts=25, cooling_rate=0.99, seed=NULL) {

    method <- match.arg(method)
    res    <- list()
    plots  <- list()
    plots2 <- list()
    k      <- 1

    ## set random number seed if provided
    if( !is.null(seed) ){ set.seed(seed); }

    ## degree distribution of network
    d = as.numeric(degree(gg))

    ## find alpha & xmin for d using poweRlaw's MLE
    ## create new discrete power-law distribution
    X=poweRlaw::displ$new(d)

    ## network's degree values
    x.range = X$internal$values

    ## estimate best xmin and alpha from sequence of xmin values
    X.est    = estimate_xmin(X, seq(1,max(d),1))
    xmin     = X.est$xmin
    zl       = sort(X$dat)
    zl       = zl[zl>=xmin]
    alpha    = X.est$pars

    ## store node probability distribution for original network
    gg.dat  = data.frame(x=X$internal$values,
                         y=data_cdf_probs(X$internal$values))

    ## generate an initial time
    t     = t.lower[1] + runif(1) * (t.upper - t.lower)

    cont = FALSE
    loss  = NA

    while ( !cont ){

        ## Perform coarse graining at time: t
        lrg   = real.LRG(e=e, v=v, vinv=vinv, L=L, t=t, gg=gg,
                         complex=complex, method=method, expm_method=expm_method,
                         tol=tol, order=order)


        ## degree distribution of coarse-grain network
        lrg.d = as.numeric(degree(lrg$gg))

        ## Try to find alpha & xmin for d using poweRlaw's MLE
        output <- tryCatch(
            {
                ## create new discrete power-law distribution
                ## This is the code that might fail
                lrg.X = poweRlaw::displ$new(lrg.d)
                ## If successful, return TRUE
                TRUE
            },
            error = function(e) {
                ## If there's an error, return FALSE
                FALSE
            }
        )

        ## If error occurred, generate a new time t and try again
        if( !output && restarts > 0 ){
            t        = t.lower[1] + runif(1) * (t.upper - t.lower)
            restarts = restarts - 1
        } else {
            ## If no error, exit the loop
            cont = TRUE
        }
    }

    ## save time
    tau = t

    if( cont && restarts >= 0 ){

        ## estimate best xmin and alpha from sequence of xmin values
        lrg.est    = estimate_xmin(lrg.X, seq(1,max(lrg.d),1))
        lrg.xmin   = lrg.est$xmin
        lrg.alpha  = lrg.est$pars

        ## store node probability distribution for coarse-grain network
        lrg.dat  = data.frame(x=lrg.X$internal$values,
                              y=data_cdf_probs(lrg.X$internal$values))
        if( dim(gg.dat)[1]  > 1 &
            dim(lrg.dat)[1] > 1 ){

            ## preform powerlaw extrapolation of the two degree distributions
            gg_pdf  = pareto_powerlaw(x.range, xmin, alpha)
            lrg_pdf = pareto_powerlaw(x.range, lrg.xmin, lrg.alpha)

            gg_cdf  = cdf_powerlaw(x.range, xmin, alpha)
            lrg_cdf = cdf_powerlaw(x.range, lrg.xmin, lrg.alpha)

            ks_dist = ks_dist(x=zl, xmin=lrg.xmin, alpha=lrg.alpha)

            loss  = js.distance(we=c(1/2,1/2), x=cbind(p=ks_dist[[2]], q=ks_dist[[3]]))[[4]]
            curr_t      = t
            curr_loss   = loss
            deltaL      = NA

            res[[k]]    = c("Iter"=1, "Loss(old)"=NA, "Loss(new)"=loss,
                            #"Loss(delta)"=NA, "Prob:"=NA,
                            "alpha_gg"=alpha, "alpha_lrg"=lrg.alpha,
                            "xmin_gg"=xmin, "xmin_lrg"=lrg.xmin,
                            "t"=NA, "new_t"=t)#, "dt"=dt)

            plots[[k]]  = log_plot(df.x=gg.dat,
                                       df.y=lrg.dat,
                                       x.lab="K", y.lab="LRG.K")

            plots2[[k]] = log_plot(df.x=data.frame(x=x.range,
                                                       y=gg_pdf),
                                       df.y=data.frame(x=x.range,
                                                       y=lrg_pdf),
                                       x.lab="K", y.lab="LRG.K")

            k = k + 1

        }

        ## Print progress
        if(print){cat(sprintf("Iteration: %d, Loss(old): %f, Loss(new): %f, Loss(delta): %f, Prob: %f, t: %f, t(new): %f, dt: %f\n",
                              0, NA, loss, NA, NA, NA, t, dt))}



        for (iter in 1:max_iter) {

            ## Generate a new candidate state
            new_t = curr_t + rnorm(n=1, mean=dt)
            while( new_t < t.lower || new_t > t.upper ){
                new_t = curr_t + rnorm(n=1, mean=dt)
            }

            ## add time to set of tau
            tau = append(tau, new_t)

            ## Perform coarse graining at time: new_t
            lrg   = real.LRG(e=e, v=v, vinv=vinv, L=L, t=new_t, gg=gg,
                             complex=complex, method=method, expm_method=expm_method,
                             tol=tol, order=order)

            lrg.sn.sz.min = 1
            lrg.sn.sz.max = igraph::vcount(lrg$gg)
            lrg.sn        = table(lrg$supernodes)
            lrg.sn.sz     = as.numeric(names(lrg.sn))
            lrg.sn.sz     = lrg.sn.sz[-c(lrg.sn.sz.min, lrg.sn.sz.max)]
            lrg.d         = as.numeric(degree(lrg$gg))

            if( length(lrg.sn.sz) >= 1 && sum(lrg.d!=1) > 1 ){

                ## find alpha & xmin for d using poweRlaw's MLE
                ## create new discrete power-law distribution
                rm(lrg.X,lrg.est)
                lrg.X=poweRlaw::displ$new(lrg.d)

                ## estimate best xmin and alpha from sequence of xmin values
                lrg.est    = estimate_xmin(lrg.X, seq(1,max(lrg.d),1))
                lrg.xmin   = lrg.est$xmin
                lrg.alpha  = lrg.est$pars

                ## store node probability distribution for coarse-grain network
                lrg.dat  = data.frame(x=lrg.X$internal$values,
                                      y=data_cdf_probs(lrg.X$internal$values))


                if( dim(gg.dat)[1]  > 1 &
                    dim(lrg.dat)[1] > 1 ){


                    ## preform powerlaw extrapolation of the two degree distributions
                    gg_pdf  = pareto_powerlaw(x.range, xmin, alpha)
                    lrg_pdf = pareto_powerlaw(x.range, lrg.xmin, lrg.alpha)

                    gg_cdf  = cdf_powerlaw(x.range, xmin, alpha)
                    lrg_cdf = cdf_powerlaw(x.range, lrg.xmin, lrg.alpha)

                    ks_dist = ks_dist(x=zl, xmin=lrg.xmin, alpha=lrg.alpha)


                    ## Calculate cross-entropy loss
                    #loss     = old_loss
                    #new_loss = test_dist(x=zl, xmin=lrg.xmin, alpha=lrg.alpha)[[1]]
                    #new_loss = min(abs(log(gg_pdf+1e-8)-log(lrg_pdf+1e-8)))
                    new_loss    = js.distance(we=c(1/2,1/2), x=cbind(p=ks_dist[[2]], q=ks_dist[[3]]))[[4]]
                    #new_loss = js.distance(we=c(1/2,1/2), x=cbind(gg_pdf, lrg_pdf))[[4]]
                    #new_loss = as.vector(ks.test(gg_cdf, lrg_cdf)[[1]])
                    #new_loss = cross.entropy.loss(gg.dat$y, lrg.dat$y)
                    #new_loss = cross.entropy.loss(gg_pdf, lrg_pdf)

                    if( new_loss < loss ){
                        ## store best results
                        res[[k]]    = c("Iter"=iter, "Loss(old)"=loss, "Loss(new)"=new_loss,
                                        #"Loss(delta)"=deltaL,
                                        #"Prob:"=acceptance_prob,
                                        "alpha_gg"=alpha, "alpha_lrg"=lrg.alpha,
                                        "xmin_gg"=xmin, "xmin_lrg"=lrg.xmin,
                                        "t"=t, "new_t"=new_t)#, "dt"=dt)

                        plots[[k]]  = log_plot(df.x=gg.dat,
                                                   df.y=lrg.dat,
                                                   x.lab="K", y.lab="LRG.K")

                        plots2[[k]] = log_plot(df.x=data.frame(x=x.range,
                                                                   y=gg_pdf),
                                                   df.y=data.frame(x=x.range,
                                                                   y=lrg_pdf),
                                                   x.lab="K", y.lab="LRG.K")
                        ## store best time
                        loss = new_loss;
                        t    = new_t;
                        k    = k + 1
                    }

                    ## difference between new and old losses
                    deltaL = (new_loss-curr_loss)

                    ## Cooling schedule: Update t
                    dt = dt * cooling_rate

                    ## Calculate the acceptance probability
                    acceptance_prob <- exp(-abs(deltaL)/dt)

                    ## Accept the new state with a certain probability
                    if ( (deltaL < 0) || (runif(1) < acceptance_prob)) {

                        ## update current time and loss
                        curr_t = new_t; curr_loss = new_loss;

                    }

                    ## Print progress
                    if(print){cat(sprintf("Iteration: %d, Loss(old): %f, Loss(new): %f, Loss(delta): %f, Prob: %f, t: %f, t(new): %f, dt: %f\n",
                                          iter, curr_loss, new_loss, deltaL, acceptance_prob, curr_t, new_t, dt))}

                }
            } else {
                ## Print progress
                if(print){
                 cat(sprintf("Iteration: %d, Loss(old): %f, Loss(new): %f, Loss(delta): %f, Prob: %f, t: %f, t(new): %f, dt: %f\n",
                             iter, curr_loss, NA, NA, NA, curr_t, new_t, dt))}
            }
        }
    }

    if( !isempty(res) ){
        df = data.frame(do.call(rbind, lapply(res, unlist)))
        colnames(df) = c("Iter", "Loss.old", "Loss.new",
                         "alpha_gg", "alpha_lrg", "xmin_gg", "xmin_lrg", "t", "new_t")##, "dt")
    }

    return(list(df=df, tau=tau, plots=plots, plots2=plots2, restarts=restarts))

}



#### Rho calculation functions ####
cal.rho.balanced <- function(L, t, method="Higham08.b", order=1, tol=1e-5){
    #methods=c("Higham08.b", "Higham08",
    #          "AlMohy-Hi09",
    #          "Ward77", "PadeRBS", "Pade", "Taylor", "PadeO", "TaylorO",
    #          "R_Eigen", "R_Pade", "R_Ward77", "hybrid_Eigen_Ward")
    n   = nrow(L)
    rho = expm::expm(x=(-L*t), method=method, order=order, tol=tol)
    rho = rho/n
    rho = rho/sum(diag(rho))
    rho
}

taylor.approx <- function(L, t, order=3){
    n           <- nrow(L)
    I           <- diag(1, n)   ## Identity matrix
    term        <- I            ## Start with the first term (I)
    rho         <- I            ## Initialize the result with I
    factorial_k <- 1            ## k!

    for (k in 1:order) {
        factorial_k <- factorial_k * k             ## Compute k!
        term        <- term %*% (-L * t)           ## Compute (L t)^k
        rho         <- rho + term / factorial_k    ## Add the k-th term
    }
    rho
}

cal.rho.approx <- function(L, t, order=3) {
    n   <- nrow(L)
    rho <- taylor.approx(L=L, t=t, order=order)
    rho <- rho/n
    rho <- rho/sum(diag(rho))
    return(rho)
}


cal.rho.eigen <- function(e, v, vinv, t){
    n = length(e)

    if( is.null(vinv) ){ vinv = solve(v); }

    exp_nte = diag(exp(-t*e))
    St      = v %*% exp_nte %*% vinv
    St      = St/n
    rho     = St/sum(diag(St))
    rho
}

cal.rho.square <- function(L, t, order=3){

    n      <- nrow(L)

    ## Step 1: Calculate the norm of A
    norm_L <- sqrt(sum(L^2))

    ## Step 2: Determine m, the smallest power of two for which A/m has a sufficiently small norm
    m <- 1
    while (norm_L * t / m > 1) {
        m <- 2 * m
    }

    # Step 3: Compute the matrix exponential of -A*t/m
    L_scaled     <- L/m
    exp_L_scaled <- taylor.approx(L=L_scaled, t=t, order=order)

    # Step 4: Square the result m times
    rho <- exp_L_scaled
    for (i in seq_len(log2(m))) {
        rho <- rho %*% rho
    }

    rho = rho/n
    rho = rho/sum(diag(rho))
    rho
}

cal.rho <- function(e, v, vinv=NULL, t, complex=FALSE){

    N = length(e)

    if( is.null(vinv) ){ vinv = solve(v); }

    ##if( complex ){
    ##  exp_nte = diag(exp(-(0+1i)*t*e))
    ## else {
    exp_nte = diag(exp(-t*e))
    ##}

    St      = v %*% exp_nte %*% vinv
    St      = St/N
    rho     = St/sum(diag(St))

    #if( complex ){
    #  rho = Mod(rho)
    #  #rho = Re(rho) #Im(rho) ##Mod(rho)
    #}

    rho
}

meta.edge <- function(x, mn, mn_ii){
    (x/pmin(mn,mn_ii))
}

meta.binary  <- function(x, mn, mn_ii){
    as.numeric( (meta.edge(x,mn,mn_ii))-1 >= 0)
}

meta.graph <- function(rho,complex=FALSE){

    N    = dim(rho)[1]
    adj  = NULL
    if( complex ){
        re  = diag(Re(rho)) ## Probability of being at node i at time t
        im  = Im(rho)       ## Probability flow between node i to j at time t
        adj = matrix(0,N,N)

        for( i in 1:N ){
            edges   = im[i,]/min(re[i],re)
            adj[i,] = ifelse(edges-1>0 & sign(edges),1,0)
        }

    } else {
        mn   = diag(rho)
        adj  = sapply(1:N, function(i) meta.binary(x=rho[i,],
                                                   mn=mn,
                                                   mn_ii=mn[i]) )
    }

    return(adj)
}

get.supernodes <- function(adj){

    ## build network from rho...
    ## remember to add node names to columns/row of rho.
    cc = igraph::graph_from_adjacency_matrix(adj,
                                             mode=ifelse(isSymmetric(adj),
                                                         "undirected",
                                                         "directed"))

    ## function clusters gives a bit more information than decompose
    clusters  = igraph::clusters(cc)

}

coarse.grain.graph <- function(gg, supernodes,comb='concat',simplify=TRUE,
                               remove.loops=TRUE, remove.multiple=TRUE){

    ## group nodes into supernodes
    gg2 = igraph::contract(gg, supernodes,vertex.attr.comb = comb)

    ## remove self-loops and multiple edges from
    if(simplify){
    gg2 = igraph::simplify(gg2, remove.loops=remove.loops,
                           remove.multiple=remove.multiple)
    }

}



#### Helper functions ####

check.t <- function(t){
    if( t <= 1e-10 ){ t = 0 }
    if( t >= 1e10 ) { t = 1e10 }
    return(t)
}


rho.tau <- function(e,t){
    t = check.t(t)
    exp(-1*t*e)
}


u.tau <- function(e, t){
    small = 1e-30
    num   = exp(-1*t*e)
    dem   = sum(num)
    (num/dem) + small
}

## Entropy Measure given t=tau for graph Laplacian's eigenvalues
S.tau <- function(e, t, n){
    mu = u.tau(e=e,t=t)
    (-1/log(n))*sum(mu*log(mu))
}


## 1st derivative of u.tau
du_dt = Deriv(u.tau,"t")

## dS(t)/log(t)
dS_dlogt <- function(e,t,n){
    ut     = u.tau(t=t,e=e)
    du     = du_dt(t=t,e=e)
    log_ut = log(ut)
    (-1/log(n))*sum(du*t*log_ut)
}


dS_dt   <- Deriv(S.tau, "t")
d2S_dt2 <- Deriv(dS_dt, "t")

dC_dt_test <- function(e,t,n){
    -1*(d2S_dt2(e=e, t=t, n=n) * t + dS_dt(e=e, t=t, n=n))
}

dC_dt_test_wrapper <- function(e,t,n){
    sapply(1:length(t), function(i) dC_dt_test(e=e, t=t[i], n=n) )
}


turing.points <- function(X,Y, tp.range=2, seed=NULL,verbatim=FALSE){
    ## Find turning-point using rBeast
    ## get y values from plot
    out   = NULL
    if(is.null(seed)){
        mseed=0
    }else{
            mseed=seed
    }
    reg.x = NULL
        c.max = max(Y)
        out = beast(Y, season="none",
                    print.options=FALSE,
                    print.progress=!verbatim,
                    quiet=TRUE,
                    mcmc.seed=mseed)

        ## Turning point index
        tp    = as.vector(na.omit(out$trend$cp))
        n.tp  = length(tp)

        cn = c("tp", "tp.peak.pr","tp.x_axis", "tp.range", "tp.min",
               "tp.min.x_axis", "tp.max", "tp.max.x_axis","tp.cum_pr")
        tp.stats = matrix(nrow=n.tp, ncol=length(cn))
        colnames(tp.stats) = cn

        tp.stats[,1] = tp

        ## the probabilities associated with the change points cp
        tp.pr = as.vector(na.omit(out$trend$cpPr))
        tp.stats[,2] = tp.pr

        ## plot where change, or turning, points occurs on plot
        tp.stats[,3] = X[tp]

        ## range round the turning point to accumulate probability
        tp.stats[,4] = rep(tp.range,n.tp)

        ## the probability distribution of having a change point at each point of time
        tp.oc.pr = out$tren$cpOccPr

        ## the chance of observing a change point AROUND that point, ie the summed probability,
        ## which maybe greater than probability AT that point.
        for( i in 1:n.tp ){
            tp.i   = tp.stats[i,1]
            tp.min = tp.i-tp.stats[i,4]
            tp.max = tp.i+tp.stats[i,4]
            tp.cum_pr = sum(tp.oc.pr[(tp.min):(tp.max)])
            tp.stats[i,5] = tp.min
            tp.stats[i,6] = X[tp.min]
            tp.stats[i,7] = tp.max
            tp.stats[i,8] = X[tp.max]
            tp.stats[i,9] = tp.cum_pr
            if(verbatim){ cat("tp: ", tp.i, ", Sum Pr: ", tp.cum_pr, "\n") }
        }


        ## turning point round region of interest, ordered first to last
        reg.x = tp.stats[order(tp.stats[,1]),]
        reg.x = as.data.frame(reg.x)
    return(list(pt=out, tp=reg.x))

}

log_ab <- function(a,b,SMALL=1e-20){
    ## calculate a * log(a/b)
    a = as.numeric(a)
    b = as.numeric(b)
    a = a+SMALL
    b = b+SMALL

    return(a*log(a/b))

}

## Kullback–Leibler divergence
kl.divergence <- function(p,q){
    ## Ref: https://en.wikipedia.org/wiki/Kullback%E2%80%93Leibler_divergence
    n=length(p)
    return(sum(sapply(1:n,function(i) log_ab(p[i],q[i]))))
}

## generalised entropy function
H <- function(x){
    x = x[x>0]
    return(sum(-x*log(x)))
}

## jenson-shannon divergence for more than two distributions
js.divergence <- function(x,we){
    ##Ref:https://en.wikipedia.org/wiki/Jensen%E2%80%93Shannon_divergence
    ##    https://stackoverflow.com/questions/11226627/jensen-shannon-divergence-in-r
    ## x=rep(0,10)
    ## y=rep(1,10)
    ## m=0.5*(x+y)
    ## js.divergence(we=c(1/3,1/3,1/3), x=cbind(x,y,m)
    if(sum(we)!=1){
        stop('Weights sould sum up to 1: sum(',we,')=',sum(we),'\n')
    }
    return(H(x %*% we) - apply(x,2,H) %*% we)

}


js.distance <- function(x,we){
    if(sum(we)!=1){
        stop('Weights sould sum up to 1: sum(',we,')=',sum(we),'\n')
    }
    ## Ref: https://en.wikipedia.org/wiki/Jensen%E2%80%93Shannon_divergence
    ## Example:
    ## x=rep(0,10)
    ## y=rep(1,10)
    ## m=0.5*(x+y)
    ## js.divergence(we=c(1/3,1/3,1/3), x=cbind(x,y,m)

    n       = dim(x)[2]
    norm    = log(n)
    js.div  = js.divergence(we=we, x=x)
    jsd     = sqrt(js.div)
    js.norm = js.div/norm

    return(list(n=n, norm=norm, js.div=js.div, jsd=jsd, js.norm=js.norm))

}

cross.entropy.loss <- function(p,q){ kl.divergence(p,q) + H(p) }

powerlaw.constant <- function(x_min, alpha){
    ## p(x) = C*x^-alpha
    (alpha-1) * x_min^(alpha-1)
}


data_cdf_probs <- function(x){
    ## x ==> graph node degree frequency
    ## as.vector(degree(gg))
    ## Use for plotting
    ## https://github.com/csgillespie/poweRlaw/blob/main/R/dist_data_cdf_methods.R
    p = x/sum(x)
    p = rev(cumsum(rev(p)))
    p
}

pareto_powerlaw <- function(x, xmin, alpha){
    ((alpha-1)/xmin) * (x/xmin)^(-alpha)
}

ccdf_powerlaw <- function(x, xmin, alpha){
    ## https://en.wikipedia.org/wiki/Power_law#Plotting_power-law_distributions
    (x/xmin)^(-alpha+1)
}

pdf_powerlaw <- function(x, xmin, alpha, log = FALSE){
    ## using poweRlaw::dpldis function
    ## https://github.com/csgillespie/poweRlaw/blob/main/R/def_displ.R
    xmin     = floor(xmin)
    constant = zeta(alpha)
    if (xmin > 1)
        constant = constant - sum((1:(xmin - 1))^(-alpha))
    if (log) {
        pdf = -alpha * log(x) - log(constant)
        pdf[round(x) < round(xmin)] = -Inf
    }
    else {
        pdf = x^(-alpha)/constant
        pdf[round(x) < round(xmin)] = 0
    }

    ###if(!log) pdf = exp(pdf)

    pdf
}

cdf_powerlaw <- function(x, xmin, alpha, lower.tail = TRUE){
    xmin     = floor(xmin)
    constant = zeta(alpha)
    if (xmin > 1)
        constant = constant - sum((1:(xmin - 1))^(-alpha))
    cdf = 1 - (constant - sapply(x, function(i) sum((xmin:i)^(-alpha))))/constant
    if (!lower.tail)
        cdf = 1 - (cdf - pdf_powerlaw(x, xmin, alpha))
    cdf[round(x) < round(xmin)] = 0
    cdf
}


ks_dist <- function(x, xmin, alpha){
    ## chrome-extension://efaidnbmnnnibpcajpcglclefindmkaj/https://www.stat.berkeley.edu/~aldous/Research/Ugrad/Willy_Lai.pdf
    ## Distance between the empirical cdf of original network and fitted theoretical cdf
    n       = length(x)
    cdf_emp = n:1/n
    cdf_fit = ccdf_powerlaw(x, xmin, alpha)
    ks_dist = max(abs(cdf_emp-cdf_fit))
    return(list(ks_dist=ks_dist, cdf_emp=cdf_emp, cdf_fit=cdf_fit))
}



#### ggplot wrapper functions ####

log_plot <- function(df.x, df.y, x.lab="K", y.lab="LRG.K",title=NULL){

    colours        = c("red","blue")
    names(colours) = c(x.lab, y.lab)

    ## Log-Log plot of original and coarse-grain graph's degree distributions
    gp = ggplot()+
        geom_line(data=df.x, aes(log(x), log(y), color=x.lab), linewidth=1)+
        geom_line(data=df.y, aes(log(x), log(y), color=y.lab), linewidth=1)+
        labs(x = "log(K)", y = "P(K)",title = title) +
        ##scale_y_continuous(labels = scales::dollar) +
        scale_color_manual(values=colours, name="")+
        theme_light() +
        theme(legend.position = "bottom")

    return(gp)
}


## Return: 1-S
S_wrapper <- function(e,t,n, negate=1){
    if( negate ){
        sapply(1:length(t), function(i) 1-S.tau(e=e, t=t[i], n=n) )
    } else {
        sapply(1:length(t), function(i) S.tau(e=e, t=t[i], n=n) )
    }
}

## Return: C = -dS(t)/d(log(t))
dS_dlogt_wrapper <- function(e,t,n, scale=TRUE){
    if( scale ){
        sapply(1:length(t), function(i)
            -log(n)*dS_dlogt(e=e, t=t[i], n=n) )
    } else {
        sapply(1:length(t), function(i)
            -1*dS_dlogt(e=e, t=t[i], n=n) )
    }
}


