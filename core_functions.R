## Normalise small numbers
norm <- function(x){
    
    x.max  = max(x)
    result = x.max - log(exp(x-x.max)) 
    
    return(result)
}

lcc <- function(gg, min.vertices=5){
  gg  = igraph::simplify(gg)
  dec = decompose(gg, min.vertices=min.vertices)
  n.dec = lapply(dec, function(x)length(V(x)))
  indx.dec = which.max(unlist(n.dec))
  dec[[indx.dec]]
}

read_hint_dat <- function(x){
  ed  = x[,1:2]
  gg  = graph_from_edgelist(as.matrix(ed),directed=FALSE)
  lcc(gg)
}

rm.non_connect_nodes <- function(gg){
  D = degree(gg)
  return(igraph::delete.vertices(gg,names(D[D==0])))
}

get.H <- function(gg, g=0.5){
  
  ## ref: https://arxiv.org/pdf/1606.08266 
  ## Magnetic Eigenmaps for the Visualization of Directed Networks
  ## Micha ̈el Fanuel, Carlos M. Ala ́ız, A ́ngela Fern ́andez, Johan A. K. Suykens
  ## KU Leuven, Department of Electrical Engineering (ESAT), Kasteelpark Arenberg 10, B-3001 Leuven, Belgium
  
  
  ## 0 <= g < 0.5
  ## g = 1/3, 1/4, 2/5
  
  ## Adjacency Matrix
  adj = as.matrix(get.adjacency(gg))
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


get.L <- function(gg,g=0.5){

    L = NULL
    
    ## is directed
    if( is.directed(gg) ){
      
      L = get.H(gg,g=g)
      
    } else {
        
        ## Degree Matrix
        D   = diag(degree(gg))

        ## Laplacian
        L   = D - as.matrix(get.adjacency(gg))
    }
    
    return(L)
    
}

get.eigen <- function(x, only.values=TRUE, inv.vec=FALSE){

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

Z.tau <- function(e, t){
  small = 1e-30
  sum(exp(-1*t*e))+small
}

Zp.tau <- function(e, t){
  small = 1e-30
  sum(e*exp(-1*t*e))+small
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

## ggplot wrapper functions

## calculate C wrapper
## Z(t)      = Sum exp(-e_i*t)
## lambda(t) = Sum e_i * exp(-e_i*t) / Z(t)
## C0        = lambda(t) * t
C0_wrapper <- function(e,t){
  Z_t      = sapply(t, function(i) Z.tau(e,i)) 
  Zp_t     = sapply(t, function(i) Zp.tau(e,i))
  lambda_t = Zp_t/Z_t
  return(t*lambda_t)
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


meta.edge <- function(x, mn, mn_ii){
    (x/pmin(mn,mn_ii))
}

meta.binary  <- function(x, mn, mn_ii){
    as.numeric( (meta.edge(x,mn,mn_ii))-1 >= 0)
}

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
    ##clusters  = igraph::clusters(cc)
    clusters  = igraph::components(cc)
    
}

coarse.grain.graph <- function(gg, supernodes){

    ## group nodes into supernodes
    gg2 = igraph::contract(gg, supernodes)

    ## remove self-loops and multiple edges from  
    gg2 = igraph::simplify(gg2, remove.loops=T, remove.multiple=T)
   
}


## Laplacian Renormalisation Group in real-space
real.LRG <- function(e, v, vinv=NULL, L=NULL, t, gg, complex=FALSE, 
                     method=c("eigen", "balanced", "square"), 
                     expm_method=c("Higham08.b"), tol=1e-5, order=1){
  
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
  #rho = cal.rho(e=e, v=v, vinv=vinv, t=t)
  
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


Fourier.LRG <- function(e, v, vinv=NULL, t, t_x, complex=FALSE, cal.rho=FALSE){

    ## e == eigenvalues
    ## v == eigenvectors
    ## t == tau
    ## t == tau, when dC/dt = 0
    
    rho = NA
    if(cal.rho){ 
      rho = cal.rho(e=e, v=v, vinv=vinv, t=t_x)
    }
  
    N    = length(e)
    indx = seq(1,N,1)[e < 1/t_x]
    Lr   = matrix(0,N,N)
    for( i in indx ){
        Lr = Lr + e[i] * v[,i] %*% t(v[,i])
    }

    t.scaled = t/t_x

    ## To be completed...
    ##Lr = t_x * Lr
    
    return(list(rho=rho, L=Lr, t=t.scaled, n=length(indx)))
    
}    

ll.network <- function(rho, sigma){
  ## rho   == network's One propagation matrix
  ## sigma == network's Two propagation matrix
  rho[rho==0+0i]=NA
  sigma[sigma==0+0i]=NA
  
  ll = sum(diag(rho*(log(sigma,base=2))), na.rm=TRUE)
  if(is.complex(ll)){ ll=Mod(ll); }
  ll
}

kl.div.network <- function(rho, sigma){
  
  ## rho   == network's One propagation matrix
  ## sigma == network's Two propagation matrix
  rho[rho==0+0i]=NA
  sigma[sigma==0+0i]=NA

  sum(diag(rho*(log(rho,base=2) - log(sigma,base=2))),na.rm=TRUE)
  
}

js.div.network <- function(rho, sigma){
  
  ##S.p = -sum(diag(rho*(log(rho,base=2))))
  
  mu  = 0.5*(rho + sigma)
  mu[mu==0+0i]=NA
  
  kl.pq = 0.5*kl.div.network(rho=rho,   sigma=mu)
  kl.qp = 0.5*kl.div.network(rho=sigma, sigma=mu)
  
  js.div = kl.pq + kl.qp
  
  if( is.complex(js.div) ){
    js.div = Mod(js.div)  
  }
  
  js.div 
}

## calculate a * log(a/b)
log.ab <- function(a,b,SMALL=1e-20){

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
    return(sum(sapply(1:n,function(i) log.ab(p[i],q[i]))))
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
    return(H(x %*% we) - apply(x,2,H) %*% we)
    
}


js.distance <- function(x,we){

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


alpha_mle <- function(x, xmin) {
  ## Estimate alpha for power-law using MLE
  x   = as.vector(x) 
  z   = x[x>=xmin]
  n   = length(z)
  slx = sum(log(z))
  
  ## Estimate alpha using MLE
  alpha <- 1 + n * sum(slx - log(xmin - 1/2) * n)^(-1)
  
  ## estimate sd for alpha
  sd    <- (alpha-1)/sqrt(n)
  
  return(list(alpha=alpha, sd=sd))
  
}
  

degree.fit <- function(x, nsim=500, threads=4){

  x    = x[x > 0]
  m_pl = displ$new(x)
  est  = estimate_xmin(m_pl)
  m_pl$setXmin(est)
  
  gof <- bootstrap_p(m_pl, no_of_sims=nsim, threads=threads)
  
  d = plot(m_pl,draw=F)

  x.min = median(gof$bootstraps$xmin)
  x.max = median(gof$bootstraps$ntail)
  x.max = ifelse(x.min > x.max, max(gof$bootstraps$ntail), x.max)
  alpha = median(gof$bootstraps$pars)
    
  return(list(model=m_pl, d=d, x.min=x.min, x.max=x.max, alpha=alpha))
  
}

eigen.fit <- function(x, nsim=500, threads=4){
  
  ## fit power-law to eigenvalues using the poweRlaw package
  x    = x[x > 0]
  m_pl = conpl$new(x)
  est  = estimate_xmin(m_pl)
  m_pl$setXmin(est)
  
  gof <- bootstrap_p(m_pl, no_of_sims=nsim, threads=threads)
  
  d = plot(m_pl,draw=F)
  
  x.min = median(gof$bootstraps$xmin)
  x.max = median(gof$bootstraps$ntail)
  x.max = ifelse(x.min > x.max, max(gof$bootstraps$ntail), x.max)
  alpha = median(gof$bootstraps$pars)
  
  return(list(model=m_pl, d=d, x.min=x.min, x.max=x.max, alpha=alpha))
  
}

eigen.fit2 <- function(x){

## find alpha & xmin for d using poweRlaw's MLE
## create new discrete power-law distribution
x=as.numeric(x)
x=x[x>0]
X=poweRlaw::conpl$new(x)

## estimate best xmin and alpha from sequence of xmin values
X.est    = estimate_xmin(X)#, seq(0,max(x),1))
xmin     = X.est$xmin
alpha    = X.est$pars

## store node probability distribution for original network
gg.dat  = data.frame(x=X$internal$dat,
                     y=data_cdf_probs(X$internal$dat))  

return(list(df=gg.dat, xmin=xmin, alpha=alpha))

}

compute_spectral_dim <- function(x){
## https://arxiv.org/pdf/2406.19104  
## lambda_f = N^(-2/sd)
## sd = -2*log(N)/log(lambda_f)
## sd = spectral dimension of network  
## x == network eigenvalues
N = length(x)
x = as.numeric(x)
x = sort(x)
lambda_f = x[2]
-2*log(N)/log(lambda_f)
}

# Step 4: Define C(τ_F(α)) using the Gamma function
C_tau <- function(alpha, ds) {
  gamma1 <- gamma(ds/2)
  gamma2 <- gammainc(ds/2 + 2, alpha)[1]
  gamma3 <- gammainc(ds/2 + 1, alpha)[1]
  
  (gamma1 * gamma2 - gamma3^2) / (gamma1^2)
}


log.log_plot <- function(df.x, df.y, x.lab="K", y.lab="LRG.K"){

colours        = c("red","blue")
names(colours) = c(x.lab, y.lab)
  
## Log-Log plot of original and coarse-grain graph's degree distributions
gp = ggplot()+
  geom_line(data=df.x, aes(log(x), log(y), color=x.lab), linewidth=1)+
  geom_line(data=df.y, aes(log(x), log(y), color=y.lab), linewidth=1)+
  labs(x = "log(K)", y = "P(K)") +
  ##scale_y_continuous(labels = scales::dollar) +
  scale_color_manual(values=colours, name="")+
  theme_light() +
  theme(legend.position = "bottom")

  return(gp)
}


scan.time <- function(gg, e=NULL, v=NULL, vinv=NULL, L=NULL, t.lower, t.upper, nsim=500, complex=FALSE, 
                      n.steps=10, method=c("eigen", "balanced", "square"), 
                      expm_method=c("Higham08.b"), tol=1e-5, order=1, print=1){
  
  method <- match.arg(method)
  
  ## define time
  ##n.steps = 10
  gap     = t.upper - t.lower 
  steps   = gap/n.steps
  time    = seq(t.lower, t.upper, steps)
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
  #X$setXmin(xmin)
  #X$setPars(alpha)
  
  ## store node probability distribution for original network
  gg.dat  = data.frame(x=X$internal$values,
                       y=data_cdf_probs(X$internal$values))  
  
  plots  = list()
  plots2 = list()
  
  
  for( t in 1:length(time) ){
    
    ## Perform coarse graining at time: t
    lrg   = real.LRG(e=e, v=v, vinv=vinv, L=L, t=time[t], gg=gg, 
                     complex=complex, method=method, expm_method=expm_method, 
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
    #lrg.X$setXmin(lrg.xmin)
    #lrg.X$setPars(lrg.alpha)
    
    
    ## store node probability distribution for coarse-grain network
    lrg.dat  = data.frame(x=lrg.X$internal$values,
                          y=data_cdf_probs(lrg.X$internal$values)) 
    
    
    #Cn      = table(as.numeric(lrg$mapping[,2]))
    #Cn.gt1  = sum(as.numeric(names(Cn))>1) ##length(Cn[Cn>1])
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
      
      #p = gg.dat$y
      #q = pdf_powerlaw(x=gg.dat$x, lrg.xmin, lrg.alpha)
      #js= js.distance(x=cbind(p,q),we=rep(0.5,2))
      
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
                      "t"=time[t])
      
      k=k+1
      
      ## Print progress
      if(print){cat(sprintf("Iteration: %d, ks_dist: %f, jsd: %f, loss: %f, lrg.sn: %d, alpha: %f, alpha.lrg: %f, xmin: %d, xmin.lrg: %d, t: %f\n", 
                            t, ks_dist[[1]], js$jsd, loss, length(lrg.sn), alpha, lrg.alpha, xmin, lrg.xmin, time[t]))}
      
      
      ## Log-Log plot of original and coarse-grain graph's degree distributions
      plots[[t]] = log.log_plot(df.x=gg.dat, 
                                df.y=lrg.dat, 
                                x.lab="K", y.lab="LRG.K")
      
      
      ## Log-Log plot of original and coarse-grain graph's degree distributions
      plots2[[t]] = log.log_plot(df.x=data.frame(x=x.range,
                                                 y=gg_pdf), 
                                 df.y=data.frame(x=x.range,
                                                 y=lrg_pdf), 
                                 x.lab="K", y.lab="LRG.K")
      }
    }
  }
}

  df = data.frame(do.call(rbind, lapply(res, unlist)))
    
  return(list(df=df,plots=plots, plots2=plots2))
  
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
  constant = pracma::zeta(alpha)
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
  constant = pracma::zeta(alpha)
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

#ks_dist <- function(x, x2, xmin1, alpha1, xmin2, alpha2){
#  cdf_emp = ccdf_powerlaw(x, xmin1, alpha1)
#  cdf_fit = ccdf_powerlaw(x2, xmin2, alpha2)
#  n.min   = min(length(cdf_emp), length(cdf_fit))
#  cdf_emp = cdf_emp[1:n.min]
#  cdf_fit = cdf_fit[1:n.min]
#  ks_dist = max(abs(cdf_emp-cdf_fit)) 
#  return(list(ks_dist=ks_dist, cdf_emp=cdf_emp, cdf_fit=cdf_fit))
#}

#ks_dist <- function(x, x2){
#  cdf_emp = data_cdf_probs(x)
#  cdf_fit = data_cdf_probs(x2)
#  n.min   = min(length(cdf_emp), length(cdf_fit))
# cdf_emp = cdf_emp[1:n.min]
#  cdf_fit = cdf_fit[1:n.min]
#  ks_dist = max(abs(cdf_emp-cdf_fit)) 
#  return(list(ks_dist=ks_dist, cdf_emp=cdf_emp, cdf_fit=cdf_fit))
#}

# Simulated Annealing with Scaling and Squaring
anneal.time <- function(gg, e=NULL, v=NULL, vinv=NULL, L=NULL, t.lower, t.upper, complex=FALSE, 
                        n.steps=10, method=c("eigen", "balanced", "square"), length.out=150,
                        expm_method=c("Higham08.b"), tol=1e-5, order=1, print=1,
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
    times = t
    
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
    #old_loss    = js.distance(we=c(1/2,1/2), x=cbind(gg_pdf, lrg_pdf))[[4]]
    #old_loss = as.vector(ks.test(gg_cdf, lrg_cdf)[[1]])
    #old_loss    = cross.entropy.loss(gg.dat$y, lrg.dat$y)
    #old_loss    = cross.entropy.loss(gg_pdf, lrg_pdf)
    #loss        = old_loss
      curr_t      = t
      curr_loss   = loss
      deltaL      = NA
      
    res[[k]]    = c("Iter"=1, "Loss(old)"=NA, "Loss(new)"=loss, 
                    #"Loss(delta)"=NA, "Prob:"=NA, 
                    "alpha_gg"=alpha, "alpha_lrg"=lrg.alpha, 
                    "xmin_gg"=xmin, "xmin_lrg"=lrg.xmin,
                    "t"=NA, "new_t"=t)#, "dt"=dt)
    
    plots[[k]]  = log.log_plot(df.x=gg.dat, 
                               df.y=lrg.dat, 
                               x.lab="K", y.lab="LRG.K")
    
    plots2[[k]] = log.log_plot(df.x=data.frame(x=x.range,
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
    while( new_t < 0 ){
      new_t = curr_t + rnorm(n=1, mean=dt)
    }
    
    ## add time to set of times
    times = append(times, new_t)
    
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
          
          plots[[k]]  = log.log_plot(df.x=gg.dat, 
                                     df.y=lrg.dat, 
                                     x.lab="K", y.lab="LRG.K")
          
          plots2[[k]] = log.log_plot(df.x=data.frame(x=x.range,
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
      #if(print){
      #  cat(sprintf("Iteration: %d, Loss(old): %f, Loss(new): %f, Loss(delta): %f, Prob: %f, t: %f, t(new): %f, dt: %f\n", 
      #              iter, curr_loss, NA, NA, NA, curr_t, new_t, dt))}
    }
  }
}
    
  if( !isempty(res) ){  
    df = data.frame(do.call(rbind, lapply(res, unlist)))
    colnames(df) = c("Iter", "Loss.old", "Loss.new",
                     "alpha_gg", "alpha_lrg", "xmin_gg", "xmin_lrg", "t", "new_t")##, "dt")  
  }
    
  return(list(df=df, times=times, plots=plots, plots2=plots2, restarts=restarts))
  
}


eigen_loss <- function(t, gg, L, e, v, vinv, complex=FALSE, 
                       method=c("eigen", "balanced", "square"), mle=0,
                       loss_test=c("ks_dist", "kl_div", "js_dist"),
                       expm_method=c("Higham08.b"), tol=1e-5, order=1, plots=0) {
  
  method    <- match.arg(method)
  loss_test <- match.arg(loss_test)
  plts      <- NULL
  data      <- NULL
  params    <- NULL
  ks_result <- NULL
  losses    <- NULL
  df.x      <- NULL
  df.y      <- NULL
  
   tryCatch({
     
    N = dim(L)[1]
     
    ## Perform coarse graining at time: new_t
    lrg   = real.LRG(e=e, v=v, vinv=vinv, L=L, t=t, gg=gg, 
                     complex=complex, method=method, expm_method=expm_method, 
                     tol=tol, order=order)

    clusters = as.numeric(names(table(lrg$supernodes)))
    if( sum(!clusters %in% c(1,N)) == 0 ) break;
    
    ## Eigenvalue distribution of the original network
    e.sd  <- compute_spectral_dim(e)
    e.cut <- e[e>0]
    X     <- poweRlaw::conpl$new(e.cut)
    X.est <- poweRlaw::estimate_xmin(X)
    xmin  <- X.est$xmin
    zl    <- sort(X$internal$dat)
    zl    <- zl[zl >= xmin]
    alpha <- X.est$pars
    alpha_fit <- alpha_mle(x=e, xmin=xmin)
    
    ## Degree distribution of the coarse-grain network
    loss <- tryCatch({
      lrg.L     <- get.L(lrg$gg)
      lrg.Eig   <- get.eigen(lrg.L, only.values=TRUE)
      lrg.E     <- lrg.Eig[[1]]
      lrg.sd    <- compute_spectral_dim(lrg.E)
      lrg.E.cut <- lrg.E[lrg.E>0]
      lrg.X     <- poweRlaw::conpl$new(lrg.E.cut)
      lrg.est   <- poweRlaw::estimate_xmin(lrg.X)
      lrg.xmin  <- lrg.est$xmin
      lrg.alpha <- lrg.est$pars
      lrg.alpha_fit <- alpha_mle(x=lrg.E, xmin=lrg.xmin)
      
      lrg.zl = sort(lrg.X$internal$dat)
      lrg.zl = lrg.zl[lrg.zl >= lrg.xmin]
            
      ## set alpha & xmin values to use in study
      if( mle==1 ){
        alpha_study     = alpha_fit[[1]]
        lrg.alpha_study = lrg.alpha_fit[[1]]
      } else {
        alpha_study     = alpha
        lrg.alpha_study = lrg.alpha
      }
      
      ## Compute C0 from alpha      
      C0      <- alpha_study+1
      lrg.C0  <- lrg.alpha_study+1
      
      ## Calculate losses 
      #ks_result <- ks_dist(x=zl, x2=lrg.zl)
      #ks_result <- ks_dist(x=zl, x2=lrg.zl, xmin1=xmin, alpha1=alpha_study, xmin2=lrg.xmin, alpha2=lrg.alpha_study)
      ks_result <- ks_dist(x=zl, xmin=lrg.xmin, alpha=lrg.alpha_study)
      kl_result <- kl.divergence(p=ks_result[[2]], q=ks_result[[3]])
      js_dist   <- js.distance(we = c(1/2, 1/2), x = cbind(p = ks_result[[2]], q = ks_result[[3]]))##[[4]]
      losses    <- list(ks_dist=ks_result[[1]], kl=kl_result, js_dist=js_dist)
      
      ## Select loss 
      loss      <- switch(loss_test,
                          "ks_dist"= ks_result[[1]],
                          "kl_div" = kl_result[[1]],
                          "js_dist"= js_dist[[4]],
                          Inf)
      
      if( plots ){
        
        data = data.frame(x=X$internal$dat, 
                          y=data_cdf_probs(X$internal$dat))
        
        #zl.min = length(ks_result[[2]])
        #
        #df.x = data.frame(x=zl[1:zl.min], y=ks_result[[2]])
        #df.y = data.frame(x=lrg.zl[1:zl.min], y=ks_result[[3]])
        
        df.x = data.frame(x=zl, y=ks_result[[2]])
        df.y = data.frame(x=zl, y=ks_result[[3]])
        
        plts[[1]] = log.log_plot(df.x=df.x, df.y=df.y, 
                                 x.lab="E(emp.)", y.lab="LRG.E(fit)")
        
        plts[[2]] = lambda_fit_plot2(df=data, 
                                     xmin1=xmin,     alpha1=alpha_study, 
                                     xmin2=lrg.xmin, alpha2=lrg.alpha_study)

        params = list(x.fit=X.est, lrg.fit=lrg.est, 
                      spec_dim=e.sd, spec_dim_lrg=lrg.sd,
                      C0=C0, lrg.C0=lrg.C0, mle_flag=mle,
                      alpha_fit=alpha_fit, lrg.alpha_fit=lrg.alpha_fit)
      }
    
      return(list(loss=loss, plots=plts, params=params, df.x=df.x, df.y=df.y,
                  ks_result=ks_result, data=data, losses=losses))
    }, 
   
     error = function(e){
      Inf
      }
    )
    
    if (is.na(loss) || is.nan(loss)) loss <- Inf
    return(list(loss=loss, plots=plts, params=params, df.x=df.x, df.y=df.y,
                ks_result=ks_result, data=data, losses=losses))
    
  }, error = function(e) {
    # Return NA or a default value on error
    warning(paste("Error in compute_loss at t =", t, ":", e$message))
    return(list(loss=Inf, plots=plts, params=params, df.x=df.x, df.y=df.y,
                ks_result=ks_result, data=data, losses=losses))
  })
}

## Define a wrapper for parallel evaluation of the loss function
parallel_eigen_loss <- function(t, gg, L, e, v, vinv, complex=FALSE, 
                              method=c("eigen", "balanced", "square"), 
                              loss_test=c("ks_dist", "kl_div", "js_dist"),
                              expm_method="Higham08.b", mle=0,
                              tol=1e-5, order=1, plots=0, cl) {
  
  ## Load packages on each worker
  clusterEvalQ(cl, {
    library(poweRlaw)
    library(igraph)
    library(expm)
    library(pracma)
  })
  
  clusterExport(cl, list("real.LRG", "ks_dist", "js.distance", "js.divergence",
                         "H", "cal.rho.approx", "taylor.approx", "kl.divergence", 
                         "cal.rho.eigen", "cal.rho.balanced", "log.ab",
                         "cal.rho.square","meta.graph", "meta.binary",                                  "meta.edge", "get.supernodes", 
                         "coarse.grain.graph", "ccdf_powerlaw",
                         "cdf_powerlaw", "pdf_powerlaw", "data_cdf_probs",
                         "get.L", "get.eigen", "compute_spectral_dim", "loss_test",
                         "t", "gg", "L", "e", "v", "vinv", "mle", "alpha_mle",
                         "complex", "method", "expm_method", "tol", "plots",
                         "order", "eigen_loss"), envir = environment())
  
  # Define the computation for each particle's position
  compute_loss <- function(x) {
    eigen_loss(t=x, gg=gg, L=L, e=e, v=v, vinv=vinv, 
               complex=complex, method=method, loss_test=loss_test, mle=mle,
               expm_method=expm_method, tol=tol, order=order, plots=plots)[[1]]
  }
  
  # Parallel computation of the loss for `t`
  losses <- parSapply(cl, t, compute_loss)
  return(losses)
}


iterative_eigen_optimization <- function(t_lower, t_upper, 
                                         gg, L, e, v, vinv, 
                                         complex=FALSE, method=c("eigen", "balanced", "square"),
                                         loss_test=c("ks_dist", "kl_div", "js_dist"),
                                         expm_method="Higham08.b", tol=1e-5, 
                                         order=1, cl, max_iter=10, n_points=10, 
                                         dt=1, cooling_rate=0.99,
                                         refine_factor=0.5, mle=0,
                                         print=1, plots=0) {
  method <- match.arg(method)
  best_t <- NULL
  best_loss <- Inf  # Start with a very large loss value
  res    <- list()
  k      <- 1
  
  ## Load packages on each worker
  clusterEvalQ(cl, {
    library(poweRlaw)
    library(igraph)
    library(expm)
    library(pracma)
  })
  
  clusterExport(cl, list("real.LRG", "ks_dist", "js.distance", "js.divergence",
                         "H", "cal.rho.approx", "taylor.approx", "kl.divergence",
                         "cal.rho.eigen", "cal.rho.balanced", "log.ab",
                         "cal.rho.square","meta.graph", "meta.binary",                                  "meta.edge", "get.supernodes", 
                         "coarse.grain.graph", "ccdf_powerlaw",
                         "cdf_powerlaw", "pdf_powerlaw", "data_cdf_probs",
                         "get.L", "get.eigen", "compute_spectral_dim", "loss_test",
                         "t", "gg", "L", "e", "v", "vinv", "mle", "alpha_mle",
                         "complex", "method", "expm_method", "tol", "plots", 
                         "order", "eigen_loss"), envir = environment())
  
  
  compute_loss <- function(x) {
    eigen_loss(t=x, gg=gg, L=L, e=e, v=v, vinv=vinv, 
               complex=complex, method=method, loss_test=loss_test, mle=mle,
               expm_method=expm_method, tol=tol, order=order, plots=plots)[[1]]
  }
  
  ## Generate current state
  n_states = n_points + floor(n_points*0.5)
  t        = t_lower + runif(n_states)*(t_upper - t_lower)
  cont     = FALSE 
  
  while ( !cont ){
    
    ## Compute losses in parallel
    losses <- parSapply(cl, t, compute_loss)
    
    ## Remove NA values
    valid_losses <- losses[!is.na(losses)]
    valid_t      <- t[!is.na(losses)]
    
    ## Remove Inf values
    valid_losses <- valid_losses[!is.infinite(valid_losses)]
    valid_t      <- valid_t[!is.infinite(valid_losses)]
    
    if (length(valid_losses) > 0) {
      cont = TRUE
    } else {
      t    = t_lower + runif(n_states)*(t_upper - t_lower)
    }
    
  }
  
  ## Find the best `t` in this iteration
  curr_t    <- valid_t[which.min(valid_losses)]
  curr_loss <- min(valid_losses)
  deltaL    <- NA
  rm(losses, valid_losses, valid_t)
  
  if( print ){ cat("Starting t:", curr_t, ", and loss", curr_loss,"\n") }
  
  for (iter in seq_len(max_iter)) {
    if( print ){ cat("Iteration:", iter, "\n") }
    
    # Generate a vector of `t` values within the current range
    ##t     = t.lower[1] + runif(1) * (t.upper - t.lower)
    new_t = t_lower + rnorm(n=n_points, sd=curr_t) * dt
    for( i in 1:n_points ){
      while( new_t[i] < t_lower | new_t[i] > t_upper ){ 
        new_t[i] = t_lower + rnorm(n=1, sd=curr_t) * dt 
      }
    }
    
    if( print ){ cat("new_t:", new_t, "\n") }
    
    # Compute losses in parallel
    losses <- parSapply(cl, new_t, compute_loss)
    
    if( print ){ cat("losses:", losses, "\n") }
    
    ## Remove NA values
    valid_losses <- losses[!is.na(losses)]
    valid_t      <- new_t[!is.na(losses)]
    
    ## Remove Inf values
    valid_losses <- valid_losses[!is.infinite(valid_losses)]
    valid_t      <- valid_t[!is.infinite(valid_losses)]
    
    if (length(valid_losses) == 0) {
      if( print ){ cat("No valid losses found in iteration", iter, "\n") }
      break
    }
    
    ## Find the best `t` in this iteration
    iter_best_loss <- min(valid_losses)
    iter_best_t    <- valid_t[which.min(valid_losses)]
    
    if( print ){ cat("Best t in this iteration:", iter_best_t, "with loss:", iter_best_loss, "\n") }
    
    # Update the global best `t` and loss if the current iteration is better
    if (iter_best_loss < best_loss) {
      
      res[[k]]  <- c("Iter"=iter, "t_old"=best_t, "t_new"=iter_best_t, 
                     "loss_old"=best_loss, "loss_new"=iter_best_loss,
                     "t_lower"=t_lower, "t_upper"=t_upper)
      
      best_loss <- iter_best_loss
      best_t    <- iter_best_t
      k         <- k + 1
    }
    
    rm(losses, valid_losses, valid_t)
    
    ## difference between new and old losses
    deltaL = (iter_best_loss-curr_loss)
    
    ## Cooling schedule: Update t
    dt = dt * cooling_rate
    
    ## Calculate the acceptance probability
    acceptance_prob <- exp(-abs(deltaL)/dt)
    
    ## Accept the new state with a certain probability
    if ( (deltaL < 0) || (runif(1) < acceptance_prob)) {
      
      ## update current time and loss
      curr_t = iter_best_t; curr_loss = iter_best_loss;
      
    }
    
    # Refine the search range
    #range_width <- (t_upper - t_lower) * refine_factor
    #t_lower <- max(t_lower, iter_best_t - range_width / 2)
    #t_upper <- min(t_upper, iter_best_t + range_width / 2)
    
    # Check for convergence (range too small)
    if ( abs(deltaL) < tol) {
      if( print ){ cat("Converged with deltaL", deltaL, "\n") }#range [", t_lower, ",", t_upper, "]\n")
      break
    }
  }
  
  df = data.frame(do.call(rbind, lapply(res, unlist)))
  
  # Return the best `t` and loss
  list(best_t = best_t, best_loss = best_loss, df=df)
}



## General Degree loss function
degree_loss <- function(t, gg, L, e, v, vinv, complex=FALSE, 
                     method=c("eigen", "balanced", "square"),
                     expm_method=c("Higham08.b"), tol=1e-5, order=1) {
  
  method <- match.arg(method)
    
  tryCatch({
    ## Perform coarse graining at time: new_t
    lrg   = real.LRG(e=e, v=v, vinv=vinv, L=L, t=t, gg=gg, 
                     complex=complex, method=method, expm_method=expm_method, 
                     tol=tol, order=order)
    
    # Degree distribution of the original network
    d     <- as.numeric(igraph::degree(gg))
    X     <- poweRlaw::displ$new(d)
    X.est <- poweRlaw::estimate_xmin(X, seq(1, max(d), 1))
    xmin  <- X.est$xmin
    zl    <- sort(X$dat)
    zl    <- zl[zl >= xmin]
    alpha <- X.est$pars
    
    ## Degree distribution of the coarse-grain network
    lrg.d     <- as.numeric(igraph::degree(lrg$gg))
    lrg.X     <- poweRlaw::displ$new(lrg.d)
    lrg.est   <- poweRlaw::estimate_xmin(lrg.X, seq(1, max(lrg.d), 1))
    lrg.xmin  <- lrg.est$xmin
    lrg.alpha <- lrg.est$pars
    
    ## Calculate Jensen-Shannon distance
    ks_result <- ks_dist(x = zl, xmin = lrg.xmin, alpha = lrg.alpha)
    loss      <- js.distance(we = c(1/2, 1/2), x = cbind(p = ks_result[[2]], q = ks_result[[3]]))[[4]]
    
    if (is.na(loss) || is.nan(loss)) loss <- Inf
    return(loss)
    
  }, error = function(e) {
    # Return NA or a default value on error
    warning(paste("Error in compute_loss at t =", t, ":", e$message))
    return(Inf)
  })
  }

# Define a wrapper for parallel evaluation of the loss function
parallel_degree_loss <- function(t, gg, L, e, v, vinv, complex=FALSE, 
                                 method=c("eigen", "balanced", "square"), expm_method="Higham08.b", 
                                 tol=1e-5, order=1, cl) {
  
  ## Load packages on each worker
  clusterEvalQ(cl, {
    library(poweRlaw)
    library(igraph)
    library(expm)
    library(pracma)
  })
  
  clusterExport(cl, list("real.LRG", "ks_dist", "js.distance", "js.divergence",
                        "H", "cal.rho.approx", "taylor.approx", 
                        "cal.rho.eigen", "cal.rho.balanced", 
                        "cal.rho.square","meta.graph", "meta.binary",                                  "meta.edge", "get.supernodes", 
                        "coarse.grain.graph", "ccdf_powerlaw",
                        "cdf_powerlaw", "pdf_powerlaw", "t",
                        "gg", "L", "e", "v", "vinv",
                        "complex", "method", "expm_method", "tol",
                        "order", "pso_loss"), envir = environment())
  
  # Define the computation for each particle's position
  compute_loss <- function(x) {
    degree_loss(t=x, gg=gg, L=L, e=e, v=v, vinv=vinv, 
                complex=complex, method=method, 
                expm_method=expm_method, tol=tol, order=order)
  }
  
  # Parallel computation of the loss for `t`
  losses <- parSapply(cl, t, compute_loss)
  return(losses)
}


iterative_degree_optimization <- function(t_lower, t_upper, 
                                          gg, L, e, v, vinv, 
                                          complex=FALSE, method=c("eigen", "balanced", "square"),
                                          expm_method="Higham08.b", tol=1e-5, 
                                          order=1, cl, max_iter=10, n_points=10, 
                                          dt=1, cooling_rate=0.99,
                                          refine_factor=0.5) {
  method <- match.arg(method)
  best_t <- NULL
  best_loss <- Inf  # Start with a very large loss value
  res    <- list()
  k      <- 1
  
  ## Load packages on each worker
  clusterEvalQ(cl, {
    library(poweRlaw)
    library(igraph)
    library(expm)
    library(pracma)
  })
  
  clusterExport(cl, list("real.LRG", "ks_dist", "js.distance", "js.divergence",
                         "H", "cal.rho.approx", "taylor.approx", 
                         "cal.rho.eigen", "cal.rho.balanced", 
                         "cal.rho.square","meta.graph", "meta.binary",                                  "meta.edge", "get.supernodes", 
                         "coarse.grain.graph", "ccdf_powerlaw",
                         "cdf_powerlaw", "pdf_powerlaw", "t",
                         "gg", "L", "e", "v", "vinv",
                         "complex", "method", "expm_method", "tol",
                         "order", "pso_loss"), envir = environment())
  
  
  compute_loss <- function(x) {
    degree_loss(t=x, gg=gg, L=L, e=e, v=v, vinv=vinv, 
                complex=complex, method=method, 
                expm_method=expm_method, tol=tol, order=order)
  }
  
  ## Generate current state
  n_states = n_points + floor(n_points*0.5)
  t    = t_lower + runif(n_states)*(t_upper - t_lower)
  cont = FALSE 
  
  while ( !cont ){
  
    ## Compute losses in parallel
    losses <- parSapply(cl, t, compute_loss)
    
    ## Remove NA values
    valid_losses <- losses[!is.na(losses)]
    valid_t      <- t[!is.na(losses)]
    
    ## Remove Inf values
    valid_losses <- valid_losses[!is.infinite(valid_losses)]
    valid_t      <- valid_t[!is.infinite(valid_losses)]
    
    if (length(valid_losses) > 0) {
      cont = TRUE
    } else {
      t    = t_lower + runif(n_states)*(t_upper - t_lower)
    }

  }
        
  ## Find the best `t` in this iteration
  curr_t    <- valid_t[which.min(valid_losses)]
  curr_loss <- min(valid_losses)
  deltaL    <- NA
  rm(losses, valid_losses, valid_t)
  
  cat("Starting t:", curr_t, ", and loss", curr_loss,"\n")
  
  for (iter in seq_len(max_iter)) {
    cat("Iteration:", iter, "\n")
    
    # Generate a vector of `t` values within the current range
    ##t     = t.lower[1] + runif(1) * (t.upper - t.lower)
    new_t = t_lower + rnorm(n=n_points, sd=curr_t) * dt
    for( i in 1:n_points ){
      while( new_t[i] < t_lower | new_t[i] > t_upper ){ 
        new_t[i] = t_lower + rnorm(n=1, sd=curr_t) * dt 
        }
    }
    
    cat("new_t:", new_t, "\n")
    
    # Compute losses in parallel
    losses <- parSapply(cl, new_t, compute_loss)
  
    cat("losses:", losses, "\n")
    
    ## Remove NA values
    valid_losses <- losses[!is.na(losses)]
    valid_t      <- new_t[!is.na(losses)]
    
    ## Remove Inf values
    valid_losses <- valid_losses[!is.infinite(valid_losses)]
    valid_t      <- valid_t[!is.infinite(valid_losses)]
    
    if (length(valid_losses) == 0) {
      cat("No valid losses found in iteration", iter, "\n")
      break
    }
    
    ## Find the best `t` in this iteration
    iter_best_loss <- min(valid_losses)
    iter_best_t    <- valid_t[which.min(valid_losses)]
    
    cat("Best t in this iteration:", iter_best_t, "with loss:", iter_best_loss, "\n")
    
    # Update the global best `t` and loss if the current iteration is better
    if (iter_best_loss < best_loss) {
      
      res[[k]]  <- c("Iter"=iter, "t(old)"=best_t, "t(new)"=iter_best_t, 
                      "loss(old)"=best_loss, "loss(new)"=iter_best_loss,
                      "t_lower"=t_lower, "t_upper"=t_upper)
      
      best_loss <- iter_best_loss
      best_t    <- iter_best_t
      k         <- k + 1
    }
    
    rm(losses, valid_losses, valid_t)
    
    ## difference between new and old losses
    deltaL = (iter_best_loss-curr_loss)
    
    ## Cooling schedule: Update t
    dt = dt * cooling_rate
    
    ## Calculate the acceptance probability
    acceptance_prob <- exp(-abs(deltaL)/dt)
    
    ## Accept the new state with a certain probability
    if ( (deltaL < 0) || (runif(1) < acceptance_prob)) {
      
      ## update current time and loss
      curr_t = iter_best_t; curr_loss = iter_best_loss;
      
    }
    
    # Refine the search range
    #range_width <- (t_upper - t_lower) * refine_factor
    #t_lower <- max(t_lower, iter_best_t - range_width / 2)
    #t_upper <- min(t_upper, iter_best_t + range_width / 2)
    
    # Check for convergence (range too small)
    if ( abs(deltaL) < tol) {
      cat("Converged with deltaL", deltaL, "\n")#range [", t_lower, ",", t_upper, "]\n")
      break
    }
  }
  
  df = data.frame(do.call(rbind, lapply(res, unlist)))
  
  # Return the best `t` and loss
  list(best_t = best_t, best_loss = best_loss, df=df)
}


## Function to aggregate nodes based on rho values
aggregate_nodes <- function(rho, N, n) {
  
  ## Step 1: Calculate absolute values of rho
  abs_rho_values <- abs(rho)
  
  ## Step 2: Order absolute values in descending order
  ordered_values <- sort(abs_rho_values, decreasing = TRUE)
  
  ## Step 3: Start with each node as its own cluster
  clusters <- lapply(1:N, function(i) { i })
  
  ## Step 4: Merge clusters until desired number of clusters is reached
  for (value in ordered_values) {
    if (length(clusters) == N - n) {
      break  ## Stop when desired number of clusters is reached
    }
    for (i in 1:N) {
      for (j in i:N) {
        if (abs_rho_values[i, j] == value) {
          ## Merge clusters containing nodes i and j
          for (k in seq_along(clusters)) {
            if (i %in% clusters[[k]]) {
              if (j <= length(clusters)) {
              clusters[[k]] <- union(clusters[[k]], clusters[[j]])
              clusters <- clusters[-j]
              break
            }
          }
        }
      }
    }
    }
  }
  
  return(clusters)
}

bin.rule <- function(x){
 ## Freedman-Diaconis rule
 ## Step 2: Calculate the IQR
  iqr <- IQR(x)
  
  # Step 3: Determine the Bin Width
  n <- length(x)
  bin_width <- 2 * iqr / (n^(1/3))
  
  # Step 4: Calculate the Number of Bins
  range_data <- range(x)
  num_bins   <- ceiling((range_data[2] - range_data[1]) / bin_width)
  return(num_bins)
}


bin.eigen <- function(x, steps=0.05, x.min=1e-4, x.max=1e2, group="A"){
  x  = x[x>0]
  x  = x[order(x)]
 
  indx = which(x<1e-3)
  if( length(indx) > 1 ){
    x.mid = x[max(indx)]
    mid.steps = steps*1e-1
    s0 = seq(x.min, x.mid, mid.steps)
    s1 = seq(x.mid, x.max, steps)
    s  = c(s0,s1)  
  } else {
    s  = seq(x.min, x.max, steps)
  }

  ##s  = seq(x.min, x.max, steps)
  
  c  = cut(x, s, labels=FALSE, include.lowest=TRUE)
  names(s)=seq(1,length(s),1)
  p  = table(c)+10^-30
  x.new = as.vector(s[match(names(p), names(s))])
  y  = as.vector(p/sum(p))
  ##y  = data_cdf_probs(p) 
  n  = length(y)
  df = data.frame(x=x.new,y=y,group=rep(group,n))
}

smooth.bin.eigen <- function(df, method="box", bandwidth=7){
  group  = levels(factor(df$group))
  method = match.arg(method, c("box","normal"))
  if( method == "box" ){
    fit <- with(df, ksmooth(x, y, kernel = "box", bandwidth = bandwidth))
  } else {
    fit <- with(df, ksmooth(x, y, kernel = "normal", bandwidth = bandwidth))
  }
  
  n = length(fit$x)
  fit$group = rep(group,n)
  
  return(list2DF(fit))
  
}

## x = vector of eigenvalues of Laplacian
## return plot of log(P(x)) Versus log(x)
lambda_plot <- function(df, steps=0.05, bandwidth=7){

  #df  = bin.eigen(x,steps=steps)  
  #fit = smooth.bin.eigen(df, bandwidth=bandwidth)
  
  #fit1 <- with(df, ksmooth(x, y, kernel = "box", bandwidth = bandwidth))
  #fit2 <- with(df, ksmooth(x, y, kernel = "normal", bandwidth = bandwidth))
  
  pt = ggplot(df, aes(x,y))+
    geom_point(color="red")+
    geom_line(linewidth=1, color="red")+
    scale_x_log10(
      name = TeX("$\\lambda$"),
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x)))+
    scale_y_log10(
      name = TeX("P($\\lambda$)"),
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x)))+
    annotation_logticks(sides="b", outside=TRUE)+
    coord_cartesian(clip = "off")+
    theme_light()+
    theme(legend.position = "right")
  
  return(pt)
  
}

compute_fit_values <- function(df, xmin, alpha){

  ## Filter data to get points for the straight line
  df_line <- df %>%
    filter(x >= xmin) %>%
    arrange(x)
  
  ## Calculate the straight line values for all x >= xmin
  df_line <- df_line %>%
    mutate(
      y_line = ccdf_powerlaw(x, xmin, alpha) 
    )

}

lambda_fit_plot <- function(df, xmin, alpha, fit.color="blue"){
  
  ## Calculate the straight line values for all x >= xmin
  fit_line = compute_fit_values(df, xmin, alpha)
  
  pt = ggplot(df, aes(x,y))+
    geom_point(color="red")+
    geom_line(linewidth=1, color="red")+
    geom_line(data = fit_line, aes(x, y_line), color = fit.color, linewidth = 1) +
    scale_x_log10(
      name = TeX("$\\lambda$"),
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x)))+
    scale_y_log10(
      name = TeX("P($\\lambda$)"),
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x)))+
    annotation_logticks(sides="b", outside=TRUE)+
    coord_cartesian(clip = "off")+
    theme_light()+
    theme(legend.position = "right")
  
  return(pt)
  
}

lambda_fit_plot2 <- function(df, xmin1, alpha1, fit.color1="red", fit.name1="E",
                                 xmin2, alpha2, fit.color2="blue", fit.name2="LRG.E"){
  
  colours        = c(fit.color1, fit.color2)
  names(colours) = c(fit.name1, fit.name2)
  
  ## Calculate the straight line values for all x >= xmin
  fit1 = compute_fit_values(df, xmin1, alpha1)
  fit2 = compute_fit_values(df, xmin2, alpha2)
  
  #browser()
  
  pt = ggplot(df, aes(x,y))+
    geom_point(color="black")+
    geom_line(linewidth=1, color="black")+
    geom_line(data = fit1, aes(x, y_line, color = fit.name1), linewidth = 1) +
    geom_line(data = fit2, aes(x, y_line, color = fit.name2), linewidth = 1) +
    scale_x_log10(
      name = TeX("$\\lambda$"),
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x)))+
    scale_y_log10(
      name = TeX("P($\\lambda$)"),
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x)))+
    scale_color_manual(values=colours, name="")+
    annotation_logticks(sides="b", outside=TRUE)+
    coord_cartesian(clip = "off")+
    theme_light()+
    theme(legend.position = "bottom")
  
  return(pt)
  
}


turing.points <- function(x, tp.range=2, print=1){
  ## Find turning-point using rBeast
  ## get y values from plot
  out   = NULL
  reg.x = NULL
  if( is.ggplot(x) ){
    pt  = ggplot_build(x)
    X   = pt$data[[1]]$x
    Y   = pt$data[[1]]$y
    c.max = max(Y)
    out = beast(Y, season="none", 
                print.options=FALSE, 
                print.progress=FALSE, 
                quiet=TRUE) 
  
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
    tp.stats[,3] = 10^(X[tp])
  
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
      tp.stats[i,6] = 10^(X[tp.min])
      tp.stats[i,7] = tp.max
      tp.stats[i,8] = 10^(X[tp.max])
      tp.stats[i,9] = tp.cum_pr
      if(print){ cat("tp: ", tp.i, ", Sum Pr: ", tp.cum_pr, "\n") }
    }
  
  
    ## turning point round region of interest, ordered first to last
    reg.x = tp.stats[order(tp.stats[,1]),]
    reg.x = as.data.frame(reg.x)
  }
    
  return(list(pt=out, tp=reg.x))
  
}

supernodes_from_rho <- function(rho, n=1){ 
 
  # Compute absolute values of rho_ij(tau*)
  rho_values <- abs(rho)
  
  # Get number of nodes (N)
  N <- nrow(rho_values)
  
  # Flatten the upper triangle of rho_values (excluding diagonal)
  rho_flat <- as.vector(rho_values[upper.tri(rho_values, diag = FALSE)])
  
  # Order indices based on descending values of rho_ij(tau*)
  ordered_indices <- order(rho_flat, decreasing = TRUE)
  
  # Initialize clusters/supernodes
  supernodes <- list()
  for (i in 1:N) {
    supernodes[[i]] <- i  # Each node starts as its own supernode
  }
  
  # Merge nodes based on ordered values until N - n clusters are obtained
  while (length(supernodes) > N - n) {
    # Get the indices of the next highest rho_ij value
    next_index <- ordered_indices[[1]]
    
    # Find the nodes i and j corresponding to next_index
    i <- (next_index - 1) %/% N + 1
    j <- (next_index - 1) %% N + 1
    
    # Find which supernodes contain nodes i and j
    supernode_i <- which(sapply(supernodes, function(s) i %in% s))
    supernode_j <- which(sapply(supernodes, function(s) j %in% s))
    
    # Merge supernodes if they are different
    if (supernode_i != supernode_j) {
      supernodes[[supernode_i]] <- c(supernodes[[supernode_i]], supernodes[[supernode_j]])
      supernodes <- supernodes[-supernode_j]
    }
    
    # Remove merged index from ordered list
    ordered_indices <- ordered_indices[-1]
  }
  
  
  return(supernodes)
  
  }
  
  reduce_base_L <- function(L, supernodes){
  
  ## K-space L
   
  ## Initialize reduced Laplacian L'
  N       <- length(supernodes)  # Number of supernodes
  L_prime <- matrix(0, nrow = N, ncol = N)
  
  ## Construct L'
  for (alpha in 1:N) {
    for (beta in 1:N) {
      if (alpha != beta) {
        # L'_alpha_alpha: sum of connections within supernode alpha
        #nodes_alpha <- supernodes[[alpha]]
        #L_prime[alpha, alpha] <- sum(L[c(nodes_alpha), c(nodes_alpha)])
     # } else {
        # L'_alpha_beta: negative of connections between supernodes alpha and beta
        nodes_alpha          <- supernodes[[alpha]]
        nodes_beta           <- supernodes[[beta]]
        L_prime[alpha, beta] <- sum(L[c(nodes_alpha), c(nodes_beta)])/(length(nodes_alpha)*length(nodes_beta))
      }
    }
  }
  
  ## Construct A'
  A_prime <- -L_prime  # A'_alpha_beta = -L'_alpha_beta
  
  ## Remove negative edges 
  A_prime = ifelse(A_prime>0, A_prime, 0)
  
  # Set A'_alpha_alpha = 0 and L'_alpha_alpha = sum_beta A'_alpha_beta
  #diag(A_prime) <- 0
  diag(L_prime) <- colSums(A_prime)
  
  return(list(L=L_prime, Adj=A_prime))
  
  }
  
quantile_normalize <- function(x) {
    ranks <- rank(x, ties.method = "average")
    quantiles <- qnorm((ranks - 0.5) / length(x))
    return(quantiles)
}
  
## Calculate effective resistances
effective_resistance <- function(gg, L) {
  
  ## Compute the Moore-Penrose pseudoinverse of the Laplacian matrix
  L_pseudo_inverse <- MASS::ginv(L)
  
  ## Precompute the edge list
  edge_list <- get.edgelist(gg)
  
  ##Initialize the vector to store effective resistance values
  resistances <- numeric(nrow(edge_list))
  
  ## Calculate the effective resistance for each edge
  for (i in seq_len(nrow(edge_list))) {
    u <- edge_list[i, 1]
    v <- edge_list[i, 2]
    
    ## Effective resistance formula: R_ij = L_pseudo_inverse[i, i] + L_pseudo_inverse[j, j] - 2 * L_pseudo_inverse[i, j]
    resistances[i] <- L_pseudo_inverse[u, u] + L_pseudo_inverse[v, v] - 2 * L_pseudo_inverse[u, v]
  }
  
  return(list(nodes=diag(L_pseudo_inverse), edges=resistances))
}

decimation_A <- function(A,L){

  ## any negative edge weights replace by zero 
  #A=ifelse(A>0,A,0)
  
  ## construct Lnew
  #Lnew       = -A
  #diag(Lnew) = colSums(A)

  n  = nrow(A) 
  
  ## construct new graph
  gg = graph_from_adjacency_matrix(A, mode="undirected", weighted=TRUE)

  ## Get effective resistances for all edges
  resistances = effective_resistance(gg=gg, L=L)   
   
  ## quantile resistance
  q_res_nodes = quantile_normalize(resistances$nodes)
  q_res_edges = quantile_normalize(resistances$edges)
    
  gg = set.vertex.attribute(gg, name="name",          index=V(gg), value=1:n)
  gg = set.vertex.attribute(gg, name="degree",        index=E(gg), value=degree(gg))
  gg = set.vertex.attribute(gg, name="resistances",   index=V(gg), value=resistances$nodes)
  gg = set.vertex.attribute(gg, name="q_resistances", index=V(gg), value=q_res_nodes)
  gg = set.vertex.attribute(gg, name="p_resistances", index=V(gg), value=1-pnorm(q_res_nodes))
  
  gg = set.edge.attribute(gg,   name="resistances",   index=E(gg), value=resistances$edges)
  gg = set.edge.attribute(gg,   name="q_resistances", index=E(gg), value=q_res_edges)
  gg = set.edge.attribute(gg,   name="p_resistances", index=E(gg), value=1-pnorm(q_res_edges))
  
  ## Get the edge list and their weights
  edge_list <- get.edgelist(gg)
  ##weights   <- E(gg)$weight
  
  ## Create a data frame of edges with their resistances
  edges_df <- data.frame(edge_list, 
                         resistances=E(gg)$resistances, 
                         q_resistances = E(gg)$q_resistances,
                         p_resistances = E(gg)$p_resistances,
                         weight=E(gg)$weight)
  
  ## Create a data frame of edges with their resistances
  nodes_df <- data.frame(name=V(gg)$name,
                         degree=V(gg)$degree,
                         resistances=V(gg)$resistances, 
                         q_resistances = V(gg)$q_resistances,
                         p_resistances = V(gg)$p_resistances)
  
  ## Sort edges by effective resistance in ascending order
  #edges_df <- edges_df[order(edges_df$resistances), ]
  
  return(list(gg=gg, edges_df=edges_df, nodes_df=nodes_df))
  
  run=0
  if(run){
  binary_A = ifelse(Adj>0,1,0)
  
  ## Apply majority rule decimation to binary_A
  num_nodes <- nrow(binary_A)
  
  decim_A = matrix(0,num_nodes,num_nodes)
  
  neighbors = binary_A %*% t(binary_A)
  diag(neighbors) = 0
  
  for (i in 1:num_nodes) {
    for (j in 1:num_nodes) {  # Only iterate over upper triangle to avoid duplication
      
        ## Apply majority rule: if more than half of the macronodes have 1, set to 1
        if (neighbors[i,j] > num_nodes / 2 ) {
          decim_A[i, j] <- 1
        } 

      }
    }

  return(decim_A)
  }
  
}


## template preferential attachment graph
## should be similar to sample_pa
template_pa <- function(m=1, gg, prior=NULL){##, directed=FALSE, prior=NULL){
  
  n        = vcount(gg)
  A        = get.adjacency(gg)
  A_names  = as.numeric(colnames(A))
  directed = is.directed(gg)
  
  if( is.null(V(gg)$name) ){
    vnames = 1:n
  } else {
    vnames = V(gg)$name
  }
  
  prior = V(gg)$p_resistances
  names(prior) = vnames
  
  ## Create an empty graph with N nodes
  seed_graph <- graph(edges = NULL, n = n, directed = directed)
  V(seed_graph)$name = vnames
  
  ## set names
  if( is.null(prior) ){ 
    V(seed_graph)$prior = rep(1,n);
  } else { 
      V(seed_graph)$prior = prior[match(V(seed_graph)$name, names(prior), nomatch=NA)];
  }
  
  ## Generate a set of all nodes
  node_set <- 1:vcount(seed_graph)
  
  ## Randomly select two starting nodes from the set
  start_nodes <- sample(node_set, 2)
  
  ## Add the selected starting nodes to the seed graph and remove them from the set
  node_set   <- node_set[-start_nodes]
  
  ## shuffle node_set, so we iterate through nodes in a random order    
  node_set   <- sample(node_set)
  
  ## Loop to add new vertices and edges until reaching N vertices
  for (i in node_set) {
    
    if (ecount(seed_graph) == 0) {
      ## Set uniform probabilities if no edges exist in the seed graph
      probs <- rep(1/vcount(seed_graph), vcount(seed_graph))
    } else {
      ## Calculate probabilities based on degrees if edges exist
      probs <- degree(seed_graph) / sum(degree(seed_graph), na.rm=TRUE)
    }
    
    ## include prior probability
    probs = V(seed_graph)$prior * probs 
    probs = probs / sum(probs, na.rm=TRUE)
    
    ## Select existing vertices
    mask       = A[i,]>0
    node_mask  = A_names[mask]
    probs_mask = probs[mask]
    selected_vertices <- sample(node_mask, m, prob = probs_mask, replace = TRUE)
 
    
    ## Connect the new vertex to selected vertices
    seed_graph <- add_edges(seed_graph, c(rep(i, length(selected_vertices)),
                                          selected_vertices))
    
  }
  
  ## Check for zero-degree nodes and add edges using preferential attachment
  zero_degree_nodes <- which(degree(seed_graph) == 0)
  for (node in zero_degree_nodes) {
    ## Select existing vertices for preferential attachment
    probs <- degree(seed_graph) / sum(degree(seed_graph), na.rm=TRUE)
    
    probs <- V(seed_graph)$prior * probs
    probs <- probs/sum(probs, na.rm=TRUE)
    
    ## Select existing vertices
    mask       = A[node,]>0
    node_mask  = A_names[mask]
    probs_mask = probs[mask]
    selected_vertex <- sample(node_mask, 1, prob = probs_mask, replace = TRUE)
    
    
    ## Connect the zero-degree node to the selected vertex
    seed_graph <- add_edges(seed_graph, c(node, selected_vertex))
  }
  
  return(seed_graph)
  
}



