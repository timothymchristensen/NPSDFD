#
# population quantities for CCAPM design in univariate case
cc_linear_gaussian <- function(beta,gamma,pars,g.grid){
  #
  kappa <- pars$kappa
  mu <- pars$mu
  sigma <- pars$sigma
  #
  rho <- beta*exp( -gamma*mu + 0.5*((gamma*sigma)/(1-kappa))^2  )
  y <- -log(rho)
  L <- log(rho) - log(beta) + gamma*mu
  #
  a1.star <- -gamma/(1-kappa)
  a1 <- a1.star*kappa
  #
  a0 <- mu*gamma*kappa/(1-kappa) - (gamma*kappa*sigma)^2/((1-kappa^2)*(1-kappa)^2)
  a0.star <- mu*gamma/(1-kappa) + ((gamma*sigma)^2/((1-kappa^2)*(1-kappa)^2))*(kappa^2-0.5*(1+kappa)^2)
  #
  phi.sim <- exp(a0 + a1*g.grid)
  phi.star.sim <- exp(a0.star + a1.star*g.grid)
  #
  return(list(rho=rho,y=y,L=L,phi.sim=phi.sim,phi.star.sim=phi.star.sim))
}
#
# population quantities for recursive/robust design in univariate case
ez_linear_gaussian <- function(beta,gamma,pars,g.grid){
  #
  kappa <- pars$kappa
  mu <- pars$mu
  sigma <- pars$sigma
  #
  theta <- 1/(gamma-1)
  #
  b1 <- -kappa/(theta*(1-kappa*beta))
  b2 <- (1-gamma+beta*b1)
  b0 <- (b2*(1-kappa)*mu + 0.5*b2^2*sigma^2)/(1-beta)
  c0 <- -(b1*mu + b1^2*sigma^2/(1-kappa^2))
  #
  v.sim <- exp(b0 + b1*g.grid)
  efn.sim <- exp(c0 + b1*g.grid)
  lambda <- exp((b0-c0)*(1-beta))
  #
  a1 <- (-b1+kappa*(beta*b1-gamma))/(1-kappa)
  a1.star <- ((beta-kappa)*b1-gamma)/(1-kappa)
  #
  a0 <- -a1*mu - a1^2*sigma^2/(1-kappa^2)
  a0.star <- - (a0 + (a1+a1.star)*mu + 0.5*(a1+a1.star)^2*sigma^2/(1-kappa^2) )
  #
  phi.sim <- exp(a0 + a1*g.grid)
  phi.star.sim <- exp(a0.star + a1.star*g.grid)
  #
  nu <- beta*b1-gamma+a1
  rho <- (beta/lambda)*exp( (beta-1)*c0+nu*(1-kappa)*mu + 0.5*nu^2*sigma^2 )
  y <- -log(rho)
  L <- log(rho)-log(beta/lambda)-(beta-1)*c0-(beta*b1-gamma-b1)*mu
  #
  V <- sigma^2/(1-kappa^2)
  V1 <- rbind(c(V,kappa*V),c(kappa*V,kappa^2*V+sigma^2))
  v.p <- rbind(-(a1+b1),a1+beta*b1-gamma)
  v.t <- rbind(a1,-a1)
  cov.pt <- t(v.p)%*%V1%*%v.t
  var.p <- t(v.p)%*%V1%*%v.p
  var.t <- t(v.t)%*%V1%*%v.t
  cor.pt <- cov.pt/sqrt(var.p*var.t)
  #
  return(list(rho=rho,y=y,L=L,lambda=lambda,phi.sim=phi.sim,phi.star.sim=phi.star.sim,v.sim=v.sim,efn.sim=efn.sim,cor.pt=cor.pt,cov.pt=cov.pt))
}
#
# population quantities for recursive/robust design in multivariate case
ez_linear_gaussian_mv <- function(beta,gamma,pars){
  #
  A <- pars$A
  mu <- pars$mu
  Sigma <- pars$Sigma
  #
  theta <- 1/(gamma-1)
  #
  dimx <- length(mu)
  #
  Gamma <- matrix(0,nrow=dimx,ncol=1)
  Gamma[1] <- 1
  #
  b1 <- -solve(diag(dimx)-beta*t(A))%*%t(A)%*%Gamma/theta
  h0 <- beta*b1-Gamma/theta
  b0 <- (t(h0)%*%(diag(dimx)-A)%*%mu+0.5*t(h0)%*%Sigma%*%h0)/(1-beta)
  #
  V <- matrix(0,nrow=dimx,ncol=dimx)
  for(i in 1:100){
    V <- A%*%V%*%t(A) + Sigma
  }
  #
  c0 <- -t(b1)%*%mu - t(b1)%*%V%*%b1
  lambda <- exp((b0-c0)*(1-beta))
  #
  a1 <- solve(diag(dimx)-t(A))%*%(t(A)%*%(beta*b1-gamma*Gamma)-b1)
  # 
  g0 <- a1 + beta*b1 - gamma*Gamma
  #
  rho <- (beta/lambda)*exp( (beta-1)*c0 + t(g0)%*%(diag(dimx)-A)%*%mu + 0.5%*%t(g0)%*%Sigma%*%g0 )
  y <- -log(rho)
  L <- log(rho)-log(beta/lambda)-(beta-1)*c0-t(beta*b1-gamma*Gamma-b1)%*%mu
  #
  V1 <- rbind(cbind(V,V%*%t(A)),cbind(A%*%V,V))
  v.p <- rbind(-(b1+a1),g0)
  v.t <- rbind(a1,-a1)
  cov.pt <- t(v.p)%*%V1%*%v.t
  var.p <- t(v.p)%*%V1%*%v.p
  var.t <- t(v.t)%*%V1%*%v.t
  cor.pt <- cov.pt/sqrt(var.p*var.t)
  #
  return(list(rho=rho,y=y,L=L,lambda=lambda,cor.pt=cor.pt,cov.pt=cov.pt))
}
#
# population quantities for recursive/robust design with stoch. volatility
ez_linear_gaussian_sv <- function(beta,gamma,pars){
  #
  # unpack parametrs
  mu <- pars$mu
  kappa <- pars$kappa
  phi <- pars$phi
  c <- pars$c
  delta <- pars$delta
  #
  # solve analytically for fixed point
  a1 <- (1-gamma)*kappa/(1-beta*kappa)
  qa <- beta*c
  nu <- 0.5*(beta*a1+1-gamma)^2
  qb <- phi*beta-1-beta*c*nu
  qc <- nu
  a2 <- (-qb - sqrt(qb^2-4*qa*qc))/(2*qa)
  a0 <- (( (1-gamma+beta*a1)*(1-kappa)*mu-delta*log(1-beta*a2*c) )/(1-beta))
  #
  # convert to SDF parameters
  d0 <- beta*exp(-(1-beta)*a0)
  d1 <- beta*a1-gamma
  d2 <- beta*a2
  d3 <- -a1
  d4 <- -a2
  #
  # solve for eigenfunction/value parameters
  c1 <- (d1*kappa+d3)/(1-kappa)
  xi <- 0.5*(d1+c1)^2+d4+d2
  qa <- c
  qb <- phi-1-c*xi
  qc <- xi
  c3 <- (-qb - sqrt(qb^2-4*qa*qc))/(2*qa)
  c2 <- c3-d2
  rho <- d0*exp( (d1+c1)*(1-kappa)*mu - delta*log(1-c*(d2+c2)) )
  #
  L <- log(rho) - (log(d0) + (d1+d3)*mu + (d2+d4)*delta*c/(1-phi))
  #
  y <- -log(rho)
  #
  v.p <- matrix(c(d3-c1,d4-c2,d1+c1,d2+c2))
  v.t <- matrix(c(c1,c2,-c1,-c2))
  #
  V <- matrix(0,4,4)
  vg <- c*delta/((1-phi)*(1-kappa^2))
  vv <- c^2*delta/((1-phi)^2)
  V[1,1] <- V[3,3] <- vg
  V[2,2] <- V[4,4] <- vv
  V[1,3] <- V[3,1] <- kappa*vg
  V[2,4] <- V[4,2] <- phi*vv
  #
  cov.pt <- t(v.p)%*%V%*%v.t
  var.p <- t(v.p)%*%V%*%v.p
  var.t <- t(v.t)%*%V%*%v.t
  cor.pt <- cov.pt/sqrt(var.p*var.t)
  #
  return(list(rho=rho,y=y,cor.pt=cor.pt,cov.pt=cov.pt,L=L))
}
#
# wrapper function for recursive/robust design with stoch. volatility
ez_sv_emp <- function(beta,gamma_grid,pars){
  #
  pars <- list(mu = 0.00549,kappa = 0.276, phi = 0.5, c = 7e-5, delta = 0.18)
  #
  # population quantities
  y.lgsv.ez <- L.lgsv.ez <- cov.lgsv.ez <- cor.lgsv.ez <- numeric(length(gamma_grid))
  #
  for(i in 1:length(gamma_grid)){
    #
    ez_lgsv_out <- ez_linear_gaussian_sv(beta,gamma_grid[i],pars)
    #
    y.lgsv.ez[i] <- ez_lgsv_out$y
    cov.lgsv.ez[i] <- ez_lgsv_out$cov.pt
    cor.lgsv.ez[i] <- ez_lgsv_out$cor.pt
    L.lgsv.ez[i] <- ez_lgsv_out$L
  }
  #
  return(list(y.lgsv.ez=y.lgsv.ez,cov.lgsv.ez=cov.lgsv.ez,cor.lgsv.ez=cor.lgsv.ez,L.lgsv.ez=L.lgsv.ez,pars=pars))
}
#
# estimate VAR(1) by OLS
var1.est <- function(X){
  #
  if(nrow(X)>1){
    X1 <- X[,-1]
    X0 <- X[,-ncol(X)]
  }else{
    X1 <- t(X[-1])
    X0 <- t(X[-length(X)])
  }
  n <- ncol(X1)
  Y1 <- X1
  Y0 <- rbind(matrix(1,1,n),X0)
  #
  theta_hat <- Y1%*%t(Y0)%*%solve(Y0%*%t(Y0))
  #
  A <- theta_hat[,-1]
  mu <- solve(diag(nrow(Y1))-A)%*%theta_hat[,1]
  #
  resid <- Y1 - theta_hat%*%Y0
  #
  Sigma <- cov(t(resid))
  #
  return(list(A=A,mu=mu,Sigma=Sigma,residuals=resid))
}
#
# estimate parameters of AR(1) with ARG(1) stoch. volatility via indirect inference
ar1arg1.ii.est <- function(x){
  #
  # specify and fit the GARCH(1,1) auxiliary model
  ar1.pars <- var1.est(t(x))
  ar1.residuals <- as.numeric(ar1.pars$residuals)
  ar1.pars <- c(ar1.pars$mu,ar1.pars$A)
  spec <- ugarchspec(variance.model = list(model = "sGARCH", 
                                           garchOrder = c(1, 1), 
                                           submodel = NULL, 
                                           external.regressors = NULL, 
                                           variance.targeting = FALSE), 
                     mean.model = list(armaOrder = c(0,0), include.mean = FALSE),
                     distribution.model = "norm")
  garch <- ugarchfit(spec = spec, data = ar1.residuals, solver="hybrid", solver.control = list(trace=0))
  beta0 <- as.numeric(c(ar1.pars,garch@fit$coef))
  #
  # initialize
  pars_0 <- list(mu=ar1.pars[1],kappa=ar1.pars[2],phi = 0.5, c = 7e-5, delta = 0.18)
  theta <- c(pars_0$mu,log(pars_0$kappa/(1-pars_0$kappa)),log(pars_0$phi/(1-pars_0$phi)),log(pars_0$c),log(pars_0$delta))
  #
  # optimize
  opt <- optim(theta,ar1arg1.ii.obj,gr=NULL,spec=spec,beta0=beta0,n=length(x),method="Nelder-Mead",control=list(maxit=10000,reltol=1e-12),hessian=FALSE)
  #
  theta <- opt$par
  pars <- list()
  pars$mu <- theta[1]
  pars$kappa <- exp(theta[2])/(1+exp(theta[2]))
  pars$phi <- exp(theta[3])/(1+exp(theta[3]))
  pars$c <- exp(theta[4])
  pars$delta <- exp(theta[5])
  #
  return(pars)
}
#
# indirect inference objective function
ar1arg1.ii.obj <- function(theta,beta0,spec,n){
  #
  # fix the seed
  set.seed(1988)
  #
  # reformat argument
  pars <- list()
  pars$mu <- theta[1]
  pars$kappa <- exp(theta[2])/(1+exp(theta[2]))
  pars$phi <- exp(theta[3])/(1+exp(theta[3]))
  pars$c <- exp(theta[4])
  pars$delta <- exp(theta[5])
  #
  H <- 500
  BB <- matrix(NA,nrow=H,ncol=5)
  #
  # initialize objective function
  obj <- 0
  #
  for(h in 1:H){
    # simulate data
    sd <- ar1arg1.sim(pars,n)
    #
    # fit auxiliary model
    ar1.pars <- var1.est(t(sd$g))
    ar1.residuals <- as.numeric(ar1.pars$residuals)
    ar1.pars <- c(ar1.pars$mu,ar1.pars$A)
    garch <- ugarchfit(spec = spec, data = ar1.residuals, solver="hybrid", solver.control = list(trace=0))
    #
    if(garch@fit$convergence!=0){
      beta <- numeric(5)
    }else{
      beta <- as.numeric(c(ar1.pars,garch@fit$coef))
    }
    #
    BB[h,] <- beta
  }
  #
  obj <- n*sum((colMeans(BB)-beta0)^2)
  #
  cat("Current value of objective function: ",obj,"\n")
  #
  return(obj)
}
#
# Simulate AR(1) with ARG(1) stoch. vol
ar1arg1.sim <- function(pars,n){
  #
  # burn in
  b <- 1000
  #
  # unpack parameters
  mu <- pars$mu
  kappa <- pars$kappa
  phi <- pars$phi
  c <- pars$c
  delta <- pars$delta
  #
  arg1pars <- list(phi=phi,c=c,delta=delta)
  #
  # simulate ARG(1) processes
  v <- arg1.sim(arg1pars,n+b)
  #
  # generate Gaussian innovations
  u <- rnorm(n+b)*sqrt(v)
  #
  x <- numeric(n+b)
  x[1] <- u[1]
  for(i in 2:(n+b)){
    x[i] <- kappa*x[i-1] + u[i] 
  }
  #
  g <- x[-(1:b)]+mu
  v <- v[-(1:b)]
  #
  return(list(g=g,v=v))
}
#
# Simulate ARG(1) process
arg1.sim <- function(pars,n){
  #
  phi <- pars$phi
  c <- pars$c
  delta <- pars$delta
  #
  x <- numeric(n)
  #
  # initial draw
  x[1] <- rgamma(1,shape=delta,scale=c/(1-phi))
  #
  for(i in 2:n){
    z <- rpois(1,phi*x[i-1]/c)
    x[i] <- rgamma(1,shape=(delta+z),scale=c)
  }
  #
  return(x)
}
#


