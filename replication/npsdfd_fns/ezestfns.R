# estimate (beta,gamma) for robust preference model from state variables and returns
est.ez.fn <- function(beta0,gamma0,G,g,B,R){
  #
  b0 <- log(beta0/(1-beta0))
  g0 <- log(gamma0-1)
  theta <- c(b0,g0)
  #
  G1 <- solve(G)
  #
  opt <- optim(theta,est.ez.obj,gr=NULL,G=G,g=g,B=B,G1=G1,R=R,method="Nelder-Mead",control=list(maxit=10000,reltol=1e-8),hessian=FALSE)
  #
  beta <- exp(opt$par[1])/(1+exp(opt$par[1]))
  gamma <- 1+exp(opt$par[2])
  #
  return(list(beta=beta,gamma=gamma))
}
#
# objective function for esitmation of (beta,gamma)
est.ez.obj <- function(theta,G,g,B,G1,R){
  #
  beta <- exp(theta[1])/(1+exp(theta[1]))
  gamma <- exp(theta[2])+1
  #
  # extract fixed point
  fpout <- nlfp(g,B,G,beta,1/(gamma-1),niter=250)
  lambda <- fpout$lambda
  efn <- fpout$efn
  #
  # form SDF
  mvec <- (beta/lambda)*exp(-gamma*g[-1])*(pmax(efn[-1],0)^beta)/pmax(efn[-n],0)
  #
  # pricing errors
  err <- (matrix(mvec)%*%matrix(1,nrow=1,ncol=ncol(R)))*R-1
  #
  # conditional GMM objective function
  n <- length(g)
  B0 <- t(B[-n,])
  B.err <- t(B0%*%err)/n
  Q <- n*sum(diag(B.err%*%G1%*%t(B.err)))
  #
  return(Q)
}
# as above, but more flexible specification allowing for more instruments
est.ez.fn.flex <- function(beta0,gamma0,G,g,B,B.inst,R){
  #
  b0 <- log(beta0/(1-beta0))
  g0 <- log(gamma0-1)
  theta <- c(b0,g0)
  #
  n <- length(g)
  G.inst <- t(B.inst[1:n-1,])%*%B.inst[1:n-1,]/n
  G1.inst <- solve(G.inst)
  #
  opt <- optim(theta,est.ez.flex.obj,gr=NULL,G=G,g=g,B=B,B.inst=B.inst,G1.inst=G1.inst,R=R,method="Nelder-Mead",control=list(maxit=10000,reltol=1e-8),hessian=FALSE)
  #
  beta <- exp(opt$par[1])/(1+exp(opt$par[1]))
  gamma <- 1+exp(opt$par[2])
  #
  return(list(beta=beta,gamma=gamma))
}
#
# objective function for esitmation of (beta,gamma)
est.ez.flex.obj <- function(theta,G,g,B,B.inst,G1.inst,R){
  #
  beta <- exp(theta[1])/(1+exp(theta[1]))
  gamma <- exp(theta[2])+1
  #
  # extract fixed point
  fpout <- nlfp(g,B,G,beta,1/(gamma-1),niter=250)
  lambda <- fpout$lambda
  efn <- fpout$efn
  #
  # form SDF
  mvec <- (beta/lambda)*exp(-gamma*g[-1])*(pmax(efn[-1],0)^beta)/pmax(efn[-n],0)
  #
  # pricing errors
  err <- (matrix(mvec)%*%matrix(1,nrow=1,ncol=ncol(R)))*R-1
  #
  # conditional GMM objective function
  n <- length(g)
  B0 <- t(B.inst[-n,])
  B.err <- t(B0%*%err)/n
  Q <- n*sum(diag(B.err%*%G1.inst%*%t(B.err)))
  #
  return(Q)
}
