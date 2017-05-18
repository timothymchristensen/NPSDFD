# Recursion to compute fixed point via nonlinear Perron-Frobenius
nlfp <- function(g,B,G,beta,theta,niter=100,B.sim=0){
  #
  n <- length(g)
  #
  G1 <- solve(G)
  B0 <- t(B[-n,])
  B1 <- B[-1,]
  #
  # initialize at projection of constant function
  c0     <- G1%*%(B0%*%matrix(1,n-1,1)/n)
  c0norm <- c0/as.numeric(sqrt(t(c0)%*%G%*%c0))
  #
  i <- 1
  while(i <= niter){
    #
    tvec <- matrix(exp(-g[-1]/theta)*(abs(B1%*%c0norm))^beta)
    Tmat <- B0%*%tvec/n
    #
    c1 <- G1%*%Tmat
    #
    c0norm <- c1/as.numeric(sqrt(t(c1)%*%G%*%c1))
    i <- i + 1
  }
  #
  if(min(B1%*%c0norm)<0){
    warning("value function truncated at zero")
  }
  efn <- abs(B%*%c0norm)
  #
  # rescaling via eigenvalue
  lambda <- as.numeric(sqrt(t(c1)%*%G%*%c1))
  scale <- lambda^(1/(1-beta))
  vfn <- efn*scale
  #
  if(identical(B.sim,0)==TRUE){
    return(list(lambda=lambda,vfn=vfn,efn=efn))
  }else{
    v.sim <- B.sim%*%c0norm*scale
    efn.sim <- B.sim%*%c0norm
    return(list(lambda=lambda,vfn=vfn,efn=efn,v.sim=v.sim,efn.sim=efn.sim))
  }
}
#
