#
# sdf series and M matrix for power utility
mmat_ccapm <- function(g,B.mat,beta,gamma){
  #
  n.g <- length(g)
  n.B <- nrow(B.mat)
  if(n.g!=n.B){
    stop("g and B.mat are not of compatible length")
  }
  #
  n <- n.g
  #
  mvec <- beta*exp(-gamma*g[-1])
  #
  B.0 <- B.mat[-n,]
  B.1 <- B.mat[-1,]
  M.1 <- B.1*(mvec%*%matrix(1,1,ncol(B.mat)))
  Mhat <- t(B.0)%*%M.1/n
  #
  return(list(mvec=mvec,Mhat=Mhat))
}
#
# sdf series and M matrix for recursive preferences with EIS=1
mmat_ez <- function(g,B.mat,G,beta,gamma,niter=100,B.sim=0){
  #
  n.g <- length(g)
  n.B <- nrow(B.mat)
  if(n.g!=n.B){
    stop("g and B.mat are not of compatible length")
  }
  #
  n <- n.g
  #
  theta <- 1/(gamma-1)
  #
  if(identical(B.sim,0)==TRUE){
    fpout <- nlfp(g,B.mat,G,beta,theta,niter)
  }else{
    fpout <- nlfp(g,B.mat,G,beta,theta,niter,B.sim)
  }
  lambda <- fpout$lambda
  efn <- fpout$efn
  #
  mvec <- (beta/lambda)*exp(-gamma*g[-1])*(pmax(efn[-1],0)^beta)/pmax(efn[-n],0)
  #
  B.0 <- B.mat[-n,]
  B.1 <- B.mat[-1,]
  M.1 <- B.1*(mvec%*%matrix(1,1,ncol(B.mat)))
  Mhat <- t(B.0)%*%M.1/n
  #
  if(identical(B.sim,0)==TRUE){
    return(list(mvec=mvec,Mhat=Mhat,lambda=fpout$lambda))
  }else{
    return(list(mvec=mvec,Mhat=Mhat,lambda=fpout$lambda,efn.sim=fpout$efn.sim))
  }
}
#
# wraper function for recursive preferences
ez.decomp <- function(g,B.mat,G,beta,gamma,niter=250,B.sim){
  #
  # solve eigenvector problem
  mmat_out <- mmat_ez(g,B.mat,G,beta,gamma,niter,B.sim)
  Mhat.ez <- mmat_out$Mhat
  est.ez <- npeigest(Mhat.ez,G,B.sim,B.mat)
  #
  # eigenfunctions
  phi <- est.ez$phi
  phi.star <- est.ez$phi.star
  c <- est.ez$c.hat
  c.star <- est.ez$c.star.hat
  #
  # functionals
  rho <- est.ez$rho
  y <- -log(est.ez$rho)
  L <- log(est.ez$rho) - mean(log(mmat_out$mvec))
  #
  lambda <- mmat_out$lambda
  chi <- mmat_out$efn.sim
  #
  return(list(rho=rho,y=y,L=L,lambda=lambda,phi=phi,phi.star=phi.star,chi=chi,c=c,c.star=c.star,mvec=mmat_out$mvec))
}

