# estimate rho and eigenfunctions
npeigest <- function(M,G,B.sim,B.mat){
  #
  eigout <- eig(M,G)
  #
  rho <- eigout$rho
  c <- eigout$c
  c.star <- eigout$c.star
  #
  n.sim <- nrow(B.sim)
  #
  if(identical(c,matrix(0,nrow=nrow(eigout$c),1))==TRUE){
    #
    phi <- phi.star <- matrix(1,nrow=n.sim,ncol=1)
    c <- c.star <- solve(G)%*%(t(B.mat[-n,])%*%matrix(1,n-1,1)/n)
    c.norm <- c/as.numeric(sqrt(t(c)%*%G%*%c))
    c.star.norm <- c.star/as.numeric(t(c.star)%*%G%*%c.norm)
    #
  }else{
    #
    c.norm <- c/as.numeric(sqrt(t(c)%*%G%*%c))
    phi <- abs(B.sim%*%c.norm)
    #
    c.star.norm <- c.star/as.numeric(t(c.star)%*%G%*%c.norm)
    phi.star <- abs(B.sim%*%c.star.norm)
    #
  }
  #
  return(list(rho=rho,phi=phi,phi.star=phi.star,c.hat=c.norm,c.star.hat=c.star.norm))
}
#
# generalized eigenvector computation wrapper function
eig <- function(M,G){
  #
  qzout <- qz.dggev(A=M,B=G,vl=TRUE,vr=TRUE)
  #
  # check stability
  if(qzout$INFO!=0){
    warning("eigenvector problem numerically unstable")
  }
  #
  eig.real <- qzout$ALPHAR/qzout$BETA
  #
  # check that rho is from real eigenvalue
  if(max(abs(qzout$ALPHAI[(eig.real==max(eig.real))]))!=0){
    #
    warning("largest eigenvalue has complex components")
    #
    # set to default
    rho <- 1
    c <- c.star <- matrix(0,nrow=ncol(G),ncol=1)
    #
  }else{
    if(sum((eig.real==max(eig.real))*1)>1){
      #
      warning("largest eigenvalue not uniquely defined")
      #
      # set to default
      rho <- 1
      c <- c.star <- matrix(0,nrow=ncol(G),ncol=1)
      #
    }else{
      #
      rho <- max(eig.real)
      #
      c <- matrix(qzout$VR[,eig.real==max(eig.real)])
      c.star <- matrix(qzout$VL[,eig.real==max(eig.real)])
      #
    }
  }
  #
  # if rho > 2 then set to defaults
  if(rho > 2){
    #
    warning("largest eigenvalue larger than 2, setting to default")
    #
    rho <- 1
    c <- c.star <- matrix(0,nrow=ncol(G),ncol=1)
    #
  }
  #
  return(list(rho=rho,c=c,c.star=c.star))
}
#
# tensor product basis
tp <- function(x1,x2){
  #
  n <- nrow(x1)
  k1<- ncol(x1)
  k2<- ncol(x2)
  uvec <- matrix(1,1,k1)
  X <- matrix(0,n,k1*k2)
  for(i in 1:k2){
   X[,seq(k1*(i-1)+1,i*k1)] <- x1*(x2[,i]%*%uvec)
  }
  return(X)
}
#
# sparse product basis
sp <- function(x1,x2){
  #
  n <- nrow(x1)
  k1<- ncol(x1)
  k2<- ncol(x2)
  #
  if(k1!=k2){stop("dimensions of univariate bases aren't equal")}
  #
  X <- matrix(0,n,k1*(k1+1)/2)
  ktemp <- 1
  for(i in 1:k2){
    for(j in 1:k1){
      if(i+j<=k1+1){
        X[,ktemp] <- x1[,j]*x2[,i]
        ktemp <- ktemp+1
      }
    }
  }
  return(X)
}
#
# sparse kronecker product
sp.kronecker <- function(x1,x2){
  #
  n1 <- nrow(x1)
  n2 <- nrow(x2)
  k1<- ncol(x1)
  k2<- ncol(x2)
  #
  if(k1!=k2){stop("dimensions of univariate bases aren't equal")}
  #
  uvec <- matrix(1,1,k1)
  X <- matrix(0,n1*n2,k1*(k1+1)/2)
  ktemp <- 1
  for(i in 1:k2){
    for(j in 1:k1){
      if(i+j<=k1+1){
        X[,ktemp] <- kronecker(x1[,j],x2[,i])
        ktemp <- ktemp+1
      }
    }
  }
  return(X)
}
#
# row-wise quantiles
rowQuantile <- function(X,q){
  #
  n <- nrow(X)
  qvec <- matrix(NA,n,1)
  for(i in 1:n){
    qvec[i] <- quantile(X[i,],q,names=FALSE)
  }
  return(qvec)
}
#
# compute RMSE of n by N matrix array of functions
rmse_mat <- function(fn.est,fn.true,weight=matrix(1,length(fn.true),1)){
  #
  n <- length(fn.true)
  N <- ncol(fn.est)
  mse.vec <- numeric(N)
  w.norm <- weight/sum(weight)
  for(i in 1:N){
    mse.vec[i] <- sum(w.norm*(fn.est[,i]-fn.true)^2)
  }
  return(sqrt(mean(mse.vec)))
}
#
# R version of find() from matlab
rfind <- function(x){
  seq(along=x)[as.logical(x)]
}
#
# construct PC and TC series form eigenvalue estimates
pctc.series <- function(m.vec,B.mat,rho,c.hat){
  #
  n <- length(m.vec)
  #
  phi.ts <- abs(B.mat%*%c.hat)
  phi0.ts <- phi.ts[-(n+1)]
  phi1.ts <- phi.ts[-1]
  #
  pc.ts <- (1/rho)*m.vec*phi1.ts/phi0.ts
  tc.ts <- rho*phi0.ts/phi1.ts
  #
  return(list(pc=pc.ts,tc=tc.ts))
}
#
# Gaussian VAR(1)
var1.est <- function(x1,x2){
  #
  fm <- VAR(cbind(x1,x2),p=1,type="const")
  #
  c1 <- matrix(fm$varresult$x1$coefficients)
  c2 <- matrix(fm$varresult$x2$coefficients)
  #
  A <- rbind(c1[-3],c2[-3])
  #
  mu <- rbind(c1[3],c2[3])
  #
  Sigma <- cov(cbind(fm$varresult$x1$residuals,fm$varresult$x2$residuals))
  #
  return(list(A=A,mu=mu,Sigma=Sigma))
}
