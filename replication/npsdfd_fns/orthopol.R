# generalized hermite basis
genhermite <- function(x,k,mu=0,sigma=1){
  x <- (x-mu)/sigma
  n <- length(x)
  X <- matrix(0,n,k)
  X[,1] <- 1
  if(k > 1){
    X[,2] <- x
  }
  if(k > 2){
    for(m in 3:k){
      X[,m] <- x*X[,m-1]-(m-2)*X[,m-2]
    }
  }
  return(X)
}
#
# B-spline basis
bspline <- function(x,r=4,knots=NULL){
  n <- length(x)
  if(is.null(knots)){
    warning("no knots specified, using default")
    knots <- c(min(x),max(x))
  }
  m <- length(knots)
  if(knots[1]>min(x)){
    knots[1] <- min(x)
    warning("decreasing minimum point")
  }
  if(knots[m]<max(x)){
    knots[m] <- max(x)
    warning("increasing maximum point")
  }
  kts <- matrix(0,1,m+2*r-2)
  kts[seq(1,r-1)] <- knots[1]
  kts[seq(m+r,m+2*r-2)] <- knots[m]
  kts[seq(r,m+r-1)] <- knots
  BB <- array(0,c(n,m+2*r-2,r))
  #
  # initialization with Haar function
  for(i in 1:n){
    ix0 <- rfind(x[i] >= kts[seq(r,r+m)] & x[i] <= kts[seq(r+1,r+m+1)])
    ix1 <- ix0[1]
    BB[i,ix1+r-1,1] <- 1
  }
  #
  # recursion
  for(j in 2:r){
    for(i in 1:(m+2*r-2-j)){
      if(i+j<=m+2*r-2){
        if(kts[i+j-1]-kts[i]!=0){
          a1 <- (x-kts[i])/(kts[i+j-1]-kts[i])
        }else{
          a1 <- matrix(0,n,1)
        }
        if(kts[i+j]-kts[i+1]!=0){
          a2 <- (x-kts[i+1])/(kts[i+j]-kts[i+1])
        }else{
          a2 <- matrix(0,n,1)
        }
      BB[,i,j] <- a1*BB[,i,j-1]+(1-a2)*BB[,i+1,j-1]
      }
    }
  }
  #
  return(BB[,seq(1:(m+r-2)),r])
}



