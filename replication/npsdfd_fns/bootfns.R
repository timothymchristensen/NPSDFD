#
# Stationary bootstrap indices
bsix <- function(n,b){
  #
  ix <- numeric()
  #
  while(length(ix) < n){
    #
    ixa <- floor(n*runif(1))
    inc <- rgeom(1,1/b)
    # inc <- b
    ixb <- (ixa:(ixa+inc))%%n
    #
    ix <- c(ix,ixb)
  }
  #
  ix <- ix[(1:n)]+1
}
#
# Bootstrap the estimators and report bootstrap CIs
# set est_pref_pars=TRUE to estimate preference parameters, otherwise does resampling for fixed (beta0,gamma0)
ez.bs <- function(g,B.mat,nboot=100,blocklength=1,beta0,gamma0,est_pref_pars=FALSE,B.inst=NULL,R=NULL){
  #
  # set the seed so reproducible 
  set.seed(1988)
  #
  rho_boot <- y_boot <- L_boot <- beta_boot <- gamma_boot <- lambda_boot <- numeric(nboot)
  #
  n <- length(g)
  #
  for(i in 1:nboot){
    #
    if(i%%50==0){
      cat("Bootstrap iteration:",i,"out of",nboot,"\n")
    }
    #
    # resample data
    ix <- bsix(n,blocklength)
    g_boot <- g[ix]
    B.mat_boot <- B.mat[ix,]
    Ghat_boot <- t(B.mat_boot[1:n-1,])%*%B.mat_boot[1:n-1,]/n
    #
    # estimate preference parameters if required
    if(est_pref_pars==TRUE){
      #
      R_boot <- R[ix,]
      B.inst_boot <- B.inst[ix,]
      #
      theta_hat_boot <- est.ez.fn.flex(beta0,gamma0,Ghat_boot,g_boot,B.mat_boot,B.inst_boot,R_boot[-1,])
      beta_hat_boot <- theta_hat_boot$beta
      gamma_hat_boot <- theta_hat_boot$gamma
      #
    }else{
      #
      beta_hat_boot <- beta0
      gamma_hat_boot <- gamma0
      #
    }
    #
    # decomposition of estimated SDF
    mmat.ez.out_boot <- mmat_ez(g_boot,B.mat_boot,Ghat_boot,beta_hat_boot,gamma_hat_boot,niter=250,B.sim)
    Mhat.ez_boot <- mmat.ez.out_boot$Mhat
    est.ez_boot <- npeigest(Mhat.ez_boot,Ghat_boot,B.mat_boot,B.mat_boot)
    #
    rho_boot[i] <- est.ez_boot$rho
    y_boot[i] <- -log(est.ez_boot$rho)
    L_boot[i] <- log(est.ez_boot$rho)-mean(log(mmat.ez.out_boot$mvec))
    #
    beta_boot[i] <- beta_hat_boot
    gamma_boot[i] <- gamma_hat_boot
    lambda_boot[i] <- mmat.ez.out_boot$lambda
  }
  #
  return(list(rho_boot=rho_boot,y_boot=y_boot,L_boot=L_boot,beta_boot=beta_boot,gamma_boot=gamma_boot,lambda_boot=lambda_boot))
}
#
# equal-tailed bootstrap confidence interval
bsci <- function(par,par_bs,alpha){
  #
  cl <- 2*par-as.numeric(quantile(par_bs,1-alpha/2))
  cu <- 2*par-as.numeric(quantile(par_bs,alpha/2))
  #
  return(list(cl=cl,cu=cu))
}
#
