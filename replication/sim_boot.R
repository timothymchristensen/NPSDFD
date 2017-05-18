# Run this file to perform simulations 
#
rm(list=ls())
#
# load packages and codes
library("MASS")
library("QZ")
library("ggplot2")
source("./npsdfd_fns/orthopol.R")
source("./npsdfd_fns/matfns.R")
source("./npsdfd_fns/Mest.R")
source("./npsdfd_fns/lgsoln.R")
source("./npsdfd_fns/nlfp.R")
#
#   Parameters
mu    <- 0.005  # AR(1) process mean
kappa <- 0.6    # AR(1) ar parameter
sigma <- 0.01   # conditional std deviation
beta  <- 0.994
gamma <- 15
#
# Configuration
N     <- 5e4    # number of simulations 
k     <- 8      # smoothing parameter for phi
n.vec <- c(400,800,1600,3200) # sample sizes for simulations
#
# Basis
basis <- 1      # 1 for Hermite polynomials, 2 for B-splines
#
# grid for plotting
g.grid <- seq(-3*sigma+mu,3*sigma+mu,length.out=2e2)
#
# Preallocate matrices for results
r.cc <- y.cc <- L.cc <- r.ez <- y.ez <- L.ez <- lambda.ez <- matrix(NA,nrow=N,ncol=length(n.vec))
phi.cc <- phi.star.cc <- phi.ez <- phi.star.ez <- efn.ez <- array(NA,c(length(g.grid),N,length(n.vec)))
#
# simulate data and store estimates
for(j in 1:length(n.vec)){
  #
  # initialize
  set.seed(1988)
  #
  n <- n.vec[j]
  #
  for(i in 1:N){
    #
    # simulate data
    g <- arima.sim(model=list(ar =kappa,order=c(1,0,0)),n=n,sd=sigma,n.start=200)+mu
    # 
    # generate bases for eigenfunction estimation and G matrix
    if(basis==1){
      B.mat <- genhermite(g,k,mu=mu,sigma=sigma/sqrt(1-kappa^2))
      B.sim <- genhermite(g.grid,k,mu=mu,sigma=sigma/sqrt(1-kappa^2))
    }
    if(basis==2){
      kts <- qnorm(seq(0.001,0.999,length.out=k-2),mean=mu,sd=sigma/sqrt(1-kappa^2))
      bspline.out <- bspline(rbind(matrix(g),matrix(g.grid)),r=4,knots=kts)
      B.mat <- bspline.out[seq(1,n),]
      B.sim <- bspline.out[seq((n+1),nrow(bspline.out)),]
    }
    if(basis!=1&basis!=2){
      stop("wrong basis flag, use 1 or 2")
    }
    Ghat <- t(B.mat[1:n-1,])%*%B.mat[1:n-1,]/n
    #
    # generate M matrices 
    mmat.cc.out <- mmat_ccapm(g,B.mat,beta,gamma)
    Mhat.cc <- mmat.cc.out$Mhat
    mmat.ez.out <- mmat_ez(g,B.mat,Ghat,beta,gamma,niter=100,B.sim=B.sim)
    Mhat.ez <- mmat.ez.out$Mhat
    #
    # compute estimators
    est.cc <- npeigest(Mhat.cc,Ghat,B.sim,B.mat)
    est.ez <- npeigest(Mhat.ez,Ghat,B.sim,B.mat)
    #
    # store estimates of parameters
    r.cc[i,j] <- est.cc$rho
    y.cc[i,j] <- -log(est.cc$rho)
    L.cc[i,j] <- log(est.cc$rho)-mean(log(mmat.cc.out$mvec))
    #
    r.ez[i,j] <- est.ez$rho
    y.ez[i,j] <- -log(est.ez$rho)
    L.ez[i,j] <- log(est.ez$rho)-mean(log(mmat.ez.out$mvec))
    #
    # store estimates of eigenfunctions
    phi.cc[,i,j] <- est.cc$phi
    phi.star.cc[,i,j] <- est.cc$phi.star
    phi.ez[,i,j] <- est.ez$phi
    phi.star.ez[,i,j] <- est.ez$phi.star
    #
    # eigenfunction and eigenvalue for continuation value recursion
    efn.ez[,i,j] <- mmat.ez.out$efn.sim
    lambda.ez[i,j] <- mmat.ez.out$lambda
    #
  }
}
#
# population quantities
pars <- list()
pars$kappa <- kappa
pars$mu <- mu
pars$sigma <- sigma
pop.cc <- cc_linear_gaussian(beta,gamma,pars,g.grid)
pop.ez <- ez_linear_gaussian(beta,gamma,pars,g.grid)
#
phi.mean.cc <- phi.l.cc <- phi.u.cc <- phi.star.mean.cc <- phi.star.l.cc <- phi.star.u.cc <- phi.mean.ez <- phi.l.ez <- phi.u.ez <- phi.star.mean.ez <- phi.star.l.ez <- phi.star.u.ez <- efn.mean.ez <- efn.l.ez <- efn.u.ez <- matrix(NA,nrow=length(g.grid),ncol=length(n.vec))
#
# format estimates for computing bias and plotting
for(i in 1:length(n.vec)){
  phi.mean.cc[,i] <- matrix(rowMeans(phi.cc[,,i]))
  phi.l.cc[,i] <- matrix(rowQuantile(phi.cc[,,i],0.05))
  phi.u.cc[,i] <- matrix(rowQuantile(phi.cc[,,i],0.95))
  phi.star.mean.cc[,i] <- matrix(rowMeans(phi.star.cc[,,i]))
  phi.star.l.cc[,i] <- matrix(rowQuantile(phi.star.cc[,,i],0.05))
  phi.star.u.cc[,i] <- matrix(rowQuantile(phi.star.cc[,,i],0.95))
  phi.mean.ez[,i] <- matrix(rowMeans(phi.ez[,,i]))
  phi.l.ez[,i] <- matrix(rowQuantile(phi.ez[,,i],0.05))
  phi.u.ez[,i] <- matrix(rowQuantile(phi.ez[,,i],0.95))
  phi.star.mean.ez[,i] <- matrix(rowMeans(phi.star.ez[,,i]))
  phi.star.l.ez[,i] <- matrix(rowQuantile(phi.star.ez[,,i],0.05))
  phi.star.u.ez[,i] <- matrix(rowQuantile(phi.star.ez[,,i],0.95))
  efn.mean.ez[,i] <- matrix(rowMeans(efn.ez[,,i]))
  efn.l.ez[,i] <- matrix(rowQuantile(efn.ez[,,i],0.05))
  efn.u.ez[,i] <- matrix(rowQuantile(efn.ez[,,i],0.95))
}
#
# Bias and rmse
bias.func <- rmse.func <- matrix(NA,nrow=length(n.vec),ncol=5)
bias.pars <- rmse.pars <- matrix(NA,nrow=length(n.vec),ncol=7)
#
w.grid <- dnorm(g.grid,mean=mu,sd=sigma/sqrt(1-kappa^2))
w.grid <- w.grid/sum(w.grid)
#
for(i in 1:length(n.vec)){
  #
  # bias for parameters
  bias.pars[i,1] <- mean(r.cc[,i])-pop.cc$rho
  bias.pars[i,2] <- mean(y.cc[,i])-pop.cc$y
  bias.pars[i,3] <- mean(L.cc[,i])-pop.cc$L
  bias.pars[i,4] <- mean(r.ez[,i])-pop.ez$rho
  bias.pars[i,5] <- mean(y.ez[,i])-pop.ez$y
  bias.pars[i,6] <- mean(L.ez[,i])-pop.ez$L
  bias.pars[i,7] <- mean(lambda.ez[,i])-pop.ez$lambda
  #
  # rmse for parameters
  rmse.pars[i,1] <- sqrt(mean((r.cc[,i]-pop.cc$rho)^2))
  rmse.pars[i,2] <- sqrt(mean((y.cc[,i]-pop.cc$y)^2))
  rmse.pars[i,3] <- sqrt(mean((L.cc[,i]-pop.cc$L)^2))
  rmse.pars[i,4] <- sqrt(mean((r.ez[,i]-pop.ez$rho)^2))
  rmse.pars[i,5] <- sqrt(mean((y.ez[,i]-pop.ez$y)^2))
  rmse.pars[i,6] <- sqrt(mean((L.ez[,i]-pop.ez$L)^2))
  rmse.pars[i,7] <- sqrt(mean((lambda.ez[,i]-pop.ez$lambda)^2))
  #
  # bias for functions
  bias.func[i,1] <- sqrt(sum(w.grid*(phi.mean.cc[,i]-pop.cc$phi.sim)^2))
  bias.func[i,2] <- sqrt(sum(w.grid*(phi.star.mean.cc[,i]-pop.cc$phi.star.sim)^2))
  bias.func[i,3] <- sqrt(sum(w.grid*(phi.mean.ez[,i]-pop.ez$phi.sim)^2))
  bias.func[i,4] <- sqrt(sum(w.grid*(phi.star.mean.ez[,i]-pop.ez$phi.star.sim)^2))
  bias.func[i,5] <- sqrt(sum(w.grid*(efn.mean.ez[,i]-pop.ez$efn.sim)^2))
  #
  # rmse for functions
  rmse.func[i,1] <- rmse_mat(phi.cc[,,i],pop.cc$phi.sim,w.grid)
  rmse.func[i,2] <- rmse_mat(phi.star.cc[,,i],pop.cc$phi.star.sim,w.grid)
  rmse.func[i,3] <- rmse_mat(phi.ez[,,i],pop.ez$phi.sim,w.grid)
  rmse.func[i,4] <- rmse_mat(phi.star.ez[,,i],pop.ez$phi.star.sim,w.grid)
  rmse.func[i,5] <- rmse_mat(efn.ez[,,i],pop.ez$efn.sim,w.grid)
}
#
# save bias and rmse to file
if(basis==1){
  sink('npsdfd_sim_hp.txt')
}else{
  sink('npsdfd_sim_bs.txt')
}
cat("\n")
cat("Bias for phi and phistar (power utility), and phi, phistar, and chi (recursive)\n")
cat("\n")
for(i in 1:length(n.vec)){
  cat("n=",n.vec[i],"  ",format(round(bias.func[i,],digits=4),nsmall=4),"\n")
}
cat("\n")
cat("RMSE for phi and phistar (power utility), and phi, phistar, and chi (recursive)\n")
cat("\n")
for(i in 1:length(n.vec)){
  cat("n=",n.vec[i],"  ",format(round(rmse.func[i,],digits=4),nsmall=4),"\n")
}
cat("\n")
cat("Bias for rho, y, L (power utility), and rho, y, L, lambda (recursive)")
cat("\n")
for(i in 1:length(n.vec)){
  cat("n=",n.vec[i],"  ",format(round(bias.pars[i,],digits=4),nsmall=4),"\n")
}
cat("\n")
cat("RMSE for rho, y, L (power utility), and rho, y, L, lambda (recursive)")
cat("\n")
for(i in 1:length(n.vec)){
  cat("n=",n.vec[i],"  ",format(round(rmse.pars[i,],digits=4),nsmall=4),"\n")
}
cat("\n")
cat("\n")
sink()
#
### Plots ###
#
ix <- 1
# plot phi for power utility
phi <- pop.cc$phi.sim
phibar <- phi.mean.cc[,ix]
phil1 <- phi.l.cc[,1]
phiu1 <- phi.u.cc[,1]
phil2 <- phi.l.cc[,2]
phiu2 <- phi.u.cc[,2]
phil3 <- phi.l.cc[,3]
phiu3 <- phi.u.cc[,3]
pframe <- data.frame(g.grid,phi,phibar,phil1,phiu1,phil2,phiu2,phil3,phiu3) 
p.cc <- ggplot(pframe,aes(g.grid,phi))+
  geom_line(data=pframe,lwd=1.2,linetype=1)+
  geom_ribbon(data=pframe,aes(ymin=phil1,ymax=phiu1),alpha=0.15)+
  geom_ribbon(data=pframe,aes(ymin=phil2,ymax=phiu2),alpha=0.25)+
  geom_ribbon(data=pframe,aes(ymin=phil3,ymax=phiu3),alpha=0.35)+
  coord_cartesian(ylim=c(0.25,2))+
  theme_classic(base_size=20)+
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))
saveG <- TRUE
if(class(dev.list()) != "NULL"){dev.off()}
if(basis==1){
  if(saveG==TRUE){cairo_pdf(file="phi_cc_hp.pdf",height=9,width=13)}
}else{
  if(saveG==TRUE){cairo_pdf(file="phi_cc_bs.pdf",height=9,width=13)}
}
p.cc
if(saveG==TRUE){dev.off()}
#
# plot phistar for power utility
phi <- pop.cc$phi.star.sim
phibar <- phi.star.mean.cc[,ix]
phil1 <- phi.star.l.cc[,1]
phiu1 <- phi.star.u.cc[,1]
phil2 <- phi.star.l.cc[,2]
phiu2 <- phi.star.u.cc[,2]
phil3 <- phi.star.l.cc[,3]
phiu3 <- phi.star.u.cc[,3]
pframe <- data.frame(g.grid,phi,phibar,phil1,phiu1,phil2,phiu2,phil3,phiu3) 
p.star.cc <- ggplot(pframe,aes(g.grid,phi))+
  geom_line(data=pframe,lwd=1.2,linetype=1)+
  geom_ribbon(data=pframe,aes(ymin=phil1,ymax=phiu1),alpha=0.15)+
  geom_ribbon(data=pframe,aes(ymin=phil2,ymax=phiu2),alpha=0.25)+
  geom_ribbon(data=pframe,aes(ymin=phil3,ymax=phiu3),alpha=0.35)+
  coord_cartesian(ylim=c(0.25,3))+
  theme_classic(base_size=20)+
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))
saveG <- TRUE
if(class(dev.list()) != "NULL"){dev.off()}
if(basis==1){
  if(saveG==TRUE){cairo_pdf(file="phistar_cc_hp.pdf",height=9,width=13)}
}else{
  if(saveG==TRUE){cairo_pdf(file="phistar_cc_bs.pdf",height=9,width=13)}
}
p.star.cc
if(saveG==TRUE){dev.off()}
#
# plot phi for recursive preferences
phi <- pop.ez$phi.sim
phibar <- phi.mean.cc[,ix]
phil1 <- phi.l.ez[,1]
phiu1 <- phi.u.ez[,1]
phil2 <- phi.l.ez[,2]
phiu2 <- phi.u.ez[,2]
phil3 <- phi.l.ez[,3]
phiu3 <- phi.u.ez[,3]
pframe <- data.frame(g.grid,phi,phibar,phil1,phiu1,phil2,phiu2,phil3,phiu3) 
p.ez <- ggplot(pframe,aes(g.grid,phi))+
  geom_line(data=pframe,lwd=1.2,linetype=1)+
  geom_ribbon(data=pframe,aes(ymin=phil1,ymax=phiu1),alpha=0.15)+
  geom_ribbon(data=pframe,aes(ymin=phil2,ymax=phiu2),alpha=0.25)+
  geom_ribbon(data=pframe,aes(ymin=phil3,ymax=phiu3),alpha=0.35)+
  coord_cartesian(ylim=c(0.8,1.2))+
  theme_classic(base_size=20)+
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))
saveG <- TRUE
if(class(dev.list()) != "NULL"){dev.off()}
if(basis==1){
  if(saveG==TRUE){cairo_pdf(file="phi_ez_hp.pdf",height=9,width=13)}
}else{
  if(saveG==TRUE){cairo_pdf(file="phi_ez_bs.pdf",height=9,width=13)}
}
p.ez
if(saveG==TRUE){dev.off()}
#
# plot phistar for recursive preferences
phi <- pop.ez$phi.star.sim
phibar <- phi.star.mean.cc[,ix]
phil1 <- phi.star.l.ez[,1]
phiu1 <- phi.star.u.ez[,1]
phil2 <- phi.star.l.ez[,2]
phiu2 <- phi.star.u.ez[,2]
phil3 <- phi.star.l.ez[,3]
phiu3 <- phi.star.u.ez[,3]
pframe <- data.frame(g.grid,phi,phibar,phil1,phiu1,phil2,phiu2,phil3,phiu3) 
p.star.ez <- ggplot(pframe,aes(g.grid,phi))+
  geom_line(data=pframe,lwd=1.2,linetype=1)+
  geom_ribbon(data=pframe,aes(ymin=phil1,ymax=phiu1),alpha=0.15)+
  geom_ribbon(data=pframe,aes(ymin=phil2,ymax=phiu2),alpha=0.25)+
  geom_ribbon(data=pframe,aes(ymin=phil3,ymax=phiu3),alpha=0.35)+
  coord_cartesian(ylim=c(0,4.5))+
  theme_classic(base_size=20)+
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))
saveG <- TRUE
if(class(dev.list()) != "NULL"){dev.off()}
if(basis==1){
  if(saveG==TRUE){cairo_pdf(file="phistar_ez_hp.pdf",height=9,width=13)}
}else{
  if(saveG==TRUE){cairo_pdf(file="phistar_ez_bs.pdf",height=9,width=13)}
}
p.star.ez
if(saveG==TRUE){dev.off()}
#
# plot continuation value function for recursive preferences
efn <- pop.ez$efn.sim
efnbar <- efn.mean.ez[,ix]
efnl1 <- efn.l.ez[,1]
efnu1 <- efn.u.ez[,1]
efnl2 <- efn.l.ez[,2]
efnu2 <- efn.u.ez[,2]
efnl3 <- efn.l.ez[,3]
efnu3 <- efn.u.ez[,3]
pframe <- data.frame(g.grid,efn,efnbar,efnl1,efnu1,efnl2,efnu2,efnl3,efnu3) 
efn.plot.ez <- ggplot(pframe,aes(g.grid,efn))+
  geom_line(data=pframe,lwd=1.2,linetype=1)+
  geom_ribbon(data=pframe,aes(ymin=efnl1,ymax=efnu1),alpha=0.15)+
  geom_ribbon(data=pframe,aes(ymin=efnl2,ymax=efnu2),alpha=0.25)+
  geom_ribbon(data=pframe,aes(ymin=efnl3,ymax=efnu3),alpha=0.35)+
  coord_cartesian(ylim=c(0.25,2.0))+
  theme_classic(base_size=20)+
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))
saveG <- TRUE
if(class(dev.list()) != "NULL"){dev.off()}
if(basis==1){
  if(saveG==TRUE){cairo_pdf(file="chi_ez_hp.pdf",height=9,width=13)}
}else{
  if(saveG==TRUE){cairo_pdf(file="chi_ez_bs.pdf",height=9,width=13)}
}
efn.plot.ez
if(saveG==TRUE){dev.off()}
