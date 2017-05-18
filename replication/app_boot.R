# Run this file to perform application 
#
rm(list=ls())
#
# load packages and codes
library("MASS")
library("QZ")
library("ggplot2")
library("scales")
library("reshape2")
library("grid")
library("tis")
source("./npsdfd_fns/orthopol.R")
source("./npsdfd_fns/matfns.R")
source("./npsdfd_fns/Mest.R")
source("./npsdfd_fns/lgsoln.R")
source("./npsdfd_fns/nlfp.R")
source("./npsdfd_fns/datacleaning.R")
source("./npsdfd_fns/ezestfns.R")
source("./npsdfd_fns/bootfns.R")
#
# import NIPA data
nipa_data <- import.nipa()
g <- nipa_data$g_c
#
### univariate case ###
#
# grid for plotting
n <- length(g)
mu <- mean(g)
sigma <- sd(g)
g.grid <- seq(-2*sigma+mu,2*sigma+mu,length.out=2e2)
#
# generate basis
k <- 8
B.mat <- genhermite(g,k,mu=mu,sigma=sigma)
B.inst <- genhermite(g,k-2,mu=mu,sigma=sigma)
B.sim <- genhermite(g.grid,k,mu=mu,sigma=sigma)
Ghat <- t(B.mat[1:n-1,])%*%B.mat[1:n-1,]/n
#
# estimate preference parameters
# import returns data
R_t90d <- import.t90d()   # 90 day t-bill rate
R_ff <- import.ff.trim()  # 6 Fama-French portfolios
R <- cbind(R_ff,R_t90d)
R <- R[-1,]
# convert to real returns using PCE deflator
g_ipd <- exp(nipa_data$g_ipd)
for(i in 1:ncol(R)){
  R[,i] <- R[,i]/g_ipd
}
#
# initial parameters for optimization
beta0 <- 0.99; gamma0 <- 30
#
theta_hat <- est.ez.fn.flex(beta0,gamma0,Ghat,g,B.mat,B.inst,R[-1,])
beta_hat <- theta_hat$beta
gamma_hat <- theta_hat$gamma
#
# decomposition of estimated SDF
ez.out <- ez.decomp(g,B.mat,Ghat,beta_hat,gamma_hat,niter=250,B.sim)
#
# plots of eigenfunction estimates
pframe <- data.frame(g.grid,ez.out$phi,ez.out$phi.star,ez.out$chi,ez.out$phi*ez.out$phi.star)
colnames(pframe) <- c("g","phi_hat","phi_star_hat","chi_hat","rn_deriv")
#
phi.plot.1 <- ggplot(pframe)+
  geom_line(data=pframe,aes(g,phi_hat),lwd=1.2,linetype=1)+
  ylab(expression(hat(phi)(x)))+
  xlab(expression(x))+
  theme_classic()+
  theme(legend.position="none",
        axis.text.x=element_text(size=20),
        axis.text.y=element_text(size=20),
        axis.title.x=element_text(size=20),
        axis.title.y=element_text(size=20))
saveG <- TRUE
if(class(dev.list()) != "NULL"){dev.off()}
if(saveG==TRUE){cairo_pdf(file="phi_emp_1.pdf",height=9,width=13)}
phi.plot.1
if(saveG==TRUE){dev.off()}
#
phi.star.plot.1 <- ggplot(pframe)+
  geom_line(data=pframe,aes(g,phi_star_hat),lwd=1.2,linetype=1)+
  ylab(c(expression(paste(hat(phi),"*",(x)))))+
  xlab(expression(x))+
  theme_classic()+
  theme(legend.position="none",
        axis.text.x=element_text(size=20),
        axis.text.y=element_text(size=20),
        axis.title.x=element_text(size=20),
        axis.title.y=element_text(size=20))
saveG <- TRUE
if(class(dev.list()) != "NULL"){dev.off()}
if(saveG==TRUE){cairo_pdf(file="phistar_emp_1.pdf",height=9,width=13)}
phi.star.plot.1
if(saveG==TRUE){dev.off()}
#
rn.plot.1 <- ggplot(pframe)+
  geom_line(data=pframe,aes(g,rn_deriv),lwd=1.2,linetype=1)+
  ylab(c(expression(paste(hat(phi)(x),hat(phi),"*",(x)))))+
  xlab(expression(x))+
  theme_classic()+
  theme(legend.position="none",
        axis.text.x=element_text(size=20),
        axis.text.y=element_text(size=20),
        axis.title.x=element_text(size=20),
        axis.title.y=element_text(size=20))
saveG <- TRUE
if(class(dev.list()) != "NULL"){dev.off()}
if(saveG==TRUE){cairo_pdf(file="rn_emp_1.pdf",height=9,width=13)}
rn.plot.1
if(saveG==TRUE){dev.off()}
#
# bootstrap for inference
# initialize
alpha <- 0.10
nboot <- 1000
blocklength <- 6
#
# bootstrap, estimating preference parameters for each iteration 
# may take a few minutes
bs.est <- ez.bs(g,B.mat,nboot,blocklength,beta0,gamma0,est_pref_pars=TRUE,B.inst,R)
#
# save output to file
sink('npsdfd_app_output.txt',append=FALSE)
cat("estimates for the univariate case: X_t = g_t\n")
cat("\n")
cat("With estimated preference parameters \n")
cat("\n")
cat("rho_hat:",ez.out$rho,"\n")
cat("y_hat:",ez.out$y,"\n")
cat("L_hat:",ez.out$L,"\n")
cat("beta_hat:",beta_hat,"\n")
cat("gamma_hat:",gamma_hat,"\n")
cat("lambda_hat",ez.out$lambda,"\n")
cat("\n")
cat("SDF entropy:",log(mean(ez.out$mvec)) - mean(log(ez.out$mvec)),"\n")
cat("\n")
cat("\\underset{(",bsci(ez.out$rho,bs.est$rho_boot[bs.est$gamma_boot>1.001],alpha)$cl,",",bsci(ez.out$rho,bs.est$rho_boot[bs.est$gamma_boot>1.001],alpha)$cu,")}{",ez.out$rho,"}\n")
cat("\\underset{(",bsci(ez.out$y,bs.est$y_boot[bs.est$gamma_boot>1.001],alpha)$cl,",",bsci(ez.out$y,bs.est$y_boot[bs.est$gamma_boot>1.001],alpha)$cu,")}{",ez.out$y,"}\n")
cat("\\underset{(",bsci(ez.out$L,bs.est$L_boot[bs.est$gamma_boot>1.001],alpha)$cl,",",bsci(ez.out$L,bs.est$L_boot[bs.est$gamma_boot>1.001],alpha)$cu,")}{",ez.out$L,"}\n")
cat("\\underset{(",bsci(beta_hat,bs.est$beta_boot[bs.est$gamma_boot>1.001],alpha)$cl,",",bsci(beta_hat,bs.est$beta_boot[bs.est$gamma_boot>1.001],alpha)$cu,")}{",beta_hat,"}\n")
cat("\\underset{(",bsci(gamma_hat,bs.est$gamma_boot[bs.est$gamma_boot>1.001],alpha)$cl,",",bsci(gamma_hat,bs.est$gamma_boot[bs.est$gamma_boot>1.001],alpha)$cu,")}{",gamma_hat,"}\n")
cat("\\underset{(",bsci(ez.out$lambda,bs.est$lambda_boot[bs.est$gamma_boot>1.001],alpha)$cl,",",bsci(ez.out$lambda,bs.est$lambda_boot[bs.est$gamma_boot>1.001],alpha)$cu,")}{",ez.out$lambda,"}\n")
cat("\n")
sink()
#
### bivariate case ###
#
# initialize
d <- nipa_data$g_div
#
k1 <- 5
k2 <- 5
d.grid <- seq(-2*sd(d)+mean(d),2*sd(d)+mean(d),length.out=3e1)
g.grid <- seq(-2*sd(g)+mean(g),2*sd(g)+mean(g),length.out=2e1)
B.mat.1 <- genhermite(g,k1,mu=mu,sigma=sigma)
B.inst.1 <- genhermite(g,k1-2,mu=mu,sigma=sigma)
B.sim.1 <- genhermite(g.grid,k1,mu=mu,sigma=sigma)
B.mat.2 <- genhermite(d,k2,mu=mean(d),sigma=sd(d))
B.inst.2 <- genhermite(d,k2-2,mu=mean(d),sigma=sd(d))
B.sim.2 <- genhermite(d.grid,k2,mu=mean(d),sigma=sd(d))
#
# sparse tensor-product basis
B.mat <- sp(B.mat.1,B.mat.2)
B.inst <- sp(B.inst.1,B.inst.2)
B.sim <- sp.kronecker(B.sim.1,B.sim.2)
Ghat <- t(B.mat[1:n-1,])%*%B.mat[1:n-1,]/n
#
# estimate preference parameters
theta_hat <- est.ez.fn.flex(beta0,gamma0,Ghat,g,B.mat,B.inst,R[-1,])
beta_hat <- theta_hat$beta
gamma_hat <- theta_hat$gamma
#
# decomposition of estimated SDF
ez.out <- ez.decomp(g,B.mat,Ghat,beta_hat,gamma_hat,niter=250,B.sim)
#
# contour plots of eigenfunction estimates
pframe <- data.frame(kronecker(g.grid,matrix(1,nrow=length(d.grid),ncol=1)),kronecker(matrix(1,nrow=length(g.grid),ncol=1),d.grid),ez.out$phi)
colnames(pframe) <- c("g","d","z")
phi.plot.2 <- ggplot(pframe,aes(g,d,z=z))+
  geom_tile(aes(fill = z))+
  scale_fill_gradient2(low = "white", high = "black", midpoint=1)+
  stat_contour(bins=25,color="black",lwd=1.0)+
  xlab(expression(g))+
  ylab(expression(d))+
  theme_classic()+
  theme(legend.title=element_blank(),
        legend.text = element_text(size = 28),
        axis.text.x=element_text(size=20),
        axis.text.y=element_text(size=20),
        axis.title.x=element_text(size=20),
        axis.title.y=element_text(size=20))+
  guides(fill = guide_colorbar(barwidth = 2, barheight = 10))
saveG <- TRUE
if(class(dev.list()) != "NULL"){dev.off()}
if(saveG==TRUE){cairo_pdf(file="phi_emp_2.pdf",height=9,width=13)}
phi.plot.2
if(saveG==TRUE){dev.off()}
#
pframe <- data.frame(kronecker(g.grid,matrix(1,nrow=length(d.grid),ncol=1)),kronecker(matrix(1,nrow=length(g.grid),ncol=1),d.grid),ez.out$phi.star)
colnames(pframe) <- c("g","d","z")
phi.star.plot.2 <- ggplot(pframe,aes(g,d,z=z))+
  geom_tile(aes(fill = z))+
  scale_fill_gradient2(low = "white", high = "black", midpoint=1)+
  stat_contour(bins=25,color="black",lwd=1.0)+
  xlab(expression(g))+
  ylab(expression(d))+
  theme_classic()+
  theme(legend.title=element_blank(),
        legend.text = element_text(size = 28),
        axis.text.x=element_text(size=20),
        axis.text.y=element_text(size=20),
        axis.title.x=element_text(size=20),
        axis.title.y=element_text(size=20))+
  guides(fill = guide_colorbar(barwidth = 2, barheight = 10))
saveG <- TRUE
if(class(dev.list()) != "NULL"){dev.off()}
if(saveG==TRUE){cairo_pdf(file="phistar_emp_2.pdf",height=9,width=13)}
phi.star.plot.2
  if(saveG==TRUE){dev.off()}
#
pframe <- data.frame(kronecker(g.grid,matrix(1,nrow=length(d.grid),ncol=1)),kronecker(matrix(1,nrow=length(g.grid),ncol=1),d.grid),ez.out$phi*ez.out$phi.star)
colnames(pframe) <- c("g","d","z")
rn.plot.2 <- ggplot(pframe,aes(g,d,z=z))+
  geom_tile(aes(fill = z))+
  scale_fill_gradient2(low = "white", high = "black", midpoint=1)+
  stat_contour(bins=25,color="black",lwd=1.0)+
  xlab(expression(g))+
  ylab(expression(d))+
  theme_classic()+
  theme(legend.title=element_blank(),
        legend.text = element_text(size = 28),
        axis.text.x=element_text(size=20),
        axis.text.y=element_text(size=20),
        axis.title.x=element_text(size=20),
        axis.title.y=element_text(size=20))+
  guides(fill = guide_colorbar(barwidth = 2, barheight = 10))
saveG <- TRUE
if(class(dev.list()) != "NULL"){dev.off()}
if(saveG==TRUE){cairo_pdf(file="rn_emp_2.pdf",height=9,width=13)}
rn.plot.2
if(saveG==TRUE){dev.off()}
#
# time series plots
# permanent and transitory components for recursive preferences
pctc <- pctc.series(ez.out$mvec,B.mat,ez.out$rho,ez.out$c)
pc.ez.est.2 <- pctc$pc
tc.ez.est.2 <- pctc$tc
#
# do decomposition with power utility, using same estimated preference parameters
mmat.cc.out <- mmat_ccapm(g,B.mat,beta_hat,gamma_hat)
Mhat.cc <- mmat.cc.out$Mhat
est.cc <- npeigest(Mhat.cc,Ghat,B.sim,B.mat)
#
# permanent and transitory components for power utility
pctc <- pctc.series(mmat.cc.out$mvec,B.mat,est.cc$rho,est.cc$c.hat)
pc.cc.est.2 <- pctc$pc
tc.cc.est.2 <- pctc$tc
#
df <- data.frame(ez.out$mvec,pc.ez.est.2,mmat.cc.out$mvec,pc.cc.est.2)
colnames(df) <- c("SDF\n(EZ)","PC\n(EZ)","SDF\n(PU)","PC\n(PU)")
df$dates <- seq(as.Date('1947-07-01'), as.Date('2016-01-01'), by = "quarter")
df.melt <- melt(df, id.vars = "dates")
#
ts0 <- ggplot(df.melt,aes(x=dates,y=value))+ 
  geom_line() +
  facet_grid(variable ~ .) 
yrng <- extendrange(eval(ts0$mapping$y, ts0$data))
rm(ts0)
nber.dates <- data.frame(as.data.frame(nberDates()),ymin=yrng[1],ymax=yrng[2])
nber.dates$Start <- as.Date(as.character(nber.dates$Start),format="%Y%m%d")
nber.dates$End <- as.Date(as.character(nber.dates$End),format="%Y%m%d")
nber.dates <- subset(nber.dates, Start >= as.Date('1947-10-01'))
ts1 <- ggplot(df.melt,aes(x=dates,y=value))+ 
  geom_line(lwd=1)+
  facet_grid(variable ~ .)+
  geom_rect(aes(xmin=Start,xmax=End,x=NULL,y=NULL,ymin=ymin,ymax=ymax),data=nber.dates,alpha = 0.25)+
  theme_classic()+
  scale_x_date(breaks=c(as.Date("1950-01-01"),
                        as.Date("1960-01-01"),
                        as.Date("1970-01-01"),
                        as.Date("1980-01-01"),
                        as.Date("1990-01-01"),
                        as.Date("2000-01-01"),
                        as.Date("2010-01-01")),labels = date_format("%Y"))+
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_text(size=14),
        axis.text.y=element_text(size=12),
        strip.text.y = element_text(size = 16,angle=0),
        strip.background = element_rect(fill = "white", color = "white", size = 1))
saveG <- TRUE
if(class(dev.list()) != "NULL"){dev.off()}
if(saveG==TRUE){cairo_pdf(file="tsplots.pdf",height=9,width=13)}
ts1
if(saveG==TRUE){dev.off()}
#
# bootstrap, estimating preference parameters for each iteration (may take a few minutes)
bs.est <- ez.bs(g,B.mat,nboot,blocklength,beta0,gamma0,est_pref_pars=TRUE,B.inst,R)
#
# save output to file
sink('npsdfd_app_output.txt',append=TRUE)
cat("estimates for the multivariate case: X_t = (g_t,d_t) \n")
cat("\n")
cat("With estimated preference parameters \n")
cat("\n")
cat("rho_hat:",ez.out$rho,"\n")
cat("y_hat:",ez.out$y,"\n")
cat("L_hat:",ez.out$L,"\n")
cat("beta_hat:",beta_hat,"\n")
cat("gamma_hat:",gamma_hat,"\n")
cat("lambda_hat",ez.out$lambda,"\n")
cat("\n")
cat("SDF entropy:",log(mean(ez.out$mvec)) - mean(log(ez.out$mvec)),"\n")
cat("\n")
cat("correlation of PC and GDP growth:",cor(nipa_data$g_gdp[-1],pc.ez.est.2),"\n")
cat("correlation of TC and GDP growth:",cor(nipa_data$g_gdp[-1],tc.ez.est.2),"\n")
cat("\n")
cat("\\underset{(",bsci(ez.out$rho,bs.est$rho_boot[bs.est$gamma_boot>1.001],alpha)$cl,",",bsci(ez.out$rho,bs.est$rho_boot[bs.est$gamma_boot>1.001],alpha)$cu,")}{",ez.out$rho,"}\n")
cat("\\underset{(",bsci(ez.out$y,bs.est$y_boot[bs.est$gamma_boot>1.001],alpha)$cl,",",bsci(ez.out$y,bs.est$y_boot[bs.est$gamma_boot>1.001],alpha)$cu,")}{",ez.out$y,"}\n")
cat("\\underset{(",bsci(ez.out$L,bs.est$L_boot[bs.est$gamma_boot>1.001],alpha)$cl,",",bsci(ez.out$L,bs.est$L_boot[bs.est$gamma_boot>1.001],alpha)$cu,")}{",ez.out$L,"}\n")
cat("\\underset{(",bsci(beta_hat,bs.est$beta_boot[bs.est$gamma_boot>1.001],alpha)$cl,",",bsci(beta_hat,bs.est$beta_boot[bs.est$gamma_boot>1.001],alpha)$cu,")}{",beta_hat,"}\n")
cat("\\underset{(",bsci(gamma_hat,bs.est$gamma_boot[bs.est$gamma_boot>1.001],alpha)$cl,",",bsci(gamma_hat,bs.est$gamma_boot[bs.est$gamma_boot>1.001],alpha)$cu,")}{",gamma_hat,"}\n")
cat("\\underset{(",bsci(ez.out$lambda,bs.est$lambda_boot[bs.est$gamma_boot>1.001],alpha)$cl,",",bsci(ez.out$lambda,bs.est$lambda_boot[bs.est$gamma_boot>1.001],alpha)$cu,")}{",ez.out$lambda,"}\n")
cat("\n")
sink()
#
# bootstrap with fixed preference parameters
beta_star <- 0.99
gamma_star <- c(20,25,30)
#
for(i in 1:length(gamma_star)){
  #
  cat("\n")
  cat("bootstrapping for gamma = ",gamma_star[i],"\n")
  cat("\n")
  #
  # bootstrap
  bs.fix <- ez.bs(g,B.mat,nboot,blocklength,beta_star,gamma_star[i],est_pref_pars=FALSE)
  #
  # SDF decomp at fixed preference parameters
  ez.out.fix <- ez.decomp(g,B.mat,Ghat,beta_star,gamma_star[i],niter=250,B.sim)
  #
  # print output
  sink('npsdfd_app_output.txt',append=TRUE)
  cat("\n")
  cat("With fixed preference parameters: beta_star =",beta_star,", gamma_star =",gamma_star[i],"\n")
  cat("rho_hat:",ez.out.fix$rho,"\n")
  cat("y_hat:",ez.out.fix$y,"\n")
  cat("L_hat:",ez.out.fix$L,"\n")
  cat("lambda_hat",ez.out.fix$lambda,"\n")
  cat("\n")
  cat("\\underset{(",bsci(ez.out.fix$rho,bs.fix$rho_boot,alpha)$cl,",",bsci(ez.out.fix$rho,bs.fix$rho_boot,alpha)$cu,")}{",ez.out.fix$rho,"}\n")
  cat("\\underset{(",bsci(ez.out.fix$y,bs.fix$y_boot,alpha)$cl,",",bsci(ez.out.fix$y,bs.fix$y_boot,alpha)$cu,")}{",ez.out.fix$y,"}\n")
  cat("\\underset{(",bsci(ez.out.fix$L,bs.fix$L_boot,alpha)$cl,",",bsci(ez.out.fix$L,bs.fix$L_boot,alpha)$cu,")}{",ez.out.fix$L,"}\n")
  cat("\\underset{(",bsci(ez.out.fix$lambda,bs.fix$lambda_boot,alpha)$cl,",",bsci(ez.out.fix$lambda,bs.fix$lambda_boot,alpha)$cu,")}{",ez.out.fix$lambda,"}\n")
  cat("\n")
  sink()
}
#
# generate figures 4a and 4b
# parameterize
beta <- 0.994
gamma_grid <- seq(1,35,length.out=35)
#
# estimate for different values of gamma
rho.ez <- y.ez <- L.ez <- cor.ez <- numeric(length(gamma_grid))
for(i in 1:length(gamma_grid)){
  #
  ez.par.out <- ez.decomp(g,B.mat,Ghat,beta,gamma_grid[i],niter=250,B.sim)
  #
  rho.ez[i] <- ez.par.out$rho
  y.ez[i] <- ez.par.out$y
  L.ez[i] <- ez.par.out$L
  #
  # correlation between recovered PC and TC
  pctc.ez <- pctc.series(ez.par.out$mvec,B.mat,ez.par.out$rho,ez.par.out$c)
  cor.ez[i] <- cor(log(pctc.ez$pc),log(pctc.ez$tc))
}
#
# compare with linear Gaussian AR(1) case
rho.lg.ez <- y.lg.ez <- L.lg.ez <- cor.lg.ez <- numeric(length(gamma_grid))
pars_var <- var1.est(rbind(t(g),t(d)))
for(i in 1:length(gamma_grid)){
  ez_lg_out <- ez_linear_gaussian_mv(beta,gamma_grid[i],pars_var)
  rho.lg.ez[i] <- ez_lg_out$rho
  y.lg.ez[i] <- ez_lg_out$y
  L.lg.ez[i] <- ez_lg_out$L
  cor.lg.ez[i] <- ez_lg_out$cor.pt
}
#
# final graph to compare with stoch. vol case.
# note this will take some time to run as using simulation-based estimation
# library("rugarch")
# pars.ar1arg <- ar1arg1.ii.est(g)
# to skip this estimation step, use 
pars.ar1arg <- list(mu=0.004140928,kappa=0.2850159,phi=0.5667164,c=9.448459e-05,delta=0.1729404)
ez.sv.out <- ez_sv_emp(beta,gamma_grid,pars.ar1arg)
y.ez.sv <- ez.sv.out$y.lgsv
cor.ez.sv <- ez.sv.out$cor.lgsv.ez
L.ez.sv <- ez.sv.out$L.lgsv.ez
#
pframe <- data.frame(gamma_grid,y.ez,y.lg.ez,y.ez.sv,cor.ez,cor.lg.ez,cor.ez.sv) 
#
yield.plot.ez <- ggplot(pframe)+
  geom_line(aes(gamma_grid,y.ez),lwd=1.2,linetype=1)+
  geom_line(aes(gamma_grid,y.lg.ez),lwd=1.2,linetype=2)+
  geom_line(aes(gamma_grid,y.ez.sv),lwd=1.2,linetype=3)+
  ylab(expression(y))+
  xlab(expression(gamma))+
  theme_classic()+
  theme(axis.title.x=element_text(size=20),
        axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=20),
        axis.text.y=element_text(size=20))
if(class(dev.list()) != "NULL"){dev.off()}
if(saveG==TRUE){cairo_pdf(file="emp_yld_sv.pdf",height=9,width=13)}
yield.plot.ez
if(saveG==TRUE){dev.off()}
#
cor.plot.ez <- ggplot(pframe)+
  geom_line(aes(gamma_grid,cor.ez),lwd=1.2,linetype=1)+
  geom_line(aes(gamma_grid,cor.lg.ez),lwd=1.2,linetype=2)+
  geom_line(aes(gamma_grid,cor.ez.sv),lwd=1.2,linetype=3)+
  ylab(expression(correlation))+
  xlab(expression(gamma))+
  theme_classic()+
  theme(axis.title.x=element_text(size=20),
        axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=20),
        axis.text.y=element_text(size=20))
if(class(dev.list()) != "NULL"){dev.off()}
if(saveG==TRUE){cairo_pdf(file="emp_cor_sv.pdf",height=9,width=13)}
cor.plot.ez
if(saveG==TRUE){dev.off()}









