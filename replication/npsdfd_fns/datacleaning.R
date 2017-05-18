# import monthly returns on fama french portfolios and convert to quarterly
import.ff <- function(){
  #
  data <- read.csv("./npsdfd_data/6_Portfolios_2x3.csv",header=FALSE,skip=20,nrows=1079)
  m.ret <- (data[-nrow(data),(2:ncol(data))])
  m.ret <- 1+m.ret/100
  #
  cumulative.m.ret <- matrix(NA,nrow=nrow(m.ret),ncol=ncol(m.ret))
  for(i in 1:ncol(m.ret)){
    cumulative.m.ret[,i] <- cumprod(m.ret[,i])
  }
  #
  date.seq <- seq(from=3,to=nrow(m.ret),by=3)
  date.seq.0 <- date.seq[-length(date.seq)]
  date.seq.1 <- date.seq[-1]
  q.ret <- matrix(NA,nrow=length(date.seq)-1,ncol=ncol(m.ret))
  for(i in 1:ncol(m.ret)){
    q.ret[,i] <- cumulative.m.ret[date.seq.1,i]/cumulative.m.ret[date.seq.0,i]
  }
  return(q.ret)
}
#
# import returns on fama french portfolios and trim to sample 1947:Q1-2016:Q1
import.ff.trim <- function(){
  #
  q.ret <- import.ff()
  q.ret <- q.ret[(82:358),]
  return(q.ret)
}
#
# import 90 day t-bill rate
import.t90d <- function(){
  #
  data <- read.csv("./npsdfd_data/TB3MS.csv",header=FALSE,skip=1)
  #
  # trim to 1947:Q1-2016:Q1 (note that the FRED aggregation method was set to the end of the period)
  t90d <- as.numeric(data[(52:328),2])
  #
  # convert to quarterly gross return
  t90d <- exp(t90d/400)
  #
  return(t90d)
}
#
# import nipa consumption data and deflate
import.nipa <- function(){
  #
  # pre-1970 nipa data
  data_nd_0 <- read.csv("./npsdfd_data/Section2All_Hist.csv",header=FALSE,skip=664,nrows=1,colClasses="character")
  data_s_0 <- read.csv("./npsdfd_data/Section2All_Hist.csv",header=FALSE,skip=669,nrows=1,colClasses="character")
  data_ipd_0 <- read.csv("./npsdfd_data/Section2All_Hist.csv",header=FALSE,skip=552,nrows=1,colClasses="character")
  data_pop_0 <- read.csv("./npsdfd_data/Section2All_Hist.csv",header=FALSE,skip=117,nrows=1,colClasses="character")
  data_ce_0 <- read.csv("./npsdfd_data/Section1All_Hist.csv",header=FALSE,skip=4232,nrows=1,colClasses="character")
  data_div_0 <- read.csv("./npsdfd_data/Section1All_Hist.csv",header=FALSE,skip=4233,nrows=1,colClasses="character")
  data_gdp_0 <- read.csv("./npsdfd_data/Section1All_Hist.csv",header=FALSE,skip=345,nrows=1,colClasses="character")
  #
  # post-1969 nipa data
  data_nd_1 <- read.csv("./npsdfd_data/Section2all_csv.csv",header=FALSE,skip=715,nrows=1,colClasses="character")
  data_s_1 <- read.csv("./npsdfd_data/Section2all_csv.csv",header=FALSE,skip=720,nrows=1,colClasses="character")
  data_ipd_1 <- read.csv("./npsdfd_data/Section2all_csv.csv",header=FALSE,skip=603,nrows=1,colClasses="character")
  data_pop_1 <- read.csv("./npsdfd_data/Section2all_csv.csv",header=FALSE,skip=117,nrows=1,colClasses="character")
  data_ce_1 <- read.csv("./npsdfd_data/Section1all_csv.csv",header=FALSE,skip=4399,nrows=1,colClasses="character")
  data_div_1 <- read.csv("./npsdfd_data/Section1all_csv.csv",header=FALSE,skip=4400,nrows=1,colClasses="character")
  data_gdp_1 <- read.csv("./npsdfd_data/Section1all_csv.csv",header=FALSE,skip=345,nrows=1,colClasses="character")
  #
  # trim
  data_nd_0 <- as.numeric(gsub(",", "", as.character(data_nd_0[,(4:95)])))
  data_s_0 <- as.numeric(gsub(",", "", as.character(data_s_0[,(4:95)])))
  data_ipd_0 <- as.numeric(gsub(",", "", as.character(data_ipd_0[,(4:95)])))
  data_pop_0 <- as.numeric(gsub(",", "", as.character(data_pop_0[,(4:95)])))
  data_ce_0 <- as.numeric(gsub(",", "", as.character(data_ce_0[,(4:95)])))
  data_div_0 <- as.numeric(gsub(",", "", as.character(data_div_0[,(4:95)])))
  data_gdp_0 <- as.numeric(gsub(",", "", as.character(data_gdp_0[,(4:95)])))
  #
  data_nd_1 <- as.numeric(gsub(",", "", as.character(data_nd_1[,(8:192)])))
  data_s_1 <- as.numeric(gsub(",", "", as.character(data_s_1[,(8:192)])))
  data_ipd_1 <- as.numeric(gsub(",", "", as.character(data_ipd_1[,(8:192)])))
  data_pop_1 <- as.numeric(gsub(",", "", as.character(data_pop_1[,(8:192)])))
  data_ce_1 <- as.numeric(gsub(",", "", as.character(data_ce_1[,(8:192)])))
  data_div_1 <- as.numeric(gsub(",", "", as.character(data_div_1[,(8:192)])))
  data_gdp_1 <- as.numeric(gsub(",", "", as.character(data_gdp_1[,(8:192)])))
  #
  # concatenate and convert to appropriate format
  data_nd <- c(data_nd_0,data_nd_1)
  data_s <- c(data_s_0,data_s_1)
  data_ipd <- c(data_ipd_0,data_ipd_1)
  data_pop <- c(data_pop_0,data_pop_1)
  data_ce <- c(data_ce_0,data_ce_1)
  data_div <- c(data_div_0,data_div_1)
  data_gdp <- c(data_gdp_0,data_gdp_1)
  data_div_adj <- data_div
  for(i in 4:length(data_div_adj)){
    data_div_adj[i] <- prod(data_div[(i-3):i])^(1/4)
  }
  #
  # convert to growth rates
  data_c <- data_nd+data_s
  g_c_n <- log(data_c[-1])-log(data_c[-length(data_c)])
  g_ce_n <- log(data_ce[-1])-log(data_ce[-length(data_ce)])
  g_div_n <- log(data_div[-1])-log(data_div[-length(data_div)])
  g_div_adj_n <- log(data_div_adj[-1])-log(data_div_adj[-length(data_div_adj)])
  g_gdp_n <- log(data_gdp[-1])-log(data_gdp[-length(data_gdp)])
  g_ipd <- log(data_ipd[-1])-log(data_ipd[-length(data_ipd)])   # PCE implicit price deflator growth
  g_pop <- log(data_pop[-1])-log(data_pop[-length(data_pop)])   # population growth
  g_c <- g_c_n-g_ipd-g_pop        # real per capita nondurables+services consumption growth 
  g_ce <- g_ce_n-g_ipd-g_pop      # real per capita corporate earnings growth
  g_div <- g_div_n-g_ipd-g_pop    # real per capita dividend growth
  g_div_adj <- g_div_adj_n-g_ipd-g_pop    # real per capita dividend growth adjusted by 4 quarter geom. average of dividends
  g_gdp <- g_gdp_n-g_ipd-g_pop    # real per capita GDP growth
  #
  return(list(g_c=g_c,g_ipd=g_ipd,g_ce=g_ce,g_div=g_div,g_gdp=g_gdp,g_div_adj=g_div_adj))
}


