library(ggplot2)
library(data.table)
library(tseries)
library(forecast)
library(lmtest)
library(urca)

rm(list=ls())
gc()

# the following function selects the autoregressive order based on the 'SIC rule' as follows: 
# (1) find the order, p (Schwarz Information Criterion); 
# (2) test for serial correlation in residuals (Breusch-Godfrey test);
# (3) increase the order by 1, if we reject the null of no serial correlation;
# (4) repeat steps 2-3 until we fail to reject the null or if the autoregressive order, lag_opt, equals the predetermined maximum order, lag_max (set to 12).
ar_order <- function(x,max_lag=12,method="SIC",bg=TRUE){
  ic_func <- function(i){# optimal lag based on (A/S)IC
    fmla <- as.formula(paste("x~",paste0("x",c(1:i),collapse="+")))
    reg <- lm(fmla,data=d)
    if(method=="SIC"){
      log(crossprod(reg$resid))+log(nrow(d))*length(reg$coefficients)/nrow(d)
    }else{
      log(crossprod(reg$resid))+2*length(reg$coefficients)/nrow(d)
    }
  }
  bg_func <- function(i){# Breusch-Godfrey
    fmla <- as.formula(paste("x~",paste0("x",c(1:i),collapse="+")))
    as.numeric(bgtest(lm(fmla,data=d))[4])
  }
  d <- data.table(x=as.numeric(x))
  d <- d[,paste0("x",c(1:max_lag)) := shift(x,n=1:max_lag,fill=NA,type="lag")]
  d <- d[complete.cases(d)]
  ic_vec <- sapply(c(1:max_lag),ic_func)
  opt_lag <- which.min(ic_vec)
  if(bg){
    bg_vec1 <- sapply(c(opt_lag:1),bg_func)
    if(!is.na(which(bg_vec1 < .05)[1])){
      opt_lag <- opt_lag-which(bg_vec1 < .05)[1]+2
    }
    bg_vec2 <- sapply(c(opt_lag:max_lag),bg_func)
    if(!is.na(which(bg_vec2 < .05)[1])){
      if(which(bg_vec2 < .05)[1]==1){
        opt_lag <- opt_lag+sort(unlist(lapply(split(which(bg_vec2 < .05),cumsum(c(1,diff(which(bg_vec2 < .05)))!=1)),max),use.names=F))[1]
      }
    }
  }
  p <- ifelse(p<=max_lag,opt_lag,max_lag)
  return(p)
}

# load the function to test linearity against star-type alternatives
source("startest.r")

# load the price data
load("dataset.RData")

c.num <- ncol(prices_dt)-1

p <- 12
h <- 12

frac <- .8

n <- nrow(prices_dt)
R <- round(frac*n)
P <- n-R-h+1

c.names <- colnames(prices_dt)[-1]

select_ls <- list()

for(r in 1:P){
  
  select_dt <- data.table(c=c.names,i=as.numeric(NA),p=as.numeric(NA),d=as.numeric(NA),n=as.numeric(NA),t=as.numeric(NA),l=as.numeric(NA))
  
  for(k in 1:c.num){
    
    cmd <- c.names[k]
    
    y <- as.matrix(log(prices_dt[1:(R+r-1),..cmd]))
    
    ord <- ar_order(y,max_lag=p,method="SIC",bg=F)
    
    adf.l <- ur.df(y,type="drift",lags=ord)
    
    adf.ts <- adf.l@teststat["statistic","tau2"]
    adf.cv <- adf.l@cval["tau2","5pct"]
    
    select_dt[,2][k] <- ifelse(adf.ts < adf.cv,0,1)
    
    dy <- c(0,diff(y))
    ord <- ord-1
    
    dmat <- embed(dy,p+h+1)
    lmat <- embed(y,p+h+1)
    
    if(ord==0){
      test.sc <- startest(y=dmat[,1],x=rep(1,nrow(dmat)),x.n=rep(1,nrow(dmat)),tvar=dmat[,2],tvar.lag=0,contemp=F,ascending=T)
    }else{
      test.sc <- startest(y=dmat[,1],x=cbind(1,dmat[,2:(2+ord-1)]),x.n=cbind(1,dmat[,2:(2+ord-1)]),tvar=dmat[,2:(2+ord-1)],tvar.lag=0,contemp=F,ascending=T)
    }

    d <- as.numeric(substr(rownames(test.sc$star)[1],3,4))
    
    p_vec <- test.sc$star[1,c("p.3","p.2","p.1")]
    t_min <- p_vec[which.min(p_vec)]
    if(which.min(p_vec)==2){
      tf <- "e"
      select_dt[,7][k] <- 0
    }else{
      tf <- "l"
      select_dt[,7][k] <- 1
    }
    
    select_dt[,3][k] <- ord+1
    select_dt[,4][k] <- d
    select_dt[,5][k] <- ifelse(as.character(test.sc$star[[5]][1]) != "",1,0)
    select_dt[,6][k] <- ifelse(as.character(test.sc$tvar[[5]]) != "",1,0)
    
  }
  
  select_ls[[r]] <- select_dt
  
  print(r)
  
}

drop_ls <- lapply(select_ls,function(x)x[,-1])

drop_ar <- array(as.numeric(unlist(drop_ls)),dim=c(c.num,6,P))

average_mat <- round(apply(drop_ar,1:2,mean),2)
average_dt <- data.table(c.names,average_mat)
colnames(average_dt) <- colnames(select_dt)
average_dt <- average_dt[order(c)]

# select commodities that present evidence of nonlinearity in each 
# and every expanding estimation window, wherein the autoregressive
# order is selected to minimize SIC subject to no serial correlation
sub_dt <- average_dt[n==1]
sub_dt <- sub_dt[order(c)]

cols <- sub_dt$c
cols <- cols[order(cols)]

subset_dt <- prices_dt[,..cols]
subset_dt$DATE <- prices_dt$DATE

save(subset_dt,sub_dt,file="subset.RData")
