library(ggplot2)
library(data.table)
library(tseries)
library(forecast)
library(lmtest)
library(urca)

rm(list=ls())
gc()

# this function selects the autoregressive order based on the SIC or the AIC rule
# SIC rule: 
# (1) find the order based on SIC; 
# (2) test for serial correlation in residuals;
# (3) increase the order by 1, if the null of no serial correlation is rejected;
# (4) repeat steps 2-3 until the autoregressive order equals the predetermined maximum order.
# AIC rule:
# (1) find the order based on AIC; 
# (2) test for serial correlation in residuals;
# (3) increase the order by 1, if the null of no serial correlation is rejected;
# (4) repeat steps 2-3 until the autoregressive order equals the predetermined maximum order;
# (5) test for serial correlation in residuals of a model with lower order;
# (6) decrease the order by 1, if the null of no serial correlation is failed to be rejected;
# (7) repeat steps 5-6 until the autoregressive order equals 1.
ar_order <- function(x,max_p=12,method="SIC"){
  d <- data.table(x=as.numeric(x))
  d <- d[,paste0("x",c(1:max_p)) := shift(y,n=1:max_p,fill=NA,type="lag")]
  d <- d[complete.cases(d)]
  ic_vec <- matrix(nrow=max_p,ncol=1)
  for(i in 1:max_p){
    fmla <- as.formula(paste("x~",paste0("x",c(1:i),collapse="+")))
    if(method=="SIC"){
      ic_vec[i,] <- log(crossprod(lm(fmla,data=d)$resid))+log(nrow(d))*(1+i)/nrow(d)
    }else{
      ic_vec[i,] <- log(crossprod(lm(fmla,data=d)$resid))+2*(1+i)/nrow(d)
    }
  }
  opt_p <- which.min(ic_vec)
  if(method=="AIC" & opt_p != 1){
    p <- opt_p
    for(j in (opt_p-1):1){
      fmla <- as.formula(paste("x~",paste0("x",c(1:j),collapse="+")))
      if(bgtest(lm(fmla,data=d))[4] < 0.05){
        p <- p
        break
      }else{
        p <- p-1
      }
    }
    opt_p <- p
  }
  p <- opt_p
  for(j in opt_p:max_p){
    fmla <- as.formula(paste("x~",paste0("x",c(1:j),collapse="+")))
    if(bgtest(lm(fmla,data=d))[4] >= 0.05){
      p <- p
      break
    }else{
      p <- p+1
    }
  }
  p <- ifelse(p<=max_p,p,max_p)
  return(p)
}

source("startest.r")

## load the data
load("dataset.RData")

c.num <- ncol(prices_dt)-1

p <- 12
h <- 12

c.names <- colnames(prices_dt)[-1]

select_dt <- data.table(c=c.names,i_sic=-99,p_sic=-99,d_sic=-99,n_sic=-99,t_sic=-99,i_aic=-99,p_aic=-99,d_aic=-99,n_aic=-99,t_aic=-99)

for(k in 1:c.num){
  
  cmd <- c.names[k]
  
  y <- as.matrix(log(prices_dt[,..cmd]))
  
  ord <- ar_order(y,max_p=p,method="SIC")
  
  adf.l <- ur.df(y,type="drift",lags=ord)
  
  adf.ts <- adf.l@teststat["statistic","tau2"]
  adf.cv <- adf.l@cval["tau2","5pct"]
  
  select_dt[,2][k] <- ifelse(adf.ts < adf.cv,0,1)
  
  # if(adf.ts < adf.cv){
  #   y <- y
  # }else{
    y <- c(0,diff(y))
    ord <- ord-1
  # }
  
  ymat <- embed(y,p+h+1)
  
  repeat{
    test.sc <- startest(y=ymat[,1],x=cbind(1,ymat[,2:(2+ord-1)]),x.n=cbind(1,ymat[,2:(2+ord-1)]),tvar=ymat[,2:(2+ord-1)],tvar.lag=0,contemp=F,ascending=T)
    if(test.sc$ac[1]>=.05 | ord==p-1){
      break
    }else{
      ord <- ord+1
    }
  }
  
  d_sic <- as.numeric(substr(rownames(test.sc$star)[1],3,4))
  
  select_dt[,3][k] <- ord+1#ifelse(adf.ts < adf.cv,ord,ord+1)
  select_dt[,4][k] <- d_sic
  select_dt[,5][k] <- ifelse(as.character(test.sc$star[[5]][1]) != "",1,0)
  select_dt[,6][k] <- ifelse(as.character(test.sc$tvar[[5]]) != "",1,0)
  
  
  y <- as.matrix(log(prices_dt[,..cmd]))
  
  ord <- ar_order(y,max_p=p,method="AIC")
  
  adf.l <- ur.df(y,type="drift",lags=ord)
  
  adf.ts <- adf.l@teststat["statistic","tau2"]
  adf.cv <- adf.l@cval["tau2","5pct"]
  
  select_dt[,7][k] <- ifelse(adf.ts < adf.cv,0,1)
  
  # if(adf.ts < adf.cv){
  #   y <- y
  # }else{
    y <- c(0,diff(y))
    ord <- ord-1
  # }
  
  ymat <- embed(y,p+h+1)
  
  repeat{
    test.sc <- startest(y=ymat[,1],x=cbind(1,ymat[,2:(2+ord-1)]),x.n=cbind(1,ymat[,2:(2+ord-1)]),tvar=ymat[,2:(2+ord-1)],tvar.lag=0,contemp=F,ascending=T)
    if(test.sc$ac[1]>=.05 | ord==p-1){
      break
    }else{
      ord <- ord+1
    }
  }
  
  d_sic <- as.numeric(substr(rownames(test.sc$star)[1],3,4))
  
  select_dt[,8][k] <- ord+1#ifelse(adf.ts < adf.cv,ord,ord+1)
  select_dt[,9][k] <- d_sic
  select_dt[,10][k] <- ifelse(as.character(test.sc$star[[5]][1]) != "",1,0)
  select_dt[,11][k] <- ifelse(as.character(test.sc$tvar[[5]]) != "",1,0)
  
}

# select commodities that present evidence of nonlinearity, 
# where the autoregressive order is first selected based on 
# the 'SIC rule' or the 'AIC rule'; omit commodities with
# auoregressive order of 10 or more.
sub_dt <- select_dt[(n_sic==1 & n_aic==1) & (p_sic<10 & p_aic<10)]

sub_dt$i <- sub_dt$i_sic
sub_dt$p <- sub_dt$p_sic
sub_dt$d <- sub_dt$d_sic

for(i in 1:nrow(sub_dt)){
  if(sub_dt$n_sic[i]!=1){
    sub_dt$i[i] <- sub_dt$i_aic[i]
    sub_dt$p[i] <- sub_dt$p_aic[i]
    sub_dt$d[i] <- sub_dt$d_aic[i]
  }
}
sub_dt <- sub_dt[order(sub_dt$c)]

cols <- sub_dt$c
cols <- cols[order(cols)]

subset_dt <- prices_dt[,..cols]
subset_dt$DATE <- prices_dt$DATE

save(subset_dt,sub_dt,file="subset.RData")
