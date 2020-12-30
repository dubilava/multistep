library(ggplot2)
library(data.table)
library(tseries)
library(forecast)
library(lmtest)

rm(list=ls())
gc()

source("startest.r")

## load the data
load("dataset.RData")

c.num <- ncol(prices_dt)-1

p <- 12
h <- 12

frac <- 1

## in-sample and out-of-sample periods
R <- round(frac*nrow(prices_dt))
P <- nrow(prices_dt)-R-h

c.names <- colnames(prices_dt)[-1]

select_dt <- data.table(c=c.names,d=-99,p=-99,n=-99,t=-99,s=-99)

for(k in 1:c.num){
  
  cmd <- c.names[k]
  
  y <- as.matrix(log(prices_dt[,..cmd])) # select the time series

  y.lvl <- y # save it in levels
  
  y <- y.lvl[1:R]
  
  adf.l <- adf.test(y)
  
  if(adf.l[4] < .05){
    y <- y
    y.all <- y.lvl 
  }else{
    y <- c(0,diff(y))
    y.all <- c(0,diff(y.lvl)) 
  }
  
  ymat <- embed(y,p+h)
  
  dum <- seasonaldummy(ts(1:nrow(ymat),f=12))
  
  
  if(adf.l[4] < .05){ # if the series are in levels
    
    ic.vec <- matrix(nrow=p,ncol=2)
    
    for(i in 1:p){
      for(j in 1:2){
        if(j==1){
          ic.ij <- log(crossprod(lm(ymat[,1]~ymat[,2:(2+i-1)])$resid))+log(nrow(ymat))*(1+i)/nrow(ymat)
        }else{
          ic.ij <- log(crossprod(lm(ymat[,1]~ymat[,2:(2+i-1)]+dum)$resid))+log(nrow(ymat))*(1+i+11)/nrow(ymat)
        }
        ic.vec[i,j] <- ic.ij
      }
    }
    
    ord <- which.min(ic.vec)
    
    if(ord > p){
      ssn <- 1
      ord <- ord-p
    }else{
      ssn <- 0
      dum <- matrix(nrow=nrow(dum),ncol=0)
    }
    
    if(ssn==1){
      if(ord < p){
        for(i in ord:(p-1)){
          if(bgtest(lm(ymat[,1]~ymat[,2:(2+i-1)]+dum))[4] >= 0.05){
            ord <- ord
            break
          }else{
            ord <- ord+1
          }
        }
      }
    }else{
      if(ord < p){
        for(i in ord:(p-1)){
          if(bgtest(lm(ymat[,1]~ymat[,2:(2+i-1)]))[4] >= 0.05){
            ord <- ord
            break
          }else{
            ord <- ord+1
          }
        }
      }
    }
    
  }else{
    
    ic.vec <- matrix(nrow=(p-1),ncol=2)
    
    for(i in 1:(p-1)){
      for(j in 1:2){
        if(j==1){
          ic.ij <- log(crossprod(lm(ymat[,1]~ymat[,2:(2+i-1)])$resid))+log(nrow(ymat))*(1+i)/nrow(ymat)
        }else{
          ic.ij <- log(crossprod(lm(ymat[,1]~ymat[,2:(2+i-1)]+dum)$resid))+log(nrow(ymat))*(1+i+11)/nrow(ymat)
        }
        ic.vec[i,j] <- ic.ij
      }
    }
    
    ord <- which.min(ic.vec)
    if(ord>(p-1)){
      ssn <- 1
      ord <- ord-(p-1)
    }else{
      ssn <- 0
      dum <- matrix(nrow=nrow(dum),ncol=0)
    }
    
    if(ssn==1){
      if(ord < (p-1)){
        for(i in ord:(p-2)){
          if(bgtest(lm(ymat[,1]~ymat[,2:(2+i-1)]+dum))[4] >= 0.05){
            ord <- ord
            break
          }else{
            ord <- ord+1
          }
        }
      }
    }else{
      if(ord < (p-1)){
        for(i in ord:(p-2)){
          if(bgtest(lm(ymat[,1]~ymat[,2:(2+i-1)]))[4] >= 0.05){
            ord <- ord
            break
          }else{
            ord <- ord+1
          }
        }
      }
    }
    
  }
  
  test.sc <- startest(y=ymat[,1],x=cbind(1,ymat[,2:(2+ord-1)],dum),x.n=cbind(1,ymat[,2:(2+ord-1)]),tvar=ymat[,2:(2+ord-1)],tvar.lag=0,contemp=F,ascending=T)

  select_dt[,2][k] <- ifelse(adf.l[4] < .05,0,1)
  select_dt[,3][k] <- ifelse(adf.l[4] < .05,ord,ord+1)
  select_dt[,4][k] <- ifelse(as.character(test.sc$star[[5]][1]) != "",1,0)
  select_dt[,5][k] <- ifelse(as.character(test.sc$tvar[[5]]) != "",1,0)
  select_dt[,6][k] <- ifelse(ssn==1,1,0)
  
}

sub_dt <- select_dt[n==1 & t==0 & s==0 & p<=6]
sub_dt <- sub_dt[order(sub_dt$c)]

cols <- sub_dt$c
cols <- cols[order(cols)]

subset_dt <- prices_dt[,..cols]
subset_dt$DATE <- prices_dt$DATE

save(subset_dt,sub_dt,file="subset.RData")
