library(ggplot2)
library(reshape2)
library(sandwich)
library(lmtest)
library(tseries)
library(dplyr)
library(urca)
library(forecast)
library(data.table)

rm(list=ls())
gc()

source("startest.r")
source("lstar.r")
source("estar.r")

## load the data
load("subset.RData")
load("select.RData")

c.num <- ncol(subset_dt)-1

p <- 12
h <- 12

B <- 5000

frac <- .85

# in-sample and out-of-sample periods
R <- round(frac*nrow(subset_dt))
P <- nrow(subset_dt)-R-h

# empty array to store: realised values, naive forecasts, iterated and direct
# forecasts from STAR models, iterated and direct forecasts from AR models
array.forecast <- array(dim=c(h,6,P,c.num))

c.names <- as.character(diagnostics_dt$c)

for(k in 1:c.num){
  
  cmd <- c.names[k]
  
  y <- as.matrix(log(subset_dt[,..cmd])) 
  
  y.lvl <- y 
  
  res <- array(dim=c(h,6,P))
  
  # if(diagnostics_dt[c==cmd]$i ==0){
  #   y.all <- y.lvl
  #   ord <- diagnostics_dt[c==cmd]$p
  # }else{
    y.all <- c(0,diff(y.lvl))
    ord <- diagnostics_dt[c==cmd]$p-1
  # }
  
  
  for(r in 1:P){
    
    y <- y.all[r:(R+r-1)] # the rolling subset for estimation
    
    y.realised <- y.lvl[(R+r):(R+r+h-1)] # the out-of-sample realized values
    y.last <- y.lvl[(R+r-1)]
    
    ## random walk forecast
    
    y.naive <- rep(y.last,h)
    
    #####################
    ## iterated method ##
    #####################
    
    ymat <- embed(y,p+h+1)
    
    y.f <- c(y,rep(NA,h))
    
    ymat.f <- embed(y.f,p+h+1)
    
    # p.vec[r,] <- ord
    
    ## LINEAR MODEL
    
    est.ar <- lm(ymat[,1]~cbind(ymat[,2:(2+ord-1)]))
    b <- est.ar$coef
    
    e <- est.ar$resid
    
    y.boot <- matrix(nrow=h,ncol=B)
    
    set.seed(k)
    e.boot <- matrix(sample(e,h*B,replace=T),nrow=h,ncol=B)
    
    for(itnum in 1:B){
      
      e.f <- e.boot[,itnum]
      
      for(i in 1:h){
        
        ymat.f[nrow(ymat)+i,2:(p+1)] <- ymat.f[nrow(ymat)+i-1,1:p]
        
        ymat.f[nrow(ymat)+i,1] <- b[1:(ord+1)]%*%c(1,ymat.f[nrow(ymat)+i,2:(ord+1)])+e.f[i]
        
      }
      
      y.boot[,itnum] <- ymat.f[(nrow(ymat)+1):(nrow(ymat)+h),1]
      
    }
    
    # if(diagnostics_dt[c==cmd]$i==1){
      y.boot <- apply(rbind(y.last,y.boot),2,cumsum)[-1,]
    # }
    
    ## point forecasts
    y.il_mean <- rowMeans(y.boot)
    
    
    ## NONLINEAR MODEL
    
    d <- star_lags[k,1]
    
    s.t <- ymat[,(d+1)]
    
    if(star_func[k,1]=="l"){
      est.star <- lstar(y=ymat[,1],x.0=cbind(1,ymat[,2:(2+ord-1)]),x.1=cbind(1,ymat[,2:(2+ord-1)]),tv=s.t)
    }else{
      est.star <- estar(y=ymat[,1],x.0=cbind(1,ymat[,2:(2+ord-1)]),x.1=cbind(1,ymat[,2:(2+ord-1)]),tv=s.t)
    }
    
    b <- est.star$coef
    e <- est.star$resid
    
    y.boot <- matrix(nrow=h,ncol=B)
    
    set.seed(k)
    e.boot <- matrix(sample(e,h*B,replace=T),nrow=h,ncol=B)
    
    for(itnum in 1:B){
      
      e.f <- e.boot[,itnum]
      
      for(i in 1:h){
        
        ymat.f[nrow(ymat)+i,2:(p+1)] <- ymat.f[nrow(ymat)+i-1,1:p]
        
        s.f <- ymat.f[nrow(ymat)+i-d,1]
        
        if(star_func[k,1]=="l"){
          g.f <- (1+exp(-b["g1"]*(s.f-b["c1"])/sd(s.t)))^(-1)
        }else{
          g.f <- (1-exp(-b["g1"]*((s.f-b["c1"])/sd(s.t))^2))
        }
        
        ymat.f[nrow(ymat)+i,1] <- b[1:(ord+1)]%*%c(1,ymat.f[nrow(ymat)+i,2:(ord+1)])+(b[(ord+1+1):((ord+1)*2)]%*%c(1,ymat.f[nrow(ymat)+i,2:(ord+1)]))*g.f+e.f[i]
        
      }
      
      y.boot[,itnum] <- ymat.f[(nrow(ymat)+1):(nrow(ymat)+h),1]
      
    }
    
    # if(diagnostics_dt[c==cmd]$i==1){
      y.boot <- apply(rbind(y.last,y.boot),2,cumsum)[-1,]
    # }
    
    ## point forecasts
    y.in_mean <- rowMeans(y.boot)
    
    
    ###################
    ## direct method ##
    ###################
    
    ymat <- embed(y,p+h+1)
    
    y.f <- c(y,rep(NA,h))
    
    ymat.f <- embed(y.f,p+h+1)
    
    y.boot <- matrix(nrow=h,ncol=B)
    y.dlin <- matrix(nrow=h,ncol=B)
    
    # if(diagnostics_dt[c==cmd]$i==0){
    #   
    #   for(j in 0:(h-1)){
    #     
    #     est.ar <- lm(ymat.f[,1]~cbind(ymat.f[,(j+2):(j+2+ord-1)]))
    #     b <- est.ar$coef
    #     e <- est.ar$resid
    #     
    #     set.seed(k)
    #     e.boot <- matrix(sample(e,h*B,replace=T),nrow=h,ncol=B)
    #     e.f <- e.boot[j+1,]
    #     
    #     # y.boot[j+1,] <- as.numeric(b[1:(ord+1)]%*%c(1,ymat.f[nrow(ymat),1:ord]))+as.matrix(e.f)
    #     y.dlin[j+1,] <- as.numeric(b[1:(ord+1)]%*%c(1,ymat.f[nrow(ymat),1:ord]))+as.matrix(e.f)
    #     
    #     dd <- star_lags[k,(j+1)]
    #     
    #     s.t <- ymat[,2+j+dd-1]
    #     s.f <- ymat.f[nrow(ymat)+1-dd,1]
    #     
    #     if(star_func[k,(j+1)]=="l"){
    #       est.star <- lstar(y=ymat[,1],x.0=cbind(1,ymat[,(2+j):(2+j+ord-1)]),x.1=cbind(1,ymat[,(2+j):(2+j+ord-1)]),tv=s.t)
    #       b <- est.star$coef
    #       g.f <- (1+exp(-b["g1"]*(s.f-b["c1"])/sd(s.t)))^(-1)
    #     }else{
    #       est.star <- estar(y=ymat[,1],x.0=cbind(1,ymat[,(2+j):(2+j+ord-1)]),x.1=cbind(1,ymat[,(2+j):(2+j+ord-1)]),tv=s.t)
    #       b <- est.star$coef
    #       g.f <- (1-exp(-b["g1"]*((s.f-b["c1"])/sd(s.t))^2))
    #     }
    #     
    #     e <- est.star$resid
    #     
    #     set.seed(k)
    #     e.boot <- matrix(sample(e,h*B,replace=T),nrow=h,ncol=B)
    #     e.f <- e.boot[j+1,]
    #     
    #     y.boot[j+1,] <- as.numeric(b[1:(ord+1)]%*%c(1,ymat.f[nrow(ymat),1:ord])+(b[(ord+1+1):((ord+1)*2)]%*%c(1,ymat.f[nrow(ymat),1:ord]))*g.f)+as.matrix(e.f)
    #     
    #   }
    #   
    # }else{
      
      for(j in 0:(h-1)){
        
        est.ar <- lm(rowSums(as.matrix(ymat.f[,1:(j+1)]))~cbind(ymat.f[,(j+2):(j+2+ord-1)]))
        b <- est.ar$coef
        e <- est.ar$resid
        
        set.seed(k)
        e.boot <- matrix(sample(e,h*B,replace=T),nrow=h,ncol=B)
        e.f <- e.boot[j+1,]
        
        # y.boot[j+1,] <- y.last+as.numeric(b[1:(ord+1)]%*%c(1,ymat.f[nrow(ymat),1:ord]))+as.matrix(e.f)
        y.dlin[j+1,] <- y.last+as.numeric(b[1:(ord+1)]%*%c(1,ymat.f[nrow(ymat),1:ord]))+as.matrix(e.f)
        
        dd <- star_lags[k,(j+1)]
        
        s.t <- ymat[,2+j+dd-1]
        s.f <- ymat.f[nrow(ymat)+1-dd,1]
        
        if(star_func[k,(j+1)]=="l"){
          est.star <- lstar(y=rowSums(as.matrix(ymat[,1:(j+1)])),x.0=cbind(1,ymat[,(2+j):(2+j+ord-1)]),x.1=cbind(1,ymat[,(2+j):(2+j+ord-1)]),tv=s.t)
          b <- est.star$coef
          g.f <- (1+exp(-b["g1"]*(s.f-b["c1"])/sd(s.t)))^(-1)
        }else{
          est.star <- estar(y=rowSums(as.matrix(ymat[,1:(j+1)])),x.0=cbind(1,ymat[,(2+j):(2+j+ord-1)]),x.1=cbind(1,ymat[,(2+j):(2+j+ord-1)]),tv=s.t)
          b <- est.star$coef
          g.f <- (1-exp(-b["g1"]*((s.f-b["c1"])/sd(s.t))^2))
        }
        
        e <- est.star$resid
        
        set.seed(k)
        e.boot <- matrix(sample(e,h*B,replace=T),nrow=h,ncol=B)
        e.f <- e.boot[j+1,]
        
        y.boot[j+1,] <- y.last+as.numeric(b[1:(ord+1)]%*%c(1,ymat.f[nrow(ymat),1:ord])+(b[(ord+1+1):((ord+1)*2)]%*%c(1,ymat.f[nrow(ymat),1:ord]))*g.f)+as.matrix(e.f)
        
      }
      
    # }
    
    ## point forecasts
    y.dl_mean <- rowMeans(y.dlin)
    y.dn_mean <- rowMeans(y.boot)
    
    # if(y.in_mean[1]!=y.dn_mean[1]){
    #   break
    # }
    
    # store the forecasts
    res[,,r] <- cbind(y.realised,y.naive,y.il_mean,y.in_mean,y.dl_mean,y.dn_mean)
    
  }
  
  array.forecast[,,,k] <- res
  print(k)
  
}


save(array.forecast,file="forecast85.RData")

