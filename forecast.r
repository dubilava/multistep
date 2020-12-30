library(ggplot2)
library(reshape2)
library(sandwich)
library(lmtest)
library(tseries)
library(hdrcde)
library(dplyr)
library(urca)
library(forecast)
library(data.table)

rm(list=ls())
gc()

my_hdr <- function(x){
  res <- hdr(x,prob=c(95))
  return(res$mode)
}

source("startest.r")
source("lstar.r")
source("estar.r")

## load the data
load("subset.RData")

c.num <- ncol(subset_dt)-1

p <- 12
h <- 12

B <- 2000

frac <- .75

## in-sample and out-of-sample periods
R <- round(frac*nrow(subset_dt))
P <- nrow(subset_dt)-R-h#-1+p


array.forecast <- array(dim=c(h,14,P,c.num))

# array.rmsfe <- array(dim=c(h,5,c.num))
# 
# array.ni <- array(dim=c(h,4,c.num))
# array.nd <- array(dim=c(h,4,c.num))
# array.id <- array(dim=c(h,4,c.num))
# array.il <- array(dim=c(h,4,c.num))
# array.dl <- array(dim=c(h,4,c.num))
# array.idl <- array(dim=c(h,4,c.num))
# 
# nl.i.vec <- matrix(nrow=c.num,ncol=1)
# nl.d.mat <- matrix(nrow=c.num,ncol=h)
# 
# array.nl.i <- array(dim=c(P,c.num))
# array.nl.d <- array(dim=c(P,h,c.num))
# 
# array.dd <- array(dim=c(P,h,c.num))

c.names <- as.character(sub_dt$c)

for(k in 1:c.num){
  
  cmd <- c.names[k]
  
  y <- as.matrix(log(subset_dt[,..cmd])) 

  y.lvl <- y 
  
  res <- array(dim=c(h,14,P))

  nl.i <- matrix(0,nrow=P,ncol=1)
  nl.d <- matrix(0,nrow=P,ncol=h)

  p.vec <- matrix(nrow=P,ncol=1)
  d.mat <- matrix(nrow=P,ncol=h)
  
  if(sub_dt[c==cmd]$d ==0){
    y.all <- y.lvl
  }else{
    y.all <- c(0,diff(y.lvl))
  }
  
  ord <- sub_dt[c==cmd]$p

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
    
    p.vec[r,] <- ord
    
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
    
    if(sub_dt[c==cmd]$d==1){
      y.boot <- apply(rbind(y.last,y.boot),2,cumsum)[-1,]
    }
    
    ## point forecasts
    y.il_mean <- rowMeans(y.boot)
    y.il_medi <- apply(y.boot,1,quantile,.5)
    y.il_mode <- apply(y.boot,1,my_hdr)
    
    
    ## NONLINEAR MODEL
    
    test.itr <- startest(y=ymat[,1],x=cbind(1,ymat[,2:(2+ord-1)]),x.n=cbind(1,ymat[,2:(2+ord-1)]),tvar=ymat[,2:(ord+1)],tvar.lag=0,contemp=F,ascending=T)
    
    d <- as.numeric(substr(rownames(test.itr$star)[1],3,4))
      
    s.t <- ymat[,(d+1)]
    
    if(as.character(test.itr$star[[1,5]])=="l"){
      est.star <- lstar(y=ymat[,1],x.0=cbind(1,ymat[,2:(2+ord-1)]),x.1=cbind(1,ymat[,2:(2+ord-1)]),tv=s.t)
    }else if(as.character(test.itr$star[[1,5]])=="e"){
      est.star <- estar(y=ymat[,1],x.0=cbind(1,ymat[,2:(2+ord-1)]),x.1=cbind(1,ymat[,2:(2+ord-1)]),tv=s.t)
    }else{# if(as.character(test.itr$star[[1,5]])=="b"){
      l.star <- tryCatch(lstar(y=ymat[,1],x.0=cbind(1,ymat[,2:(2+ord-1)]),x.1=cbind(1,ymat[,2:(2+ord-1)]),tv=s.t),error=function(e) NULL)
      e.star <- tryCatch(estar(y=ymat[,1],x.0=cbind(1,ymat[,2:(2+ord-1)]),x.1=cbind(1,ymat[,2:(2+ord-1)]),tv=s.t),error=function(e) NULL)
      if(is.null(l.star)){
        est.star <- e.star
        test.itr$star[[1,5]] <- "e"
      }else if(is.null(e.star)){
        est.star <- l.star
        test.itr$star[[1,5]] <- "l"
      }else if(crossprod(l.star$resid) <= crossprod(e.star$resid)){
        est.star <- l.star
        test.itr$star[[1,5]] <- "l"
      }else{
        est.star <- e.star
        test.itr$star[[1,5]] <- "e"
      }
    }
    
    b <- est.star$coef
    e <- est.star$resid
    
    if(as.character(test.itr$star[[1,5]])=="l"){
      G <- (1+exp(-b["g1"]*(s.t-b["c1"])/sd(s.t)))^(-1)
    }else{
      G <- (1-exp(-b["g1"]*((s.t-b["c1"])/sd(s.t))^2))
    }
    
    G_dt <- data.table(s.t,G)
    
    gg <- ggplot(G_dt,aes(x=s.t,y=G))+
      geom_point(size=1)+
      theme_classic()+
      theme(legend.title = element_blank())
    
    # print(gg)
    
    y.boot <- matrix(nrow=h,ncol=B)
    
    set.seed(k)
    e.boot <- matrix(sample(e,h*B,replace=T),nrow=h,ncol=B)
    
    for(itnum in 1:B){
      
      e.f <- e.boot[,itnum]
      
      for(i in 1:h){
        
        ymat.f[nrow(ymat)+i,2:(p+1)] <- ymat.f[nrow(ymat)+i-1,1:p]
        
        s.f <- ymat.f[nrow(ymat)+i-d,1]
        if(as.character(test.itr$star[[1,5]])=="l"){
          g.f <- (1+exp(-b["g1"]*(s.f-b["c1"])/sd(s.t)))^(-1)
        }else{
          g.f <- (1-exp(-b["g1"]*((s.f-b["c1"])/sd(s.t))^2))
        }
        
        ymat.f[nrow(ymat)+i,1] <- b[1:(ord+1)]%*%c(1,ymat.f[nrow(ymat)+i,2:(ord+1)])+(b[(ord+1+1):((ord+1)*2)]%*%c(1,ymat.f[nrow(ymat)+i,2:(ord+1)]))*g.f+e.f[i]
        
      }
      
      y.boot[,itnum] <- ymat.f[(nrow(ymat)+1):(nrow(ymat)+h),1]
      
    }
    
    if(sub_dt[c==cmd]$d==1){
      y.boot <- apply(rbind(y.last,y.boot),2,cumsum)[-1,]
    }
    
    ## point forecasts
    y.in_mean <- rowMeans(y.boot)
    y.in_medi <- apply(y.boot,1,quantile,.5)
    y.in_mode <- apply(y.boot,1,my_hdr)
    
    
    ###################
    ## direct method ##
    ###################
    
    ymat <- embed(y,p+h+1)
    
    y.f <- c(y,rep(NA,h))
    
    ymat.f <- embed(y.f,p+h+1)
    
    y.boot <- matrix(nrow=h,ncol=B)
    y.dlin <- matrix(nrow=h,ncol=B)
    
    if(sub_dt[c==cmd]$d==0){
      
      for(j in 0:(h-1)){
        
        est.ar <- lm(ymat.f[,1]~cbind(ymat.f[,(j+2):(j+2+ord-1)]))
        b <- est.ar$coef
        e <- est.ar$resid
        
        set.seed(k)
        e.boot <- matrix(sample(e,h*B,replace=T),nrow=h,ncol=B)
        e.f <- e.boot[j+1,]
        
        # y.boot[j+1,] <- as.numeric(b[1:(ord+1)]%*%c(1,ymat.f[nrow(ymat),1:ord]))+as.matrix(e.f)
        y.dlin[j+1,] <- as.numeric(b[1:(ord+1)]%*%c(1,ymat.f[nrow(ymat),1:ord]))+as.matrix(e.f)
        
        test.dir <- startest(y=ymat[,1],x=cbind(1,ymat[,(2+j):(2+j+ord-1)]),x.n=cbind(1,ymat[,(2+j):(2+j+ord-1)]),tvar=ymat[,(2+j):(2+j+ord-1)],tvar.lag=0,contemp=F,ascending=T)
        
        dd <- as.numeric(substr(rownames(test.dir$star)[1],3,4))

        s.t <- ymat[,2+j+dd-1]
        s.f <- ymat.f[nrow(ymat)+1-dd,1]
        
        if(as.character(test.dir$star[[1,5]])=="l"){
          est.star <- lstar(y=ymat[,1],x.0=cbind(1,ymat[,(2+j):(2+j+ord-1)]),x.1=cbind(1,ymat[,(2+j):(2+j+ord-1)]),tv=s.t)
          b <- est.star$coef
          g.f <- (1+exp(-b["g1"]*(s.f-b["c1"])/sd(s.t)))^(-1)
        }else if(as.character(test.dir$star[[1,5]])=="e"){
          est.star <- estar(y=ymat[,1],x.0=cbind(1,ymat[,(2+j):(2+j+ord-1)]),x.1=cbind(1,ymat[,(2+j):(2+j+ord-1)]),tv=s.t)
          b <- est.star$coef
          g.f <- (1-exp(-b["g1"]*((s.f-b["c1"])/sd(s.t))^2))
        }else{# if(as.character(test.dir$star[[1,5]])=="b"){
          l.star <- tryCatch(lstar(y=ymat[,1],x.0=cbind(1,ymat[,(2+j):(2+j+ord-1)]),x.1=cbind(1,ymat[,(2+j):(2+j+ord-1)]),tv=s.t),error=function(e) NULL)
          e.star <- tryCatch(estar(y=ymat[,1],x.0=cbind(1,ymat[,(2+j):(2+j+ord-1)]),x.1=cbind(1,ymat[,(2+j):(2+j+ord-1)]),tv=s.t),error=function(e) NULL)
          if(is.null(l.star)){
            est.star <- e.star
            b <- est.star$coef
            g.f <- (1-exp(-b["g1"]*((s.f-b["c1"])/sd(s.t))^2))
            test.dir$star[[1,5]] <- "e"
          }else if(is.null(e.star)){
            est.star <- l.star
            b <- est.star$coef
            g.f <- (1+exp(-b["g1"]*(s.f-b["c1"])/sd(s.t)))^(-1)
            test.dir$star[[1,5]] <- "l"
          }else if(crossprod(l.star$resid) <= crossprod(e.star$resid)){
            est.star <- l.star
            b <- est.star$coef
            g.f <- (1+exp(-b["g1"]*(s.f-b["c1"])/sd(s.t)))^(-1)
            test.dir$star[[1,5]] <- "l"
          }else{
            est.star <- e.star
            b <- est.star$coef
            g.f <- (1-exp(-b["g1"]*((s.f-b["c1"])/sd(s.t))^2))
            test.dir$star[[1,5]] <- "e"
          }
        }
        
        e <- est.star$resid
        
        set.seed(k)
        e.boot <- matrix(sample(e,h*B,replace=T),nrow=h,ncol=B)
        e.f <- e.boot[j+1,]
        
        y.boot[j+1,] <- as.numeric(b[1:(ord+1)]%*%c(1,ymat.f[nrow(ymat),1:ord])+(b[(ord+1+1):((ord+1)*2)]%*%c(1,ymat.f[nrow(ymat),1:ord]))*g.f)+as.matrix(e.f)
        
      }
      
    }else{
      
      for(j in 0:(h-1)){
        
        est.ar <- lm(rowSums(as.matrix(ymat.f[,1:(j+1)]))~cbind(ymat.f[,(j+2):(j+2+ord-1)]))
        b <- est.ar$coef
        e <- est.ar$resid
        
        set.seed(k)
        e.boot <- matrix(sample(e,h*B,replace=T),nrow=h,ncol=B)
        e.f <- e.boot[j+1,]
        
        y.boot[j+1,] <- y.last+as.numeric(b[1:(ord+1)]%*%c(1,ymat.f[nrow(ymat),1:ord]))+as.matrix(e.f)
        y.dlin[j+1,] <- y.last+as.numeric(b[1:(ord+1)]%*%c(1,ymat.f[nrow(ymat),1:ord]))+as.matrix(e.f)
        
        test.dir <- startest(y=rowSums(as.matrix(ymat[,1:(j+1)])),x=cbind(1,ymat[,(j+2):(j+2+ord-1)]),x.n=cbind(1,ymat[,(j+2):(j+2+ord-1)]),tvar=ymat[,(2+j):(2+j+ord-1)],tvar.lag=0,contemp=F,ascending=T)
        
        dd <- as.numeric(substr(rownames(test.dir$star)[1],3,4))
          
          s.t <- ymat[,2+j+dd-1]
          s.f <- ymat.f[nrow(ymat)+1-dd,1]
          
          if(as.character(test.dir$star[[1,5]])=="l"){
            est.star <- lstar(y=rowSums(as.matrix(ymat[,1:(j+1)])),x.0=cbind(1,ymat[,(2+j):(2+j+ord-1)]),x.1=cbind(1,ymat[,(2+j):(2+j+ord-1)]),tv=s.t)
            b <- est.star$coef
            g.f <- (1+exp(-b["g1"]*(s.f-b["c1"])/sd(s.t)))^(-1)
          }else if(as.character(test.dir$star[[1,5]])=="e"){
            est.star <- estar(y=rowSums(as.matrix(ymat[,1:(j+1)])),x.0=cbind(1,ymat[,(2+j):(2+j+ord-1)]),x.1=cbind(1,ymat[,(2+j):(2+j+ord-1)]),tv=s.t)
            b <- est.star$coef
            g.f <- (1-exp(-b["g1"]*((s.f-b["c1"])/sd(s.t))^2))
          }else{# if(as.character(test.dir$star[[1,5]])=="b"){
            tryCatch(l.star <- lstar(y=rowSums(as.matrix(ymat[,1:(j+1)])),x.0=cbind(1,ymat[,(2+j):(2+j+ord-1)]),x.1=cbind(1,ymat[,(2+j):(2+j+ord-1)]),tv=s.t),error=function(e) NULL)
            tryCatch(e.star <- estar(y=rowSums(as.matrix(ymat[,1:(j+1)])),x.0=cbind(1,ymat[,(2+j):(2+j+ord-1)]),x.1=cbind(1,ymat[,(2+j):(2+j+ord-1)]),tv=s.t),error=function(e) NULL)
            if(is.null(l.star)){
              est.star <- e.star
              b <- est.star$coef
              g.f <- (1-exp(-b["g1"]*((s.f-b["c1"])/sd(s.t))^2))
              test.dir$star[[1,5]] <- "e"
            }else if(is.null(e.star)){
              est.star <- l.star
              b <- est.star$coef
              g.f <- (1+exp(-b["g1"]*(s.f-b["c1"])/sd(s.t)))^(-1)
              test.dir$star[[1,5]] <- "l"
            }else if(crossprod(l.star$resid) <= crossprod(e.star$resid)){
              est.star <- l.star
              b <- est.star$coef
              g.f <- (1+exp(-b["g1"]*(s.f-b["c1"])/sd(s.t)))^(-1)
              test.dir$star[[1,5]] <- "l"
            }else{
              est.star <- e.star
              b <- est.star$coef
              g.f <- (1-exp(-b["g1"]*((s.f-b["c1"])/sd(s.t))^2))
              test.dir$star[[1,5]] <- "e"
            }
          }
          
          e <- est.star$resid
          
          set.seed(k)
          e.boot <- matrix(sample(e,h*B,replace=T),nrow=h,ncol=B)
          e.f <- e.boot[j+1,]
          
          y.boot[j+1,] <- y.last+as.numeric(b[1:(ord+1)]%*%c(1,ymat.f[nrow(ymat),1:ord])+(b[(ord+1+1):((ord+1)*2)]%*%c(1,ymat.f[nrow(ymat),1:ord]))*g.f)+as.matrix(e.f)
        
      }
      
    }

    ## point forecasts
    y.dl_mean <- rowMeans(y.dlin)
    y.dn_mean <- rowMeans(y.boot)
    y.dl_medi <- apply(y.dlin,1,quantile,.5)
    y.dn_medi <- apply(y.boot,1,quantile,.5)
    y.dl_mode <- apply(y.dlin,1,my_hdr)
    y.dn_mode <- apply(y.boot,1,my_hdr)
    
    # if(y.in_mean[1]!=y.dn_mean[1]){
    #   break
    # }
    
    # store the forecasts
    res[,,r] <- cbind(y.realised,y.naive,y.il_mean,y.in_mean,y.dl_mean,y.dn_mean,y.il_medi,y.in_medi,y.dl_medi,y.dn_medi,y.il_mode,y.in_mode,y.dl_mode,y.dn_mode)
    
  }
  
  array.forecast[,,,k] <- res
  print(k)
  
}




save(array.forecast,file="forecast.RData")

