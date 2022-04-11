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


arh_order <- function(x,h,max_lag=12,method="SIC",bg=TRUE){
  ic_func <- function(i){# optimal lag based on (A/S)IC
    fmla <- as.formula(paste("x~",paste0("x",c((j+1):(j+i)),collapse="+")))
    reg <- lm(fmla,data=d)
    if(method=="SIC"){
      log(crossprod(reg$resid))+log(nrow(d))*length(reg$coefficients)/nrow(d)
    }else{
      log(crossprod(reg$resid))+2*length(reg$coefficients)/nrow(d)
    }
  }
  bg_func <- function(i){# Breusch-Godfrey
    fmla <- as.formula(paste("x~",paste0("x",c((j+1):(j+i)),collapse="+")))
    as.numeric(bgtest(lm(fmla,data=d))[4])
  }
  j <- h-1
  d <- data.table(x=as.numeric(x))
  d <- d[,paste0("x",c((j+1):(j+max_lag))) := shift(x,n=(j+1):(j+max_lag),fill=NA,type="lag")]
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


# load the data
load("subset.RData")

c_num <- ncol(subset_dt)-1

p <- 12
h <- 12

B <- 5000

frac <- .80

# in-sample and out-of-sample periods
n <- nrow(subset_dt)
R <- round(frac*n)
P <- n-R-h+1

# empty array to store: realized values, naive forecasts, iterated and direct
# forecasts from STAR models, iterated and direct forecasts from AR models
array_forecast <- array(dim=c(h,6,P,c_num))
array_pval <- array(dim=c(h,P,c_num))
array_order <- array(dim=c(h,P,c_num))

c_names <- as.character(names(subset_dt)[1:c_num])

for(k in 1:c_num){
  
  cmd <- c_names[k]
  
  y <- as.matrix(log(subset_dt[,..cmd])) 
  
  res <- array(dim=c(h,6,P))
  pval <- matrix(nrow=h,ncol=P)
  lag_order <- matrix(nrow=h,ncol=P)
  
  dy <- c(0,diff(y))
  
  
  for(r in 1:P){
    
    # the expanding window for estimation
    w <- y[1:(R-1+r)]
    dw <- dy[1:(R-1+r)]
    
    ord_w <- ar_order(w,max_lag=p,method="SIC",bg=F)
    
    y_realised <- y[(R+r):(R+r+h-1)] # the out-of-sample realized values
    y_last <- y[(R-1+r)]
    
    # random walk forecast
    y_naive <- rep(y_last,h)
    
    # iterated method
    
    dmat <- embed(dw,p+h+1)
    lmat <- embed(w,p+h+1)
    
    dy_f <- c(dw,rep(NA,h))
    y_f <- c(w,rep(NA,h))
    
    dmat_f <- embed(dy_f,p+h+1)
    lmat_f <- embed(y_f,p+h+1)
    
    # linear iterated
    
    if(ord_w>1){
      est_ar <- lm(dmat[,1]~dmat[,2:ord_w])
    }else{
      est_ar <- lm(dmat[,1]~1)
    }
    
    b <- est_ar$coef
    e <- est_ar$resid
    
    y_boot <- matrix(nrow=h,ncol=B)
    
    set.seed(k)
    e_boot <- matrix(sample(e,h*B,replace=T),nrow=h,ncol=B)
    
    for(itnum in 1:B){
      
      e_f <- e_boot[,itnum]
      
      for(i in 1:h){
        
        if(ord_w>1){
          dmat_f[nrow(dmat)+i,2:ord_w] <- dmat_f[nrow(dmat)+i-1,1:(ord_w-1)]
          dmat_f[nrow(dmat)+i,1] <- b[1:ord_w]%*%c(1,dmat_f[nrow(dmat)+i,2:ord_w])+e_f[i]
        }else{
          dmat_f[nrow(dmat)+i,1] <- b+e_f[i]
        }
        
      }
      
      y_boot[,itnum] <- dmat_f[(nrow(dmat)+1):(nrow(dmat)+h),1]
      
    }
    
    y_boot <- apply(rbind(y_last,y_boot),2,cumsum)[-1,]

    y_lin_iter <- rowMeans(y_boot)
    
    
    # nonlinear iterated
    if(ord_w>1){
      test_nl <- startest(y=dmat[,1],x=cbind(1,dmat[,2:ord_w]),x.n=cbind(1,dmat[,2:ord_w]),tvar=dmat[,2:ord_w],tvar.lag=0,contemp=F,ascending=T)
    }else{
      test_nl <- startest(y=dmat[,1],x=rep(1,nrow(dmat)),x.n=rep(1,nrow(dmat)),tvar=dmat[,2],tvar.lag=0,contemp=F,ascending=T)
    }
    
    d <- as.numeric(substr(rownames(test_nl$star)[1],3,4))
    
    s_t <- dmat[,(1+d)]
    
    p_vec <- test_nl$star[1,c("p.3","p.2","p.1")]
    t_min <- p_vec[which.min(p_vec)]
    if(which.min(p_vec)==2){
      tf <- "e"
    }else{
      tf <- "l"
    }
    
    if(ord_w>1){
      if(tf=="l"){
        est_star <- lstar(y=dmat[,1],x.0=cbind(1,dmat[,2:ord_w]),x.1=cbind(1,dmat[,2:ord_w]),tv=s_t)
      }else{
        est_star <- estar(y=dmat[,1],x.0=cbind(1,dmat[,2:ord_w]),x.1=cbind(1,dmat[,2:ord_w]),tv=s_t)
      }
    }else{
      if(tf=="l"){
        est_star <- lstar(y=dmat[,1],x.0=rep(1,nrow(dmat)),x.1=rep(1,nrow(dmat)),tv=s_t)
      }else{
        est_star <- estar(y=dmat[,1],x.0=rep(1,nrow(dmat)),x.1=rep(1,nrow(dmat)),tv=s_t)
      }
    }
    
    
    b <- est_star$coef
    e <- est_star$resid
    
    y_boot <- matrix(nrow=h,ncol=B)
    
    set.seed(k)
    e_boot <- matrix(sample(e,h*B,replace=T),nrow=h,ncol=B)
    
    for(itnum in 1:B){
      
      e_f <- e_boot[,itnum]
      
      for(i in 1:h){
        
        if(ord_w>1){
          dmat_f[nrow(dmat)+i,2:ord_w] <- dmat_f[nrow(dmat)+i-1,1:(ord_w-1)]
          s_f <- dmat_f[nrow(dmat)+i-d,1]
          if(tf=="l"){
            g_f <- (1+exp(-b["g1"]*(s_f-b["c1"])/sd(s_t)))^(-1)
          }else{
            g_f <- (1-exp(-b["g1"]*((s_f-b["c1"])/sd(s_t))^2))
          }
          dmat_f[nrow(dmat)+i,1] <- b[1:ord_w]%*%c(1,dmat_f[nrow(dmat)+i,2:ord_w])+(b[(ord_w+1):(ord_w*2)]%*%c(1,dmat_f[nrow(dmat)+i,2:ord_w]))*g_f+e_f[i]
        }else{
          s_f <- dmat_f[nrow(dmat)+i-d,1]
          if(tf=="l"){
            g_f <- (1+exp(-b["g1"]*(s_f-b["c1"])/sd(s_t)))^(-1)
          }else{
            g_f <- (1-exp(-b["g1"]*((s_f-b["c1"])/sd(s_t))^2))
          }
          dmat_f[nrow(dmat)+i,1] <- b[1]+b[2]*g_f+e_f[i]
        }
        
      }
      
      y_boot[,itnum] <- dmat_f[(nrow(dmat)+1):(nrow(dmat)+h),1]
      
    }
    
    y_boot <- apply(rbind(y_last,y_boot),2,cumsum)[-1,]
    
    y_nln_iter <- rowMeans(y_boot)
    

    # direct method
    
    y_dlin <- matrix(nrow=h,ncol=B)
    y_dnln <- matrix(nrow=h,ncol=B)
    
    for(j in 0:(h-1)){
      
      dj <- c(rep(0,j+1),diff(w,j+1))
      
      dmat_j <- embed(dj,p+h+1)
      
      ord_w <- arh_order(w,j+1,max_lag=p,method="SIC",bg=F)
      
      # linear direct
      if(ord_w>1){
        est_ar <- lm(dmat_j[,1]~dmat[,(j+2):(j+ord_w)])
      }else{
        est_ar <- lm(dmat_j[,1]~1)
      }
      b <- est_ar$coef
      e <- est_ar$resid
      
      set.seed(k)
      e_boot <- matrix(sample(e,h*B,replace=T),nrow=h,ncol=B)
      e_f <- e_boot[j+1,]
      
      if(ord_w>1){
        y_dlin[j+1,] <- y_last+as.numeric(b%*%c(1,dmat_f[nrow(dmat),1:(ord_w-1)]))+as.matrix(e_f)
      }else{
        y_dlin[j+1,] <- y_last+as.numeric(b)+as.matrix(e_f)
      }
      
      # nonlinear direct
      if(ord_w>1){
        test_nl <- startest(y=dmat_j[,1],x=cbind(1,dmat[,(j+2):(j+ord_w)]),x.n=cbind(1,dmat[,(j+2):(j+ord_w)]),tvar=dmat[,(j+2):(j+ord_w)],tvar.lag=0,contemp=F,ascending=T)
      }else{
        test_nl <- startest(y=dmat_j[,1],x=rep(1,nrow(dmat)),x.n=rep(1,nrow(dmat)),tvar=dmat[,(j+2)],tvar.lag=0,contemp=F,ascending=T)
      }
      
      dd <- as.numeric(substr(rownames(test_nl$star)[1],3,4))
      
      s_t <- dmat[,(j+1+dd)]
      s_f <- dmat[nrow(dmat),dd]
      
      p_vec <- test_nl$star[1,c("p.3","p.2","p.1")]
      t_min <- p_vec[which.min(p_vec)]
      if(which.min(p_vec)==2){
        tf <- "e"
      }else{
        tf <- "l"
      }
      
      if(ord_w>1){
        if(tf=="l"){
          est_star <- lstar(y=dmat_j[,1],x.0=cbind(1,dmat[,(j+2):(j+ord_w)]),x.1=cbind(1,dmat[,(j+2):(j+ord_w)]),tv=s_t)
          b <- est_star$coef
          g_f <- (1+exp(-b["g1"]*(s_f-b["c1"])/sd(s_t)))^(-1)
        }else{
          est_star <- estar(y=dmat_j[,1],x.0=cbind(1,dmat[,(j+2):(j+ord_w)]),x.1=cbind(1,dmat[,(j+2):(j+ord_w)]),tv=s_t)
          b <- est_star$coef
          g_f <- (1-exp(-b["g1"]*((s_f-b["c1"])/sd(s_t))^2))
        }
      }else{
        if(tf=="l"){
          est_star <- lstar(y=dmat_j[,1],x.0=rep(1,nrow(dmat)),x.1=rep(1,nrow(dmat)),tv=s_t)
          b <- est_star$coef
          g_f <- (1+exp(-b["g1"]*(s_f-b["c1"])/sd(s_t)))^(-1)
        }else{
          est_star <- estar(y=dmat_j[,1],x.0=rep(1,nrow(dmat)),x.1=rep(1,nrow(dmat)),tv=s_t)
          b <- est_star$coef
          g_f <- (1-exp(-b["g1"]*((s_f-b["c1"])/sd(s_t))^2))
        }
      }
      
      e <- est_star$resid
      
      set.seed(k)
      e_boot <- matrix(sample(e,h*B,replace=T),nrow=h,ncol=B)
      e_f <- e_boot[j+1,]
      
      if(ord_w>1){
        y_dnln[j+1,] <- y_last+as.numeric(b[1:ord_w]%*%c(1,dmat_f[nrow(dmat),1:(ord_w-1)])+(b[(ord_w+1):(ord_w*2)]%*%c(1,dmat_f[nrow(dmat),1:(ord_w-1)]))*g_f)+as.matrix(e_f)
      }else{
        y_dnln[j+1,] <- y_last+as.numeric(b[1]+b[2]*g_f)+as.matrix(e_f)
      }
      
      pval[j+1,r] <- as.numeric(t_min)
      
      lag_order[j+1,r] <- ord_w
      
    }
    
    y_lin_dir <- rowMeans(y_dlin)
    y_nln_dir <- rowMeans(y_dnln)
    
    res[,,r] <- cbind(y_realised,y_naive,y_lin_iter,y_nln_iter,y_lin_dir,y_nln_dir)
    
  }
  
  array_forecast[,,,k] <- res
  array_pval[,,k] <- pval
  array_order[,,k] <- lag_order
  print(k)
  
}

save(array_forecast,array_pval,file="forecast80.RData")
