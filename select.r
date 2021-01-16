library(ggplot2)
library(reshape2)
library(sandwich)
library(lmtest)
library(tseries)
library(dplyr)
library(urca)
library(moments)
library(forecast)
library(data.table)

rm(list=ls())
gc()

source("startest.r")
source("lstar.r")
source("estar.r")

## load the data
load("subset.RData")
# sub_dt[c %in% c("LGST","SWST")]$p <- sub_dt[c %in% c("LGST","SWST")]$p+2
# sub_dt[c %in% c("SRMP")]$p <- sub_dt[c %in% c("SRMP")]$p+1

c.num <- ncol(subset_dt)-1

p <- 12
h <- 12

c.names <- as.character(sub_dt$c)

diagnose_dt <- data.table(c=c.names,sig=-99,eta=-99,kap=-99,JB=-99,ET=-99,Gf=-99,g.est=-99,g.se=-99,c.est=-99,c.se=-99)

star_func <- matrix(nrow=length(c.names),ncol=h)
star_lags <- matrix(nrow=length(c.names),ncol=h)
star_pval <- matrix(nrow=length(c.names),ncol=h)

star_ls <- list()

for(k in 1:c.num){
  
  cmd <- c.names[k]
  
  y <- as.matrix(log(subset_dt[,..cmd])) 
  
  y.lvl <- y 
  
  # if(sub_dt[c==cmd]$i ==0){
  #   y.all <- y.lvl
  #   ord <- sub_dt[c==cmd]$p
  # }else{
    y.all <- c(0,diff(y.lvl))
    ord <- sub_dt[c==cmd]$p-1
  # }
  
  y <- y.all
  
  ymat <- embed(y,p+h+1)
  
  ## NONLINEAR MODEL
  
  # ord <- ord+1
  
  test.itr <- startest(y=ymat[,1],x=cbind(1,ymat[,2:(2+ord-1)]),x.n=cbind(1,ymat[,2:(2+ord-1)]),tvar=ymat[,2:(ord+1)],tvar.lag=0,contemp=F,ascending=T)
  
  test.itr
  
  d <- as.numeric(substr(rownames(test.itr$star)[1],3,4))
  star_lags[k,1] <- d
  
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
  
  diagnose_dt[k,2] <- round(sd(e),3)
  diagnose_dt[k,3] <- round(skewness(e),3)
  diagnose_dt[k,4] <- round(kurtosis(e),3)-3
  diagnose_dt[k,5] <- round(jarque.bera.test(e)$p.value,3)
  diagnose_dt[k,6] <- round(test.itr$ac[1],3)
  diagnose_dt[k,7] <- ifelse(test.itr$star[[1,5]]=="l",1,2)
  diagnose_dt[k,8] <- round(est.star$coef["g1"],3)
  diagnose_dt[k,9] <- round(est.star$se["g1"],3)
  diagnose_dt[k,10] <- round(est.star$coef["c1"],3)
  diagnose_dt[k,11] <- round(est.star$se["c1"],3)
  
  star_pval[k,1] <- as.numeric(test.itr$star[[1,1]])
  
  if(as.character(test.itr$star[[1,5]])=="l"){
    star_func[k,1] <- "l"
    G <- (1+exp(-b["g1"]*(s.t-b["c1"])/sd(s.t)))^(-1)
  }else{
    star_func[k,1] <- "e"
    G <- (1-exp(-b["g1"]*((s.t-b["c1"])/sd(s.t))^2))
  }
  
  G_dt <- data.table(Commodity=c.names[k],s.t,G)
  
  gg <- ggplot(G_dt,aes(x=s.t,y=G))+
    geom_point(size=1)+
    theme_classic()+
    theme(legend.title = element_blank())
  
  print(gg)
  
  star_ls[[k]] <- G_dt
  
  
  ###################
  ## direct method ##
  ###################
  
  ymat <- embed(y,p+h+1)
  
  # if(sub_dt[c==cmd]$i==0){
  #   
  #   for(j in 0:(h-1)){
  #     
  #     
  #     test.dir <- startest(y=ymat[,1],x=cbind(1,ymat[,(2+j):(2+j+ord-1)]),x.n=cbind(1,ymat[,(2+j):(2+j+ord-1)]),tvar=ymat[,(2+j):(2+j+ord-1)],tvar.lag=0,contemp=F,ascending=T)
  #     
  #     dd <- as.numeric(substr(rownames(test.dir$star)[1],3,4))
  #     star_lags[k,(j+1)] <- dd
  #     
  #     s.t <- ymat[,2+j+dd-1]
  #     
  #     if(as.character(test.dir$star[[1,5]])=="l"){
  #       est.star <- lstar(y=ymat[,1],x.0=cbind(1,ymat[,(2+j):(2+j+ord-1)]),x.1=cbind(1,ymat[,(2+j):(2+j+ord-1)]),tv=s.t)
  #     }else if(as.character(test.dir$star[[1,5]])=="e"){
  #       est.star <- estar(y=ymat[,1],x.0=cbind(1,ymat[,(2+j):(2+j+ord-1)]),x.1=cbind(1,ymat[,(2+j):(2+j+ord-1)]),tv=s.t)
  #     }else{# if(as.character(test.dir$star[[1,5]])=="b"){
  #       l.star <- tryCatch(lstar(y=ymat[,1],x.0=cbind(1,ymat[,(2+j):(2+j+ord-1)]),x.1=cbind(1,ymat[,(2+j):(2+j+ord-1)]),tv=s.t),error=function(e) NULL)
  #       e.star <- tryCatch(estar(y=ymat[,1],x.0=cbind(1,ymat[,(2+j):(2+j+ord-1)]),x.1=cbind(1,ymat[,(2+j):(2+j+ord-1)]),tv=s.t),error=function(e) NULL)
  #       if(is.null(l.star)){
  #         est.star <- e.star
  #         test.dir$star[[1,5]] <- "e"
  #       }else if(is.null(e.star)){
  #         est.star <- l.star
  #         test.dir$star[[1,5]] <- "l"
  #       }else if(crossprod(l.star$resid) <= crossprod(e.star$resid)){
  #         est.star <- l.star
  #         test.dir$star[[1,5]] <- "l"
  #       }else{
  #         est.star <- e.star
  #         test.dir$star[[1,5]] <- "e"
  #       }
  #     }
  #     
  #     star_pval[k,(j+1)] <- as.numeric(test.dir$star[[1,1]])
  #     
  #     if(as.character(test.dir$star[[1,5]])=="l"){
  #       star_func[k,(j+1)] <- "l"
  #     }else{
  #       star_func[k,(j+1)] <- "e"
  #     }
  #     
  #   }
  #   
  # }else{
    
    for(j in 0:(h-1)){

      test.dir <- startest(y=rowSums(as.matrix(ymat[,1:(j+1)])),x=cbind(1,ymat[,(j+2):(j+2+ord-1)]),x.n=cbind(1,ymat[,(j+2):(j+2+ord-1)]),tvar=ymat[,(2+j):(2+j+ord-1)],tvar.lag=0,contemp=F,ascending=T)
      
      dd <- as.numeric(substr(rownames(test.dir$star)[1],3,4))
      star_lags[k,(j+1)] <- dd
      
      s.t <- ymat[,2+j+dd-1]
      
      if(as.character(test.dir$star[[1,5]])=="l"){
        est.star <- lstar(y=rowSums(as.matrix(ymat[,1:(j+1)])),x.0=cbind(1,ymat[,(2+j):(2+j+ord-1)]),x.1=cbind(1,ymat[,(2+j):(2+j+ord-1)]),tv=s.t)
      }else if(as.character(test.dir$star[[1,5]])=="e"){
        est.star <- estar(y=rowSums(as.matrix(ymat[,1:(j+1)])),x.0=cbind(1,ymat[,(2+j):(2+j+ord-1)]),x.1=cbind(1,ymat[,(2+j):(2+j+ord-1)]),tv=s.t)
      }else{# if(as.character(test.dir$star[[1,5]])=="b"){
        tryCatch(l.star <- lstar(y=rowSums(as.matrix(ymat[,1:(j+1)])),x.0=cbind(1,ymat[,(2+j):(2+j+ord-1)]),x.1=cbind(1,ymat[,(2+j):(2+j+ord-1)]),tv=s.t),error=function(e) NULL)
        tryCatch(e.star <- estar(y=rowSums(as.matrix(ymat[,1:(j+1)])),x.0=cbind(1,ymat[,(2+j):(2+j+ord-1)]),x.1=cbind(1,ymat[,(2+j):(2+j+ord-1)]),tv=s.t),error=function(e) NULL)
        if(is.null(l.star)){
          est.star <- e.star
          test.dir$star[[1,5]] <- "e"
        }else if(is.null(e.star)){
          est.star <- l.star
          test.dir$star[[1,5]] <- "l"
        }else if(crossprod(l.star$resid) <= crossprod(e.star$resid)){
          est.star <- l.star
          test.dir$star[[1,5]] <- "l"
        }else{
          est.star <- e.star
          test.dir$star[[1,5]] <- "e"
        }
      }
      
      star_pval[k,(j+1)] <- as.numeric(test.dir$star[[1,1]])
      
      if(as.character(test.dir$star[[1,5]])=="l"){
        star_func[k,(j+1)] <- "l"
      }else{
        star_func[k,(j+1)] <- "e"
      }
      
    }
    
  # }
  
  print(k)
  
}

diagnostics_dt <- merge(sub_dt[,.(c,i,p,d)],diagnose_dt,by="c")

save(star_func,star_lags,star_pval,star_ls,diagnostics_dt,file="select.RData")

