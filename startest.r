startest <- function(y,x,x.n=NULL,tvar,tvar.lag=6,ascending=F,diag.lag=6,crit=.05,contemp=T){
  
  y <- as.matrix(y)
  n.y <- nrow(y)
  
  x <- as.matrix(x)
  n.l <- nrow(x)
  m.l <- ncol(x)
  
  if(is.null(x.n)){
    x.n <- x
  }
  x.n <- as.matrix(x.n) 
  n.n <- nrow(x.n)
  m.n <- ncol(x.n)
  
  tv     <- as.matrix(tvar)
  tvar.num <- ncol(tv)
  if(tvar.num==1){
    tv.mat <- as.matrix(embed(tv,tvar.lag+1))
  }else{
    tv.mat <- tv
  }
  n.t    <- nrow(tv.mat)
  
  if(contemp==T){
    if(tvar.lag==0 & tvar.num==1){
      colnames(tv.mat)   <- "s.0"
    }else if(tvar.lag!=0 & tvar.num==1){
      colnames(tv.mat)   <- c(paste0("s.0"),colnames(tv.mat[,-1],do.NULL=F,prefix="s."))
    }else{
      colnames(tv.mat)   <- colnames(tv.mat,do.NULL=F,prefix="s.")
    }
  }else{
    if(tvar.lag==0 & tvar.num==1){
      colnames(tv.mat)   <- "s.1"
    }else{
      colnames(tv.mat)   <- colnames(tv.mat,do.NULL=F,prefix="s.")
    }
  }
  
  n <- min(c(n.y,n.l,n.n,n.t))
  
  if(n.y > n){y      <- as.matrix(y[-c(1:(n.y-n)),])}
  if(n.l > n){x      <- as.matrix(x[-c(1:(n.l-n)),])}
  if(n.n > n){x.n    <- as.matrix(x.n[-c(1:(n.n-n)),])}
  if(n.t > n){
    tv.mat <- as.matrix(tv.mat[-c(1:(n.t-n)),])
    if(tvar.lag==0){
      colnames(tv.mat)   <- "s.0"  
    }
  }
  
  e <- as.matrix(lm(y~x-1)$resid)
  
  e.mat <- matrix(0,n,diag.lag)
  for(i in 1:diag.lag){
    e.mat[(i+1):n,i] <- e[1:(n-i),]
  }
  e2.mat <- e.mat^2
  ac.vec <- matrix(NA,diag.lag,1)
  arch.vec <- matrix(NA,diag.lag,1)
  for(i in 1:diag.lag){
    x.e          <- cbind(e.mat[,1:i],x) 
    ac.vec[i,]   <- round(anova(lm(e~x.e-1),lm(e~x-1))[2,6],4)
    e2           <- as.matrix(e^2)
    x.e2         <- cbind(1,e2.mat[,1:i])
    x.r2         <- matrix(1,n,1)
    arch.vec[i,] <- round(anova(lm(e2~x.e2-1),lm(e2~x.r2-1))[2,6],4)
  }
  rownames(ac.vec)   <- rownames(ac.vec,do.NULL=F,prefix="AC.")
  rownames(arch.vec) <- rownames(arch.vec,do.NULL=F,prefix="ARCH.")
  colnames(ac.vec)   <- 'AC'
  colnames(arch.vec) <- 'ARCH'
  
  
  nl.out <- matrix(NA,ncol(tv.mat),4)
  
  for(i in 1:ncol(tv.mat)){
    
    s.1 <- as.matrix(tv.mat[,i])
    s.2 <- s.1^2
    s.3 <- s.1^3
    
    t.1 <- as.matrix(c(1:n)/n)
    t.2 <- t.1^2
    t.3 <- t.1^3
    
    zs.1 <- as.numeric(s.1)*x.n
    zs.2 <- as.numeric(s.2)*x.n
    zs.3 <- as.numeric(s.3)*x.n
    
    zt.1 <- as.numeric(t.1)*x.n
    zt.2 <- as.numeric(t.2)*x.n
    zt.3 <- as.numeric(t.3)*x.n
    
    x.0  <- cbind(x)
    
    xs.1 <- cbind(x,zs.1)
    xs.2 <- cbind(x,zs.1,zs.2)
    xs.3 <- cbind(x,zs.1,zs.2,zs.3)
    
    xt.1 <- cbind(x,zt.1)
    xt.2 <- cbind(x,zt.1,zt.2)
    xt.3 <- cbind(x,zt.1,zt.2,zt.3)
    
    xs.1 <- xs.1[,!duplicated(t(xs.1))]
    xs.2 <- xs.2[,!duplicated(t(xs.2))]
    xs.3 <- xs.3[,!duplicated(t(xs.3))]
    
    xt.1 <- xt.1[,!duplicated(t(xt.1))]
    xt.2 <- xt.2[,!duplicated(t(xt.2))]
    xt.3 <- xt.3[,!duplicated(t(xt.3))]
    
    reg.0 <- lm(y~x.0-1)
    
    regs.1 <- lm(y~xs.1-1)
    regs.2 <- lm(y~xs.2-1)
    regs.3 <- lm(y~xs.3-1)
    
    regt.1 <- lm(y~xt.1-1)
    regt.2 <- lm(y~xt.2-1)
    regt.3 <- lm(y~xt.3-1)
    
    ps.n <- round(anova(regs.3,reg.0)[2,6],4)
    ps.3 <- round(anova(regs.3,regs.2)[2,6],4)
    ps.2 <- round(anova(regs.2,regs.1)[2,6],4)
    ps.1 <- round(anova(regs.1,reg.0)[2,6],4)
    
    pt.n <- round(anova(regt.3,reg.0)[2,6],4)
    pt.3 <- round(anova(regt.3,regt.2)[2,6],4)
    pt.2 <- round(anova(regt.2,regt.1)[2,6],4)
    pt.1 <- round(anova(regt.1,reg.0)[2,6],4)
    
    nl.out[i,] <- c(ps.n,ps.3,ps.2,ps.1)
    sc.out <- matrix(c(pt.n,pt.3,pt.2,pt.1),1,4)
    
  }
  
  rownames(nl.out) <- colnames(tv.mat)
  rownames(sc.out) <- 't*'
  
  # if(contemp == FALSE){
  #   nl.out <- nl.out[-1,]
  # }
  
  if(ascending == TRUE){
    if(ncol(tv.mat)!=0){
      nl.out <- nl.out[order(nl.out[,1]),]
    }
  }
  
  nl.rows <- min(ncol(tv.mat),nrow(nl.out))
  
  if(nl.rows==1){
    nl.out <- matrix(nl.out,1,4)
    rownames(nl.out) <- colnames(tv.mat)
  }
  
  nl.lab <- matrix(NA,nl.rows,1)
  for(k in 1:nl.rows){
    if(nl.out[k,1] >= crit){
      nl.lab[k,] <- paste('')
    }else if(((nl.out[k,2] < crit) | (nl.out[k,4] < crit)) & (nl.out[k,3] < crit)){
      nl.lab[k,] <- paste('b')
    }else if(((nl.out[k,2] == min(nl.out[k,2:4])) & (nl.out[k,2] < crit)) | ((nl.out[k,4] == min(nl.out[k,2:4])) & (nl.out[k,4] < crit))){
      nl.lab[k,] <- paste('l')
    }else if((nl.out[k,3] == min(nl.out[k,2:4])) & (nl.out[k,3] < crit)){
      nl.lab[k,] <- paste('e')
    }else{
      nl.lab[k,] <- paste('')
    }
  }


  if(sc.out[,1] >= crit){
    sc.lab <- paste('')
  }else if(((sc.out[,2] == min(sc.out[,2:4])) & (sc.out[,2] < crit)) | ((sc.out[,4] == min(sc.out[,2:4])) & (sc.out[,4] < crit))){
    sc.lab <- paste('l')
  }else if((sc.out[,3] == min(sc.out[,2:4])) & (sc.out[,3] < crit)){
    sc.lab <- paste('e')
  }else{
    sc.lab <- paste('')
  }

  
  nl.out <- formatC(nl.out,format="e",digits=3)
  sc.out <- formatC(sc.out,format="e",digits=3)
  
  nl.mat <- as.data.frame(cbind(nl.out,nl.lab),stringsAsFactors=F)
  sc.mat <- as.data.frame(cbind(sc.out,sc.lab),stringsAsFactors=F)
  
  colnames(nl.mat) <- c('p.0','p.3','p.2','p.1','MODEL')
  colnames(sc.mat) <- c('p.0','p.3','p.2','p.1','MODEL')
  
  res.list <- return(list(star=nl.mat,tvar=sc.mat,ac=ac.vec,arch=arch.vec))
}
