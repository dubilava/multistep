estar <- function(y,x.0,x.1=x.0,tv,crit=.25,g.opt=NULL,c.opt=NULL){
  
  require(minpack.lm)
  
  y   <- as.matrix(y)
  
  x.0 <- as.matrix(x.0)
  x.1 <- as.matrix(x.1)
  m.x0 <- ncol(x.0)
  m.x1 <- ncol(x.1)
  
  tv  <- as.matrix(tv)
  sig <- sd(tv)
  
  c.rng <- quantile(tv,c(crit,1-crit))
  
  if(is.null(g.opt)){
    g.opt <- 0.5
  }
  
  if(is.null(c.opt)){
    c.opt <- quantile(tv,.5)
  }
  
  glab <- "g1"
  clab <- "c1" 
  
  trans <- as.numeric((1-exp(-(g.opt/(sig^2))*((tv-c.opt)^2))))
  xb <- cbind(x.0,x.1*trans)
  
  res.mat <- matrix(nrow=3,ncol=3)
  
  g.vec <- c(.5,1,2)
  c.vec <- quantile(tv,c(.4,.5,.6))
  
  for(gg in 1:3){
    for(cc in 1:3){
      
      g.opt <- g.vec[gg]
      c.opt <- c.vec[cc]
      
      b <- as.numeric(c(lm(y~xb-1)$coef,g.opt,c.opt))
      set.seed(1981)
      b <- 0.9*b+rnorm(length(b),0,0.01)
      names(b) <- c(sprintf("a%d",1:m.x0),sprintf("b%d",1:m.x1),glab,clab)
      
      x0 <- paste0("x.0[,",1:m.x0,"]*","a",1:m.x0,collapse="+")
      x1 <- paste0("x.1[,",1:m.x1,"]*","b",1:m.x1,collapse="+")
      
      func <- paste0("((1-exp(-",glab,"/(sig^2)*((tv-(",clab,"))^2))))")
      fmla <- as.formula(paste0("y ~ ",x0,"+(",x1,")*",func))
      
      reg <- tryCatch(nlsLM(fmla,start=b,lower=c(rep(-Inf,length(b)-2),.2,c.rng[1]),upper=c(rep(Inf,length(b)-2),5,c.rng[2]),control=nls.control(maxiter=1000,warnOnly=T)),error=function(e) NULL)
      
      if(is.null(reg)){
        res.mat[gg,cc] <- 1e+32
      }else{
        res.mat[gg,cc] <- crossprod(resid(reg))
      }
      
    }
  }
  
  g.opt <- g.vec[which(res.mat==min(res.mat),arr.ind=TRUE)[1]]
  c.opt <- c.vec[which(res.mat==min(res.mat),arr.ind=TRUE)[2]]
  
  b <- as.numeric(c(lm(y~xb-1)$coef,g.opt,c.opt))
  set.seed(1981)
  b <- 0.9*b+rnorm(length(b),0,0.01)
  names(b) <- c(sprintf("a%d",1:m.x0),sprintf("b%d",1:m.x1),glab,clab)
  
  x0 <- paste0("x.0[,",1:m.x0,"]*","a",1:m.x0,collapse="+")
  x1 <- paste0("x.1[,",1:m.x1,"]*","b",1:m.x1,collapse="+")
  
  func <- paste0("((1-exp(-",glab,"/(sig^2)*((tv-(",clab,"))^2))))")
  fmla <- as.formula(paste0("y ~ ",x0,"+(",x1,")*",func))
  
  reg <- tryCatch(nlsLM(fmla,start=b,lower=c(rep(-Inf,length(b)-2),.2,c.rng[1]),upper=c(rep(Inf,length(b)-2),5,c.rng[2]),control=nls.control(maxiter=1000,warnOnly=T)),error=function(e) NULL)
  
  return(list(coef=coef(reg),se=summary(reg)$coefficients[,2],resid=resid(reg),vcov=vcov(reg)))
  
}
