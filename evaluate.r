library(ggplot2)
library(reshape2)
library(sandwich)
library(lmtest)
library(tseries)
library(dplyr)
library(urca)
library(forecast)
library(data.table)
library(stargazer)
library(xtable)

rm(list=ls())
gc()

w_frac <- 80

## load the data
load("forecast80.RData")
load("subset.RData")

p <- 12
h <- 12

frac <- w_frac/100

# in-sample and out-of-sample periods
n <- nrow(subset_dt)
R <- round(frac*n)
P <- n-R-h+1

subset_dt$DATE[c(h+p+1,R)]

# P <- dim(array_forecast)[3]

c_num <- ncol(subset_dt)-1

c_names <- names(subset_dt)[1:c_num]
  
commodities <- c("Aluminum","Bananas","Cotton","Diammonium Phosphate","Fishmeal",
                 "Groundnut Oil","Groundnuts","Gold","Hides","Lamb",
                 "Lead","Soft Logs","Oranges","Olive Oil","Platinum",
                 "Plywood","Palm Oil","Rice","Rapeseed Oil","Rubber",
                 "Soybean Oil","Hard Sawnwood","Soft Sawnwood","Uran","Urea")


multistep_dt <- data.table(commodities,t(round(apply(array_pval,c(1,3),mean),3)))
colnames(multistep_dt)[-1] <- paste0("h",c(1:12))


rmsfe_mean <- array(dim=c(h,8,c_num))

mean.nl_i <- array(dim=c(h,4,c_num))
mean.nl_d <- array(dim=c(h,4,c_num))
mean.id_l <- array(dim=c(h,4,c_num))
mean.id_n <- array(dim=c(h,4,c_num))
mean.indl <- array(dim=c(h,4,c_num))

mean_i <- array(dim=c(h,4,c_num))
mean_n <- array(dim=c(h,4,c_num))

mean.nl_i_a <- array(dim=c(h,4,c_num))
mean.nl_d_a <- array(dim=c(h,4,c_num))
mean.id_l_a <- array(dim=c(h,4,c_num))
mean.id_n_a <- array(dim=c(h,4,c_num))
mean.indl_a <- array(dim=c(h,4,c_num))

mean_i_a <- array(dim=c(h,4,c_num))
mean_n_a <- array(dim=c(h,4,c_num))

for(k in 1:c_num){
  
  # forecast errors
  e.il_mean <- array_forecast[,1,,k]-array_forecast[,3,,k]
  e.in_mean <- array_forecast[,1,,k]-array_forecast[,4,,k]
  e.dl_mean <- array_forecast[,1,,k]-array_forecast[,5,,k]
  e.dn_mean <- array_forecast[,1,,k]-array_forecast[,6,,k]
  
  # combined forecast errors
  e.i_mean <- .5*e.il_mean+.5*e.in_mean
  e.d_mean <- .5*e.dl_mean+.5*e.dn_mean
  e.l_mean <- .5*e.il_mean+.5*e.dl_mean
  e.n_mean <- .5*e.in_mean+.5*e.dn_mean

  # root mean square forecast errors
  rmsfe.il_mean <- apply(e.il_mean,1,function(x) sqrt(mean(x^2)))
  rmsfe.in_mean <- apply(e.in_mean,1,function(x) sqrt(mean(x^2)))
  rmsfe.dl_mean <- apply(e.dl_mean,1,function(x) sqrt(mean(x^2)))
  rmsfe.dn_mean <- apply(e.dn_mean,1,function(x) sqrt(mean(x^2)))
  rmsfe.i_mean <- apply(e.i_mean,1,function(x) sqrt(mean(x^2)))
  rmsfe.d_mean <- apply(e.d_mean,1,function(x) sqrt(mean(x^2)))
  rmsfe.l_mean <- apply(e.l_mean,1,function(x) sqrt(mean(x^2)))
  rmsfe.n_mean <- apply(e.n_mean,1,function(x) sqrt(mean(x^2)))

  rmsfe_mean[,,k] <- cbind(rmsfe.il_mean,rmsfe.in_mean,rmsfe.dl_mean,rmsfe.dn_mean,rmsfe.i_mean,rmsfe.d_mean,rmsfe.l_mean,rmsfe.n_mean)

  # (1) iterated STAR vs iterated AR
  d.nl_i <- e.in_mean^2-e.il_mean^2
  # (2) direct STAR vs direct AR
  d.nl_d <- e.dn_mean^2-e.dl_mean^2
  # (3) iterated AR vs direct AR
  d.id_l <- e.il_mean^2-e.dl_mean^2
  # (4) iterated STAR vs direct STAR
  d.id_n <- e.in_mean^2-e.dn_mean^2
  # (5) iterated STAR vs direct AR
  d.indl <- e.in_mean^2-e.dl_mean^2
  
  # (C1) iterated combined vs direct combined
  d_i <- e.i_mean^2-e.d_mean^2
  # (C2) nonlinear combined vs linear combined
  d_n <- e.n_mean^2-e.l_mean^2
  
  
  # (1) iterated STAR vs iterated AR
  mat.nl_i <- matrix(nrow=h,ncol=4)
  mat.nl_i_a <- matrix(nrow=h,ncol=4)
  for(j in 1:h){
    d.reg <- lm(d.nl_i[j,]~1)
    mat.nl_i[j,] <- coeftest(d.reg,vcov.=NeweyWest(d.reg,lag=h))
    
    harvey <- sqrt((P+1-2*j+j*(j-1)/P)/P)
    mat.nl_i[j,3] <- mat.nl_i[j,3]*harvey
    mat.nl_i[j,4] <- 2*pt(abs(mat.nl_i[j,3]),df=P-1,lower.tail=F)
    
    if(j>=2){
      wd.nl_i <- colMeans(d.nl_i[1:j,])
      wd.reg <- lm(wd.nl_i~1)
      mat.nl_i_a[j,] <- coeftest(wd.reg,vcov.=NeweyWest(wd.reg,lag=h))
      mat.nl_i_a[j,3] <- mat.nl_i_a[j,3]*harvey
      mat.nl_i_a[j,4] <- 2*pt(abs(mat.nl_i_a[j,3]),df=P-1,lower.tail=F)
    }else{
      mat.nl_i_a[j,] <-  mat.nl_i[j,]
    }
    
  }
  mean.nl_i[,,k] <- mat.nl_i
  mean.nl_i_a[,,k] <- mat.nl_i_a
  
  # (2) direct STAR vs direct AR
  mat.nl_d <- matrix(nrow=h,ncol=4)
  mat.nl_d_a <- matrix(nrow=h,ncol=4)
  for(j in 1:h){
    d.reg <- lm(d.nl_d[j,]~1)
    mat.nl_d[j,] <- coeftest(d.reg,vcov.=NeweyWest(d.reg,lag=h))
    
    harvey <- sqrt((P+1-2*j+j*(j-1)/P)/P)
    mat.nl_d[j,3] <- mat.nl_d[j,3]*harvey
    mat.nl_d[j,4] <- 2*pt(abs(mat.nl_d[j,3]),df=P-1,lower.tail=F)
    
    if(j>=2){
      wd.nl_d <- colMeans(d.nl_d[1:j,])
      wd.reg <- lm(wd.nl_d~1)
      mat.nl_d_a[j,] <- coeftest(wd.reg,vcov.=NeweyWest(wd.reg,lag=h))
      mat.nl_d_a[j,3] <- mat.nl_d_a[j,3]*harvey
      mat.nl_d_a[j,4] <- 2*pt(abs(mat.nl_d_a[j,3]),df=P-1,lower.tail=F)
    }else{
      mat.nl_d_a[j,] <-  mat.nl_d[j,]
    }
    
  }
  mean.nl_d[,,k] <- mat.nl_d
  mean.nl_d_a[,,k] <- mat.nl_d_a
  
  # (3) iterated AR vs direct AR
  mat.id_l <- matrix(nrow=h,ncol=4)
  mat.id_l_a <- matrix(nrow=h,ncol=4)
  for(j in 2:h){
    d.reg <- lm(d.id_l[j,]~1)
    mat.id_l[j,] <- coeftest(d.reg,vcov.=NeweyWest(d.reg,lag=h))
    
    harvey <- sqrt((P+1-2*j+j*(j-1)/P)/P)
    mat.id_l[j,3] <- mat.id_l[j,3]*harvey
    mat.id_l[j,4] <- 2*pt(abs(mat.id_l[j,3]),df=P-1,lower.tail=F)
    
    if(j>2){
      wd.id_l <- colMeans(d.id_l[2:j,])
      wd.reg <- lm(wd.id_l~1)
      mat.id_l_a[j,] <- coeftest(wd.reg,vcov.=NeweyWest(wd.reg,lag=h))
      mat.id_l_a[j,3] <- mat.id_l_a[j,3]*harvey
      mat.id_l_a[j,4] <- 2*pt(abs(mat.id_l_a[j,3]),df=P-1,lower.tail=F)
    }else{
      mat.id_l_a[j,] <-  mat.id_l[j,]
    }
    
  }
  mean.id_l[,,k] <- mat.id_l
  mean.id_l_a[,,k] <- mat.id_l_a

  # (4) iterated STAR vs direct STAR
  mat.id_n <- matrix(nrow=h,ncol=4)
  mat.id_n_a <- matrix(nrow=h,ncol=4)
  for(j in 2:h){
    d.reg <- lm(d.id_n[j,]~1)
    mat.id_n[j,] <- coeftest(d.reg,vcov.=NeweyWest(d.reg,lag=h))
    
    harvey <- sqrt((P+1-2*j+j*(j-1)/P)/P)
    mat.id_n[j,3] <- mat.id_n[j,3]*harvey
    mat.id_n[j,4] <- 2*pt(abs(mat.id_n[j,3]),df=P-1,lower.tail=F)
    
    if(j>2){
      wd.id_n <- colMeans(d.id_n[2:j,])
      wd.reg <- lm(wd.id_n~1)
      mat.id_n_a[j,] <- coeftest(wd.reg,vcov.=NeweyWest(wd.reg,lag=h))
      mat.id_n_a[j,3] <- mat.id_n_a[j,3]*harvey
      mat.id_n_a[j,4] <- 2*pt(abs(mat.id_n_a[j,3]),df=P-1,lower.tail=F)
    }else{
      mat.id_n_a[j,] <-  mat.id_n[j,]
    }
    
  }
  mean.id_n[,,k] <- mat.id_n
  mean.id_n_a[,,k] <- mat.id_n_a

  # (5) iterated STAR vs direct AR
  mat.indl <- matrix(nrow=h,ncol=4)
  mat.indl_a <- matrix(nrow=h,ncol=4)
  for(j in 1:h){
    d.reg <- lm(d.indl[j,]~1)
    mat.indl[j,] <- coeftest(d.reg,vcov.=NeweyWest(d.reg,lag=h))
    
    harvey <- sqrt((P+1-2*j+j*(j-1)/P)/P)
    mat.indl[j,3] <- mat.indl[j,3]*harvey
    mat.indl[j,4] <- 2*pt(abs(mat.indl[j,3]),df=P-1,lower.tail=F)
    
    if(j>=2){
      wd.indl <- colMeans(d.indl[1:j,])
      wd.reg <- lm(wd.indl~1)
      mat.indl_a[j,] <- coeftest(wd.reg,vcov.=NeweyWest(wd.reg,lag=h))
      mat.indl_a[j,3] <- mat.indl_a[j,3]*harvey
      mat.indl_a[j,4] <- 2*pt(abs(mat.indl_a[j,3]),df=P-1,lower.tail=F)
    }else{
      mat.indl_a[j,] <-  mat.indl[j,]
    }
    
  }
  mean.indl[,,k] <- mat.indl
  mean.indl_a[,,k] <- mat.indl_a
  
  # (C1) iterated combined vs direct combined
  mat_i <- matrix(nrow=h,ncol=4)
  mat_i_a <- matrix(nrow=h,ncol=4)
  for(j in 1:h){
    d.reg <- lm(d_i[j,]~1)
    mat_i[j,] <- coeftest(d.reg,vcov.=NeweyWest(d.reg,lag=h))
    
    harvey <- sqrt((P+1-2*j+j*(j-1)/P)/P)
    mat_i[j,3] <- mat_i[j,3]*harvey
    mat_i[j,4] <- 2*pt(abs(mat_i[j,3]),df=P-1,lower.tail=F)
    
    if(j>=2){
      wd_i <- colMeans(d_i[1:j,])
      wd.reg <- lm(wd_i~1)
      mat_i_a[j,] <- coeftest(wd.reg,vcov.=NeweyWest(wd.reg,lag=h))
      mat_i_a[j,3] <- mat_i_a[j,3]*harvey
      mat_i_a[j,4] <- 2*pt(abs(mat_i_a[j,3]),df=P-1,lower.tail=F)
    }else{
      mat_i_a[j,] <-  mat_i[j,]
    }
    
  }
  mean_i[,,k] <- mat_i
  mean_i_a[,,k] <- mat_i_a
  
  # (C2) nonlinear combined vs linear combined
  mat_n <- matrix(nrow=h,ncol=4)
  mat_n_a <- matrix(nrow=h,ncol=4)
  for(j in 1:h){
    d.reg <- lm(d_n[j,]~1)
    mat_n[j,] <- coeftest(d.reg,vcov.=NeweyWest(d.reg,lag=h))
    
    harvey <- sqrt((P+1-2*j+j*(j-1)/P)/P)
    mat_n[j,3] <- mat_n[j,3]*harvey
    mat_n[j,4] <- 2*pt(abs(mat_n[j,3]),df=P-1,lower.tail=F)
    
    if(j>=2){
      wd_n <- colMeans(d_n[1:j,])
      wd.reg <- lm(wd_n~1)
      mat_n_a[j,] <- coeftest(wd.reg,vcov.=NeweyWest(wd.reg,lag=h))
      mat_n_a[j,3] <- mat_n_a[j,3]*harvey
      mat_n_a[j,4] <- 2*pt(abs(mat_n_a[j,3]),df=P-1,lower.tail=F)
    }else{
      mat_n_a[j,] <-  mat_n[j,]
    }
    
  }
  mean_n[,,k] <- mat_n
  mean_n_a[,,k] <- mat_n_a
  
}

# (1) iterated STAR vs iterated AR
rat.nl_i_mean <- round(t(rmsfe_mean[,2,]/rmsfe_mean[,1,]),3)
dm.nl_i_mean <- round(t(mean.nl_i[,4,]),3)
dm.nl_i_mean_a <- round(t(mean.nl_i_a[,4,]),3)

# (2) direct STAR vs direct AR
rat.nl_d_mean <- round(t(rmsfe_mean[,4,]/rmsfe_mean[,3,]),3)
dm.nl_d_mean <- round(t(mean.nl_d[,4,]),3)
dm.nl_d_mean_a <- round(t(mean.nl_d_a[,4,]),3)

# (3) iterated AR vs direct AR
rat.id_l_mean <- round(t(rmsfe_mean[,1,]/rmsfe_mean[,3,]),3)
dm.id_l_mean <- round(t(mean.id_l[,4,]),3)
dm.id_l_mean_a <- round(t(mean.id_l_a[,4,]),3)

# (4) iterated STAR vs direct STAR
rat.id_n_mean <- round(t(rmsfe_mean[,2,]/rmsfe_mean[,4,]),3)
dm.id_n_mean <- round(t(mean.id_n[,4,]),3)
dm.id_n_mean_a <- round(t(mean.id_n_a[,4,]),3)

# (5) iterated STAR vs direct AR
rat.indl_mean <- round(t(rmsfe_mean[,2,]/rmsfe_mean[,3,]),3)
dm.indl_mean <- round(t(mean.indl[,4,]),3)
dm.indl_mean_a <- round(t(mean.indl_a[,4,]),3)

# (C1) iterated combined vs direct combined
rat.i_mean <- round(t(rmsfe_mean[,5,]/rmsfe_mean[,6,]),3)
dm.i_mean <- round(t(mean_i[,4,]),3)
dm.i_mean_a <- round(t(mean_i_a[,4,]),3)

# (C2) nonlinear combined vs linear combined
rat.n_mean <- round(t(rmsfe_mean[,8,]/rmsfe_mean[,7,]),3)
dm.n_mean <- round(t(mean_n[,4,]),3)
dm.n_mean_a <- round(t(mean_n_a[,4,]),3)

## commodity ordering
## this is sort of arbitrary, but in general: 
## agricultural--forestry--metals and minerals
cmd_sequence <- c(5,.9,1,4,2,
                  2,2,5,1.9,1.8,
                  5,3,.9,2,5,
                  3,2,1,2,2.5,
                  2,3,3,5,4)


################
###  TABLES  ###
################

# (1) iterated STAR vs iterated AR
tab1_dt <- data.table(Commodity=commodities,rat.nl_i_mean)
tab1_dt$grp <- cmd_sequence
tab1_dt <- tab1_dt[order(grp,Commodity)]
tab1_dt$grp <- NULL
colnames(tab1_dt) <- c("Commodity",paste0("h",c(1:12)))

dm1_dt <- data.table(Commodity=commodities,dm.nl_i_mean)
dm1_dt$grp <- cmd_sequence
dm1_dt <- dm1_dt[order(grp,Commodity)]
dm1_dt$grp <- NULL
colnames(dm1_dt) <- c("Commodity",paste0("h",c(1:12)))

dma1_dt <- data.table(Commodity=commodities,dm.nl_i_mean_a)
dma1_dt$grp <- cmd_sequence
dma1_dt <- dma1_dt[order(grp,Commodity)]
dma1_dt$grp <- NULL
colnames(dma1_dt) <- c("Commodity",paste0("h",c(1:12)))

tab1_dt$Methods <- "iSTAR/iAR"
tab1_lg <- melt(tab1_dt,id.vars=c("Commodity","Methods"),variable.name="Horizon",value.name="Ratio")

tab1_dt$Methods <- NULL
print(xtable(tab1_dt),include.rownames=FALSE)
write.table(tab1_dt,file="Paper/table_2.txt",sep=",",quote=F,row.names=F)


# (2) direct STAR vs direct AR
tab2_dt <- data.table(Commodity=commodities,rat.nl_d_mean)
tab2_dt$grp <- cmd_sequence
tab2_dt <- tab2_dt[order(grp,Commodity)]
tab2_dt$grp <- NULL
colnames(tab2_dt) <- c("Commodity",paste0("h",c(1:12)))

dm2_dt <- data.table(Commodity=commodities,dm.nl_d_mean)
dm2_dt$grp <- cmd_sequence
dm2_dt <- dm2_dt[order(grp,Commodity)]
dm2_dt$grp <- NULL
colnames(dm2_dt) <- c("Commodity",paste0("h",c(1:12)))

dma2_dt <- data.table(Commodity=commodities,dm.nl_d_mean_a)
dma2_dt$grp <- cmd_sequence
dma2_dt <- dma2_dt[order(grp,Commodity)]
dma2_dt$grp <- NULL
colnames(dma2_dt) <- c("Commodity",paste0("h",c(1:12)))

tab2_dt$Methods <- "dSTAR/dAR"
tab2_lg <- melt(tab2_dt,id.vars=c("Commodity","Methods"),variable.name="Horizon",value.name="Ratio")

tab2_dt$Methods <- NULL
print(xtable(tab2_dt),include.rownames=FALSE)
write.table(tab2_dt,file="Paper/table_3.txt",sep=",",quote=F,row.names=F)


# (3) iterated AR vs direct AR
tab3_dt <- data.table(Commodity=commodities,rat.id_l_mean)
tab3_dt$grp <- cmd_sequence
tab3_dt <- tab3_dt[order(grp,Commodity)]
tab3_dt$grp <- NULL
colnames(tab3_dt) <- c("Commodity",paste0("h",c(1:12)))

dm3_dt <- data.table(Commodity=commodities,dm.id_l_mean)
dm3_dt$grp <- cmd_sequence
dm3_dt <- dm3_dt[order(grp,Commodity)]
dm3_dt$grp <- NULL
colnames(dm3_dt) <- c("Commodity",paste0("h",c(1:12)))

dma3_dt <- data.table(Commodity=commodities,dm.id_l_mean_a)
dma3_dt$grp <- cmd_sequence
dma3_dt <- dma3_dt[order(grp,Commodity)]
dma3_dt$grp <- NULL
colnames(dma3_dt) <- c("Commodity",paste0("h",c(1:12)))

tab3_dt$Methods <- "iAR/dAR"
tab3_lg <- melt(tab3_dt,id.vars=c("Commodity","Methods"),variable.name="Horizon",value.name="Ratio")

tab3_dt$Methods <- NULL
print(xtable(tab3_dt),include.rownames=FALSE)
write.table(tab3_dt,file="Paper/table_6.txt",sep=",",quote=F,row.names=F)


# (4) iterated STAR vs direct STAR
tab4_dt <- data.table(Commodity=commodities,rat.id_n_mean)
tab4_dt$grp <- cmd_sequence
tab4_dt <- tab4_dt[order(grp,Commodity)]
tab4_dt$grp <- NULL
colnames(tab4_dt) <- c("Commodity",paste0("h",c(1:12)))

dm4_dt <- data.table(Commodity=commodities,dm.id_n_mean)
dm4_dt$grp <- cmd_sequence
dm4_dt <- dm4_dt[order(grp,Commodity)]
dm4_dt$grp <- NULL
colnames(dm4_dt) <- c("Commodity",paste0("h",c(1:12)))

dma4_dt <- data.table(Commodity=commodities,dm.id_n_mean_a)
dma4_dt$grp <- cmd_sequence
dma4_dt <- dma4_dt[order(grp,Commodity)]
dma4_dt$grp <- NULL
colnames(dma4_dt) <- c("Commodity",paste0("h",c(1:12)))

tab4_dt$Methods <- "iSTAR/dSTAR"
tab4_lg <- melt(tab4_dt,id.vars=c("Commodity","Methods"),variable.name="Horizon",value.name="Ratio")

tab4_dt$Methods <- NULL
print(xtable(tab4_dt),include.rownames=FALSE)
write.table(tab4_dt,file="Paper/table_4.txt",sep=",",quote=F,row.names=F)


# (5) iterated STAR vs direct AR
tab5_dt <- data.table(Commodity=commodities,rat.indl_mean)
tab5_dt$grp <- cmd_sequence
tab5_dt <- tab5_dt[order(grp,Commodity)]
tab5_dt$grp <- NULL
colnames(tab5_dt) <- c("Commodity",paste0("h",c(1:12)))

dm5_dt <- data.table(Commodity=commodities,dm.indl_mean)
dm5_dt$grp <- cmd_sequence
dm5_dt <- dm5_dt[order(grp,Commodity)]
dm5_dt$grp <- NULL
colnames(dm5_dt) <- c("Commodity",paste0("h",c(1:12)))

dma5_dt <- data.table(Commodity=commodities,dm.indl_mean_a)
dma5_dt$grp <- cmd_sequence
dma5_dt <- dma5_dt[order(grp,Commodity)]
dma5_dt$grp <- NULL
colnames(dma5_dt) <- c("Commodity",paste0("h",c(1:12)))

tab5_dt$Methods <- "iSTAR/dAR"
tab5_lg <- melt(tab5_dt,id.vars=c("Commodity","Methods"),variable.name="Horizon",value.name="Ratio")

tab5_dt$Methods <- NULL
print(xtable(tab5_dt),include.rownames=FALSE)
write.table(tab5_dt,file="Paper/table_5.txt",sep=",",quote=F,row.names=F)


################
###  FIGURE  ###
################

tab_lg <- rbind(tab1_lg,tab2_lg,tab4_lg,tab5_lg,tab3_lg)

tab_lg$Commodity <- factor(tab_lg$Commodity,levels=unique(tab_lg$Commodity))
tab_lg$Methods <- factor(tab_lg$Methods,levels=unique(tab_lg$Methods))

gg_ratios <- ggplot(tab_lg,aes(x=Horizon,y=Ratio,color=Methods,linetype=Methods,group=Methods))+
  geom_line(size=0.4)+
  geom_hline(yintercept = 1,color="gray70",linetype=5,size=.2)+
  facet_wrap(~Commodity,ncol=5,scale="free")+
  scale_color_manual(values=c("gray30","steelblue","indianred","seagreen","goldenrod"))+
  scale_x_discrete(breaks=c("h3","h6","h9","h12"),labels=c(3,6,9,12))+
  scale_y_continuous(labels = scales::number_format(accuracy=0.01))+
  labs(x="Forecast Horizon",y="RMSFE Ratio")+
  theme_classic()+
  theme(axis.title = element_text(size=10),axis.text.x=element_text(size=6,hjust=.8),axis.text.y=element_text(size=6),strip.text=element_text(size=6),strip.background=element_blank(),legend.position = "top",legend.title=element_blank())


ggsave(paste("Paper/Forecasts",w_frac,".eps",sep=""),gg_ratios,width=6.5,height=6.5,units="in",device="eps")
ggsave(paste("Paper/Forecasts",w_frac,".png",sep=""),gg_ratios,width=6.5,height=6.5,units="in",device="png")
ggsave(paste("Presentation/Forecasts",w_frac,".png",sep=""),gg_ratios,width=6.5,height=4.5,units="in",device="png")

