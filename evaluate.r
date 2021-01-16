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

rm(list=ls())
gc()

## load the data
load("forecast80.RData")
load("subset.RData")
load("select.RData")

h=12

P <- dim(array.forecast)[3]

c.num <- ncol(subset_dt)-1

c.names <- as.character(diagnostics_dt$c)
commodities <- c("Aluminum","Beef","Barley","Coffee (Arabicas)","Cocoa","Crude Oil","Copper","DAP","Fish","Fishmeal","Groundnut Oil","Groundnuts","Gold","Hides","Lamb","Lead","Soft Logs","Olive Oil","Platinum","Plywood","Palm Oil","Rice","Rubber","Soybean Oil","Shrimp","Hard Sawnwood","Tin","TSP","Uranium","Urea")

rmsfe_mean <- array(dim=c(h,5,c.num))

mean.nl_i <- array(dim=c(h,4,c.num))
mean.nl_d <- array(dim=c(h,4,c.num))
mean.id_l <- array(dim=c(h,4,c.num))
mean.id_n <- array(dim=c(h,4,c.num))
mean.indl <- array(dim=c(h,4,c.num))
mean.nrw <- array(dim=c(h,4,c.num))
mean.lrw <- array(dim=c(h,4,c.num))

mean.nl_i_a <- array(dim=c(h,4,c.num))
mean.nl_d_a <- array(dim=c(h,4,c.num))
mean.id_l_a <- array(dim=c(h,4,c.num))
mean.id_n_a <- array(dim=c(h,4,c.num))
mean.indl_a <- array(dim=c(h,4,c.num))
mean.nrw_a <- array(dim=c(h,4,c.num))
mean.lrw_a <- array(dim=c(h,4,c.num))

for(k in 1:c.num){
  
  e.rw      <- array.forecast[,1,,k]-array.forecast[,2,,k]
  e.il_mean <- array.forecast[,1,,k]-array.forecast[,3,,k]
  e.in_mean <- array.forecast[,1,,k]-array.forecast[,4,,k]
  e.dl_mean <- array.forecast[,1,,k]-array.forecast[,5,,k]
  e.dn_mean <- array.forecast[,1,,k]-array.forecast[,6,,k]

  rmsfe.rw <- apply(e.rw,1,function(x) sqrt(mean(x^2)))
  rmsfe.il_mean <- apply(e.il_mean,1,function(x) sqrt(mean(x^2)))
  rmsfe.in_mean <- apply(e.in_mean,1,function(x) sqrt(mean(x^2)))
  rmsfe.dl_mean <- apply(e.dl_mean,1,function(x) sqrt(mean(x^2)))
  rmsfe.dn_mean <- apply(e.dn_mean,1,function(x) sqrt(mean(x^2)))

  rmsfe_mean[,,k] <- cbind(rmsfe.il_mean,rmsfe.in_mean,rmsfe.dl_mean,rmsfe.dn_mean,rmsfe.rw)

  
  d.nl_i <- e.in_mean^2-e.il_mean^2
  d.nl_d <- e.dn_mean^2-e.dl_mean^2
  d.id_l <- e.il_mean^2-e.dl_mean^2
  d.id_n <- e.in_mean^2-e.dn_mean^2
  d.indl <- e.in_mean^2-e.dl_mean^2
  d.lrw <- e.il_mean^2-e.rw^2
  d.nrw <- e.in_mean^2-e.rw^2
  
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

  
  mat.indl <- matrix(nrow=h,ncol=4)
  mat.indl_a <- matrix(nrow=h,ncol=4)
  for(j in 2:h){
    d.reg <- lm(d.indl[j,]~1)
    mat.indl[j,] <- coeftest(d.reg,vcov.=NeweyWest(d.reg,lag=h))
    
    harvey <- sqrt((P+1-2*j+j*(j-1)/P)/P)
    mat.indl[j,3] <- mat.indl[j,3]*harvey
    mat.indl[j,4] <- 2*pt(abs(mat.indl[j,3]),df=P-1,lower.tail=F)
    
    if(j>2){
      wd.indl <- colMeans(d.indl[2:j,])
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
  
  
  mat.lrw <- matrix(nrow=h,ncol=4)
  mat.lrw_a <- matrix(nrow=h,ncol=4)
  for(j in 1:h){
    d.reg <- lm(d.lrw[j,]~1)
    mat.lrw[j,] <- coeftest(d.reg,vcov.=NeweyWest(d.reg,lag=h))
    
    harvey <- sqrt((P+1-2*j+j*(j-1)/P)/P)
    mat.lrw[j,3] <- mat.lrw[j,3]*harvey
    mat.lrw[j,4] <- 2*pt(abs(mat.lrw[j,3]),df=P-1,lower.tail=F)
    
    if(j>=2){
      wd.lrw <- colMeans(d.lrw[1:j,])
      wd.reg <- lm(wd.lrw~1)
      mat.lrw_a[j,] <- coeftest(wd.reg,vcov.=NeweyWest(wd.reg,lag=h))
      mat.lrw_a[j,3] <- mat.lrw_a[j,3]*harvey
      mat.lrw_a[j,4] <- 2*pt(abs(mat.lrw_a[j,3]),df=P-1,lower.tail=F)
    }else{
      mat.lrw_a[j,] <-  mat.lrw[j,]
    }
    
  }
  mean.lrw[,,k] <- mat.lrw
  mean.lrw_a[,,k] <- mat.lrw_a
  
  
  mat.nrw <- matrix(nrow=h,ncol=4)
  mat.nrw_a <- matrix(nrow=h,ncol=4)
  for(j in 1:h){
    d.reg <- lm(d.nrw[j,]~1)
    mat.nrw[j,] <- coeftest(d.reg,vcov.=NeweyWest(d.reg,lag=h))
    
    harvey <- sqrt((P+1-2*j+j*(j-1)/P)/P)
    mat.nrw[j,3] <- mat.nrw[j,3]*harvey
    mat.nrw[j,4] <- 2*pt(abs(mat.nrw[j,3]),df=P-1,lower.tail=F)
    
    if(j>=2){
      wd.nrw <- colMeans(d.nrw[1:j,])
      wd.reg <- lm(wd.nrw~1)
      mat.nrw_a[j,] <- coeftest(wd.reg,vcov.=NeweyWest(wd.reg,lag=h))
      mat.nrw_a[j,3] <- mat.nrw_a[j,3]*harvey
      mat.nrw_a[j,4] <- 2*pt(abs(mat.nrw_a[j,3]),df=P-1,lower.tail=F)
    }else{
      mat.nrw_a[j,] <-  mat.nrw[j,]
    }
    
  }
  mean.nrw[,,k] <- mat.nrw
  mean.nrw_a[,,k] <- mat.nrw_a
  
}


## rmsfe: il, in, dl, dn


rat.nl_i_mean <- round(t(rmsfe_mean[,2,]/rmsfe_mean[,1,]),3)
dm.nl_i_mean <- round(t(mean.nl_i[,4,]),3)
dm.nl_i_mean_a <- round(t(mean.nl_i_a[,4,]),3)

rat.nl_d_mean <- round(t(rmsfe_mean[,4,]/rmsfe_mean[,3,]),3)
dm.nl_d_mean <- round(t(mean.nl_d[,4,]),3)
dm.nl_d_mean_a <- round(t(mean.nl_d_a[,4,]),3)

rat.id_n_mean <- round(t(rmsfe_mean[,2,]/rmsfe_mean[,4,]),3)
dm.id_n_mean <- round(t(mean.id_n[,4,]),3)
dm.id_n_mean_a <- round(t(mean.id_n_a[,4,]),3)

rat.id_l_mean <- round(t(rmsfe_mean[,1,]/rmsfe_mean[,3,]),3)
dm.id_l_mean <- round(t(mean.id_l[,4,]),3)
dm.id_l_mean_a <- round(t(mean.id_l_a[,4,]),3)

rat.indl_mean <- round(t(rmsfe_mean[,2,]/rmsfe_mean[,3,]),3)
dm.indl_mean <- round(t(mean.indl[,4,]),3)
dm.indl_mean_a <- round(t(mean.indl_a[,4,]),3)

rat.lrw_mean <- round(t(rmsfe_mean[,1,]/rmsfe_mean[,5,]),3)
dm.lrw_mean <- round(t(mean.lrw[,4,]),3)
dm.lrw_mean_a <- round(t(mean.lrw_a[,4,]),3)

rat.nrw_mean <- round(t(rmsfe_mean[,2,]/rmsfe_mean[,5,]),3)
dm.nrw_mean <- round(t(mean.nrw[,4,]),3)
dm.nrw_mean_a <- round(t(mean.nrw_a[,4,]),3)


cmd_sequence <- c(5,3,0,2,2,6,5,6,3,1,1,1,5,3,3,5,4,1,5,4,1,0,4,1,3,4,5,6,6,6)

tab1_dt <- data.table(Commodity=commodities,nrw=rat.nrw_mean[,1],rat.id_n_mean[,-1])
tab1_dt$grp <- cmd_sequence
tab1_dt <- tab1_dt[order(grp,Commodity)]
tab1_dt$grp <- NULL
colnames(tab1_dt) <- c("Commodity",paste0("h",c(1:12)))
# stargazer(as.matrix(tab1_dt),summary=F)
tab1_dt$Methods <- "iSTAR/dSTAR"
tab1_dt$h1 <- NA
tab1_lg <- melt(tab1_dt,id.vars=c("Commodity","Methods"),variable.name="Horizon",value.name="Ratio")

# tab1_dt <- data.table(Commodity=commodities,nrw=dm.nrw_mean[,1],dm.id_n_mean[,-1])
# tab1_dt$grp <- cmd_sequence
# tab1_dt <- tab1_dt[order(grp,Commodity)]
# tab1_dt$grp <- NULL
# cmdty <- tab1_dt$Commodity
# tab1_dt[tab1_dt >= 0.05] <- NA
# tab1_dt$Commodity <- cmdty
# 
# tab1_dt <- data.table(Commodity=commodities,nrw=dm.nrw_mean_a[,1],dm.id_n_mean_a[,-1])
# tab1_dt$grp <- cmd_sequence
# tab1_dt <- tab1_dt[order(grp,Commodity)]
# tab1_dt$grp <- NULL
# cmdty <- tab1_dt$Commodity
# tab1_dt[tab1_dt >= 0.05] <- NA
# tab1_dt$Commodity <- cmdty


tab2_dt <- data.table(Commodity=commodities,rat.indl_mean)
tab2_dt$grp <- cmd_sequence
tab2_dt <- tab2_dt[order(grp,Commodity)]
tab2_dt$grp <- NULL
colnames(tab2_dt) <- c("Commodity",paste0("h",c(1:12)))
# stargazer(as.matrix(tab2_dt),summary=F)
tab2_dt$Methods <- "iSTAR/dAR"
tab2_dt$h1 <- NA
tab2_lg <- melt(tab2_dt,id.vars=c("Commodity","Methods"),variable.name="Horizon",value.name="Ratio")

# tab2_dt <- data.table(Commodity=commodities,dm.indl_mean)
# tab2_dt$grp <- cmd_sequence
# tab2_dt <- tab2_dt[order(grp,Commodity)]
# tab2_dt$grp <- NULL
# cmdty <- tab2_dt$Commodity
# tab2_dt[tab2_dt >= 0.05] <- NA
# tab2_dt$Commodity <- cmdty
# 
# tab2_dt <- data.table(Commodity=commodities,dm.indl_mean_a)
# tab2_dt$grp <- cmd_sequence
# tab2_dt <- tab2_dt[order(grp,Commodity)]
# tab2_dt$grp <- NULL
# cmdty <- tab2_dt$Commodity
# tab2_dt[tab2_dt >= 0.05] <- NA
# tab2_dt$Commodity <- cmdty


tab3_dt <- data.table(Commodity=commodities,rat.nl_i_mean)
tab3_dt$grp <- cmd_sequence
tab3_dt <- tab3_dt[order(grp,Commodity)]
tab3_dt$grp <- NULL
colnames(tab3_dt) <- c("Commodity",paste0("h",c(1:12)))
# stargazer(as.matrix(tab3_dt),summary=F)
tab3_dt$Methods <- "iSTAR/iAR"
tab3_dt$h1 <- NA
tab3_lg <- melt(tab3_dt,id.vars=c("Commodity","Methods"),variable.name="Horizon",value.name="Ratio")

# tab3_dt <- data.table(Commodity=commodities,dm.nl_i_mean)
# tab3_dt$grp <- cmd_sequence
# tab3_dt <- tab3_dt[order(grp,Commodity)]
# tab3_dt$grp <- NULL
# cmdty <- tab3_dt$Commodity
# tab3_dt[tab3_dt >= 0.05] <- NA
# tab3_dt$Commodity <- cmdty
# 
# tab3_dt <- data.table(Commodity=commodities,dm.nl_i_mean_a)
# tab3_dt$grp <- cmd_sequence
# tab3_dt <- tab3_dt[order(grp,Commodity)]
# tab3_dt$grp <- NULL
# cmdty <- tab3_dt$Commodity
# tab3_dt[tab3_dt >= 0.05] <- NA
# tab3_dt$Commodity <- cmdty


tab4_dt <- data.table(Commodity=commodities,lrw=rat.lrw_mean[,1],rat.id_l_mean[,-1])
tab4_dt$grp <- cmd_sequence
tab4_dt <- tab4_dt[order(grp,Commodity)]
tab4_dt$grp <- NULL
colnames(tab4_dt) <- c("Commodity",paste0("h",c(1:12)))
# stargazer(as.matrix(tab4_dt),summary=F)
tab4_dt$Methods <- "iAR/dAR"
tab4_dt$h1 <- NA
tab4_lg <- melt(tab4_dt,id.vars=c("Commodity","Methods"),variable.name="Horizon",value.name="Ratio")

# tab4_dt <- data.table(Commodity=commodities,lrw=dm.lrw_mean[,1],dm.id_l_mean[,-1])
# tab4_dt$grp <- cmd_sequence
# tab4_dt <- tab4_dt[order(grp,Commodity)]
# tab4_dt$grp <- NULL
# cmdty <- tab4_dt$Commodity
# tab4_dt[tab4_dt >= 0.05] <- NA
# tab4_dt$Commodity <- cmdty
# 
# tab4_dt <- data.table(Commodity=commodities,lrw=dm.lrw_mean_a[,1],dm.id_l_mean_a[,-1])
# tab4_dt$grp <- cmd_sequence
# tab4_dt <- tab4_dt[order(grp,Commodity)]
# tab4_dt$grp <- NULL
# cmdty <- tab4_dt$Commodity
# tab4_dt[tab4_dt >= 0.05] <- NA
# tab4_dt$Commodity <- cmdty


tab_lg <- rbind(tab1_lg,tab2_lg,tab3_lg)

tab_lg$Commodity <- factor(tab_lg$Commodity,levels=unique(tab_lg$Commodity))
tab_lg$Methods <- factor(tab_lg$Methods,levels=unique(tab_lg$Methods))

gg_ratios <- ggplot(tab_lg,aes(x=Horizon,y=Ratio,color=Methods,linetype=Methods,group=Methods))+
  geom_line(size=0.4)+
  geom_hline(yintercept = 1,color="gray70",linetype=5,size=.2)+
  facet_wrap(~Commodity,ncol=5)+
  coord_cartesian(ylim=c(.85,1.15))+
  scale_color_manual(values=c("gray30","steelblue","indianred"))+
  scale_x_discrete(breaks=c("h3","h6","h9","h12"),labels=c(3,6,9,12))+
  labs(x="Forecast Horizon",y="RMSFE Ratio")+
  theme_classic()+
  theme(axis.title = element_text(size=10),axis.text.x=element_text(size=6,hjust=.8),axis.text.y=element_text(size=6),strip.text=element_text(size=7),strip.background=element_blank(),legend.position = "bottom",legend.title=element_blank())

ggsave("Forecasts80.png",gg_ratios,width=6.5,height=6.5,units="in")

# star_dt <- data.table(Commodity=commodities,star_pval)
# star_dt$grp <- cmd_sequence
# star_dt <- star_dt[order(grp,Commodity)]
# star_dt$grp <- NULL
# cmdty <- star_dt$Commodity
# star_dt[star_dt < 0.001] <- "0.001$^*$"
# star_dt$Commodity <- cmdty

