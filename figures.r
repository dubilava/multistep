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
library(stringr)

rm(list=ls())
gc()

source("startest.r")
source("lstar.r")
source("estar.r")

## load the data
load("forecast80.RData")
load("subset.RData")
load("select.RData")


c.names <- sub_dt$c

commodities <- c("Aluminum","Beef","Barley","Coffee (Arabicas)","Cocoa","Coconut Oil","Crude Oil","Copper","DAP","Fish","Fishmeal","Groundnut Oil","Gold","Hides","Lamb","Lead","Hard Logs","Soft Logs","Platinum","Plywood","Palm Oil","Rice","Rapeseed Oil","Rubber","Soybeans","Soybean Meal","Silver","Sorghum","Shrimp","Hard Sawnwood","Soft Sawnwood","Tin","TSP","Uranium","Wool (Coarse)","Wool (Fine)")

sub_dt$cmd <- commodities

# groups: grains=1, oils=2, other ag=3, forest=4, metals=5, other=6
grp_order <- c(5,3,1,3,3,2,6,5,6,3,2,2,4,3,3,5,4,4,5,4,2,1,2,4,1,2,5,1,3,4,4,5,6,6,3,3) 
sub_dt$grp <- grp_order
sub_dt <- sub_dt[order(grp,cmd)]

colnames(subset_dt) <- c(commodities,"Date")

subset_lg <- melt(subset_dt,id.vars = "Date", variable.name = "Commodity",value.name = "Price")

subset_lg$Commodity <- factor(subset_lg$Commodity,levels=sub_dt$cmd)

gg_panel <- ggplot(subset_lg,aes(x=Date,y=log(Price)))+
  geom_line(size=0.5)+
  facet_wrap(~Commodity,scales="free",ncol=6)+
  labs(x="Year",y="Price")+
  theme_classic()+
  theme(axis.title = element_text(size=10),axis.text.x=element_text(size=6,angle=45,hjust=.8),axis.text.y=element_text(size=6),strip.text=element_text(size=7),strip.background=element_blank())

ggsave("Prices.png",gg_panel,width=6.5,height=6.5,units="in")


transfunc_lg <- Reduce(rbind,star_ls)

for(i in 1:nrow(sub_dt)){
  transfunc_lg$Commodity <- gsub(sub_dt$c[i],sub_dt$cmd[i],transfunc_lg$Commodity)
}

transfunc_lg$Commodity <- factor(transfunc_lg$Commodity,levels=sub_dt$cmd)

gg_tf <- ggplot(transfunc_lg,aes(x=s.t,y=G))+
  geom_point(size=0.2)+
  facet_wrap(~Commodity,scales="free_x",ncol=6)+
  labs(x=bquote(s[t]),y=bquote(G(s[t])))+
  theme_classic()+
  theme(axis.title = element_text(size=10),axis.text.x=element_text(size=6,angle=45,hjust=.8),axis.text.y=element_text(size=6),strip.text=element_text(size=7),strip.background=element_blank())

ggsave("Transition_Functions.png",gg_tf,width=6.5,height=6.5,units="in")



