library(ggplot2)
library(reshape2)
library(sandwich)
library(lmtest)
library(tseries)
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


c_num <- ncol(subset_dt)-1

c_names <- names(subset_dt)[1:c_num]

commodities <- c("Aluminum","Bananas","Cotton","Diammonium Phosphate","Fishmeal",
                 "Groundnut Oil","Groundnuts","Gold","Hides","Lamb",
                 "Lead","Soft Logs","Oranges","Olive Oil","Platinum",
                 "Plywood","Palm Oil","Rice","Rapeseed Oil","Rubber",
                 "Soybean Oil","Hard Sawnwood","Soft Sawnwood","Uran","Urea")

sub_dt$cmd <- commodities

## this is sort of arbitrary, but 
# in general: agricultural(1)--forestry(2,3)--metals & minerals(4,5)
grp_order <- c(5,.9,1,4,2,
               2,2,5,1.9,1.8,
               5,3,.9,2,5,
               3,2,1,2,2.5,
               2,3,3,5,4)
sub_dt$grp <- grp_order
sub_dt <- sub_dt[order(grp,cmd)]

colnames(subset_dt) <- c(commodities,"Date")

subset_lg <- melt(subset_dt,id.vars = "Date", variable.name = "Commodity",value.name = "Price")

subset_lg$Commodity <- factor(subset_lg$Commodity,levels=sub_dt$cmd)

gg_panel <- ggplot(subset_lg,aes(x=Date,y=log(Price)))+
  geom_line(size=0.4,color="gray30")+
  facet_wrap(~Commodity,scales="free",ncol=5)+
  scale_y_continuous(labels = scales::number_format(accuracy=0.1))+
  labs(x="Year",y="Price")+
  theme_classic()+
  theme(axis.title = element_text(size=10),axis.text.x=element_text(size=6,angle=45,hjust=.8),axis.text.y=element_text(size=6),strip.text=element_text(size=6),strip.background=element_blank())

# ggsave("Paper/Prices.eps",gg_panel,width=6.5,height=6.5,units="in",device="eps")
# ggsave("Paper/Prices.png",gg_panel,width=6.5,height=6.5,units="in",device="png")
# ggsave("Presentation/Prices.png",gg_panel,width=6.5,height=4.5,units="in",device="png")


forecast_ar <- array_forecast[,,63,]

dimnames(forecast_ar) <- list(month.abb,c("realised","rw","iAR","iSTAR","dAR","dSTAR"),commodities)

forecast_dt <- as.data.table(forecast_ar)

colnames(forecast_dt) <- c("horizon","series","commodity","value")
forecast_dt <- forecast_dt[series!="rw"]

forecast_dt$commodity <- factor(forecast_dt$commodity,levels=sub_dt$cmd)
forecast_dt$horizon <- factor(forecast_dt$horizon,levels=month.abb)
forecast_dt$series <- factor(forecast_dt$series,levels=unique(forecast_dt$series)[c(5,4,3,2,1)])

forecast_dt$value <- as.numeric(forecast_dt$value)

gg_forecast <- ggplot(forecast_dt,aes(x=horizon,y=value,color=series,linetype=series,group=series))+
  geom_line(size=0.4)+
  facet_wrap(~commodity,scales="free",ncol=5)+
  scale_color_manual(values=c("darkgray","indianred","steelblue","goldenrod","seagreen"))+
  scale_x_discrete(breaks=month.abb[c(3,6,9,12)],labels=month.abb[c(3,6,9,12)])+
  scale_y_continuous(labels = scales::number_format(accuracy=0.01))+
  labs(x="Month",y="Price")+
  theme_classic()+
  theme(axis.title = element_text(size=10),axis.text.x=element_text(size=6),axis.text.y=element_text(size=6),strip.text=element_text(size=6),strip.background=element_blank(),legend.position = "top",legend.title=element_blank())

# ggsave("Paper/Forecasts2018.eps",gg_forecast,width=6.5,height=6.5,units="in",device="eps")
# ggsave("Paper/Forecasts2018.png",gg_forecast,width=6.5,height=6.5,units="in",device="png")
# ggsave("Presentation/Forecasts2018.png",gg_forecast,width=6.5,height=4.5,units="in",device="png")


