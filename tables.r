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
library(xtable)
library(fixest)

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

count_fn <- function(x,crit=.05){
  sum(x < crit)/length(x)
}

# table 1
tab1_dt <- as.data.table(t(round(apply(array_pval,c(1,3),count_fn),2)))
colnames(tab1_dt) <- paste0("h",1:12)
tab1_dt$cmd <- commodities
tab1_dt$grp <- grp_order
tab1_dt <- tab1_dt[order(grp,cmd)]

tab1_dt$grp <- NULL
tab1_dt <- setcolorder(tab1_dt,c("cmd",paste0("h",1:12)))
print(xtable(tab1_dt),include.rownames=FALSE)
write.table(tab1_dt,file="Paper/table_1.txt",sep=",",quote=F,row.names=F)

# appendix table 1
taba1_dt <- sub_dt[,.(Commodity=cmd,p,l,d,t)]
print(xtable(taba1_dt),include.rownames=FALSE)
write.table(taba1_dt,file="Paper/table_a1.txt",sep=",",quote=F,row.names=F)





