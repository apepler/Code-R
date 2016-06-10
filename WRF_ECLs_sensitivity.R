rm(list=ls())

library(RNetCDF)
library(fields)
library(akima)
library("R.matlab")
readMat('~/Documents/Data/Useful.mat')->Useful
mask<-t(Useful$mask)
mask[is.na(mask)]=0
source('~/R/x86_64-pc-linux-gnu-library/2.15/gplots/R/rich.colors.R')

setwd("/srv/ccrc/data36/z3478332/WRF/output/ERAI_R2_ensemble_notopo/out/")

makePDF = function(data1,data2,xlabel,labloc="topleft") {
  a=density(data1,na.rm=T)
  b=density(data2,na.rm=T)
  
  lims=range(data1,data2,na.rm=T)
  if((lims[2]-lims[1])<10)
  {
    lims[1]=floor(lims[1])
    lims[2]=ceiling(lims[2])
  } else {
    lims[1]=floor(lims[1]/5)*5
    lims[2]=ceiling(lims[2]/5)*5
  }
  
  plot(a,col=rgb(0,0,1,1/4),xlim=lims,ylim=range(0,a$y,b$y),
       xlab=xlabel,ylab="Frequency",cex.main=1.2,
       main="June 2007 ECL statistics for 20 runs")
  polygon(a,col=rgb(0,0,1,1/4),density=-1)
  polygon(b,col=rgb(1,0,0,1/4),density=-1)
  legend(labloc,legend=c("Control","NoTopo"),
         col=c(rgb(0,0,1,1/4),rgb(1,0,0,1/4)),lwd=4,cex=1,bty="n")   
}

events<-fixes<-events_notopo<-fixes_notopo <- list()
dir="/home/nfs/z3478332/output/outputUM_wrf_2007-06/ERAI_R2_ensemble_notopo/"
names=c("2007052700_v2","2007052700_d02_v2")
for(n in 1:2)
  {
    events[[n]]=read.csv(paste(dir,"ECLevents_",names[n],".csv",sep=""))
    fixes[[n]]=read.csv(paste(dir,"ECLfixes_",names[n],".csv",sep=""))
    fixes[[n]]$Date2=fixes[[n]]$Date+(as.numeric(fixes[[n]]$Time)-1)/4
  }

keydates=c(20070608,20070616,20070619,20070628)

for(t in 1:4)
  {
    b=fixes[[1]]
    I=which(b$Date==keydates[t] & b$Time=="00:00")
    a=unique(b$ID[I])
    for(j in 1:length(a)) 
    {
      c=b[b$ID==a[j],]
      if(t==1 & j==1) fixesshort=c else fixesshort=rbind(fixesshort,c)
    }
    
    b=fixes[[2]]
    I=which(b$Date==keydates[t] & b$Time=="00:00")
    a=unique(b$ID[I])
    for(j in 1:length(a)) 
    {
      c=b[b$ID==a[j],]
      if(t==1 & j==1) fixesshort_NT=c else fixesshort_NT=rbind(fixesshort_NT,c)
    }
  }

dates=data.frame(Date=seq(20070601,20070630.75,0.25),
                 CV=rep(NaN,120),MSLP=rep(NaN,120),CV_NT=rep(NaN,120),MSLP_NT=rep(NaN,120))
for(i in 1:120)
{
  I=which(fixesshort$Date2==dates$Date[i])
  if(length(I)>0){
    dates[i,2]=mean(fixesshort$CV[I])
    dates[i,3]=mean(fixesshort$MSLP[I])
  } 
  
  I=which(fixesshort_NT$Date2==dates$Date[i])
  if(length(I)>0){
    dates[i,4]=mean(fixesshort_NT$CV[I])
    dates[i,5]=mean(fixesshort_NT$MSLP[I])
  } 
}
plot(dates$Date,dates[,2],type="l",lwd=3,col="blue",
     ylim=c(0,3),xlab="Date",ylab="CV",main="Mean ECL curvature across all runs")
lines(dates$Date,dates[,4],lwd=3,col="red")
abline(v=keydates,col="lightgray")
legend("topleft",legend=c("D01","D02"),col=c("blue","red"),lwd=3)


