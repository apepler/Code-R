rm(list=ls())

library(RNetCDF)
library(fields)
library(akima)
library("R.matlab")
readMat('~/Documents/Data/Useful_ECL.mat')->Useful
readMat('~/Documents/Data/Useful.mat')->Useful2
mask<-t(Useful$mask)
mask[is.na(mask)]=0
source('~/R/x86_64-pc-linux-gnu-library/2.15/gplots/R/rich.colors.R')

load('~/Documents/Data/Useful_WRF.RData')

totrain_default_d01<-totrain_notopo_d01<-maxrain_default_d01<-maxrain_notopo_d01<-matrix(0,720,20)
totrain_default_d02<-totrain_notopo_d02<-maxrain_default_d02<-maxrain_notopo_d02<-matrix(0,720,20)

n=1
for(day in 27:31)
  for(hour in c("00","06","12","18"))
  {
    dir1=paste("/srv/ccrc/data36/z3478332/WRF/output/ERAI_R2_ensemble/out/200705",day,hour,"/",sep="")
    dir2=paste("/srv/ccrc/data36/z3478332/WRF/output/ERAI_R2_ensemble_notopo/out/200705",day,hour,"/",sep="")

    fname="wrfhrly_d01_2007-06-01_00:00:00"    
    a=open.nc(paste(dir1,fname,sep=""))
    latW=var.get.nc(a,"XLAT",c(1,1,1),c(215,144,1))
    lonW=var.get.nc(a,"XLONG",c(1,1,1),c(215,144,1))
    rain1=var.get.nc(a,"PREC_ACC_C")+var.get.nc(a,"PREC_ACC_NC")
    a=open.nc(paste(dir2,fname,sep=""))
    rain2=var.get.nc(a,"PREC_ACC_C")+var.get.nc(a,"PREC_ACC_NC")
    
    for(i in 1:720)
    {
      totrain_default_d01[i,n]=mean(rain1[,,i]*mask50a,na.rm=T)
      totrain_notopo_d01[i,n]=mean(rain2[,,i]*mask50a,na.rm=T)
      maxrain_default_d01[i,n]=max(rain1[,,i]*mask50a,na.rm=T)
      maxrain_notopo_d01[i,n]=max(rain2[,,i]*mask50a,na.rm=T)
    }
    
    fname="wrfhrly_d02_2007-06-01_00:00:00"    
    a=open.nc(paste(dir1,fname,sep=""))
    latW=var.get.nc(a,"XLAT",c(1,1,1),c(325,200,1))
    lonW=var.get.nc(a,"XLONG",c(1,1,1),c(325,200,1))
    rain1=var.get.nc(a,"PREC_ACC_C")+var.get.nc(a,"PREC_ACC_NC")
    a=open.nc(paste(dir2,fname,sep=""))
    rain2=var.get.nc(a,"PREC_ACC_C")+var.get.nc(a,"PREC_ACC_NC")
    
    for(i in 1:720)
    {
      totrain_default_d02[i,n]=mean(rain1[,,i]*mask10a,na.rm=T)
      totrain_notopo_d02[i,n]=mean(rain2[,,i]*mask10a,na.rm=T)
      maxrain_default_d02[i,n]=max(rain1[,,i]*mask10a,na.rm=T)
      maxrain_notopo_d02[i,n]=max(rain2[,,i]*mask10a,na.rm=T)
    }
    
    n=n+1
  }


setwd("/srv/ccrc/data36/z3478332/WRF/output/ERAI_R2_ensemble_notopo/out/")
load("WRF_ESB_dailyrain.RData")

data=cbind(apply(totrain_default_d01,2,sum),apply(totrain_notopo_d01,2,sum),apply(totrain_default_d02,2,sum),apply(totrain_notopo_d02,2,sum))
data2=cbind(apply(totrain_default_d01,2,max),apply(totrain_notopo_d01,2,max),apply(totrain_default_d02,2,max),apply(totrain_notopo_d02,2,max))
data3=cbind(apply(maxrain_default_d01,2,max),apply(maxrain_notopo_d01,2,max),apply(maxrain_default_d02,2,max),apply(maxrain_notopo_d02,2,max))

threshcount = function(data,thresh) {
  I=which(data>=thresh)
  return(length(I))
}
data4=cbind(apply(maxrain_default_d01,2,threshcount,25),apply(maxrain_notopo_d01,2,threshcount,25),apply(maxrain_default_d02,2,threshcount,25),apply(maxrain_notopo_d02,2,threshcount,25))

a=density(data2[,1],na.rm=T)
b=density(data2[,2],na.rm=T)

lims=range(data2[,1:2])
lims[1]=floor(lims[1]/10)*10
lims[2]=ceiling(lims[2]/10)*10

pdf(file="/srv/ccrc/data36/z3478332/WRF/output/ERAI_R2_ensemble_notopo/out/D01_maxESBdailyrain_notopo_allruns_v2.pdf",
    width=6,height=4)
plot(a,col=rgb(0,0,1,1/4),xlim=lims,ylim=range(0,a$y,b$y),
     xlab="Rain (mm)",ylab="Frequency",cex.main=1.2,
     main="Max daily June 2007 ESB rainfall for 20 runs")
polygon(a,col=rgb(0,0,1,1/4),density=-1)
polygon(b,col=rgb(1,0,0,1/4),density=-1)
legend("topright",legend=c("Control","NoTopo"),
       col=c(rgb(0,0,1,1/4),col=rgb(1,0,0,1/4),"red"),lwd=4,cex=1,bty="n")   
dev.off()

a=density(data2[,3],na.rm=T)
b=density(data2[,4],na.rm=T)

lims=range(data2[,3:4])
lims[1]=floor(lims[1]/10)*10
lims[2]=ceiling(lims[2]/10)*10

pdf(file="/srv/ccrc/data36/z3478332/WRF/output/ERAI_R2_ensemble_notopo/out/D02_maxESBdailyrain_notopo_allruns_v2.pdf",
    width=6,height=4)
plot(a,col=rgb(0,0,1,1/4),xlim=lims,ylim=range(0,a$y,b$y),
     xlab="Rain (mm)",ylab="Frequency",cex.main=1.2,
     main="Max daily June 2007 ESB rainfall for 20 runs")
polygon(a,col=rgb(0,0,1,1/4),density=-1)
polygon(b,col=rgb(1,0,0,1/4),density=-1)
legend("topright",legend=c("Control","NoTopo"),
       col=c(rgb(0,0,1,1/4),col=rgb(1,0,0,1/4),"red"),lwd=4,cex=1,bty="n")   
dev.off()



daily=cbind(apply(totrain_default_d01,1,mean),apply(totrain_notopo_d01,1,mean),apply(totrain_default_d02,1,mean),apply(totrain_notopo_d02,1,mean))

plot(seq(1:30),daily[,1],col=1,type="l",lwd=3,xlab="Day (June 2007)",ylab="Mean rain (mm)",main="Daily mean ESB rain during June 2007 (20 runs)")
for(i in 2:4) lines(seq(1:30),daily[,i],col=i,type="l",lwd=3)
legend("topright",c("D01 default","D01 notopo","D02 default","D02 notopo"),lwd=3,col=1:4,ncol=2,bty="n")
