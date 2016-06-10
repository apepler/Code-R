rm(list=ls())
setwd('/srv/ccrc/data36/z3478332/WRF/output/ERAI_R2_ensemble_notopo/out/')

library(RNetCDF)
library(fields)
library(akima)
library("R.matlab")
readMat('~/Documents/Data/Useful_ECL.mat')->Useful
readMat('~/Documents/Data/Useful.mat')->Useful2
mask<-t(Useful$mask)
mask[is.na(mask)]=0
source('~/R/x86_64-pc-linux-gnu-library/2.15/gplots/R/rich.colors.R')

# meanR=matrix(0,20,4)
# n=1
# for(day in 27:31)
#   for(hour in c("00","06","12","18"))
#   {
#     dir1=paste("/srv/ccrc/data36/z3478332/WRF/output/ERAI_R2_ensemble/out/200705",day,hour,"/",sep="")
#     dir2=paste("/srv/ccrc/data36/z3478332/WRF/output/ERAI_R2_ensemble_notopo/out/200705",day,hour,"/",sep="")
#     fname="wrfhrly_d01_2007-06-01_00:00:00"
#     
#     a=open.nc(paste(dir1,fname,sep=""))
#     latW=var.get.nc(a,"XLAT",c(1,1,1),c(215,144,1))
#     lonW=var.get.nc(a,"XLONG",c(1,1,1),c(215,144,1))
#     rain1=apply(var.get.nc(a,"PREC_ACC_C"),c(1,2),sum)+apply(var.get.nc(a,"PREC_ACC_NC"),c(1,2),sum)
#     a=open.nc(paste(dir2,fname,sep=""))
#     rain2=apply(var.get.nc(a,"PREC_ACC_C"),c(1,2),sum)+apply(var.get.nc(a,"PREC_ACC_NC"),c(1,2),sum)
#     
#     
#     lont=as.vector(lonW)
#     latt=as.vector(latW)
#     rr=as.vector(rain1)
#     rr[which(is.na(rr))]=0
#     rr2=interp(lont,latt,rr,Useful$x,Useful$y)
#     meanR[n,1]=mean(rr2$z*t(Useful$mask),na.rm=T)
#     
#     rr=as.vector(rain2)
#     rr[which(is.na(rr))]=0
#     rr2=interp(lont,latt,rr,Useful$x,Useful$y)
#     meanR[n,2]=mean(rr2$z*t(Useful$mask),na.rm=T)
#       
#     rr=as.vector(rain2-rain1)
#     rr[which(is.na(rr))]=0
#     rr2=interp(lont,latt,rr,Useful$x,Useful$y)
#     meanR[n,3]=mean(rr2$z*t(Useful$mask),na.rm=T) 
#     n=n+1
#   }
# 
# meanR[,4]=meanR[,2]-meanR[,1]
# write.csv(meanR,file="/srv/ccrc/data36/z3478332/WRF/output/ERAI_R2_ensemble_notopo/out/D01_totalESBrain.csv")

read.csv("/srv/ccrc/data36/z3478332/WRF/output/ERAI_R2_ensemble_notopo/out/D01_totalESBrain.csv")->meanR

##Make a fun density plot

a=density(meanR[,1],na.rm=T)
b=density(meanR[,2],na.rm=T)

lims=range(meanR[,1:2])
lims[1]=floor(lims[1]/10)*10
lims[2]=ceiling(lims[2]/10)*10

pdf(file="/srv/ccrc/data36/z3478332/WRF/output/ERAI_R2_ensemble_notopo/D01_totalESBrain_notopo_allruns.pdf",
    width=6,height=4)
plot(a,col=rgb(0,0,1,1/4),xlim=lims,ylim=range(0,a$y,b$y),
     xlab="Total rain (mm)",ylab="Frequency",cex.main=1.2,
     main="Total June 2007 ESB rainfall for 20 runs")
polygon(a,col=rgb(0,0,1,1/4),density=-1)
polygon(b,col=rgb(1,0,0,1/4),density=-1)
legend("topleft",legend=c("Control","NoTopo"),
       col=c(rgb(0,0,1,1/4),col=rgb(1,0,0,1/4),"red"),lwd=4,cex=1,bty="n")   
dev.off()





##########  Trying with daily prop > value

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
# dailyrain<-dailyrain_nt<-dailymaxrain<-dailymaxrain_nt<-array(0,c(215,144,30,20))
# 
# n=1
# for(day in 27:31)
#   for(hour in c("00","06","12","18"))
#   {
#     dir1=paste("/srv/ccrc/data36/z3478332/WRF/output/ERAI_R2_ensemble/out/200705",day,hour,"/",sep="")
#     dir2=paste("/srv/ccrc/data36/z3478332/WRF/output/ERAI_R2_ensemble_notopo/out/200705",day,hour,"/",sep="")
#     
#     fname="wrfhrly_d01_2007-06-01_00:00:00"    
#     a=open.nc(paste(dir1,fname,sep=""))
#     latW=var.get.nc(a,"XLAT",c(1,1,1),c(215,144,1))
#     lonW=var.get.nc(a,"XLONG",c(1,1,1),c(215,144,1))
#     rain1=var.get.nc(a,"PREC_ACC_C")+var.get.nc(a,"PREC_ACC_NC")
#     a=open.nc(paste(dir2,fname,sep=""))
#     rain2=var.get.nc(a,"PREC_ACC_C")+var.get.nc(a,"PREC_ACC_NC")
#     
#     for(i in 1:30)
#     {
#       dailyrain[,,i,n]=apply(rain1[,,((i-1)*24+1):(i*24)],c(1,2),sum)
#       dailyrain_nt[,,i,n]=apply(rain2[,,((i-1)*24+1):(i*24)],c(1,2),sum)
#       dailymaxrain[,,i,n]=apply(rain1[,,((i-1)*24+1):(i*24)],c(1,2),max)
#       dailymaxrain_nt[,,i,n]=apply(rain2[,,((i-1)*24+1):(i*24)],c(1,2),max)
#     }
#     n=n+1
#   }
# 
# save(dailyrain,dailymaxrain,dailyrain_nt,dailymaxrain_nt,file="/srv/ccrc/data36/z3478332/WRF/output/ERAI_R2_ensemble_notopo/WRFdailyrain_ensemble_d01.RData")

load('WRFdailyrain_ensemble_d01.RData')

rain25<-rain50<-rain100<-rain25_nt<-rain50_nt<-rain100_nt<-matrix(0,30,20)
ncells=length(which(is.na(mask50a)==0))

for(n in 1:20)
for(i in 1:30)
{
  I=which(dailyrain[,,i,n]*mask50a>=25)
  if(length(I)>0) rain25[i,n]=length(I)/ncells
  I=which(dailyrain[,,i,n]*mask50a>=50)
  if(length(I)>0) rain50[i,n]=length(I)/ncells
  I=which(dailyrain[,,i,n]*mask50a>=100)
  if(length(I)>0) rain100[i,n]=length(I)/ncells
  
  I=which(dailyrain_nt[,,i,n]*mask50a>=25)
  if(length(I)>0) rain25_nt[i,n]=length(I)/ncells
  I=which(dailyrain_nt[,,i,n]*mask50a>=50)
  if(length(I)>0) rain50_nt[i,n]=length(I)/ncells
  I=which(dailyrain_nt[,,i,n]*mask50a>=100)
  if(length(I)>0) rain100_nt[i,n]=length(I)/ncells
}

daily=cbind(apply(rain25,1,mean),apply(rain50,1,mean),apply(rain100,1,mean),
            apply(rain25_nt,1,mean),apply(rain50_nt,1,mean),apply(rain100_nt,1,mean))*100

pdf(file="/srv/ccrc/data36/z3478332/WRF/output/ERAI_R2_ensemble_notopo/D01_dailyrainthresh_notopo.pdf",
    width=6,height=4)
plot(seq(20070601,20070630),daily[,1],type="l",lwd=3,col=1,
     ylim=c(0,30),xlab="Date",ylab="Area (%)",main="Daily mean area of ESB exceeding rain thresholds")
for(i in 2:3) lines(seq(20070601,20070630),daily[,i],lwd=3,col=i)
for(i in 1:3) lines(seq(20070601,20070630),daily[,i+3],lwd=3,col=i,lty=2)
legend("topleft",legend=c("25mm","50mm","100mm"),col=1:3,lwd=3,bty="n")
legend("topright",legend=c("Control","NoTopo"),col=1,lwd=3,bty="n",lty=c(1,2))
dev.off()


a=density(apply(rain25,2,max)*100,na.rm=T)
b=density(apply(rain25_nt,2,max)*100,na.rm=T)
lims=range(apply(rain25,2,max),apply(rain25_nt,2,max))*100 

lims[1]=floor(lims[1]/5)*5
lims[2]=ceiling(lims[2]/5)*5

pdf(file="/srv/ccrc/data36/z3478332/WRF/output/ERAI_R2_ensemble_notopo/D01_maxarea25mm_notopo.pdf",
    width=6,height=4)
plot(a,col=rgb(0,0,1,1/4),xlim=lims,ylim=range(0,a$y,b$y),
     xlab="Rain area (%)",ylab="Frequency",cex.main=1.2,
     main=paste("Percentage area above 25 mm"))
polygon(a,col=rgb(0,0,1,1/4),density=-1)
polygon(b,col=rgb(1,0,0,1/4),density=-1)
legend("topleft",legend=c("Control","NoTopo"),
       col=c(rgb(0,0,1,1/4),col=rgb(1,0,0,1/4),"red"),lwd=4,cex=1,bty="n")   
dev.off()

a=density(apply(rain50,2,max)*100,na.rm=T)
b=density(apply(rain50_nt,2,max)*100,na.rm=T)
lims=range(apply(rain50,2,max),apply(rain50_nt,2,max))*100 

lims[1]=floor(lims[1]/5)*5
lims[2]=ceiling(lims[2]/5)*5

pdf(file="/srv/ccrc/data36/z3478332/WRF/output/ERAI_R2_ensemble_notopo/D01_maxarea50mm_notopo.pdf",
    width=6,height=4)
plot(a,col=rgb(0,0,1,1/4),xlim=lims,ylim=range(0,a$y,b$y),
     xlab="Rain area (%)",ylab="Frequency",cex.main=1.2,
     main=paste("Percentage area above 50 mm"))
polygon(a,col=rgb(0,0,1,1/4),density=-1)
polygon(b,col=rgb(1,0,0,1/4),density=-1)
legend("topright",legend=c("Control","NoTopo"),
       col=c(rgb(0,0,1,1/4),col=rgb(1,0,0,1/4),"red"),lwd=4,cex=1,bty="n")   
dev.off()

  a=density(apply(rain100,2,max)*100,na.rm=T)
  b=density(apply(rain100_nt,2,max)*100,na.rm=T)
  lims=range(apply(rain100,2,max),apply(rain100_nt,2,max))*100 

lims[1]=floor(lims[1]/5)*5
  lims[2]=ceiling(lims[2]/5)*5
  
  pdf(file="/srv/ccrc/data36/z3478332/WRF/output/ERAI_R2_ensemble_notopo/D01_maxarea100mm_notopo.pdf",
      width=6,height=4)
  plot(a,col=rgb(0,0,1,1/4),xlim=lims,ylim=range(0,a$y,b$y),
       xlab="Rain area (%)",ylab="Frequency",cex.main=1.2,
       main=paste("Percentage area above 100 mm"))
  polygon(a,col=rgb(0,0,1,1/4),density=-1)
  polygon(b,col=rgb(1,0,0,1/4),density=-1)
  legend("topright",legend=c("Control","NoTopo"),
         col=c(rgb(0,0,1,1/4),col=rgb(1,0,0,1/4),"red"),lwd=4,cex=1,bty="n")   
  dev.off()



for(day in c(6,7,8,15,19,26))
{
 a=density(rain25[day,]*100,na.rm=T)
 b=density(rain25_nt[day,]*100,na.rm=T)
 lims=range(rain25[day,],rain25_nt[day,])*100

lims[1]=floor(lims[1]/5)*5
lims[2]=ceiling(lims[2]/5)*5

pdf(file=paste("/srv/ccrc/data36/z3478332/WRF/output/ERAI_R2_ensemble_notopo/D01_Jun",day,"_area25mm_notopo.pdf",sep=""),
    width=6,height=4)
plot(a,col=rgb(0,0,1,1/4),xlim=lims,ylim=range(0,a$y,b$y),
     xlab="Rain area (%)",ylab="Frequency",cex.main=1.2,
     main=paste("Percentage area above 25 mm"))
polygon(a,col=rgb(0,0,1,1/4),density=-1)
polygon(b,col=rgb(1,0,0,1/4),density=-1)
legend("topright",legend=c("Control","NoTopo"),
       col=c(rgb(0,0,1,1/4),col=rgb(1,0,0,1/4),"red"),lwd=4,cex=1,bty="n")   
dev.off()
}


meanR<-meanR_nt<-matrix(0,30,20)
for(i in 1:30)
  for(n in 1:20)
  {
    meanR[i,n]=mean(dailyrain[,,i,n]*mask50a,na.rm=T)
    meanR_nt[i,n]=mean(dailyrain_nt[,,i,n]*mask50a,na.rm=T)
  }

signif<-matrix(0,30,4)
for(i in 1:30)
{
  a=ks.test(meanR[i,],meanR_nt[i,])
  signif[i,1]=a$p.value
  a=ks.test(rain25[i,],rain25_nt[i,])
  signif[i,2]=a$p.value
  a=ks.test(rain50[i,],rain50_nt[i,])
  signif[i,3]=a$p.value
  a=ks.test(rain100[i,],rain100_nt[i,])
  signif[i,4]=a$p.value
}

pdf(file="/srv/ccrc/data36/z3478332/WRF/output/ERAI_R2_ensemble_notopo/D01_dailyrainsig_notopo.pdf",
    width=6,height=4)
plot(seq(20070601,20070630),apply(meanR,1,mean),type="l",lwd=3,col="blue",
     ylim=c(0,30),xlab="Date",ylab="Rain (mm)",main="Daily mean ESB rainfall")
lines(seq(20070601,20070630),apply(meanR_nt,1,mean),lwd=3,col="red")
I=which(signif[,1]<=0.05)
a=seq(20070601,20070630)
points(a[I],rep(0,length(I)),pch=4,col="black",lwd=2)
legend("topright",legend=c("Control","NoTopo"),col=c("blue","red"),lwd=3,bty="n")
dev.off()

####################

rm(list=ls())
setwd('/srv/ccrc/data36/z3478332/WRF/output/ERAI_R2_ensemble_notopo/out/')

library(RNetCDF)
library(fields)
library(akima)
library("R.matlab")
readMat('~/Documents/Data/Useful_ECL.mat')->Useful
readMat('~/Documents/Data/Useful.mat')->Useful2
mask<-t(Useful$mask)
mask[is.na(mask)]=0
source('~/R/x86_64-pc-linux-gnu-library/2.15/gplots/R/rich.colors.R')

meanR=matrix(0,20,4)
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
    rain1=apply(var.get.nc(a,"PREC_ACC_C"),c(1,2),sum)+apply(var.get.nc(a,"PREC_ACC_NC"),c(1,2),sum)
    a=open.nc(paste(dir2,fname,sep=""))
    rain2=apply(var.get.nc(a,"PREC_ACC_C"),c(1,2),sum)+apply(var.get.nc(a,"PREC_ACC_NC"),c(1,2),sum)
    
    meanR[n,1]=mean(rain1,na.rm=T)
    meanR[n,2]=mean(rain2,na.rm=T)
    meanR[n,3]=mean(rain2-rain1,na.rm=T) 
    n=n+1
  }

meanR[,4]=meanR[,2]-meanR[,1]

a=density(meanR[,1],na.rm=T)
b=density(meanR[,2],na.rm=T)

lims=range(meanR[,1:2])
lims[1]=floor(lims[1]/10)*10
lims[2]=ceiling(lims[2]/10)*10

pdf(file="/srv/ccrc/data36/z3478332/WRF/output/ERAI_R2_ensemble_notopo/D01_totald01rain_notopo_allruns.pdf",
    width=6,height=4)
plot(a,col=rgb(0,0,1,1/4),xlim=lims,ylim=range(0,a$y,b$y),
     xlab="Total rain (mm)",ylab="Frequency",cex.main=1.2,
     main="Total June 2007 ESB rainfall for 20 runs")
polygon(a,col=rgb(0,0,1,1/4),density=-1)
polygon(b,col=rgb(1,0,0,1/4),density=-1)
legend("topleft",legend=c("Control","NoTopo"),
       col=c(rgb(0,0,1,1/4),col=rgb(1,0,0,1/4),"red"),lwd=4,cex=1,bty="n")   
dev.off()
