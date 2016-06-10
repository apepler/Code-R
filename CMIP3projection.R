##############33

rm(list=ls())
setwd("~/output")

now<-future<-matrix(0,12,12)
rownames(now)<-rownames(future)<-rep("aaa",12)

cmip=c("echam5","csiromk3","miroc","cccma")
wrf=c("R1","R2","R3")


n=1
for(i in 1:4)
  for(j in 1:3)
  {
    rownames(now)[n]<-rownames(future)[n]<-paste(cmip[i],wrf[j])
    dir=paste("outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv1/",sep="")
    data=read.csv(paste(dir,"ECLevents_umelb_",cmip[i],"_wrf",wrf[j],"_rad2cv1.csv",sep=""))
    mm=floor((data$Date1%%10000)/100)
    for(m in 1:12)
    {
      I=which(mm==m)
      now[n,m]=length(I)/20
    }
    data=read.csv(paste(dir,"ECLevents_umelb_",cmip[i],"_wrf",wrf[j],"_rad2cv1_6079.csv",sep=""))
    mm=floor((data$Date1%%10000)/100)
    for(m in 1:12)
    {
      I=which(mm==m)
      future[n,m]=length(I)/20
    }
    
    n=n+1
  }

now2=now*20
plot(seq(1:12),now2[1,],ylim=c(0,90),xlab="Month",ylab="Number of ECLs")
for(i in 2:12) lines(seq(1:12),now2[i,])
colnames(now2)=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
boxplot(now2/20,main="Number of ECLs p.a.",ylim=c(0,4))

change=((future/now)-1)*100
colnames(change)=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")

boxplot(change,main="% Change in ECL numbers between 1990-2009 and 2060-2079")
abline(h=0,col="red")

bymodel=matrix(0,4,12)
rownames(bymodel)=cmip
for(i in 1:4) bymodel[i,]=apply(change[((i-1)*3+1):(i*3),],2,mean)

bywrf=matrix(0,3,12)
rownames(bywrf)=wrf
for(i in 1:3) bywrf[i,]=apply(change[seq(i,12,3),],2,mean)

colnames(bywrf)<-colnames(bymodel)<-colnames(change)

plot(seq(1:12),bymodel[1,],col=1,type="p",pch=19,ylim=c(-60,110),
     xlab="Month",ylab="",main="% Change in ECL numbers between 1990-2009 and 2060-2079 - by model")
for(i in 2:4) points(seq(1:12),bymodel[i,],col=i,pch=19)
abline(h=0,col="red")
legend("topright",pch=19,col=1:4,legend=cmip)

plot(seq(1:12),bywrf[1,],col=1,type="p",pch=19,ylim=c(-60,110),
     xlab="Month",ylab="",main="% Change in ECL numbers between 1990-2009 and 2060-2079 - by WRF version")
for(i in 2:3) points(seq(1:12),bywrf[i,],col=i,pch=19)
abline(h=0,col="red")
legend("topright",pch=19,col=1:3,legend=wrf)

################ Version 2 

rm(list=ls())
setwd("~/output")

now<-future<-matrix(0,12,12)
rownames(now)<-rownames(future)<-rep("aaa",12)

cmip=c("echam5","csiromk3","miroc","cccma")
wrf=c("R1","R2","R3")


n=1
for(i in 1:4)
  for(j in 1:3)
  {
    rownames(now)[n]<-rownames(future)[n]<-paste(cmip[i],wrf[j])
    dir=paste("outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad5cv25/",sep="")
    data=read.csv(paste(dir,"ECLevents_umelb_",cmip[i],"_wrf",wrf[j],"_rad5cv25_9009.csv",sep=""))
    mm=floor((data$Date1%%10000)/100)
    for(m in 1:12)
    {
      I=which(mm==m)
      now[n,m]=length(I)/20
    }
    data=read.csv(paste(dir,"ECLevents_umelb_",cmip[i],"_wrf",wrf[j],"_rad5cv25_6079.csv",sep=""))
    mm=floor((data$Date1%%10000)/100)
    for(m in 1:12)
    {
      I=which(mm==m)
      future[n,m]=length(I)/20
    }
    
    n=n+1
  }

now2=now*20
plot(seq(1:12),now2[1,],ylim=c(0,110),type="l",xlab="Month",ylab="Number of ECLs")
for(i in 2:12) lines(seq(1:12),now2[i,])
colours=c("red","orange","green","blue")
#for(i in 1:4) lines(seq(1:12),apply(now2[((i-1)*3+1):(i*3),],2,mean),lwd=4,col=colours[i])
#legend("bottomright",cmip,ncol=4,lwd=4,col=colours)
for(i in 1:3) lines(seq(1:12),apply(now2[seq(i,12,3),],2,mean),lwd=4,col=colours[i])
legend("bottomright",wrf,ncol=4,lwd=4,col=colours)


change=((future/now)-1)*100
colnames(change)=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")

boxplot(change,main="% Change in ECL numbers between 1990-2009 and 2060-2079",ylim=c(-50,200))
abline(h=0,col="red")

bymodel=matrix(0,4,12)
rownames(bymodel)=cmip
for(i in 1:4) bymodel[i,]=apply(change[((i-1)*3+1):(i*3),],2,mean)

bywrf=matrix(0,3,12)
rownames(bywrf)=wrf
for(i in 1:3) bywrf[i,]=apply(change[seq(i,12,3),],2,mean)

colnames(bywrf)<-colnames(bymodel)<-colnames(change)

plot(seq(1:12),bymodel[1,],col=1,type="p",pch=19,ylim=c(-60,110),
     xlab="Month",ylab="",main="% Change in ECL numbers between 1990-2009 and 2060-2079 - by model")
for(i in 2:4) points(seq(1:12),bymodel[i,],col=i,pch=19)
abline(h=0,col="red")
legend("topright",pch=19,col=1:4,legend=cmip)

plot(seq(1:12),bywrf[1,],col=1,type="p",pch=19,ylim=c(-60,110),
     xlab="Month",ylab="",main="% Change in ECL numbers between 1990-2009 and 2060-2079 - by WRF version")
for(i in 2:3) points(seq(1:12),bywrf[i,],col=i,pch=19)
abline(h=0,col="red")
legend("topright",pch=19,col=1:3,legend=wrf)


##########Now, restrict to the event CV threshold that gives ~ 5 or ~ 2 p.a. in historical

rm(list=ls())
setwd("~/output")

now<-future<-array(0,c(12,12,5))
dimnames(now)[[1]]<-dimnames(future)[[1]]<-rep("aaa",12)

cmip=c("echam5","csiromk3","miroc","cccma")
wrf=c("R1","R2","R3")
thresh=c(20,10,5,2,1)
n=1
for(t in 1:5)
{
  n=1
for(i in 1:4)
  for(j in 1:3)
  {
    dimnames(now)[[1]][n]<-dimnames(future)[[1]][n]<-paste(cmip[i],wrf[j])
    dir=paste("outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad5cv25/",sep="")
    data=read.csv(paste(dir,"ECLevents_umelb_",cmip[i],"_wrf",wrf[j],"_rad5cv25_9009.csv",sep=""))
    
    b=order(data$CV2,decreasing=T)
    thresh2=data$CV2[b[20*thresh[t]]]
    data2=data[data$CV2>=thresh2,]
    
    mm=floor((data2$Date1%%10000)/100)
    for(m in 1:12)
    {
      I=which(mm==m)
      now[n,m,t]=length(I)
    }
    
    data=read.csv(paste(dir,"ECLevents_umelb_",cmip[i],"_wrf",wrf[j],"_rad5cv25_6079.csv",sep=""))
    data2=data[data$CV2>=thresh2,]
    mm=floor((data2$Date1%%10000)/100)
    for(m in 1:12)
    {
      I=which(mm==m)
      future[n,m,t]=length(I)
    }
    
    n=n+1
  }
}

dimnames(now)[[2]]=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
lims=c(4,2,1,0.6,0.4)
for(i in 1:5)
{
  plot(seq(1:12),now[1,,i]/20,ylim=c(0,lims[i]),xlab="Month",ylab="Number of ECLs p.a.",type="l")
  for(m in 2:12) lines(seq(1:12),now[m,,i]/20)
  lines(seq(1:12),apply(now[,,i],2,mean)/20,lwd=3,col="red")
  #boxplot(now[,,i]/20,main="Number of ECLs p.a.")
}

coolprop=matrix(0,12,5)
rownames(coolprop)<-dimnames(now)[[1]]
colnames(coolprop)<-paste(thresh,"pa")
for(i in 1:5) coolprop[,i]=(apply(now[,5:10,i],1,sum)/apply(now[,,i],1,sum))
boxplot(coolprop*100,main="Proportion of ECLs during cool season")
abline(h=50,col="red")


change=((future/now)-1)*100
dimnames(change)[[2]]=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")

for(i in 1:5)
{
boxplot(change[,,i],main=paste("% Change in ECLs for thresh of",thresh[i],"events p.a. between 1990-2009 and 2060-2079"),ylim=c(-200,200))
abline(h=0,col="red")
}

nowY=apply(now,c(1,3),sum)
futureY=apply(future,c(1,3),sum)
changeY=((futureY/nowY)-1)*100
dimnames(changeY)[[2]]=thresh
boxplot(changeY,main="% Change in ECL events p.a. between 1990-2009 and 2060-2079",cex.main=1,ylim=c(-50,50),xlab="Threshold: ECLs p.a.")
abline(h=0,col="red")


season=c("MAM","JJA","SON","DJF")
library(abind)
nowS1=abind(now,now[,1:2,],along=2)
futureS1=abind(future,future[,1:2,],along=2)

changeS2<-array(0,c(12,4,5))
dimnames(changeS2)[[1]]=dimnames(change[[1]])
dimnames(changeS2)[[2]]=season
dimnames(changeS2)[[3]]=thresh
nowS2<-futureS2<-changeS2
for(i in 1:4)
{
  nowS2[,i,]=apply(nowS1[,(i*3):(i*3+2),],c(1,3),sum)
  futureS2[,i,]=apply(futureS1[,(i*3):(i*3+2),],c(1,3),sum)
  changeS2[,i,]=((futureS2[,i,]/nowS2[,i,])-1)*100
  boxplot(changeS2[,i,],main=paste("% Change in",season[i],"ECL events p.a. between 1990-2009 and 2060-2079"),cex.main=1,ylim=c(-100,100),xlab="Threshold: ECLs p.a.")
  abline(h=0,col="red")
}


nowS3=apply(now[,c(1:4,11:12),],c(1,3),sum)
futureS3=apply(future[,c(1:4,11:12),],c(1,3),sum)
changeS3=((futureS3/nowS3)-1)*100
boxplot(changeS3,main=paste("% Change in warm season ECL events p.a. between 1990-2009 and 2060-2079"),cex.main=1,ylim=c(-100,100),xlab="Threshold: ECLs p.a.")
abline(h=0,col="red")

nowSA=apply(nowS2,c(2,3),mean)
futureSA=apply(futureS2,c(2,3),mean)
changeSA=((futureSA/nowSA)-1)*100

meannow=apply(now,c(2,3),mean)/20
meanfuture=apply(future,c(2,3),mean)/20
change2=(meanfuture/meannow-1)*100
change3=apply(change,c(2,3),mean,na.rm=T)

plot(seq(1:12),meannow[,1],type="l",col=1,lwd=2,xlab="Month",ylab="Average ECLs p.a.",ylim=c(0,3))
for(i in 2:5) lines(seq(1:12),meannow[,i],col=i,lwd=2)
for(i in 1:5) lines(seq(1:12),meanfuture[,i],col=i,lwd=2,lty=2)
abline(h=0,col="grey")
legend("topright",legend=paste("Thresh:",thresh,"events p.a."),col=1:5,lwd=2)

plot(seq(1:12),change2[,1],type="l",col=1,lwd=2,xlab="Month",ylab="% change",ylim=c(-150,150))
for(i in 2:5) lines(seq(1:12),change2[,i],col=i,lwd=2)
legend("topright",legend=paste("Thresh:",thresh,"events p.a."),col=1:5,lwd=2)
abline(h=0,col="grey")



###Oooh, try to make fun density plots?

rm(list=ls())
setwd("~/output")
library(RNetCDF)
library(fields)
library(akima)
library("R.matlab")
readMat('~/Documents/Data/Useful.mat')->Useful
mask<-t(Useful$mask)
mask[is.na(mask)]=0
source('~/R/x86_64-pc-linux-gnu-library/2.15/gplots/R/rich.colors.R')

cmip=c("echam5","csiromk3","miroc","cccma")
wrf=c("R1","R2","R3")

for(thresh in c(10,5,2))
{
n=1
  for(i in 1:4)
    for(j in 1:3)
    {
      dir=paste("outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad5cv25/",sep="")
      data=read.csv(paste(dir,"ECLevents_umelb_",cmip[i],"_wrf",wrf[j],"_rad5cv25_9009.csv",sep=""))     
      
      b=order(data$CV2,decreasing=T)
      thresh2=data$CV2[b[20*thresh]]
      data1=data[data$CV2>=thresh2,]
      if(n==1) now=data[data$CV2>=thresh2,] else now=rbind(now,data[data$CV2>=thresh2,])   
      
      data=read.csv(paste(dir,"ECLevents_umelb_",cmip[i],"_wrf",wrf[j],"_rad5cv25_6079.csv",sep=""))
      b=order(data$CV2,decreasing=T)
      thresh2=data$CV2[b[20*thresh]]    
      if(n==1) future=data[data$CV2>=thresh2,]  else future=rbind(now,data[data$CV2>=thresh2,] )   
      data2=data[data$CV2>=thresh2,]
      
      a=density(data1$CV2,na.rm=T)
      b=density(data2$CV2,na.rm=T)
      
      lims=c(0,3)
      
      pdf(file=paste("~/Documents/ECLs/WRFCMIP/ECL_CVchange_pdf_",cmip[i],"_",wrf[j],"_thresh",thresh,".pdf",sep=""),width=5,height=4)

      plot(a,col=rgb(0,0,1,1/4),xlim=lims,ylim=range(0,a$y,b$y),
           xlab="Intensity",ylab="Frequency",cex.main=1.2,
           main="ECL intensity distribution")
      polygon(a,col=rgb(0,0,1,1/4),density=-1)
      polygon(b,col=rgb(1,0,0,1/4),density=-1)
      legend("topright",legend=c("1990-2009","2060-2079"),
             col=c(rgb(0,0,1,1/4),rgb(1,0,0,1/4)),lwd=4,cex=1,bty="n")   
      dev.off()
      
      n=n+1
    }

a=density(now$CV2,na.rm=T)
b=density(future$CV2,na.rm=T)

pdf(file=paste("~/Documents/ECLs/WRFCMIP/ECL_CVchange_pdf_thresh",thresh,".pdf",sep=""),width=5,height=4)
plot(a,col=rgb(0,0,1,1/4),xlim=lims,ylim=range(0,a$y,b$y),
     xlab="Intensity",ylab="Frequency",cex.main=1.2,
     main="ECL intensity distribution")
polygon(a,col=rgb(0,0,1,1/4),density=-1)
polygon(b,col=rgb(1,0,0,1/4),density=-1)
legend("topright",legend=c("1990-2009","2060-2079"),
       col=c(rgb(0,0,1,1/4),rgb(1,0,0,1/4)),lwd=4,cex=1,bty="n")   
dev.off()
}




