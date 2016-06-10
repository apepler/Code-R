rm(list=ls())

library(RNetCDF)
library(fields)
library(akima)
library("R.matlab")
readMat('~/Documents/Data/Useful.mat')->Useful
mask<-t(Useful$mask)
mask[is.na(mask)]=0
source('~/R/x86_64-pc-linux-gnu-library/2.15/gplots/R/rich.colors.R')

setwd("/srv/ccrc/data36/z3478332/WRF/output/ERAI_R2_ensemble_notopo/out")

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

n=1
for(day in 27:31)
  for(hour in c("00","06","12","18"))
    {
     dir1="/home/nfs/z3478332/output/outputUM_wrf_2007-06/ERAI_R2_ensemble/d01_p60/"
     dir2="/home/nfs/z3478332/output/outputUM_wrf_2007-06/ERAI_R2_ensemble_notopo/d01_p60/"
     events[[n]]=read.csv(paste(dir1,"ECLevents_200705",day,hour,".csv",sep=""))
     fixes[[n]]=read.csv(paste(dir1,"ECLfixes_200705",day,hour,".csv",sep=""))
     events_notopo[[n]]=read.csv(paste(dir2,"ECLevents_200705",day,hour,".csv",sep=""))
     fixes_notopo[[n]]=read.csv(paste(dir2,"ECLfixes_200705",day,hour,".csv",sep=""))
     
     n=n+1
}

events2=events[[1]]
fixes2=fixes[[1]]
events2_notopo=events_notopo[[1]]
fixes2_notopo=fixes_notopo[[1]]
for(i in 2:20)
{
  events2=rbind(events2,events[[i]])
  fixes2=rbind(fixes2,fixes[[i]])
  events2_notopo=rbind(events2_notopo,events_notopo[[i]])
  fixes2_notopo=rbind(fixes2_notopo,fixes_notopo[[i]])
}


####Start: All events/all fixes joined together in a density plot

pdf(file="D01_ECLs_eventCV.pdf",width=6,height=4)
makePDF(events2$CV2,events2_notopo$CV2,"Event maximum curvature for all lows in the ESB",labloc="topright")   
dev.off()
pdf(file="D01_ECLs_eventMSLP.pdf",width=6,height=4)
makePDF(events2$MSLP2,events2_notopo$MSLP2,"Event minimum MSLP for all lows in the ESB")
dev.off()
pdf(file="D01_ECLs_eventlength.pdf",width=6,height=4)
makePDF(events2$Length2,events2_notopo$Length2,"Event duration for all lows in the ESB",labloc="topright")
dev.off()
pdf(file="D01_ECLs_eventrad.pdf",width=6,height=4)
makePDF(events2$Rad2,events2_notopo$Rad2,"Event average radius for all lows in the ESB")
dev.off()
pdf(file="D01_ECLs_fixCV.pdf",width=6,height=4)
makePDF(fixes2$CV,fixes2_notopo$CV,"Curvature for all lows in the ESB",labloc="topright") 
dev.off()
pdf(file="D01_ECLs_fixMSLP.pdf",width=6,height=4)
makePDF(fixes2$MSLP,fixes2_notopo$MSLP,"MSLP for all lows in the ESB")
dev.off()

### Next, try to compare event CV/MSLP throughout time

fixes2$Date2=fixes2$Date+(as.numeric(fixes2$Time)-1)/4
fixes2_notopo$Date2=fixes2_notopo$Date+(as.numeric(fixes2_notopo$Time)-1)/4
dates=data.frame(Date=seq(20070601,20070630.75,0.25),
                 CV=rep(NaN,120),MSLP=rep(NaN,120),CV_NT=rep(NaN,120),MSLP_NT=rep(NaN,120))
for(i in 1:120)
{
  I=which(fixes2$Date2==dates$Date[i] & fixes2$Location==1)
  if(length(I)>0){
    dates[i,2]=mean(fixes2$CV[I])
    dates[i,3]=mean(fixes2$MSLP[I])
  } 
  
  I=which(fixes2_notopo$Date2==dates$Date[i] & fixes2_notopo$Location==1)
  if(length(I)>0){
    dates[i,4]=mean(fixes2_notopo$CV[I])
    dates[i,5]=mean(fixes2_notopo$MSLP[I])
  } 
}

pdf(file="D01_ECLs_dailyCV.pdf",width=7,height=4)
plot(dates$Date,dates[,2],type="l",lwd=3,col="blue",
     ylim=c(0,3),xlab="Date",ylab="CV",main="Mean ECL curvature across all runs")
lines(dates$Date,dates[,4],lwd=3,col="red")
legend("topleft",legend=c("Control","NoTopo"),col=c("blue","red"),lwd=3)
dev.off()
pdf(file="D01_ECLs_dailyMSLP.pdf",width=7,height=4)
plot(dates$Date,dates[,3],type="l",lwd=3,col="blue",
     ylim=c(970,1030),xlab="Date",ylab="MSLP",main="Mean ECL MSLP across all runs")
lines(dates$Date,dates[,5],lwd=3,col="red")
legend("bottomleft",legend=c("Control","NoTopo"),col=c("blue","red"),lwd=3)
dev.off()

### Next, look at changes in location

lat=seq(-45,-20,2.5)
lon=seq(140,170,2.5)

ECLloc<-ECLloc_notopo<-CVloc<-CVloc_notopo<-matrix(NaN,length(lon),length(lat))

for(i in 1:length(lon))
  for(j in 1:length(lat))
  {
    I=which(fixes2$Lat>=lat[j]-1.25 & fixes2$Lat<lat[j]+1.25 & 
              fixes2$Lon>=lon[i]-1.25 & fixes2$Lon<lon[i]+1.25)
    if(length(I)>0)
    {
      ECLloc[i,j]=length(I)
      CVloc[i,j]=mean(fixes2$CV[I])      
    }
    
    I=which(fixes2_notopo$Lat>=lat[j]-1.25 & fixes2_notopo$Lat<lat[j]+1.25 & 
              fixes2_notopo$Lon>=lon[i]-1.25 & fixes2_notopo$Lon<lon[i]+1.25)
    if(length(I)>0)
    {
      ECLloc_notopo[i,j]=length(I)
      CVloc_notopo[i,j]=mean(fixes2_notopo$CV[I])      
    }
  }

ECLloc[ECLloc>120]=120
ECLloc_notopo[ECLloc_notopo>120]=120
CVloc[CVloc>60]=60
CVloc_notopo[CVloc_notopo>60]=60

ECLpcchange = (ECLloc_notopo/ECLloc -1)*100
CVpcchange = (CVloc_notopo/CVloc -1)*100
cols=gray(seq(0.9,0.1,-0.15))
pdf(file="D01_ECLloc_control.pdf",width=5,height=6)
image.plot(lon,lat,ECLloc/20,xlab="",ylab="",col=cols,xlim=c(145,165),breaks=0:6,zlim=c(0,6),main="Average number of lows")
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
dev.off()
pdf(file="D01_ECLloc_notopo.pdf",width=5,height=6)
image.plot(lon,lat,ECLloc_notopo/20,xlab="",ylab="",col=cols,xlim=c(145,165),breaks=0:6,zlim=c(0,6),main="Average number of lows")
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
dev.off()
pdf(file="D01_ECLcv.pdf",width=5,height=6)
image.plot(lon,lat,CVloc,xlab="",ylab="",col=cols,xlim=c(145,165),breaks=seq(0,3,0.5),zlim=c(0,3),main="Average low CV")
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
dev.off()
pdf(file="D01_ECLcv_notopo.pdf",width=5,height=6)
image.plot(lon,lat,CVloc_notopo,xlab="",ylab="",col=cols,xlim=c(145,165),breaks=seq(0,3,0.5),zlim=c(0,3),main="Average low CV")
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
dev.off()

### Next, look at changes in specific events?


for(day in c(8,16,19,28))
{
  date=20070600+day
  eventmaxCV=matrix(NaN,20,2)
  
  for(i in 1:20)
{
  a=fixes[[i]]
  a$Date2=a$Date+(as.numeric(a$Time)-1)/4
  I=which(a$Date2==date)
  J=which((a$ID %in% unique(a$ID[I])) & a$Location==1)
  if(length(J)>0) eventmaxCV[i,1]=max(a$CV[J])
  
  a=fixes_notopo[[i]]
  a$Date2=a$Date+(as.numeric(a$Time)-1)/4
  I=which(a$Date2==date)
  J=which((a$ID %in% unique(a$ID[I])) & a$Location==1)
  if(length(J)>0) eventmaxCV[i,2]=max(a$CV[J])
}

pdf(file=paste("D01_eventmaxCV_",day,"Jun.pdf",sep=""),width=6,height=4)
makePDF(eventmaxCV[,1],eventmaxCV[,2],paste("Event maximum CV -",day,"June 2007"))
dev.off()
}



##Oooh, how many have an ECL?

fixes2$Date2=fixes2$Date+(as.numeric(fixes2$Time)-1)/4
fixes2_notopo$Date2=fixes2_notopo$Date+(as.numeric(fixes2_notopo$Time)-1)/4
dates=data.frame(Date=seq(20070601,20070630.75,0.25),
                 Count=rep(0,120),Count_nt=rep(0,120))
for(i in 1:120)
  for(j in 1:20)
{
    a=fixes[[j]]
    a$Date2=a$Date+(as.numeric(a$Time)-1)/4
    I=which(a$Date2==dates$Date[i] & a$Location==1 & a$CV>=1)
  if(length(I)>0) dates[i,2]=dates[i,2]+1
  
  a=fixes_notopo[[j]]
  a$Date2=a$Date+(as.numeric(a$Time)-1)/4
  I=which(a$Date2==dates$Date[i] & a$Location==1 & a$CV>=1)
  if(length(I)>0) dates[i,3]=dates[i,3]+1
}

dates[,2:3]=dates[,2:3]*5

pdf(file="D01_ECLs_dailyprob_cv1.pdf",width=7,height=4)
plot(dates$Date,dates[,2],type="l",lwd=3,col="blue",
     ylim=c(0,100),xlab="Date",ylab="Proportions (%)",main="Proportion of members with an ECL of CV>=1")
lines(dates$Date,dates[,3],lwd=3,col="red")
legend("topleft",legend=c("Control","NoTopo"),col=c("blue","red"),lwd=3)
dev.off()

###############
## Mean ECL track

##Tracks for ECLs - June 2007
library(RNetCDF)
f1<-open.nc("/srv/ccrc/data34/z3478332/ERAI/ERAI_mslp_2005-01_2010-12.nc")
lonP=var.get.nc(f1,'longitude')
latP=var.get.nc(f1,'latitude')
time=var.get.nc(f1,'time')
hh=time%%24
time=as.Date(time/24,origin="1900-01-01")
date=as.numeric(substr(time,1,4))*10000+as.numeric(substr(time,6,7))*100+
  as.numeric(substr(time,9,10))

source('~/Documents/R/filled.contour3.R')
source('~/Documents/R/filled.legend.R')
library("R.matlab")
readMat('~/Documents/Data/Useful.mat')->Useful
mask<-t(Useful$mask)
mask[is.na(mask)]=0
																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																		
keydates=c(20070608,20070616,20070619,20070628)

for(t in 1:4)
{
pdf(file=paste("ECL_meantrack_",keydates[t],"_d01.pdf",sep=""),width=5,height=5,pointsize=12)
par(mar=c(4,4,2,2))
I=which(date==keydates[t] & hh==0)
MSLP=var.get.nc(f1,"msl",c(1,1,I),c(length(lonP),length(latP),1),unpack=T)/100
#pdf(file=paste("Fig8_forposter_",keydates[t],".pdf",sep=""),width=9,height=10,pointsize=36)
filled.contour3(Useful$x,Useful$y,mask,col=c("white","lightgreen"),xlim=c(145,170),ylim=c(-45,-20))
contour(lonP,latP[length(latP):1],MSLP[,length(latP):1],levels=seq(992,1050,4),add=T,lwd=2,labcex=1) #19980808 = 29 (1 + 7*4)
  
for(i in 1:20)
  {
    b=fixes[[i]]
    I=which(b$Date==keydates[t] & b$Time=="00:00")
    a=unique(b$ID[I])
    for(j in 1:length(a)) 
    {
      c=b[b$ID==a[j],]
      lines(c$Lon,c$Lat,col=rgb(0,0,1,1/4))
      if(i==1 & j==1) fixesshort=c else fixesshort=rbind(fixesshort,c)
    }
    
    b=fixes_notopo[[i]]
    I=which(b$Date==keydates[t] & b$Time=="00:00")
    a=unique(b$ID[I])
    for(j in 1:length(a)) 
    {
      c=b[b$ID==a[j],]
      lines(c$Lon,c$Lat,col=rgb(1,0,0,1/4))
      if(i==1 & j==1) fixesshort_NT=c else fixesshort_NT=rbind(fixesshort_NT,c)
    }
  }

fixesshort$Date2=fixesshort$Date+(as.numeric(fixesshort$Time)-1)/4
fixesshort_NT$Date2=fixesshort_NT$Date+(as.numeric(fixesshort_NT$Time)-1)/4
dates=data.frame(Date=seq(keydates[t]-3,keydates[t]+3,0.25),count=rep(NaN,25),
                 lat=rep(NaN,25),lon=rep(NaN,25),count_NT=rep(NaN,25),lat_NT=rep(NaN,25),lon_NT=rep(NaN,25))

for(i in 1:length(dates[,1]))
{
  I=which(fixesshort$Date2==dates$Date[i])
  dates$count[i]=length(I)
  if(length(I)>10)
  {
    dates$lat[i]=mean(fixesshort$Lat[I])
    dates$lon[i]=mean(fixesshort$Lon[I])
  }
  
  I=which(fixesshort_NT$Date2==dates$Date[i])
  dates$count_NT[i]=length(I)
  if(length(I)>10)
  {
    dates$lat_NT[i]=mean(fixesshort_NT$Lat[I])
    dates$lon_NT[i]=mean(fixesshort_NT$Lon[I])
  }
}

lines(dates$lon,dates$lat,col="blue",lwd=3)
lines(dates$lon_NT,dates$lat_NT,col="red",lwd=3)
I=which(dates$Date==keydates[t])
points(dates$lon[I],dates$lat[I],col="blue",lwd=3,cex=2,pch=4)
points(dates$lon_NT[I],dates$lat_NT[I],col="red",lwd=3,cex=2,pch=4)
legend("topright",legend=c("Control","NoTopo"),
             col=c("blue","red"),lwd=3) 
dev.off()
}

###Daily CV v2

for(t in 1:4)
  for(i in 1:20)
  {
    b=fixes[[i]]
    I=which(b$Date==keydates[t] & b$Time=="00:00")
    a=unique(b$ID[I])
    for(j in 1:length(a)) 
    {
      c=b[b$ID==a[j],]
      if(i==1 & j==1 & t==1) fixesshort=c else fixesshort=rbind(fixesshort,c)
    }
    
    b=fixes_notopo[[i]]
    I=which(b$Date==keydates[t] & b$Time=="00:00")
    a=unique(b$ID[I])
    for(j in 1:length(a)) 
    {
      c=b[b$ID==a[j],]
      if(i==1 & j==1 & t==1) fixesshort_NT=c else fixesshort_NT=rbind(fixesshort_NT,c)
    }
  }
  
  fixesshort$Date2=fixesshort$Date+(as.numeric(fixesshort$Time)-1)/4
  fixesshort_NT$Date2=fixesshort_NT$Date+(as.numeric(fixesshort_NT$Time)-1)/4
  dates=data.frame(Date=seq(20070601,20070630.75,0.25),
                 CV=rep(NaN,120),MSLP=rep(NaN,120),CV_NT=rep(NaN,120),MSLP_NT=rep(NaN,120))
for(i in 1:120)
{
  I=which(fixesshort$Date2==dates$Date[i])
  if(length(I)>10){
    dates[i,2]=mean(fixesshort$CV[I])
    dates[i,3]=mean(fixesshort$MSLP[I])
  } 
  
  I=which(fixesshort_NT$Date2==dates$Date[i])
  if(length(I)>10){
    dates[i,4]=mean(fixesshort_NT$CV[I])
    dates[i,5]=mean(fixesshort_NT$MSLP[I])
  } 
}

pdf(file="D01_ECLs_dailyCV_events.pdf",width=7,height=4)
plot(dates$Date,dates[,2],type="l",lwd=3,col="blue",
     ylim=c(0,3),xlab="Date",ylab="CV",main="Mean ECL curvature across all runs")
lines(dates$Date,dates[,4],lwd=3,col="red")
abline(v=keydates,col="lightgray")
legend("topleft",legend=c("Control","NoTopo"),col=c("blue","red"),lwd=3)
dev.off()
pdf(file="D01_ECLs_dailyMSLP_events.pdf",width=7,height=4)
plot(dates$Date,dates[,3],type="l",lwd=3,col="blue",
     ylim=c(970,1030),xlab="Date",ylab="MSLP",main="Mean ECL MSLP across all runs")
lines(dates$Date,dates[,5],lwd=3,col="red")
abline(v=keydates,col="lightgray")
legend("bottomleft",legend=c("Control","NoTopo"),col=c("blue","red"),lwd=3)
dev.off()

##########
#####
# R2 only - look at distribution of max ECL intensity on 27-29 June 2007
  
rm(list=ls())
setwd("/srv/ccrc/data36/z3478332/WRF/output/ERAI_R2_ensemble_notopo/out")
fixes<- list()

event26=matrix(0,20,2)

n=1
for(day in 27:31)
  for(hour in c("00","06","12","18"))
  {
    dir1="/srv/ccrc/data34/z3478332/ECLtracks/outputUM_wrf_2007-06/ERAI_R2_ensemble/d02-as-d01_p240/"
    fixes[[n]]=read.csv(paste(dir1,"ECLfixes_200705",day,hour,"_d02.csv",sep=""))
    fixes[[n]]$Date2=as.POSIXct(paste(as.character(fixes[[n]]$Date),substr(fixes[[n]]$Time,1,2),sep=""),format="%Y%m%d%H",tz="GMT")
    
    I=which(fixes[[n]]$Date2>=as.POSIXct("2007062700",format="%Y%m%d%H",tz="GMT") &
              fixes[[n]]$Date2<=as.POSIXct("2007062918",format="%Y%m%d%H",tz="GMT") &
              fixes[[n]]$Location==1)
    
    event26[n,1]=max(fixes[[n]]$CV[I])
    event26[n,2]=min(fixes[[n]]$MSLP[I])
    
    n=n+1
  }
apply(event26,2,mean)