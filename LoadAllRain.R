###Still need to do average cool-season and warm-season rain across versions...
##AWAP
rm(list=ls())
source('~/Documents/R/CMIP5_play.R')
setwd("~/Documents/Data/CMIP5")
library(RNetCDF)
library(fields)
library("R.matlab")

load('~/Documents/Data/AWAPrain.RData')
yearsA=seq(1900,2012)
I=which(yearsA>=1950 & yearsA<=2005)
coolAWAP=apply(cool[,,I],c(1,2),mean)
coolAWAPa=apply(cool[,,I],c(1,2),median)
yearsB=seq(1901,2012)
I=which(yearsB>=1950 & yearsB<=2005) #So has the year ending in 1950 - year ending in 2005
warmAWAP=apply(warm[,,I],c(1,2),mean) #So has the year ending in 1950 - year ending in 2005
warmAWAPa=apply(warm[,,I],c(1,2),median) #So has the year ending in 1950 - year ending in 2005
rm(cool,warm)

##CMIP5
read.csv("Modellist.csv")->list
list=list[c(-24,-14),]
coolA<-warmA<-array(0,c(144,73,37))
for(i in 1:37)
{
  a=CMIP_rain_v2(as.character(list[i,1]),as.character(list[i,2]),as.character(list[i,3]),"N")
  coolA[,,i]=unlist(a[2])
  warmA[,,i]=unlist(a[3])
}
coolCMIP5=apply(coolA,c(1,2),mean,na.rm=T)
warmCMIP5=apply(warmA,c(1,2),mean,na.rm=T)
coolCMIP5a=apply(coolA,c(1,2),median,na.rm=T)
warmCMIP5a=apply(warmA,c(1,2),median,na.rm=T)

###CMIP3

a=23
b=47
yy=seq(1850,2012)
timeREF=matrix(0,length(yy)*12,3)
n=1
for(i in 1:163)
  for(j in 1:12)
  {
    timeREF[n,1]=yy[i]
    timeREF[n,2]=j
    timeREF[n,3]=yy[i]*100+j
    n=n+1
  }

name1=c("mpi-echam5","csiromk3","miroc3.2","cgcm3.1")
date1=c(186001,187101,185001,185001)
date2=c(210012,200012,200012,200012)
coolA<-warmA<-array(NaN,c(144,73,4))

for(n in 1:4)
{
  Rfile=paste("/srv/ccrc/data34/z3478332/CMIP3/precip/prcp_20c3m_",name1[n],"_2.5deg.nc",sep="")
  f1=open.nc(Rfile)
  lat=var.get.nc(f1,"lat")
  lon=var.get.nc(f1,"lon")
  rain=var.get.nc(f1,"pr")
  time=timeREF[(timeREF[,3]>=date1[n] & timeREF[,3]<=date2[n]),1:2]
  I=which(time[,1]>=1950 & time[,1]<=2005)
  rain=rain[,,I]
  time=time[I,]
  
  I=which(time[,2] %in% c(1,3,5,7,8,10,12))
  rain[,,I]=rain[,,I]*31*60*60*24
  I=which(time[,2] %in% c(4,6,9,11))
  rain[,,I]=rain[,,I]*30*60*60*24
  I=which(time[,2]==2)
  rain[,,I]=rain[,,I]*28*60*60*24
  
  years=seq(min(time[,1]),max(time[,1]))
  coolR<-warmR<-array(0,c(144,73,length(years)))
  for(i in 1:length(years))
  {
    I=which(time[,1]==years[i] & time[,2]>=5 & time[,2]<=10)
    coolR[,,i]=apply(rain[,,I],c(1,2),sum)
    I=which((time[,1]==years[i] & time[,2]<=4) | (time[,1]==years[i]-1 & time[,2]>=11))
    warmR[,,i]=apply(rain[,,I],c(1,2),sum)
  }
  coolA[,,n]=apply(coolR,c(1,2),mean)
  warmA[,,n]=apply(warmR,c(1,2),mean)
}
coolCMIP3=apply(coolA,c(1,2),mean,na.rm=T)
warmCMIP3=apply(warmA,c(1,2),mean,na.rm=T)
coolCMIP3a=apply(coolA,c(1,2),median,na.rm=T)
warmCMIP3a=apply(warmA,c(1,2),median,na.rm=T)

##WRF50
name2=c("echam5","csiromk3","miroc3.2","cccma")
coolA<-warmA<-array(NaN,dim=c(151,101,12))
count=1
for(r in 1:3)
  for(n in 1:4)
  {
    W_pfile=paste("/srv/ccrc/data36/z3478332/CMIP3/WRF_R",r,"_",name2[n],"_precip.nc",sep="")    
    f1=open.nc(W_pfile)
    lat50=var.get.nc(f1,"lat50")
    lon50=var.get.nc(f1,"lon50")
    rain=var.get.nc(f1,"prcp50")
    time=var.get.nc(f1,"time")
    year=floor(time/100)
    month=time%%100
    time=cbind(year,month)
    
    #Finally, let's do the two correlation plots
    years=seq(min(time[,1]),max(time[,1]))
    coolR<-warmR<-array(0,c(151,101,length(years)))
    for(i in 1:length(years))
    {
      I=which(time[,1]==years[i] & time[,2]>=5 & time[,2]<=10)
      coolR[,,i]=apply(rain[,,I],c(1,2),sum)
      I=which((time[,1]==years[i] & time[,2]<=4) | (time[,1]==years[i]-1 & time[,2]>=11))
      warmR[,,i]=apply(rain[,,I],c(1,2),sum)
    }
    
    coolA[,,count]=apply(coolR,c(1,2),mean)
    warmA[,,count]=apply(warmR,c(1,2),mean)
    count=count+1    
  }

cool50=apply(coolA,c(1,2),mean,na.rm=T)
warm50=apply(warmA,c(1,2),mean,na.rm=T)
cool50a=apply(coolA,c(1,2),median,na.rm=T)
warm50a=apply(warmA,c(1,2),median,na.rm=T)

##WRF50N
coolA<-warmA<-coolA1<-warmA1<-array(NaN,dim=c(151,101,3))
count=1
for(r in 1:3)
  {
    W_pfile=paste("/srv/ccrc/data34/z3478332/WRF/WRF_R",r,"_ncep1_precip.nc",sep="")    
    f1=open.nc(W_pfile)
    lat50=var.get.nc(f1,"lat50")
    lon50=var.get.nc(f1,"lon50")
    rain=var.get.nc(f1,"prcp50")
    time=var.get.nc(f1,"time")
    year=floor(time/100)
    month=time%%100
    time=cbind(year,month)
    
    #Finally, let's do the two correlation plots
    years=seq(min(time[,1]),max(time[,1]))
    coolM=matrix(0,length(years),2)
    coolR<-warmR<-array(0,c(151,101,length(years)))
    for(i in 1:length(years))
    {
      I=which(time[,1]==years[i] & time[,2]>=5 & time[,2]<=10)
      coolR[,,i]=apply(rain[,,I],c(1,2),sum)
      I=which((time[,1]==years[i] & time[,2]<=4) | (time[,1]==years[i]-1 & time[,2]>=11))
      warmR[,,i]=apply(rain[,,I],c(1,2),sum)
    }
    
    I=which(years>=1950 & years<=2005)
    coolA[,,count]=apply(coolR[,,I],c(1,2),mean)
    warmA[,,count]=apply(warmR[,,I],c(1,2),mean)
    I=which(years>=1990 & years<=2009)
    coolA1[,,count]=apply(coolR[,,I],c(1,2),mean)
    warmA1[,,count]=apply(warmR[,,I],c(1,2),mean)
    count=count+1    
  }

cool50N=apply(coolA,c(1,2),mean,na.rm=T)
warm50N=apply(warmA,c(1,2),mean,na.rm=T)
cool50Na=apply(coolA,c(1,2),median,na.rm=T)
warm50Na=apply(warmA,c(1,2),median,na.rm=T)
cool50N_99=apply(coolA1,c(1,2),mean,na.rm=T)
warm50N_99=apply(warmA1,c(1,2),mean,na.rm=T)
cool50Na_99=apply(coolA1,c(1,2),median,na.rm=T)
warm50Na_99=apply(warmA1,c(1,2),median,na.rm=T)

##WRF10
coolA<-warmA<-array(NaN,dim=c(301,251,12))
count=1
for(r in 1:3)
  for(n in 1:4)
  {
    W_pfile=paste("/srv/ccrc/data36/z3478332/CMIP3/WRF10_R",r,"_",name2[n],"_precip.nc",sep="")
    f1=open.nc(W_pfile)
    lat10=var.get.nc(f1,"lat10")
    lon10=var.get.nc(f1,"lon10")
    rain=var.get.nc(f1,"prcp10")
    time=var.get.nc(f1,"time")
    year=floor(time/100)
    month=time%%100
    time=cbind(year,month)
    
    #Finally, let's do the two correlation plots
    years=seq(min(time[,1]),max(time[,1]))
    years=cbind(seq(min(time[,1]),max(time[,1])),matrix(0,(length(years)),2))
    coolM=matrix(0,length(years[,1]),2)
    coolR<-warmR<-array(0,c(301,251,length(years[,1])))
    for(i in 1:length(years[,1]))
    {
      I=which(time[,1]==years[i,1] & time[,2]>=5 & time[,2]<=10)
      coolR[,,i]=apply(rain[,,I],c(1,2),sum)
      I=which((time[,1]==years[i,1] & time[,2]<=4) | (time[,1]==years[i,1]-1 & time[,2]>=11))
      warmR[,,i]=apply(rain[,,I],c(1,2),sum)
    }
    
    coolA[,,count]=apply(coolR,c(1,2),mean)
    warmA[,,count]=apply(warmR,c(1,2),mean)
    count=count+1    
  }

cool10=apply(coolA,c(1,2),mean,na.rm=T)
warm10=apply(warmA,c(1,2),mean,na.rm=T)
cool10a=apply(coolA,c(1,2),median,na.rm=T)
warm10a=apply(warmA,c(1,2),median,na.rm=T)
save.image(file="RAIN.RData")

library(abind)
AWAP=abind(coolAWAP,warmAWAP,coolAWAPa,warmAWAPa,along=3)
CMIP5=abind(coolCMIP5,warmCMIP5,coolCMIP5a,warmCMIP5a,along=3)
CMIP3=abind(coolCMIP3,warmCMIP3,coolCMIP3a,warmCMIP3a,along=3)
WRF50=abind(cool50,warm50,cool50a,warm50a,along=3)
WRF50N=abind(cool50N,warm50N,cool50Na,warm50Na,along=3)
WRF50N_9009=abind(cool50N_99,warm50N_99,cool50Na_99,warm50Na_99,along=3)
WRF10=abind(cool10,warm10,cool10a,warm10a,along=3)
names=c("cool_mean","warm_mean","cool_median","warm_median")

AWAP250=array(NaN,dim(CMIP5))
AWAP50=array(NaN,dim(WRF50))
AWAP10=array(NaN,dim(WRF10))

library(fields)
library("R.matlab")
readMat('~/Documents/Data/Useful.mat')->Useful
mask<-t(Useful$mask)
mask[is.na(mask)]=0
source('~/Documents/R/filled.contour3.R')
for(k in 1:4)
{
  a=interp.surface.grid(list(x=Useful$x,y=Useful$y,z=t(AWAP[,,k])),list(x=lon,y=lat))
  AWAP250[,,k]=a$z
  a=interp.surface.grid(list(x=Useful$x,y=Useful$y,z=t(AWAP[,,k])),list(x=lon50,y=lat50))
  AWAP50[,,k]=a$z
  a=interp.surface.grid(list(x=Useful$x,y=Useful$y,z=t(AWAP[,,k])),list(x=lon10,y=lat10))
  AWAP10[,,k]=a$z
}

for(k in 1:4)
{
  fname=paste("Rain_AWAP250_",names[k],"_smaller.tiff",sep="")
  tiff(file=fname, height=550, width=500)
  par(plt = c(0.1,0.85,0.15,0.85),las = 1,cex.axis = 1)
  image(lon,lat,AWAP250[,,k],xlab="",ylab="",breaks=c(0,10,25,50,100,200,300,400,600,800,1200,1600,2000,10000),
        col=rich.colors(13),zlim=c(0,10000),xlim=c(130,155),ylim=c(-45,-23))
  contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
  par(new = "TRUE",plt = c(0.85,0.9,0.15,0.85),las = 1,cex.axis = 1)
  image.plot(legend.only=TRUE,breaks=seq(0,13),lab.breaks=c("0","10","25","50","100","200","300","400","600","800","1200","1600","2000",""),col=rich.colors(13),zlim=c(0,13))
  dev.off()
  
  fname=paste("Rain_CMIP5_",names[k],"_smaller.tiff",sep="")
  tiff(file=fname, height=550, width=500)
  par(plt = c(0.1,0.85,0.15,0.85),las = 1,cex.axis = 1)
  image(lon,lat,CMIP5[,,k],xlab="",ylab="",breaks=c(0,10,25,50,100,200,300,400,600,800,1200,1600,2000,10000),
        col=rich.colors(13),zlim=c(0,10000),xlim=c(130,155),ylim=c(-45,-23))
  contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
  par(new = "TRUE",plt = c(0.85,0.9,0.15,0.85),las = 1,cex.axis = 1)
  image.plot(legend.only=TRUE,breaks=seq(0,13),lab.breaks=c("0","10","25","50","100","200","300","400","600","800","1200","1600","2000",""),col=rich.colors(13),zlim=c(0,13))
  dev.off()
  
  fname=paste("Rain_CMIP3_",names[k],"_smaller.tiff",sep="")
  tiff(file=fname, height=550, width=500)
  par(plt = c(0.1,0.85,0.15,0.85),las = 1,cex.axis = 1)
  image(lon,lat,CMIP3[,,k],xlab="",ylab="",breaks=c(0,10,25,50,100,200,300,400,600,800,1200,1600,2000,10000),
        col=rich.colors(13),zlim=c(0,10000),xlim=c(130,155),ylim=c(-45,-23))
  contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
  par(new = "TRUE",plt = c(0.85,0.9,0.15,0.85),las = 1,cex.axis = 1)
  image.plot(legend.only=TRUE,breaks=seq(0,13),lab.breaks=c("0","10","25","50","100","200","300","400","600","800","1200","1600","2000",""),col=rich.colors(13),zlim=c(0,13))
  dev.off()
  
  fname=paste("Rain_AWAP50_",names[k],"_smaller.tiff",sep="")
  tiff(file=fname, height=550, width=500)
  par(plt = c(0.1,0.85,0.15,0.85),las = 1,cex.axis = 1)
  image(lon50,lat50,AWAP50[,,k],xlab="",ylab="",breaks=c(0,10,25,50,100,200,300,400,600,800,1200,1600,2000,10000),
        col=rich.colors(13),zlim=c(0,10000),xlim=c(130,155),ylim=c(-45,-23))
  contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
  par(new = "TRUE",plt = c(0.85,0.9,0.15,0.85),las = 1,cex.axis = 1)
  image.plot(legend.only=TRUE,breaks=seq(0,13),lab.breaks=c("0","10","25","50","100","200","300","400","600","800","1200","1600","2000",""),col=rich.colors(13),zlim=c(0,13))
  dev.off()
  
  fname=paste("Rain_WRF50_",names[k],"_smaller.tiff",sep="")
  tiff(file=fname, height=550, width=500)
  par(plt = c(0.1,0.85,0.15,0.85),las = 1,cex.axis = 1)
  image(lon50,lat50,WRF50[,,k],xlab="",ylab="",breaks=c(0,10,25,50,100,200,300,400,600,800,1200,1600,2000,10000),
        col=rich.colors(13),zlim=c(0,10000),xlim=c(130,155),ylim=c(-45,-23))
  contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
  par(new = "TRUE",plt = c(0.85,0.9,0.15,0.85),las = 1,cex.axis = 1)
  image.plot(legend.only=TRUE,breaks=seq(0,13),lab.breaks=c("0","10","25","50","100","200","300","400","600","800","1200","1600","2000",""),col=rich.colors(13),zlim=c(0,13))
  dev.off()
  
  fname=paste("Rain_WRF50_NCEP5005_",names[k],"_smaller.tiff",sep="")
  tiff(file=fname, height=550, width=500)
  par(plt = c(0.1,0.85,0.15,0.85),las = 1,cex.axis = 1)
  image(lon50,lat50,WRF50N[,,k],xlab="",ylab="",breaks=c(0,10,25,50,100,200,300,400,600,800,1200,1600,2000,10000),
        col=rich.colors(13),zlim=c(0,10000),xlim=c(130,155),ylim=c(-45,-23))
  contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
  par(new = "TRUE",plt = c(0.85,0.9,0.15,0.85),las = 1,cex.axis = 1)
  image.plot(legend.only=TRUE,breaks=seq(0,13),lab.breaks=c("0","10","25","50","100","200","300","400","600","800","1200","1600","2000",""),col=rich.colors(13),zlim=c(0,13))
  dev.off()
  
  fname=paste("Rain_WRF50_NCEP9009_",names[k],"_smaller.tiff",sep="")
  tiff(file=fname, height=550, width=500)
  par(plt = c(0.1,0.85,0.15,0.85),las = 1,cex.axis = 1)
  image(lon50,lat50,WRF50N_9009[,,k],xlab="",ylab="",breaks=c(0,10,25,50,100,200,300,400,600,800,1200,1600,2000,10000),
        col=rich.colors(13),zlim=c(0,10000),xlim=c(130,155),ylim=c(-45,-23))
  contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
  par(new = "TRUE",plt = c(0.85,0.9,0.15,0.85),las = 1,cex.axis = 1)
  image.plot(legend.only=TRUE,breaks=seq(0,13),lab.breaks=c("0","10","25","50","100","200","300","400","600","800","1200","1600","2000",""),col=rich.colors(13),zlim=c(0,13))
  dev.off()
  
  fname=paste("Rain_AWAP10_",names[k],"_smaller.tiff",sep="")
  tiff(file=fname, height=550, width=500)
  par(plt = c(0.1,0.85,0.15,0.85),las = 1,cex.axis = 1)
  image(lon10,lat10,AWAP10[,,k],xlab="",ylab="",breaks=c(0,10,25,50,100,200,300,400,600,800,1200,1600,2000,10000),
        col=rich.colors(13),zlim=c(0,10000),xlim=c(130,155),ylim=c(-45,-23))
  contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
  par(new = "TRUE",plt = c(0.85,0.9,0.15,0.85),las = 1,cex.axis = 1)
  image.plot(legend.only=TRUE,breaks=seq(0,13),lab.breaks=c("0","10","25","50","100","200","300","400","600","800","1200","1600","2000",""),col=rich.colors(13),zlim=c(0,13))
  dev.off()
  
  fname=paste("Rain_WRF10_",names[k],"_smaller.tiff",sep="")
  tiff(file=fname, height=550, width=500)
  par(plt = c(0.1,0.85,0.15,0.85),las = 1,cex.axis = 1)
  image(lon10,lat10,WRF10[,,k],xlab="",ylab="",breaks=c(0,10,25,50,100,200,300,400,600,800,1200,1600,2000,10000),
        col=rich.colors(13),zlim=c(0,10000),xlim=c(130,155),ylim=c(-45,-23))
  contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
  par(new = "TRUE",plt = c(0.85,0.9,0.15,0.85),las = 1,cex.axis = 1)
  image.plot(legend.only=TRUE,breaks=seq(0,13),lab.breaks=c("0","10","25","50","100","200","300","400","600","800","1200","1600","2000",""),col=rich.colors(13),zlim=c(0,13))
  dev.off()
  
  fname=paste("Rain_AWAP_",names[k],"_smaller.tiff",sep="")
  tiff(file=fname, height=550, width=500)
  par(plt = c(0.1,0.85,0.15,0.85),las = 1,cex.axis = 1)
  filled.contour3(Useful$x,Useful$y,t(AWAP[,,k]*Useful$mask),xlab="",
                  c(0,10,25,50,100,200,300,400,600,800,1200,1600,2000,10000),col=rich.colors(13)
                  ,zlim=c(0,10000),xlim=c(130,155),ylim=c(-45,-23))
  par(xpd = NA)
  contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
  par(new = "TRUE",plt = c(0.85,0.9,0.15,0.85),las = 1,cex.axis = 1)
  image.plot(legend.only=TRUE,breaks=seq(0,13),lab.breaks=c("0","10","25","50","100","200","300","400","600","800","1200","1600","2000",""),col=rich.colors(13),zlim=c(0,13))
  dev.off()
}

save(names,AWAP,AWAP250,AWAP50,AWAP10,CMIP5,CMIP3,WRF50,WRF50N,WRF50N_9009,WRF10,file="RAIN.RData")
